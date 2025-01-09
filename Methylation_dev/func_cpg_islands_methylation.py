import pandas as pd
import pyranges as pr # type: ignore
from config import CONFIG, PATHS, logger
import numpy as np
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import os
from tqdm import tqdm
import pysam
from typing import Dict, List, Tuple, Any, Union
import pyBigWig
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing

os.chdir("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Methylation_dev")

def load_cpg_islands(cpg_file: str) -> pr.PyRanges:
    """Load CpG islands from BED file into PyRanges object
    
    Format:
    chr1    3531624    3531843    611    CpG:    27
    """
    # Read the BED file
    cpg_df = pd.read_csv(cpg_file, sep='\t', header=None,
                        names=['Chromosome', 'Start', 'End', 'ID', 'Type', 'cpg_count'])
    
    # Convert CpG count to numeric (it's the last column)
    cpg_df['cpg_count'] = pd.to_numeric(cpg_df['cpg_count'])
    
    # Debug info
    logger.debug(f"Loaded CpG islands data shape: {cpg_df.shape}")
    logger.debug(f"Columns: {cpg_df.columns.tolist()}")
    logger.debug(f"First few rows:\n{cpg_df.head()}")
    
    # Convert to PyRanges
    cpg_islands = pr.PyRanges(cpg_df)
    logger.info(f"Loaded {len(cpg_islands)} CpG islands")
    return cpg_islands

def calculate_methylation_levels_parallel(region_df: pd.DataFrame,
                                        medip_dir: str,
                                        cell_type: str,
                                        genome_fasta: str,
                                        n_processes: int = None) -> pd.DataFrame:
    """Calculate methylation levels for genomic regions in parallel."""
    
    if n_processes is None:
        n_processes = max(1, multiprocessing.cpu_count() - 1)

    # Split data into chunks
    chunk_size = max(1, len(region_df) // (n_processes * 4))
    chunks = np.array_split(region_df, len(region_df) // chunk_size + 1)

    # Process chunks in parallel
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = [
            executor.submit(calculate_methylation_for_regions,
                          chunk, medip_dir, cell_type, genome_fasta)
            for chunk in chunks
        ]
        
        results = []
        for future in tqdm(futures, desc="Calculating methylation levels"):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                logger.error(f"Error processing chunk: {str(e)}")
                continue

    if not results:
        raise RuntimeError("No valid results obtained from parallel processing")
    
    return pd.concat(results, ignore_index=True)

def calculate_methylation_for_regions(region_df: pd.DataFrame,
                                    medip_dir: str,
                                    cell_type: str,
                                    genome_fasta: str) -> pd.DataFrame:
    """Calculate methylation for a set of genomic regions."""
    
    # Map cell types to file prefixes
    prefix_map = {'NEU': 'N', 'NSC': 'PP'}
    prefix = prefix_map.get(cell_type)
    if not prefix:
        raise ValueError(f"Unknown cell type: {cell_type}")

    # Get replicate file paths
    ip_files = [os.path.join(medip_dir, f"Medip_{prefix}_output_r{i}.bw") 
                for i in range(1, 4)]
    input_files = [os.path.join(medip_dir, f"Medip_{prefix}_input_r{i}.bw") 
                  for i in range(1, 4)]

    # Check files exist
    for f in ip_files + input_files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Missing required file: {f}")

    results = []
    with pysam.FastaFile(genome_fasta) as fasta:
        for _, region in region_df.iterrows():
            try:
                methylation, cpg_count = calculate_region_methylation(
                    region, ip_files, input_files, fasta
                )
                results.append({
                    'methylation': methylation,
                    'cpg_count': cpg_count
                })
            except Exception as e:
                logger.debug(f"Error processing region {region['chr']}:{region['start']}-{region['end']}: {str(e)}")
                results.append({
                    'methylation': 0,
                    'cpg_count': 0
                })

    return pd.DataFrame(results)

def calculate_region_methylation(region: pd.Series,
                               ip_files: List[str],
                               input_files: List[str],
                               fasta: pysam.FastaFile) -> Tuple[float, int]:
    """Calculate methylation with better handling of missing data"""
    
    replicate_values = []
    
    # Get sequence and validate CpG content
    sequence = fasta.fetch(region['chr'], int(region['start']), int(region['end']))
    cpg_count = sequence.upper().count('CG')
    
    if cpg_count == 0:
        logger.warning(f"No CpGs in region {region['chr']}:{region['start']}-{region['end']}")
        return None, cpg_count  # Mark as missing rather than zero

    valid_replicates = 0
    for ip_file, input_file in zip(ip_files, input_files):
        with pyBigWig.open(ip_file) as bw_ip, pyBigWig.open(input_file) as bw_input:
            try:
                ip_values = bw_ip.values(region['chr'], int(region['start']), int(region['end']))
                input_values = bw_input.values(region['chr'], int(region['start']), int(region['end']))
                
                if ip_values and input_values:
                    valid_pairs = [(ip, inp) for ip, inp in zip(ip_values, input_values)
                                 if ip is not None and inp is not None]
                    
                    if valid_pairs:
                        ip_vals, input_vals = zip(*valid_pairs)
                        methylation = calculate_normalized_methylation(ip_vals, input_vals, sequence)
                        if methylation is not None:
                            replicate_values.append(methylation)
                            valid_replicates += 1
                        
            except Exception as e:
                logger.debug(f"Error processing replicate: {str(e)}")
                continue

    # Require minimum number of valid replicates
    MIN_REPLICATES = 2
    if valid_replicates >= MIN_REPLICATES:
        return np.mean(replicate_values), cpg_count
    else:
        logger.warning(f"Insufficient valid replicates ({valid_replicates}) for region")
        return None, cpg_count

def calculate_normalized_methylation(ip_values: List[float],
                                  input_values: List[float],
                                  sequence: str) -> float:
    """Calculate methylation with better handling of edge cases"""
    
    # 1. Signal validation
    if not ip_values or not input_values:
        logger.warning("Missing IP or input values")
        return None  # Better to mark as missing than assume zero
        
    # 2. Calculate mean signals with minimum threshold
    MIN_SIGNAL = 0.001  # Prevent division by very small numbers
    ip_mean = max(MIN_SIGNAL, np.mean(ip_values))
    input_mean = max(MIN_SIGNAL, np.mean(input_values))
    
    # 3. Calculate enrichment
    enrichment = ip_mean / input_mean
    
    # 4. Normalize by CpG density with minimum threshold
    cpg_count = sequence.upper().count('CG')
    region_length = len(sequence)
    MIN_DENSITY = 0.001  # Minimum density threshold
    cpg_density = max(MIN_DENSITY, cpg_count / region_length)
    
    normalized_enrichment = enrichment / cpg_density
    
    # 5. Convert to methylation percentage with baseline
    # Assuming minimum ~20% methylation in mammalian cells
    BASE_METHYLATION = 20
    methylation = BASE_METHYLATION + min(80, normalized_enrichment * 15)
    
    return methylation

def analyze_cpg_methylation(merged_df: pd.DataFrame,
                          cpg_islands: pr.PyRanges,
                          medip_dir: str,
                          cell_type: str,
                          genome_fasta: str,
                          n_processes: int = None) -> pd.DataFrame:
    """Calculate methylation levels for CpG islands and analyze their relationship with gene expression
    
    Args:
        merged_df: DataFrame containing gene info and expression data
        cpg_islands: PyRanges object containing CpG island coordinates
        medip_dir: Directory containing MeDIP bigWig files
        cell_type: Cell type code ('NSC' or 'NEU')
        genome_fasta: Path to genome FASTA file
        n_processes: Number of parallel processes
    
    Returns:
        DataFrame with CpG island methylation levels and analysis results
    """
    # Create PyRanges objects for promoters and gene bodies
    promoters = pr.PyRanges(
        pd.DataFrame({
            'Chromosome': merged_df['chr'],
            'Start': merged_df['promoter_start'],
            'End': merged_df['promoter_end'],
            'gene_id': merged_df['gene_id'],
            'binding_type': merged_df['binding_type'],
            'expression_status': merged_df['expression_status'],
            'region_type': 'promoter'
        })
    )
    
    gene_bodies = pr.PyRanges(
        pd.DataFrame({
            'Chromosome': merged_df['chr'],
            'Start': merged_df['gene_body_start'],
            'End': merged_df['gene_body_end'],
            'gene_id': merged_df['gene_id'],
            'binding_type': merged_df['binding_type'],
            'expression_status': merged_df['expression_status'],
            'region_type': 'gene_body'
        })
    )
    
    # Find overlaps between CpG islands and gene regions
    promoter_cpgs = promoters.join(cpg_islands)
    gene_body_cpgs = gene_bodies.join(cpg_islands)
    
    # Combine overlaps
    all_cpgs = pd.concat([
        promoter_cpgs.as_df(),
        gene_body_cpgs.as_df()
    ])
    
    if all_cpgs.empty:
        logger.warning("No CpG islands found overlapping with genes")
        return pd.DataFrame()
    
    # Calculate methylation levels
    cpg_regions = pd.DataFrame({
        'chr': all_cpgs['Chromosome'],
        'start': all_cpgs['Start'],
        'end': all_cpgs['End'],
        'gene_id': all_cpgs['gene_id'],
        'binding_type': all_cpgs['binding_type'],
        'expression_status': all_cpgs['expression_status'],
        'region_type': all_cpgs['region_type'],
        'cpg_count': all_cpgs['cpg_count']
    })
    
    # Calculate methylation levels in parallel
    methylation_results = calculate_methylation_levels_parallel(
        cpg_regions,
        medip_dir,
        cell_type,
        genome_fasta,
        n_processes
    )
    
    # Merge results
    result_df = cpg_regions.merge(methylation_results,
                                left_index=True,
                                right_index=True)
    
    return result_df

def plot_cpg_methylation_by_category(data: pd.DataFrame, cell_type: str):
    """Create plots analyzing methylation levels in CpG islands by binding category and location
    
    Args:
        data: DataFrame containing methylation data for CpG islands with columns:
            - binding_type: type of binding ('both', 'exo_only', etc.)
            - region_type: location ('promoter' or 'gene_body')
            - expression_status: gene expression status
            - methylation: methylation percentage
        cell_type: Cell type code ('NSC' or 'NEU')
    """
    # Create output directory
    output_dir = os.path.join('plots', 'cpg_methylation', cell_type)
    os.makedirs(output_dir, exist_ok=True)
    
    # Set plot style
    plt.rcParams.update({
        'figure.figsize': (15, 6),
        'axes.grid': True,
        'grid.alpha': 0.3,
        'axes.spines.top': False,
        'axes.spines.right': False
    })
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Define binding categories
    binding_categories = {
        'exo_enriched': data[data['binding_type'].isin(['exo_only', 'both'])],
        'no_exo_enriched': data[~data['binding_type'].isin(['exo_only', 'both'])]
    }
    
    # Plot for promoter CpG islands
    plot_data = []
    for cat_name, cat_data in binding_categories.items():
        promoter_data = cat_data[cat_data['region_type'] == 'promoter']
        if not promoter_data.empty:
            for status in ['not_deregulated', 'upregulated', 'downregulated']:
                status_data = promoter_data[promoter_data['expression_status'] == status]
                if not status_data.empty:
                    plot_data.append({
                        'Category': cat_name,
                        'Expression': status,
                        'Methylation': status_data['methylation'].mean(),
                        'Count': len(status_data)
                    })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create bar plot for promoter CpGs
    sns.barplot(data=plot_df,
                x='Category',
                y='Methylation',
                hue='Expression',
                ax=ax1)
    
    ax1.set_title(f'{cell_type} - CpG Island Methylation in Promoters')
    ax1.set_xlabel('Binding Category')
    ax1.set_ylabel('Methylation Level (%)')
    
    # Add sample sizes
    for i, cat in enumerate(plot_df['Category'].unique()):
        cat_data = plot_df[plot_df['Category'] == cat]
        ax1.text(i, ax1.get_ylim()[1], f'n={cat_data["Count"].sum()}',
                horizontalalignment='center', verticalalignment='bottom')
    
    # Plot for gene body CpG islands
    plot_data = []
    for cat_name, cat_data in binding_categories.items():
        body_data = cat_data[cat_data['region_type'] == 'gene_body']
        if not body_data.empty:
            for status in ['not_deregulated', 'upregulated', 'downregulated']:
                status_data = body_data[body_data['expression_status'] == status]
                if not status_data.empty:
                    plot_data.append({
                        'Category': cat_name,
                        'Expression': status,
                        'Methylation': status_data['methylation'].mean(),
                        'Count': len(status_data)
                    })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create bar plot for gene body CpGs
    sns.barplot(data=plot_df,
                x='Category',
                y='Methylation',
                hue='Expression',
                ax=ax2)
    
    ax2.set_title(f'{cell_type} - CpG Island Methylation in Gene Bodies')
    ax2.set_xlabel('Binding Category')
    ax2.set_ylabel('Methylation Level (%)')
    
    # Add sample sizes
    for i, cat in enumerate(plot_df['Category'].unique()):
        cat_data = plot_df[plot_df['Category'] == cat]
        ax2.text(i, ax2.get_ylim()[1], f'n={cat_data["Count"].sum()}',
                horizontalalignment='center', verticalalignment='bottom')
    
    # Adjust layout and save
    plt.suptitle(f'{cell_type} - CpG Island Methylation Analysis by Binding Category', y=1.05)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{cell_type}_cpg_methylation_by_category.pdf'),
                bbox_inches='tight')
    plt.close()
    
    return plot_df

