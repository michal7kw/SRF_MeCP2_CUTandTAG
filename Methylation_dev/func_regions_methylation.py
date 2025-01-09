import pandas as pd
import numpy as np
import pyranges as pr
import pysam
import pyBigWig
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Any
import os
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from config import CONFIG, PATHS, logger
import multiprocessing
from functools import partial

os.chdir("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Methylation_dev")

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

def calculate_region_methylation(region: pd.Series,
                               ip_files: List[str],
                               input_files: List[str],
                               fasta: pysam.FastaFile) -> Tuple[float, int]:
    """Calculate methylation with better handling of missing data"""
    
    replicate_values = []
    
    # Get sequence and validate CpG content
    try:
        sequence = fasta.fetch(region['chr'], int(region['start']), int(region['end']))
        cpg_count = sequence.upper().count('CG')
    except Exception as e:
        logger.debug(f"Error fetching sequence: {str(e)}")
        return None, 0
    
    if cpg_count == 0:
        logger.debug(f"No CpGs in region {region['chr']}:{region['start']}-{region['end']}")
        return None, 0

    valid_replicates = 0
    for ip_file, input_file in zip(ip_files, input_files):
        try:
            with pyBigWig.open(ip_file) as bw_ip, pyBigWig.open(input_file) as bw_input:
                ip_values = bw_ip.values(region['chr'], int(region['start']), int(region['end']))
                input_values = bw_input.values(region['chr'], int(region['start']), int(region['end']))
                
                if ip_values and input_values:
                    valid_pairs = [(ip, inp) for ip, inp in zip(ip_values, input_values)
                                 if ip is not None and inp is not None and inp > 0]
                    
                    if valid_pairs:
                        ip_vals, input_vals = zip(*valid_pairs)
                        ip_mean = np.mean(ip_vals)
                        input_mean = np.mean(input_vals)
                        
                        # Calculate enrichment and normalize by CpG density
                        enrichment = ip_mean / input_mean
                        region_length = region['end'] - region['start']
                        cpg_density = (cpg_count / region_length) * 1000  # CpGs per kb
                        
                        # Convert to methylation percentage
                        methylation = min(100, max(0, 20 + enrichment * cpg_density))
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
        logger.debug(f"Insufficient valid replicates ({valid_replicates}) for region")
        return None, cpg_count

def plot_peak_methylation_by_category(data: pd.DataFrame, cell_type: str):
    """Create plots analyzing methylation levels in peaks by binding category and location
    
    Args:
        data: DataFrame containing methylation data for peaks with columns:
            - binding_type: type of binding ('both', 'exo_only', etc.)
            - binding_place: location of binding ('promoter' or 'gene_body')
            - expression_status: gene expression status
            - methylation: methylation percentage
        cell_type: Cell type code ('NSC' or 'NEU')
    """
    # Create output directory if it doesn't exist
    output_dir = os.path.join('plots', 'peak_methylation', cell_type)
    os.makedirs(output_dir, exist_ok=True)
    
    # Set default style parameters
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
    
    # Plot for promoter regions
    plot_data = []
    for cat_name, cat_data in binding_categories.items():
        promoter_data = cat_data[cat_data['binding_place'] == 'promoter']
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
    
    # Create bar plot for promoter
    sns.barplot(data=plot_df, 
                x='Category', 
                y='Methylation',
                hue='Expression',
                ax=ax1)
    
    ax1.set_title(f'{cell_type} - Peak Methylation in Promoters')
    ax1.set_xlabel('Binding Category')
    ax1.set_ylabel('Methylation Level (%)')
    
    # Add sample size annotations for promoter
    for i, cat in enumerate(plot_df['Category'].unique()):
        cat_data = plot_df[plot_df['Category'] == cat]
        ax1.text(i, ax1.get_ylim()[1], f'n={cat_data["Count"].sum()}',
                horizontalalignment='center', verticalalignment='bottom')
    
    # Plot for gene body regions
    plot_data = []
    for cat_name, cat_data in binding_categories.items():
        body_data = cat_data[cat_data['binding_place'] == 'gene_body']
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
    
    # Create bar plot for gene body
    sns.barplot(data=plot_df, 
                x='Category', 
                y='Methylation',
                hue='Expression',
                ax=ax2)
    
    ax2.set_title(f'{cell_type} - Peak Methylation in Gene Bodies')
    ax2.set_xlabel('Binding Category')
    ax2.set_ylabel('Methylation Level (%)')
    
    # Add sample size annotations for gene body
    for i, cat in enumerate(plot_df['Category'].unique()):
        cat_data = plot_df[plot_df['Category'] == cat]
        ax2.text(i, ax2.get_ylim()[1], f'n={cat_data["Count"].sum()}',
                horizontalalignment='center', verticalalignment='bottom')
    
    # Adjust layout and save
    plt.suptitle(f'{cell_type} - Peak Methylation Analysis by Binding Category', y=1.05)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{cell_type}_peak_methylation_by_category.pdf'),
                bbox_inches='tight')
    plt.close()
    
    return plot_df

def analyze_peak_methylation(merged_df: pd.DataFrame,
                           medip_dir: str,
                           cell_type: str,
                           genome_fasta: str,
                           n_processes: int = None) -> pd.DataFrame:
    """Calculate methylation levels for MeCP2 binding peaks and analyze their impact on gene expression"""
    
    # Filter for genes with MeCP2 binding
    bound_genes = merged_df[merged_df['mecp2_bound']].copy()
    logger.info(f"Analyzing methylation for {len(bound_genes)} MeCP2-bound genes")
    
    # Determine binding place for each peak
    bound_genes['binding_place'] = bound_genes.apply(
        lambda row: 'promoter' if (
            row['peak_start'] >= row['promoter_start'] and 
            row['peak_end'] <= row['promoter_end']
        ) else 'gene_body',
        axis=1
    )
    
    # Create regions DataFrame for methylation calculation
    peak_regions = pd.DataFrame({
        'chr': bound_genes['chr'],
        'start': bound_genes['peak_start'],
        'end': bound_genes['peak_end'],
        'gene_id': bound_genes['gene_id'],
        'binding_type': bound_genes['binding_type'],
        'binding_place': bound_genes['binding_place'],
        'expression_status': bound_genes['expression_status'],
        'log2FoldChange': bound_genes['log2FoldChange']
    }).dropna(subset=['start', 'end'])  # Remove any rows with NaN coordinates
    
    # Convert coordinates to integers
    peak_regions['start'] = peak_regions['start'].astype(int)
    peak_regions['end'] = peak_regions['end'].astype(int)
    
    logger.info(f"Processing {len(peak_regions)} peaks")
    logger.info("\nPeak distribution:")
    logger.info(peak_regions.groupby(['binding_type', 'binding_place']).size())
    
    # Calculate methylation levels
    methylation_results = calculate_methylation_levels_parallel(
        peak_regions,
        medip_dir,
        cell_type,
        genome_fasta,
        n_processes
    )
    
    # Merge results back
    result_df = peak_regions.merge(methylation_results, 
                                 left_index=True, 
                                 right_index=True)
    
    # Remove rows with NaN methylation values
    result_df = result_df.dropna(subset=['methylation'])
    
    logger.info("\nMethylation calculation summary:")
    logger.info(f"Total peaks with valid methylation: {len(result_df)}")
    logger.info("\nMethylation statistics by binding type:")
    for bt in result_df['binding_type'].unique():
        bt_data = result_df[result_df['binding_type'] == bt]
        logger.info(f"\n{bt}:")
        logger.info(f"Mean methylation: {bt_data['methylation'].mean():.2f}%")
        logger.info(f"Median methylation: {bt_data['methylation'].median():.2f}%")
    
    return result_df