#%% Import required libraries
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import os
import pyBigWig
import pysam
import logging
from typing import Dict, List, Tuple, Any, Union
import pyranges as pr
from functools import partial
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from tqdm import tqdm
import time
import argparse
from config import CONFIG, PATHS, logger
from cache_utils import save_to_cache, load_from_cache, clear_cache
import multiprocessing


# Set default plotting style
sns.set_theme(style="whitegrid")  # Use seaborn's built-in style
plt.rcParams['figure.figsize'] = [10, 6]
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300

# set working directory
os.chdir("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Methylation_dev")

from typing import Tuple

def analyze_methylation_patterns_parallel(results: Dict[str, pd.DataFrame], 
                                        output_dir: str,
                                        n_processes: int = None):
    """Parallel version of methylation pattern analysis"""
    if n_processes is None:
        n_processes = max(1, multiprocessing.cpu_count() - 1)

    # Create tasks for parallel processing
    tasks = []
    for cell_type, df in results.items():
        tasks.extend([
            (df, cell_type, output_dir, 'distribution'),
            (df, cell_type, output_dir, 'expression'),
            (df, cell_type, output_dir, 'binding')
        ])

    # Process tasks in parallel
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = []
        for task in tasks:
            if task[3] == 'distribution':
                futures.append(executor.submit(analyze_methylation_distribution, *task[:-1]))
            elif task[3] == 'expression':
                futures.append(executor.submit(analyze_methylation_expression_relationship, *task[:-1]))
            else:
                futures.append(executor.submit(analyze_binding_patterns, *task[:-1]))

        # Wait for all tasks to complete
        for future in tqdm(futures, desc="Analyzing methylation patterns"):
            try:
                future.result()
            except Exception as e:
                logger.error(f"Error in parallel analysis: {str(e)}")

def calculate_coverage_uniformity(df: pd.DataFrame) -> pd.Series:
    """Calculate coverage uniformity score"""
    try:
        # For now, return a simple placeholder score
        # This should be improved with actual coverage calculations
        return pd.Series(0.8, index=df.index)
    except Exception as e:
        logger.warning(f"Error calculating coverage uniformity: {str(e)}")
        return pd.Series(1.0, index=df.index)

def calculate_signal_to_noise(df: pd.DataFrame) -> pd.Series:
    """Calculate signal-to-noise ratio"""
    try:
        # For now, return a simple placeholder score
        # This should be improved with actual signal calculations
        return pd.Series(2.5, index=df.index)
    except Exception as e:
        logger.warning(f"Error calculating signal-to-noise ratio: {str(e)}")
        return pd.Series(2.0, index=df.index)
    
#Debug data merging
def debug_merge_issues():
    """Debug the data merging issues"""
    # Load data
    genes_df, gene_name_to_id = load_gene_annotations(PATHS['gtf_file'])  # Get both return values
    expression_data = load_expression_data(PATHS['rnaseq'], gene_name_to_id)  # Pass the mapping
    mecp2_binding = pd.read_csv(os.path.join(PATHS['mecp2_dir'], PATHS['mecp2_file']))
    
    # Print sample data from each DataFrame
    print("\nGene annotations sample:")
    print("Shape:", genes_df.shape)
    print("Sample gene_ids:", genes_df['gene_id'].head().tolist())
    print("\nColumns:", genes_df.columns.tolist())
    
    # Print sample of gene name mapping
    print("\nSample of gene name mapping:")
    sample_mapping = dict(list(gene_name_to_id.items())[:5])
    print(sample_mapping)
    
    for cell_type, expr_df in expression_data.items():
        print(f"\n{cell_type} expression data:")
        print("Shape:", expr_df.shape)
        print("Sample gene_ids:", expr_df['gene_id'].head().tolist())
        print("Columns:", expr_df.columns.tolist())
    
    print("\nMeCP2 binding data:")
    print("Shape:", mecp2_binding.shape)
    print("Columns:", mecp2_binding.columns.tolist())
    
    # Check for common gene_ids
    for cell_type, expr_df in expression_data.items():
        common_genes = set(genes_df['gene_id']) & set(expr_df['gene_id'])
        print(f"\nNumber of common genes between annotations and {cell_type}:", len(common_genes))
        if len(common_genes) > 0:
            print("Sample common genes:", list(common_genes)[:5])
    
    # Test merge with a single cell type
    cell_type = 'NEU'
    expr_df = expression_data[cell_type]
    
    # Try merge and print intermediate results
    merged_df = genes_df.merge(expr_df, on='gene_id', how='inner')
    print(f"\nMerged DataFrame for {cell_type}:")
    print("Shape:", merged_df.shape)
    if merged_df.shape[0] > 0:
        print("Sample rows:")
        print(merged_df[['gene_id', 'gene_name', 'log2FoldChange', 'padj']].head())
    
    return genes_df, expression_data, mecp2_binding

def identify_cpg_associated_regions(genes_df: pd.DataFrame, 
                                  cpg_islands: pr.PyRanges) -> pd.DataFrame:
    """Identify the specific CpG island regions within promoters and gene bodies"""
    
    # Create PyRanges objects for promoters and gene bodies
    promoters = pr.PyRanges(
        pd.DataFrame({
            'Chromosome': genes_df['chr'],
            'Start': genes_df['promoter_start'],
            'End': genes_df['promoter_end'],
            'gene_id': genes_df['gene_id'],
            'region_type': 'promoter'
        })
    )
    
    gene_bodies = pr.PyRanges(
        pd.DataFrame({
            'Chromosome': genes_df['chr'],
            'Start': genes_df['gene_body_start'],
            'End': genes_df['gene_body_end'],
            'gene_id': genes_df['gene_id'],
            'region_type': 'gene_body'
        })
    )
    
    # Find overlaps using intersection method
    promoter_cpg = promoters.join(cpg_islands)  # Changed from intersect to join
    gene_body_cpg = gene_bodies.join(cpg_islands)  # Changed from intersect to join
    
    cpg_regions = pd.DataFrame()
    
    if len(promoter_cpg) > 0:
        promoter_regions = promoter_cpg.df
        promoter_regions = pd.DataFrame({
            'gene_id': promoter_regions['gene_id'],
            'chr': promoter_regions['Chromosome'],
            'start': promoter_regions['Start'],
            'end': promoter_regions['End'],
            'region_type': 'promoter',
            'cpg_count': promoter_regions['cpg_count'].astype(int)  # Convert to int
        })
        cpg_regions = pd.concat([cpg_regions, promoter_regions])
    
    if len(gene_body_cpg) > 0:
        gene_body_regions = gene_body_cpg.df
        gene_body_regions = pd.DataFrame({
            'gene_id': gene_body_regions['gene_id'],
            'chr': gene_body_regions['Chromosome'],
            'start': gene_body_regions['Start'],
            'end': gene_body_regions['End'],
            'region_type': 'gene_body',
            'cpg_count': gene_body_regions['cpg_count'].astype(int)  # Convert to int
        })
        cpg_regions = pd.concat([cpg_regions, gene_body_regions])
    
    if cpg_regions.empty:
        logger.warning("No CpG islands found in promoter or gene body regions")
        return pd.DataFrame(columns=['gene_id', 'chr', 'start', 'end', 
                                   'region_type', 'cpg_count'])
    
    # Reset index and sort
    cpg_regions = cpg_regions.reset_index(drop=True)
    cpg_regions = cpg_regions.sort_values(['chr', 'start'])
    
    # Log summary
    logger.info(f"Found {len(cpg_regions)} CpG island regions "
                f"({len(cpg_regions[cpg_regions['region_type'] == 'promoter'])} in promoters, "
                f"{len(cpg_regions[cpg_regions['region_type'] == 'gene_body'])} in gene bodies)")
    logger.info(f"Average CpG count: {cpg_regions['cpg_count'].mean():.1f}")
    
    return cpg_regions

def calculate_contextual_methylation_cpg_only(merged_df: pd.DataFrame,
                                            cpg_regions: pd.DataFrame,
                                            medip_dir: str,
                                            cell_type: str,
                                            genome_fasta: str,
                                            n_processes: int = None) -> pd.DataFrame:
    """Calculate methylation levels only for CpG island regions with CpG density weighting
    
    Args:
        merged_df: DataFrame with gene information and binding data
        cpg_regions: DataFrame with CpG island regions and their CpG counts
        medip_dir: Directory containing MeDIP bigWig files
        cell_type: Cell type code ('NSC' or 'NEU')
        genome_fasta: Path to genome FASTA file
        n_processes: Number of parallel processes
    
    Returns:
        DataFrame with methylation levels weighted by CpG density
    """
    
    logger.info(f"Calculating CpG-weighted methylation for {cell_type}")
    logger.info(f"Number of input regions: {len(cpg_regions)}")
    
    # Calculate methylation for CpG regions
    methylation_results = calculate_methylation_levels_parallel(
        cpg_regions[['chr', 'start', 'end']],
        medip_dir,
        cell_type,
        genome_fasta,
        n_processes
    )
    
    # Combine results by region type
    results_by_gene = {}
    for gene_id in cpg_regions['gene_id'].unique():
        gene_regions = cpg_regions[cpg_regions['gene_id'] == gene_id]
        gene_results = {
            'gene_id': gene_id,
            'promoter_methylation': None,
            'promoter_cpg_count': 0,
            'promoter_cpg_density': 0,
            'gene_body_methylation': None,
            'gene_body_cpg_count': 0,
            'gene_body_cpg_density': 0
        }
        
        for region_type in ['promoter', 'gene_body']:
            type_regions = gene_regions[gene_regions['region_type'] == region_type]
            if not type_regions.empty:
                weighted_methylation = 0
                total_cpgs = 0
                total_length = 0
                
                for idx, region in type_regions.iterrows():
                    region_methylation = methylation_results.iloc[idx]['methylation']
                    region_length = region['end'] - region['start']
                    region_cpgs = region['cpg_count']
                    
                    if region_methylation is not None and region_cpgs > 0:
                        # Weight methylation by number of CpGs in the region
                        weighted_methylation += region_methylation * region_cpgs
                        total_cpgs += region_cpgs
                        total_length += region_length
                
                if total_cpgs > 0:
                    # Calculate CpG-density-weighted methylation
                    gene_results[f'{region_type}_methylation'] = weighted_methylation / total_cpgs
                    gene_results[f'{region_type}_cpg_count'] = total_cpgs
                    # Calculate CpG density per kb
                    gene_results[f'{region_type}_cpg_density'] = (total_cpgs * 1000) / total_length
                    
                    # Log extreme values for debugging
                    if gene_results[f'{region_type}_methylation'] > 90:
                        logger.debug(f"High methylation detected for {gene_id} {region_type}: "
                                   f"{gene_results[f'{region_type}_methylation']:.1f}%")
        
        results_by_gene[gene_id] = gene_results
    
    # Create final DataFrame
    result_df = pd.DataFrame(results_by_gene.values())
    
    # Add summary statistics
    logger.info("Methylation summary statistics:")
    for region_type in ['promoter', 'gene_body']:
        meth_col = f'{region_type}_methylation'
        density_col = f'{region_type}_cpg_density'
        valid_data = result_df[result_df[meth_col].notna()]
        
        if not valid_data.empty:
            logger.info(f"\n{region_type.title()} statistics:")
            logger.info(f"Mean methylation: {valid_data[meth_col].mean():.1f}%")
            logger.info(f"Median methylation: {valid_data[meth_col].median():.1f}%")
            logger.info(f"Mean CpG density: {valid_data[density_col].mean():.1f} CpGs/kb")
            logger.info(f"Number of regions: {len(valid_data)}")
    
    # Merge with original data
    final_df = merged_df[['gene_id', 'gene_name', 'binding_type', 'mecp2_bound', 
                         'log2FoldChange', 'padj']].merge(result_df, on='gene_id', how='right')
    
    # Add quality control flags
    final_df['low_cpg_promoter'] = final_df['promoter_cpg_density'] < 10  # Flag low CpG density
    final_df['high_meth_promoter'] = final_df['promoter_methylation'] > 90  # Flag suspiciously high methylation
    
    logger.info(f"Final output shape: {final_df.shape}")
    logger.info(f"Genes with low CpG density in promoter: {final_df['low_cpg_promoter'].sum()}")
    logger.info(f"Genes with high methylation in promoter: {final_df['high_meth_promoter'].sum()}")
    
    return final_df

def identify_peak_regions(merged_df: pd.DataFrame, 
                         mecp2_binding: pd.DataFrame) -> pd.DataFrame:
    """Identify MeCP2 binding peaks associated with each gene"""
    
    # Debug input data
    logger.info("\nIdentifying peak regions:")
    logger.info(f"Merged data shape: {merged_df.shape}")
    logger.info(f"MeCP2 binding data shape: {mecp2_binding.shape}")
    
    # Create PyRanges objects for genes and peaks
    genes_df = pd.DataFrame({
        'Chromosome': merged_df['chr'],
        'Start': merged_df['promoter_start'],
        'End': merged_df['gene_body_end'],
        'gene_id': merged_df['gene_id']
    })
    
    peaks_df = pd.DataFrame({
        'Chromosome': mecp2_binding['chr'],
        'Start': mecp2_binding['start'],
        'End': mecp2_binding['end'],
        'peak_id': range(len(mecp2_binding)),
        'binding_type': mecp2_binding['binding_type']
    })
    
    logger.info("\nGenes DataFrame:")
    logger.info(f"Shape: {genes_df.shape}")
    logger.info(f"Columns: {genes_df.columns.tolist()}")
    logger.info("First few rows:")
    logger.info(genes_df.head())
    
    logger.info("\nPeaks DataFrame:")
    logger.info(f"Shape: {peaks_df.shape}")
    logger.info(f"Columns: {peaks_df.columns.tolist()}")
    logger.info("First few rows:")
    logger.info(peaks_df.head())
    
    # Create PyRanges objects
    genes_pr = pr.PyRanges(genes_df)
    peaks_pr = pr.PyRanges(peaks_df)
    
    # Find overlaps between genes and peaks
    overlaps = peaks_pr.join(genes_pr)
    
    if overlaps is None or len(overlaps) == 0:
        logger.warning("No overlaps found between peaks and genes")
        logger.warning("Check chromosome formats and coordinates")
        return pd.DataFrame()
        
    # Convert to DataFrame and add peak information
    peak_regions = overlaps.as_df()
    
    logger.info("\nOverlaps found:")
    logger.info(f"Shape: {peak_regions.shape}")
    logger.info(f"Columns: {peak_regions.columns.tolist()}")
    logger.info("First few rows:")
    logger.info(peak_regions.head())
    
    # Merge with original peak information
    peak_regions = peak_regions.merge(
        mecp2_binding,
        left_on=['Chromosome', 'Start', 'End'],
        right_on=['chr', 'start', 'end']
    )
    
    logger.info(f"\nFinal peak regions:")
    logger.info(f"Shape: {peak_regions.shape}")
    logger.info(f"Columns: {peak_regions.columns.tolist()}")
    logger.info("First few rows:")
    logger.info(peak_regions.head())
    
    return peak_regions

def calculate_peak_methylation(peak_regions: pd.DataFrame,
                             medip_dir: str,
                             cell_type: str,
                             genome_fasta: str,
                             n_processes: int = None) -> pd.DataFrame:
    """Calculate methylation levels in MeCP2 binding peaks"""
    
    # Debug input
    logger.info(f"\nCalculating methylation for {len(peak_regions)} peaks")
    logger.info(f"Input peak regions columns: {peak_regions.columns.tolist()}")
    
    # Check input files
    prefix_map = {'NEU': 'N', 'NSC': 'PP'}
    prefix = prefix_map.get(cell_type)
    if not prefix:
        raise ValueError(f"Unknown cell type: {cell_type}")
    
    ip_files = [os.path.join(medip_dir, f"Medip_{prefix}_output_r{i}.bw") 
                for i in range(1, 4)]
    input_files = [os.path.join(medip_dir, f"Medip_{prefix}_input_r{i}.bw") 
                  for i in range(1, 4)]
    
    logger.info("\nChecking input files:")
    for f in ip_files + input_files:
        exists = os.path.exists(f)
        logger.info(f"{f}: {'exists' if exists else 'MISSING'}")
    
    # Check genome file
    logger.info(f"\nGenome file: {genome_fasta}")
    if not os.path.exists(genome_fasta):
        raise FileNotFoundError(f"Genome file not found: {genome_fasta}")
    
    if n_processes is None:
        n_processes = os.cpu_count() - 1
    
    # Split peaks into chunks for parallel processing
    chunk_size = max(1, len(peak_regions) // (n_processes * 4))
    chunks = np.array_split(peak_regions, len(peak_regions) // chunk_size + 1)
    
    # Process chunks in parallel
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = [
            executor.submit(calculate_methylation_for_peaks,
                          chunk, medip_dir, cell_type, genome_fasta)
            for chunk in chunks
        ]
        
        results = []
        for future in tqdm(futures, desc="Calculating peak methylation"):
            try:
                result = future.result()
                if not result.empty:  # Only append non-empty results
                    results.append(result)
                    logger.debug(f"Processed chunk with {len(result)} peaks")
            except Exception as e:
                logger.error(f"Error processing chunk: {str(e)}")
                continue
    
    if not results:
        logger.warning("No valid methylation results obtained")
        # Return empty DataFrame with expected columns
        return pd.DataFrame(columns=['peak_id', 'gene_id', 'chr', 'start', 'end', 
                                   'binding_type', 'methylation', 'cpg_count', 'cpg_density'])
    
    final_df = pd.concat(results, ignore_index=True)
    logger.info(f"\nFinal methylation results:")
    logger.info(f"Shape: {final_df.shape}")
    logger.info(f"Columns: {final_df.columns.tolist()}")
    logger.info("First few rows:")
    logger.info(final_df.head())
    
    return final_df

def calculate_methylation_for_peaks(peaks_df: pd.DataFrame,
                                  medip_dir: str,
                                  cell_type: str,
                                  genome_fasta: str) -> pd.DataFrame:
    """Calculate methylation for a set of peaks with accurate CpG methylation estimation"""
    
    prefix_map = {'NEU': 'N', 'NSC': 'PP'}
    prefix = prefix_map.get(cell_type)
    
    ip_files = [os.path.join(medip_dir, f"Medip_{prefix}_output_r{i}.bw") 
                for i in range(1, 4)]
    input_files = [os.path.join(medip_dir, f"Medip_{prefix}_input_r{i}.bw") 
                  for i in range(1, 4)]
    
    results = []
    with pysam.FastaFile(genome_fasta) as fasta:
        for _, peak in peaks_df.iterrows():
            try:
                # Debug info
                logger.debug(f"Processing peak: {peak['Chromosome']}:{peak['Start']}-{peak['End']}")
                
                # Get sequence and count CpGs
                sequence = fasta.fetch(str(peak['Chromosome']), 
                                     int(peak['Start']), 
                                     int(peak['End']))
                cpg_count = sequence.upper().count('CG')
                
                if cpg_count == 0:
                    logger.debug(f"No CpGs in peak {peak['peak_id']}")
                    continue
                
                # Calculate methylation for each replicate
                replicate_values = []
                for ip_file, input_file in zip(ip_files, input_files):
                    with pyBigWig.open(ip_file) as bw_ip, pyBigWig.open(input_file) as bw_input:
                        try:
                            # Get signal values
                            ip_values = bw_ip.values(str(peak['Chromosome']), 
                                                   int(peak['Start']), 
                                                   int(peak['End']))
                            input_values = bw_input.values(str(peak['Chromosome']), 
                                                         int(peak['Start']), 
                                                         int(peak['End']))
                            
                            if ip_values and input_values:
                                # Calculate enrichment
                                ip_mean = np.nanmean([x for x in ip_values if x is not None])
                                input_mean = np.nanmean([x for x in input_values if x is not None])
                                
                                if input_mean > 0:  # Avoid division by zero
                                    enrichment = ip_mean / input_mean
                                    # Convert enrichment to methylation percentage
                                    methylation = min(100, max(0, enrichment * 20))  # Scale factor of 20
                                    replicate_values.append(methylation)
                                    
                        except Exception as e:
                            logger.debug(f"Error processing replicate for peak {peak['peak_id']}: {str(e)}")
                            continue
                
                # Require at least 2 valid replicates
                if len(replicate_values) >= 2:
                    methylation = np.mean(replicate_values)
                    results.append({
                        'peak_id': peak['peak_id'],
                        'gene_id': peak['gene_id'],
                        'chr': peak['Chromosome'],
                        'start': peak['Start'],
                        'end': peak['End'],
                        'binding_type': peak['binding_type_x'],
                        'methylation': methylation,
                        'cpg_count': cpg_count,
                        'cpg_density': (cpg_count * 1000) / (peak['End'] - peak['Start'])
                    })
                    logger.debug(f"Successfully processed peak {peak['peak_id']}: "
                               f"methylation={methylation:.1f}%, CpGs={cpg_count}")
                
            except Exception as e:
                logger.error(f"Error processing peak {peak['peak_id']}: {str(e)}")
                continue
    
    if not results:
        logger.warning("No valid methylation results obtained")
        return pd.DataFrame(columns=['peak_id', 'gene_id', 'chr', 'start', 'end', 
                                   'binding_type', 'methylation', 'cpg_count', 'cpg_density'])
    
    result_df = pd.DataFrame(results)
    logger.info(f"Successfully calculated methylation for {len(result_df)} peaks")
    logger.info(f"Mean methylation: {result_df['methylation'].mean():.1f}%")
    logger.info(f"Mean CpG density: {result_df['cpg_density'].mean():.1f} CpGs/kb")
    
    return result_df
