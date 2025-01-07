import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import os
import pyBigWig
import pysam
import logging
from typing import Dict, List, Tuple, Any, Union, Optional
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


def analyze_tss_binding_patterns(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze binding patterns specifically around TSS and their relationship with gene regulation"""
    
    # Check if we have methylation data
    has_methylation = all(col in df.columns for col in ['cpg_island_count', 'cpg_island_mean_methylation'])
    logger.info(f"Methylation data available: {has_methylation}")
    
    # Create output directories
    tss_dir = os.path.join(output_dir, 'tss_analysis', cell_type)
    os.makedirs(tss_dir, exist_ok=True)
    
    # Define binding categories
    binding_categories = {
        'exo_enriched': df[df['binding_type'].isin(['exo', 'both'])],
        'exo_only': df[df['binding_type'] == 'exo'],
        'endo_only': df[df['binding_type'] == 'endo'],
        'non_enriched': df[~df['mecp2_bound']]
    }
    
    # Analyze each category
    results = {}
    for category, category_df in binding_categories.items():
        # Skip empty categories
        if len(category_df) == 0:
            logger.warning(f"No genes found in group: {category} for {cell_type}")
            continue
            
        # Split by regulation status
        regulation_groups = {
            'not_deregulated': category_df[category_df['expression_status'] == 'unchanged'],
            'upregulated': category_df[category_df['expression_status'] == 'upregulated'],
            'downregulated': category_df[category_df['expression_status'] == 'downregulated']
        }
        
        # Calculate statistics for each group
        group_stats = {}
        for reg_status, group_df in regulation_groups.items():
            if len(group_df) > 0:
                stats = {
                    'count': len(group_df)
                }
                
                # Add methylation stats only if available
                if has_methylation:
                    stats.update({
                        'methylation': {
                            'mean': group_df['cpg_island_mean_methylation'].mean(),
                            'std': group_df['cpg_island_mean_methylation'].std(),
                            'median': group_df['cpg_island_mean_methylation'].median()
                        },
                        'cpg_count': {
                            'mean': group_df['cpg_island_count'].mean(),
                            'std': group_df['cpg_island_count'].std(),
                            'median': group_df['cpg_island_count'].median()
                        }
                    })
                
                group_stats[reg_status] = stats
                
                # Save gene lists
                output_file = os.path.join(tss_dir, f'{category}_{reg_status}_genes.csv')
                group_df.to_csv(output_file, index=False)
        
        results[category] = {
            'regulation_counts': pd.Series({k: len(v) for k, v in regulation_groups.items()}),
            'group_stats': group_stats
        }
    
    return results

def run_analysis_pipeline(force_recompute: bool = False, n_processes: int = None,
                        cell_type: str = None, chromosome: str = None) -> Dict[str, Any]:
    """Run the analysis pipeline with improved error handling"""
    try:
        if force_recompute:
            clear_cache()
            logger.info("Cache cleared due to force recompute flag")
        
        # Verify input files first
        if not verify_input_files():
            raise FileNotFoundError("Required input files are missing")
        
        # Create base output directory structure
        base_output_dir = os.path.join(PATHS['output_dir'], f'analyze_mecp2_cpg_enrichment_{CONFIG["experiment"]}')
        try:
            os.makedirs(base_output_dir, exist_ok=True)
            logger.info(f"Created base output directory: {base_output_dir}")
            
            # Store the base output directory in PATHS
            PATHS['base_output_dir'] = base_output_dir
            
            # Create required subdirectories
            subdirs = ['cpg_analysis', 'detailed_analysis']
            for subdir in subdirs:
                subdir_path = os.path.join(base_output_dir, subdir)
                os.makedirs(subdir_path, exist_ok=True)
                logger.info(f"Created subdirectory: {subdir_path}")
            
            # Create cell-type specific directories if cell_type is specified
            if cell_type:
                cell_type_dir = os.path.join(base_output_dir, 'cpg_analysis', cell_type)
                os.makedirs(cell_type_dir, exist_ok=True)
                logger.info(f"Created cell-type directory: {cell_type_dir}")
            
        except Exception as e:
            logger.error(f"Error creating directory structure: {str(e)}")
            raise
        
        # Load data with validation
        try:
            if not verify_input_files():
                raise FileNotFoundError("Required input files are missing")
                
            mecp2_path = os.path.join(PATHS['mecp2_dir'], PATHS['mecp2_file'])
            mecp2_binding = validate_mecp2_data(mecp2_path)
            logger.info(f"Successfully loaded MeCP2 binding data: {len(mecp2_binding)} rows")
            logger.debug(f"MeCP2 columns: {mecp2_binding.columns.tolist()}")
        except Exception as e:
            logger.error(f"Failed to load MeCP2 data: {str(e)}")
            raise
        
        # Load all data first
        genes_df, gene_name_to_id = load_gene_annotations(PATHS['gtf_file'])
        
        # Filter by chromosome if specified
        if chromosome:
            mecp2_binding = mecp2_binding[mecp2_binding['chr'] == chromosome]
            genes_df = genes_df[genes_df['chr'] == chromosome]
            logger.info(f"Filtered data for chromosome {chromosome}")
        
        # Load only specified cell type if provided
        if cell_type:
            expression_data = {
                cell_type: load_expression_data(PATHS['rnaseq'][cell_type], gene_name_to_id, single_file=True)
            }
            logger.info(f"Analyzing cell type: {cell_type}")
        else:
            expression_data = load_expression_data(PATHS['rnaseq'], gene_name_to_id)
        
        # Run analysis
        results = {}
        tss_results = {}
        
        for ct, expr_df in expression_data.items():
            logger.info(f"\nProcessing {ct}...")
            logger.info(f"Expression data shape: {expr_df.shape}")
            logger.debug(f"Expression columns: {expr_df.columns.tolist()}")
            
            # Use base_output_dir for output paths
            output_dir = os.path.join(base_output_dir, 'cpg_analysis', ct)
            os.makedirs(output_dir, exist_ok=True)
            logger.info(f"Using output directory: {output_dir}")
            
            # Handle SettingWithCopyWarning
            pd.options.mode.chained_assignment = None
            
            # Create detailed analysis directory
            detailed_dir = os.path.join(base_output_dir, 'detailed_analysis')
            os.makedirs(detailed_dir, exist_ok=True)
            logger.info(f"Using detailed analysis directory: {detailed_dir}")
            
            # Merge data
            merged_df = validate_and_merge_data(genes_df, expr_df, mecp2_binding)
            
            # Ensure column names match PyRanges requirements
            if 'chr' in merged_df.columns:
                merged_df = merged_df.rename(columns={
                    'chr': 'Chromosome',
                    'start': 'Start',
                    'end': 'End'
                })
            
            logger.info(f"Merged data columns: {merged_df.columns.tolist()}")
            
            # Calculate methylation levels
            methylation_data = calculate_methylation_levels(
                merged_df,
                PATHS['medip_dir'],
                ct,
                PATHS['genome_fasta']
            )
            
            # Merge methylation data with other data if available
            if not methylation_data.empty and len(methylation_data.columns) > 0:
                results[ct] = pd.concat([merged_df, methylation_data], axis=1)
                logger.info("Methylation data added successfully")
            else:
                logger.warning("No methylation data available - proceeding with basic analysis")
                results[ct] = merged_df
            
            # Save intermediate results
            output_dir = os.path.join(PATHS['output_dir'], 'cpg_analysis', ct)
            os.makedirs(output_dir, exist_ok=True)
            
            logger.info(f"Saving results to {output_dir}")
            results[ct].to_csv(os.path.join(output_dir, f'{ct}_full_results.csv'), index=False)
            
            # Log data shape and columns for debugging
            logger.info(f"Results shape for {ct}: {results[ct].shape}")
            logger.info(f"Available columns: {results[ct].columns.tolist()}")
            
            # Perform TSS binding analysis with correct output path
            logger.info(f"Performing TSS binding analysis for {ct}...")
            tss_results[ct] = analyze_tss_binding_patterns(results[ct], ct, base_output_dir)
            
            # Add detailed methylation analysis with correct output path
            methylation_stats = analyze_methylation_patterns(results[ct], ct, base_output_dir)
            
            # Add new detailed binding-methylation analysis
            binding_results = analyze_binding_methylation_expression(results[ct], ct, base_output_dir)
            
            # Run other analyses with correct output paths
            exo_enriched_results = analyze_exo_enriched_methylation(results[ct], ct, base_output_dir)
            binding_patterns = analyze_mecp2_detailed_binding(results[ct], ct, base_output_dir)
            occupancy_patterns = analyze_cpg_occupancy(results[ct], ct, base_output_dir)
            regulatory_features = analyze_regulatory_features(results[ct], ct, base_output_dir)
        
        return {
            'methylation_results': results,
            'tss_results': tss_results,
            'genes_df': genes_df,
            'expression_data': expression_data,
            'mecp2_binding': mecp2_binding
        }
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        logger.error("Stack trace:", exc_info=True)
        raise
    finally:
        # Reset pandas warning settings
        pd.options.mode.chained_assignment = 'warn'

def verify_input_files():
    """Verify that all required input files exist with detailed logging"""
    all_files_found = True
    
    # Check mecp2 file
    mecp2_path = os.path.join(PATHS['mecp2_dir'], PATHS['mecp2_file'])
    logger.info(f"Checking for MeCP2 file at: {mecp2_path}")
    
    if os.path.exists(mecp2_path):
        logger.info(f"Found MeCP2 file: {mecp2_path}")
        return True  # If we found the file, return immediately
    else:
        logger.warning(f"MeCP2 file not found at: {mecp2_path}")
        # Try to find the file in subdirectories
        for root, dirs, files in os.walk(PATHS['mecp2_dir']):
            for file in files:
                if file == PATHS['mecp2_file']:
                    found_path = os.path.join(root, file)
                    logger.info(f"Found MeCP2 file at: {found_path}")
                    # Update the paths to use the correct location
                    PATHS['mecp2_dir'] = root
                    logger.info(f"Updated mecp2_dir to: {root}")
                    return True
    
    # If we get here, we didn't find the file
    logger.error(f"Could not find MeCP2 file {PATHS['mecp2_file']} in any subdirectory of {PATHS['mecp2_dir']}")
    logger.error("Directory structure:")
    os.system(f"tree {PATHS['mecp2_dir']}")
    return False

def analyze_cpg_binding_correlation(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze correlation between CpG methylation and MeCP2 binding"""
    
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df, x='binding_type', y='cpg_island_mean_methylation')
    plt.title(f'{cell_type}: CpG Methylation by MeCP2 Binding Type')
    plt.xlabel('Binding Type')
    plt.ylabel('CpG Island Methylation (%)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cpg_binding_correlation.pdf'))
    plt.close()
    
    # Save statistics
    stats_file = os.path.join(output_dir, 'cpg_binding_stats.txt')
    with open(stats_file, 'w') as f:
        f.write(f"CpG Methylation by Binding Type - {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        for binding_type in df['binding_type'].unique():
            subset = df[df['binding_type'] == binding_type]
            f.write(f"\n{binding_type}:\n")
            f.write(f"Count: {len(subset)}\n")
            f.write(f"Mean methylation: {subset['cpg_island_mean_methylation'].mean():.2f}%\n")
            f.write(f"Median methylation: {subset['cpg_island_mean_methylation'].median():.2f}%\n")
            f.write(f"Std Dev: {subset['cpg_island_mean_methylation'].std():.2f}%\n")

def analyze_cpg_expression_correlation(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze correlation between CpG methylation and gene expression"""
    
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df, x='cpg_island_mean_methylation', y='log2FoldChange', 
                    hue='binding_type', alpha=0.6)
    plt.title(f'{cell_type}: CpG Methylation vs Expression')
    plt.xlabel('CpG Island Methylation (%)')
    plt.ylabel('log2 Fold Change')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cpg_expression_correlation.pdf'))
    plt.close()

#Helper functions for genomic regions
def get_promoter_region(gene_start: int, gene_end: int, strand: str) -> Tuple[int, int]:
    """Get promoter coordinates based on gene location and strand"""
    if strand == '+':
        return (gene_start + CONFIG['genomic_regions']['promoter'][0], 
                gene_start + CONFIG['genomic_regions']['promoter'][1])
    else:
        return (gene_end - CONFIG['genomic_regions']['promoter'][1], 
                gene_end - CONFIG['genomic_regions']['promoter'][0])

def get_gene_body_region(gene_start: int, gene_end: int, strand: str) -> Tuple[int, int]:
    """Get gene body coordinates based on gene location and strand"""
    if strand == '+':
        return (gene_start + CONFIG['genomic_regions']['gene_body_start'], gene_end)
    else:
        return (gene_start, gene_end - CONFIG['genomic_regions']['gene_body_start']) 

#Analyze methylation patterns
def analyze_methylation_patterns(genes_df: pd.DataFrame, 
                               expression_data: Dict[str, pd.DataFrame],
                               mecp2_binding: pd.DataFrame,
                               medip_dir: str,
                               genome_fasta: str,
                               use_cache: bool = True) -> Dict[str, pd.DataFrame]:
    """Improved methylation pattern analysis"""
    results = {}
    
    for cell_type, expr_df in expression_data.items():
        logger.info(f"Analyzing {cell_type} methylation patterns...")
        
        # Try to load from cache if enabled
        cache_key = f"methylation_patterns_{cell_type}"
        if use_cache and CONFIG['cache']['enabled']:
            cached_result = load_from_cache(cache_key)
            if cached_result is not None:
                results[cell_type] = cached_result
                continue
        
        # Merge data with proper validation
        merged_df = validate_and_merge_data(genes_df, expr_df, mecp2_binding)
        
        # Calculate methylation with biological context
        methylation_results = calculate_contextual_methylation(
            merged_df,
            medip_dir,
            cell_type,
            genome_fasta
        )
        
        # Add genomic context analysis
        methylation_results = add_genomic_context(methylation_results)
        
        # Validate results against known biological patterns
        methylation_results = validate_methylation_patterns(methylation_results)
        
        results[cell_type] = methylation_results
        
        # Cache results if enabled
        if CONFIG['cache']['enabled']:
            save_to_cache(methylation_results, cache_key)
    
    return results


#Modified print_regulated_summary function
def print_regulated_summary_v2(results: Dict[str, pd.DataFrame]):
    """Print summary statistics for regulated MeCP2-bound genes with improved analysis"""
    for cell_type, df in results.items():
        print(f"\nSummary for {cell_type}:")
        print("="*50)
        
        # First, check regulated genes
        regulated_df = df[df['expression_status'].isin(['upregulated', 'downregulated'])]
        print(f"\nTotal regulated genes: {len(regulated_df)}")
        
        if len(regulated_df) == 0:
            print("No regulated genes found.")
            continue
        
        # Check MeCP2 binding
        mecp2_bound_df = regulated_df[regulated_df['mecp2_bound']]
        print(f"MeCP2-bound regulated genes: {len(mecp2_bound_df)}")
        
        if len(mecp2_bound_df) == 0:
            print("No MeCP2-bound regulated genes found.")
            continue
        
        # Check binding types present
        binding_types = mecp2_bound_df['binding_type'].value_counts()
        print("\nBinding types present:")
        print(binding_types)
        
        # Analyze each binding type
        for binding_type in binding_types.index:
            if pd.isna(binding_type):
                continue
                
            subset_df = mecp2_bound_df[mecp2_bound_df['binding_type'] == binding_type]
            
            print(f"\n{binding_type.upper()} BINDING:")
            print("-"*20)
            
            # Gene counts and basic statistics
            expr_counts = subset_df['expression_status'].value_counts()
            print("\nGene counts by expression status:")
            print(expr_counts)
            
            # Calculate and display statistics
            metrics = {
                'promoter_methylation': 'Promoter Methylation',
                'gene_body_methylation': 'Gene Body Methylation',
                'promoter_cpg_count': 'Promoter CpG Count'
            }
            
            for metric, metric_name in metrics.items():
                print(f"\n{metric_name}:")
                
                # Calculate statistics by expression status
                stats_df = subset_df.groupby('expression_status')[metric].agg([
                    'count', 'mean', 'std', 'min', 'max'
                ]).round(3)
                print(stats_df)
                
                # Perform statistical tests
                up_genes = subset_df[subset_df['expression_status'] == 'upregulated']
                down_genes = subset_df[subset_df['expression_status'] == 'downregulated']
                
                if len(up_genes) > 0 and len(down_genes) > 0:
                    stats_results = perform_statistical_tests(
                        up_genes[metric],
                        down_genes[metric],
                        metric_name
                    )
                    
                    if stats_results.get('mannwhitney'):
                        print(f"\nMann-Whitney U test for {metric_name}:")
                        print(f"Statistic: {stats_results['mannwhitney']['statistic']:.2f}")
                        print(f"P-value: {stats_results['mannwhitney']['pvalue']:.4f}")
                        
                        if stats_results.get('cohens_d') is not None:
                            print(f"Cohen's d: {stats_results['cohens_d']:.4f}")
                        else:
                            print("Cohen's d: Not available")
                    
                    if stats_results.get('error'):
                        print(f"Error in statistical tests: {stats_results['error']}")
            
            # Calculate correlation
            if len(subset_df) > 0:
                corr_result = calculate_correlation(
                    subset_df['promoter_methylation'],
                    subset_df['gene_body_methylation']
                )
                if corr_result:
                    print("\nCorrelation between promoter and gene body methylation:")
                    print(f"Spearman correlation: rho={corr_result['statistic']:.4f}, p={corr_result['pvalue']:.4f}")
                else:
                    print("\nCorrelation analysis: Not available (insufficient data or no variation)")

#Calculate methylation levels for multiple replicates
def calculate_methylation_levels_with_replicates(region_df: pd.DataFrame,
                                               medip_dir: str,
                                               cell_type_prefix: str,
                                               genome_fasta: str) -> pd.DataFrame:
    """Calculate methylation levels using improved normalization"""
    cell_type_map = {
        'NEU': 'N',
        'NSC': 'PP'
    }
    
    prefix = cell_type_map.get(cell_type_prefix)
    if not prefix:
        raise ValueError(f"Unknown cell type: {cell_type_prefix}")
    
    # Get replicate and input files
    replicate_files = [
        os.path.join(medip_dir, f"Medip_{prefix}_output_r{i}.bw")
        for i in range(1, 4)
    ]
    input_files = [
        os.path.join(medip_dir, f"Medip_{prefix}_input_r{i}.bw")
        for i in range(1, 4)
    ]
    
    logger.info(f"Processing replicates for {cell_type_prefix} ({prefix})")
    
    # Validate and clean region coordinates
    region_df = region_df.copy()
    region_df['start'] = region_df['start'].clip(lower=0)
    
    # Get chromosome sizes and validate coordinates
    first_bw = pyBigWig.open(replicate_files[0])
    chrom_sizes = dict(first_bw.chroms().items())
    first_bw.close()
    
    for chrom in region_df['chr'].unique():
        if chrom in chrom_sizes:
            mask = (region_df['chr'] == chrom)
            region_df.loc[mask, 'end'] = region_df.loc[mask, 'end'].clip(upper=chrom_sizes[chrom])
    
    # Calculate methylation for each replicate
    replicate_results = []
    for rep_file, input_file in zip(replicate_files, input_files):
        if not (os.path.exists(rep_file) and os.path.exists(input_file)):
            logger.error(f"File not found: {rep_file} or {input_file}")
            continue
            
        logger.info(f"Processing {os.path.basename(rep_file)}")
        try:
            bw_ip = pyBigWig.open(rep_file)
            bw_input = pyBigWig.open(input_file)
            fasta = pysam.FastaFile(genome_fasta)  # Open genome fasta file
            
            methylation_values = []
            cpg_counts = []
            
            # Process regions
            for _, row in region_df.iterrows():
                try:
                    # Get IP and input values
                    ip_values = bw_ip.values(row['chr'], int(row['start']), int(row['end']))
                    input_values = bw_input.values(row['chr'], int(row['start']), int(row['end']))
                    
                    if ip_values is None or input_values is None:
                        methylation = 0
                        cpg_count = 0
                    else:
                        # Filter out None values
                        valid_pairs = [
                            (ip, inp) for ip, inp in zip(ip_values, input_values)
                            if ip is not None and inp is not None
                        ]
                        
                        if valid_pairs:
                            ip_vals, input_vals = zip(*valid_pairs)
                            
                            # Get sequence and calculate methylation
                            sequence = fasta.fetch(row['chr'], int(row['start']), int(row['end']))
                            methylation = calculate_normalized_methylation(ip_vals, input_vals, sequence)
                            cpg_count = sequence.upper().count('CG')
                        else:
                            methylation = 0
                            cpg_count = 0
                    
                    methylation_values.append(methylation)
                    cpg_counts.append(cpg_count)
                    
                except Exception as e:
                    logger.debug(f"Error processing region {row['chr']}:{row['start']}-{row['end']}: {str(e)}")
                    methylation_values.append(0)
                    cpg_counts.append(0)
            
            bw_ip.close()
            bw_input.close()
            fasta.close()
            
            replicate_results.append(methylation_values)
            
        except Exception as e:
            logger.error(f"Error processing {rep_file}: {str(e)}")
            continue
    
    if not replicate_results:
        raise RuntimeError(f"No valid data found for {cell_type_prefix}")
    
    # Calculate average methylation across replicates
    avg_methylation = np.mean(replicate_results, axis=0)
    
    return pd.DataFrame({
        'methylation': avg_methylation,
        'cpg_count': cpg_counts
    })

def calculate_normalized_methylation(ip_values: List[float], 
                                  input_values: List[float],
                                  window_size: Optional[int] = None) -> float:
    """
    Calculate normalized methylation level from IP and input values
    
    Args:
        ip_values: List of IP signal values
        input_values: List of input signal values
        window_size: Optional window size for smoothing
    
    Returns:
        Normalized methylation percentage
    """
    try:
        # Convert to numpy arrays and remove any None values
        ip_array = np.array([x for x in ip_values if x is not None])
        input_array = np.array([x for x in input_values if x is not None])
        
        if len(ip_array) == 0 or len(input_array) == 0:
            return 0.0
            
        # Apply smoothing if window size is provided
        if window_size:
            ip_array = np.convolve(ip_array, np.ones(window_size)/window_size, mode='valid')
            input_array = np.convolve(input_array, np.ones(window_size)/window_size, mode='valid')
        
        # Calculate mean values
        ip_mean = np.mean(ip_array)
        input_mean = np.mean(input_array)
        
        # Avoid division by zero
        if input_mean == 0:
            if ip_mean == 0:
                return 0.0
    else:
                return 100.0
                
        # Calculate normalized methylation
        methylation = (ip_mean / input_mean) * 100
        
        # Clip to valid range
        methylation = np.clip(methylation, 0, 100)
        
        logger.debug(f"Methylation calculation:")
        logger.debug(f"IP mean: {ip_mean:.2f}")
        logger.debug(f"Input mean: {input_mean:.2f}")
        logger.debug(f"Methylation: {methylation:.2f}%")
        
        return float(methylation)
        
    except Exception as e:
        logger.error(f"Error calculating methylation: {str(e)}")
        return 0.0

#Statistical analysis
def perform_statistical_analysis(results: Dict[str, pd.DataFrame]) -> Dict[str, Dict]:
    """Perform statistical tests on methylation patterns"""
    stats_results = {}
    
    for cell_type, df in results.items():
        cell_stats = {}
        
        # Analyze differences between up and down regulated genes
        for region in ['promoter', 'gene_body']:
            # Methylation differences
            up_meth = df[df['expression_status'] == 'upregulated'][f'{region}_methylation']
            down_meth = df[df['expression_status'] == 'downregulated'][f'{region}_methylation']
            
            stat, pval = stats.mannwhitneyu(up_meth, down_meth)
            
            # CpG density differences
            up_cpg = df[df['expression_status'] == 'upregulated'][f'{region}_cpg_count']
            down_cpg = df[df['expression_status'] == 'downregulated'][f'{region}_cpg_count']
            
            cpg_stat, cpg_pval = stats.mannwhitneyu(up_cpg, down_cpg)
            
            cell_stats[region] = {
                'methylation_comparison': {
                    'statistic': stat,
                    'pvalue': pval
                },
                'cpg_density_comparison': {
                    'statistic': cpg_stat,
                    'pvalue': cpg_pval
                }
            }
        
        stats_results[cell_type] = cell_stats
    
    return stats_results

#Create visualization functions
def create_methylation_patterns_plot(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Create detailed methylation pattern plots"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 15))
    
    # Promoter methylation by binding type and expression
    sns.boxplot(data=df, x='binding_type', y='promoter_methylation',
                hue='expression_status', ax=ax1)
    ax1.set_title('Promoter Methylation')
    ax1.set_xlabel('Binding Type')
    ax1.set_ylabel('Methylation (%)')
    
    # Gene body methylation
    sns.boxplot(data=df, x='binding_type', y='gene_body_methylation',
                hue='expression_status', ax=ax2)
    ax2.set_title('Gene Body Methylation')
    ax2.set_xlabel('Binding Type')
    ax2.set_ylabel('Methylation (%)')
    
    # Add individual points
    sns.stripplot(data=df, x='binding_type', y='promoter_methylation',
                 hue='expression_status', dodge=True, alpha=0.3, ax=ax1)
    sns.stripplot(data=df, x='binding_type', y='gene_body_methylation',
                 hue='expression_status', dodge=True, alpha=0.3, ax=ax2)
    
    # Methylation correlation
    sns.scatterplot(data=df, x='promoter_methylation', y='gene_body_methylation',
                   hue='expression_status', style='binding_type', ax=ax3)
    ax3.set_title('Methylation Correlation')
    
    # CpG density
    sns.boxplot(data=df, x='binding_type', y='promoter_cpg_count',
                hue='expression_status', ax=ax4)
    ax4.set_title('CpG Density')
    
    plt.suptitle(f'{cell_type}: Methylation Patterns', y=1.02, size=16)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'regulatory_patterns',
                            f'{cell_type}_methylation_patterns.pdf'))
    plt.close()

def create_binding_correlation_plot(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Create correlation plots for binding strength and methylation"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 15))
    
    # Exo signal correlations
    sns.scatterplot(data=df, x='exo_signal', y='promoter_methylation',
                   hue='expression_status', ax=ax1)
    ax1.set_title('Exo Signal vs Promoter Methylation')
    
    sns.scatterplot(data=df, x='exo_signal', y='gene_body_methylation',
                   hue='expression_status', ax=ax2)
    ax2.set_title('Exo Signal vs Gene Body Methylation')
    
    # Endo signal correlations
    sns.scatterplot(data=df, x='endo_signal', y='promoter_methylation',
                   hue='expression_status', ax=ax3)
    ax3.set_title('Endo Signal vs Promoter Methylation')
    
    sns.scatterplot(data=df, x='endo_signal', y='gene_body_methylation',
                   hue='expression_status', ax=ax4)
    ax4.set_title('Endo Signal vs Gene Body Methylation')
    
    plt.suptitle(f'{cell_type}: Binding Strength Correlations', y=1.02, size=16)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'regulatory_patterns',
                            f'{cell_type}_binding_correlations.pdf'))
    plt.close()

def create_methylation_distribution_plot(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Create distribution plots for methylation and binding signals"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 15))
    
    # Function to safely create KDE plot
    def safe_kdeplot(data, x, hue, ax, title):
        try:
            if data[x].var() > 0:  # Check if there's variance in the data
                sns.kdeplot(data=data, x=x, hue=hue, ax=ax)
            else:
                logger.warning(f"No variance in {x} data for {title}")
                ax.text(0.5, 0.5, 'No variation in data',
                       horizontalalignment='center',
                       verticalalignment='center',
                       transform=ax.transAxes)
        except Exception as e:
            logger.warning(f"Error creating KDE plot for {title}: {str(e)}")
            ax.text(0.5, 0.5, 'Error creating plot',
                   horizontalalignment='center',
                   verticalalignment='center',
                   transform=ax.transAxes)
        ax.set_title(title)
    
    # Create plots with error handling
    safe_kdeplot(df, 'promoter_methylation', 'expression_status', ax1, 'Promoter Methylation Distribution')
    safe_kdeplot(df, 'gene_body_methylation', 'expression_status', ax2, 'Gene Body Methylation Distribution')
    safe_kdeplot(df, 'exo_signal', 'expression_status', ax3, 'Exo Signal Distribution')
    safe_kdeplot(df, 'endo_signal', 'expression_status', ax4, 'Endo Signal Distribution')
    
    plt.suptitle(f'{cell_type}: Feature Distributions', y=1.02, size=16)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'regulatory_patterns',
                            f'{cell_type}_distributions.pdf'))
    plt.close()

#Debug the filtering steps
def debug_regulated_genes_filtering(results: Dict[str, pd.DataFrame]):
    """Debug the filtering of regulated genes with MeCP2 binding"""
    for cell_type, df in results.items():
        print(f"\nDebugging {cell_type}:")
        print("="*50)
        
        # Check total genes
        print(f"\nTotal genes: {len(df)}")
        
        # Check expression status distribution
        print("\nExpression status distribution:")
        print(df['expression_status'].value_counts())
        
        # Check MeCP2 binding
        print("\nMeCP2 binding distribution:")
        print(df['mecp2_bound'].value_counts())
        
        # Check binding types
        print("\nBinding type distribution:")
        print(df['binding_type'].value_counts())
        
        # Check regulated genes
        regulated = df[df['expression_status'].isin(['upregulated', 'downregulated'])]
        print(f"\nRegulated genes: {len(regulated)}")
        
        # Check MeCP2-bound regulated genes
        bound_regulated = regulated[regulated['mecp2_bound']]
        print(f"\nMeCP2-bound regulated genes: {len(bound_regulated)}")
        
        if len(bound_regulated) > 0:
            print("\nBinding type distribution in regulated genes:")
            print(bound_regulated['binding_type'].value_counts())
            print("\nExpression status in MeCP2-bound regulated genes:")
            print(bound_regulated['expression_status'].value_counts())
            
            # Sample of the data
            print("\nSample of MeCP2-bound regulated genes:")
            sample_cols = ['gene_name', 'expression_status', 'binding_type', 
                         'promoter_methylation', 'gene_body_methylation']
            print(bound_regulated[sample_cols].head())

#Modified statistical test functions
def perform_statistical_tests(up_data: pd.Series, down_data: pd.Series, metric_name: str) -> Dict:
    """Perform statistical tests with proper error handling"""
    results = {}
    
    try:
        # Check if we have enough samples
        if len(up_data) < 2 or len(down_data) < 2:
            logger.warning(f"Not enough samples for {metric_name} statistical tests")
            return {
                'mannwhitney': None,
                'cohens_d': None,
                'error': 'Insufficient sample size'
            }
        
        # Remove NaN values
        up_data = up_data.dropna()
        down_data = down_data.dropna()
        
        # Mann-Whitney U test
        if len(up_data) >= 2 and len(down_data) >= 2:
            stat, pval = stats.mannwhitneyu(up_data, down_data, alternative='two-sided')
            results['mannwhitney'] = {
                'statistic': stat,
                'pvalue': pval
            }
            
            # Calculate Cohen's d only if we have valid data
            if up_data.std() > 0 or down_data.std() > 0:
                pooled_std = np.sqrt((up_data.std() ** 2 + down_data.std() ** 2) / 2)
                if pooled_std > 0:
                    d = (up_data.mean() - down_data.mean()) / pooled_std
                    results['cohens_d'] = d
                else:
                    results['cohens_d'] = None
                    logger.warning(f"Zero pooled standard deviation for {metric_name}")
            else:
                results['cohens_d'] = None
                logger.warning(f"No variation in data for {metric_name}")
        else:
            results['mannwhitney'] = None
            results['cohens_d'] = None
            
    except Exception as e:
        logger.error(f"Error in statistical tests for {metric_name}: {str(e)}")
        results['error'] = str(e)
        results['mannwhitney'] = None
        results['cohens_d'] = None
        
    return results

def calculate_correlation(x: pd.Series, y: pd.Series, method: str = 'spearman') -> Dict:
    """Calculate correlation between two series with improved error handling"""
    try:
        # Remove NaN values
        mask = ~(x.isna() | y.isna())
        x_clean = x[mask]
        y_clean = y[mask]
        
        # Check if we have enough data points
        if len(x_clean) < 3 or len(y_clean) < 3:
            logger.warning("Not enough valid data points for correlation")
            return None
            
        # Check for variance in both variables
        if x_clean.var() == 0 or y_clean.var() == 0:
            logger.warning("No variance in one or both variables")
            return None
            
        if method == 'spearman':
            correlation = stats.spearmanr(x_clean, y_clean)
        else:
            correlation = stats.pearsonr(x_clean, y_clean)
            
        return {
            'statistic': correlation.statistic,
            'pvalue': correlation.pvalue,
            'method': method,
            'n': len(x_clean)
        }
    except Exception as e:
        logger.warning(f"Could not calculate {method} correlation: {str(e)}")
        return None

def save_analysis_summary(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Save statistical summary of the analysis"""
    os.makedirs(os.path.join(output_dir, 'summaries'), exist_ok=True)
    
    with open(os.path.join(output_dir, 'summaries', f'{cell_type}_analysis_summary.txt'), 'w') as f:
        f.write(f"Analysis Summary for {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        # Overall statistics
        f.write("Overall Statistics:\n")
        f.write("-"*20 + "\n")
        f.write(f"Total genes analyzed: {len(df)}\n")
        f.write(f"Upregulated genes: {len(df[df['expression_status'] == 'upregulated'])}\n")
        f.write(f"Downregulated genes: {len(df[df['expression_status'] == 'downregulated'])}\n\n")
        
        # Binding type distribution
        f.write("Binding Type Distribution:\n")
        f.write("-"*20 + "\n")
        binding_counts = df['binding_type'].value_counts()
        for binding_type, count in binding_counts.items():
            f.write(f"{binding_type}: {count}\n")
        f.write("\n")
        
        # Methylation statistics by binding type and expression status
        for binding_type in df['binding_type'].unique():
            subset = df[df['binding_type'] == binding_type]
            f.write(f"\n{binding_type.upper()} BINDING:\n")
            f.write("-"*20 + "\n")
            
            for status in ['upregulated', 'downregulated']:
                status_subset = subset[subset['expression_status'] == status]
                if len(status_subset) > 0:
                    f.write(f"\n{status.capitalize()} Genes:\n")
                    f.write(f"Count: {len(status_subset)}\n")
                    
                    # Methylation statistics
                    for metric in ['promoter_methylation', 'gene_body_methylation']:
                        stats_dict = status_subset[metric].describe()
                        f.write(f"\n{metric.replace('_', ' ').title()}:\n")
                        f.write(f"Mean: {stats_dict['mean']:.2f}\n")
                        f.write(f"Std: {stats_dict['std']:.2f}\n")
                        f.write(f"Min: {stats_dict['min']:.2f}\n")
                        f.write(f"Max: {stats_dict['max']:.2f}\n")
            
            # Correlation analysis
            if len(subset) > 0:
                corr_result = calculate_correlation(
                    subset['promoter_methylation'],
                    subset['gene_body_methylation']
                )
                if corr_result:
                    print("\nCorrelation between promoter and gene body methylation:")
                    print(f"Spearman correlation: rho={corr_result['statistic']:.4f}, p={corr_result['pvalue']:.4f}")
                else:
                    print("\nCorrelation analysis: Not available (insufficient data or no variation)")

def calculate_effect_size(group1: pd.Series, group2: pd.Series) -> float:
    """Calculate Cohen's d effect size between two groups"""
    try:
        # Remove NaN values
        group1_clean = group1.dropna()
        group2_clean = group2.dropna()
        
        if len(group1_clean) < 2 or len(group2_clean) < 2:
            return None
            
        # Calculate Cohen's d
        d = (group1_clean.mean() - group2_clean.mean()) / np.sqrt(
            ((group1_clean.std() ** 2 + group2_clean.std() ** 2) / 2)
        )
        return d
    except Exception as e:
        logger.warning(f"Could not calculate effect size: {str(e)}")
        return None

#Create focused visualization of significant findings
def plot_significant_differences(results: Dict[str, pd.DataFrame], output_dir: str):
    """Create plots focusing on the significant differences found in the analysis"""
    os.makedirs(os.path.join(output_dir, 'significant_findings'), exist_ok=True)
    
    # Plot settings
    colors = {'upregulated': '#2ecc71', 'downregulated': '#e74c3c'}
    
    for cell_type, df in results.items():
        mecp2_bound_df = df[df['mecp2_bound']]
        
        # 1. NEU: Gene body methylation in both binding
        if cell_type == 'NEU':
            both_binding = mecp2_bound_df[
                (mecp2_bound_df['binding_type'] == 'both') &
                (mecp2_bound_df['expression_status'].isin(['upregulated', 'downregulated']))
            ]
            
            plt.figure(figsize=(8, 6))
            sns.boxplot(data=both_binding,
                       x='expression_status',
                       y='gene_body_methylation',
                       hue='expression_status',
                       palette=colors,
                       legend=False)
            plt.title('NEU: Gene Body Methylation\nBoth Binding (p=0.0101)')
            plt.ylabel('Gene Body Methylation (%)')
            plt.xlabel('Expression Status')
            
            # Add individual points
            sns.stripplot(data=both_binding,
                         x='expression_status',
                         y='gene_body_methylation',
                         color='black',
                         alpha=0.3,
                         jitter=0.2,
                         size=4)
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'significant_findings',
                                    'NEU_both_binding_gene_body.pdf'))
            plt.close()
        
        # 2. NSC: Gene body methylation comparisons
        elif cell_type == 'NSC':
            # Both binding
            both_binding = mecp2_bound_df[
                (mecp2_bound_df['binding_type'] == 'both') &
                (mecp2_bound_df['expression_status'].isin(['upregulated', 'downregulated']))
            ]
            
            # Endo only
            endo_only = mecp2_bound_df[
                (mecp2_bound_df['binding_type'] == 'endo_only') &
                (mecp2_bound_df['expression_status'].isin(['upregulated', 'downregulated']))
            ]
            
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # Updated boxplot calls
            sns.boxplot(data=both_binding,
                       x='expression_status',
                       y='gene_body_methylation',
                       hue='expression_status',
                       palette=colors,
                       legend=False,
                       ax=ax1)
            ax1.set_title('Both Binding (p=0.0225)')
            ax1.set_ylabel('Gene Body Methylation (%)')
            ax1.set_xlabel('Expression Status')
            
            sns.boxplot(data=endo_only,
                       x='expression_status',
                       y='gene_body_methylation',
                       hue='expression_status',
                       palette=colors,
                       legend=False,
                       ax=ax2)
            ax2.set_title('Endo-only Binding (p=0.0347)')
            ax2.set_ylabel('Gene Body Methylation (%)')
            ax2.set_xlabel('Expression Status')
            
            plt.suptitle('NSC: Gene Body Methylation Patterns', y=1.02)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'significant_findings',
                                    'NSC_gene_body_methylation.pdf'))
            plt.close()

#Modified plot_mecp2_regulatory_patterns function
def plot_mecp2_regulatory_patterns(results: Dict[str, pd.DataFrame], output_dir: str):
    """Create comprehensive visualization of MeCP2's regulatory patterns"""
    os.makedirs(os.path.join(output_dir, 'regulatory_patterns'), exist_ok=True)
    
    for cell_type, df in results.items():
        # Filter for regulated genes
        regulated_df = df[df['expression_status'].isin(['upregulated', 'downregulated'])]
        if len(regulated_df) == 0:
            logger.warning(f"No regulated genes found for {cell_type}")
            continue
            
        mecp2_bound_df = regulated_df[regulated_df['mecp2_bound']]
        if len(mecp2_bound_df) == 0:
            logger.warning(f"No MeCP2-bound regulated genes found for {cell_type}")
            continue
        
        # Create visualizations
        create_methylation_patterns_plot(mecp2_bound_df, cell_type, output_dir)
        create_binding_correlation_plot(mecp2_bound_df, cell_type, output_dir)
        create_methylation_distribution_plot(mecp2_bound_df, cell_type, output_dir)

#Create focused visualization of significant findings
def plot_significant_differences(results: Dict[str, pd.DataFrame], output_dir: str):
    """Create plots focusing on the significant differences found in the analysis"""
    os.makedirs(os.path.join(output_dir, 'significant_findings'), exist_ok=True)
    
    # Plot settings
    colors = {'upregulated': '#2ecc71', 'downregulated': '#e74c3c'}
    
    for cell_type, df in results.items():
        mecp2_bound_df = df[df['mecp2_bound']]
        
        # 1. NEU: Gene body methylation in both binding
        if cell_type == 'NEU':
            both_binding = mecp2_bound_df[
                (mecp2_bound_df['binding_type'] == 'both') &
                (mecp2_bound_df['expression_status'].isin(['upregulated', 'downregulated']))
            ]
            
            plt.figure(figsize=(8, 6))
            sns.boxplot(data=both_binding,
                       x='expression_status',
                       y='gene_body_methylation',
                       hue='expression_status',
                       palette=colors,
                       legend=False)
            plt.title('NEU: Gene Body Methylation\nBoth Binding (p=0.0101)')
            plt.ylabel('Gene Body Methylation (%)')
            plt.xlabel('Expression Status')
            
            # Add individual points
            sns.stripplot(data=both_binding,
                         x='expression_status',
                         y='gene_body_methylation',
                         color='black',
                         alpha=0.3,
                         jitter=0.2,
                         size=4)
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'significant_findings',
                                    'NEU_both_binding_gene_body.pdf'))
            plt.close()
        
        # 2. NSC: Gene body methylation comparisons
        elif cell_type == 'NSC':
            # Both binding
            both_binding = mecp2_bound_df[
                (mecp2_bound_df['binding_type'] == 'both') &
                (mecp2_bound_df['expression_status'].isin(['upregulated', 'downregulated']))
            ]
            
            # Endo only
            endo_only = mecp2_bound_df[
                (mecp2_bound_df['binding_type'] == 'endo_only') &
                (mecp2_bound_df['expression_status'].isin(['upregulated', 'downregulated']))
            ]
            
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # Updated boxplot calls
            sns.boxplot(data=both_binding,
                       x='expression_status',
                       y='gene_body_methylation',
                       hue='expression_status',
                       palette=colors,
                       legend=False,
                       ax=ax1)
            ax1.set_title('Both Binding (p=0.0225)')
            ax1.set_ylabel('Gene Body Methylation (%)')
            ax1.set_xlabel('Expression Status')
            
            sns.boxplot(data=endo_only,
                       x='expression_status',
                       y='gene_body_methylation',
                       hue='expression_status',
                       palette=colors,
                       legend=False,
                       ax=ax2)
            ax2.set_title('Endo-only Binding (p=0.0347)')
            ax2.set_ylabel('Gene Body Methylation (%)')
            ax2.set_xlabel('Expression Status')
            
            plt.suptitle('NSC: Gene Body Methylation Patterns', y=1.02)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'significant_findings',
                                    'NSC_gene_body_methylation.pdf'))
            plt.close()

#Modified print_regulated_summary function
def print_regulated_summary_v2(results: Dict[str, pd.DataFrame]):
    """Print summary statistics for regulated MeCP2-bound genes with improved analysis"""
    for cell_type, df in results.items():
        print(f"\nSummary for {cell_type}:")
        print("="*50)
        
        # First, check regulated genes
        regulated_df = df[df['expression_status'].isin(['upregulated', 'downregulated'])]
        print(f"\nTotal regulated genes: {len(regulated_df)}")
        
        if len(regulated_df) == 0:
            print("No regulated genes found.")
            continue
        
        # Check MeCP2 binding
        mecp2_bound_df = regulated_df[regulated_df['mecp2_bound']]
        print(f"MeCP2-bound regulated genes: {len(mecp2_bound_df)}")
        
        if len(mecp2_bound_df) == 0:
            print("No MeCP2-bound regulated genes found.")
            continue
        
        # Check binding types present
        binding_types = mecp2_bound_df['binding_type'].value_counts()
        print("\nBinding types present:")
        print(binding_types)
        
        # Analyze each binding type
        for binding_type in binding_types.index:
            if pd.isna(binding_type):
                continue
                
            subset_df = mecp2_bound_df[mecp2_bound_df['binding_type'] == binding_type]
            
            print(f"\n{binding_type.upper()} BINDING:")
            print("-"*20)
            
            # Gene counts and basic statistics
            expr_counts = subset_df['expression_status'].value_counts()
            print("\nGene counts by expression status:")
            print(expr_counts)
            
            # Calculate and display statistics
            metrics = {
                'promoter_methylation': 'Promoter Methylation',
                'gene_body_methylation': 'Gene Body Methylation',
                'promoter_cpg_count': 'Promoter CpG Count'
            }
            
            for metric, metric_name in metrics.items():
                print(f"\n{metric_name}:")
                
                # Calculate statistics by expression status
                stats_df = subset_df.groupby('expression_status')[metric].agg([
                    'count', 'mean', 'std', 'min', 'max'
                ]).round(3)
                print(stats_df)
                
                # Perform statistical tests
                up_genes = subset_df[subset_df['expression_status'] == 'upregulated']
                down_genes = subset_df[subset_df['expression_status'] == 'downregulated']
                
                if len(up_genes) > 0 and len(down_genes) > 0:
                    stats_results = perform_statistical_tests(
                        up_genes[metric],
                        down_genes[metric],
                        metric_name
                    )
                    
                    if stats_results.get('mannwhitney'):
                        print(f"\nMann-Whitney U test for {metric_name}:")
                        print(f"Statistic: {stats_results['mannwhitney']['statistic']:.2f}")
                        print(f"P-value: {stats_results['mannwhitney']['pvalue']:.4f}")
                        
                        if stats_results.get('cohens_d') is not None:
                            print(f"Cohen's d: {stats_results['cohens_d']:.4f}")
                        else:
                            print("Cohen's d: Not available")
                    
                    if stats_results.get('error'):
                        print(f"Error in statistical tests: {stats_results['error']}")
            
            # Calculate correlation
            if len(subset_df) > 0:
                corr_result = calculate_correlation(
                    subset_df['promoter_methylation'],
                    subset_df['gene_body_methylation']
                )
                if corr_result:
                    print("\nCorrelation between promoter and gene body methylation:")
                    print(f"Spearman correlation: rho={corr_result['statistic']:.4f}, p={corr_result['pvalue']:.4f}")
                else:
                    print("\nCorrelation analysis: Not available (insufficient data or no variation)")

#Load and process gene annotations with gene name mapping
def load_gene_annotations(gtf_file: str) -> pd.DataFrame:
    """Load gene annotations with debug mode support"""
    logger.info("Loading gene annotations...")
    
    genes = pr.read_gtf(gtf_file)
    genes = genes[genes.Feature == 'gene']
    
    # Convert to DataFrame
    genes_df = pd.DataFrame({
        'gene_id': genes.gene_id.str.split('.').str[0],
        'gene_name': genes.gene_name,
        'chr': genes.Chromosome,
        'start': genes.Start,
        'end': genes.End,
        'strand': genes.Strand
    })
    
    # In debug mode, sample a subset of genes
    if CONFIG['debug']['enabled']:
        genes_df = genes_df.sample(n=min(len(genes_df), CONFIG['debug']['sample_size']), 
                                 random_state=42)
        logger.info(f"Debug mode: sampled {len(genes_df)} genes")
    
    # Create mapping and add coordinates
    gene_name_to_id = dict(zip(genes_df['gene_name'], genes_df['gene_id']))
    
    genes_df[['promoter_start', 'promoter_end']] = genes_df.apply(
        lambda x: get_promoter_region(x.start, x.end, x.strand), 
        axis=1, result_type='expand'
    )
    
    genes_df[['gene_body_start', 'gene_body_end']] = genes_df.apply(
        lambda x: get_gene_body_region(x.start, x.end, x.strand),
        axis=1, result_type='expand'
    )
    
    return genes_df, gene_name_to_id

#Load and process differential expression data with gene ID mapping
def load_expression_data(rnaseq_files: Union[Dict[str, str], str], 
                        gene_name_to_id: dict,
                        single_file: bool = False) -> Union[Dict[str, pd.DataFrame], pd.DataFrame]:
    """Load expression data with debug mode support"""
    logger.info("Loading differential expression data...")
    
    if single_file:
        # Handle single file case
        df = pd.read_csv(rnaseq_files)
        
        # Map gene symbols to Ensembl IDs
        if 'gene' in df.columns:
            df['gene_id'] = df['gene'].map(gene_name_to_id)
        else:
            df['gene_id'] = df['gene_id'].map(gene_name_to_id)
        
        # Remove rows where mapping failed
        df = df.dropna(subset=['gene_id'])
        
        # In debug mode, sample a subset of genes
        if CONFIG['debug']['enabled']:
            df = df.sample(n=min(len(df), CONFIG['debug']['sample_size']), 
                         random_state=42)
            logger.info(f"Debug mode: sampled {len(df)} genes")
        
        # Classify genes
        df['expression_status'] = 'unchanged'
        df.loc[(df['log2FoldChange'] > CONFIG['expression_thresholds']['log2fc']) & 
               (df['padj'] < CONFIG['expression_thresholds']['padj']), 
               'expression_status'] = 'upregulated'
        df.loc[(df['log2FoldChange'] < -CONFIG['expression_thresholds']['log2fc']) & 
               (df['padj'] < CONFIG['expression_thresholds']['padj']), 
               'expression_status'] = 'downregulated'
        
        return df
    
    # Handle dictionary of files case
    expression_data = {}
    for cell_type, file_path in rnaseq_files.items():
        expression_data[cell_type] = load_expression_data(file_path, gene_name_to_id, single_file=True)
    
    return expression_data

def analyze_methylation_distribution(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze the distribution of methylation levels in different regions"""
    try:
        # Create violin plots for methylation distribution
        plt.figure(figsize=(12, 6))
        
        # Plot methylation levels for promoters and gene bodies
        data_to_plot = pd.melt(df, 
                              value_vars=['promoter_methylation', 'gene_body_methylation'],
                              var_name='Region', value_name='Methylation Level')
        
        sns.violinplot(data=data_to_plot, x='Region', y='Methylation Level')
        plt.title(f'{cell_type}: Distribution of Methylation Levels')
        plt.savefig(os.path.join(output_dir, f'{cell_type}_methylation_distribution.pdf'))
        plt.close()
        
        # Calculate and save statistics
        stats_dict = {
            'promoter': df['promoter_methylation'].describe(),
            'gene_body': df['gene_body_methylation'].describe()
        }
        
        # Calculate correlation between promoter and gene body methylation
        try:
            valid_data = df[['promoter_methylation', 'gene_body_methylation']].dropna()
            if len(valid_data) > 1:
                correlation = stats.spearmanr(valid_data['promoter_methylation'],
                                           valid_data['gene_body_methylation'])
                stats_dict['correlation'] = {
                    'rho': correlation.statistic,
                    'pvalue': correlation.pvalue,
                    'n': len(valid_data)
                }
            else:
                stats_dict['correlation'] = None
                logger.warning(f"{cell_type}: Insufficient data for correlation calculation")
        except Exception as e:
            stats_dict['correlation'] = None
            logger.warning(f"{cell_type}: Error calculating correlation: {str(e)}")
        
        return stats_dict
    
    except Exception as e:
        logger.error(f"Error in methylation distribution analysis for {cell_type}: {str(e)}")
        return None

def analyze_methylation_expression_relationship(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze how methylation levels relate to gene expression changes"""
    try:
        # Create separate plots for promoter and gene body methylation
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # 1. Promoter methylation vs expression
        sns.scatterplot(data=df, 
                       x='promoter_methylation', 
                       y='log2FoldChange',
                       hue='expression_status',
                       alpha=0.6,
                       ax=ax1)
        ax1.set_title('Promoter Methylation vs Expression Change')
        
        # 2. Gene body methylation vs expression
        sns.scatterplot(data=df,
                       x='gene_body_methylation',
                       y='log2FoldChange',
                       hue='expression_status',
                       alpha=0.6,
                       ax=ax2)
        ax2.set_title('Gene Body Methylation vs Expression Change')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{cell_type}_methylation_expression.pdf'))
        plt.close()
        
        # Calculate statistics for different expression groups
        stats = {}
        for status in ['upregulated', 'downregulated', 'unchanged']:
            subset = df[df['expression_status'] == status]
            if len(subset) > 0:
                stats[status] = {
                    'promoter_methylation': subset['promoter_methylation'].describe(),
                    'gene_body_methylation': subset['gene_body_methylation'].describe(),
                    'count': len(subset)
                }
        
        return stats
    
    except Exception as e:
        logger.error(f"Error in methylation-expression analysis for {cell_type}: {str(e)}")
        return None

def analyze_binding_patterns(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze the relationship between methylation and MeCP2 binding"""
    try:
        # Create plots for binding analysis
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. Methylation levels by binding type
        sns.boxplot(data=df, x='binding_type', y='promoter_methylation', ax=ax1)
        ax1.set_title('Promoter Methylation by Binding Type')
        
        sns.boxplot(data=df, x='binding_type', y='gene_body_methylation', ax=ax2)
        ax2.set_title('Gene Body Methylation by Binding Type')
        
        # 2. Expression changes by binding type
        sns.boxplot(data=df, x='binding_type', y='log2FoldChange', ax=ax3)
        ax3.set_title('Expression Changes by Binding Type')
        
        # 3. CpG density by binding type
        sns.boxplot(data=df, x='binding_type', y='promoter_cpg_count', ax=ax4)
        ax4.set_title('CpG Density by Binding Type')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{cell_type}_binding_patterns.pdf'))
        plt.close()
        
        # Calculate statistics
        stats = {}
        for binding_type in df['binding_type'].dropna().unique():
            if pd.isna(binding_type):
                continue
            subset = df[df['binding_type'] == binding_type]
            stats[binding_type] = {
                'count': len(subset),
                'methylation': {
                    'promoter': subset['promoter_methylation'].describe(),
                    'gene_body': subset['gene_body_methylation'].describe()
                },
                'expression': subset['log2FoldChange'].describe(),
                'cpg_density': {
                    'promoter': subset['promoter_cpg_count'].describe(),
                    'gene_body': subset['gene_body_cpg_count'].describe()
                }
            }
        
        return stats
    
    except Exception as e:
        logger.error(f"Error in binding pattern analysis for {cell_type}: {str(e)}")
        return None

def analyze_methylation_patterns_detailed(results: Dict[str, pd.DataFrame], output_dir: str):
    """Perform detailed analysis of methylation patterns and their relationship with gene regulation"""
    
    for cell_type, df in results.items():
        logger.info(f"\nAnalyzing {cell_type} methylation patterns...")
        
        try:
            # Create output directory for this analysis
            analysis_dir = os.path.join(output_dir, 'methylation_analysis')
            os.makedirs(analysis_dir, exist_ok=True)
            
            # 1. Analyze methylation levels distribution
            dist_stats = analyze_methylation_distribution(df, cell_type, analysis_dir)
            if dist_stats:
                logger.info(f"Completed methylation distribution analysis for {cell_type}")
            
            # 2. Analyze relationship between methylation and expression
            expr_stats = analyze_methylation_expression_relationship(df, cell_type, analysis_dir)
            if expr_stats:
                logger.info(f"Completed methylation-expression analysis for {cell_type}")
            
            # 3. Analyze binding patterns
            binding_stats = analyze_binding_patterns(df, cell_type, analysis_dir)
            if binding_stats:
                logger.info(f"Completed binding pattern analysis for {cell_type}")
            
            # 4. Generate comprehensive report
            generate_methylation_report(df, cell_type, analysis_dir)
            logger.info(f"Generated methylation report for {cell_type}")
            
        except Exception as e:
            logger.error(f"Error analyzing {cell_type}: {str(e)}")
            continue

def generate_methylation_report(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Generate a comprehensive report of the methylation analysis"""
    try:
        report_file = os.path.join(output_dir, f'{cell_type}_methylation_report.txt')
        
        with open(report_file, 'w') as f:
            f.write(f"Methylation Analysis Report for {cell_type}\n")
            f.write("="*50 + "\n\n")
            
            # 1. Overall statistics
            f.write("1. Overall Statistics\n")
            f.write("-"*20 + "\n")
            f.write(f"Total genes analyzed: {len(df)}\n")
            f.write(f"Genes with MeCP2 binding: {df['mecp2_bound'].sum()}\n")
            f.write("\nExpression status distribution:\n")
            f.write(str(df['expression_status'].value_counts()) + "\n\n")
            
            # 2. Methylation patterns
            f.write("2. Methylation Patterns\n")
            f.write("-"*20 + "\n")
            for region in ['promoter_methylation', 'gene_body_methylation']:
                f.write(f"\n{region.replace('_', ' ').title()}:\n")
                f.write(str(df[region].describe()) + "\n")
            
            # 3. Binding analysis
            f.write("\n3. Binding Analysis\n")
            f.write("-"*20 + "\n")
            for binding_type in df['binding_type'].dropna().unique():
                subset = df[df['binding_type'] == binding_type]
                f.write(f"\n{binding_type.upper()}:\n")
                f.write(f"Number of genes: {len(subset)}\n")
                f.write("Expression status distribution:\n")
                f.write(str(subset['expression_status'].value_counts()) + "\n")
                
                try:
                    # Calculate correlations
                    corr = stats.spearmanr(subset['promoter_methylation'],
                                         subset['log2FoldChange'])
                    f.write(f"\nPromoter methylation vs expression correlation:\n")
                    f.write(f"Spearman rho: {corr.statistic:.3f}, p-value: {corr.pvalue:.3e}\n")
                    
                    corr = stats.spearmanr(subset['gene_body_methylation'],
                                         subset['log2FoldChange'])
                    f.write(f"\nGene body methylation vs expression correlation:\n")
                    f.write(f"Spearman rho: {corr.statistic:.3f}, p-value: {corr.pvalue:.3e}\n")
                except Exception as e:
                    logger.warning(f"Could not calculate correlations for {binding_type}: {str(e)}")
            
            # 4. Additional correlation analysis
            f.write("\n4. Additional Correlations\n")
            f.write("-"*20 + "\n")
            try:
                if len(df) > 0:
                    # Promoter methylation vs expression
                    corr_result = calculate_correlation(
                        df['promoter_methylation'],
                        df['log2FoldChange']
                    )
                    if corr_result:
                        f.write("\nPromoter Methylation vs Expression:\n")
                        f.write(f"rho: {corr_result['statistic']:.4f}\n")
                        f.write(f"p-value: {corr_result['pvalue']:.4f}\n")
                        f.write(f"n: {corr_result['n']}\n")
                    
                    # Gene body methylation vs expression
                    corr_result = calculate_correlation(
                        df['gene_body_methylation'],
                        df['log2FoldChange']
                    )
                    if corr_result:
                        f.write("\nGene Body Methylation vs Expression:\n")
                        f.write(f"rho: {corr_result['statistic']:.4f}\n")
                        f.write(f"p-value: {corr_result['pvalue']:.4f}\n")
                        f.write(f"n: {corr_result['n']}\n")
            except Exception as e:
                logger.warning(f"Could not calculate additional correlations: {str(e)}")
                
    except Exception as e:
        logger.error(f"Error generating methylation report for {cell_type}: {str(e)}")
        return None

def calculate_methylation_levels(region_df: pd.DataFrame,
                               medip_dir: str,
                               cell_type: str,
                               genome_fasta: str) -> pd.DataFrame:
    """Calculate methylation levels for each region"""
    logger.info(f"Starting methylation calculation for {len(region_df)} regions")
    
    # Map cell type to file prefix
    cell_type_map = {
        'NSC': 'N',
        'NEU': 'N',
        'PP': 'PP',
        'DP': 'DP'
    }
    prefix = cell_type_map.get(cell_type)
    if not prefix:
        raise ValueError(f"Unknown cell type: {cell_type}")
    
    logger.info(f"Using file prefix '{prefix}' for cell type '{cell_type}'")
    
    # Check medip directory contents
    logger.info(f"Checking contents of medip_dir: {medip_dir}")
    files = os.listdir(medip_dir)
    logger.info(f"Files in directory:\n{files}")
    
    # Get BigWig files
    ip_file = os.path.join(medip_dir, f"Medip_{prefix}_output_r1.bw")
    input_file = os.path.join(medip_dir, f"Medip_{prefix}_input_r1.bw")
    
    logger.info("Looking for BigWig files:")
    logger.info(f"IP file: {ip_file}")
    logger.info(f"Input file: {input_file}")
    
    # Load CpG islands
    logger.info("Loading CpG islands from ../DATA/cpg_islands.bed")
    cpg_islands = pd.read_csv("../DATA/cpg_islands.bed", sep='\t',
                             names=['Chromosome', 'Start', 'End', 'id', 'type', 'cpg_count'])
    
    logger.info(f"Loaded {len(cpg_islands)} CpG islands")
    logger.info(f"Columns: {cpg_islands.columns.tolist()}")
    logger.info("Sample data:")
    print(cpg_islands.head())
    
    # Convert to PyRanges
    logger.info(f"Converting {len(cpg_islands)} CpG islands to PyRanges")
    cpg_pr = pr.PyRanges(cpg_islands)
    
    # Open BigWig files
    bw_ip = pyBigWig.open(ip_file)
    bw_input = pyBigWig.open(input_file)
    
    # Get valid chromosomes
    valid_chroms = set(bw_ip.chroms().keys())
    logger.info(f"Valid chromosomes in BigWig: {valid_chroms}")
    
    # Initialize results
    results = []
    
    # Process each region with progress bar
    for idx, row in tqdm(region_df.iterrows(), total=len(region_df), desc="Analyzing regions"):
        try:
            # Ensure chromosome name is correct format
            chrom = row['chr'] if 'chr' in row else row['Chromosome']  # Handle both column names
            if not chrom.startswith('chr'):
                chrom = f"chr{chrom}"
                
            if chrom not in valid_chroms:
                continue
                
            # Get region coordinates
            start = int(row['Start'])
            end = int(row['End'])
            
            # Find overlapping CpG islands
            region_pr = pr.PyRanges(chromosomes=[chrom], starts=[start], ends=[end])
            overlaps = region_pr.join(cpg_pr)
            
            if overlaps is None or len(overlaps) == 0:
                results.append({
                    'cpg_island_count': 0,
                    'cpg_island_mean_methylation': 0.0,
                    'cpg_island_methylation_std': 0.0
                })
                continue
            
            # Calculate methylation for each overlapping CpG island
            methylation_values = []
            for _, cpg in overlaps.as_df().iterrows():
                try:
                    ip_values = bw_ip.values(chrom, cpg.Start, cpg.End)
                    input_values = bw_input.values(chrom, cpg.Start, cpg.End)
                    
                    # Calculate methylation
                    methylation = calculate_normalized_methylation(ip_values, input_values)
                    methylation_values.append(methylation)
                except Exception as e:
                    logger.error(f"Error calculating methylation for CpG island: {str(e)}")
                    continue
            
            # Store results
            results.append({
                'cpg_island_count': len(overlaps),
                'cpg_island_mean_methylation': np.mean(methylation_values) if methylation_values else 0.0,
                'cpg_island_methylation_std': np.std(methylation_values) if len(methylation_values) > 1 else 0.0
            })
            
        except Exception as e:
            logger.error(f"Error processing region {idx}: {str(e)}")
            results.append({
                'cpg_island_count': 0,
                'cpg_island_mean_methylation': 0.0,
                'cpg_island_methylation_std': 0.0
            })
    
    # Close BigWig files
    bw_ip.close()
    bw_input.close()
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    # Log completion
    logger.info("Completed methylation calculation:")
    logger.info(f"- Processed {len(region_df)} regions")
    logger.info(f"- Found {results_df['cpg_island_count'].sum()} CpG islands")
    logger.info(f"- Average methylation: {results_df['cpg_island_mean_methylation'].mean():.2f}%")
    logger.info(f"- Data shape: {results_df.shape}")
    logger.info(f"- Columns: {results_df.columns.tolist()}")
    
    # Verify data
    logger.info("Verifying methylation data...")
    logger.info(f"Methylation data columns: {results_df.columns.tolist()}")
    logger.info("Sample of methylation data:")
    print(results_df.head())
    
    return results_df

def validate_methylation_levels(df: pd.DataFrame) -> pd.DataFrame:
    """
    Validate and clean methylation level calculations
    
    Args:
        df: DataFrame with methylation data
    
    Returns:
        Validated and cleaned DataFrame
    """
    # Create copy to avoid modifying original
    df = df.copy()
    
    # Ensure methylation values are within valid range
    df['methylation'] = df['methylation'].clip(0, 100)
    
    # Flag low-confidence measurements
    df['low_confidence'] = (
        (df['cpg_count'] < CONFIG['quality_thresholds']['min_cpgs']) |
        (df['methylation_confidence'] < CONFIG['quality_thresholds']['coverage'])
    )
    
    # Log validation results
    total_regions = len(df)
    low_conf_regions = df['low_confidence'].sum()
    logger.info(f"Methylation validation results:")
    logger.info(f"Total regions: {total_regions}")
    logger.info(f"Low confidence regions: {low_conf_regions} ({low_conf_regions/total_regions*100:.1f}%)")
    
    return df

def add_methylation_qc_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Enhanced quality control metrics"""
    df = df.copy()
    
    try:
        # Calculate region lengths and CpG density
        if 'start' in df.columns and 'end' in df.columns:
            region_lengths = df['end'] - df['start']
        else:
            logger.warning("No start/end columns found. Using default region length.")
            region_lengths = 1000  # default length
            
        # Calculate improved metrics
        df['region_length'] = region_lengths
        df['cpg_density'] = df['cpg_count'] / df['region_length']
        
        # Add coverage uniformity if available
        try:
            df['coverage_uniformity'] = calculate_coverage_uniformity(df)
        except Exception as e:
            logger.warning(f"Could not calculate coverage uniformity: {str(e)}")
            df['coverage_uniformity'] = 1.0  # default value
            
        # Add signal-to-noise ratio if available
        try:
            df['signal_to_noise'] = calculate_signal_to_noise(df)
        except Exception as e:
            logger.warning(f"Could not calculate signal-to-noise ratio: {str(e)}")
            df['signal_to_noise'] = 2.0  # default value
        
        # Calculate composite quality score
        df['methylation_confidence'] = df.apply(
            lambda x: calculate_enhanced_confidence_score(
                x['methylation'],
                x['cpg_density'],
                x['cpg_count'],
                x.get('coverage_uniformity', 1.0),
                x.get('signal_to_noise', 2.0)
            ),
            axis=1
        )
        
        return df
        
    except Exception as e:
        logger.error(f"Error in add_methylation_qc_metrics: {str(e)}")
        logger.error(f"DataFrame columns: {df.columns.tolist()}")
        logger.error(f"DataFrame head:\n{df.head().to_string()}")
        raise

def calculate_enhanced_confidence_score(methylation: float, 
                                     cpg_density: float, 
                                     cpg_count: int,
                                     coverage_uniformity: float,
                                     signal_to_noise: float) -> float:
    """Calculate confidence score for methylation measurements"""
    # Minimum requirements
    min_cpgs = 3
    min_density = 0.001  # 1 CpG per 1000bp
    
    # Calculate component scores
    cpg_score = min(1.0, cpg_count / min_cpgs)
    density_score = min(1.0, cpg_density / min_density)
    coverage_score = coverage_uniformity
    signal_score = min(1.0, (signal_to_noise - 1) / 2)  # Normalize to [0,1]
    
    # Weight the components
    confidence = (
        0.4 * cpg_score +          # CpG count is important
        0.3 * density_score +      # Density matters
        0.2 * coverage_score +     # Coverage uniformity
        0.1 * signal_score        # Signal-to-noise ratio
    )
    
    return max(0.0, min(1.0, confidence))

def map_binding_to_genes(mecp2_binding: pd.DataFrame, genes_df: pd.DataFrame) -> pd.DataFrame:
    """Map MeCP2 binding regions to genes based on genomic coordinates"""
    try:
        # Standardize chromosome format
        mecp2_binding = mecp2_binding.copy()
        genes_df = genes_df.copy()
        
        # Ensure chr prefix consistency
        mecp2_binding['chr'] = mecp2_binding['chr'].str.replace('^chr', '', regex=True)
        genes_df['chr'] = genes_df['chr'].str.replace('^chr', '', regex=True)
        
        # Create window around genes
        window_size = 5000  # 5kb window
        genes_df['promoter_start'] = genes_df['promoter_start'] - window_size
        genes_df['gene_body_end'] = genes_df['gene_body_end'] + window_size
        
        # Debug coordinate ranges
        logger.info("\nCoordinate ranges after standardization:")
        for chrom in sorted(set(mecp2_binding['chr'].unique()) | set(genes_df['chr'].unique())):
            binding_regions = mecp2_binding[mecp2_binding['chr'] == chrom]
            genes = genes_df[genes_df['chr'] == chrom]
            if not binding_regions.empty and not genes.empty:
                logger.info(f"\nChromosome {chrom}:")
                logger.info(f"Binding regions: {len(binding_regions)}, range: {binding_regions['start'].min()}-{binding_regions['end'].max()}")
                logger.info(f"Genes: {len(genes)}, range: {genes['promoter_start'].min()}-{genes['gene_body_end'].max()}")
        
        # Create PyRanges objects
        binding_pr = pr.PyRanges(
            chromosomes=mecp2_binding['chr'],
            starts=mecp2_binding['start'].astype(int),
            ends=mecp2_binding['end'].astype(int)
        )
        
        genes_pr = pr.PyRanges(
            chromosomes=genes_df['chr'],
            starts=genes_df['promoter_start'].astype(int),
            ends=genes_df['gene_body_end'].astype(int)
        )
        
        # Add metadata
        for col in ['binding_type', 'exo_signal', 'endo_signal']:
            if col in mecp2_binding.columns:
                setattr(binding_pr, col, mecp2_binding[col].values)
        
        genes_pr.gene_id = genes_df['gene_id'].values
        
        # Find overlaps
        logger.info("\nFinding overlaps...")
        overlaps = binding_pr.join(genes_pr)
        
        if overlaps is None or len(overlaps) == 0:
            logger.warning("No overlaps found between binding regions and genes")
            return pd.DataFrame(columns=['gene_id', 'binding_type', 'exo_signal', 'endo_signal'])
        
        # Convert to DataFrame
        binding_genes = overlaps.as_df()
        logger.info(f"Found {len(binding_genes)} overlaps")
        
        # Group by gene_id and aggregate binding info
        gene_binding = binding_genes.groupby('gene_id').agg({
            'binding_type': lambda x: x.iloc[0] if len(set(x)) == 1 else 'both',
            'exo_signal': 'max',
            'endo_signal': 'max'
        }).reset_index()
        
        return gene_binding
        
    except Exception as e:
        logger.error(f"Error in map_binding_to_genes: {str(e)}")
        raise

def validate_and_merge_data(genes_df: pd.DataFrame, expr_df: pd.DataFrame, 
                          mecp2_binding: pd.DataFrame) -> pd.DataFrame:
    """Validate and merge gene annotations with expression data"""
    try:
        # Merge gene annotations with expression data
        merged_df = genes_df.merge(expr_df, on='gene_id', how='inner')
        logger.info(f"Merged gene annotations with expression data. Shape: {merged_df.shape}")
        
        # Debug chromosome formats
        logger.info("\nChromosome format check:")
        logger.info(f"Genes chromosomes: {merged_df['chr'].head().tolist()}")
        logger.info(f"Binding chromosomes: {mecp2_binding['chr'].head().tolist()}")
        
        # Standardize chromosome format
        merged_df = merged_df.copy()
        mecp2_binding = mecp2_binding.copy()
        
        # Ensure 'chr' prefix consistency
        merged_df['chr'] = merged_df['chr'].astype(str)
        mecp2_binding['chr'] = mecp2_binding['chr'].astype(str)
        
        if not merged_df['chr'].str.startswith('chr').all():
            merged_df['chr'] = 'chr' + merged_df['chr']
        if not mecp2_binding['chr'].str.startswith('chr').all():
            mecp2_binding['chr'] = 'chr' + mecp2_binding['chr']
            
        logger.info("\nAfter standardization:")
        logger.info(f"Genes chromosomes: {merged_df['chr'].head().tolist()}")
        logger.info(f"Binding chromosomes: {mecp2_binding['chr'].head().tolist()}")
        
        # Map binding regions to genes
        binding_data = map_binding_to_genes(mecp2_binding, merged_df)
        
        # Add binding information to merged data
        if len(binding_data) > 0:
            merged_df = merged_df.merge(binding_data, on='gene_id', how='left')
            merged_df['mecp2_bound'] = ~merged_df['binding_type'].isna()
        else:
            merged_df['mecp2_bound'] = False
            merged_df['binding_type'] = None
            merged_df['exo_signal'] = None
            merged_df['endo_signal'] = None
        
        # Debug final data
        logger.info("\nFinal data check:")
        logger.info(f"Total genes: {len(merged_df)}")
        logger.info(f"MeCP2-bound genes: {merged_df['mecp2_bound'].sum()}")
        if merged_df['mecp2_bound'].any():
            logger.info("Sample of bound genes:")
            logger.info(merged_df[merged_df['mecp2_bound']][
                ['gene_id', 'chr', 'binding_type']
            ].head())
        
        return merged_df
        
    except Exception as e:
        logger.error(f"Error in validate_and_merge_data: {str(e)}")
        logger.error(f"MeCP2 binding data shape: {mecp2_binding.shape}")
        logger.error(f"MeCP2 binding columns: {mecp2_binding.columns.tolist()}")
        raise

def calculate_contextual_methylation(merged_df: pd.DataFrame,
                                   medip_dir: str,
                                   cell_type: str,
                                   genome_fasta: str,
                                   n_processes: int = None) -> pd.DataFrame:
    """Parallel version of contextual methylation calculation"""
    logger.info("Calculating methylation levels in parallel...")

    # Prepare promoter and gene body data
    promoter_data = merged_df[['chr', 'promoter_start', 'promoter_end']].rename(
        columns={'promoter_start': 'start', 'promoter_end': 'end'}
    )
    gene_body_data = merged_df[['chr', 'gene_body_start', 'gene_body_end']].rename(
        columns={'gene_body_start': 'start', 'gene_body_end': 'end'}
    )

    # Calculate both in parallel using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=2) as executor:
        promoter_future = executor.submit(
            calculate_methylation_levels_parallel,
            promoter_data, medip_dir, cell_type, genome_fasta, n_processes
        )
        gene_body_future = executor.submit(
            calculate_methylation_levels_parallel,
            gene_body_data, medip_dir, cell_type, genome_fasta, n_processes
        )

        promoter_methylation = promoter_future.result()
        gene_body_methylation = gene_body_future.result()

    # Create result DataFrame
    result_df = pd.DataFrame({
        'gene_id': merged_df['gene_id'],
        'gene_name': merged_df['gene_name'],
        'expression_status': merged_df['expression_status'],
        'log2FoldChange': merged_df['log2FoldChange'],
        'padj': merged_df['padj'],
        'promoter_methylation': promoter_methylation['methylation'],
        'promoter_cpg_count': promoter_methylation['cpg_count'],
        'gene_body_methylation': gene_body_methylation['methylation'],
        'gene_body_cpg_count': gene_body_methylation['cpg_count']
    })

    # Add MeCP2 binding information
    for col in ['mecp2_bound', 'binding_type', 'exo_signal', 'endo_signal']:
        if col in merged_df.columns:
            result_df[col] = merged_df[col]
        else:
            result_df[col] = None if col != 'mecp2_bound' else False

    return result_df

def add_genomic_context(df: pd.DataFrame) -> pd.DataFrame:
    """Add genomic context information to methylation data"""
    df = df.copy()
    
    # Add CpG density categories
    df['promoter_cpg_density'] = df['promoter_cpg_count'] / (
        CONFIG['genomic_regions']['promoter'][1] - CONFIG['genomic_regions']['promoter'][0]
    )
    
    # Classify CpG density
    df['promoter_cpg_class'] = pd.cut(
        df['promoter_cpg_density'],
        bins=[0, 0.02, 0.05, np.inf],
        labels=['Low', 'Medium', 'High']
    )
    
    return df

def validate_methylation_patterns(df: pd.DataFrame) -> pd.DataFrame:
    """Validate methylation patterns against known biological expectations"""
    df = df.copy()
    
    # Flag unexpected patterns
    df['unexpected_pattern'] = (
        ((df['expression_status'] == 'upregulated') & 
         (df['promoter_methylation'] > CONFIG['methylation_thresholds']['hyper'])) |
        ((df['expression_status'] == 'downregulated') & 
         (df['promoter_methylation'] < CONFIG['methylation_thresholds']['hypo']))
    )
    
    if df['unexpected_pattern'].any():
        logger.warning(
            f"Found {df['unexpected_pattern'].sum()} genes with unexpected "
            "methylation patterns"
        )
    
    return df

def create_binding_pattern_plots(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Create plots for binding patterns with empty data handling"""
    try:
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Debug info about input data
        logger.info(f"\nBinding pattern analysis debug for {cell_type}:")
        logger.info(f"Total genes: {len(df)}")
        logger.info(f"Columns available: {df.columns.tolist()}")
        logger.info(f"MeCP2 bound genes: {df['mecp2_bound'].sum()}")
        logger.info(f"Expression status distribution:\n{df['expression_status'].value_counts()}")
        if 'binding_type' in df.columns:
            logger.info(f"Binding type distribution:\n{df['binding_type'].value_counts(dropna=False)}")
        
        if not df['mecp2_bound'].any():
            logger.warning(f"No MeCP2-bound genes found for {cell_type}")
            
            # Create empty plots with message
            for plot_type in ['binding_distribution', 'methylation_vs_binding']:
                plt.figure(figsize=(10, 6))
                plt.text(0.5, 0.5, 'No MeCP2-bound genes found', 
                        horizontalalignment='center',
                        verticalalignment='center')
                plt.savefig(os.path.join(output_dir, f'{cell_type}_{plot_type}.pdf'))
                plt.close()
            return
        
        # Get bound genes
        bound_genes = df[df['mecp2_bound']]
        logger.info(f"\nBound genes data:")
        logger.info(f"Number of bound genes: {len(bound_genes)}")
        logger.info(f"Expression status of bound genes:\n{bound_genes['expression_status'].value_counts()}")
        logger.info(f"Binding types of bound genes:\n{bound_genes['binding_type'].value_counts(dropna=False)}")
        
        # 1. Binding distribution by expression status
        plt.figure(figsize=(10, 6))
        if len(bound_genes) > 0:
            # Debug info before plotting
            logger.info("\nPlotting data summary:")
            for status in bound_genes['expression_status'].unique():
                for btype in bound_genes['binding_type'].unique():
                    count = len(bound_genes[(bound_genes['expression_status'] == status) & 
                                         (bound_genes['binding_type'] == btype)])
                    logger.info(f"Status: {status}, Type: {btype}, Count: {count}")
            
            try:
                sns.countplot(data=bound_genes, x='expression_status', hue='binding_type')
                plt.title(f'{cell_type}: MeCP2 Binding Distribution')
            except Exception as e:
                logger.error(f"Error in countplot: {str(e)}")
                logger.error("Data sample for debugging:")
                logger.error(bound_genes[['expression_status', 'binding_type']].head())
        else:
            plt.text(0.5, 0.5, 'No binding distribution data available')
        
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{cell_type}_binding_distribution.pdf'))
        plt.close()
        
        # 2. Methylation vs Binding Signal
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        if len(bound_genes) > 0:
            # Debug methylation and signal data
            logger.info("\nMethylation and signal data summary:")
            for col in ['promoter_methylation', 'gene_body_methylation', 'exo_signal', 'endo_signal']:
                if col in bound_genes.columns:
                    logger.info(f"{col} stats:")
                    logger.info(f"Mean: {bound_genes[col].mean():.2f}")
                    logger.info(f"Range: {bound_genes[col].min():.2f} - {bound_genes[col].max():.2f}")
                    logger.info(f"Null values: {bound_genes[col].isnull().sum()}")
            
            try:
                sns.scatterplot(data=bound_genes, x='promoter_methylation', 
                              y='exo_signal', ax=ax1)
                ax1.set_title('Promoter Methylation vs Exo Signal')
                
                sns.scatterplot(data=bound_genes, x='gene_body_methylation',
                              y='endo_signal', ax=ax2)
                ax2.set_title('Gene Body Methylation vs Endo Signal')
            except Exception as e:
                logger.error(f"Error in scatterplot: {str(e)}")
                logger.error("Data sample for debugging:")
                logger.error(bound_genes[['promoter_methylation', 'gene_body_methylation', 
                                       'exo_signal', 'endo_signal']].head())
        else:
            ax1.text(0.5, 0.5, 'No methylation-binding data available',
                    horizontalalignment='center')
            ax2.text(0.5, 0.5, 'No methylation-binding data available',
                    horizontalalignment='center')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{cell_type}_methylation_vs_binding.pdf'))
        plt.close()
        
    except Exception as e:
        logger.error(f"Error in binding pattern analysis for {cell_type}: {str(e)}")
        logger.error("Full error info:", exc_info=True)

def process_chunk(chunk: pd.DataFrame, medip_dir: str, cell_type_prefix: str, genome_fasta: str) -> pd.DataFrame:
    """Process a chunk of regions for methylation calculation"""
    return calculate_methylation_levels_with_replicates(
        chunk, medip_dir, cell_type_prefix, genome_fasta
    )

def calculate_methylation_levels_parallel(region_df: pd.DataFrame,
                                        medip_dir: str,
                                        cell_type_prefix: str,
                                        genome_fasta: str,
                                        n_processes: int = None) -> pd.DataFrame:
    """Parallel version of calculate_methylation_levels"""
    if n_processes is None:
        n_processes = max(1, multiprocessing.cpu_count() - 1)

    # Split the data into chunks
    chunk_size = max(1, len(region_df) // (n_processes * 4))  # Create more chunks than processes
    chunks = np.array_split(region_df, len(region_df) // chunk_size + 1)

    # Process chunks in parallel using the dedicated function
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = [
            executor.submit(process_chunk, chunk, medip_dir, cell_type_prefix, genome_fasta)
            for chunk in chunks
        ]
        
        # Process results with progress bar
        results = []
        for future in tqdm(futures, desc="Calculating methylation levels"):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                logger.error(f"Error processing chunk: {str(e)}")
                continue

    # Combine results
    if not results:
        raise RuntimeError("No valid results obtained from parallel processing")
    
    return pd.concat(results, ignore_index=True)

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

def analyze_binding_enrichment_groups(results, output_dir):
    """Analyze binding enrichment patterns across different groups"""
    
    # Handle the case where results is a dictionary of cell types
    if isinstance(results, dict):
        all_results = {}
        for cell_type, df in results.items():
            cell_results = analyze_single_cell_type(df, cell_type, output_dir)
            if cell_results:
                all_results[cell_type] = cell_results
        return all_results
    else:
        # Handle single DataFrame case
        return analyze_single_cell_type(results, 'combined', output_dir)

def analyze_single_cell_type(df, cell_type, output_dir):
    """Analyze binding enrichment for a single cell type"""
    results = {}
    
    # Check required columns - updated to match available columns
    required_columns = ['binding_type', 'expression_status', 
                       'promoter_methylation', 'gene_body_methylation']
    
    if not all(col in df.columns for col in required_columns):
        logger.warning(f"Missing required columns for {cell_type}. "
                      f"Available columns: {df.columns.tolist()}")
        return None
    
    # Define binding groups using binding_type column
    binding_groups = {
        'all_genes': df,
        'mecp2_bound': df[df['mecp2_bound']],
        'exo_enriched': df[df['binding_type'].isin(['exo', 'both'])],
        'endo_only': df[df['binding_type'] == 'endo'],
        'non_enriched': df[~df['mecp2_bound']]
    }
    
    # Analyze each group
    for group_name, group_df in binding_groups.items():
        if len(group_df) == 0:
            logger.warning(f"No genes found in group: {group_name} for {cell_type}")
            continue
            
        # Calculate regulation counts
        regulation_counts = group_df['expression_status'].value_counts()
        
        # Calculate methylation statistics
        methylation_stats = group_df.groupby('expression_status').agg({
            'promoter_methylation': ['mean', 'std', 'count'],
            'gene_body_methylation': ['mean', 'std', 'count']
        })
        
        # Store results
        results[group_name] = {
            'regulation_counts': regulation_counts,
            'methylation_stats': methylation_stats
        }
        
        # Create visualizations
        create_group_visualizations(group_df, group_name, cell_type, output_dir)
    
    # Save summary to file
    summary_file = os.path.join(output_dir, 'binding_enrichment', cell_type,
                               f'binding_enrichment_summary.txt')
    os.makedirs(os.path.dirname(summary_file), exist_ok=True)
    
    with open(summary_file, 'w') as f:
        f.write(f"Binding Enrichment Analysis Summary - {cell_type}\n")
        f.write("="*40 + "\n\n")
        
        for group_name, group_results in results.items():
            f.write(f"\n{group_name.upper()}\n")
            f.write("-"*30 + "\n")
            f.write("Regulation distribution:\n")
            f.write(group_results['regulation_counts'].to_string())
            f.write("\n\nMethylation statistics:\n")
            f.write(group_results['methylation_stats'].to_string())
            f.write("\n" + "="*40 + "\n")
    
    return results

def create_group_visualizations(df, group_name, cell_type, output_dir):
    """Create visualizations for binding enrichment group analysis"""
    plot_dir = os.path.join(output_dir, 'binding_enrichment', cell_type, 'plots')
    os.makedirs(plot_dir, exist_ok=True)
    
    # Set style for better readability
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.2)
    
    # 1. Combined methylation visualization
    plt.figure(figsize=(12, 6))
    
    # Reshape data for paired visualization
    plot_data = pd.melt(
        df,
        value_vars=['promoter_methylation', 'gene_body_methylation'],
        id_vars=['expression_status'],
        var_name='Region',
        value_name='Methylation'
    )
    
    # Create violin plot with individual points
    sns.violinplot(
        data=plot_data,
        x='expression_status',
        y='Methylation',
        hue='Region',
        split=True,
        inner='quartile',
        cut=0
    )
    
    plt.title(f'{cell_type} - {group_name}\nMethylation Patterns')
    plt.xlabel('Expression Status')
    plt.ylabel('Methylation Level (%)')
    plt.xticks(rotation=45)
    
    # Adjust legend
    plt.legend(title='Region', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, f'{cell_type}_{group_name}_methylation.pdf'),
                bbox_inches='tight', dpi=300)
    plt.close()
    
    # 2. Regulation distribution - horizontal bars for better readability
    plt.figure(figsize=(8, 6))
    regulation_counts = df['expression_status'].value_counts()
    
    # Create horizontal bar plot with custom colors
    colors = {'upregulated': '#2ecc71', 'downregulated': '#e74c3c', 'unchanged': '#3498db'}
    bars = plt.barh(regulation_counts.index, regulation_counts.values,
                   color=[colors.get(x, '#95a5a6') for x in regulation_counts.index])
    
    # Add value labels with offset instead of padding
    for bar in bars:
        width = bar.get_width()
        plt.text(width + width * 0.02,  # Add 2% of width as offset
                bar.get_y() + bar.get_height()/2,
                f'{int(width):,}',
                ha='left', va='center',
                fontsize=10,
                fontweight='bold')
    
    plt.title(f'{cell_type} - {group_name}\nGene Regulation Distribution')
    plt.xlabel('Number of Genes')
    
    # Add some padding to the x-axis to prevent label cutoff
    plt.margins(x=0.2)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, f'{cell_type}_{group_name}_regulation.pdf'),
                bbox_inches='tight', dpi=300)
    plt.close()

def save_gene_lists(df: pd.DataFrame, group_name: str, cell_type: str, output_dir: str):
    """Save gene lists for each regulation category"""
    list_dir = os.path.join(output_dir, 'binding_enrichment', 'gene_lists')
    os.makedirs(list_dir, exist_ok=True)
    
    for status in ['unchanged', 'upregulated', 'downregulated']:
        status_df = df[df['expression_status'] == status]
        if len(status_df) > 0:
            output_file = os.path.join(list_dir, 
                                     f'{cell_type}_{group_name}_{status}_genes.csv')
            status_df.to_csv(output_file, index=False)

def save_summary_stats(summary: Dict, cell_type: str, output_dir: str):
    """Save summary statistics for each binding group"""
    summary_dir = os.path.join(output_dir, 'binding_enrichment', 'summaries')
    os.makedirs(summary_dir, exist_ok=True)
    
    with open(os.path.join(summary_dir, f'{cell_type}_binding_summary.txt'), 'w') as f:
        f.write(f"Binding Enrichment Analysis Summary for {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        for group_name, stats in summary.items():
            f.write(f"\n{group_name.upper()}\n")
            f.write("-"*30 + "\n")
            
            f.write("\nRegulation Distribution:\n")
            for status, count in stats['regulation_counts'].items():
                f.write(f"{status}: {count}\n")
            
            f.write("\nMethylation Statistics:\n")
            for status, meth_stats in stats['methylation_stats'].items():
                f.write(f"\n{status.capitalize()}:\n")
                f.write(f"Count: {meth_stats['count']}\n")
                f.write("Promoter methylation: "
                       f"mean={meth_stats['promoter_methylation']['mean']:.2f} ± "
                       f"{meth_stats['promoter_methylation']['std']:.2f}\n")
                f.write("Gene body methylation: "
                       f"mean={meth_stats['gene_body_methylation']['mean']:.2f} ± "
                       f"{meth_stats['gene_body_methylation']['std']:.2f}\n")

def create_exo_up_visualizations(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Create detailed visualizations for exo-bound upregulated genes"""
    plot_dir = os.path.join(output_dir, 'exo_upregulated', 'plots')
    os.makedirs(plot_dir, exist_ok=True)
    
    # 1. Methylation distribution
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Promoter methylation distribution
    sns.histplot(data=df, x='promoter_methylation', bins=30, ax=ax1)
    ax1.set_title(f'{cell_type} - Promoter Methylation\nExo-bound Upregulated Genes')
    ax1.set_xlabel('Methylation Level (%)')
    ax1.set_ylabel('Count')
    
    # Gene body methylation distribution
    sns.histplot(data=df, x='gene_body_methylation', bins=30, ax=ax2)
    ax2.set_title(f'{cell_type} - Gene Body Methylation\nExo-bound Upregulated Genes')
    ax2.set_xlabel('Methylation Level (%)')
    ax2.set_ylabel('Count')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, f'{cell_type}_methylation_distribution.pdf'))
    plt.close()
    
    # 2. Methylation vs Expression
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df, x='promoter_methylation', y='log2FoldChange')
    plt.title(f'{cell_type} - Methylation vs Expression\nExo-bound Upregulated Genes')
    plt.xlabel('Promoter Methylation (%)')
    plt.ylabel('log2 Fold Change')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, f'{cell_type}_methylation_vs_expression.pdf'))
    plt.close()
    
    # 3. CpG density analysis
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Promoter CpG density vs methylation
    sns.scatterplot(data=df, x='promoter_cpg_count', y='promoter_methylation', ax=ax1)
    ax1.set_title('Promoter CpG Count vs Methylation')
    ax1.set_xlabel('CpG Count')
    ax1.set_ylabel('Methylation Level (%)')
    
    # Gene body CpG density vs methylation
    sns.scatterplot(data=df, x='gene_body_cpg_count', y='gene_body_methylation', ax=ax2)
    ax2.set_title('Gene Body CpG Count vs Methylation')
    ax2.set_xlabel('CpG Count')
    ax2.set_ylabel('Methylation Level (%)')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, f'{cell_type}_cpg_density_analysis.pdf'))
    plt.close()

def analyze_binding_enrichment_detailed(results: Dict[str, pd.DataFrame], output_dir: str):
    """Perform detailed analysis of binding enrichment patterns"""
    for cell_type, df in results.items():
        # Create binding groups
        binding_groups = {
            'exo_enriched': df[df['binding_type'].isin(['exo', 'both'])],
            'exo_only': df[df['binding_type'] == 'exo'],
            'endo_only': df[df['binding_type'] == 'endo'],
            'non_enriched': df[~df['mecp2_bound']]
        }
        
        # Analyze each group
        for group_name, group_df in binding_groups.items():
            if len(group_df) == 0:
                continue
                
            # Split by regulation status
            regulation_groups = {
                status: group_df[group_df['expression_status'] == status]
                for status in ['unchanged', 'upregulated', 'downregulated']
            }
            
            # Calculate statistics
            stats = calculate_group_statistics(regulation_groups)
            
            # Save detailed results
            save_detailed_results(
                stats, 
                regulation_groups,
                cell_type,
                group_name,
                output_dir
            )

def calculate_group_statistics(regulation_groups: Dict[str, pd.DataFrame]) -> Dict:
    """Calculate detailed statistics for each regulation group"""
    stats = {}
    
    for status, group_df in regulation_groups.items():
        if len(group_df) == 0:
            continue
            
        stats[status] = {
            'count': len(group_df),
            'promoter_methylation': {
                'mean': group_df['promoter_methylation'].mean(),
                'std': group_df['promoter_methylation'].std(),
                'median': group_df['promoter_methylation'].median(),
                'q1': group_df['promoter_methylation'].quantile(0.25),
                'q3': group_df['promoter_methylation'].quantile(0.75)
            },
            'gene_body_methylation': {
                'mean': group_df['gene_body_methylation'].mean(),
                'std': group_df['gene_body_methylation'].std(),
                'median': group_df['gene_body_methylation'].median(),
                'q1': group_df['gene_body_methylation'].quantile(0.25),
                'q3': group_df['gene_body_methylation'].quantile(0.75)
            },
            'cpg_density': {
                'promoter_mean': group_df['promoter_cpg_count'].mean(),
                'gene_body_mean': group_df['gene_body_cpg_count'].mean()
            },
            'expression': {
                'mean_log2fc': group_df['log2FoldChange'].mean(),
                'median_log2fc': group_df['log2FoldChange'].median()
            }
        }
    
    return stats

def save_detailed_results(stats: Dict,
                         regulation_groups: Dict[str, pd.DataFrame],
                         cell_type: str,
                         group_name: str,
                         output_dir: str):
    """Save detailed analysis results"""
    # Create output directory
    results_dir = os.path.join(output_dir, 'binding_enrichment', 'detailed_results')
    os.makedirs(results_dir, exist_ok=True)
    
    # Save statistics
    with open(os.path.join(results_dir, f'{cell_type}_{group_name}_stats.txt'), 'w') as f:
        f.write(f"Detailed Analysis Results - {cell_type} - {group_name}\n")
        f.write("="*50 + "\n\n")
        
        for status, stat in stats.items():
            f.write(f"\n{status.upper()}\n")
            f.write("-"*30 + "\n")
            f.write(f"Number of genes: {stat['count']}\n\n")
            
            f.write("Promoter Methylation:\n")
            f.write(f"  Mean ± SD: {stat['promoter_methylation']['mean']:.2f} ± {stat['promoter_methylation']['std']:.2f}\n")
            f.write(f"  Median (Q1-Q3): {stat['promoter_methylation']['median']:.2f} ")
            f.write(f"({stat['promoter_methylation']['q1']:.2f}-{stat['promoter_methylation']['q3']:.2f})\n\n")
            
            f.write("Gene Body Methylation:\n")
            f.write(f"  Mean ± SD: {stat['gene_body_methylation']['mean']:.2f} ± {stat['gene_body_methylation']['std']:.2f}\n")
            f.write(f"  Median (Q1-Q3): {stat['gene_body_methylation']['median']:.2f} ")
            f.write(f"({stat['gene_body_methylation']['q1']:.2f}-{stat['gene_body_methylation']['q3']:.2f})\n\n")
            
            f.write("CpG Density:\n")
            f.write(f"  Promoter mean: {stat['cpg_density']['promoter_mean']:.2f}\n")
            f.write(f"  Gene body mean: {stat['cpg_density']['gene_body_mean']:.2f}\n\n")
            
            f.write("Expression:\n")
            f.write(f"  Mean log2FC: {stat['expression']['mean_log2fc']:.2f}\n")
            f.write(f"  Median log2FC: {stat['expression']['median_log2fc']:.2f}\n\n")
    
    # Save gene lists
    for status, group_df in regulation_groups.items():
        if len(group_df) > 0:
            output_file = os.path.join(results_dir, 
                                     f'{cell_type}_{group_name}_{status}_genes_detailed.csv')
            group_df.to_csv(output_file, index=False)

def create_additional_visualizations(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Create additional visualizations for binding enrichment analysis"""
    plot_dir = os.path.join(output_dir, 'binding_enrichment', 'plots', 'additional')
    os.makedirs(plot_dir, exist_ok=True)
    
    # 1. Combined binding and methylation heatmap
    plt.figure(figsize=(12, 8))
    pivot_data = df.pivot_table(
        values=['promoter_methylation', 'gene_body_methylation'],
        index='binding_type',
        columns='expression_status',
        aggfunc='mean'
    )
    sns.heatmap(pivot_data, annot=True, fmt='.2f', cmap='YlOrRd')
    plt.title(f'{cell_type} - Methylation Patterns by Binding and Expression')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, f'{cell_type}_methylation_heatmap.pdf'))
    plt.close()
    
    # 2. Binding type proportions in different expression categories
    plt.figure(figsize=(10, 6))
    binding_props = pd.crosstab(
        df['expression_status'], 
        df['binding_type'],
        normalize='index'
    )
    binding_props.plot(kind='bar', stacked=True)
    plt.title(f'{cell_type} - Binding Type Distribution')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, f'{cell_type}_binding_proportions.pdf'))
    plt.close()
    
    # 3. Statistical comparison plots
    stats_dir = os.path.join(output_dir, 'binding_enrichment', 'statistics')
    os.makedirs(stats_dir, exist_ok=True)
    
    # Perform statistical tests
    binding_types = df['binding_type'].unique()
    expr_status = df['expression_status'].unique()
    
    stats_results = {
        'promoter': {},
        'gene_body': {}
    }
    
    for region in ['promoter', 'gene_body']:
        meth_col = f'{region}_methylation'
        for binding in binding_types:
            for status in expr_status:
                group_data = df[
                    (df['binding_type'] == binding) & 
                    (df['expression_status'] == status)
                ][meth_col]
                
                if len(group_data) > 0:
                    # Perform statistical tests
                    stats_results[region][f'{binding}_{status}'] = {
                        'mean': group_data.mean(),
                        'std': group_data.std(),
                        'n': len(group_data)
                    }
    
    # Save statistical results
    with open(os.path.join(stats_dir, f'{cell_type}_statistical_tests.txt'), 'w') as f:
        f.write(f"Statistical Analysis Results - {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        for region, results in stats_results.items():
            f.write(f"\n{region.upper()} METHYLATION\n")
            f.write("-"*30 + "\n")
            
            for group, stats in results.items():
                f.write(f"\n{group}:\n")
                f.write(f"  N = {stats['n']}\n")
                f.write(f"  Mean ± SD: {stats['mean']:.2f} ± {stats['std']:.2f}\n")

def create_exo_enriched_comparisons(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Create visualizations comparing exo-enriched binding patterns"""
    plot_dir = os.path.join(output_dir, 'exo_enriched_analysis', cell_type, 'plots')
    os.makedirs(plot_dir, exist_ok=True)
    
    # Check if required columns exist
    required_columns = ['binding_type', 'promoter_methylation', 'gene_body_methylation']
    if not all(col in df.columns for col in required_columns):
        logger.warning(f"Missing required columns for exo-enriched comparisons in {cell_type}")
        return
    
    # Filter for exo-enriched genes
    exo_enriched = df[df['binding_type'] == 'exo_only'].copy()
    
    if len(exo_enriched) == 0:
        logger.warning(f"No exo-enriched genes found for {cell_type}")
        return
    
    # Create methylation comparison plot
    plt.figure(figsize=(10, 6))
    
    # Create boxplot with individual points
    sns.boxplot(data=pd.melt(exo_enriched,
                            value_vars=['promoter_methylation', 'gene_body_methylation'],
                            var_name='Region',
                            value_name='Methylation Level'),
                x='Region',
                y='Methylation Level',
                showfliers=False)
    
    sns.stripplot(data=pd.melt(exo_enriched,
                              value_vars=['promoter_methylation', 'gene_body_methylation'],
                              var_name='Region',
                              value_name='Methylation Level'),
                  x='Region',
                  y='Methylation Level',
                  color='red',
                  alpha=0.3,
                  jitter=0.2,
                  size=4)
    
    plt.title(f'{cell_type} - Methylation Patterns in Exo-enriched Genes')
    plt.xlabel('Region')
    plt.ylabel('Methylation Level (%)')
    
    # Add statistical annotation
    stat_result = stats.wilcoxon(exo_enriched['promoter_methylation'],
                                exo_enriched['gene_body_methylation'])
    plt.text(0.5, plt.ylim()[1], f'p = {stat_result.pvalue:.2e}',
             horizontalalignment='center',
             verticalalignment='bottom')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'methylation_comparison.pdf'),
                bbox_inches='tight', dpi=300)
    plt.close()
    
    # Add additional visualizations as needed
    logger.info(f"Created exo-enriched comparisons for {cell_type}")

def analyze_mecp2_binding_patterns(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze MeCP2 binding patterns and their relationship with methylation and gene regulation"""
    
    # First, let's identify the binding regions
    # Assuming we have columns like 'location' or similar that indicates binding regions
    # Let's print the columns we have
    logger.info(f"Available columns: {df.columns.tolist()}")
    
    # 1. Analyze CpG island targeting patterns
    def analyze_cpg_island_targeting():
        # Identify differential binding patterns
        exo_only = df[
            (df['binding_type'] == 'exo') & 
            (df['location'] == 'promoter')  # Assuming 'location' is the correct column
        ]
        
        common_targets = df[
            (df['binding_type'] == 'both') & 
            (df['location'] == 'promoter')
        ]
        
        # Calculate binding strength ratio for common targets
        if 'exo_binding_strength' in df.columns and 'endo_binding_strength' in df.columns:
            common_targets['binding_ratio'] = common_targets['exo_binding_strength'] / \
                                            common_targets['endo_binding_strength']
        else:
            # If we don't have separate binding strength columns, use what we have
            logger.warning("Missing binding strength columns - using alternative metrics")
            # Add alternative calculation here based on available columns
        
        exo_enriched = common_targets[common_targets['binding_ratio'] > 1] \
            if 'binding_ratio' in common_targets.columns else pd.DataFrame()
        
        return {
            'exo_only': exo_only,
            'common_targets': common_targets,
            'exo_enriched': exo_enriched
        }
    
    # 2. Analyze methylation patterns in different genomic contexts
    def analyze_methylation_patterns():
        # Define genomic contexts based on available columns
        contexts = {
            'promoter': df[df['location'] == 'promoter'],
            'gene_body': df[df['location'] == 'gene_body'],
            'intergenic': df[df['location'] == 'intergenic']
        }
        
        methylation_stats = {}
        for context, context_df in contexts.items():
            methylation_stats[context] = {
                'methylation': {
                    'mean': context_df['methylation_level'].mean() if 'methylation_level' in df.columns else None,
                    'std': context_df['methylation_level'].std() if 'methylation_level' in df.columns else None
                }
            }
        
        return methylation_stats
    
    # 3. Analyze MeCP2 occupancy patterns
    def analyze_mecp2_occupancy():
        promoter_regions = df[df['location'] == 'promoter']
        
        occupancy_stats = {
            'endo_occupancy': len(promoter_regions[promoter_regions['binding_type'].isin(['endo', 'both'])]) / len(promoter_regions),
            'exo_occupancy': len(promoter_regions[promoter_regions['binding_type'].isin(['exo', 'both'])]) / len(promoter_regions),
            'binding_competition': {
                'endo_only': len(promoter_regions[promoter_regions['binding_type'] == 'endo']),
                'exo_only': len(promoter_regions[promoter_regions['binding_type'] == 'exo']),
                'both': len(promoter_regions[promoter_regions['binding_type'] == 'both'])
            }
        }
        
        return occupancy_stats
    
    # Run analyses and log results
    logger.info(f"Starting analysis for {cell_type}")
    results = {
        'cpg_targeting': analyze_cpg_island_targeting(),
        'methylation_patterns': analyze_methylation_patterns(),
        'occupancy': analyze_mecp2_occupancy()
    }
    
    # Create visualizations
    create_analysis_visualizations(results, cell_type, output_dir)
    
    # Save detailed results
    save_analysis_results(results, cell_type, output_dir)
    
    return results

def create_analysis_visualizations(results, cell_type, output_dir):
    """Create comprehensive visualizations for the analysis"""
    plot_dir = os.path.join(output_dir, 'mecp2_analysis', cell_type, 'plots')
    os.makedirs(plot_dir, exist_ok=True)
    
    # Set style
    plt.style.use('seaborn')
    sns.set_theme(style="whitegrid")
    sns.set_context("paper", font_scale=1.2)
    
    # 1. Methylation by binding type and expression status
    plt.figure(figsize=(12, 7))
    
    # Use violin plots for better distribution visualization
    sns.violinplot(
        data=results['methylation_data'],
        x='binding_type',
        y='gene_body_methylation',
        hue='expression_status',
        split=True,
        inner='quartile',
        cut=0
    )
    
    plt.title(f'{cell_type} - Gene Body Methylation by Binding Type')
    plt.xlabel('Binding Type')
    plt.ylabel('Gene Body Methylation (%)')
    plt.xticks(rotation=0)
    
    # Adjust legend
    plt.legend(title='Expression Status', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'methylation_by_binding.pdf'),
                bbox_inches='tight', dpi=300)
    plt.close()
    
    # 2. Binding patterns - use horizontal bars for better label readability
    plt.figure(figsize=(10, 6))
    occupancy = results['occupancy']['binding_competition']
    
    # Sort values for better visualization
    categories = ['Both', 'Exo Only', 'Endo Only']
    values = [occupancy['both'], occupancy['exo_only'], occupancy['endo_only']]
    
    # Create horizontal bar plot
    bars = plt.barh(categories, values)
    
    # Add value labels on the bars
    for bar in bars:
        width = bar.get_width()
        plt.text(width, bar.get_y() + bar.get_height()/2,
                f'{int(width):,}',
                ha='left', va='center', fontsize=10)
    
    plt.title(f'{cell_type} - MeCP2 Binding Patterns')
    plt.xlabel('Number of Regions')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'binding_patterns.pdf'))
    plt.close()
    
    # 3. Methylation distribution comparison
    if 'methylation_comparison' in results:
        plt.figure(figsize=(10, 6))
        
        # Use kernel density estimation for smooth distribution visualization
        for binding_type in ['endo_only', 'exo_only', 'both']:
            data = results['methylation_comparison'][binding_type]
            sns.kdeplot(
                data=data,
                label=binding_type.replace('_', ' ').title(),
                common_norm=False
            )
        
        plt.title(f'{cell_type} - Methylation Distribution by Binding Type')
        plt.xlabel('Methylation Level (%)')
        plt.ylabel('Density')
        plt.legend(title='Binding Type')
        
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'methylation_distribution.pdf'))
        plt.close()

def save_analysis_results(results, cell_type, output_dir):
    """Save detailed analysis results"""
    results_dir = os.path.join(output_dir, 'mecp2_analysis', cell_type, 'results')
    os.makedirs(results_dir, exist_ok=True)
    
    with open(os.path.join(results_dir, 'analysis_summary.txt'), 'w') as f:
        f.write(f"MeCP2 Binding Analysis - {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        # CpG island targeting
        f.write("CpG ISLAND TARGETING\n")
        f.write("-"*30 + "\n")
        cpg_results = results['cpg_targeting']
        f.write(f"Exo-only targets: {len(cpg_results['exo_only'])}\n")
        f.write(f"Common targets: {len(cpg_results['common_targets'])}\n")
        f.write(f"Exo-enriched targets: {len(cpg_results['exo_enriched'])}\n\n")
        
        # Occupancy patterns
        f.write("OCCUPANCY PATTERNS\n")
        f.write("-"*30 + "\n")
        occ = results['occupancy']
        f.write(f"Endo occupancy: {occ['endo_occupancy']:.2%}\n")
        f.write(f"Exo occupancy: {occ['exo_occupancy']:.2%}\n\n")
        
        # Save detailed gene lists and additional statistics...

def create_tss_analysis_plots(results: dict, cell_type: str, output_dir: str):
    """Create visualizations for TSS binding analysis"""
    
    # 1. Regulation distribution plot
    plt.figure(figsize=(12, 6))
    
    # Prepare data for plotting
    categories = []
    not_dereg = []
    up = []
    down = []
    
    for category, stats in results.items():
        categories.append(category)
        not_dereg.append(stats.get('not_deregulated', {}).get('count', 0))
        up.append(stats.get('upregulated', {}).get('count', 0))
        down.append(stats.get('downregulated', {}).get('count', 0))
    
    # Create stacked bar plot
    width = 0.35
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.bar(categories, not_dereg, width, label='Not Deregulated')
    ax.bar(categories, up, width, bottom=not_dereg, label='Upregulated')
    ax.bar(categories, down, width, bottom=[i+j for i,j in zip(not_dereg, up)], 
           label='Downregulated')
    
    plt.title(f'{cell_type} - Gene Regulation by Binding Category')
    plt.xlabel('Binding Category')
    plt.ylabel('Number of Genes')
    plt.legend()
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'regulation_distribution.pdf'))
    plt.close()
    
    # 2. Methylation comparison plot
    create_methylation_comparison_plot(results, cell_type, output_dir)

def create_methylation_comparison_plot(results: dict, cell_type: str, output_dir: str):
    """Create methylation comparison plots for different binding categories with statistics"""
    
    # Prepare data for plotting
    plot_data = {
        'Category': [],
        'Regulation': [],
        'Region': [],
        'Methylation': [],
        'StdDev': []
    }
    
    for category, stats in results.items():
        for reg_status, reg_stats in stats.items():
            # Promoter methylation
            plot_data['Category'].append(category)
            plot_data['Regulation'].append(reg_status)
            plot_data['Region'].append('Promoter')
            plot_data['Methylation'].append(
                reg_stats['promoter_methylation']['mean']
            )
            plot_data['StdDev'].append(
                reg_stats['promoter_methylation']['std']
            )
            
            # Gene body methylation
            plot_data['Category'].append(category)
            plot_data['Regulation'].append(reg_status)
            plot_data['Region'].append('Gene Body')
            plot_data['Methylation'].append(
                reg_stats['gene_body_methylation']['mean']
            )
            plot_data['StdDev'].append(
                reg_stats['gene_body_methylation']['std']
            )
    
    df = pd.DataFrame(plot_data)
    
    # Create separate plots for Promoter and Gene Body
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # Promoter methylation plot
    promoter_data = df[df['Region'] == 'Promoter']
    sns.barplot(data=promoter_data, x='Category', y='Methylation',
                hue='Regulation', ax=ax1)
    ax1.set_title(f'{cell_type} - Promoter Methylation')
    ax1.set_xlabel('Binding Category')
    ax1.set_ylabel('Methylation Level (%)')
    ax1.tick_params(axis='x', rotation=45)
    
    # Gene body methylation plot
    gene_body_data = df[df['Region'] == 'Gene Body']
    sns.barplot(data=gene_body_data, x='Category', y='Methylation',
                hue='Regulation', ax=ax2)
    ax2.set_title(f'{cell_type} - Gene Body Methylation')
    ax2.set_xlabel('Binding Category')
    ax2.set_ylabel('Methylation Level (%)')
    ax2.tick_params(axis='x', rotation=45)
    
    # Add statistical significance
    for ax, data in [(ax1, promoter_data), (ax2, gene_body_data)]:
        for i, category in enumerate(data['Category'].unique()):
            category_data = data[data['Category'] == category]
            
            # Compare upregulated vs downregulated
            up_data = category_data[category_data['Regulation'] == 'upregulated']
            down_data = category_data[category_data['Regulation'] == 'downregulated']
            
            if len(up_data) > 0 and len(down_data) > 0:
                stats = calculate_group_significance(
                    up_data['Methylation'],
                    down_data['Methylation'],
                    'upregulated',
                    'downregulated'
                )
                
                if stats and stats['pvalue'] < 0.05:
                    y_max = max(up_data['Methylation'].max(), down_data['Methylation'].max())
                    ax.text(i, y_max + 2,
                           f"p={stats['pvalue']:.2e}",
                           ha='center')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'methylation_comparison.pdf'))
    plt.close()

def calculate_group_significance(group1: pd.Series, group2: pd.Series, 
                               label1: str, label2: str) -> Dict:
    """Calculate statistical significance between two groups with improved handling of small samples"""
    try:
        # Remove NaN values
        group1_clean = group1.dropna()
        group2_clean = group2.dropna()
        
        # Check if we have enough samples (at least 3 per group for statistical validity)
        if len(group1_clean) < 3 or len(group2_clean) < 3:
            logger.debug(f"Insufficient samples for statistical test: {label1}={len(group1_clean)}, {label2}={len(group2_clean)}")
            return None
            
        # Perform Mann-Whitney U test
        stat, pval = stats.mannwhitneyu(group1_clean, group2_clean, alternative='two-sided')
        
        # Calculate effect size only if we have enough samples
        if len(group1_clean) >= 3 and len(group2_clean) >= 3:
            d = (group1_clean.mean() - group2_clean.mean()) / np.sqrt(
                ((group1_clean.std() ** 2 + group2_clean.std() ** 2) / 2)
            )
        else:
            d = None
        
        return {
            'test': 'Mann-Whitney U',
            'statistic': stat,
            'pvalue': pval,
            'cohens_d': d,
            'group1_mean': group1_clean.mean(),
            'group2_mean': group2_clean.mean(),
            'group1_n': len(group1_clean),
            'group2_n': len(group2_clean)
        }
    except Exception as e:
        logger.debug(f"Error calculating significance between {label1} and {label2}: {str(e)}")
        return None

def summarize_binding_analysis(results: Dict[str, pd.DataFrame], output_dir: str):
    """Create a summary of the binding analysis results"""
    summary_file = os.path.join(output_dir, 'binding_analysis_summary.txt')
    
    with open(summary_file, 'w') as f:
        f.write("MeCP2 Binding Analysis Summary\n")
        f.write("="*50 + "\n\n")
        
        for cell_type, df in results.items():
            f.write(f"\n{cell_type} Analysis\n")
            f.write("-"*30 + "\n")
            
            # Overall statistics
            f.write(f"Total genes analyzed: {len(df)}\n")
            f.write(f"MeCP2-bound genes: {df['mecp2_bound'].sum()}\n")
            
            # Binding type distribution
            if 'binding_type' in df.columns:
                f.write("\nBinding type distribution:\n")
                f.write(df['binding_type'].value_counts().to_string())
                
            # Expression status distribution
            f.write("\n\nExpression status distribution:\n")
            f.write(df['expression_status'].value_counts().to_string())
            
            # Methylation statistics
            f.write("\n\nMethylation statistics:\n")
            for region in ['promoter', 'gene_body']:
                f.write(f"\n{region.title()} methylation:\n")
                f.write(f"Mean: {df[f'{region}_methylation'].mean():.2f}\n")
                f.write(f"Std: {df[f'{region}_methylation'].std():.2f}\n")
            
            f.write("\n" + "="*50 + "\n")

def identify_cpg_islands(sequence: str, min_length: int = 200, min_gc: float = 0.5, min_oe: float = 0.6) -> List[Tuple[int, int]]:
    """
    Identify CpG islands using standard criteria:
    - Length > 200bp
    - GC content > 50%
    - Observed/Expected CpG ratio > 0.6
    
    Args:
        sequence: DNA sequence to analyze
        min_length: Minimum length for CpG island (default 200bp)
        min_gc: Minimum GC content (default 50%)
        min_oe: Minimum observed/expected CpG ratio (default 0.6)
    
    Returns:
        List of tuples containing (start, end) positions of CpG islands
    """
    cpg_islands = []
    sequence = sequence.upper()
    
    # Sliding window approach
    window_size = 200  # Start with minimum CpG island size
    
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        
        # Calculate GC content
        gc_content = (window.count('G') + window.count('C')) / len(window)
        
        # Calculate observed/expected CpG ratio
        cpg_count = window.count('CG')
        c_count = window.count('C')
        g_count = window.count('G')
        expected_cpg = (c_count * g_count) / len(window)
        if expected_cpg > 0:
            oe_ratio = cpg_count / expected_cpg
        else:
            oe_ratio = 0
            
        # Check if window meets CpG island criteria
        if (gc_content >= min_gc and 
            oe_ratio >= min_oe and 
            cpg_count >= 4):  # Minimum 4 CpGs
            
            # Extend window to find full CpG island
            start = i
            end = i + window_size
            
            # Extend forward
            while (end < len(sequence) and
                   is_cpg_island_region(sequence[start:end+1], min_gc, min_oe)):
                end += 1
                
            # Extend backward
            while (start > 0 and
                   is_cpg_island_region(sequence[start-1:end], min_gc, min_oe)):
                start -= 1
            
            # Add if meets minimum length requirement
            if end - start >= min_length:
                cpg_islands.append((start, end))
                i = end  # Skip to end of found CpG island
                
    # Merge overlapping islands
    if cpg_islands:
        cpg_islands.sort()
        merged = []
        current_start, current_end = cpg_islands[0]
        
        for start, end in cpg_islands[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                merged.append((current_start, current_end))
                current_start, current_end = start, end
        merged.append((current_start, current_end))
        
        return merged
    
    return cpg_islands

def is_cpg_island_region(sequence: str, min_gc: float, min_oe: float) -> bool:
    """Helper function to check if a sequence meets CpG island criteria"""
    sequence = sequence.upper()
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
    
    cpg_count = sequence.count('CG')
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    expected_cpg = (c_count * g_count) / len(sequence)
    
    if expected_cpg > 0:
        oe_ratio = cpg_count / expected_cpg
    else:
        oe_ratio = 0
        
    return gc_content >= min_gc and oe_ratio >= min_oe

def calculate_cpg_island_methylation(sequence: str, ip_values: List[float], 
                                   input_values: List[float]) -> Dict[str, Any]:
    """Calculate methylation specifically for CpG islands"""
    try:
        # Identify CpG islands
        cpg_islands = identify_cpg_islands(sequence)
        
        # Calculate methylation for each island
        methylation_levels = []
        for start, end in cpg_islands:
            # Convert to numpy arrays for memory efficiency
            island_ip = np.array(ip_values[start:end], dtype=np.float32)
            island_input = np.array(input_values[start:end], dtype=np.float32)
            island_seq = sequence[start:end]
            
            methylation = calculate_normalized_methylation(
                island_ip,
                island_input,
                island_seq
            )
            methylation_levels.append(methylation)
            
            # Clear memory
            del island_ip
            del island_input
        
        # Calculate statistics
        if methylation_levels:
            stats = {
                'mean_methylation': float(np.mean(methylation_levels)),
                'std_methylation': float(np.std(methylation_levels)),
                'n_islands': len(cpg_islands),
                'total_length': sum(end - start for start, end in cpg_islands)
            }
        else:
            stats = {
                'mean_methylation': 0.0,
                'std_methylation': 0.0,
                'n_islands': 0,
                'total_length': 0
            }
        
        return {
            'cpg_islands': cpg_islands,
            'methylation_levels': methylation_levels,
            'statistics': stats
        }
        
    except Exception as e:
        logger.error(f"Error in CpG island methylation calculation: {str(e)}")
        return {
            'cpg_islands': [],
            'methylation_levels': [],
            'statistics': {
                'mean_methylation': 0.0,
                'std_methylation': 0.0,
                'n_islands': 0,
                'total_length': 0
            }
        }

def analyze_cpg_island_patterns(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze CpG island methylation patterns"""
    
    logger.info(f"Starting CpG island analysis for {cell_type}")
    
    # Create output directory
    cpg_dir = os.path.join(output_dir, 'cpg_analysis', cell_type)
    os.makedirs(cpg_dir, exist_ok=True)
    logger.info(f"Created output directory: {cpg_dir}")
    
    # Check if we have CpG island data
    if not all(col in df.columns for col in ['cpg_island_count', 'cpg_island_mean_methylation']):
        logger.error(f"Missing CpG island data columns for {cell_type}")
        return None
    
    logger.info(f"Found {df['cpg_island_count'].sum()} CpG islands in {len(df)} regions")
    
    try:
        # 1. Compare CpG island vs non-island methylation
        plt.figure(figsize=(10, 6))
        data_to_plot = pd.DataFrame({
            'CpG Islands': df['cpg_island_mean_methylation'],
            'Overall': df['methylation']
        })
        sns.boxplot(data=pd.melt(data_to_plot), x='variable', y='value')
        plt.title(f'{cell_type}: CpG Island vs Overall Methylation')
        plt.ylabel('Methylation Level (%)')
        plt.savefig(os.path.join(cpg_dir, 'cpg_vs_overall_methylation.pdf'))
        plt.close()
        logger.info("Created methylation comparison plot")
        
        # 2. CpG island count distribution
        plt.figure(figsize=(10, 6))
        sns.histplot(data=df, x='cpg_island_count', bins=30)
        plt.title(f'{cell_type}: CpG Island Count Distribution')
        plt.xlabel('Number of CpG Islands')
        plt.savefig(os.path.join(cpg_dir, 'cpg_island_count_distribution.pdf'))
        plt.close()
        logger.info("Created CpG island count distribution plot")
        
        # 3. Save detailed statistics
        stats_file = os.path.join(cpg_dir, 'cpg_analysis_stats.txt')
        with open(stats_file, 'w') as f:
            f.write(f"CpG Island Analysis - {cell_type}\n")
            f.write("="*50 + "\n\n")
            
            # Basic statistics
            f.write("Overall Statistics:\n")
            f.write("-"*20 + "\n")
            total_islands = df['cpg_island_count'].sum()
            regions_with_islands = (df['cpg_island_count'] > 0).sum()
            f.write(f"Total CpG islands: {total_islands}\n")
            f.write(f"Regions with CpG islands: {regions_with_islands} ({regions_with_islands/len(df)*100:.1f}%)\n")
            f.write(f"Average CpG islands per region: {df['cpg_island_count'].mean():.2f}\n\n")
            
            # Methylation statistics
            f.write("Methylation Statistics:\n")
            f.write("-"*20 + "\n")
            f.write("CpG Island Methylation:\n")
            f.write(f"Mean: {df['cpg_island_mean_methylation'].mean():.2f}%\n")
            f.write(f"Median: {df['cpg_island_mean_methylation'].median():.2f}%\n")
            f.write(f"Std Dev: {df['cpg_island_mean_methylation'].std():.2f}%\n\n")
            
            # By binding type if available
            if 'binding_type' in df.columns:
                f.write("CpG Islands by Binding Type:\n")
                f.write("-"*20 + "\n")
                for binding_type in df['binding_type'].unique():
                    subset = df[df['binding_type'] == binding_type]
                    if len(subset) > 0:
                        f.write(f"\n{binding_type}:\n")
                        f.write(f"Count: {subset['cpg_island_count'].sum()}\n")
                        f.write(f"Mean methylation: {subset['cpg_island_mean_methylation'].mean():.2f}%\n")
        
        logger.info(f"Saved CpG island statistics to {stats_file}")
        
        return {
            'total_islands': total_islands,
            'regions_with_islands': regions_with_islands,
            'mean_methylation': df['cpg_island_mean_methylation'].mean(),
            'std_methylation': df['cpg_island_mean_methylation'].std()
        }
        
    except Exception as e:
        logger.error(f"Error in CpG island analysis for {cell_type}: {str(e)}")
        return None

def analyze_cpg_islands_detailed(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Detailed analysis of CpG islands and their relationship with MeCP2 binding and expression"""
    
    logger.info(f"Starting detailed CpG analysis for {cell_type}")
    logger.info(f"DataFrame shape: {df.shape}")
    logger.info(f"Available columns: {df.columns.tolist()}")
    
    # Check if we have CpG data
    required_columns = ['cpg_island_count', 'cpg_island_mean_methylation', 'binding_type']
    if not all(col in df.columns for col in required_columns):
        missing = [col for col in required_columns if col not in df.columns]
        logger.error(f"Missing required columns for CpG analysis: {missing}")
        logger.error(f"Data preview:\n{df.head()}")
        return None
    
    # Verify we have valid data
    if df['cpg_island_count'].sum() == 0:
        logger.error("No CpG islands found in any regions")
        return None
    
    # Create specific output directory for CpG analysis
    cpg_dir = os.path.join(output_dir, cell_type)
    os.makedirs(cpg_dir, exist_ok=True)
    logger.info(f"Created CpG analysis directory: {cpg_dir}")
    
    # 1. CpG Islands vs MeCP2 Binding
    binding_dir = os.path.join(cpg_dir, 'binding_analysis')
    os.makedirs(binding_dir, exist_ok=True)
    
    # Plot CpG methylation by binding type
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=df, x='binding_type', y='cpg_island_mean_methylation')
    plt.title(f'{cell_type}: CpG Island Methylation by MeCP2 Binding Type')
    plt.xlabel('Binding Type')
    plt.ylabel('CpG Island Methylation (%)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(binding_dir, 'cpg_methylation_by_binding.pdf'))
    plt.close()
    
    # 2. CpG Islands vs Expression
    expression_dir = os.path.join(cpg_dir, 'expression_analysis')
    os.makedirs(expression_dir, exist_ok=True)
    
    if 'log2FoldChange' in df.columns:
        # Scatter plot of CpG methylation vs expression
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=df, x='cpg_island_mean_methylation', y='log2FoldChange', 
                       hue='binding_type', alpha=0.6)
        plt.title(f'{cell_type}: CpG Methylation vs Expression')
        plt.xlabel('CpG Island Methylation (%)')
        plt.ylabel('log2 Fold Change')
        plt.tight_layout()
        plt.savefig(os.path.join(expression_dir, 'cpg_methylation_vs_expression.pdf'))
        plt.close()
    
    # 3. Save detailed statistics
    stats_file = os.path.join(cpg_dir, 'detailed_cpg_analysis.txt')
    with open(stats_file, 'w') as f:
        f.write(f"Detailed CpG Island Analysis - {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        # Overall statistics
        f.write("1. Overall Statistics\n")
        f.write("-"*20 + "\n")
        total_islands = df['cpg_island_count'].sum()
        regions_with_islands = (df['cpg_island_count'] > 0).sum()
        f.write(f"Total CpG islands: {total_islands}\n")
        f.write(f"Regions with islands: {regions_with_islands} ({regions_with_islands/len(df)*100:.1f}%)\n\n")
        
        # Binding analysis
        f.write("2. CpG Islands by Binding Type\n")
        f.write("-"*20 + "\n")
        for binding_type in df['binding_type'].unique():
            subset = df[df['binding_type'] == binding_type]
            if len(subset) > 0:
                f.write(f"\n{binding_type}:\n")
                f.write(f"Number of regions: {len(subset)}\n")
                f.write(f"Total CpG islands: {subset['cpg_island_count'].sum()}\n")
                f.write(f"Mean methylation: {subset['cpg_island_mean_methylation'].mean():.2f}%\n")
                f.write(f"Median methylation: {subset['cpg_island_mean_methylation'].median():.2f}%\n")
                
                # Expression analysis if available
                if 'log2FoldChange' in subset.columns:
                    corr = calculate_correlation(
                        subset['cpg_island_mean_methylation'],
                        subset['log2FoldChange']
                    )
                    if corr:
                        f.write(f"Correlation with expression: rho={corr['statistic']:.3f}, p={corr['pvalue']:.2e}\n")
        
        # Save summary table
        summary_df = pd.DataFrame({
            'binding_type': df['binding_type'].unique(),
            'n_regions': [len(df[df['binding_type'] == bt]) for bt in df['binding_type'].unique()],
            'total_cpg_islands': [df[df['binding_type'] == bt]['cpg_island_count'].sum() 
                                for bt in df['binding_type'].unique()],
            'mean_methylation': [df[df['binding_type'] == bt]['cpg_island_mean_methylation'].mean() 
                               for bt in df['binding_type'].unique()]
        })
        summary_file = os.path.join(cpg_dir, 'cpg_summary.csv')
        summary_df.to_csv(summary_file, index=False)
    
    return {
        'total_islands': total_islands,
        'regions_with_islands': regions_with_islands,
        'binding_stats': summary_df.to_dict('records')
    }

def load_cpg_islands(bed_file: str) -> pd.DataFrame:
    """Load CpG islands from BED file"""
    logger.info(f"Loading CpG islands from {bed_file}")
    try:
        # Read with correct column names for PyRanges
        cpg_islands = pd.read_csv(bed_file, sep='\t', 
                                 names=['Chromosome', 'Start', 'End', 'id', 'type', 'cpg_count'])
        
        # Ensure chromosome column format matches your data
        if not cpg_islands['Chromosome'].str.contains('chr').all():
            cpg_islands['Chromosome'] = 'chr' + cpg_islands['Chromosome'].astype(str)
        
        logger.info(f"Loaded {len(cpg_islands)} CpG islands")
        logger.info(f"Columns: {cpg_islands.columns.tolist()}")
        logger.info(f"Sample data:\n{cpg_islands.head()}")
        
        return cpg_islands
    except Exception as e:
        logger.error(f"Error loading CpG islands: {str(e)}")
        return pd.DataFrame()

def create_cpg_methylation_plots(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Create visualizations for CpG island methylation analysis"""
    
    try:
        # 1. CpG methylation distribution
        plt.figure(figsize=(10, 6))
        sns.histplot(data=df, x='cpg_island_mean_methylation', bins=50)
        plt.title(f'{cell_type}: CpG Island Methylation Distribution')
        plt.xlabel('Methylation Level (%)')
        plt.ylabel('Count')
        plt.savefig(os.path.join(output_dir, 'cpg_methylation_distribution.pdf'))
        plt.close()
        
        # 2. CpG count vs methylation
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=df, x='cpg_island_count', y='cpg_island_mean_methylation')
        plt.title(f'{cell_type}: CpG Count vs Methylation')
        plt.xlabel('Number of CpG Islands')
        plt.ylabel('Mean Methylation Level (%)')
        plt.savefig(os.path.join(output_dir, 'cpg_count_vs_methylation.pdf'))
        plt.close()
        
        # 3. CpG methylation by binding type
        if 'binding_type' in df.columns:
            plt.figure(figsize=(12, 6))
            sns.boxplot(data=df, x='binding_type', y='cpg_island_mean_methylation')
            plt.title(f'{cell_type}: CpG Methylation by Binding Type')
            plt.xlabel('Binding Type')
            plt.ylabel('Methylation Level (%)')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'cpg_methylation_by_binding.pdf'))
            plt.close()
            
    except Exception as e:
        logger.error(f"Error creating CpG methylation plots: {str(e)}")

def generate_cpg_report(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Generate detailed report of CpG island analysis"""
    
    report_file = os.path.join(output_dir, f'{cell_type}_cpg_analysis_report.txt')
    logger.info(f"Generating CpG analysis report: {report_file}")
    
    try:
        with open(report_file, 'w') as f:
            f.write(f"CpG Island Analysis Report - {cell_type}\n")
            f.write("="*50 + "\n\n")
            
            # Basic statistics
            f.write("1. Overall Statistics\n")
            f.write("-"*20 + "\n")
            total_regions = len(df)
            regions_with_cpg = (df['cpg_island_count'] > 0).sum()
            total_cpgs = df['cpg_island_count'].sum()
            
            f.write(f"Total regions analyzed: {total_regions}\n")
            f.write(f"Regions with CpG islands: {regions_with_cpg} ({regions_with_cpg/total_regions*100:.1f}%)\n")
            f.write(f"Total CpG islands found: {total_cpgs}\n")
            f.write(f"Average CpG islands per region: {total_cpgs/total_regions:.2f}\n\n")
            
            # Methylation statistics
            f.write("2. Methylation Statistics\n")
            f.write("-"*20 + "\n")
            f.write("CpG Island Methylation:\n")
            f.write(f"Mean: {df['cpg_island_mean_methylation'].mean():.2f}%\n")
            f.write(f"Median: {df['cpg_island_mean_methylation'].median():.2f}%\n")
            f.write(f"Std Dev: {df['cpg_island_mean_methylation'].std():.2f}%\n\n")
            
            # Save summary to CSV
            summary_file = os.path.join(output_dir, f'{cell_type}_cpg_summary.csv')
            summary_df = pd.DataFrame({
                'metric': ['total_regions', 'regions_with_cpg', 'total_cpgs', 'mean_methylation'],
                'value': [total_regions, regions_with_cpg, total_cpgs, 
                         df['cpg_island_mean_methylation'].mean()]
            })
            summary_df.to_csv(summary_file, index=False)
            
        logger.info(f"Generated CpG analysis report and summary")
        return True
        
    except Exception as e:
        logger.error(f"Error generating CpG report: {str(e)}")
        return False

def analyze_methylation_patterns(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze methylation patterns in relation to MeCP2 binding and gene expression"""
    
    logger.info(f"\nAnalyzing methylation patterns for {cell_type}...")
    
    # Create output directory for detailed analysis
    analysis_dir = os.path.join(output_dir, 'methylation_analysis')
    os.makedirs(analysis_dir, exist_ok=True)
    
    # 1. Calculate summary statistics
    binding_stats = df.groupby('binding_type').agg({
        'cpg_island_count': ['mean', 'std', 'count'],
        'cpg_island_mean_methylation': ['mean', 'std', 'median']
    }).round(2)
    
    # Save statistics
    stats_file = os.path.join(analysis_dir, f'{cell_type}_methylation_stats.txt')
    with open(stats_file, 'w') as f:
        f.write(f"Methylation Analysis Report - {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        # Overall statistics
        f.write("1. Overall Statistics\n")
        f.write("-"*20 + "\n")
        f.write(f"Total regions analyzed: {len(df)}\n")
        f.write(f"Regions with CpG islands: {(df['cpg_island_count'] > 0).sum()}\n")
        f.write(f"Average methylation: {df['cpg_island_mean_methylation'].mean():.2f}%\n\n")
        
        # Binding type statistics
        f.write("2. Statistics by Binding Type\n")
        f.write("-"*20 + "\n")
        f.write(binding_stats.to_string())
        f.write("\n\n")
        
        # Expression correlation
        if 'expression_status' in df.columns:
            expr_stats = df.groupby('expression_status')['cpg_island_mean_methylation'].agg(['mean', 'std', 'count'])
            f.write("3. Methylation by Expression Status\n")
            f.write("-"*20 + "\n")
            f.write(expr_stats.to_string())
    
    # 2. Create visualizations
    # Methylation distribution by binding type
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df[df['cpg_island_count'] > 0], 
                x='binding_type', 
                y='cpg_island_mean_methylation')
    plt.title(f'{cell_type}: Methylation Levels by MeCP2 Binding')
    plt.xlabel('Binding Type')
    plt.ylabel('Methylation Level (%)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_dir, f'{cell_type}_methylation_by_binding.pdf'))
    plt.close()
    
    # Methylation vs Expression
    if 'expression_status' in df.columns:
        plt.figure(figsize=(12, 6))
        sns.boxplot(data=df[df['cpg_island_count'] > 0],
                   x='expression_status',
                   y='cpg_island_mean_methylation',
                   hue='binding_type')
        plt.title(f'{cell_type}: Methylation by Expression and Binding')
        plt.xlabel('Expression Status')
        plt.ylabel('Methylation Level (%)')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(analysis_dir, f'{cell_type}_methylation_expression.pdf'))
        plt.close()
    
    logger.info(f"Analysis completed. Results saved to {analysis_dir}")
    return binding_stats

def analyze_binding_methylation_expression(df: pd.DataFrame, cell_type: str, output_dir: str):
    """
    Analyze the relationship between MeCP2 binding, methylation, and gene expression
    with focus on Exo-enriched and control groups
    """
    logger.info(f"\nAnalyzing binding-methylation-expression patterns for {cell_type}...")
    
    # Create output directory
    analysis_dir = os.path.join(output_dir, 'binding_methylation_analysis')
    os.makedirs(analysis_dir, exist_ok=True)
    
    # Define binding categories of interest
    binding_categories = {
        'Exo-enriched': df[df['binding_type'].isin(['exo', 'both'])],
        'Exo-only': df[df['binding_type'] == 'exo'],
        'Non-enriched': df[~df['mecp2_bound']],
        'Endo-only': df[df['binding_type'] == 'endo']
    }
    
    # Initialize results dictionary
    results = {}
    
    # Generate summary report
    report_file = os.path.join(analysis_dir, f'{cell_type}_binding_analysis_report.txt')
    with open(report_file, 'w') as f:
        f.write(f"MeCP2 Binding Analysis Report - {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        for category, category_df in binding_categories.items():
            if len(category_df) == 0:
                logger.warning(f"No genes found in category: {category}")
                continue
                
            # Get expression groups
            expr_groups = {
                'not_deregulated': category_df[category_df['expression_status'] == 'unchanged'],
                'upregulated': category_df[category_df['expression_status'] == 'upregulated'],
                'downregulated': category_df[category_df['expression_status'] == 'downregulated']
            }
            
            # Save gene lists
            for status, group_df in expr_groups.items():
                output_file = os.path.join(analysis_dir, f'{category}_{status}_genes.csv')
                if len(group_df) > 0:
                    group_df.to_csv(output_file, index=False)
            
            # Write statistics to report
            f.write(f"\n{category} Analysis\n")
            f.write("-"*30 + "\n")
            f.write(f"Total genes: {len(category_df)}\n")
            
            for status, group_df in expr_groups.items():
                f.write(f"\n{status.title()}:\n")
                f.write(f"  Count: {len(group_df)} ({len(group_df)/len(category_df)*100:.1f}%)\n")
                if len(group_df) > 0:
                    f.write(f"  Mean promoter methylation: {group_df['cpg_island_mean_methylation'].mean():.2f}%\n")
                    
            results[category] = expr_groups
    
    # Create visualizations
    
    # 1. Expression status distribution by binding category
    plt.figure(figsize=(12, 6))
    data = []
    categories = []
    statuses = []
    for category, expr_groups in results.items():
        for status, group_df in expr_groups.items():
            data.append(len(group_df))
            categories.append(category)
            statuses.append(status)
    
    plot_df = pd.DataFrame({
        'Category': categories,
        'Status': statuses,
        'Count': data
    })
    
    sns.barplot(data=plot_df, x='Category', y='Count', hue='Status')
    plt.title(f'{cell_type}: Gene Expression Status by Binding Category')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_dir, f'{cell_type}_expression_by_binding.pdf'))
    plt.close()
    
    # 2. Methylation levels by binding category and expression status
    plt.figure(figsize=(15, 6))
    methylation_data = []
    
    for category, expr_groups in results.items():
        for status, group_df in expr_groups.items():
            if len(group_df) > 0:
                methylation_data.append(pd.DataFrame({
                    'Category': category,
                    'Status': status,
                    'Methylation': group_df['cpg_island_mean_methylation'],
                }))
    
    if methylation_data:
        methylation_df = pd.concat(methylation_data)
        sns.boxplot(data=methylation_df, x='Category', y='Methylation', hue='Status')
        plt.title(f'{cell_type}: Methylation Levels by Binding and Expression')
        plt.xlabel('Binding Category')
        plt.ylabel('Methylation Level (%)')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(analysis_dir, f'{cell_type}_methylation_patterns.pdf'))
        plt.close()
    
    # 3. Methylation comparison between promoter regions
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=df[df['cpg_island_count'] > 0], 
                x='binding_type', 
                y='cpg_island_mean_methylation',
                hue='expression_status')
    plt.title(f'{cell_type}: Promoter Methylation by Binding and Expression')
    plt.xlabel('Binding Type')
    plt.ylabel('Methylation Level (%)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_dir, f'{cell_type}_promoter_methylation.pdf'))
    plt.close()
    
    logger.info(f"Analysis completed. Results saved to {analysis_dir}")
    return results

def analyze_mecp2_detailed_binding(df: pd.DataFrame, cell_type: str, output_dir: str):
    """
    Detailed analysis of MeCP2 binding patterns focusing on:
    1. Exo vs Endo binding comparison
    2. CpG island enrichment
    3. Methylation status of bound regions
    4. Cell-type specific occupancy patterns
    """
    # Create output directory
    analysis_dir = os.path.join(output_dir, 'detailed_mecp2_analysis')
    os.makedirs(analysis_dir, exist_ok=True)
    
    # 1. Analyze Exo vs Endo binding
    exo_only = df[df['binding_type'] == 'exo']
    endo_only = df[df['binding_type'] == 'endo']
    common_targets = df[df['binding_type'] == 'both']
    
    # Compare Exo vs Endo signal in common targets
    common_targets['exo_endo_ratio'] = common_targets['exo_signal'] / common_targets['endo_signal']
    exo_enriched = common_targets[common_targets['exo_endo_ratio'] > 1]
    
    # Save results
    binding_report = os.path.join(analysis_dir, f'{cell_type}_binding_comparison.txt')
    with open(binding_report, 'w') as f:
        f.write(f"MeCP2 Binding Pattern Analysis - {cell_type}\n")
        f.write("="*50 + "\n\n")
        f.write(f"Exo-only targets: {len(exo_only)}\n")
        f.write(f"Endo-only targets: {len(endo_only)}\n")
        f.write(f"Common targets: {len(common_targets)}\n")
        f.write(f"Targets with higher Exo signal: {len(exo_enriched)}\n")
    
    # 2. Analyze gene body methylation for downregulated genes
    downreg_genes = df[df['expression_status'] == 'downregulated']
    
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=downreg_genes, x='binding_type', 
                y='cpg_island_mean_methylation', hue='expression_status')
    plt.title(f'{cell_type}: Gene Body Methylation in Downregulated Genes')
    plt.savefig(os.path.join(analysis_dir, f'{cell_type}_downreg_methylation.pdf'))
    plt.close()
    
    return {
        'exo_only': exo_only,
        'endo_only': endo_only,
        'common_targets': common_targets,
        'exo_enriched': exo_enriched
    }

def analyze_cpg_occupancy(df: pd.DataFrame, cell_type: str, output_dir: str):
    """
    Analyze CpG island occupancy patterns and their relationship with gene regulation
    """
    analysis_dir = os.path.join(output_dir, 'cpg_occupancy_analysis')
    os.makedirs(analysis_dir, exist_ok=True)
    
    # Calculate occupancy metrics
    cpg_targets = df[df['cpg_island_count'] > 0]
    total_cpgs = len(cpg_targets)
    
    occupancy_stats = {
        'endo': len(cpg_targets[cpg_targets['endo_signal'] > 0]) / total_cpgs * 100,
        'exo': len(cpg_targets[cpg_targets['exo_signal'] > 0]) / total_cpgs * 100,
        'both': len(cpg_targets[cpg_targets['binding_type'] == 'both']) / total_cpgs * 100
    }
    
    # Save occupancy statistics
    with open(os.path.join(analysis_dir, f'{cell_type}_occupancy_stats.txt'), 'w') as f:
        f.write(f"CpG Island Occupancy Analysis - {cell_type}\n")
        f.write("="*50 + "\n\n")
        f.write(f"Total CpG islands: {total_cpgs}\n")
        f.write(f"Endogenous MeCP2 occupancy: {occupancy_stats['endo']:.1f}%\n")
        f.write(f"Exogenous MeCP2 occupancy: {occupancy_stats['exo']:.1f}%\n")
        f.write(f"Dual occupancy: {occupancy_stats['both']:.1f}%\n")
    
    # Create visualizations
    plt.figure(figsize=(10, 6))
    sns.barplot(x=['Endogenous', 'Exogenous', 'Dual'], 
                y=[occupancy_stats['endo'], occupancy_stats['exo'], occupancy_stats['both']])
    plt.title(f'{cell_type}: CpG Island Occupancy by MeCP2')
    plt.ylabel('Percentage of CpG Islands')
    plt.savefig(os.path.join(analysis_dir, f'{cell_type}_occupancy_patterns.pdf'))
    plt.close()
    
    return occupancy_stats

def analyze_regulatory_features(df: pd.DataFrame, cell_type: str, output_dir: str):
    """
    Analyze the relationship between MeCP2 binding, gene regulation, and genomic features
    """
    analysis_dir = os.path.join(output_dir, 'regulatory_analysis')
    os.makedirs(analysis_dir, exist_ok=True)
    
    # Focus on genes with differential MeCP2 binding
    diff_bound = df[df['binding_type'] == 'both'].copy()
    diff_bound['exo_endo_ratio'] = diff_bound['exo_signal'] / diff_bound['endo_signal']
    
    # Analyze expression patterns
    high_exo = diff_bound[diff_bound['exo_endo_ratio'] > 1]
    
    # Create summary plots
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=high_exo, x='expression_status', 
                y='cpg_island_mean_methylation', hue='binding_type')
    plt.title(f'{cell_type}: Methylation in Genes with High Exo/Endo Ratio')
    plt.savefig(os.path.join(analysis_dir, f'{cell_type}_high_exo_methylation.pdf'))
    plt.close()
    
    # Save detailed results
    high_exo.to_csv(os.path.join(analysis_dir, f'{cell_type}_high_exo_genes.csv'))
    
    return {
        'high_exo_genes': high_exo,
        'diff_bound_genes': diff_bound
    }

def analyze_exo_enriched_methylation(df: pd.DataFrame, cell_type: str, output_dir: str):
    """
    Analyze methylation patterns in TSS-associated CpG islands and gene bodies
    specifically for exo-enriched samples (including exo-only), comparing up/down regulated genes.
    """
    logger.info(f"\nAnalyzing Exo-enriched methylation patterns for {cell_type}...")
    
    # Create output directory with explicit path
    analysis_dir = os.path.join(output_dir, 'exo_enriched_analysis')
    os.makedirs(analysis_dir, exist_ok=True)
    logger.info(f"Created analysis directory: {analysis_dir}")
    
    # Define exo-enriched samples:
    # 1. Get exo-only samples
    exo_only = df[df['binding_type'] == 'exo'].copy()
    logger.info(f"Found {len(exo_only)} exo-only regions")
    
    # 2. Select regions with both exo and endo binding and calculate ratio
    common_targets = df[df['binding_type'] == 'both'].copy()
    common_targets['exo_endo_ratio'] = common_targets['exo_signal'] / common_targets['endo_signal']
    
    # 3. Define exo-enriched from common targets as those with ratio > 1
    exo_enriched_common = common_targets[common_targets['exo_endo_ratio'] > 1]
    logger.info(f"Found {len(exo_enriched_common)} exo-enriched regions from common targets")
    
    # 4. Combine exo-only and exo-enriched common targets
    exo_enriched = pd.concat([exo_only, exo_enriched_common])
    logger.info(f"Total exo-enriched regions: {len(exo_enriched)}")
    
    # Save the exo-enriched genes list with source information
    exo_enriched['source'] = exo_enriched.apply(
        lambda x: 'exo_only' if x['binding_type'] == 'exo' else 'exo_enriched_common',
        axis=1
    )
    exo_enriched.to_csv(os.path.join(analysis_dir, f'{cell_type}_exo_enriched_genes.csv'), index=False)
    
    # Group by expression status
    expression_groups = {
        'upregulated': exo_enriched[exo_enriched['expression_status'] == 'upregulated'],
        'downregulated': exo_enriched[exo_enriched['expression_status'] == 'downregulated'],
        'unchanged': exo_enriched[exo_enriched['expression_status'] == 'unchanged']
    }
    
    # Calculate statistics and create report
    stats_file = os.path.join(analysis_dir, f'{cell_type}_methylation_analysis.txt')
    with open(stats_file, 'w') as f:
        f.write(f"Exo-enriched Methylation Analysis - {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        f.write("1. Overall Statistics\n")
        f.write("-"*20 + "\n")
        f.write(f"Total exo-enriched regions: {len(exo_enriched)}\n")
        f.write(f"- Exo-only regions: {len(exo_only)}\n")
        f.write(f"- Exo-enriched common regions: {len(exo_enriched_common)}\n")
        if len(exo_enriched_common) > 0:
            f.write(f"Mean exo/endo ratio in common regions: {exo_enriched_common['exo_endo_ratio'].mean():.2f}\n")
        f.write("\n")
        
        f.write("2. Expression Groups\n")
        f.write("-"*20 + "\n")
        for status, group in expression_groups.items():
            f.write(f"\n{status.title()}:\n")
            f.write(f"Count: {len(group)}\n")
            if len(group) > 0:
                # TSS methylation stats
                f.write("TSS CpG islands:\n")
                f.write(f"- Mean methylation: {group['cpg_island_mean_methylation'].mean():.2f}%\n")
                f.write(f"- Median methylation: {group['cpg_island_mean_methylation'].median():.2f}%\n")
                f.write(f"- Std Dev: {group['cpg_island_mean_methylation'].std():.2f}%\n")
                
                # Source breakdown
                exo_only_count = len(group[group['source'] == 'exo_only'])
                common_count = len(group[group['source'] == 'exo_enriched_common'])
                f.write(f"\nSource breakdown:\n")
                f.write(f"- Exo-only: {exo_only_count} ({exo_only_count/len(group)*100:.1f}%)\n")
                f.write(f"- Exo-enriched common: {common_count} ({common_count/len(group)*100:.1f}%)\n")
    
    # Create visualizations
    
    # 1. TSS CpG island methylation comparison by expression and source
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=exo_enriched, 
                x='expression_status', 
                y='cpg_island_mean_methylation',
                hue='source')
    plt.title(f'{cell_type}: TSS CpG Island Methylation by Expression and Source')
    plt.xlabel('Expression Status')
    plt.ylabel('Methylation Level (%)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_dir, f'{cell_type}_tss_methylation_by_source.pdf'))
    plt.close()
    
    # 2. Methylation distribution
    plt.figure(figsize=(10, 6))
    sns.kdeplot(data=exo_enriched, 
                x='cpg_island_mean_methylation',
                hue='expression_status',
                common_norm=False)
    plt.title(f'{cell_type}: Distribution of TSS Methylation Levels')
    plt.xlabel('Methylation Level (%)')
    plt.ylabel('Density')
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_dir, f'{cell_type}_methylation_distribution.pdf'))
    plt.close()
    
    # 3. Exo signal vs Methylation colored by expression
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=exo_enriched, 
                    x='exo_signal', 
                    y='cpg_island_mean_methylation',
                    hue='expression_status',
                    style='source',
                    alpha=0.6)
    plt.title(f'{cell_type}: Exo Signal vs Methylation')
    plt.xlabel('Exo Signal')
    plt.ylabel('Methylation Level (%)')
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_dir, f'{cell_type}_signal_vs_methylation.pdf'))
    plt.close()
    
    return {
        'exo_enriched': exo_enriched,
        'expression_groups': expression_groups,
        'source_breakdown': {
            'exo_only': exo_only,
            'exo_enriched_common': exo_enriched_common
        }
    }

def verify_directory_structure(base_dir: str):
    """Verify the analysis directory structure exists and is writable"""
    required_dirs = [
        'cpg_analysis',
        'detailed_analysis',
        os.path.join('cpg_analysis', 'NSC'),
        os.path.join('detailed_analysis', 'exo_enriched_analysis')
    ]
    
    logger.info("Verifying directory structure...")
    for dir_path in required_dirs:
        full_path = os.path.join(base_dir, dir_path)
        try:
            # Create directory if it doesn't exist
            os.makedirs(full_path, exist_ok=True)
            
            # Test if directory is writable by creating and removing a test file
            test_file = os.path.join(full_path, '.test')
            try:
                with open(test_file, 'w') as f:
                    f.write('test')
                if os.path.exists(test_file):  # Verify file was created
                    os.remove(test_file)  # Clean up
                    logger.info(f"Verified directory is writable: {full_path}")
                else:
                    raise OSError(f"Could not create test file in {full_path}")
            except Exception as e:
                logger.error(f"Directory {full_path} is not writable: {str(e)}")
                raise OSError(f"Directory {full_path} is not writable: {str(e)}")
                
        except Exception as e:
            logger.error(f"Error with directory {full_path}: {str(e)}")
            # Try to create parent directories
            parent_dir = os.path.dirname(full_path)
            if not os.path.exists(parent_dir):
                try:
                    os.makedirs(parent_dir, exist_ok=True)
                    logger.info(f"Created parent directory: {parent_dir}")
                except Exception as pe:
                    logger.error(f"Could not create parent directory {parent_dir}: {str(pe)}")
            raise

def validate_mecp2_data(file_path: str) -> pd.DataFrame:
    """
    Validate and load MeCP2 data with detailed error checking
    """
    logger.info(f"Attempting to load MeCP2 data from: {file_path}")
    
    try:
        df = pd.read_csv(file_path)
        logger.info(f"Successfully loaded MeCP2 data: {len(df)} rows")
        logger.info(f"Columns found: {df.columns.tolist()}")
        return df
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
        # Check parent directory
        parent_dir = os.path.dirname(file_path)
        if os.path.exists(parent_dir):
            logger.info(f"Contents of {parent_dir}:")
            for item in os.listdir(parent_dir):
                logger.info(f"  - {item}")
        raise
    except Exception as e:
        logger.error(f"Error loading MeCP2 data: {str(e)}")
        raise

def save_tss_analysis_results(results: dict, cell_type: str, output_dir: str):
    """Save detailed statistics from TSS binding analysis"""
    
    # Create output file
    stats_file = os.path.join(output_dir, f'{cell_type}_tss_analysis_stats.txt')
    
    with open(stats_file, 'w') as f:
        f.write(f"TSS Binding Analysis Results - {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        for category, stats in results.items():
            f.write(f"\n{category.upper()}\n")
            f.write("-"*30 + "\n")
            
            total_genes = sum(stats.get(status, {}).get('count', 0) 
                            for status in ['not_deregulated', 'upregulated', 'downregulated'])
            
            f.write(f"Total genes: {total_genes}\n\n")
            
            for status, measurements in stats.items():
                f.write(f"{status}:\n")
                f.write(f"  Count: {measurements['count']}\n")
                f.write(f"  Percentage: {(measurements['count']/total_genes)*100:.2f}%\n")
                
                f.write("\n  Promoter Methylation:\n")
                f.write(f"    Mean ± SD: {measurements['promoter_methylation']['mean']:.2f} ± "
                       f"{measurements['promoter_methylation']['std']:.2f}\n")
                f.write(f"    Median: {measurements['promoter_methylation']['median']:.2f}\n")
                
                f.write("\n  Gene Body Methylation:\n")
                f.write(f"    Mean ± SD: {measurements['gene_body_methylation']['mean']:.2f} ± "
                       f"{measurements['gene_body_methylation']['std']:.2f}\n")
                f.write(f"    Median: {measurements['gene_body_methylation']['median']:.2f}\n\n")

def analyze_tss_binding_patterns(df: pd.DataFrame, cell_type: str, base_output_dir: str):
    """Analyze binding patterns specifically around TSS"""
    tss_dir = os.path.join(base_output_dir, 'tss_analysis', cell_type)
    os.makedirs(tss_dir, exist_ok=True)
    logger.info(f"Using TSS analysis directory: {tss_dir}")
    
    # Check if we have methylation data
    has_methylation = all(col in df.columns for col in ['cpg_island_count', 'cpg_island_mean_methylation'])
    logger.info(f"Methylation data available: {has_methylation}")
    
    # Define binding categories
    binding_categories = {
        'exo_enriched': df[df['binding_type'].isin(['exo', 'both'])],
        'exo_only': df[df['binding_type'] == 'exo'],
        'endo_only': df[df['binding_type'] == 'endo'],
        'non_enriched': df[~df['mecp2_bound']]
    }
    
    # Analyze each category
    results = {}
    for category, category_df in binding_categories.items():
        # Skip empty categories
        if len(category_df) == 0:
            logger.warning(f"No genes found in group: {category} for {cell_type}")
            continue
            
        # Split by regulation status
        regulation_groups = {
            'not_deregulated': category_df[category_df['expression_status'] == 'unchanged'],
            'upregulated': category_df[category_df['expression_status'] == 'upregulated'],
            'downregulated': category_df[category_df['expression_status'] == 'downregulated']
        }
        
        # Calculate statistics for each group
        group_stats = {}
        for status, group_df in regulation_groups.items():
            if len(group_df) == 0:
                continue
                
            stats = {
                'count': len(group_df),
                'promoter_methylation': {
                    'mean': group_df['cpg_island_mean_methylation'].mean(),
                    'median': group_df['cpg_island_mean_methylation'].median(),
                    'std': group_df['cpg_island_mean_methylation'].std()
                }
            }
            group_stats[status] = stats
            
        results[category] = group_stats
    
    # Save results
    save_tss_analysis_results(results, cell_type, tss_dir)
    
    return results

def analyze_methylation_patterns(df: pd.DataFrame, cell_type: str, base_output_dir: str):
    """Analyze methylation patterns"""
    analysis_dir = os.path.join(base_output_dir, 'methylation_analysis', cell_type)
    os.makedirs(analysis_dir, exist_ok=True)
    logger.info(f"Using methylation analysis directory: {analysis_dir}")
    
    logger.info(f"\nAnalyzing methylation patterns for {cell_type}...")
    
    # 1. Calculate summary statistics for CpG islands only
    methylation_stats = {
        'cpg_island': {
            'mean': df['cpg_island_mean_methylation'].mean(),
            'median': df['cpg_island_mean_methylation'].median(),
            'std': df['cpg_island_mean_methylation'].std()
        }
    }
    
    # 2. Create visualizations
    
    # CpG Island methylation distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(data=df, x='cpg_island_mean_methylation', bins=50)
    plt.title(f'{cell_type}: CpG Island Methylation Distribution')
    plt.xlabel('Methylation Level (%)')
    plt.ylabel('Count')
    plt.savefig(os.path.join(analysis_dir, 'cpg_methylation_dist.pdf'))
    plt.close()
    
    # Save statistics to file
    stats_file = os.path.join(analysis_dir, 'methylation_statistics.txt')
    with open(stats_file, 'w') as f:
        f.write(f"Methylation Analysis Results - {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        f.write("CpG Island Methylation:\n")
        f.write("-"*30 + "\n")
        f.write(f"Mean: {methylation_stats['cpg_island']['mean']:.2f}%\n")
        f.write(f"Median: {methylation_stats['cpg_island']['median']:.2f}%\n")
        f.write(f"Std Dev: {methylation_stats['cpg_island']['std']:.2f}%\n\n")
    
    return methylation_stats


