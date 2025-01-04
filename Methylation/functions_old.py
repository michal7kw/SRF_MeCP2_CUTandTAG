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
from typing import Dict, List, Tuple, Any
import pyranges as pr
from functools import partial
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import time
import argparse
from config import CONFIG, PATHS, logger
from cache_utils import save_to_cache, load_from_cache, clear_cache


# Set default plotting style
sns.set_theme(style="whitegrid")  # Use seaborn's built-in style
plt.rcParams['figure.figsize'] = [10, 6]
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300

# set working directory
os.chdir("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Methylation")


#%% Helper functions for genomic regions
def get_promoter_region(gene_start: int, gene_end: int, strand: str) -> Tuple[int, int]:
    """Get promoter coordinates based on gene location and strand"""
    if strand == '+':
        return (gene_start + CONFIG['promoter_region'][0], 
                gene_start + CONFIG['promoter_region'][1])
    else:
        return (gene_end - CONFIG['promoter_region'][1], 
                gene_end - CONFIG['promoter_region'][0])

def get_gene_body_region(gene_start: int, gene_end: int, strand: str) -> Tuple[int, int]:
    """Get gene body coordinates based on gene location and strand"""
    if strand == '+':
        return (gene_start + CONFIG['gene_body_start'], gene_end)
    else:
        return (gene_start, gene_end - CONFIG['gene_body_start'])

#%% Analyze methylation patterns
def analyze_methylation_patterns(genes_df: pd.DataFrame, 
                               expression_data: Dict[str, pd.DataFrame],
                               mecp2_binding: pd.DataFrame,
                               medip_dir: str,
                               genome_fasta: str,
                               use_cache: bool = True) -> Dict[str, pd.DataFrame]:
    """Analyze methylation patterns with caching support"""
    
    # Try to load from cache first
    if use_cache and CONFIG['cache']['enabled']:
        cached_results = load_from_cache(CONFIG['cache']['methylation_file'])
        if cached_results is not None:
            return cached_results
    
    # If not in cache or cache disabled, perform analysis
    results = {}
    start_time = time.time()
    
    logger.info("Starting methylation pattern analysis...")
    logger.info(f"Number of genes to process: {len(genes_df)}")
    
    if CONFIG['debug']['enabled']:
        logger.info("DEBUG MODE: Processing reduced dataset")
    
    # Process MeCP2 binding data
    logger.info("Processing MeCP2 binding data...")
    mecp2_df_for_pr = mecp2_binding.copy()
    mecp2_df_for_pr = mecp2_df_for_pr.rename(columns={
        'chr': 'Chromosome',
        'start': 'Start',
        'end': 'End'
    })
    
    logger.info("Converting to PyRanges...")
    mecp2_pr = pr.PyRanges(mecp2_df_for_pr)
    
    genes_df_for_pr = genes_df.copy()
    genes_df_for_pr = genes_df_for_pr.rename(columns={
        'chr': 'Chromosome',
        'start': 'Start',
        'end': 'End'
    })
    genes_pr = pr.PyRanges(genes_df_for_pr)
    
    logger.info("Finding overlaps...")
    overlaps = mecp2_pr.join(genes_pr)
    overlaps_df = overlaps.as_df()
    
    logger.info("Creating gene mappings...")
    mecp2_gene_mapping = overlaps_df.groupby(
        ['Start', 'End', 'Chromosome'],
        observed=True
    )['gene_id'].first().reset_index()
    
    mecp2_binding = mecp2_binding.merge(
        mecp2_gene_mapping.rename(columns={
            'Chromosome': 'chr',
            'Start': 'start',
            'End': 'end'
        }),
        on=['chr', 'start', 'end'],
        how='left'
    )
    
    for cell_type, expr_df in expression_data.items():
        logger.info(f"\nProcessing {cell_type}...")
        process_start = time.time()
        
        # Merge expression data with gene annotations
        try:
            merged_df = genes_df.merge(expr_df, on='gene_id', how='inner')
            logger.info(f"Merged gene annotations with expression data. Shape: {merged_df.shape}")
        except Exception as e:
            logger.error(f"Error merging gene annotations with expression data: {str(e)}")
            continue
        
        # Add MeCP2 binding information
        try:
            merged_df = merged_df.merge(
                mecp2_binding[['gene_id', 'binding_type', 'exo_signal', 'endo_signal']].dropna(), 
                on='gene_id', 
                how='left'
            )
            logger.info(f"Merged MeCP2 binding data. Shape: {merged_df.shape}")
        except Exception as e:
            logger.error(f"Error merging MeCP2 binding data: {str(e)}")
            continue
        
        merged_df['mecp2_bound'] = merged_df['binding_type'].notna()
        
        # Calculate methylation
        logger.info("Calculating promoter methylation...")
        promoter_methylation = calculate_methylation_levels_with_replicates(
            merged_df[['chr', 'promoter_start', 'promoter_end']].rename(
                columns={'promoter_start': 'start', 'promoter_end': 'end'}
            ),
            medip_dir,
            cell_type,
            genome_fasta
        )
        
        logger.info("Calculating gene body methylation...")
        gene_body_methylation = calculate_methylation_levels_with_replicates(
            merged_df[['chr', 'gene_body_start', 'gene_body_end']].rename(
                columns={'gene_body_start': 'start', 'gene_body_end': 'end'}
            ),
            medip_dir,
            cell_type,
            genome_fasta
        )
        
        # Combine results
        results[cell_type] = pd.DataFrame({
            'gene_id': merged_df['gene_id'],
            'gene_name': merged_df['gene_name'],
            'expression_status': merged_df['expression_status'],
            'log2FoldChange': merged_df['log2FoldChange'],
            'padj': merged_df['padj'],
            'mecp2_bound': merged_df['mecp2_bound'],
            'binding_type': merged_df['binding_type'],
            'exo_signal': merged_df['exo_signal'],
            'endo_signal': merged_df['endo_signal'],
            'promoter_methylation': promoter_methylation['methylation'],
            'promoter_cpg_count': promoter_methylation['cpg_count'],
            'gene_body_methylation': gene_body_methylation['methylation'],
            'gene_body_cpg_count': gene_body_methylation['cpg_count']
        })
        
        process_time = time.time() - process_start
        logger.info(f"Completed {cell_type} processing in {process_time:.2f} seconds")
    
    total_time = time.time() - start_time
    logger.info(f"\nTotal analysis time: {total_time:.2f} seconds")
    
    # Cache the results if enabled
    if CONFIG['cache']['enabled']:
        save_to_cache(results, CONFIG['cache']['methylation_file'])
    
    return results

# Add a function to run the complete analysis pipeline
def run_analysis_pipeline(force_recompute: bool = False) -> Dict[str, Any]:
    """
    Run the complete analysis pipeline with caching support
    
    Args:
        force_recompute: If True, ignore cache and recompute everything
    """
    if force_recompute:
        clear_cache()
        logger.info("Cache cleared due to force recompute flag")
    
    # Check if we're in debug mode
    if CONFIG['debug']['enabled']:
        logger.info(f"Running pipeline in DEBUG mode with sample size: {CONFIG['debug']['sample_size']}")
    
    # Load basic data
    genes_df, gene_name_to_id = load_gene_annotations(PATHS['gtf_file'])
    expression_data = load_expression_data(PATHS['rnaseq'], gene_name_to_id)
    
    # Load MeCP2 binding data
    mecp2_binding = pd.read_csv(os.path.join(PATHS['mecp2_dir'], PATHS['mecp2_file']))
    
    # If in debug mode, sample MeCP2 binding data
    if CONFIG['debug']['enabled']:
        sample_size = CONFIG['debug']['sample_size']
        mecp2_binding = mecp2_binding.sample(n=min(len(mecp2_binding), sample_size), 
                                           random_state=42)
        logger.info(f"Debug mode: sampled {len(mecp2_binding)} MeCP2 binding regions")
    
    # Run methylation analysis
    results = analyze_methylation_patterns(
        genes_df,
        expression_data,
        mecp2_binding,
        PATHS['medip_dir'],
        PATHS['genome_fasta'],
        use_cache=not force_recompute
    )
    
    return {
        'methylation_results': results,
        'genes_df': genes_df,
        'expression_data': expression_data,
        'mecp2_binding': mecp2_binding
    }

#%% Modified print_regulated_summary function
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

#%% Calculate methylation levels for multiple replicates
def calculate_methylation_levels_with_replicates(region_df: pd.DataFrame,
                                               medip_dir: str,
                                               cell_type_prefix: str,
                                               fasta_file: str,
                                               batch_size: int = 1000) -> pd.DataFrame:
    """Calculate average methylation levels across replicates with proper normalization"""
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
            fasta = pysam.FastaFile(fasta_file)
            
            methylation_values = []
            cpg_counts = []
            
            # Process regions in batches
            total_regions = len(region_df)
            batches = range(0, total_regions, batch_size)
            
            with tqdm(total=total_regions, desc=f"Calculating methylation for {os.path.basename(rep_file)}") as pbar:
                for start_idx in batches:
                    end_idx = min(start_idx + batch_size, total_regions)
                    batch = region_df.iloc[start_idx:end_idx]
                    
                    batch_methylation = []
                    batch_cpg = []
                    
                    for _, row in batch.iterrows():
                        try:
                            # Get IP and input values
                            ip_values = bw_ip.values(row['chr'], int(row['start']), int(row['end']))
                            input_values = bw_input.values(row['chr'], int(row['start']), int(row['end']))
                            
                            if ip_values is None or input_values is None:
                                methylation = 0
                            else:
                                # Filter out None values
                                valid_pairs = [
                                    (ip, inp) for ip, inp in zip(ip_values, input_values)
                                    if ip is not None and inp is not None
                                ]
                                
                                if valid_pairs:
                                    ip_vals, input_vals = zip(*valid_pairs)
                                    
                                    # Calculate normalized methylation
                                    # Use log2 ratio and scale to 0-100 range
                                    ip_mean = np.mean(ip_vals)
                                    input_mean = np.mean(input_vals)
                                    
                                    if input_mean > 0 and ip_mean > 0:
                                        log2_ratio = np.log2(ip_mean / input_mean)
                                        # Scale log2 ratio to 0-100 range
                                        # Typically, log2 ratios range from -3 to 3
                                        methylation = max(0, min(100, (log2_ratio + 3) * 16.67))
                                    else:
                                        methylation = 0
                                else:
                                    methylation = 0
                            
                            batch_methylation.append(methylation)
                            
                            # Count CpGs (only need to do this once)
                            if len(replicate_results) == 0:
                                sequence = fasta.fetch(row['chr'], int(row['start']), int(row['end']))
                                cpg_count = sequence.upper().count('CG')
                                batch_cpg.append(cpg_count)
                                
                        except Exception as e:
                            logger.debug(f"Error processing region {row['chr']}:{row['start']}-{row['end']}: {str(e)}")
                            batch_methylation.append(0)
                            if len(replicate_results) == 0:
                                batch_cpg.append(0)
                    
                    methylation_values.extend(batch_methylation)
                    if len(replicate_results) == 0:
                        cpg_counts.extend(batch_cpg)
                    pbar.update(len(batch))
            
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
        'cpg_count': cpg_counts if cpg_counts else [0] * len(avg_methylation)
    })

    #%% Visualization functions
def create_methylation_plots(results: Dict[str, pd.DataFrame], output_dir: str):
    """Create comprehensive visualization of methylation patterns"""
    if CONFIG['debug']['enabled']:
        output_dir = os.path.join(output_dir, 'debug_output')
        logger.info(f"Debug mode: saving plots to {output_dir}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    for cell_type, df in results.items():
        # 1. Promoter methylation by expression status
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=df, x='expression_status', y='promoter_methylation', 
                   hue='mecp2_bound')
        plt.title(f'{cell_type}: Promoter Methylation by Expression Status')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{cell_type}_promoter_methylation.pdf')
        plt.close()
        
        # 2. Gene body methylation by expression status
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=df, x='expression_status', y='gene_body_methylation',
                   hue='mecp2_bound')
        plt.title(f'{cell_type}: Gene Body Methylation by Expression Status')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{cell_type}_gene_body_methylation.pdf')
        plt.close()
        
        # 3. CpG density analysis
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        sns.boxplot(data=df, x='expression_status', y='promoter_cpg_count',
                   hue='mecp2_bound', ax=ax1)
        ax1.set_title('Promoter CpG Count')
        ax1.tick_params(axis='x', rotation=45)
        
        sns.boxplot(data=df, x='expression_status', y='gene_body_cpg_count',
                   hue='mecp2_bound', ax=ax2)
        ax2.set_title('Gene Body CpG Count')
        ax2.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{cell_type}_cpg_density.pdf')
        plt.close()


#%% Statistical analysis
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


#%% Create visualization functions
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


#%% Debug the filtering steps
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


#%% Debug data merging
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


#%% Modified statistical test functions
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
                    f.write("\nCorrelation Analysis:\n")
                    f.write(f"Promoter-Gene Body Correlation:\n")
                    f.write(f"rho: {corr_result['statistic']:.4f}\n")
                    f.write(f"p-value: {corr_result['pvalue']:.4f}\n")
                    f.write(f"n: {corr_result['n']}\n")

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

#%% Create focused visualization of significant findings
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


#%% Modified plot_mecp2_regulatory_patterns function
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


#%% Create focused visualization of significant findings
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

#%% Modified print_regulated_summary function
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

#%% Print summary for regulated genes
def print_regulated_summary(results: Dict[str, pd.DataFrame]):
    """Print summary statistics for regulated MeCP2-bound genes"""
    for cell_type, df in results.items():
        print(f"\nSummary for {cell_type}:")
        print("="*50)
        
        # Filter for regulated genes
        regulated_df = df[
            df['mecp2_bound'] & 
            (df['expression_status'].isin(['upregulated', 'downregulated']))
        ]
        
        # Separate by binding type
        for binding_type in ['exo', 'endo']:
            subset_df = regulated_df[regulated_df['binding_type'] == binding_type]
            
            print(f"\n{binding_type.upper()} BINDING:")
            print("-"*20)
            
            # Count genes by expression status
            print("\nGene counts by expression status:")
            print(subset_df['expression_status'].value_counts())
            
            # Mean methylation by expression status
            print("\nMean promoter methylation:")
            print(subset_df.groupby('expression_status')['promoter_methylation'].mean())
            
            print("\nMean gene body methylation:")
            print(subset_df.groupby('expression_status')['gene_body_methylation'].mean())
            
            # Mean CpG density
            print("\nMean promoter CpG count:")
            print(subset_df.groupby('expression_status')['promoter_cpg_count'].mean())

#%% Create detailed visualizations for MeCP2-bound regulated genes
def create_mecp2_regulated_analysis(results: Dict[str, pd.DataFrame], output_dir: str):
    """Create detailed analysis of methylation patterns for MeCP2-bound regulated genes"""
    os.makedirs(os.path.join(output_dir, 'mecp2_regulated_analysis'), exist_ok=True)
    
    for cell_type, df in results.items():
        # Filter for MeCP2-bound and regulated genes only
        mecp2_bound_df = df[
            df['mecp2_bound'] & 
            (df['expression_status'].isin(['upregulated', 'downregulated']))
        ]
        
        # Separate by binding type
        exo_df = mecp2_bound_df[mecp2_bound_df['binding_type'] == 'exo']
        endo_df = mecp2_bound_df[mecp2_bound_df['binding_type'] == 'endo']
        
        # 1. Methylation distribution for exo binding
        if len(exo_df) > 0:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # Promoter methylation
            sns.boxplot(data=exo_df, 
                       x='expression_status', 
                       y='promoter_methylation',
                       ax=ax1)
            ax1.set_title(f'{cell_type}: Exo Binding\nPromoter Methylation')
            ax1.tick_params(axis='x', rotation=45)
            
            # Gene body methylation
            sns.boxplot(data=exo_df,
                       x='expression_status',
                       y='gene_body_methylation',
                       ax=ax2)
            ax2.set_title(f'{cell_type}: Exo Binding\nGene Body Methylation')
            ax2.tick_params(axis='x', rotation=45)
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'mecp2_regulated_analysis', 
                                    f'{cell_type}_exo_methylation.pdf'))
            plt.close()
        
        # 2. Methylation distribution for endo binding
        if len(endo_df) > 0:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # Promoter methylation
            sns.boxplot(data=endo_df, 
                       x='expression_status', 
                       y='promoter_methylation',
                       ax=ax1)
            ax1.set_title(f'{cell_type}: Endo Binding\nPromoter Methylation')
            ax1.tick_params(axis='x', rotation=45)
            
            # Gene body methylation
            sns.boxplot(data=endo_df,
                       x='expression_status',
                       y='gene_body_methylation',
                       ax=ax2)
            ax2.set_title(f'{cell_type}: Endo Binding\nGene Body Methylation')
            ax2.tick_params(axis='x', rotation=45)
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'mecp2_regulated_analysis', 
                                    f'{cell_type}_endo_methylation.pdf'))
            plt.close()
        
        # 3. CpG density analysis
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        if len(exo_df) > 0:
            sns.boxplot(data=exo_df,
                       x='expression_status',
                       y='promoter_cpg_count',
                       ax=ax1)
            ax1.set_title('Exo Binding: CpG Density')
            ax1.tick_params(axis='x', rotation=45)
        
        if len(endo_df) > 0:
            sns.boxplot(data=endo_df,
                       x='expression_status',
                       y='promoter_cpg_count',
                       ax=ax2)
            ax2.set_title('Endo Binding: CpG Density')
            ax2.tick_params(axis='x', rotation=45)
        
        plt.suptitle(f'{cell_type}: CpG Density in Regulated Genes', y=1.02)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'mecp2_regulated_analysis',
                                f'{cell_type}_cpg_density.pdf'))
        plt.close()
        
        # 4. Statistical analysis
        stats_summary = {
            'exo': {'promoter': {}, 'gene_body': {}},
            'endo': {'promoter': {}, 'gene_body': {}}
        }
        
        # Analyze each binding type separately
        for binding_type, subset_df in [('exo', exo_df), ('endo', endo_df)]:
            if len(subset_df) > 0:
                for region in ['promoter', 'gene_body']:
                    # Compare methylation between up and down regulated genes
                    up_meth = subset_df[subset_df['expression_status'] == 'upregulated'][f'{region}_methylation']
                    down_meth = subset_df[subset_df['expression_status'] == 'downregulated'][f'{region}_methylation']
                    
                    if len(up_meth) > 0 and len(down_meth) > 0:
                        stat, pval = stats.mannwhitneyu(up_meth, down_meth)
                        stats_summary[binding_type][region]['methylation_diff'] = {
                            'statistic': stat,
                            'pvalue': pval,
                            'up_mean': up_meth.mean(),
                            'down_mean': down_meth.mean(),
                            'up_count': len(up_meth),
                            'down_count': len(down_meth)
                        }
                    
                    # Compare CpG density
                    up_cpg = subset_df[subset_df['expression_status'] == 'upregulated'][f'{region}_cpg_count']
                    down_cpg = subset_df[subset_df['expression_status'] == 'downregulated'][f'{region}_cpg_count']
                    
                    if len(up_cpg) > 0 and len(down_cpg) > 0:
                        stat, pval = stats.mannwhitneyu(up_cpg, down_cpg)
                        stats_summary[binding_type][region]['cpg_density_diff'] = {
                            'statistic': stat,
                            'pvalue': pval,
                            'up_mean': up_cpg.mean(),
                            'down_mean': down_cpg.mean()
                        }
        
        # Save statistical summary
        with open(os.path.join(output_dir, 'mecp2_regulated_analysis',
                              f'{cell_type}_regulated_stats.txt'), 'w') as f:
            f.write(f"Statistical Analysis for {cell_type} MeCP2-Bound Regulated Genes\n")
            f.write("="*60 + "\n\n")
            
            for binding_type in ['exo', 'endo']:
                f.write(f"\n{binding_type.upper()} BINDING:\n")
                f.write("="*20 + "\n")
                
                for region, stats_dict in stats_summary[binding_type].items():
                    f.write(f"\n{region.upper()}:\n")
                    f.write("-"*20 + "\n")
                    
                    for analysis, values in stats_dict.items():
                        f.write(f"{analysis}:\n")
                        for key, value in values.items():
                            f.write(f"{key}: {value}\n")

#%% Load and process gene annotations with gene name mapping
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

#%% Load and process differential expression data with gene ID mapping
def load_expression_data(rnaseq_files: Dict[str, str], gene_name_to_id: dict) -> Dict[str, pd.DataFrame]:
    """Load expression data with debug mode support"""
    logger.info("Loading differential expression data...")
    
    expression_data = {}
    for cell_type, file_path in rnaseq_files.items():
        df = pd.read_csv(file_path)
        
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
            logger.info(f"Debug mode: sampled {len(df)} genes for {cell_type}")
        
        # Classify genes
        df['expression_status'] = 'unchanged'
        df.loc[(df['log2FoldChange'] > CONFIG['expression_thresholds']['log2fc']) & 
               (df['padj'] < CONFIG['expression_thresholds']['padj']), 
               'expression_status'] = 'upregulated'
        df.loc[(df['log2FoldChange'] < -CONFIG['expression_thresholds']['log2fc']) & 
               (df['padj'] < CONFIG['expression_thresholds']['padj']), 
               'expression_status'] = 'downregulated'
        
        expression_data[cell_type] = df
    
    return expression_data

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
        
        # Safely calculate correlation
        try:
            # Remove NaN values and calculate correlation
            valid_data = df[['promoter_methylation', 'gene_body_methylation']].dropna()
            if len(valid_data) > 1:  # Need at least 2 points for correlation
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
