#Import required libraries
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
from typing import Tuple

# Set default plotting style
sns.set_theme(style="whitegrid")  # Use seaborn's built-in style
plt.rcParams['figure.figsize'] = [10, 6]
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300

# set working directory
os.chdir("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Methylation_dev")

def get_promoter_region(gene_start: int, gene_end: int, strand: str) -> Tuple[int, int]:
    """Get promoter coordinates based on gene location and strand
    
    For genes on the + strand:
    - Takes CONFIG['genomic_regions']['promoter'][0] bases upstream of TSS as start
    - Takes CONFIG['genomic_regions']['promoter'][1] bases upstream of TSS as end
    
    For genes on the - strand:
    - Takes CONFIG['genomic_regions']['promoter'][1] bases downstream of TTS as start 
    - Takes CONFIG['genomic_regions']['promoter'][0] bases downstream of TTS as end
    
    Args:
        gene_start: Start coordinate of the gene
        gene_end: End coordinate of the gene  
        strand: DNA strand ('+' or '-')
        
    Returns:
        Tuple of (promoter_start, promoter_end) coordinates
    """
    if strand == '+':
        # For + strand, promoter is upstream of gene start (TSS)
        return (gene_start + CONFIG['genomic_regions']['promoter'][0], 
                gene_start + CONFIG['genomic_regions']['promoter'][1])
    else:
        # For - strand, promoter is downstream of gene end (TSS) 
        return (gene_end - CONFIG['genomic_regions']['promoter'][1], 
                gene_end - CONFIG['genomic_regions']['promoter'][0])

def get_gene_body_region(gene_start: int, gene_end: int, strand: str) -> Tuple[int, int]:
    """Get gene body coordinates based on gene location and strand
    
    For genes on the + strand:
    - Starts CONFIG['genomic_regions']['gene_body_start'] bases downstream of TSS
    - Ends at gene end coordinate
    
    For genes on the - strand:  
    - Starts at gene start coordinate
    - Ends CONFIG['genomic_regions']['gene_body_start'] bases upstream of TTS
    
    Args:
        gene_start: Start coordinate of the gene
        gene_end: End coordinate of the gene
        strand: DNA strand ('+' or '-')
        
    Returns:
        Tuple of (gene_body_start, gene_body_end) coordinates
    """
    if strand == '+':
        # For + strand, gene body starts downstream of TSS
        return (gene_start + CONFIG['genomic_regions']['gene_body_start'], gene_end)
    else:
        # For - strand, gene body ends upstream of TTS
        return (gene_start, gene_end - CONFIG['genomic_regions']['gene_body_start'])

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

# Add a function to run the complete analysis pipeline
def run_analysis_pipeline(force_recompute: bool = False) -> Dict[str, Any]:
    """Run the complete analysis pipeline with caching support"""
    if force_recompute:
        clear_cache()
        logger.info("Cache cleared due to force recompute flag")
    
    # Load all data first
    mecp2_binding = pd.read_csv(os.path.join(PATHS['mecp2_dir'], PATHS['mecp2_file']))
    genes_df, gene_name_to_id = load_gene_annotations(PATHS['gtf_file'])
    expression_data = load_expression_data(PATHS['rnaseq'], gene_name_to_id)
    
    # If in debug mode, sample smartly
    if CONFIG['debug']['enabled']:
        sample_size = CONFIG['debug']['sample_size']
        
        # First, get genes with expression data
        expressed_genes = set()
        for cell_type, expr_df in expression_data.items():
            expressed_genes.update(expr_df['gene_id'])
        
        # Filter genes to those with expression data
        genes_df = genes_df[genes_df['gene_id'].isin(expressed_genes)]
        logger.info(f"Found {len(genes_df)} genes with expression data")
        
        # Find genes with binding among expressed genes
        binding_regions = pr.PyRanges(
            chromosomes=mecp2_binding['chr'],
            starts=mecp2_binding['start'].astype(int),
            ends=mecp2_binding['end'].astype(int)
        )
        
        genes_pr = pr.PyRanges(
            chromosomes=genes_df['chr'],
            starts=genes_df['promoter_start'].astype(int),
            ends=genes_df['gene_body_end'].astype(int)
        )
        genes_pr.gene_id = genes_df['gene_id'].values
        
        # Find overlaps
        overlaps = binding_regions.join(genes_pr)
        if overlaps is not None and len(overlaps) > 0:
            bound_genes = pd.unique(overlaps.as_df()['gene_id'])
            logger.info(f"Found {len(bound_genes)} expressed genes with MeCP2 binding")
            
            # Sample bound genes by chromosome
            bound_genes_df = genes_df[genes_df['gene_id'].isin(bound_genes)]
            chrs_with_binding = bound_genes_df['chr'].unique()
            
            # Calculate target number of bound genes
            target_n_bound = min(sample_size // 2, len(bound_genes))
            genes_per_chr = max(1, target_n_bound // len(chrs_with_binding))
            
            # Sample genes from each chromosome with binding
            sampled_bound = []
            for chr_name in chrs_with_binding:
                chr_bound = bound_genes_df[bound_genes_df['chr'] == chr_name]['gene_id'].values
                if len(chr_bound) > 0:
                    n_from_chr = min(genes_per_chr, len(chr_bound))
                    sampled_bound.extend(np.random.choice(chr_bound, size=n_from_chr, replace=False))
            
            n_bound = len(sampled_bound)
            logger.info(f"Sampled {n_bound} bound genes from {len(chrs_with_binding)} chromosomes")
            
            # Sample unbound expressed genes from same chromosomes
            unbound_genes = genes_df[~genes_df['gene_id'].isin(bound_genes)]
            n_unbound = sample_size - n_bound
            genes_per_chr = max(1, n_unbound // len(chrs_with_binding))
            
            sampled_unbound = []
            for chr_name in chrs_with_binding:
                chr_unbound = unbound_genes[unbound_genes['chr'] == chr_name]['gene_id'].values
                if len(chr_unbound) > 0:
                    n_from_chr = min(genes_per_chr, len(chr_unbound))
                    sampled_unbound.extend(np.random.choice(chr_unbound, size=n_from_chr, replace=False))
            
            # Combine samples
            sampled_genes = np.concatenate([sampled_bound, sampled_unbound])
            genes_df = genes_df[genes_df['gene_id'].isin(sampled_genes)]
            
            # Get binding regions that overlap with sampled genes
            overlapping_regions = overlaps.as_df()
            overlapping_regions = overlapping_regions[
                overlapping_regions['gene_id'].isin(sampled_genes)
            ]
            binding_indices = overlapping_regions.index.unique()
            mecp2_binding = mecp2_binding.loc[binding_indices]
            
            logger.info(f"Final sample: {len(genes_df)} genes ({n_bound} bound, {len(sampled_unbound)} unbound)")
            logger.info(f"Binding regions: {len(mecp2_binding)}")
            logger.info(f"Chromosomes: {sorted(genes_df['chr'].unique())}")
        else:
            logger.warning("No overlaps found, falling back to random sampling")
            genes_df = genes_df.sample(n=sample_size, random_state=42)
            mecp2_binding = mecp2_binding.sample(n=sample_size, random_state=42)
    
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

#Modified print_regulated_summary function
def print_regulated_summary_v2(results: Dict[str, pd.DataFrame]):
    """Print summary statistics for regulated MeCP2-bound genes with improved analysis"""
    
    # Create output directory if it doesn't exist
    output_dir = "analysis_output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Open file for writing
    with open(os.path.join(output_dir, "regulated_genes_summary.txt"), "w") as f:
        for cell_type, df in results.items():
            summary = f"\nSummary for {cell_type}:\n"
            summary += "="*50 + "\n"
            
            f.write(summary)
            print(summary, end="")
            
            # First, check regulated genes
            regulated_df = df[df['expression_status'].isin(['upregulated', 'downregulated'])]
            regulated_summary = f"\nTotal regulated genes: {len(regulated_df)}\n"
            
            f.write(regulated_summary)
            print(regulated_summary, end="")
            
            if len(regulated_df) == 0:
                msg = "No regulated genes found.\n"
                f.write(msg)
                print(msg, end="")
                continue
            
            # Check MeCP2 binding
            mecp2_bound_df = regulated_df[regulated_df['mecp2_bound']]
            binding_summary = f"MeCP2-bound regulated genes: {len(mecp2_bound_df)}\n"
            
            f.write(binding_summary)
            print(binding_summary, end="")
            
            if len(mecp2_bound_df) == 0:
                msg = "No MeCP2-bound regulated genes found.\n"
                f.write(msg)
                print(msg, end="")
                continue
            
            # Check binding types present
            binding_types = mecp2_bound_df['binding_type'].value_counts()
            binding_type_summary = "\nBinding types present:\n" + binding_types.to_string() + "\n"
            
            f.write(binding_type_summary)
            print(binding_type_summary, end="")
            
            # Analyze each binding type
            for binding_type in binding_types.index:
                if pd.isna(binding_type):
                    continue
                    
                subset_df = mecp2_bound_df[mecp2_bound_df['binding_type'] == binding_type]
                
                binding_header = f"\n{binding_type.upper()} BINDING:\n" + "-"*20 + "\n"
                f.write(binding_header)
                print(binding_header, end="")
                
                # Gene counts and basic statistics
                expr_counts = subset_df['expression_status'].value_counts()
                counts_summary = "\nGene counts by expression status:\n" + expr_counts.to_string() + "\n"
                
                f.write(counts_summary)
                print(counts_summary, end="")
                
                # Calculate and display statistics
                metrics = {
                    'promoter_methylation': 'Promoter Methylation',
                    'gene_body_methylation': 'Gene Body Methylation',
                    'promoter_cpg_count': 'Promoter CpG Count'
                }
                
                for metric, metric_name in metrics.items():
                    metric_header = f"\n{metric_name}:\n"
                    f.write(metric_header)
                    print(metric_header, end="")
                    
                    # Calculate statistics by expression status
                    stats_df = subset_df.groupby('expression_status')[metric].agg([
                        'count', 'mean', 'std', 'min', 'max'
                    ]).round(3)
                    stats_summary = stats_df.to_string() + "\n"
                    
                    f.write(stats_summary)
                    print(stats_summary, end="")
                    
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
                            stats_summary = (f"\nMann-Whitney U test for {metric_name}:\n"
                                          f"Statistic: {stats_results['mannwhitney']['statistic']:.2f}\n"
                                          f"P-value: {stats_results['mannwhitney']['pvalue']:.4f}\n")
                            
                            if stats_results.get('cohens_d') is not None:
                                stats_summary += f"Cohen's d: {stats_results['cohens_d']:.4f}\n"
                            else:
                                stats_summary += "Cohen's d: Not available\n"
                                
                            f.write(stats_summary)
                            print(stats_summary, end="")
                        
                        if stats_results.get('error'):
                            error_msg = f"Error in statistical tests: {stats_results['error']}\n"
                            f.write(error_msg)
                            print(error_msg, end="")
                
                # Calculate correlation
                if len(subset_df) > 0:
                    corr_result = calculate_correlation(
                        subset_df['promoter_methylation'],
                        subset_df['gene_body_methylation']
                    )
                    if corr_result:
                        corr_summary = ("\nCorrelation between promoter and gene body methylation:\n"
                                      f"Spearman correlation: rho={corr_result['statistic']:.4f}, "
                                      f"p={corr_result['pvalue']:.4f}\n")
                        f.write(corr_summary)
                        print(corr_summary, end="")
                    else:
                        no_corr_msg = "\nCorrelation analysis: Not available (insufficient data or no variation)\n"
                        f.write(no_corr_msg)
                        print(no_corr_msg, end="")


#########################################################################################

def calculate_contextual_methylation(merged_df: pd.DataFrame,
                                   medip_dir: str,
                                   cell_type: str,
                                   genome_fasta: str,
                                   n_processes: int = None) -> pd.DataFrame:
    """Calculate methylation levels for promoter and gene body regions.
    
    Args:
        merged_df: DataFrame with columns from validate_and_merge_data():
            - chr, gene_start, gene_end, strand (from GTF)
            - promoter_start, promoter_end (calculated)
            - gene_body_start, gene_body_end (calculated)
            - gene_id, gene_name (identifiers)
            - binding_type (from mecp2_binding)
            - mecp2_bound (boolean)
            - exo_signal, endo_signal (from binding data)
            - log2FoldChange, padj (from expression data)
        medip_dir: Directory containing MeDIP bigWig files
        cell_type: Cell type code ('NSC' or 'NEU')
        genome_fasta: Path to genome FASTA file
        n_processes: Number of parallel processes
    
    Returns:
        DataFrame with all columns from merged_df plus:
        - promoter_methylation: float (0-100)
        - promoter_cpg_count: int
        - gene_body_methylation: float (0-100)
        - gene_body_cpg_count: int
    """
    logger.info(f"Starting methylation calculation for {cell_type}")
    logger.info(f"Input data shape: {merged_df.shape}")
    
    # Verify required columns
    required_columns = ['chr', 'promoter_start', 'promoter_end', 
                       'gene_body_start', 'gene_body_end']
    missing_columns = [col for col in required_columns if col not in merged_df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    # Prepare region data
    promoter_data = merged_df[['chr', 'promoter_start', 'promoter_end']].rename(
        columns={'promoter_start': 'start', 'promoter_end': 'end'}
    )
    gene_body_data = merged_df[['chr', 'gene_body_start', 'gene_body_end']].rename(
        columns={'gene_body_start': 'start', 'gene_body_end': 'end'}
    )

    logger.info("Calculating methylation for promoters and gene bodies")
    # Calculate methylation for both regions in parallel
    with ThreadPoolExecutor(max_workers=2) as executor:
        promoter_future = executor.submit(
            calculate_methylation_levels_parallel,
            promoter_data, medip_dir, cell_type, genome_fasta, n_processes
        )
        gene_body_future = executor.submit(
            calculate_methylation_levels_parallel,
            gene_body_data, medip_dir, cell_type, genome_fasta, n_processes
        )
        
        try:
            promoter_methylation = promoter_future.result()
            gene_body_methylation = gene_body_future.result()
        except Exception as e:
            logger.error(f"Error in parallel methylation calculation: {str(e)}")
            raise

    # Create result DataFrame preserving all original columns
    result_df = merged_df.copy()
    
    # Add methylation data
    result_df['promoter_methylation'] = promoter_methylation['methylation']
    result_df['promoter_cpg_count'] = promoter_methylation['cpg_count']
    result_df['gene_body_methylation'] = gene_body_methylation['methylation']
    result_df['gene_body_cpg_count'] = gene_body_methylation['cpg_count']

    # Fill NaN values if any
    result_df['promoter_methylation'] = result_df['promoter_methylation'].fillna(0)
    result_df['gene_body_methylation'] = result_df['gene_body_methylation'].fillna(0)
    result_df['promoter_cpg_count'] = result_df['promoter_cpg_count'].fillna(0)
    result_df['gene_body_cpg_count'] = result_df['gene_body_cpg_count'].fillna(0)

    # Ensure binding columns exist
    if 'binding_type' not in result_df.columns:
        result_df['binding_type'] = 'none'
    if 'mecp2_bound' not in result_df.columns:
        result_df['mecp2_bound'] = False

    logger.info(f"Completed methylation calculation. Output shape: {result_df.shape}")
    logger.info("Sample of results:")
    logger.info(result_df[['gene_id', 'promoter_methylation', 
                          'gene_body_methylation', 'binding_type']].head())

    return result_df


#########################################################################################

#Visualization functions
def create_methylation_plots(results: Dict[str, pd.DataFrame], output_dir: str):
    """Create comprehensive visualization of methylation patterns
    
    Args:
        results: Dictionary containing DataFrames with methylation data
        output_dir: Base output directory for saving plots
    """
    plot_dir = os.path.join(output_dir, 'methylation_plots')
    os.makedirs(plot_dir, exist_ok=True)
    
    for cell_type, df in results.items():
        # Add expression status based on log2FoldChange and padj
        df = df.copy()  # Create a copy to avoid modifying original
        df['expression_status'] = 'unchanged'
        df.loc[(df['log2FoldChange'] > 1) & (df['padj'] < 0.05), 'expression_status'] = 'upregulated'
        df.loc[(df['log2FoldChange'] < -1) & (df['padj'] < 0.05), 'expression_status'] = 'downregulated'
        
        # Filter out low quality data
        high_quality_df = df[
            ~df['low_cpg_promoter'] & 
            ~df['high_meth_promoter']
        ]
        
        # Log data quality info
        logger.info(f"\nMethylation plot data for {cell_type}:")
        logger.info(f"Total genes: {len(df)}")
        logger.info(f"High quality genes: {len(high_quality_df)}")
        logger.info("Expression status distribution:")
        logger.info(high_quality_df['expression_status'].value_counts().to_string())
        
        # 1. Promoter methylation by expression status
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=high_quality_df, 
                   x='binding_type', 
                   y='promoter_methylation',
                   hue='expression_status')
        plt.title(f'{cell_type}: Promoter Methylation by Binding Type')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, f'{cell_type}_promoter_methylation.pdf'))
        plt.close()
        
        # 2. Gene body methylation by expression status
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=high_quality_df,
                   x='binding_type',
                   y='gene_body_methylation',
                   hue='expression_status')
        plt.title(f'{cell_type}: Gene Body Methylation by Binding Type')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, f'{cell_type}_gene_body_methylation.pdf'))
        plt.close()
        
        # 3. CpG density analysis
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Promoter CpG metrics
        sns.boxplot(data=high_quality_df,
                   x='binding_type',
                   y='promoter_cpg_count',
                   hue='expression_status',
                   ax=ax1)
        ax1.set_title('Promoter CpG Count')
        ax1.tick_params(axis='x', rotation=45)
        
        sns.boxplot(data=high_quality_df,
                   x='binding_type',
                   y='promoter_cpg_density',
                   hue='expression_status',
                   ax=ax2)
        ax2.set_title('Promoter CpG Density')
        ax2.tick_params(axis='x', rotation=45)
        
        # Gene body CpG metrics
        sns.boxplot(data=high_quality_df,
                   x='binding_type',
                   y='gene_body_cpg_count',
                   hue='expression_status',
                   ax=ax3)
        ax3.set_title('Gene Body CpG Count')
        ax3.tick_params(axis='x', rotation=45)
        
        sns.boxplot(data=high_quality_df,
                   x='binding_type',
                   y='gene_body_cpg_density',
                   hue='expression_status',
                   ax=ax4)
        ax4.set_title('Gene Body CpG Density')
        ax4.tick_params(axis='x', rotation=45)
        
        plt.suptitle(f'{cell_type}: CpG Metrics by Region and Binding Type', y=1.02)
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, f'{cell_type}_cpg_metrics.pdf'))
        plt.close()
        
        # 4. Methylation vs CpG density scatter plots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        sns.scatterplot(data=high_quality_df,
                       x='promoter_cpg_density',
                       y='promoter_methylation',
                       hue='expression_status',
                       style='binding_type',
                       ax=ax1)
        ax1.set_title('Promoter Methylation vs CpG Density')
        
        sns.scatterplot(data=high_quality_df,
                       x='gene_body_cpg_density',
                       y='gene_body_methylation',
                       hue='expression_status',
                       style='binding_type',
                       ax=ax2)
        ax2.set_title('Gene Body Methylation vs CpG Density')
        
        plt.suptitle(f'{cell_type}: Methylation vs CpG Density', y=1.02)
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, f'{cell_type}_methylation_vs_density.pdf'))
        plt.close()

#Statistical analysis
def perform_statistical_analysis(results: Dict[str, pd.DataFrame]) -> Dict[str, Dict]:
    """Perform statistical tests on methylation patterns
    
    Args:
        results: Dictionary containing DataFrames with methylation data
        
    Returns:
        Dictionary containing statistical test results for each cell type and region
    """
    stats_results = {}
    
    for cell_type, df in results.items():
        # Add expression status based on log2FoldChange and padj
        df = df.copy()  # Create a copy to avoid modifying original
        df['expression_status'] = 'unchanged'
        df.loc[(df['log2FoldChange'] > 1) & (df['padj'] < 0.05), 'expression_status'] = 'upregulated'
        df.loc[(df['log2FoldChange'] < -1) & (df['padj'] < 0.05), 'expression_status'] = 'downregulated'
        
        # Filter out low quality data
        high_quality_df = df[
            ~df['low_cpg_promoter'] & 
            ~df['high_meth_promoter']
        ]
        
        cell_stats = {}
        
        # Log sample sizes
        logger.info(f"\nStatistical analysis for {cell_type}:")
        logger.info(f"Total samples: {len(df)}")
        logger.info(f"High quality samples: {len(high_quality_df)}")
        logger.info("Expression status distribution:")
        logger.info(high_quality_df['expression_status'].value_counts().to_string())
        
        # Analyze differences between up and down regulated genes
        for region in ['promoter', 'gene_body']:
            region_stats = {}
            
            # Get data for each comparison
            up_data = high_quality_df[high_quality_df['expression_status'] == 'upregulated']
            down_data = high_quality_df[high_quality_df['expression_status'] == 'downregulated']
            
            # Only perform tests if we have enough samples
            if len(up_data) >= 3 and len(down_data) >= 3:
                # Methylation differences
                try:
                    stat, pval = stats.mannwhitneyu(
                        up_data[f'{region}_methylation'].dropna(),
                        down_data[f'{region}_methylation'].dropna(),
                        alternative='two-sided'
                    )
                    region_stats['methylation_comparison'] = {
                        'statistic': stat,
                        'pvalue': pval,
                        'up_mean': up_data[f'{region}_methylation'].mean(),
                        'down_mean': down_data[f'{region}_methylation'].mean(),
                        'up_n': len(up_data),
                        'down_n': len(down_data)
                    }
                except Exception as e:
                    logger.warning(f"Could not perform methylation test for {region}: {str(e)}")
                    region_stats['methylation_comparison'] = None
                
                # CpG density differences
                try:
                    cpg_stat, cpg_pval = stats.mannwhitneyu(
                        up_data[f'{region}_cpg_density'].dropna(),
                        down_data[f'{region}_cpg_density'].dropna(),
                        alternative='two-sided'
                    )
                    region_stats['cpg_density_comparison'] = {
                        'statistic': cpg_stat,
                        'pvalue': cpg_pval,
                        'up_mean': up_data[f'{region}_cpg_density'].mean(),
                        'down_mean': down_data[f'{region}_cpg_density'].mean(),
                        'up_n': len(up_data),
                        'down_n': len(down_data)
                    }
                except Exception as e:
                    logger.warning(f"Could not perform CpG density test for {region}: {str(e)}")
                    region_stats['cpg_density_comparison'] = None
            else:
                logger.warning(f"Insufficient samples for {region} statistical tests")
                region_stats['methylation_comparison'] = None
                region_stats['cpg_density_comparison'] = None
            
            cell_stats[region] = region_stats
        
        stats_results[cell_type] = cell_stats
    
    return stats_results

#Create visualization functions
def create_methylation_patterns_plot(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Create detailed methylation pattern plots with CpG density consideration
    
    Args:
        df: DataFrame containing methylation data with CpG density metrics
        cell_type: Cell type being analyzed
        output_dir: Base output directory
    """
    plot_dir = os.path.join(output_dir, 'methylation_patterns')
    os.makedirs(plot_dir, exist_ok=True)
    
    # Filter out low quality data
    high_quality_df = df[
        ~df['low_cpg_promoter'] & 
        ~df['high_meth_promoter']
    ]
    
    # Create figure with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Methylation vs CpG density
    sns.scatterplot(data=high_quality_df, 
                   x='promoter_cpg_density',
                   y='promoter_methylation',
                   hue='expression_status',
                   style='binding_type',
                   ax=ax1)
    ax1.set_title('Promoter Methylation vs CpG Density')
    ax1.set_xlabel('CpG Density (CpGs/kb)')
    ax1.set_ylabel('Methylation (%)')
    
    sns.scatterplot(data=high_quality_df,
                   x='gene_body_cpg_density',
                   y='gene_body_methylation',
                   hue='expression_status',
                   style='binding_type',
                   ax=ax2)
    ax2.set_title('Gene Body Methylation vs CpG Density')
    ax2.set_xlabel('CpG Density (CpGs/kb)')
    ax2.set_ylabel('Methylation (%)')
    
    # 2. Methylation by binding type with size reflecting CpG count
    sns.scatterplot(data=high_quality_df,
                   x='binding_type',
                   y='promoter_methylation',
                   size='promoter_cpg_count',
                   hue='expression_status',
                   ax=ax3)
    ax3.set_title('Promoter Methylation by Binding Type')
    ax3.set_xlabel('Binding Type')
    ax3.set_ylabel('Methylation (%)')
    
    sns.scatterplot(data=high_quality_df,
                   x='binding_type',
                   y='gene_body_methylation',
                   size='gene_body_cpg_count',
                   hue='expression_status',
                   ax=ax4)
    ax4.set_title('Gene Body Methylation by Binding Type')
    ax4.set_xlabel('Binding Type')
    ax4.set_ylabel('Methylation (%)')
    
    # Add sample size annotations
    for ax in [ax3, ax4]:
        for i, binding_type in enumerate(high_quality_df['binding_type'].unique()):
            n = len(high_quality_df[high_quality_df['binding_type'] == binding_type])
            ax.text(i, ax.get_ylim()[0], f'n={n}', ha='center', va='bottom')
    
    plt.suptitle(f'{cell_type} - Methylation Patterns\n(High quality data only)', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, f'{cell_type}_methylation_patterns.pdf'))
    plt.close()
    
    # Log analysis summary
    logger.info(f"\nMethylation pattern analysis for {cell_type}:")
    logger.info(f"Total samples: {len(df)}")
    logger.info(f"High quality samples: {len(high_quality_df)}")
    logger.info(f"Average promoter CpG density: {high_quality_df['promoter_cpg_density'].mean():.1f}")
    logger.info(f"Average gene body CpG density: {high_quality_df['gene_body_cpg_density'].mean():.1f}")

def create_binding_correlation_plot(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Create correlation plots between binding signals and methylation
    
    Args:
        df: DataFrame with methylation and binding data
        cell_type: Cell type being analyzed
        output_dir: Directory to save plots
    """
    # Filter for high quality data
    high_quality_df = df[
        ~df['low_cpg_promoter'] & 
        ~df['high_meth_promoter']
    ]
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Methylation vs binding type
    sns.boxplot(data=high_quality_df,
               x='binding_type',
               y='promoter_methylation',
               hue='expression_status',
               ax=ax1)
    ax1.set_title('Promoter Methylation by Binding Type')
    ax1.set_xlabel('Binding Type')
    ax1.set_ylabel('Methylation (%)')
    ax1.tick_params(axis='x', rotation=45)
    
    # 2. Gene body methylation vs binding type
    sns.boxplot(data=high_quality_df,
               x='binding_type',
               y='gene_body_methylation',
               hue='expression_status',
               ax=ax2)
    ax2.set_title('Gene Body Methylation by Binding Type')
    ax2.set_xlabel('Binding Type')
    ax2.set_ylabel('Methylation (%)')
    ax2.tick_params(axis='x', rotation=45)
    
    # 3. CpG density vs binding type
    sns.boxplot(data=high_quality_df,
               x='binding_type',
               y='promoter_cpg_density',
               hue='expression_status',
               ax=ax3)
    ax3.set_title('Promoter CpG Density by Binding Type')
    ax3.set_xlabel('Binding Type')
    ax3.set_ylabel('CpG Density')
    ax3.tick_params(axis='x', rotation=45)
    
    # 4. Gene body CpG density vs binding type
    sns.boxplot(data=high_quality_df,
               x='binding_type',
               y='gene_body_cpg_density',
               hue='expression_status',
               ax=ax4)
    ax4.set_title('Gene Body CpG Density by Binding Type')
    ax4.set_xlabel('Binding Type')
    ax4.set_ylabel('CpG Density')
    ax4.tick_params(axis='x', rotation=45)
    
    plt.suptitle(f'{cell_type}: Binding Type Correlations', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'regulatory_patterns', 
                            f'{cell_type}_binding_correlations.pdf'))
    plt.close()

def create_methylation_distribution_plot(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Create distribution plots for methylation, binding signals and CpG density
    
    Args:
        df: DataFrame containing methylation data with CpG metrics
        cell_type: Cell type being analyzed
        output_dir: Base output directory
    """
    plot_dir = os.path.join(output_dir, 'methylation_distributions')
    os.makedirs(plot_dir, exist_ok=True)
    
    # Filter out low quality data
    high_quality_df = df[
        ~df['low_cpg_promoter'] & 
        ~df['high_meth_promoter']
    ]
    
    # Create figure with 6 subplots (3x2)
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(15, 18))
    
    # Function to safely create KDE plot with error handling
    def safe_kdeplot(data, x, hue, ax, title):
        try:
            if data[x].var() > 0:  # Check if there's variance in the data
                sns.kdeplot(data=data, x=x, hue=hue, ax=ax)
                ax.set_title(title)
                
                # Add summary statistics
                for status in data[hue].unique():
                    subset = data[data[hue] == status][x]
                    if len(subset) > 0:
                        mean_val = subset.mean()
                        median_val = subset.median()
                        ax.axvline(mean_val, color='gray', linestyle='--', alpha=0.5)
                        ax.text(mean_val, ax.get_ylim()[1], 
                               f'{status}\nMean: {mean_val:.1f}\nMedian: {median_val:.1f}',
                               rotation=90, va='top', ha='right', fontsize=8)
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
    
    # Create distribution plots
    safe_kdeplot(high_quality_df, 'promoter_methylation', 'expression_status', 
                 ax1, 'Promoter Methylation Distribution')
    safe_kdeplot(high_quality_df, 'gene_body_methylation', 'expression_status', 
                 ax2, 'Gene Body Methylation Distribution')
    
    safe_kdeplot(high_quality_df, 'promoter_cpg_density', 'expression_status',
                 ax3, 'Promoter CpG Density Distribution')
    safe_kdeplot(high_quality_df, 'gene_body_cpg_density', 'expression_status',
                 ax4, 'Gene Body CpG Density Distribution')
    
    safe_kdeplot(high_quality_df, 'exo_signal', 'expression_status',
                 ax5, 'Exo Signal Distribution')
    safe_kdeplot(high_quality_df, 'endo_signal', 'expression_status',
                 ax6, 'Endo Signal Distribution')
    
    plt.suptitle(f'{cell_type} - Feature Distributions\n(High quality data only)', 
                 y=1.02, size=16)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, f'{cell_type}_distributions.pdf'))
    plt.close()
    
    # Log distribution statistics
    logger.info(f"\nDistribution analysis for {cell_type}:")
    for status in high_quality_df['expression_status'].unique():
        subset = high_quality_df[high_quality_df['expression_status'] == status]
        logger.info(f"\n{status.upper()}:")
        for metric in ['promoter_methylation', 'gene_body_methylation',
                      'promoter_cpg_density', 'gene_body_cpg_density']:
            data = subset[metric].dropna()
            if len(data) > 0:
                logger.info(f"{metric}:")
                logger.info(f"  Mean ± SD: {data.mean():.1f} ± {data.std():.1f}")
                logger.info(f"  Median (Q1-Q3): {data.median():.1f} "
                          f"({data.quantile(0.25):.1f}-{data.quantile(0.75):.1f})")
                logger.info(f"  n = {len(data)}")

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
                    f.write("\nCorrelation Analysis:\n")
                    f.write(f"Promoter-Gene Body Correlation:\n")
                    f.write(f"rho: {corr_result['statistic']:.4f}\n")
                    f.write(f"p-value: {corr_result['pvalue']:.4f}\n")
                    f.write(f"n: {corr_result['n']}\n")

#Create focused visualization of significant findings
def plot_significant_differences(results: Dict[str, pd.DataFrame], output_dir: str):
    """Create plots showing significant differences in methylation patterns
    
    Args:
        results: Dictionary containing DataFrames with methylation data
        output_dir: Directory to save output plots
    """
    # Create output directory
    plot_dir = os.path.join(output_dir, 'significant_differences')
    os.makedirs(plot_dir, exist_ok=True)
    
    # Set up colors for plotting
    colors = {'upregulated': '#2ecc71', 'downregulated': '#e74c3c'}
    
    for cell_type, df in results.items():
        # Filter for MeCP2-bound genes (any binding type)
        mecp2_bound_df = df[df['binding_type'].notna()].copy()
        
        # Filter out low quality data
        high_quality_df = mecp2_bound_df[
            ~mecp2_bound_df['low_cpg_promoter'] & 
            ~mecp2_bound_df['high_meth_promoter']
        ]
        
        # Filter for regulated genes
        regulated_df = high_quality_df[
            high_quality_df['expression_status'].isin(['upregulated', 'downregulated'])
        ]
        
        # Log sample sizes
        logger.info(f"\nSignificant differences analysis for {cell_type}:")
        logger.info(f"Total MeCP2-bound genes: {len(mecp2_bound_df)}")
        logger.info(f"High quality genes: {len(high_quality_df)}")
        logger.info(f"Regulated genes: {len(regulated_df)}")
        
        if len(regulated_df) == 0:
            logger.warning(f"No regulated genes found for {cell_type}")
            continue
        
        # Create plots for each binding type
        for binding_type in ['exo', 'endo', 'both']:
            binding_df = regulated_df[regulated_df['binding_type'] == binding_type]
            if len(binding_df) == 0:
                continue
                
            # 1. Methylation patterns
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # Promoter methylation
            sns.boxplot(data=binding_df,
                       x='expression_status',
                       y='promoter_methylation',
                       palette=colors,
                       ax=ax1)
            ax1.set_title('Promoter Methylation')
            ax1.tick_params(axis='x', rotation=45)
            
            # Gene body methylation
            sns.boxplot(data=binding_df,
                       x='expression_status',
                       y='gene_body_methylation',
                       palette=colors,
                       ax=ax2)
            ax2.set_title('Gene Body Methylation')
            ax2.tick_params(axis='x', rotation=45)
            
            plt.suptitle(f'{cell_type}: {binding_type.title()} Binding Methylation Patterns', y=1.02)
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, f'{cell_type}_{binding_type}_methylation.pdf'))
            plt.close()
            
            # 2. CpG density patterns
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # Promoter CpG density
            sns.boxplot(data=binding_df,
                       x='expression_status',
                       y='promoter_cpg_density',
                       palette=colors,
                       ax=ax1)
            ax1.set_title('Promoter CpG Density')
            ax1.tick_params(axis='x', rotation=45)
            
            # Gene body CpG density
            sns.boxplot(data=binding_df,
                       x='expression_status',
                       y='gene_body_cpg_density',
                       palette=colors,
                       ax=ax2)
            ax2.set_title('Gene Body CpG Density')
            ax2.tick_params(axis='x', rotation=45)
            
            plt.suptitle(f'{cell_type}: {binding_type.title()} Binding CpG Density', y=1.02)
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, f'{cell_type}_{binding_type}_cpg_density.pdf'))
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

#Print summary for regulated genes
def print_regulated_summary(results: Dict[str, pd.DataFrame]):
    """Print summary statistics for regulated MeCP2-bound genes
    
    Args:
        results: Dictionary containing DataFrames with methylation data
    """
    for cell_type, cell_data in results.items():
        # Skip if cell_data is not a DataFrame
        if not isinstance(cell_data, pd.DataFrame):
            logger.warning(f"Skipping {cell_type}: data is not a DataFrame")
            continue
            
        print(f"\nSummary for {cell_type}:")
        print("="*50)
        
        # Add expression status if not present
        df = cell_data.copy()
        if 'expression_status' not in df.columns:
            df['expression_status'] = 'unchanged'
            df.loc[(df['log2FoldChange'] > 1) & (df['padj'] < 0.05), 'expression_status'] = 'upregulated'
            df.loc[(df['log2FoldChange'] < -1) & (df['padj'] < 0.05), 'expression_status'] = 'downregulated'
        
        # Filter for MeCP2-bound and regulated genes
        regulated_df = df[
            df['binding_type'].notna() & 
            (df['expression_status'].isin(['upregulated', 'downregulated']))
        ]
        
        # Filter out low quality data
        high_quality_df = regulated_df[
            ~regulated_df['low_cpg_promoter'] & 
            ~regulated_df['high_meth_promoter']
        ]
        
        # Print overall statistics
        print(f"\nTotal genes: {len(df)}")
        print(f"Regulated genes: {len(regulated_df)}")
        print(f"High quality regulated genes: {len(high_quality_df)}")
        
        # Separate by binding type
        for binding_type in ['exo', 'both', 'endo']:
            subset_df = high_quality_df[high_quality_df['binding_type'] == binding_type]
            
            if len(subset_df) == 0:
                continue
                
            print(f"\n{binding_type.upper()} BINDING:")
            print("-"*20)
            
            # Count genes by expression status
            print("\nGene counts by expression status:")
            print(subset_df['expression_status'].value_counts())
            
            # Mean methylation by expression status
            print("\nMean promoter methylation:")
            methylation_stats = subset_df.groupby('expression_status')['promoter_methylation'].agg(['mean', 'std', 'count'])
            print(methylation_stats.round(2))
            
            print("\nMean gene body methylation:")
            gene_body_stats = subset_df.groupby('expression_status')['gene_body_methylation'].agg(['mean', 'std', 'count'])
            print(gene_body_stats.round(2))
            
            # CpG density statistics
            print("\nPromoter CpG density:")
            cpg_stats = subset_df.groupby('expression_status')['promoter_cpg_density'].agg(['mean', 'std', 'count'])
            print(cpg_stats.round(2))
            
            print("\nGene body CpG density:")
            gene_body_cpg_stats = subset_df.groupby('expression_status')['gene_body_cpg_density'].agg(['mean', 'std', 'count'])
            print(gene_body_cpg_stats.round(2))

#Create detailed visualizations for MeCP2-bound regulated genes
def create_mecp2_regulated_analysis(results: Dict[str, pd.DataFrame], output_dir: str):
    """Create analysis of genes regulated by MeCP2 binding
    
    Args:
        results: Dictionary containing DataFrames with methylation data
        output_dir: Base output directory for saving results
    """
    # Create output directory
    analysis_dir = os.path.join(output_dir, 'mecp2_regulated_analysis')
    os.makedirs(analysis_dir, exist_ok=True)
    
    for cell_type, df in results.items():
        # Verify we have a DataFrame
        if not isinstance(df, pd.DataFrame):
            logger.error(f"Data for {cell_type} is not a DataFrame")
            continue
            
        logger.info(f"\nAnalyzing regulated genes for {cell_type}")
        
        # Add expression status if not present
        df = df.copy()
        if 'expression_status' not in df.columns:
            df['expression_status'] = 'unchanged'
            df.loc[(df['log2FoldChange'] > 1) & (df['padj'] < 0.05), 'expression_status'] = 'upregulated'
            df.loc[(df['log2FoldChange'] < -1) & (df['padj'] < 0.05), 'expression_status'] = 'downregulated'
        
        # Rest of the function remains the same...
        # Filter for MeCP2-bound (any binding type except NaN) and regulated genes
        mecp2_bound_df = df[
            df['binding_type'].notna() & 
            (df['expression_status'].isin(['upregulated', 'downregulated']))
        ]
        
        # Filter out low quality data
        high_quality_df = mecp2_bound_df[
            ~mecp2_bound_df['low_cpg_promoter'] & 
            ~mecp2_bound_df['high_meth_promoter']
        ]
        
        # Log sample sizes
        logger.info(f"\nMeCP2 regulated analysis for {cell_type}:")
        logger.info(f"Total bound and regulated genes: {len(mecp2_bound_df)}")
        logger.info(f"High quality genes: {len(high_quality_df)}")
        logger.info("\nBinding type distribution:")
        logger.info(high_quality_df['binding_type'].value_counts().to_string())
        logger.info("\nExpression status distribution:")
        logger.info(high_quality_df['expression_status'].value_counts().to_string())
        
        # Separate by binding type
        binding_groups = {
            'exo_enriched': high_quality_df[high_quality_df['binding_type'].isin(['exo', 'both'])],
            'exo_only': high_quality_df[high_quality_df['binding_type'] == 'exo'],
            'endo_only': high_quality_df[high_quality_df['binding_type'] == 'endo']
        }
        
        # Create plots
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Plot methylation levels by binding type and regulation
        sns.boxplot(data=high_quality_df,
                   x='binding_type',
                   y='promoter_methylation',
                   hue='expression_status',
                   ax=ax1)
        ax1.set_title('Promoter Methylation')
        ax1.tick_params(axis='x', rotation=45)
        
        sns.boxplot(data=high_quality_df,
                   x='binding_type',
                   y='gene_body_methylation',
                   hue='expression_status',
                   ax=ax2)
        ax2.set_title('Gene Body Methylation')
        ax2.tick_params(axis='x', rotation=45)
        
        # Plot CpG density
        sns.boxplot(data=high_quality_df,
                   x='binding_type',
                   y='promoter_cpg_density',
                   hue='expression_status',
                   ax=ax3)
        ax3.set_title('Promoter CpG Density')
        ax3.tick_params(axis='x', rotation=45)
        
        sns.boxplot(data=high_quality_df,
                   x='binding_type',
                   y='gene_body_cpg_density',
                   hue='expression_status',
                   ax=ax4)
        ax4.set_title('Gene Body CpG Density')
        ax4.tick_params(axis='x', rotation=45)
        
        plt.suptitle(f'{cell_type}: MeCP2-Regulated Genes Analysis', y=1.02)
        plt.tight_layout()
        plt.savefig(os.path.join(analysis_dir, f'{cell_type}_regulated_analysis.pdf'))
        plt.close()
        
        # Save gene lists and statistics
        stats_file = os.path.join(analysis_dir, f'{cell_type}_regulated_stats.txt')
        with open(stats_file, 'w') as f:
            f.write(f"MeCP2 Regulated Genes Analysis - {cell_type}\n")
            f.write("="*50 + "\n\n")
            
            for binding_type, group_df in binding_groups.items():
                if len(group_df) == 0:
                    continue
                    
                f.write(f"\n{binding_type.upper()}\n")
                f.write("-"*30 + "\n")
                f.write(f"Total genes: {len(group_df)}\n")
                
                # Expression status breakdown
                f.write("\nExpression status:\n")
                expr_counts = group_df['expression_status'].value_counts()
                f.write(expr_counts.to_string())
                f.write("\n")
                
                # Save gene lists
                output_file = os.path.join(analysis_dir, 
                                         f'{cell_type}_{binding_type}_regulated_genes.csv')
                group_df.to_csv(output_file, index=False)
                
                # Calculate and write statistics
                for region in ['promoter', 'gene_body']:
                    f.write(f"\n{region.title()} statistics:\n")
                    
                    # Methylation stats
                    f.write("  Methylation:\n")
                    for status in ['upregulated', 'downregulated']:
                        status_data = group_df[group_df['expression_status'] == status]
                        if len(status_data) > 0:
                            f.write(f"    {status}:\n")
                            f.write(f"      Mean: {status_data[f'{region}_methylation'].mean():.2f}\n")
                            f.write(f"      Std: {status_data[f'{region}_methylation'].std():.2f}\n")
                            f.write(f"      n: {len(status_data)}\n")
                    
                    # CpG density stats
                    f.write("\n  CpG Density:\n")
                    for status in ['upregulated', 'downregulated']:
                        status_data = group_df[group_df['expression_status'] == status]
                        if len(status_data) > 0:
                            f.write(f"    {status}:\n")
                            f.write(f"      Mean: {status_data[f'{region}_cpg_density'].mean():.2f}\n")
                            f.write(f"      Std: {status_data[f'{region}_cpg_density'].std():.2f}\n")
                            f.write(f"      n: {len(status_data)}\n")

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

def analyze_single_cell_binding(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze binding enrichment for a single cell type"""
    
    # Create exo-enriched group (both exo-only and exo-enriched)
    exo_enriched = df[df['binding_type'].isin(['exo', 'both'])]
    endo_only = df[df['binding_type'] == 'endo']
    non_enriched = df[~df['mecp2_bound']]
    
    # Analyze each group
    groups = {
        'exo_enriched': exo_enriched,
        'endo_only': endo_only,
        'non_enriched': non_enriched
    }
    
    group_results = {}
    
    for group_name, group_df in groups.items():
        if len(group_df) == 0:
            logger.warning(f"No genes found in {group_name} group for {cell_type}")
            continue
            
        # Get regulation counts
        regulation_counts = group_df['expression_status'].value_counts()
        
        # Calculate methylation statistics for each regulation status
        stats = {}
        for status in ['not_deregulated', 'upregulated', 'downregulated']:
            subset = group_df[group_df['expression_status'] == status]
            if len(subset) > 0:
                stats[status] = {
                    'count': len(subset),
                    'promoter_methylation': {
                        'mean': subset['promoter_methylation'].mean(),
                        'std': subset['promoter_methylation'].std(),
                        'median': subset['promoter_methylation'].median()
                    },
                    'gene_body_methylation': {
                        'mean': subset['gene_body_methylation'].mean(),
                        'std': subset['gene_body_methylation'].std(),
                        'median': subset['gene_body_methylation'].median()
                    }
                }
        
        # Save results
        output_file = os.path.join(output_dir, f'{group_name}_analysis.txt')
        with open(output_file, 'w') as f:
            f.write(f"Analysis for {group_name} - {cell_type}\n")
            f.write("="*50 + "\n\n")
            
            f.write("Regulation counts:\n")
            f.write(str(regulation_counts) + "\n\n")
            
            for status, stat in stats.items():
                f.write(f"\n{status.upper()}:\n")
                f.write(f"Count: {stat['count']}\n")
                f.write("\nPromoter methylation:\n")
                f.write(f"Mean ± SD: {stat['promoter_methylation']['mean']:.2f} ± {stat['promoter_methylation']['std']:.2f}\n")
                f.write(f"Median: {stat['promoter_methylation']['median']:.2f}\n")
                f.write("\nGene body methylation:\n")
                f.write(f"Mean ± SD: {stat['gene_body_methylation']['mean']:.2f} ± {stat['gene_body_methylation']['std']:.2f}\n")
                f.write(f"Median: {stat['gene_body_methylation']['median']:.2f}\n")
        
        # Save gene lists
        for status in ['not_deregulated', 'upregulated', 'downregulated']:
            subset = group_df[group_df['expression_status'] == status]
            if len(subset) > 0:
                output_file = os.path.join(output_dir, f'{group_name}_{status}_genes.csv')
                subset.to_csv(output_file, index=False)
        
        group_results[group_name] = {
            'regulation_counts': regulation_counts,
            'methylation_stats': pd.DataFrame(stats).T
        }
    
    return group_results

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

def create_methylation_comparison_plot(results: dict, cell_type: str, output_dir: str):
    """Create methylation comparison plots for different binding categories with statistics
    
    This function generates two side-by-side bar plots comparing methylation levels:
    1. Promoter methylation across binding categories and regulation status
    2. Gene body methylation across binding categories and regulation status
    
    It also performs statistical testing between upregulated and downregulated genes
    within each binding category and annotates significant differences.
    
    Args:
        results: Dictionary containing methylation statistics for each binding category and regulation status
        cell_type: String indicating cell type being analyzed (e.g. 'NSC', 'NEU')
        output_dir: Directory path where plots should be saved
    """
    
    # Initialize dictionary to store plot data in a format suitable for pandas DataFrame
    plot_data = {
        'Category': [],      # Binding category (e.g. exo_enriched, endo_only)
        'Regulation': [],    # Regulation status (up/down/unchanged)
        'Region': [],        # Genomic region (Promoter or Gene Body)
        'Methylation': [],   # Mean methylation level
        'StdDev': []        # Standard deviation of methylation
    }
    
    # Extract methylation data for each category, regulation status and genomic region
    for category, stats in results.items():
        for reg_status, reg_stats in stats.items():
            # Add promoter methylation data
            plot_data['Category'].append(category)
            plot_data['Regulation'].append(reg_status)
            plot_data['Region'].append('Promoter')
            plot_data['Methylation'].append(
                reg_stats['promoter_metrics']['methylation']['mean']
            )
            plot_data['StdDev'].append(
                reg_stats['promoter_metrics']['methylation']['std']
            )
            
            # Add gene body methylation data
            plot_data['Category'].append(category)
            plot_data['Regulation'].append(reg_status)
            plot_data['Region'].append('Gene Body')
            plot_data['Methylation'].append(
                reg_stats['gene_body_metrics']['methylation']['mean']
            )
            plot_data['StdDev'].append(
                reg_stats['gene_body_metrics']['methylation']['std']
            )
    
    # Convert to DataFrame for easier plotting
    df = pd.DataFrame(plot_data)
    
    # Create figure with two subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # Create promoter methylation plot
    promoter_data = df[df['Region'] == 'Promoter']
    sns.barplot(data=promoter_data, x='Category', y='Methylation',
                hue='Regulation', ax=ax1)
    ax1.set_title(f'{cell_type} - Promoter Methylation')
    ax1.set_xlabel('Binding Category')
    ax1.set_ylabel('Methylation Level (%)')
    ax1.tick_params(axis='x', rotation=45)
    
    # Create gene body methylation plot
    gene_body_data = df[df['Region'] == 'Gene Body']
    sns.barplot(data=gene_body_data, x='Category', y='Methylation',
                hue='Regulation', ax=ax2)
    ax2.set_title(f'{cell_type} - Gene Body Methylation')
    ax2.set_xlabel('Binding Category')
    ax2.set_ylabel('Methylation Level (%)')
    ax2.tick_params(axis='x', rotation=45)
    
    # Add statistical significance annotations
    for ax, data in [(ax1, promoter_data), (ax2, gene_body_data)]:
        for i, category in enumerate(data['Category'].unique()):
            category_data = data[data['Category'] == category]
            
            # Compare methylation between upregulated and downregulated genes
            up_data = category_data[category_data['Regulation'] == 'upregulated']
            down_data = category_data[category_data['Regulation'] == 'downregulated']
            
            # Only perform statistical test if we have data for both groups
            if len(up_data) > 0 and len(down_data) > 0:
                stats = calculate_group_significance(
                    up_data['Methylation'],
                    down_data['Methylation'],
                    'upregulated',
                    'downregulated'
                )
                
                # Add p-value annotation if difference is significant
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

def analyze_binding_enrichment_groups(df: pd.DataFrame, output_dir: str) -> Dict[str, Dict]:
    """Analyze binding enrichment patterns for different groups
    
    Args:
        df: DataFrame containing methylation and binding data
        output_dir: Directory to save output files
        
    Returns:
        Dictionary containing analysis results for each binding group:
        {
            'exo_enriched': {
                'methylation_stats': {...},
                'cpg_density_stats': {...}
            },
            'endo_only': {...},
            'unbound': {...}
        }
    """
    # Define binding groups
    binding_groups = {
        'exo_enriched': df[df['binding_type'].isin(['exo', 'both'])],
        'endo_only': df[df['binding_type'] == 'endo'],
        'unbound': df[df['binding_type'].isna()]  # Changed from mecp2_bound to binding_type
    }
    
    results = {}
    
    for group_name, group_df in binding_groups.items():
        # Skip if group is empty
        if len(group_df) == 0:
            logger.warning(f"No genes found in {group_name} group")
            continue
            
        # Filter for high quality data
        high_quality_df = group_df[
            ~group_df['low_cpg_promoter'] & 
            ~group_df['high_meth_promoter']
        ]
        
        if len(high_quality_df) == 0:
            logger.warning(f"No high quality genes found in {group_name} group")
            continue
        
        # Calculate statistics
        results[group_name] = {
            'methylation_stats': {
                'promoter': {
                    'mean': high_quality_df['promoter_methylation'].mean(),
                    'std': high_quality_df['promoter_methylation'].std(),
                    'n': len(high_quality_df)
                },
                'gene_body': {
                    'mean': high_quality_df['gene_body_methylation'].mean(),
                    'std': high_quality_df['gene_body_methylation'].std(),
                    'n': len(high_quality_df)
                }
            },
            'cpg_density_stats': {
                'promoter': {
                    'mean': high_quality_df['promoter_cpg_density'].mean(),
                    'std': high_quality_df['promoter_cpg_density'].std(),
                    'n': len(high_quality_df)
                },
                'gene_body': {
                    'mean': high_quality_df['gene_body_cpg_density'].mean(),
                    'std': high_quality_df['gene_body_cpg_density'].std(),
                    'n': len(high_quality_df)
                }
            }
        }
        
        # Create plots
        # create_binding_group_plots(high_quality_df, group_name, output_dir)
    
    return results

def analyze_exo_upregulated(results: Dict[str, pd.DataFrame], output_dir: str):
    """Detailed analysis of upregulated genes with exo MeCP2 binding
    
    This function:
    1. Filters for genes that are upregulated and have exo or both binding
    2. Saves detailed results to CSV
    3. Creates visualizations of the patterns
    
    Args:
        results: Dictionary of DataFrames with methylation results by cell type
        output_dir: Directory to save output files
    """
    for cell_type, df in results.items():
        # Filter for exo-bound and upregulated genes
        exo_up_df = df[
            (df['binding_type'].isin(['exo', 'both'])) & 
            (df['expression_status'] == 'upregulated')
        ]
        
        if len(exo_up_df) == 0:
            logger.warning(f"No exo-bound upregulated genes found for {cell_type}")
            continue
        
        # Save detailed results
        output_file = os.path.join(output_dir, 'exo_upregulated', 
                                 f'{cell_type}_exo_upregulated_analysis.csv')
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        exo_up_df.to_csv(output_file, index=False)
        
        # Create visualizations
        create_exo_up_visualizations(exo_up_df, cell_type, output_dir)

def create_exo_enriched_analysis_plots(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Create detailed plots specifically for exo-enriched samples
    
    Args:
        df: DataFrame with methylation and binding data
        cell_type: Cell type being analyzed
        output_dir: Directory to save plots
    """
    # Create output directory
    plot_dir = os.path.join(output_dir, 'exo_enriched_analysis')
    os.makedirs(plot_dir, exist_ok=True)
    
    # Filter for exo-enriched samples (exo or both binding)
    exo_df = df[df['binding_type'].isin(['exo', 'both'])]
    
    if len(exo_df) == 0:
        logger.warning(f"No exo-enriched samples found for {cell_type}")
        return
    
    # Filter for high quality data
    high_quality_df = exo_df[
        ~exo_df['low_cpg_promoter'] & 
        ~exo_df['high_meth_promoter']
    ]
    
    if len(high_quality_df) == 0:
        logger.warning(f"No high quality exo-enriched samples found for {cell_type}")
        return
    
    # Create plots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Methylation patterns
    sns.boxplot(data=high_quality_df,
               x='binding_type',
               y='promoter_methylation',
               hue='expression_status',
               ax=ax1)
    ax1.set_title('Promoter Methylation')
    ax1.tick_params(axis='x', rotation=45)
    
    sns.boxplot(data=high_quality_df,
               x='binding_type',
               y='gene_body_methylation',
               hue='expression_status',
               ax=ax2)
    ax2.set_title('Gene Body Methylation')
    ax2.tick_params(axis='x', rotation=45)
    
    # 2. CpG density patterns
    sns.boxplot(data=high_quality_df,
               x='binding_type',
               y='promoter_cpg_density',
               hue='expression_status',
               ax=ax3)
    ax3.set_title('Promoter CpG Density')
    ax3.tick_params(axis='x', rotation=45)
    
    sns.boxplot(data=high_quality_df,
               x='binding_type',
               y='gene_body_cpg_density',
               hue='expression_status',
               ax=ax4)
    ax4.set_title('Gene Body CpG Density')
    ax4.tick_params(axis='x', rotation=45)
    
    plt.suptitle(f'{cell_type}: Exo-Enriched Gene Analysis', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, f'{cell_type}_exo_enriched_patterns.pdf'))
    plt.close()
    
    # Save statistics
    stats_file = os.path.join(plot_dir, f'{cell_type}_exo_enriched_stats.txt')
    with open(stats_file, 'w') as f:
        f.write(f"Exo-Enriched Analysis - {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        f.write("Sample sizes:\n")
        f.write(f"Total exo-enriched genes: {len(exo_df)}\n")
        f.write(f"High quality genes: {len(high_quality_df)}\n\n")
        
        for region in ['promoter', 'gene_body']:
            f.write(f"\n{region.title()} Statistics:\n")
            f.write("-"*30 + "\n")
            
            # Methylation stats
            f.write("\nMethylation:\n")
            for status in ['upregulated', 'downregulated', 'unchanged']:
                status_data = high_quality_df[high_quality_df['expression_status'] == status]
                if len(status_data) > 0:
                    f.write(f"  {status}:\n")
                    f.write(f"    Mean: {status_data[f'{region}_methylation'].mean():.2f}\n")
                    f.write(f"    Std: {status_data[f'{region}_methylation'].std():.2f}\n")
                    f.write(f"    n: {len(status_data)}\n")
            
            # CpG density stats
            f.write("\nCpG Density:\n")
            for status in ['upregulated', 'downregulated', 'unchanged']:
                status_data = high_quality_df[high_quality_df['expression_status'] == status]
                if len(status_data) > 0:
                    f.write(f"  {status}:\n")
                    f.write(f"    Mean: {status_data[f'{region}_cpg_density'].mean():.2f}\n")
                    f.write(f"    Std: {status_data[f'{region}_cpg_density'].std():.2f}\n")
                    f.write(f"    n: {len(status_data)}\n")

def analyze_exo_enriched_patterns(results: Dict[str, pd.DataFrame], output_dir: str, 
                                exo_enrichment_threshold: float = 1.5):
    """Analyze patterns in exo-enriched samples across all cell types
    
    Args:
        results: Dictionary of DataFrames with methylation results by cell type
        output_dir: Directory to save output files
        exo_enrichment_threshold: Minimum ratio of exo/endo signal to consider as enriched (default: 1.5)
    """
    for cell_type, df in results.items():
        logger.info(f"Analyzing exo-enriched patterns for {cell_type}")
        create_exo_enriched_analysis_plots(df, cell_type, output_dir)

def create_mecp2_regulated_analysis(results: Dict[str, pd.DataFrame], output_dir: str):
    """Create analysis of genes regulated by MeCP2 binding
    
    Args:
        results: Dictionary containing DataFrames with methylation data
        output_dir: Base output directory for saving results
    """
    # Create output directory
    analysis_dir = os.path.join(output_dir, 'mecp2_regulated_analysis')
    os.makedirs(analysis_dir, exist_ok=True)
    
    for cell_type, df in results.items():
        # Verify we have a DataFrame
        if not isinstance(df, pd.DataFrame):
            logger.error(f"Data for {cell_type} is not a DataFrame")
            continue
            
        logger.info(f"\nAnalyzing regulated genes for {cell_type}")
        
        # Add expression status if not present
        df = df.copy()
        if 'expression_status' not in df.columns:
            df['expression_status'] = 'unchanged'
            df.loc[(df['log2FoldChange'] > 1) & (df['padj'] < 0.05), 'expression_status'] = 'upregulated'
            df.loc[(df['log2FoldChange'] < -1) & (df['padj'] < 0.05), 'expression_status'] = 'downregulated'
        
        # Rest of the function remains the same...
        # Filter for MeCP2-bound (any binding type except NaN) and regulated genes
        mecp2_bound_df = df[
            df['binding_type'].notna() & 
            (df['expression_status'].isin(['upregulated', 'downregulated']))
        ]
        
        # Filter out low quality data
        high_quality_df = mecp2_bound_df[
            ~mecp2_bound_df['low_cpg_promoter'] & 
            ~mecp2_bound_df['high_meth_promoter']
        ]
        
        # Log sample sizes
        logger.info(f"\nMeCP2 regulated analysis for {cell_type}:")
        logger.info(f"Total bound and regulated genes: {len(mecp2_bound_df)}")
        logger.info(f"High quality genes: {len(high_quality_df)}")
        logger.info("\nBinding type distribution:")
        logger.info(high_quality_df['binding_type'].value_counts().to_string())
        logger.info("\nExpression status distribution:")
        logger.info(high_quality_df['expression_status'].value_counts().to_string())
        
        # Separate by binding type
        binding_groups = {
            'exo_enriched': high_quality_df[high_quality_df['binding_type'].isin(['exo', 'both'])],
            'exo_only': high_quality_df[high_quality_df['binding_type'] == 'exo'],
            'endo_only': high_quality_df[high_quality_df['binding_type'] == 'endo']
        }
        
        # Create plots
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Plot methylation levels by binding type and regulation
        sns.boxplot(data=high_quality_df,
                   x='binding_type',
                   y='promoter_methylation',
                   hue='expression_status',
                   ax=ax1)
        ax1.set_title('Promoter Methylation')
        ax1.tick_params(axis='x', rotation=45)
        
        sns.boxplot(data=high_quality_df,
                   x='binding_type',
                   y='gene_body_methylation',
                   hue='expression_status',
                   ax=ax2)
        ax2.set_title('Gene Body Methylation')
        ax2.tick_params(axis='x', rotation=45)
        
        # Plot CpG density
        sns.boxplot(data=high_quality_df,
                   x='binding_type',
                   y='promoter_cpg_density',
                   hue='expression_status',
                   ax=ax3)
        ax3.set_title('Promoter CpG Density')
        ax3.tick_params(axis='x', rotation=45)
        
        sns.boxplot(data=high_quality_df,
                   x='binding_type',
                   y='gene_body_cpg_density',
                   hue='expression_status',
                   ax=ax4)
        ax4.set_title('Gene Body CpG Density')
        ax4.tick_params(axis='x', rotation=45)
        
        plt.suptitle(f'{cell_type}: MeCP2-Regulated Genes Analysis', y=1.02)
        plt.tight_layout()
        plt.savefig(os.path.join(analysis_dir, f'{cell_type}_regulated_analysis.pdf'))
        plt.close()
        
        # Save gene lists and statistics
        stats_file = os.path.join(analysis_dir, f'{cell_type}_regulated_stats.txt')
        with open(stats_file, 'w') as f:
            f.write(f"MeCP2 Regulated Genes Analysis - {cell_type}\n")
            f.write("="*50 + "\n\n")
            
            for binding_type, group_df in binding_groups.items():
                if len(group_df) == 0:
                    continue
                    
                f.write(f"\n{binding_type.upper()}\n")
                f.write("-"*30 + "\n")
                f.write(f"Total genes: {len(group_df)}\n")
                
                # Expression status breakdown
                f.write("\nExpression status:\n")
                expr_counts = group_df['expression_status'].value_counts()
                f.write(expr_counts.to_string())
                f.write("\n")
                
                # Save gene lists
                output_file = os.path.join(analysis_dir, 
                                         f'{cell_type}_{binding_type}_regulated_genes.csv')
                group_df.to_csv(output_file, index=False)
                
                # Calculate and write statistics
                for region in ['promoter', 'gene_body']:
                    f.write(f"\n{region.title()} statistics:\n")
                    
                    # Methylation stats
                    f.write("  Methylation:\n")
                    for status in ['upregulated', 'downregulated']:
                        status_data = group_df[group_df['expression_status'] == status]
                        if len(status_data) > 0:
                            f.write(f"    {status}:\n")
                            f.write(f"      Mean: {status_data[f'{region}_methylation'].mean():.2f}\n")
                            f.write(f"      Std: {status_data[f'{region}_methylation'].std():.2f}\n")
                            f.write(f"      n: {len(status_data)}\n")
                    
                    # CpG density stats
                    f.write("\n  CpG Density:\n")
                    for status in ['upregulated', 'downregulated']:
                        status_data = group_df[group_df['expression_status'] == status]
                        if len(status_data) > 0:
                            f.write(f"    {status}:\n")
                            f.write(f"      Mean: {status_data[f'{region}_cpg_density'].mean():.2f}\n")
                            f.write(f"      Std: {status_data[f'{region}_cpg_density'].std():.2f}\n")
                            f.write(f"      n: {len(status_data)}\n")

def get_methylation_data():
    """Get methylation data from cache
    
    Returns:
        Dictionary containing DataFrames with methylation data
    """
    data = load_from_cache('methylation_results')
    if data is None:
        logger.error("Could not load methylation results from cache")
        return None
        
    # Verify data structure
    if not isinstance(data, dict):
        logger.error("Cached data is not a dictionary")
        return None
        
    # Verify each value is a DataFrame
    for cell_type, df in data.items():
        if not isinstance(df, pd.DataFrame):
            logger.error(f"Data for {cell_type} is not a DataFrame")
            return None
            
    return data

def analyze_binding_enrichment_detailed(results: Dict[str, pd.DataFrame], output_dir: str):
    """Perform detailed analysis of binding enrichment patterns
    
    Args:
        results: Dictionary containing DataFrames with methylation data
        output_dir: Directory to save output files
    """
    # Create output directory
    analysis_dir = os.path.join(output_dir, 'binding_enrichment')
    os.makedirs(analysis_dir, exist_ok=True)
    
    for cell_type, df in results.items():
        # Create binding groups
        binding_groups = {
            'exo_enriched': df[df['binding_type'].isin(['exo', 'both'])],
            'exo_only': df[df['binding_type'] == 'exo'],
            'endo_only': df[df['binding_type'] == 'endo'],
            'unbound': df[df['binding_type'].isna()]  # Changed from mecp2_bound
        }
        
        # Analyze each group
        for group_name, group_df in binding_groups.items():
            if len(group_df) == 0:
                logger.warning(f"No genes found in {group_name} group for {cell_type}")
                continue
                
            # Filter for high quality data
            high_quality_df = group_df[
                ~group_df['low_cpg_promoter'] & 
                ~group_df['high_meth_promoter']
            ]
            
            if len(high_quality_df) == 0:
                logger.warning(f"No high quality genes found in {group_name} group for {cell_type}")
                continue
            
            # Calculate statistics
            stats = {
                'total_genes': len(group_df),
                'high_quality_genes': len(high_quality_df),
                'methylation': {
                    'promoter': {
                        'mean': high_quality_df['promoter_methylation'].mean(),
                        'std': high_quality_df['promoter_methylation'].std()
                    },
                    'gene_body': {
                        'mean': high_quality_df['gene_body_methylation'].mean(),
                        'std': high_quality_df['gene_body_methylation'].std()
                    }
                },
                'cpg_density': {
                    'promoter': {
                        'mean': high_quality_df['promoter_cpg_density'].mean(),
                        'std': high_quality_df['promoter_cpg_density'].std()
                    },
                    'gene_body': {
                        'mean': high_quality_df['gene_body_cpg_density'].mean(),
                        'std': high_quality_df['gene_body_cpg_density'].std()
                    }
                }
            }
            
            # Save results
            output_file = os.path.join(analysis_dir, f'{cell_type}_{group_name}_analysis.txt')
            with open(output_file, 'w') as f:
                f.write(f"Binding Enrichment Analysis - {cell_type} - {group_name}\n")
                f.write("="*50 + "\n\n")
                
                f.write("Sample sizes:\n")
                f.write(f"Total genes: {stats['total_genes']}\n")
                f.write(f"High quality genes: {stats['high_quality_genes']}\n\n")
                
                for metric in ['methylation', 'cpg_density']:
                    f.write(f"\n{metric.title()} Statistics:\n")
                    f.write("-"*30 + "\n")
                    for region in ['promoter', 'gene_body']:
                        f.write(f"\n{region.title()}:\n")
                        f.write(f"  Mean: {stats[metric][region]['mean']:.2f}\n")
                        f.write(f"  Std: {stats[metric][region]['std']:.2f}\n")