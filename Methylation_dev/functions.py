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
from typing import Dict, List, Tuple, Any, Optional, Union
import pyranges as pr
from functools import partial
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
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

#%% Helper functions for genomic regions
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

#%% Analyze methylation patterns
def analyze_methylation_patterns(args):
    """
    Main analysis function
    """
    # Load binding data
    binding_df = load_binding_data(PATHS['binding_data'])
    if binding_df is None:
        logger.error("Failed to load binding data")
        return {}
    
    results = {}
    for cell_type in args.cell_types:
        try:
            # Load required data
            genes_df = load_gene_data(cell_type)
            expr_df = load_expression_data(cell_type) 
            cpg_associations = load_cpg_associations(cell_type)
            
            # Check if all data loaded successfully
            if any(x is None for x in [genes_df, expr_df, cpg_associations]):
                logger.error(f"Failed to load required data for {cell_type}")
                continue
            
            # Merge data
            merged_df = merge_analysis_data(genes_df, expr_df, cpg_associations)
            
            if merged_df is None or len(merged_df) == 0:
                logger.warning(f"No valid data for {cell_type}")
                continue
        
            # Assign binding types
            merged_df = assign_binding_types(merged_df, binding_df)
            
            # Store results
            results[cell_type] = merged_df
            
        except Exception as e:
            logger.error(f"Error analyzing {cell_type}: {str(e)}")
            continue
            
    return results

def assign_binding_types(genes_df: pd.DataFrame, binding_df: pd.DataFrame) -> pd.DataFrame:
    """
    Assign binding types to genes based on overlapping binding regions
    """
    if binding_df is None:
        logger.warning("No binding data provided - setting all binding types to 'none'")
        genes_df['binding_type'] = 'none'
        genes_df['exo_signal'] = None
        genes_df['endo_signal'] = None
        genes_df['mecp2_bound'] = False
        return genes_df

    # Set index for faster lookups
    genes_df = genes_df.set_index(['chromosome', 'start', 'end'])

    # Convert to PyRanges for efficient overlap detection
    genes_pr = pr.PyRanges(
        chromosomes=genes_df.reset_index()['chromosome'],
        starts=genes_df.reset_index()['start'],
        ends=genes_df.reset_index()['end'],
        strands=genes_df['strand']
    )
    
    binding_pr = pr.PyRanges(
        chromosomes=binding_df['chr'],
        starts=binding_df['start'],
        ends=binding_df['end'],
        strands=None
    )

    # Find overlaps
    overlaps = genes_pr.join(binding_pr).as_df()
    
    # Initialize columns
    genes_df['binding_type'] = 'none'
    genes_df['exo_signal'] = 0.0
    genes_df['endo_signal'] = 0.0
    genes_df['mecp2_bound'] = False

    # Process overlaps more efficiently
    for _, overlap in overlaps.iterrows():
        idx = (overlap['Chromosome'], overlap['Start'], overlap['End'])
        
        current_type = genes_df.at[idx, 'binding_type']
        new_type = overlap['binding_type']
        
        # Update binding type
        if current_type == 'none':
            genes_df.at[idx, 'binding_type'] = new_type
        elif current_type != new_type and new_type != 'none':
            genes_df.at[idx, 'binding_type'] = 'both'
        
        # Update signals (take maximum value)
        genes_df.at[idx, 'exo_signal'] = max(
            genes_df.at[idx, 'exo_signal'],
            overlap['exo_signal']
        )
        genes_df.at[idx, 'endo_signal'] = max(
            genes_df.at[idx, 'endo_signal'],
            overlap['endo_signal']
        )
        
        genes_df.at[idx, 'mecp2_bound'] = True

    # Reset index for consistency with rest of pipeline
    genes_df = genes_df.reset_index()

    # Log results
    logger.info("\nBinding type distribution after assignment:")
    logger.info(genes_df['binding_type'].value_counts())
    
    return genes_df

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

def calculate_normalized_methylation(ip_values, input_values, sequence):
    """Calculate biologically meaningful methylation levels"""
    if not ip_values or not input_values:
        return 0
        
    # Count CpGs in region
    cpg_count = sequence.upper().count('CG')
    if cpg_count == 0:
        return 0
        
    # Calculate average signals with proper normalization
    ip_mean = np.mean([x for x in ip_values if x is not None])
    input_mean = np.mean([x for x in input_values if x is not None])
    
    if input_mean <= 0:
        return 0
    
    # Improved biological normalization:
    # 1. CpG density normalization
    region_length = len(sequence)
    cpg_density = cpg_count / region_length
    
    # 2. Calculate enrichment with local background correction
    local_background = np.percentile([x for x in input_values if x is not None], 25)
    enrichment = (ip_mean - local_background) / (input_mean * cpg_density)
    
    # 3. Convert to methylation percentage using sigmoid transformation
    # This better reflects biological methylation levels
    methylation = 100 / (1 + np.exp(-enrichment))
    
    return max(0, min(100, methylation))

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
def analyze_methylation_distribution(df, cell_type, output_dir):
    """Analyze methylation distribution patterns"""
    try:
        if not isinstance(df, pd.DataFrame):
            logger.error(f"Expected DataFrame, got {type(df)}")
            return None
            
        if len(df) == 0:
            logger.warning(f"No data to analyze for {cell_type}")
            return None
            
        # Create output directory
        plot_dir = os.path.join(output_dir, 'methylation_analysis', cell_type)
        os.makedirs(plot_dir, exist_ok=True)
        
        # Create violin plots for methylation distribution
        plt.figure(figsize=(12, 6))
        
        # Plot methylation levels for promoters and gene bodies
        data_to_plot = pd.melt(
            df,
                              value_vars=['promoter_methylation', 'gene_body_methylation'],
            var_name='Region',
            value_name='Methylation Level'
        )
        
        sns.violinplot(data=data_to_plot, x='Region', y='Methylation Level')
        plt.title(f'{cell_type}: Distribution of Methylation Levels')
        plt.savefig(os.path.join(plot_dir, 'methylation_distribution.pdf'))
        plt.close()
        
        # Calculate and return statistics
        return {
            'promoter': df['promoter_methylation'].describe(),
            'gene_body': df['gene_body_methylation'].describe(),
            'correlation': calculate_methylation_correlation(df)
        }
        
    except Exception as e:
        logger.error(f"Error in methylation distribution analysis: {str(e)}")
        return None

def calculate_methylation_correlation(df):
    """Calculate correlation between promoter and gene body methylation"""
        try:
            valid_data = df[['promoter_methylation', 'gene_body_methylation']].dropna()
            if len(valid_data) > 1:
            correlation = stats.spearmanr(
                valid_data['promoter_methylation'],
                valid_data['gene_body_methylation']
            )
            return {
                    'rho': correlation.statistic,
                    'pvalue': correlation.pvalue,
                    'n': len(valid_data)
                }
        except Exception as e:
        logger.error(f"Error calculating correlation: {str(e)}")
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

def analyze_binding_patterns(df, cell_type, output_dir):
    """Analyze binding patterns"""
    try:
        if not isinstance(df, pd.DataFrame):
            logger.error(f"Expected DataFrame, got {type(df)}")
            return
            
        required_cols = ['binding_type', 'expression_status']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            logger.error(f"Missing required columns: {missing_cols}")
            return
            
        # Rest of the function...
    
    except Exception as e:
        logger.error(f"Error in binding pattern analysis: {str(e)}")

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
                               cell_type_prefix: str,
                               genome_fasta: str,
                               use_cache: bool = True) -> pd.DataFrame:
    """Calculate methylation levels for genomic regions"""
    # Ensure we have the required columns
    required_cols = ['chr', 'start', 'end']
    if not all(col in region_df.columns for col in required_cols):
        raise ValueError(f"Missing required columns. Need {required_cols}, got {region_df.columns.tolist()}")
    
    # Calculate methylation using replicates
    methylation_data = calculate_methylation_levels_with_replicates(
        region_df,
        medip_dir,
        cell_type_prefix,
        genome_fasta
    )
    
    # Add region information to methylation data
    methylation_data = methylation_data.copy()
    methylation_data['start'] = region_df['start'].values
    methylation_data['end'] = region_df['end'].values
    
    # Add quality metrics
    methylation_data = add_methylation_qc_metrics(methylation_data)
    
    # Validate results
    methylation_data = validate_methylation_levels(methylation_data)
    
    # Remove temporary columns if needed
    if 'start' in methylation_data.columns:
        methylation_data = methylation_data.drop(['start', 'end'], axis=1)
    
    return methylation_data

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

def validate_and_merge_data(genes_df, expr_df, cpg_associations):
    """Validate and merge data with better error handling"""
    try:
        # Validate input data
        for df, name in [(genes_df, 'genes_df'), (expr_df, 'expr_df'), 
                        (cpg_associations, 'cpg_associations')]:
            if df is None or len(df) == 0:
                raise ValueError(f"Empty or None {name}")
        
        # First merge: gene annotations with expression data
        merged_df = genes_df.merge(expr_df, on='gene_id', how='inner')
        logger.info(f"Initial merge shape: {merged_df.shape}")
        
        if len(merged_df) == 0:
            raise ValueError("No matches between genes and expression data")
            
        # Add binding information
        merged_df['binding_type'] = 'none'  # Default value
        
        # Add CpG associations
        final_df = merged_df.merge(cpg_associations, on='gene_id', how='left')
        logger.info(f"Final merge shape: {final_df.shape}")
        
        # Validate final data
        if len(final_df) == 0:
            raise ValueError("No data after all merges")
            
        return final_df
        
    except Exception as e:
        logger.error(f"Error in validate_and_merge_data: {str(e)}")
        logger.error("Available columns:")
        for df, name in [(genes_df, 'genes_df'), (expr_df, 'expr_df'),
                        (cpg_associations, 'cpg_associations')]:
            if df is not None:
                logger.error(f"{name}: {df.columns.tolist()}")
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

def process_task_wrapper(args):
    """Wrapper function to unpack arguments for process_analysis_task"""
    return process_analysis_task(*args)

def analyze_methylation_patterns_parallel(results, output_dir, n_processes=None):
    """Parallel version of methylation pattern analysis"""
    if n_processes is None:
        n_processes = max(1, multiprocessing.cpu_count() - 1)

    logger.info(f"Running analysis with {n_processes} processes")

    # Create tasks for parallel processing
    tasks = []
    for cell_type, cell_results in results.items():
        methylation_df = cell_results.get('methylation')
        if isinstance(methylation_df, pd.DataFrame) and len(methylation_df) > 0:
            tasks.extend([
                (methylation_df, cell_type, output_dir, 'distribution'),
                (methylation_df, cell_type, output_dir, 'expression'),
                (methylation_df, cell_type, output_dir, 'binding')
            ])
        else:
            logger.warning(f"No valid methylation data for {cell_type}")
    
    if not tasks:
        logger.warning("No valid tasks for parallel analysis")
        return {}
    
    # Process tasks in parallel
    analysis_results = {}
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = {executor.submit(run_analysis_task, task): task for task in tasks}
        
        for future in tqdm(
            as_completed(futures), 
            total=len(futures),
            desc="Analyzing methylation patterns"
        ):
            task = futures[future]
            try:
                result = future.result()
                if result is not None:
                    cell_type = task[1]
                    analysis_type = task[3]
                    if cell_type not in analysis_results:
                        analysis_results[cell_type] = {}
                    analysis_results[cell_type][analysis_type] = result
            except Exception as e:
                logger.error(f"Error processing task {task}: {str(e)}")
    
    return analysis_results

def determine_binding_type(row):
    """Determine binding type based on exo and endo signals"""
    try:
        exo = float(row['exo_signal'])
        endo = float(row['endo_signal'])
        
        # Define thresholds (you may need to adjust these)
        signal_threshold = 1.0
        
        if exo > signal_threshold and endo > signal_threshold:
            return 'common'
        elif exo > signal_threshold:
            return 'exo_enriched'
        elif endo > signal_threshold:
            return 'endo_only'
        else:
            return 'none'
    except (ValueError, TypeError):
        logger.warning(f"Invalid signal values: exo={row['exo_signal']}, endo={row['endo_signal']}")
        return 'none'

def analyze_binding_enrichment_groups(methylation_df, output_dir):
    """Analyze binding enrichment with improved validation and signal-based categorization"""
    try:
        if not isinstance(methylation_df, pd.DataFrame):
            logger.error(f"Expected DataFrame, got {type(methylation_df)}")
            return None
            
        # Debug information about the DataFrame
        logger.info(f"DataFrame shape: {methylation_df.shape}")
        logger.info(f"Columns available: {methylation_df.columns.tolist()}")
        
        # Determine binding types based on signal values
        if 'exo_signal' in methylation_df.columns and 'endo_signal' in methylation_df.columns:
            logger.info("Determining binding types from signal values...")
            
            # Log signal value distributions
            logger.info("\nExo signal distribution:")
            logger.info(methylation_df['exo_signal'].describe())
            logger.info("\nEndo signal distribution:")
            logger.info(methylation_df['endo_signal'].describe())
            
            # Assign binding types based on signals
            methylation_df['binding_type'] = methylation_df.apply(determine_binding_type, axis=1)
            
        binding_counts = methylation_df['binding_type'].value_counts()
        logger.info(f"\nBinding type counts:\n{binding_counts}")
        
        # Sample of data for each binding type
        for binding_type in methylation_df['binding_type'].unique():
            subset = methylation_df[methylation_df['binding_type'] == binding_type].head()
            logger.info(f"\nSample of {binding_type} genes:")
            sample_cols = ['gene_id', 'gene_name', 'binding_type', 'exo_signal', 'endo_signal']
            logger.info(subset[sample_cols])
        
        # Rest of the function remains the same...
        binding_categories = {
            'exo_enriched': methylation_df[
                methylation_df['binding_type'] == 'exo_enriched'
            ],
            'endo_only': methylation_df[
                methylation_df['binding_type'] == 'endo_only'
            ],
            'common': methylation_df[
                methylation_df['binding_type'] == 'common'
            ]
        }
        
        results = {}
        for category, subset in binding_categories.items():
            if len(subset) == 0:
                logger.warning(f"No genes in {category} (from {len(methylation_df)} total genes)")
                continue
                
            logger.info(f"Found {len(subset)} genes in {category}")
            results[category] = {
                'gene_list': subset,
                'regulation_distribution': subset['expression_status'].value_counts(),
                'methylation_stats': {
                    'promoter': subset['promoter_methylation'].describe(),
                    'gene_body': subset['gene_body_methylation'].describe()
                }
            }
            
        if not results:
            logger.warning("No binding categories had data")
            return None
            
        return results
        
    except Exception as e:
        logger.error(f"Error in binding analysis: {str(e)}")
        logger.error(f"Full error: {str(e.__class__.__name__)}: {str(e)}")
        import traceback
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        return None

def save_category_summary(category, regulation_dist, methylation_stats, subset, output_dir):
    """Save summary for a binding category"""
    summary_file = os.path.join(output_dir, 'binding_enrichment', 
                               f'{category}_analysis.txt')
    os.makedirs(os.path.dirname(summary_file), exist_ok=True)
    
    with open(summary_file, 'w') as f:
        f.write(f"{category.upper()} ANALYSIS\n")
        f.write("="*50 + "\n\n")
        
        f.write("Regulation Distribution:\n")
        f.write("-"*30 + "\n")
        f.write(regulation_dist.to_string())
        f.write("\n\n")
        
        f.write("Methylation Statistics:\n")
        f.write("-"*30 + "\n")
        f.write(methylation_stats.to_string())
        f.write("\n\n")
        
        # Save gene lists by regulation status
        for status in ['not_regulated', 'upregulated', 'downregulated']:
            status_genes = subset[subset['expression_status'] == status]
            if len(status_genes) > 0:
                f.write(f"\n{status.upper()} GENES:\n")
                f.write("-"*30 + "\n")
                for _, gene in status_genes.iterrows():
                    f.write(f"Gene: {gene['gene_name']}\n")
                    f.write(f"Promoter methylation: {gene['promoter_methylation']:.2f}\n")
                    f.write(f"Gene body methylation: {gene['gene_body_methylation']:.2f}\n")
                    f.write("-"*20 + "\n")

def create_group_visualizations(df, category, cell_type, output_dir):
    """Create visualizations for a binding category"""
    try:
        # Set default style
        plt.style.use('default')
        
        # Create output directory
        plot_dir = os.path.join(output_dir, 'visualizations', cell_type, category)
        os.makedirs(plot_dir, exist_ok=True)
        
        # Methylation patterns
        plt.figure(figsize=(12, 6))
        sns.set_style("whitegrid")
        
        # Create violin plots
        plot_data = pd.melt(
            df,
            value_vars=['promoter_methylation', 'gene_body_methylation'],
            id_vars=['expression_status'],
            var_name='Region',
            value_name='Methylation'
        )
        
        sns.violinplot(
            data=plot_data,
            x='expression_status',
            y='Methylation',
            hue='Region',
            split=True
        )
        
        plt.title(f'{cell_type} - {category}\nMethylation Patterns')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'methylation_patterns.pdf'))
        plt.close()
        
    except Exception as e:
        logger.error(f"Error creating visualizations for {category} in {cell_type}: {str(e)}")

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

def create_analysis_visualizations(df, cell_type, output_dir):
    """Create detailed analysis visualizations"""
    try:
        # Set default style
        plt.style.use('default')
        
        # Create output directory
        plot_dir = os.path.join(output_dir, 'visualizations', cell_type, 'analysis')
        os.makedirs(plot_dir, exist_ok=True)
        
        # 1. Methylation distribution
        plt.figure(figsize=(10, 6))
        sns.set_style("whitegrid")
        
        sns.boxplot(data=df, x='binding_type', y='promoter_methylation')
        plt.title(f'{cell_type} - Promoter Methylation by Binding Type')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'promoter_methylation.pdf'))
        plt.close()
        
        # 2. Expression relationship
        plt.figure(figsize=(10, 6))
        sns.scatterplot(
            data=df,
            x='promoter_methylation',
            y='log2FoldChange',
            hue='binding_type',
            alpha=0.6
        )
        plt.title(f'{cell_type} - Expression vs Methylation')
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, 'expression_methylation.pdf'))
        plt.close()
        
    except Exception as e:
        logger.error(f"Error creating analysis visualizations for {cell_type}: {str(e)}")

def save_analysis_results(results, cell_type, output_dir):
    """Save analysis results with better error handling"""
    try:
        # Create output directory
        analysis_dir = os.path.join(output_dir, 'analysis_results', cell_type)
        os.makedirs(analysis_dir, exist_ok=True)
        
        # Save methylation results
        if results.get('methylation') is not None:
            pd.DataFrame(results['methylation']).to_csv(
                os.path.join(analysis_dir, 'methylation_results.csv'),
                index=False
            )
        
        # Save binding results
        if results.get('binding') is not None:
            with open(os.path.join(analysis_dir, 'binding_results.txt'), 'w') as f:
                f.write(f"Binding Analysis Results - {cell_type}\n")
                f.write("="*50 + "\n\n")
                for category, data in results['binding'].items():
                    f.write(f"\n{category.upper()}:\n")
                    f.write("-"*30 + "\n")
                    if isinstance(data, dict):
                        for key, value in data.items():
                            f.write(f"\n{key}:\n{value}\n")
        
        # Save CpG analysis results
        if results.get('cpg') is not None:
            with open(os.path.join(analysis_dir, 'cpg_analysis.txt'), 'w') as f:
                f.write(f"CpG Analysis Results - {cell_type}\n")
                f.write("="*50 + "\n\n")
                for key, value in results['cpg'].items():
                    f.write(f"\n{key}:\n{value}\n")
        
    except Exception as e:
        logger.error(f"Error saving results for {cell_type}: {str(e)}")

def plot_methylation_comparison(results, output_dir):
    """Create comparative plots of methylation patterns"""
    if not results:  # Add check for None
        logger.warning("No results available for methylation comparison")
        return
    
    try:
        plt.figure(figsize=(12, 6))
        for category, data in results.items():
            if 'gene_list' not in data:
                continue
                
            df = data['gene_list']
            sns.boxplot(
                data=df,
                x='expression_status',
                y='promoter_methylation',
                label=category
            )
        
        plt.title('Methylation Patterns by Expression Status')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'methylation_comparison.pdf'))
        plt.close()
        
    except Exception as e:
        logger.error(f"Error in plot_methylation_comparison: {str(e)}")

def identify_promoter_cpg_islands(genes_df, cpg_islands_bed, upstream_distance=2000, downstream_distance=500):
    """
    Identify CpG islands in promoter regions of genes
    """
    try:
        # Validate input data
        required_columns = ['chr', 'start', 'end', 'strand', 'gene_id']
        if not all(col in genes_df.columns for col in required_columns):
            missing = [col for col in required_columns if col not in genes_df.columns]
            raise ValueError(f"Missing required columns in genes_df: {missing}")
        
        if not os.path.exists(cpg_islands_bed):
            raise FileNotFoundError(f"CpG islands BED file not found: {cpg_islands_bed}")
        
        # Create DataFrame in correct format for PyRanges
        pr_df = pd.DataFrame({
            'Chromosome': genes_df['chr'],
            'Start': genes_df['start'],
            'End': genes_df['end'],
            'Strand': genes_df['strand'],
            'gene_id': genes_df['gene_id']
        })
        
        # Convert to PyRanges using from_dict
        genes_pr = pr.from_dict(pr_df.to_dict('list'))
        
        logger.info(f"Created PyRanges object with {len(genes_df)} genes")
        
        # Load CpG islands
        cpg_islands = pr.read_bed(cpg_islands_bed)
        logger.info(f"Loaded {len(cpg_islands.df)} CpG islands from BED file")
        
        # Define promoter regions based on strand
        promoter_regions = []
        for _, gene in genes_df.iterrows():
            if gene['strand'] == '+':
                start = max(0, gene['start'] - upstream_distance)
                end = gene['start'] + downstream_distance
            else:
                start = gene['end'] - downstream_distance
                end = gene['end'] + upstream_distance
                
            promoter_regions.append({
                'Chromosome': gene['chr'],
                'Start': start,
                'End': end,
                'gene_id': gene['gene_id'],
                'Strand': gene['strand']
            })
        
        # Convert promoter regions to PyRanges
        promoters_pr = pr.from_dict(pd.DataFrame(promoter_regions).to_dict('list'))
        
        # Find overlaps between promoters and CpG islands
        overlaps = promoters_pr.join(cpg_islands)
        
        # Calculate overlap characteristics
        cpg_gene_associations = []
        for _, overlap in overlaps.df.iterrows():
            overlap_size = overlap.End - overlap.Start
            distance_to_tss = min(
                abs(overlap.Start - overlap.Start_b),
                abs(overlap.End - overlap.Start_b)
            )
            
            cpg_gene_associations.append({
                'gene_id': overlap.gene_id,  # Using gene_id instead of Name
                'cpg_chr': overlap.Chromosome,
                'cpg_start': overlap.Start,
                'cpg_end': overlap.End,
                'overlap_size': overlap_size,
                'distance_to_tss': distance_to_tss,
                'gc_content': calculate_gc_content(overlap.Chromosome, overlap.Start, overlap.End)
            })
        
        result_df = pd.DataFrame(cpg_gene_associations)
        logger.info(f"Found {len(result_df)} CpG island-gene associations")
        
        return result_df
        
    except Exception as e:
        logger.error(f"Error in identify_promoter_cpg_islands: {str(e)}")
        raise

def calculate_gc_content(chromosome, start, end, genome_fasta=PATHS['genome_fasta']):
    """Calculate GC content for a given genomic region"""
    with pysam.FastaFile(genome_fasta) as fasta:
        sequence = fasta.fetch(chromosome, start, end)
        gc_count = sequence.count('G') + sequence.count('C')
        total = len(sequence)
        return (gc_count / total) * 100 if total > 0 else 0

def analyze_cpg_methylation_patterns(methylation_df, cpg_associations):
    """Analyze CpG methylation patterns"""
    try:
        if not isinstance(methylation_df, pd.DataFrame):
            logger.error(f"Expected DataFrame, got {type(methylation_df)}")
            return None
            
        if len(methylation_df) == 0:
            logger.warning("Empty methylation data")
            return None
            
        # Analyze CpG patterns
        cpg_results = {
            'cpg_targeting': analyze_cpg_targeting(methylation_df),
            'methylation_levels': analyze_methylation_levels(methylation_df),
            'regulation_patterns': analyze_regulation_patterns(methylation_df)
        }
        
        return cpg_results
        
    except Exception as e:
        logger.error(f"Error in CpG methylation analysis: {str(e)}")
        return None

def analyze_cpg_targeting(methylation_df):
    """Analyze CpG island targeting patterns"""
    # Implement CpG island targeting analysis logic
    # This is a placeholder and should be replaced with actual implementation
    return {
        'exo_only': methylation_df[methylation_df['binding_type'] == 'exo_only'],
        'common_targets': methylation_df[methylation_df['binding_type'] == 'both'],
        'exo_enriched': methylation_df[methylation_df['binding_type'] == 'exo_enriched']
    }

def analyze_methylation_levels(methylation_df):
    """Analyze methylation levels"""
    # Implement methylation levels analysis logic
    # This is a placeholder and should be replaced with actual implementation
    return {
        'mean': methylation_df['methylation'].mean(),
        'std': methylation_df['methylation'].std(),
        'median': methylation_df['methylation'].median(),
        'q1': methylation_df['methylation'].quantile(0.25),
        'q3': methylation_df['methylation'].quantile(0.75)
    }

def analyze_regulation_patterns(methylation_df):
    """Analyze regulation patterns"""
    # Implement regulation patterns analysis logic
    # This is a placeholder and should be replaced with actual implementation
    return {
        'upregulated': methylation_df[methylation_df['expression_status'] == 'upregulated'],
        'downregulated': methylation_df[methylation_df['expression_status'] == 'downregulated'],
        'unchanged': methylation_df[methylation_df['expression_status'] == 'unchanged']
    }

def analyze_cpg_methylation_patterns(df, cpg_associations):
    """
    Analyze methylation patterns in relation to CpG islands and gene regulation
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with gene expression and methylation data
    cpg_associations : pd.DataFrame
        DataFrame with CpG island associations
    """
    # Merge gene data with CpG associations
    merged_data = df.merge(cpg_associations, on='gene_id', how='inner')
    
    # Categorize genes based on CpG island characteristics
    merged_data['cpg_category'] = merged_data.apply(
        lambda x: categorize_cpg_region(
            x['gc_content'], 
            x['overlap_size'], 
            x['distance_to_tss']
        ), 
        axis=1
    )
    
    # Analyze regulation patterns for each category
    categories = ['high_cpg', 'intermediate_cpg', 'low_cpg']
    results = {}
    
    for category in categories:
        subset = merged_data[merged_data['cpg_category'] == category]
        results[category] = {
            'regulation_counts': subset['expression_status'].value_counts(),
            'methylation_stats': subset.groupby('expression_status').agg({
                'promoter_methylation': ['mean', 'std', 'count'],
                'gene_body_methylation': ['mean', 'std', 'count']
            })
        }
    
    return results

def categorize_cpg_region(gc_content, size, distance_to_tss):
    """Categorize CpG regions based on their characteristics"""
    if (gc_content >= 50 and 
        size >= 200 and 
        distance_to_tss <= 1000):
        return 'high_cpg'
    elif (gc_content >= 40 and 
          size >= 100 and 
          distance_to_tss <= 1500):
        return 'intermediate_cpg'
    else:
        return 'low_cpg'

# Core analysis functions
def analyze_methylation_patterns_parallel(results, output_dir, n_processes=None):
    """Parallel version of methylation pattern analysis"""
    if n_processes is None:
        n_processes = max(1, multiprocessing.cpu_count() - 1)
    
    logger.info(f"Running analysis with {n_processes} processes")
    
    # Create tasks for parallel processing
    tasks = []
    for cell_type, cell_results in results.items():
        methylation_df = cell_results.get('methylation')
        if isinstance(methylation_df, pd.DataFrame) and len(methylation_df) > 0:
        tasks.extend([
                (methylation_df, cell_type, output_dir, 'distribution'),
                (methylation_df, cell_type, output_dir, 'expression'),
                (methylation_df, cell_type, output_dir, 'binding')
            ])
        else:
            logger.warning(f"No valid methylation data for {cell_type}")
    
    if not tasks:
        logger.warning("No valid tasks for parallel analysis")
        return {}

    # Process tasks in parallel
    analysis_results = {}
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = {executor.submit(run_analysis_task, task): task for task in tasks}
        
        for future in tqdm(
            as_completed(futures), 
            total=len(futures),
            desc="Analyzing methylation patterns"
        ):
            task = futures[future]
            try:
                result = future.result()
                if result is not None:
                    cell_type = task[1]
                    analysis_type = task[3]
                    if cell_type not in analysis_results:
                        analysis_results[cell_type] = {}
                    analysis_results[cell_type][analysis_type] = result
            except Exception as e:
                logger.error(f"Error processing task {task}: {str(e)}")
    
    return analysis_results

def process_analysis_task(df, cell_type, output_dir, analysis_type):
    """Process individual analysis tasks"""
    try:
        if not isinstance(df, pd.DataFrame):
            logger.error(f"Expected DataFrame for {analysis_type}, got {type(df)}")
            return
            
        required_cols = {
            'distribution': ['promoter_methylation', 'gene_body_methylation'],
            'expression': ['promoter_methylation', 'gene_body_methylation', 'log2FoldChange'],
            'binding': ['binding_type', 'expression_status']
        }
        
        missing = [col for col in required_cols.get(analysis_type, []) 
                  if col not in df.columns]
        if missing:
            logger.error(f"Missing required columns for {analysis_type}: {missing}")
            return
            
        if analysis_type == 'distribution':
            analyze_methylation_distribution(df, cell_type, output_dir)
        elif analysis_type == 'expression':
            analyze_methylation_expression_relationship(df, cell_type, output_dir)
        elif analysis_type == 'binding':
            analyze_binding_patterns(df, cell_type, output_dir)
            
    except Exception as e:
        logger.error(f"Error in {analysis_type} analysis for {cell_type}: {str(e)}")

def analyze_cell_types(genes_df, expression_data, cpg_associations, args):
    """Analyze methylation patterns for each cell type"""
    results = {}
    
    if not expression_data:
        logger.error("No expression data provided")
        return results
        
    for cell_type, expr_df in expression_data.items():
        logger.info(f"\nAnalyzing {cell_type}...")
        
        try:
            # Merge data
            merged_df = merge_analysis_data(genes_df, expr_df, cpg_associations)
            
            if merged_df is None or len(merged_df) == 0:
                logger.warning(f"No valid data for {cell_type}")
                continue
            
            # Calculate signal values
            bigwig_files = {
                'exo': [f for f in os.listdir(PATHS['bigwig_dir']) if 'exo' in f.lower()],
                'endo': [f for f in os.listdir(PATHS['bigwig_dir']) if 'endo' in f.lower()]
            }
            
            merged_df = calculate_signal_values(merged_df, bigwig_files)
            
            # Calculate methylation levels
            methylation_df = calculate_methylation_levels(merged_df, args)
            
            results[cell_type] = {
                'methylation': methylation_df,
                'binding': analyze_binding_enrichment_groups(methylation_df, PATHS['output_dir']),
                'cpg': analyze_cpg_methylation_patterns(methylation_df, cpg_associations)
            }
            
        except Exception as e:
            logger.error(f"Error analyzing {cell_type}: {str(e)}")
            continue
    
    return results
def print_summary_statistics(cell_type, binding_analysis, cpg_analysis):
    """Print summary statistics for a cell type analysis"""
    print(f"\n{cell_type.upper()} ANALYSIS SUMMARY")
    print("="*50)
    
    # Binding analysis summary
    print("\nBinding Analysis:")
    print("-"*30)
    if binding_analysis:
        for category, data in binding_analysis.items():
            print(f"\n{category.upper()}:")
            if 'regulation_distribution' in data:
                print("\nRegulation Distribution:")
                print(data['regulation_distribution'])
            if 'methylation_stats' in data:
                print("\nMethylation Statistics:")
                print(data['methylation_stats'])
            else:
        print("No binding analysis results available")
    
    # CpG analysis summary
    print("\nCpG Analysis:")
    print("-"*30)
    if cpg_analysis:
        for category, data in cpg_analysis.items():
            print(f"\n{category.upper()}:")
            if 'regulation_counts' in data:
                print("\nRegulation Counts:")
                print(data['regulation_counts'])
            if 'methylation_stats' in data:
                print("\nMethylation Statistics:")
                print(data['methylation_stats'])
    else:
        print("No CpG analysis results available")

def store_and_visualize_results(methylation_results, binding_analysis, 
                              cpg_analysis, cell_type):
    """Store results and create initial visualizations"""
    results = {
        'methylation': methylation_results,
        'binding': binding_analysis,
        'cpg': cpg_analysis
    }
    
    # Create initial visualizations only if we have binding analysis results
    if binding_analysis:
        plot_methylation_comparison(binding_analysis, PATHS['output_dir'])
    
    # Print summary statistics
    print_summary_statistics(cell_type, binding_analysis, cpg_analysis)
    
    return results

def analyze_category(subset):
    """Analyze a single binding category"""
    return {
        'regulation_distribution': subset['expression_status'].value_counts(),
        'methylation_stats': subset.groupby('expression_status').agg({
            'promoter_methylation': ['mean', 'std', 'count'],
            'gene_body_methylation': ['mean', 'std', 'count']
        }),
        'gene_list': subset[['gene_name', 'expression_status', 
                           'promoter_methylation', 'gene_body_methylation']]
    }

def analyze_exo_enriched_genes(results):
    """
    Analyze genes enriched in Exo binding, focusing on their regulation and methylation patterns
    
    Parameters:
    -----------
    results : dict
        Dictionary containing analysis results for each cell type
    """
    logger.info("\nAnalyzing exo-enriched genes...")
    
    if not results:
        logger.warning("No results available for exo-enriched analysis")
        return
    
    for cell_type, cell_results in results.items():
        methylation_data = cell_results.get('methylation')
        binding_data = cell_results.get('binding')
        
        try:
            if not methylation_data is not None or binding_data is None:
                logger.warning(f"Missing required data for {cell_type}")
                continue
            
            # Filter for exo-enriched genes
            exo_genes = methylation_data[
                methylation_data['binding_type'] == 'exo_only'
            ].copy()
            
            if len(exo_genes) == 0:
                logger.warning(f"No exo-enriched genes found for {cell_type}")
                continue
            
            # Analyze regulation distribution
            regulation_dist = exo_genes['expression_status'].value_counts()
            logger.info(f"\n{cell_type} - Regulation distribution of exo-enriched genes:")
            logger.info(regulation_dist)
            
            # Analyze methylation patterns by regulation status
            methylation_stats = exo_genes.groupby('expression_status').agg({
                'promoter_methylation': ['mean', 'std', 'count'],
                'gene_body_methylation': ['mean', 'std', 'count']
            })
            logger.info(f"\nMethylation statistics:")
            logger.info(methylation_stats)
            
            # Save detailed results
            save_exo_enriched_analysis(
                cell_type,
                exo_genes,
                regulation_dist,
                methylation_stats,
                PATHS['output_dir']
            )
            
            # Create visualizations
            create_exo_enriched_visualizations(
                exo_genes,
                cell_type,
                PATHS['output_dir']
            )
            
            except Exception as e:
            logger.error(f"Error analyzing exo-enriched genes for {cell_type}: {str(e)}")
            continue

def save_exo_enriched_analysis(cell_type, exo_genes, regulation_dist, 
                             methylation_stats, output_dir):
    """Save detailed analysis of exo-enriched genes"""
    output_dir = os.path.join(output_dir, 'exo_enriched_analysis')
    os.makedirs(output_dir, exist_ok=True)
    
    # Save summary file
    summary_file = os.path.join(output_dir, f'{cell_type}_exo_enriched_summary.txt')
    with open(summary_file, 'w') as f:
        f.write(f"EXO-ENRICHED GENES ANALYSIS - {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        f.write("Regulation Distribution:\n")
        f.write("-"*30 + "\n")
        f.write(regulation_dist.to_string())
        f.write("\n\n")
        
        f.write("Methylation Statistics:\n")
        f.write("-"*30 + "\n")
        f.write(methylation_stats.to_string())
        f.write("\n\n")
        
        # Write gene lists by regulation status
        for status in ['not_regulated', 'upregulated', 'downregulated']:
            status_genes = exo_genes[exo_genes['expression_status'] == status]
            if len(status_genes) > 0:
                f.write(f"\n{status.upper()} GENES:\n")
                f.write("-"*30 + "\n")
                for _, gene in status_genes.iterrows():
                    f.write(f"Gene: {gene['gene_name']}\n")
                    f.write(f"Promoter methylation: {gene['promoter_methylation']:.2f}\n")
                    f.write(f"Gene body methylation: {gene['gene_body_methylation']:.2f}\n")
                    f.write("-"*20 + "\n")
    
    # Save detailed gene list
    exo_genes.to_csv(
        os.path.join(output_dir, f'{cell_type}_exo_enriched_genes.csv'),
        index=False
    )

def create_exo_enriched_visualizations(exo_genes, cell_type, output_dir):
    """Create visualizations for exo-enriched genes analysis"""
    plot_dir = os.path.join(output_dir, 'exo_enriched_analysis', 'plots')
    os.makedirs(plot_dir, exist_ok=True)
    
    try:
        # 1. Methylation patterns by regulation status
        plt.figure(figsize=(12, 6))
        
        # Create violin plots for promoter and gene body methylation
        plot_data = pd.melt(
            exo_genes,
            value_vars=['promoter_methylation', 'gene_body_methylation'],
            id_vars=['expression_status'],
            var_name='Region',
            value_name='Methylation'
        )
        
        sns.violinplot(
            data=plot_data,
            x='expression_status',
            y='Methylation',
            hue='Region',
            split=True
        )
        
        plt.title(f'{cell_type} - Exo-enriched Genes\nMethylation Patterns')
        plt.xlabel('Expression Status')
        plt.ylabel('Methylation Level (%)')
        
        plt.tight_layout()
        plt.savefig(
            os.path.join(plot_dir, f'{cell_type}_exo_enriched_methylation.pdf'),
            bbox_inches='tight'
        )
        plt.close()
        
        # 2. Regulation distribution
        plt.figure(figsize=(8, 6))
        regulation_counts = exo_genes['expression_status'].value_counts()
        
        sns.barplot(x=regulation_counts.index, y=regulation_counts.values)
        plt.title(f'{cell_type} - Exo-enriched Genes\nRegulation Distribution')
        plt.xlabel('Expression Status')
        plt.ylabel('Number of Genes')
        
        plt.tight_layout()
        plt.savefig(
            os.path.join(plot_dir, f'{cell_type}_exo_enriched_regulation.pdf'),
            bbox_inches='tight'
        )
        plt.close()
        
    except Exception as e:
        logger.error(f"Error creating exo-enriched visualizations: {str(e)}")

def run_analysis_task(task):
    """Run a single analysis task"""
    df, cell_type, output_dir, analysis_type = task
    try:
        if not isinstance(df, pd.DataFrame):
            logger.error(f"Expected DataFrame for {analysis_type}, got {type(df)}")
            return None
            
        if analysis_type == 'distribution':
            return analyze_methylation_distribution(df, cell_type, output_dir)
        elif analysis_type == 'expression':
            return analyze_methylation_expression_relationship(df, cell_type, output_dir)
        elif analysis_type == 'binding':
            return analyze_binding_patterns(df, cell_type, output_dir)
            
    except Exception as e:
        logger.error(f"Error in {analysis_type} analysis for {cell_type}: {str(e)}")
        return None

def analyze_methylation_patterns_parallel(results, output_dir, n_processes=None):
    """Parallel version of methylation pattern analysis"""
    if n_processes is None:
        n_processes = max(1, multiprocessing.cpu_count() - 1)
    
    logger.info(f"Running analysis with {n_processes} processes")
    
    # Create tasks for parallel processing
    tasks = []
    for cell_type, cell_results in results.items():
        methylation_df = cell_results.get('methylation')
        if isinstance(methylation_df, pd.DataFrame) and len(methylation_df) > 0:
            tasks.extend([
                (methylation_df, cell_type, output_dir, 'distribution'),
                (methylation_df, cell_type, output_dir, 'expression'),
                (methylation_df, cell_type, output_dir, 'binding')
            ])
        else:
            logger.warning(f"No valid methylation data for {cell_type}")
    
    if not tasks:
        logger.warning("No valid tasks for parallel analysis")
        return {}
    
    # Process tasks in parallel
    analysis_results = {}
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = {executor.submit(run_analysis_task, task): task for task in tasks}
        
        for future in tqdm(
            as_completed(futures), 
            total=len(futures),
            desc="Analyzing methylation patterns"
        ):
            task = futures[future]
            try:
                result = future.result()
                if result is not None:
                    cell_type = task[1]
                    analysis_type = task[3]
                    if cell_type not in analysis_results:
                        analysis_results[cell_type] = {}
                    analysis_results[cell_type][analysis_type] = result
            except Exception as e:
                logger.error(f"Error processing task {task}: {str(e)}")
    
    return analysis_results

def calculate_signal_values(gene_df: pd.DataFrame, bigwig_files: Dict[str, str]) -> pd.DataFrame:
    """Calculate exo and endo signal values for each gene"""
    logger.info("Calculating signal values from bigwig files...")
    
    try:
        # Initialize signal columns
        gene_df['exo_signal'] = 0.0
        gene_df['endo_signal'] = 0.0
        
        # Get file paths
        exo_files = [f for f in bigwig_files.values() if 'exo' in f.lower()]
        endo_files = [f for f in bigwig_files.values() if 'endo' in f.lower()]
        
        if not exo_files or not endo_files:
            logger.error("Missing exo or endo bigwig files")
            return gene_df
            
        logger.info(f"Found {len(exo_files)} exo files and {len(endo_files)} endo files")
        
        # Calculate signals for each gene
        for idx, row in gene_df.iterrows():
            try:
                # Get gene coordinates
                chrom = row['chromosome']
                start = int(row['start'])
                end = int(row['end'])
                
                # Calculate exo signal
                exo_values = []
                for bw_file in exo_files:
                    with pyBigWig.open(bw_file) as bw:
                        try:
                            values = bw.values(chrom, start, end)
                            values = [v for v in values if v is not None]
                            if values:
                                exo_values.append(np.mean(values))
                        except RuntimeError:
                            continue
                
                # Calculate endo signal
                endo_values = []
                for bw_file in endo_files:
                    with pyBigWig.open(bw_file) as bw:
                        try:
                            values = bw.values(chrom, start, end)
                            values = [v for v in values if v is not None]
                            if values:
                                endo_values.append(np.mean(values))
                        except RuntimeError:
                            continue
                
                # Update signal values
                gene_df.at[idx, 'exo_signal'] = np.mean(exo_values) if exo_values else 0.0
                gene_df.at[idx, 'endo_signal'] = np.mean(endo_values) if endo_values else 0.0
                
            except Exception as e:
                logger.warning(f"Error calculating signals for gene {row['gene_id']}: {str(e)}")
                continue
        
        # Log signal value distributions
        logger.info("\nExo signal distribution:")
        logger.info(gene_df['exo_signal'].describe())
        logger.info("\nEndo signal distribution:")
        logger.info(gene_df['endo_signal'].describe())
        
        return gene_df
        
    except Exception as e:
        logger.error(f"Error in signal calculation: {str(e)}")
        return gene_df

def analyze_cell_types(genes_df, expression_data, cpg_associations, args):
    """Analyze methylation patterns for each cell type"""
    results = {}
    
    if not expression_data:
        logger.error("No expression data provided")
        return results
        
    for cell_type, expr_df in expression_data.items():
        logger.info(f"\nAnalyzing {cell_type}...")
        
        try:
            # Merge data
            merged_df = merge_analysis_data(genes_df, expr_df, cpg_associations)
            
            if merged_df is None or len(merged_df) == 0:
                logger.warning(f"No valid data for {cell_type}")
                continue
            
            # Calculate signal values
            bigwig_files = {
                'exo': [f for f in os.listdir(PATHS['bigwig_dir']) if 'exo' in f.lower()],
                'endo': [f for f in os.listdir(PATHS['bigwig_dir']) if 'endo' in f.lower()]
            }
            
            merged_df = calculate_signal_values(merged_df, bigwig_files)
            
            # Calculate methylation levels
            methylation_df = calculate_methylation_levels(merged_df, args)
            
            results[cell_type] = {
                'methylation': methylation_df,
                'binding': analyze_binding_enrichment_groups(methylation_df, PATHS['output_dir']),
                'cpg': analyze_cpg_methylation_patterns(methylation_df, cpg_associations)
            }
            
        except Exception as e:
            logger.error(f"Error analyzing {cell_type}: {str(e)}")
            continue
    
    return results

def load_mecp2_binding_data(chromosome: str = None) -> pd.DataFrame:
    """Load MeCP2 binding data from CSV file"""
    try:
        binding_file = os.path.join(PATHS['data_dir'], 'mecp2_cpg_enrichment_parallel.csv')
        binding_df = pd.read_csv(binding_file)
        
        # Filter for chromosome if specified
        if chromosome:
            binding_df = binding_df[binding_df['chr'] == chromosome].copy()
            
        logger.info(f"Loaded {len(binding_df)} MeCP2 binding regions")
        logger.info(f"Binding types distribution:\n{binding_df['binding_type'].value_counts()}")
        
        return binding_df
    except Exception as e:
        logger.error(f"Error loading MeCP2 binding data: {str(e)}")
        return None

def assign_binding_types(df: pd.DataFrame, binding_df: pd.DataFrame, thresholds: Dict[str, float]) -> pd.DataFrame:
    """
    Assign binding types based on MeCP2 signal data using calculated thresholds
    """
    logger.info("Determining binding types from signal values...")
    
    if binding_df is None or binding_df.empty:
        logger.warning("No binding data provided - setting all binding types to 'none'")
        df['binding_type'] = 'none'
        df['exo_signal'] = None 
        df['endo_signal'] = None
        df['mecp2_bound'] = False
        return df

    # Merge binding data with gene data
    merged = pd.merge(df, binding_df, on='gene_id', how='left')
    
    # Use calculated thresholds
    exo_threshold = thresholds['exo']
    endo_threshold = thresholds['endo']
    
    logger.info(f"Using thresholds - Exo: {exo_threshold:.2f}, Endo: {endo_threshold:.2f}")
    
    # Calculate binding types
    conditions = [
        (merged['exo_signal'] >= exo_threshold) & (merged['endo_signal'] >= endo_threshold),
        (merged['exo_signal'] >= exo_threshold) & (merged['endo_signal'] < endo_threshold),
        (merged['exo_signal'] < exo_threshold) & (merged['endo_signal'] >= endo_threshold),
    ]
    choices = ['both', 'exo_only', 'endo_only']
    
    merged['binding_type'] = np.select(conditions, choices, default='none')
    merged['mecp2_bound'] = merged['binding_type'] != 'none'
    
    # Log binding type distributions
    binding_counts = merged['binding_type'].value_counts()
    logger.info("\nBinding type counts:")
    logger.info(binding_counts)
    
    # Calculate percentages
    binding_pcts = (binding_counts / len(merged) * 100).round(2)
    logger.info("\nBinding type percentages:")
    logger.info(binding_pcts)
    
    return merged

def analyze_cell_types(genes_df, expression_data, cpg_associations, args):
    """Analyze methylation patterns for each cell type"""
    results = {}
    
    if not expression_data:
        logger.error("No expression data provided")
        return results
    
    # Load MeCP2 binding data
    binding_df = load_mecp2_binding_data(args.chromosome)
    if binding_df is None:
        logger.error("Failed to load MeCP2 binding data")
        return results
        
    for cell_type, expr_df in expression_data.items():
        logger.info(f"\nAnalyzing {cell_type}...")
        
        try:
            # Merge data
            merged_df = merge_analysis_data(genes_df, expr_df, cpg_associations)
            
            if merged_df is None or len(merged_df) == 0:
                logger.warning(f"No valid data for {cell_type}")
                continue
            
            # Assign binding types from MeCP2 data
            merged_df = assign_binding_types(merged_df, binding_df, CONFIG['binding_thresholds'])
            
            # Calculate methylation levels
            methylation_df = calculate_methylation_levels(merged_df, args)
            
            results[cell_type] = {
                'methylation': methylation_df,
                'binding': analyze_binding_enrichment_groups(methylation_df, PATHS['output_dir']),
                'cpg': analyze_cpg_methylation_patterns(methylation_df, cpg_associations)
            }
            
        except Exception as e:
            logger.error(f"Error analyzing {cell_type}: {str(e)}")
            continue
    
    return results

def load_binding_data(binding_file: str) -> pd.DataFrame:
    """
    Load pre-calculated MeCP2 binding data
    """
    try:
        if not os.path.exists(binding_file):
            logger.error(f"Binding file not found: {binding_file}")
            return None
            
        binding_df = pd.read_csv(binding_file)
        
        # Verify required columns exist
        required_cols = ['chr', 'start', 'end', 'exo_signal', 'endo_signal', 'binding_type']
        if not all(col in binding_df.columns for col in required_cols):
            logger.error(f"Binding file missing required columns: {required_cols}")
            return None
            
        logger.info(f"Loaded binding data with {len(binding_df)} regions")
        logger.info("\nBinding type distribution:")
        logger.info(binding_df['binding_type'].value_counts())
        
        return binding_df
        
    except Exception as e:
        logger.error(f"Error loading binding data: {str(e)}")
        return None

def assign_binding_types(genes_df: pd.DataFrame, binding_df: pd.DataFrame) -> pd.DataFrame:
    """
    Assign binding types to genes based on overlapping binding regions
    """
    if binding_df is None:
        logger.warning("No binding data provided - setting all binding types to 'none'")
        genes_df['binding_type'] = 'none'
        genes_df['exo_signal'] = None
        genes_df['endo_signal'] = None
        genes_df['mecp2_bound'] = False
        return genes_df

    # Convert to PyRanges for efficient overlap detection
    genes_pr = pr.PyRanges(
        chromosomes=genes_df['chromosome'],
        starts=genes_df['start'],
        ends=genes_df['end'],
        strands=genes_df['strand']
    )
    
    binding_pr = pr.PyRanges(
        chromosomes=binding_df['chr'],
        starts=binding_df['start'],
        ends=binding_df['end'],
        strands=None
    )

    # Find overlaps
    overlaps = genes_pr.join(binding_pr).as_df()
    
    # Initialize columns
    genes_df['binding_type'] = 'none'
    genes_df['exo_signal'] = 0.0
    genes_df['endo_signal'] = 0.0
    genes_df['mecp2_bound'] = False

    # Process overlaps
    for _, overlap in overlaps.iterrows():
        gene_idx = genes_df.index[
            (genes_df['chromosome'] == overlap['Chromosome']) & 
            (genes_df['start'] == overlap['Start']) & 
            (genes_df['end'] == overlap['End'])
        ].values[0]
        
        current_type = genes_df.at[gene_idx, 'binding_type']
        new_type = overlap['binding_type']
        
        # Update binding type
        if current_type == 'none':
            genes_df.at[gene_idx, 'binding_type'] = new_type
        elif current_type != new_type and new_type != 'none':
            genes_df.at[gene_idx, 'binding_type'] = 'both'
        
        # Update signals (take maximum value)
        genes_df.at[gene_idx, 'exo_signal'] = max(
            genes_df.at[gene_idx, 'exo_signal'],
            overlap['exo_signal']
        )
        genes_df.at[gene_idx, 'endo_signal'] = max(
            genes_df.at[gene_idx, 'endo_signal'],
            overlap['endo_signal']
        )
        
        genes_df.at[gene_idx, 'mecp2_bound'] = True

    # Log results
    logger.info("\nBinding type distribution after assignment:")
    logger.info(genes_df['binding_type'].value_counts())
    
    return genes_df

def classify_expression_changes(expr_df: pd.DataFrame, 
                             padj_threshold: float = 0.05,
                             lfc_threshold: float = 1.0) -> pd.DataFrame:
    """
    Classify genes based on their expression changes
    """
    expr_df['expression_status'] = 'unchanged'
    
    # Significant changes
    significant = expr_df['padj'] < padj_threshold
    
    # Upregulated: log2FC > 1 (2-fold increase)
    expr_df.loc[significant & (expr_df['log2FoldChange'] >= lfc_threshold), 
                'expression_status'] = 'upregulated'
    
    # Downregulated: log2FC < -1 (2-fold decrease)
    expr_df.loc[significant & (expr_df['log2FoldChange'] <= -lfc_threshold), 
                'expression_status'] = 'downregulated'
    
    return expr_df

def identify_promoter_cpgs(genes_df: pd.DataFrame, 
                         cpg_df: pd.DataFrame,
                         upstream_distance: int = 2000,
                         downstream_distance: int = 500) -> pd.DataFrame:
    """
    Identify CpG islands in promoter regions
    """
    # Convert to PyRanges for efficient overlap detection
    genes_pr = pr.PyRanges(
        chromosomes=genes_df['chromosome'],
        starts=genes_df['start'] - upstream_distance,  # Include upstream region
        ends=genes_df['start'] + downstream_distance,  # Include downstream
        strands=genes_df['strand']
    )
    
    cpg_pr = pr.PyRanges(
        chromosomes=cpg_df['chr'],
        starts=cpg_df['start'],
        ends=cpg_df['end']
    )
    
    # Find overlaps
    overlaps = genes_pr.join(cpg_pr).as_df()
    
    return overlaps

def analyze_binding_and_expression(genes_df: pd.DataFrame,
                                 binding_df: pd.DataFrame,
                                 expr_df: pd.DataFrame,
                                 cpg_df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """
    Analyze relationships between MeCP2 binding, CpG islands, and expression
    """
    # 1. Identify promoter CpG islands
    promoter_cpgs = identify_promoter_cpgs(genes_df, cpg_df)
    
    # 2. Assign binding types to genes
    genes_with_binding = assign_binding_types(genes_df, binding_df)
    
    # 3. Classify expression changes
    expr_with_status = classify_expression_changes(expr_df)
    
    # 4. Merge all data
    merged_df = pd.merge(genes_with_binding, expr_with_status, 
                        left_on='gene_id', right_on='gene')
    
    # 5. Analyze each binding category
    results = {}
    for binding_type in ['exo_only', 'endo_only', 'both', 'none']:
        subset = merged_df[merged_df['binding_type'] == binding_type]
        
        # Get expression distribution
        expr_counts = subset['expression_status'].value_counts()
        
        # Calculate methylation statistics
        methylation_stats = {
            'promoter_methylation': subset.groupby('expression_status')['promoter_methylation'].mean(),
            'gene_body_methylation': subset.groupby('expression_status')['gene_body_methylation'].mean()
        }
        
        results[binding_type] = {
            'data': subset,
            'expression_counts': expr_counts,
            'methylation_stats': methylation_stats
        }
        
        # Log results
        logger.info(f"\nResults for {binding_type}:")
        logger.info(f"Expression distribution:\n{expr_counts}")
        logger.info("\nMethylation levels by expression status:")
        logger.info(f"Promoter:\n{methylation_stats['promoter_methylation']}")
        logger.info(f"Gene body:\n{methylation_stats['gene_body_methylation']}")
    
    return results

def visualize_binding_analysis(results: Dict[str, Dict], output_dir: str):
    """
    Create visualizations for binding analysis results
    """
    # 1. Expression distribution by binding type
    plt.figure(figsize=(12, 6))
    expr_data = {bt: res['expression_counts'] for bt, res in results.items()}
    df_expr = pd.DataFrame(expr_data)
    
    sns.barplot(data=df_expr)
    plt.title('Expression Changes by Binding Type')
    plt.ylabel('Number of Genes')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/expression_by_binding.pdf')
    plt.close()
    
    # 2. Methylation levels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Promoter methylation
    methylation_data = []
    for bt, res in results.items():
        for expr_status, meth_level in res['methylation_stats']['promoter_methylation'].items():
            methylation_data.append({
                'binding_type': bt,
                'expression': expr_status,
                'methylation': meth_level,
                'region': 'promoter'
            })
    
    meth_df = pd.DataFrame(methylation_data)
    
    sns.boxplot(data=meth_df, x='binding_type', y='methylation', 
                hue='expression', ax=ax1)
    ax1.set_title('Promoter Methylation')
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45)
    
    # Similar plot for gene body methylation...
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/methylation_patterns.pdf')
    plt.close()

def load_gene_data(cell_type: str) -> pd.DataFrame:
    """
    Load gene data for a specific cell type
    """
    try:
        file_path = os.path.join(PATHS['gene_data'], f"{cell_type}_genes.csv")
        if not os.path.exists(file_path):
            logger.error(f"Gene data file not found: {file_path}")
            return None
            
        df = pd.read_csv(file_path)
        logger.info(f"Loaded gene data for {cell_type}: {len(df)} genes")
        return df
        
    except Exception as e:
        logger.error(f"Error loading gene data for {cell_type}: {str(e)}")
        return None

def load_cpg_associations(cell_type: str) -> pd.DataFrame:
    """
    Load CpG associations for a specific cell type
    """
    try:
        file_path = os.path.join(PATHS['cpg_data'], f"{cell_type}_cpg_associations.csv")
        if not os.path.exists(file_path):
            logger.error(f"CpG associations file not found: {file_path}")
            return None
            
        df = pd.read_csv(file_path)
        logger.info(f"Loaded CpG associations for {cell_type}: {len(df)} associations")
        return df
        
    except Exception as e:
        logger.error(f"Error loading CpG associations for {cell_type}: {str(e)}")
        return None

def merge_analysis_data(genes_df: pd.DataFrame, 
                       expr_df: pd.DataFrame,
                       cpg_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge gene data with expression and CpG information
    """
    try:
        if genes_df is None or expr_df is None or cpg_df is None:
            logger.error("One or more required dataframes is None")
            return None
            
        # Merge gene data with expression data
        merged = pd.merge(genes_df, expr_df, 
                         left_on='gene_id', 
                         right_on='gene',
                         how='inner')
        
        # Merge with CpG associations
        merged = pd.merge(merged, cpg_df,
                         on='gene_id',
                         how='left')
        
        # Fill NaN values for CpG-related columns
        cpg_columns = ['promoter_cpg_count', 'gene_body_cpg_count']
        merged[cpg_columns] = merged[cpg_columns].fillna(0)
        
        logger.info(f"Merged data shape: {merged.shape}")
        return merged
        
    except Exception as e:
        logger.error(f"Error merging analysis data: {str(e)}")
        return None
