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

##########################################################

def validate_and_merge_data(genes_df: pd.DataFrame, expr_df: pd.DataFrame, 
                          mecp2_binding: pd.DataFrame) -> pd.DataFrame:
    """Validate and merge gene annotations with expression data and MeCP2 binding data.
    
    This function performs several key steps:
    1. Merges gene annotations with RNA expression data
    2. Standardizes chromosome formats between datasets (ensuring 'chr' prefix consistency)
    3. Maps MeCP2 binding regions to genes using map_binding_to_genes()
    4. Adds binding information (bound/unbound status, binding type, signal strengths)
    
    Args:
        genes_df: DataFrame containing gene annotations (chr, coordinates, gene IDs)
        expr_df: DataFrame containing RNA expression data (gene IDs, expression values)
        mecp2_binding: DataFrame containing MeCP2 binding regions and signals
        
    Returns:
        DataFrame containing merged gene annotations, expression data and binding information
        with standardized chromosome formats and binding status for each gene
    """
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

##########################################################

def analyze_tss_binding_patterns(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze binding patterns specifically around TSS and their relationship with gene regulation.
    
    This function analyzes MeCP2 binding patterns around transcription start sites (TSS) and how they
    relate to gene regulation. It categorizes genes based on binding type (exo/endo enriched) and 
    regulation status (up/down/unchanged).
    
    For each category, it calculates and saves:
    - Gene lists as CSV files: {output_dir}/tss_analysis/{cell_type}/{binding_category}_{regulation_status}_genes.csv
    - Statistical metrics (counts, methylation means/medians) in: {output_dir}/tss_analysis/{cell_type}/{cell_type}_tss_analysis_stats.txt
    - Visualization plots showing the relationships between binding and regulation
    
    Args:
        df: DataFrame containing gene info, binding data and expression status
        cell_type: Cell type being analyzed (e.g. 'NSC', 'NEU') 
        output_dir: Base output directory for saving results
        
    Returns:
        Dictionary containing statistical results for each binding category and regulation status
    """
    
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
                group_stats[reg_status] = {
                    'count': len(group_df),
                    'promoter_methylation': {
                        'mean': group_df['promoter_methylation'].mean(),
                        'std': group_df['promoter_methylation'].std(),
                        'median': group_df['promoter_methylation'].median()
                    },
                    'gene_body_methylation': {
                        'mean': group_df['gene_body_methylation'].mean(),
                        'std': group_df['gene_body_methylation'].std(),
                        'median': group_df['gene_body_methylation'].median()
                    }
                }
                
                # Save gene lists
                output_file = os.path.join(tss_dir, f'{category}_{reg_status}_genes.csv')
                group_df.to_csv(output_file, index=False)
        
        results[category] = group_stats
    
    # Create visualization
    create_tss_analysis_plots(results, cell_type, tss_dir)
    
    # Save detailed statistics
    save_tss_analysis_results(results, cell_type, tss_dir)
    
    return results

def create_tss_analysis_plots(results: dict, cell_type: str, output_dir: str):
    """Create visualizations for TSS binding analysis
    
    This function generates plots to visualize the relationship between MeCP2 binding 
    and gene regulation at transcription start sites (TSS).
    
    Args:
        results: Dictionary containing analysis results for each binding category
        cell_type: String indicating the cell type being analyzed (e.g. 'NSC', 'NEU')
        output_dir: Directory path where plots should be saved
    
    Creates two main visualizations:
    1. A stacked bar plot showing distribution of regulated genes in each binding category
    2. A methylation comparison plot (created by create_methylation_comparison_plot)
    """
    
    # 1. Regulation distribution plot - Shows counts of regulated genes by binding category
    plt.figure(figsize=(12, 6))
    
    # Initialize lists to store data for each regulation category
    categories = []  # Binding categories (e.g. exo_enriched, endo_only)
    not_dereg = []  # Counts of non-deregulated genes
    up = []         # Counts of upregulated genes  
    down = []       # Counts of downregulated genes
    
    # Extract counts for each category and regulation status
    for category, stats in results.items():
        categories.append(category)
        # Use .get() with default 0 to handle missing categories
        not_dereg.append(stats.get('not_deregulated', {}).get('count', 0))
        up.append(stats.get('upregulated', {}).get('count', 0))
        down.append(stats.get('downregulated', {}).get('count', 0))
    
    # Create stacked bar plot showing all regulation categories
    width = 0.35
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot bars for each regulation category, stacking them
    ax.bar(categories, not_dereg, width, label='Not Deregulated')
    ax.bar(categories, up, width, bottom=not_dereg, label='Upregulated')
    ax.bar(categories, down, width, bottom=[i+j for i,j in zip(not_dereg, up)], 
           label='Downregulated')
    
    # Customize plot appearance
    plt.title(f'{cell_type} - Gene Regulation by Binding Category')
    plt.xlabel('Binding Category')
    plt.ylabel('Number of Genes')
    plt.legend()
    plt.xticks(rotation=45)  # Rotate category labels for better readability
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'regulation_distribution.pdf'))
    plt.close()
    
    # 2. Generate methylation comparison plot
    create_methylation_comparison_plot(results, cell_type, output_dir)

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
                reg_stats['promoter_methylation']['mean']
            )
            plot_data['StdDev'].append(
                reg_stats['promoter_methylation']['std']
            )
            
            # Add gene body methylation data
            plot_data['Category'].append(category)
            plot_data['Regulation'].append(reg_status)
            plot_data['Region'].append('Gene Body')
            plot_data['Methylation'].append(
                reg_stats['gene_body_methylation']['mean']
            )
            plot_data['StdDev'].append(
                reg_stats['gene_body_methylation']['std']
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

##########################################################

def calculate_contextual_methylation(merged_df: pd.DataFrame,
                                   medip_dir: str,
                                   cell_type: str,
                                   genome_fasta: str,
                                   n_processes: int = None) -> pd.DataFrame:
    """
    Calculate methylation levels for promoter and gene body regions in parallel.
    
    This function takes gene location data and calculates methylation levels by:
    1. Extracting promoter and gene body coordinates for each gene
    2. Running parallel methylation calculations for both regions using ThreadPoolExecutor
    3. Combining the methylation results with gene expression and binding data
    
    Args:
        merged_df: DataFrame with gene locations and expression data
        medip_dir: Directory containing MeDIP-seq data
        cell_type: Cell type being analyzed (e.g. 'NSC', 'NEU')
        genome_fasta: Path to genome FASTA file
        n_processes: Number of processes for parallel computation
        
    Returns:
        DataFrame containing:
        - Gene identifiers and names
        - Expression status and fold changes
        - Promoter and gene body methylation levels
        - CpG counts for both regions
        - MeCP2 binding information
    """
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

def calculate_methylation_levels_parallel(region_df: pd.DataFrame,
                                        medip_dir: str,
                                        cell_type_prefix: str,
                                        genome_fasta: str,
                                        n_processes: int = None) -> pd.DataFrame:
    """Parallel version of calculate_methylation_levels that distributes computation across multiple processes
    
    Args:
        region_df: DataFrame containing genomic regions to analyze
        medip_dir: Directory containing MeDIP-seq data files
        cell_type_prefix: Prefix identifying the cell type (e.g. 'NSC', 'NEU') 
        genome_fasta: Path to genome FASTA file
        n_processes: Number of processes to use. Defaults to CPU count - 1
        
    Returns:
        DataFrame with methylation levels calculated for each region
    """
    # Use all but one CPU core if n_processes not specified
    if n_processes is None:
        n_processes = max(1, multiprocessing.cpu_count() - 1)

    # Split data into chunks for parallel processing
    # Create 4x more chunks than processes for better load balancing
    chunk_size = max(1, len(region_df) // (n_processes * 4))  
    chunks = np.array_split(region_df, len(region_df) // chunk_size + 1)

    # Process chunks in parallel using ProcessPoolExecutor
    # Each chunk is processed independently by a separate process
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = [
            executor.submit(process_chunk, chunk, medip_dir, cell_type_prefix, genome_fasta)
            for chunk in chunks
        ]
        
        # Collect results with progress tracking
        # Handle any errors that occur during processing
        results = []
        for future in tqdm(futures, desc="Calculating methylation levels"):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                logger.error(f"Error processing chunk: {str(e)}")
                continue

    # Verify we got valid results before combining
    if not results:
        raise RuntimeError("No valid results obtained from parallel processing")
    
    # Combine all chunk results into single DataFrame
    return pd.concat(results, ignore_index=True)

def process_chunk(chunk: pd.DataFrame, medip_dir: str, cell_type_prefix: str, genome_fasta: str) -> pd.DataFrame:
    """Process a chunk of regions for methylation calculation"""
    return calculate_methylation_levels_with_replicates(
        chunk, medip_dir, cell_type_prefix, genome_fasta
    )

def calculate_methylation_levels_with_replicates(region_df: pd.DataFrame,
                                               medip_dir: str,
                                               cell_type_prefix: str,
                                               genome_fasta: str) -> pd.DataFrame:
    """Calculate methylation levels using improved normalization
    
    This function processes multiple replicate MeDIP-seq experiments to calculate methylation levels:
    1. Loads replicate and input control files for a given cell type
    2. For each genomic region, calculates methylation by comparing IP to input signal
    3. Averages methylation values across replicates
    4. Returns methylation levels and CpG counts for each region
    """
    # Map cell type codes to file prefix conventions
    cell_type_map = {
        'NEU': 'N',  # Neurons
        'NSC': 'PP'  # Neural stem cells
    }
    
    prefix = cell_type_map.get(cell_type_prefix)
    if not prefix:
        raise ValueError(f"Unknown cell type: {cell_type_prefix}")
    
    # Construct file paths for 3 replicates of both IP and input samples
    replicate_files = [
        os.path.join(medip_dir, f"Medip_{prefix}_output_r{i}.bw")
        for i in range(1, 4)
    ]
    input_files = [
        os.path.join(medip_dir, f"Medip_{prefix}_input_r{i}.bw")
        for i in range(1, 4)
    ]
    
    logger.info(f"Processing replicates for {cell_type_prefix} ({prefix})")
    
    # Ensure region coordinates are valid (non-negative starts)
    region_df = region_df.copy()
    region_df['start'] = region_df['start'].clip(lower=0)
    
    # Get chromosome sizes from first bigWig file to validate region coordinates
    first_bw = pyBigWig.open(replicate_files[0])
    chrom_sizes = dict(first_bw.chroms().items())
    first_bw.close()
    
    # Clip region ends to chromosome sizes
    for chrom in region_df['chr'].unique():
        if chrom in chrom_sizes:
            mask = (region_df['chr'] == chrom)
            region_df.loc[mask, 'end'] = region_df.loc[mask, 'end'].clip(upper=chrom_sizes[chrom])
    
    # Process each replicate pair (IP and input)
    replicate_results = []
    for rep_file, input_file in zip(replicate_files, input_files):
        if not (os.path.exists(rep_file) and os.path.exists(input_file)):
            logger.error(f"File not found: {rep_file} or {input_file}")
            continue
            
        logger.info(f"Processing {os.path.basename(rep_file)}")
        try:
            # Open data files
            bw_ip = pyBigWig.open(rep_file)
            bw_input = pyBigWig.open(input_file)
            fasta = pysam.FastaFile(genome_fasta)  # For CpG counting
            
            methylation_values = []
            cpg_counts = []
            
            # Calculate methylation for each region
            for _, row in region_df.iterrows():
                try:
                    # Extract signal values for the region
                    ip_values = bw_ip.values(row['chr'], int(row['start']), int(row['end']))
                    input_values = bw_input.values(row['chr'], int(row['start']), int(row['end']))
                    
                    if ip_values is None or input_values is None:
                        methylation = 0
                        cpg_count = 0
                    else:
                        # Filter out missing values and pair IP with input
                        valid_pairs = [
                            (ip, inp) for ip, inp in zip(ip_values, input_values)
                            if ip is not None and inp is not None
                        ]
                        
                        if valid_pairs:
                            ip_vals, input_vals = zip(*valid_pairs)
                            
                            # Get DNA sequence and calculate normalized methylation
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
            
            # Clean up file handles
            bw_ip.close()
            bw_input.close()
            fasta.close()
            
            replicate_results.append(methylation_values)
            
        except Exception as e:
            logger.error(f"Error processing {rep_file}: {str(e)}")
            continue
    
    if not replicate_results:
        raise RuntimeError(f"No valid data found for {cell_type_prefix}")
    
    # Average methylation across all valid replicates
    avg_methylation = np.mean(replicate_results, axis=0)
    
    return pd.DataFrame({
        'methylation': avg_methylation,
        'cpg_count': cpg_counts
    })

def calculate_normalized_methylation(ip_values, input_values, sequence):
    """Calculate biologically meaningful methylation levels
    
    In mammalian DNA, methylation occurs primarily at CpG sites and is binary 
    at each site (either methylated or unmethylated). The overall methylation
    level represents the proportion of methylated CpGs in the region.
    """
    if not ip_values or not input_values:
        return 0
        
    # Count CpGs in region
    cpg_count = sequence.upper().count('CG')
    if cpg_count == 0:
        return 0
        
    # Calculate average signals
    ip_mean = np.mean([x for x in ip_values if x is not None])
    input_mean = np.mean([x for x in input_values if x is not None])
    
    if input_mean <= 0:
        return 0
    
    # 1. Calculate enrichment relative to input
    # MeDIP specifically pulls down methylated DNA fragments
    enrichment = ip_mean / input_mean
    
    # 2. Normalize by local CpG density
    # Regions with more CpGs will naturally have more MeDIP signal
    region_length = len(sequence)
    cpg_density = cpg_count / region_length
    normalized_enrichment = enrichment / cpg_density
    
    # 3. Convert to methylation percentage
    # Based on typical MeDIP calibration curves, where:
    # - Low enrichment (~0-1x) indicates mostly unmethylated CpGs
    # - Medium enrichment (~2-3x) indicates partially methylated regions
    # - High enrichment (>4x) indicates heavily methylated regions
    if normalized_enrichment <= 1:
        methylation = 25 * normalized_enrichment  # Linear scaling for low values
    else:
        # Asymptotic approach to 100% for higher enrichment
        methylation = 100 * (1 - np.exp(-0.5 * (normalized_enrichment - 1)))
    
    return methylation

##########################################################

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
