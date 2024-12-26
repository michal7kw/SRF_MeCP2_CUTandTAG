import pandas as pd
import pyBigWig
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pysam
import pyranges as pr
import json
from functools import lru_cache
from typing import Dict, List, Tuple, Optional
import logging
from concurrent.futures import ProcessPoolExecutor
from functools import partial

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Configuration
CONFIG = {
    'methylation_thresholds': {
        'hypo': 30,
        'hyper': 70
    },
    'quality_thresholds': {
        'coverage': 0.7,
        'min_cpgs': 5
    },
    'expression_thresholds': {
        'log2fc': 1,
        'padj': 0.05
    },
    'standard_chromosomes': [f'chr{i}' for i in range(1, 20)] + ['chrX', 'chrY']
}

@lru_cache(maxsize=1000)
def count_cpgs_in_sequence(sequence: str) -> int:
    """
    Count CpG dinucleotides in a DNA sequence with caching
    """
    return sequence.upper().count('CG')

def calculate_methylation_batch(bw_file: str, regions_df: pd.DataFrame, 
                              fasta_file: str) -> pd.DataFrame:
    """
    Calculate methylation levels and coverage for multiple regions efficiently
    with proper signal normalization
    """
    bw = pyBigWig.open(bw_file)
    fasta = pysam.FastaFile(fasta_file)
    
    # Pre-filter valid chromosomes
    valid_chroms = set(fasta.references) & set(regions_df['chr'].unique())
    regions_df = regions_df[regions_df['chr'].isin(valid_chroms)].copy()
    
    # Vectorized calculations
    regions_df['sequence'] = regions_df.apply(
        lambda x: fasta.fetch(x['chr'], x['start'], x['end']), axis=1
    )
    regions_df['cpg_count'] = regions_df['sequence'].apply(count_cpgs_in_sequence)
    regions_df['region_length'] = regions_df['end'] - regions_df['start']
    regions_df['cpg_density'] = (regions_df['cpg_count'] * 100) / regions_df['region_length']
    
    # Calculate methylation and coverage in batches
    batch_size = 1000
    methylation_values = []
    coverage_values = []
    
    # Get global statistics for normalization
    all_values = []
    for _, row in regions_df.iterrows():
        vals = bw.values(row['chr'], row['start'], row['end'])
        if vals and any(v is not None for v in vals):
            all_values.extend([v for v in vals if v is not None])
    
    # Calculate normalization factors
    global_median = np.median(all_values)
    global_99th = np.percentile(all_values, 99)
    
    for i in range(0, len(regions_df), batch_size):
        batch = regions_df.iloc[i:i+batch_size]
        values = [
            bw.values(row['chr'], row['start'], row['end'])
            for _, row in batch.iterrows()
        ]
        
        # Calculate methylation and coverage for each region
        for vals in values:
            if vals is None or len(vals) == 0:
                methylation_values.append(0)
                coverage_values.append(0)
            else:
                valid_vals = [v for v in vals if v is not None]
                total_positions = len(vals)
                covered_positions = len(valid_vals)
                
                if covered_positions > 0:
                    # Normalize values to 0-100 range
                    normalized_vals = np.clip(
                        [100 * (v / global_99th) for v in valid_vals],
                        0, 100
                    )
                    methylation = np.mean(normalized_vals)
                else:
                    methylation = 0
                
                coverage = covered_positions / total_positions if total_positions > 0 else 0
                
                methylation_values.append(methylation)
                coverage_values.append(coverage)
    
    # Add calculated values to DataFrame
    regions_df['methylation'] = methylation_values
    regions_df['coverage'] = coverage_values
    
    # Cleanup
    bw.close()
    fasta.close()
    regions_df.drop('sequence', axis=1, inplace=True)
    
    return regions_df

def analyze_methylation_binding_relationship(mecp2_dir, mecp2_file, medip_dir):
    """
    Analyze relationship between MeCP2 binding and DNA methylation patterns
    
    Key analyses:
    1. Identify hypomethylated CpG islands bound by MeCP2
    2. Compare methylation levels between different binding categories
    3. Assess correlation between methylation and MeCP2 binding strength
    """
    mecp2_df = pd.read_csv(os.path.join(mecp2_dir, mecp2_file))
    
    cell_types = {
        'NEU': 'N',  # Neurons
        'NSC': 'PP'  # Neural Stem Cells
    }
    replicates = ['r1', 'r2', 'r3']
    
    results = []
    
    for mecp2_cell_type, medip_cell_type in cell_types.items():
        replicate_results = []
        
        # Process each replicate
        for rep in replicates:
            bw_file = os.path.join(medip_dir, f"Medip_{medip_cell_type}_output_{rep}.bw")
            meth_data = calculate_methylation_batch(
                bw_file=bw_file,
                regions_df=mecp2_df,
                fasta_file=genome_fasta
            )
            replicate_results.append(meth_data)
        
        # Combine replicate data
        temp_df = mecp2_df.copy()
        temp_df['methylation'] = np.mean([rep['methylation'].values for rep in replicate_results], axis=0)
        temp_df['coverage'] = np.mean([rep['coverage'].values for rep in replicate_results], axis=0)
        temp_df['quality'] = ['high' if cov > CONFIG['quality_thresholds']['coverage'] else 'low' 
                            for cov in temp_df['coverage']]
        temp_df['cpg_count'] = replicate_results[0]['cpg_count']
        temp_df['cpg_density'] = replicate_results[0]['cpg_density']
        
        # Classify methylation status using thresholds from CONFIG
        temp_df['methylation_status'] = pd.cut(
            temp_df['methylation'],
            bins=[0, CONFIG['methylation_thresholds']['hypo'],
                  CONFIG['methylation_thresholds']['hyper'], 100],
            labels=['hypomethylated', 'intermediate', 'hypermethylated']
        )
        
        temp_df['cell_type'] = mecp2_cell_type
        results.append(temp_df)
    
    final_df = pd.concat(results)
    
    # Analysis focusing on hypomethylation
    for mecp2_cell_type in cell_types.keys():
        high_quality_data = final_df[
            (final_df['cell_type'] == mecp2_cell_type) & 
            (final_df['quality'] == 'high')
        ]
        
        logger.info(f"\nAnalysis for {mecp2_cell_type} (high quality data only):")
        
        # Analyze methylation patterns by binding type
        for binding_type in ['exo_only', 'endo_only', 'both']:
            subset = high_quality_data[high_quality_data['binding_type'] == binding_type]
            
            # Calculate key metrics
            hypomethylated = (subset['methylation_status'] == 'hypomethylated').sum()
            total_sites = len(subset)
            
            if total_sites > 0:
                logger.info(f"\n{binding_type}:")
                logger.info(f"Total high-quality sites: {total_sites}")
                logger.info(f"Hypomethylated sites: {hypomethylated} ({100*hypomethylated/total_sites:.1f}%)")
                logger.info(f"Mean methylation: {subset['methylation'].mean():.2f}%")
                logger.info(f"Mean CpG density: {subset['cpg_density'].mean():.2f} CpGs/100bp")
        
        # Correlation analysis for sites with both binding types
        both_sites = high_quality_data[high_quality_data['binding_type'] == 'both']
        if len(both_sites) > 0:
            # Analyze correlation between methylation and binding strength
            corr_exo = stats.spearmanr(both_sites['methylation'], both_sites['exo_signal'])
            corr_endo = stats.spearmanr(both_sites['methylation'], both_sites['endo_signal'])
            
            logger.info("\nCorrelation Analysis (Spearman):")
            logger.info(f"Methylation vs Exogenous MeCP2: rho={corr_exo[0]:.3f}, p={corr_exo[1]:.3e}")
            logger.info(f"Methylation vs Endogenous MeCP2: rho={corr_endo[0]:.3f}, p={corr_endo[1]:.3e}")
    
    return final_df

def create_visualization_batch(data: pd.DataFrame, output_path: str, 
                             plot_type: str, **kwargs):
    """
    Create visualization with proper memory management
    """
    plt.figure(figsize=kwargs.get('figsize', (10, 6)))
    
    try:
        if plot_type == 'histogram':
            sns.histplot(data=data, **kwargs)
        elif plot_type == 'violin':
            sns.violinplot(data=data, **kwargs)
        elif plot_type == 'scatter':
            sns.scatterplot(data=data, **kwargs)
        
        plt.title(kwargs.get('title', ''))
        plt.xlabel(kwargs.get('xlabel', ''))
        plt.ylabel(kwargs.get('ylabel', ''))
        
        if kwargs.get('rotate_xticks'):
            plt.xticks(rotation=45)
        
        plt.tight_layout()
        plt.savefig(output_path)
    finally:
        plt.close()

def remove_outliers(df, column):
    Q1 = df[column].quantile(0.25)
    Q3 = df[column].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]

def create_comparison_plots(results_df, output_dir="plots"):
    plt.style.use('seaborn-v0_8')
    
    # Set a consistent color palette for binding types
    binding_colors = {
        'exo_only': '#FF9999',  # Light red
        'endo_only': '#66B2FF', # Light blue
        'both': '#99FF99'       # Light green
    }
    
    # First, create combined cell type comparisons
    for suffix, data in [('with_outliers', results_df), 
                        ('no_outliers', remove_outliers(results_df, 'methylation'))]:
        
        # 1. Combined box plot with cell types side by side
        plt.figure(figsize=(12, 7))
        sns.boxplot(data=data, x='binding_type', y='methylation', 
                   hue='cell_type', palette=['#FFA07A', '#98FB98'])
        plt.title(f'Methylation Levels Comparison - All Cell Types ({suffix})')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/methylation_boxplot_combined_{suffix}.pdf')
        plt.close()

        # 2. Combined violin plot
        plt.figure(figsize=(12, 7))
        sns.violinplot(data=data, x='binding_type', y='methylation',
                      hue='cell_type', palette=['#FFA07A', '#98FB98'],
                      split=True)
        plt.title(f'Methylation Distribution - All Cell Types ({suffix})')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/methylation_violin_combined_{suffix}.pdf')
        plt.close()

    # Then create individual cell type plots
    for cell_type in results_df['cell_type'].unique():
        cell_type_data = results_df[results_df['cell_type'] == cell_type]
        cell_type_data_no_outliers = remove_outliers(cell_type_data, 'methylation')
        
        for data, suffix in [(cell_type_data, 'with_outliers'), 
                            (cell_type_data_no_outliers, 'no_outliers')]:
            
            # 1. Enhanced box plot
            plt.figure(figsize=(10, 6))
            sns.boxplot(data=data, x='binding_type', y='methylation', 
                       hue='binding_type', palette=binding_colors, legend=False)
            plt.title(f'Methylation Levels - {cell_type} ({suffix})')
            plt.xticks(rotation=45)
            
            # Add statistical annotations
            means = data.groupby('binding_type')['methylation'].mean()
            for i, binding_type in enumerate(data['binding_type'].unique()):
                plt.text(i, means[binding_type], f'n={len(data[data.binding_type == binding_type])}',
                        ha='center', va='bottom')
            
            plt.tight_layout()
            plt.savefig(f'{output_dir}/methylation_boxplot_{cell_type}_{suffix}.pdf')
            plt.close()

            # 2. Enhanced violin plot with individual points
            plt.figure(figsize=(10, 6))
            sns.violinplot(data=data, x='binding_type', y='methylation',
                          hue='binding_type', palette=binding_colors, legend=False)
            # Replace swarmplot with stripplot
            sns.stripplot(data=data, x='binding_type', y='methylation',
                         color='black', alpha=0.5, size=2, jitter=True)
            plt.title(f'Methylation Distribution - {cell_type} ({suffix})')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(f'{output_dir}/methylation_violin_{cell_type}_{suffix}.pdf')
            plt.close()

            # 3. Enhanced density plot
            plt.figure(figsize=(10, 6))
            for bind_type in data['binding_type'].unique():
                subset = data[data['binding_type'] == bind_type]
                sns.kdeplot(data=subset, x='methylation', label=bind_type,
                           color=binding_colors[bind_type], fill=True, alpha=0.3)
            plt.title(f'Methylation Density - {cell_type} ({suffix})')
            plt.xlabel('Methylation Level')
            plt.ylabel('Density')
            plt.legend(title='Binding Type')
            plt.tight_layout()
            plt.savefig(f'{output_dir}/methylation_density_{cell_type}_{suffix}.pdf')
            plt.close()

        # 5. Statistical summary with additional metrics
        stats_summary = cell_type_data.groupby('binding_type').agg({
            'methylation': ['mean', 'std', 'median', 'count', 'min', 'max'],
            'exo_signal': ['mean', 'median'],
            'endo_signal': ['mean', 'median']
        }).round(3)
        
        stats_summary_no_outliers = cell_type_data_no_outliers.groupby('binding_type').agg({
            'methylation': ['mean', 'std', 'median', 'count', 'min', 'max'],
            'exo_signal': ['mean', 'median'],
            'endo_signal': ['mean', 'median']
        }).round(3)
        
        # Save statistics to CSV
        stats_summary.to_csv(f'{output_dir}/methylation_statistics_{cell_type}_with_outliers.csv')
        stats_summary_no_outliers.to_csv(f'{output_dir}/methylation_statistics_{cell_type}_no_outliers.csv')

    return {
        'with_outliers': stats_summary,
        'no_outliers': stats_summary_no_outliers
    }

def analyze_cpg_methylation(results_df, cpg_islands_file, medip_dir, genome_fasta, output_dir="plots"):
    """
    Compare methylation levels between MeCP2-bound and unbound CpG islands,
    with specific focus on demonstrating hypomethylation of MeCP2-bound regions
    """
    logger.info("Analyzing CpG island methylation patterns...")
    
    # Read and filter CpG islands
    cpg_islands = pd.read_csv(cpg_islands_file, sep='\t', 
                             names=['chr', 'start', 'end', 'id', 'type', 'cpg_count'])
    
    # Filter for standard chromosomes
    cpg_islands = cpg_islands[cpg_islands['chr'].isin(CONFIG['standard_chromosomes'])].copy()
    
    # Create bound/unbound classification
    bound_cpgs = set((row['chr'], row['start'], row['end']) for _, row in results_df.iterrows())
    cpg_islands['binding_status'] = cpg_islands.apply(
        lambda x: 'MeCP2-bound' if (x['chr'], x['start'], x['end']) in bound_cpgs 
        else 'Unbound', axis=1)
    
    # Calculate methylation levels
    cell_types = {'NEU': 'N', 'NSC': 'PP'}
    
    for cell_type, medip_prefix in cell_types.items():
        logger.info(f"Processing {cell_type}...")
        
        # Calculate methylation using parallel processing
        methylation_data = process_cell_type_parallel(
            cell_type, medip_prefix, cpg_islands, medip_dir, genome_fasta
        )
        
        # Add methylation data to cpg_islands
        cpg_islands[f'methylation_{cell_type}'] = methylation_data['methylation']
        
        # Classify methylation status using stricter threshold for hypomethylation
        cpg_islands[f'methylation_status_{cell_type}'] = pd.cut(
            cpg_islands[f'methylation_{cell_type}'],
            bins=[0, 20, 50, 100],  # Stricter threshold for hypomethylation
            labels=['hypomethylated', 'intermediate', 'hypermethylated']
        )
    
    # Create detailed methylation analysis plots
    create_cpg_methylation_plots(cpg_islands, output_dir)
    
    # Perform statistical analysis
    statistics = perform_cpg_methylation_analysis(cpg_islands)
    
    # Save detailed results
    save_cpg_analysis_results(cpg_islands, statistics, output_dir)
    
    return cpg_islands, statistics

def create_cpg_methylation_plots(data, output_dir):
    """
    Create comprehensive visualizations focusing on hypomethylation of MeCP2-bound regions
    """
    os.makedirs(f"{output_dir}/cpg_analysis", exist_ok=True)
    
    for cell_type in ['NEU', 'NSC']:
        # 1. Distribution of methylation levels
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=data, x='binding_status', y=f'methylation_{cell_type}')
        plt.title(f'{cell_type}: Methylation Levels in CpG Islands')
        plt.ylabel('Methylation Level (%)')
        plt.savefig(f'{output_dir}/cpg_analysis/methylation_distribution_{cell_type}.pdf')
        plt.close()
        
        # 2. Violin plot with individual points
        plt.figure(figsize=(10, 6))
        sns.violinplot(data=data, x='binding_status', y=f'methylation_{cell_type}')
        sns.stripplot(data=data, x='binding_status', y=f'methylation_{cell_type}',
                     color='black', alpha=0.3, size=2, jitter=True)
        plt.title(f'{cell_type}: Detailed Methylation Distribution')
        plt.ylabel('Methylation Level (%)')
        plt.savefig(f'{output_dir}/cpg_analysis/methylation_violin_{cell_type}.pdf')
        plt.close()
        
        # 3. Methylation status distribution
        status_data = pd.crosstab(
            data['binding_status'],
            data[f'methylation_status_{cell_type}'],
            normalize='index'
        ) * 100
        
        plt.figure(figsize=(10, 6))
        status_data.plot(kind='bar', stacked=True)
        plt.title(f'{cell_type}: Methylation Status Distribution')
        plt.xlabel('Binding Status')
        plt.ylabel('Percentage')
        plt.legend(title='Methylation Status')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/cpg_analysis/methylation_status_{cell_type}.pdf')
        plt.close()

def perform_cpg_methylation_analysis(data):
    """
    Perform statistical analysis to demonstrate hypomethylation of MeCP2-bound regions
    """
    statistics = {}
    
    for cell_type in ['NEU', 'NSC']:
        # 1. Basic statistics
        bound = data[data['binding_status'] == 'MeCP2-bound'][f'methylation_{cell_type}']
        unbound = data[data['binding_status'] == 'Unbound'][f'methylation_{cell_type}']
        
        # Mann-Whitney U test
        stat, pval = stats.mannwhitneyu(bound, unbound, alternative='less')  # Testing if bound < unbound
        
        # Effect size (Cohen's d)
        d = (bound.mean() - unbound.mean()) / np.sqrt((bound.std()**2 + unbound.std()**2) / 2)
        
        # 2. Hypomethylation enrichment
        bound_status = data[data['binding_status'] == 'MeCP2-bound'][f'methylation_status_{cell_type}']
        unbound_status = data[data['binding_status'] == 'Unbound'][f'methylation_status_{cell_type}']
        
        bound_hypo_pct = (bound_status == 'hypomethylated').mean() * 100
        unbound_hypo_pct = (unbound_status == 'hypomethylated').mean() * 100
        
        # Fisher's exact test for hypomethylation enrichment
        contingency = pd.crosstab(
            data['binding_status'],
            data[f'methylation_status_{cell_type}'] == 'hypomethylated'
        )
        fisher_stat, fisher_pval = stats.fisher_exact(contingency)
        
        statistics[cell_type] = {
            'mean_methylation': {
                'bound': bound.mean(),
                'unbound': unbound.mean()
            },
            'median_methylation': {
                'bound': bound.median(),
                'unbound': unbound.median()
            },
            'hypomethylation_percentage': {
                'bound': bound_hypo_pct,
                'unbound': unbound_hypo_pct
            },
            'statistical_tests': {
                'mannwhitney': {'statistic': stat, 'pvalue': pval},
                'effect_size': d,
                'fisher_exact': {'odds_ratio': fisher_stat, 'pvalue': fisher_pval}
            }
        }
    
    return statistics

def save_cpg_analysis_results(data, statistics, output_dir):
    """
    Save comprehensive analysis results
    """
    results_file = f'{output_dir}/cpg_analysis/methylation_analysis_results.txt'
    
    with open(results_file, 'w') as f:
        f.write("CpG Island Methylation Analysis Results\n")
        f.write("=====================================\n\n")
        
        for cell_type, stats in statistics.items():
            f.write(f"\n{cell_type} Analysis:\n")
            f.write("-----------------\n")
            
            # Mean methylation levels
            f.write("\nMean Methylation Levels:\n")
            f.write(f"MeCP2-bound: {stats['mean_methylation']['bound']:.2f}%\n")
            f.write(f"Unbound: {stats['mean_methylation']['unbound']:.2f}%\n")
            
            # Hypomethylation percentages
            f.write("\nHypomethylation Frequencies:\n")
            f.write(f"MeCP2-bound: {stats['hypomethylation_percentage']['bound']:.1f}%\n")
            f.write(f"Unbound: {stats['hypomethylation_percentage']['unbound']:.1f}%\n")
            
            # Statistical tests
            f.write("\nStatistical Analysis:\n")
            f.write("Mann-Whitney U test (testing if bound regions are hypomethylated):\n")
            f.write(f"p-value: {stats['statistical_tests']['mannwhitney']['pvalue']:.2e}\n")
            f.write(f"Effect size (Cohen's d): {stats['statistical_tests']['effect_size']:.3f}\n")
            
            f.write("\nFisher's exact test for hypomethylation enrichment:\n")
            f.write(f"Odds ratio: {stats['statistical_tests']['fisher_exact']['odds_ratio']:.3f}\n")
            f.write(f"p-value: {stats['statistical_tests']['fisher_exact']['pvalue']:.2e}\n")

def annotate_regions_with_genes(regions_df, gtf_file):
    """
    Annotate genomic regions with their nearest genes using GTF file
    
    Parameters:
    -----------
    regions_df : pandas DataFrame
        DataFrame containing genomic regions (chr, start, end)
    gtf_file : str
        Path to GTF file containing gene annotations
    
    Returns:
    --------
    pandas DataFrame
        Input DataFrame with additional 'nearest_gene' column
    """
    # Convert regions to PyRanges
    regions = pr.PyRanges(
        chromosomes=regions_df['chr'],
        starts=regions_df['start'],
        ends=regions_df['end']
    )
    
    # Read GTF file and extract genes
    genes = pr.read_gtf(gtf_file)
    genes = genes[genes.Feature == 'gene']
    
    # Find nearest genes
    nearest = regions.nearest(genes)
    
    # Add gene names to original DataFrame
    regions_df = regions_df.copy()
    regions_df['nearest_gene'] = nearest.gene_name
    
    return regions_df

def integrate_expression_data(methylation_df, rnaseq_files, gtf_file):
    """
    Integrate methylation data with RNA-seq expression data from DESeq2 results
    
    Parameters:
    -----------
    methylation_df : pandas DataFrame
        DataFrame containing methylation and MeCP2 binding data
    rnaseq_files : dict
        Dictionary with cell types as keys and paths to DESeq2 results as values
    gtf_file : str
        Path to GTF file for gene annotations
    """
    # First, annotate regions with nearest genes
    print("Annotating regions with nearest genes...")
    methylation_df = annotate_regions_with_genes(methylation_df, gtf_file)
    
    results = []
    
    for cell_type, rnaseq_file in rnaseq_files.items():
        print(f"\nProcessing {cell_type}...")
        
        # Load RNA-seq data for this cell type
        expression_df = pd.read_csv(rnaseq_file)
        
        # Filter methylation data for this cell type
        cell_methylation = methylation_df[methylation_df['cell_type'] == cell_type].copy()
        
        # Merge methylation and expression data
        integrated = cell_methylation.merge(
            expression_df,
            left_on='nearest_gene',
            right_on='gene',
            how='inner'
        )
        
        # Classify genes based on expression changes and significance
        integrated['expression_change'] = 'unchanged'
        
        # Significant upregulation: log2FC > 1 and padj < 0.05
        integrated.loc[(integrated['log2FoldChange'] > 1) & 
                      (integrated['padj'] < 0.05), 'expression_change'] = 'upregulated'
        
        # Significant downregulation: log2FC < -1 and padj < 0.05
        integrated.loc[(integrated['log2FoldChange'] < -1) & 
                      (integrated['padj'] < 0.05), 'expression_change'] = 'downregulated'
        
        results.append(integrated)
        
    if not results:
        raise ValueError("No data could be integrated. Check gene annotations.")
    
    final_df = pd.concat(results)
    return final_df

def create_methylation_expression_plots(integrated_data, output_dir="plots"):
    """
    Create visualizations focusing on the relationship between methylation,
    MeCP2 binding, and gene expression changes.
    """
    plt.style.use('seaborn-v0_8')
    
    # 1. Expression changes in hypomethylated vs other regions
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=integrated_data, 
                x='methylation_status', 
                y='log2FoldChange',
                hue='binding_type')
    plt.title('Expression Changes by Methylation Status and Binding Type')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/expression_by_methylation_status.pdf')
    plt.close()
    
    # 2. Methylation levels in up/down regulated genes
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=integrated_data,
                x='expression_change',
                y='methylation',
                hue='binding_type')
    plt.title('Methylation Levels by Expression Change and Binding Type')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/methylation_by_expression_change.pdf')
    plt.close()
    
    # 3. Scatter plot: Methylation vs Expression Change
    for cell_type in integrated_data['cell_type'].unique():
        cell_data = integrated_data[integrated_data['cell_type'] == cell_type]
        
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=cell_data,
                       x='methylation',
                       y='log2FoldChange',
                       hue='binding_type',
                       style='methylation_status',
                       alpha=0.6)
        plt.title(f'{cell_type}: Methylation vs Expression Change')
        plt.xlabel('Methylation Level (%)')
        plt.ylabel('log2 Fold Change')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/methylation_vs_expression_{cell_type}.pdf')
        plt.close()
        
        # 4. Focused analysis of hypomethylated regions
        hypo_data = cell_data[cell_data['methylation_status'] == 'hypomethylated']
        
        plt.figure(figsize=(10, 6))
        sns.violinplot(data=hypo_data,
                      x='binding_type',
                      y='log2FoldChange',
                      hue='expression_change')
        plt.title(f'{cell_type}: Expression Changes in Hypomethylated Regions')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/hypomethylated_expression_{cell_type}.pdf')
        plt.close()

def perform_statistical_tests(integrated_data):
    """
    Perform statistical tests to analyze the significance of observed patterns
    """
    results = {}
    
    for cell_type in integrated_data['cell_type'].unique():
        cell_data = integrated_data[integrated_data['cell_type'] == cell_type]
        
        # Test 1: Are hypomethylated regions enriched for expression changes?
        hypo = cell_data[cell_data['methylation_status'] == 'hypomethylated']
        other = cell_data[cell_data['methylation_status'] != 'hypomethylated']
        
        stat, pval = stats.mannwhitneyu(
            hypo['log2FoldChange'].abs(),
            other['log2FoldChange'].abs()
        )
        
        results[f'{cell_type}_hypo_enrichment'] = {
            'statistic': stat,
            'pvalue': pval
        }
        
        # Test 2: Different expression patterns in bound vs unbound hypomethylated regions
        hypo_bound = hypo[hypo['binding_type'].isin(['exo_only', 'both'])]
        hypo_unbound = hypo[hypo['binding_type'] == 'endo_only']
        
        if len(hypo_bound) > 0 and len(hypo_unbound) > 0:
            stat, pval = stats.mannwhitneyu(
                hypo_bound['log2FoldChange'],
                hypo_unbound['log2FoldChange']
            )
            
            results[f'{cell_type}_binding_effect'] = {
                'statistic': stat,
                'pvalue': pval
            }
    
    return results

def run_complete_analysis():
    """
    Run complete analysis with error handling and progress tracking
    """
    try:
        logger.info("Starting analysis pipeline...")
        
        # Validate required files exist
        required_files = {
            'genome_fasta': genome_fasta,
            'mecp2_file': os.path.join(mecp2_dir, mecp2_file),
            'cpg_islands_file': cpg_islands_file,
            'gtf_file': gtf_file
        }
        
        for name, filepath in required_files.items():
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"Required file not found: {name} ({filepath})")
        
        # Create output directories
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(f"{output_dir}/with_outliers", exist_ok=True)
        os.makedirs(f"{output_dir}/no_outliers", exist_ok=True)
        
        # Define analysis steps with their required arguments
        steps = [
            (
                "Analyzing methylation and MeCP2 binding",
                analyze_methylation_binding_relationship,
                {
                    'mecp2_dir': mecp2_dir,
                    'mecp2_file': mecp2_file,
                    'medip_dir': medip_dir
                }
            ),
            (
                "Analyzing CpG islands",
                analyze_cpg_methylation,
                {
                    'results_df': 'methylation_results',  # Will be updated with actual results
                    'cpg_islands_file': cpg_islands_file,
                    'medip_dir': medip_dir,
                    'genome_fasta': genome_fasta,
                    'output_dir': output_dir
                }
            ),
            (
                "Integrating RNA-seq data",
                integrate_expression_data,
                {
                    'methylation_df': 'methylation_results',  # Will be updated with actual results
                    'rnaseq_files': rnaseq_files,
                    'gtf_file': gtf_file
                }
            ),
            (
                "Creating visualizations",
                create_comprehensive_visualizations,
                {
                    'methylation_results': 'methylation_results',  # Will be updated
                    'cpg_islands': 'cpg_islands',  # Will be updated
                    'integrated_data': 'integrated_data',  # Will be updated
                    'output_dir': output_dir
                }
            ),
            (
                "Performing statistical analysis",
                perform_comprehensive_analysis,
                {
                    'methylation_results': 'methylation_results',  # Will be updated
                    'cpg_islands': 'cpg_islands',  # Will be updated
                    'integrated_data': 'integrated_data',  # Will be updated
                    'output_dir': output_dir
                }
            )
        ]
        
        results = {}
        for step_name, step_func, step_args in steps:
            logger.info(f"Starting {step_name}...")
            try:
                # Update arguments with actual results from previous steps
                current_args = {}
                for arg_name, arg_value in step_args.items():
                    if isinstance(arg_value, str) and arg_value in results:
                        current_args[arg_name] = results[arg_value]
                    else:
                        current_args[arg_name] = arg_value
                
                # Run the analysis step
                step_result = step_func(**current_args)
                
                # Store results
                if step_func == analyze_methylation_binding_relationship:
                    results['methylation_results'] = step_result
                elif step_func == analyze_cpg_methylation:
                    results['cpg_islands'], results['cpg_stats'] = step_result
                elif step_func == integrate_expression_data:
                    results['integrated_data'] = step_result
                
                logger.info(f"Completed {step_name}")
                
            except Exception as e:
                logger.error(f"Error in {step_name}: {str(e)}")
                raise
        
        logger.info("Analysis pipeline completed successfully")
        return results
        
    except Exception as e:
        logger.error(f"Analysis pipeline failed: {str(e)}")
        raise

def create_comprehensive_visualizations(methylation_results, cpg_islands, integrated_data, output_dir):
    """
    Create comprehensive set of biologically relevant visualizations
    """
    # Create subdirectories
    os.makedirs(f"{output_dir}/methylation", exist_ok=True)
    os.makedirs(f"{output_dir}/expression", exist_ok=True)
    os.makedirs(f"{output_dir}/integrated", exist_ok=True)
    
    # Set style
    plt.style.use('default')  # Use default matplotlib style
    sns.set_theme()  # Apply seaborn defaults
    
    # Define color palettes
    binding_colors = {
        'exo_only': '#1f77b4',
        'endo_only': '#ff7f0e',
        'both': '#2ca02c'
    }
    
    cpg_colors = {
        'MeCP2-bound': '#1f77b4',
        'Unbound': '#d62728'
    }
    
    # Set default figure parameters
    plt.rcParams.update({
        'figure.figsize': (10, 6),
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight'
    })
    
    # 1. Methylation distribution plots
    for data_version, data in [
        ('with_outliers', methylation_results),
        ('no_outliers', remove_outliers(methylation_results, 'methylation'))
    ]:
        # Global methylation distribution
        plt.figure(figsize=(10, 6))
        sns.histplot(data=data, x='methylation', hue='cell_type', bins=50)
        plt.title(f'Global Methylation Distribution ({data_version})')
        plt.xlabel('Methylation Level (%)')
        plt.ylabel('Count')
        plt.savefig(f'{output_dir}/methylation/global_distribution_{data_version}.pdf')
        plt.close()
        
        # Methylation by binding type
        for cell_type in data['cell_type'].unique():
            cell_data = data[data['cell_type'] == cell_type]
            
            plt.figure(figsize=(12, 6))
            sns.violinplot(data=cell_data, x='binding_type', y='methylation',
                          hue='binding_type', palette=binding_colors, legend=False)
            plt.title(f'{cell_type}: Methylation Levels by Binding Type ({data_version})')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(f'{output_dir}/methylation/binding_distribution_{cell_type}_{data_version}.pdf')
            plt.close()
            
            # CpG density vs methylation scatter
            plt.figure(figsize=(10, 6))
            for binding_type in cell_data['binding_type'].unique():
                subset = cell_data[cell_data['binding_type'] == binding_type]
                plt.scatter(subset['cpg_density'], subset['methylation'],
                          label=binding_type, alpha=0.6,
                          color=binding_colors[binding_type])
            plt.title(f'{cell_type}: CpG Density vs Methylation ({data_version})')
            plt.xlabel('CpG Density (CpGs/100bp)')
            plt.ylabel('Methylation Level (%)')
            plt.legend()
            plt.tight_layout()
            plt.savefig(f'{output_dir}/methylation/cpg_density_{cell_type}_{data_version}.pdf')
            plt.close()
    
    # 2. Expression-related plots
    for cell_type in integrated_data['cell_type'].unique():
        cell_data = integrated_data[integrated_data['cell_type'] == cell_type]
        
        # Expression vs methylation scatter
        plt.figure(figsize=(12, 8))
        for binding_type in cell_data['binding_type'].unique():
            subset = cell_data[cell_data['binding_type'] == binding_type]
            plt.scatter(subset['methylation'], subset['log2FoldChange'],
                       label=binding_type, alpha=0.6,
                       color=binding_colors[binding_type])
        
        plt.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        plt.axvline(x=30, color='gray', linestyle='--', alpha=0.5)  # Hypomethylation threshold
        
        plt.title(f'{cell_type}: Expression vs Methylation')
        plt.xlabel('Methylation Level (%)')
        plt.ylabel('log2 Fold Change')
        plt.legend(title='Binding Type')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/expression/methylation_vs_expression_{cell_type}.pdf')
        plt.close()
        
        # Expression in hypomethylated regions
        hypo_data = cell_data[cell_data['methylation_status'] == 'hypomethylated']
        
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=hypo_data, x='binding_type', y='log2FoldChange',
                   hue='binding_type', palette=binding_colors, legend=False)
        plt.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        plt.title(f'{cell_type}: Expression Changes in Hypomethylated Regions')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/expression/hypomethylated_expression_{cell_type}.pdf')
        plt.close()
    
    # 3. CpG island specific analysis
    if cpg_islands is not None:
        # Methylation distribution in CpG islands for each cell type
        cell_types = ['NEU', 'NSC']
        for cell_type in cell_types:
            plt.figure(figsize=(10, 6))
            sns.boxplot(data=cpg_islands, x='binding_status', 
                       y=f'methylation_{cell_type}',
                       hue='binding_status', palette=cpg_colors, legend=True)
            plt.title(f'{cell_type}: Methylation in CpG Islands')
            plt.xticks(rotation=45)
            plt.ylabel('Methylation Level (%)')
            plt.legend(title='Binding Status', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(f'{output_dir}/integrated/cpg_island_methylation_{cell_type}.pdf')
            plt.close()
    
    logger.info(f"Created visualizations in {output_dir}")

def perform_comprehensive_analysis(methylation_results, cpg_islands, integrated_data, output_dir):
    """
    Perform comprehensive statistical analysis with biological interpretation
    """
    results = {}
    
    # 1. Analyze methylation patterns
    for cell_type in methylation_results['cell_type'].unique():
        cell_data = methylation_results[methylation_results['cell_type'] == cell_type]
        
        # Methylation status distribution
        meth_dist = cell_data['methylation_status'].value_counts()
        
        # Binding type enrichment in hypomethylated regions
        hypo = cell_data[cell_data['methylation_status'] == 'hypomethylated']
        binding_dist = hypo['binding_type'].value_counts()
        
        results[f'{cell_type}_methylation'] = {
            'methylation_distribution': meth_dist.to_dict(),
            'hypomethylated_binding': binding_dist.to_dict()
        }
    
    # 2. Analyze expression patterns
    for cell_type in integrated_data['cell_type'].unique():
        cell_data = integrated_data[integrated_data['cell_type'] == cell_type]
        
        # Expression changes in hypomethylated regions
        hypo_data = cell_data[cell_data['methylation_status'] == 'hypomethylated']
        expr_dist = hypo_data['expression_change'].value_counts()
        
        # Correlation between methylation and expression
        corr = stats.spearmanr(cell_data['methylation'], cell_data['log2FoldChange'])
        
        results[f'{cell_type}_expression'] = {
            'expression_distribution': expr_dist.to_dict(),
            'methylation_expression_correlation': {
                'rho': corr.statistic,
                'pvalue': corr.pvalue
            }
        }
    
    # Save results
    with open(f'{output_dir}/analysis_results.txt', 'w') as f:
        for key, value in results.items():
            f.write(f"\n{key}:\n")
            f.write(json.dumps(value, indent=2))
            f.write("\n")
    
    return results

def process_cell_type_parallel(cell_type: str, medip_prefix: str, 
                             mecp2_df: pd.DataFrame, medip_dir: str, 
                             genome_fasta: str) -> pd.DataFrame:
    """
    Process methylation data for a cell type in parallel
    """
    replicates = ['r1', 'r2', 'r3']
    
    # Create partial function for parallel processing
    process_replicate = partial(
        calculate_methylation_batch,
        regions_df=mecp2_df,
        fasta_file=genome_fasta
    )
    
    # Process replicates in parallel
    bw_files = [
        os.path.join(medip_dir, f"Medip_{medip_prefix}_output_{rep}.bw")
        for rep in replicates
    ]
    
    with ProcessPoolExecutor() as executor:
        replicate_results = list(executor.map(process_replicate, bw_files))
    
    # Combine results
    temp_df = mecp2_df.copy()
    temp_df['methylation'] = np.mean([
        rep['methylation'].values for rep in replicate_results
    ], axis=0)
    temp_df['coverage'] = np.mean([
        rep['coverage'].values for rep in replicate_results
    ], axis=0)
    
    return temp_df

# Define input paths
mecp2_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment_align2_005/NSC/mecp2_cpg_enrichment_parallel"
mecp2_file = "mecp2_cpg_enrichment_parallel.csv"
medip_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeDIP/output_old/bigwig"
cpg_islands_file = "../DATA/cpg_islands.bed"
rnaseq_files = {'NEU': '../iterative_alternative/DATA/DEA_NEU_filtered.csv', 'NSC': '../iterative_alternative/DATA/DEA_NSC_filtered.csv'}
genome_fasta = "../DATA/mm10.fa"
output_dir = "plots"
gtf_file = "../DATA/gencode.vM10.annotation.gtf"

#set working directory
os.chdir("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Methylation")

if __name__ == "__main__":
    try:
        # Set working directory (although this is now handled in the shell script)
        os.chdir("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Methylation")
        
        # Run the analysis
        results = run_complete_analysis()
        
        # Save final summary
        with open(f"{output_dir}/analysis_summary.txt", "w") as f:
            f.write("Analysis Summary\n")
            f.write("===============\n\n")
            
            if 'methylation_results' in results:
                f.write("Methylation Analysis:\n")
                f.write(f"Total regions analyzed: {len(results['methylation_results'])}\n")
                
            if 'cpg_islands' in results:
                f.write("\nCpG Islands Analysis:\n")
                f.write(f"Total CpG islands: {len(results['cpg_islands'])}\n")
                
            if 'integrated_data' in results:
                f.write("\nIntegrated Analysis:\n")
                f.write(f"Total integrated regions: {len(results['integrated_data'])}\n")
                
    except Exception as e:
        logger.error(f"Program failed: {str(e)}")
        raise