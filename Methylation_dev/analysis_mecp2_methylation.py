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
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from tqdm import tqdm
import time
import argparse

# Set working directory
os.chdir('/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Methylation_dev')


import importlib
import functions
import config
import cache_utils

importlib.reload(functions)
importlib.reload(config) 
importlib.reload(cache_utils)

from functions import *
from config import *
from cache_utils import *

#%% Define the analysis pipeline function
def run_analysis_pipeline(force_recompute: bool = False, n_processes: int = None,
                        cell_type: str = None, chromosome: str = None) -> Dict[str, Any]:
    """Run the analysis pipeline with support for cell type and chromosome filtering
    
    This function orchestrates the full analysis pipeline:
    1. Loads MeCP2 binding data from ChIP-seq
    2. Loads gene annotations and creates gene name mapping
    3. Loads RNA-seq expression data for specified cell types
    4. Validates and merges the datasets
    5. Calculates methylation levels in promoter and gene body regions
    6. Analyzes TSS binding patterns
    
    Args:
        force_recompute: If True, clears cache and recomputes all results
        n_processes: Number of processes to use for parallel computation
        cell_type: Optional cell type to analyze (NEU or NSC)
        chromosome: Optional chromosome to analyze
        
    Returns:
        Dictionary containing:
        - methylation_results: Methylation analysis for each cell type
        - tss_results: TSS binding analysis results
        - genes_df: Gene annotations
        - expression_data: RNA-seq expression data
        - mecp2_binding: MeCP2 ChIP-seq binding data
    """
    if force_recompute:
        clear_cache()
        logger.info("Cache cleared due to force recompute flag")
    
    # Load MeCP2 binding data from ChIP-seq
    mecp2_binding = pd.read_csv(os.path.join(PATHS['mecp2_dir'], PATHS['mecp2_file']))
    print(mecp2_binding.head())

    # Load gene annotations and create gene name to ID mapping
    genes_df, gene_name_to_id = load_gene_annotations(PATHS['gtf_file'])
    
    # Filter by chromosome if specified
    if chromosome:
        mecp2_binding = mecp2_binding[mecp2_binding['chr'] == chromosome]
        genes_df = genes_df[genes_df['chr'] == chromosome]
        logger.info(f"Filtered data for chromosome {chromosome}")
    
    # Load RNA-seq expression data for specified cell type(s)
    if cell_type:
        expression_data = {
            cell_type: load_expression_data(PATHS['rnaseq'][cell_type], gene_name_to_id, single_file=True)
        }
        logger.info(f"Analyzing cell type: {cell_type}")
    else:
        expression_data = load_expression_data(PATHS['rnaseq'], gene_name_to_id)
    
    # Print mapping statistics for debugging
    for cell_type, expr_df in expression_data.items():
        print(f"\n{cell_type} statistics:")
        print(f"Total genes: {len(expr_df)}")
        print(f"Genes with valid mapping: {expr_df['gene_id'].notna().sum()}")
        print("\nSample of mapped genes:")
        sample_df = expr_df[expr_df['gene_id'].notna()].head()
        print(pd.DataFrame({
            'Gene Symbol': sample_df['gene'],
            'Ensembl ID': sample_df['gene_id'],
            'log2FC': sample_df['log2FoldChange'],
            'padj': sample_df['padj']
        }))

    # Run analysis for each cell type
    results = {}
    tss_results = {}
    
    for ct, expr_df in expression_data.items():
        # Merge gene annotations, expression data and binding data
        merged_df = validate_and_merge_data(genes_df, expr_df, mecp2_binding)
        
        # Calculate methylation levels in promoter and gene body regions
        results[ct] = calculate_contextual_methylation(
            merged_df,
            PATHS['medip_dir'],
            ct,
            PATHS['genome_fasta'],
            n_processes=n_processes
        )
        
        # Analyze binding patterns around transcription start sites
        logger.info(f"Performing TSS binding analysis for {ct}...")
        tss_results[ct] = analyze_tss_binding_patterns(results[ct], ct, PATHS['output_dir'])
    
    return {
        'methylation_results': results,
        'tss_results': tss_results,
        'genes_df': genes_df,
        'expression_data': expression_data,
        'mecp2_binding': mecp2_binding
    }

#%% Configure debug mode and parse command line arguments
try:
    # For Jupyter notebook
    get_ipython()
    in_notebook = True
except:
    in_notebook = False

#%% Configure debug mode and parse command line arguments
if not in_notebook:
    # Command line mode
    parser = argparse.ArgumentParser(description='Analyze MeCP2 methylation patterns')
    parser.add_argument('--debug', action='store_true', help='Run in debug mode with minimal dataset')
    parser.add_argument('--sample-size', type=int, default=100, help='Number of genes to use in debug mode')
    parser.add_argument('--force-recompute', action='store_true', 
                       help='Force recomputation of all analysis stages')
    parser.add_argument('--experiment', type=str, default='align2_005',
                       help='Experiment name to analyze')
    parser.add_argument('--processes', type=int, default=None,
                       help='Number of processes to use for parallel computation')
    parser.add_argument('--cell-type', type=str, choices=['NEU', 'NSC'], 
                       help='Cell type to analyze')
    parser.add_argument('--chromosome', type=str, 
                       help='Chromosome to analyze (e.g., chr1)')
    args = parser.parse_args()
else:
    # Default arguments for notebook
    class Args:
        def __init__(self):
            self.debug = False
            self.sample_size = 100
            self.force_recompute = True
            self.experiment = 'align2_005'
            self.processes = None
            self.cell_type = 'NSC'
            self.chromosome = None
    args = Args()

# Set experiment name in environment variable for config.py to use
os.environ['EXPERIMENT_NAME'] = args.experiment

# Reload config to update paths with new experiment name
import importlib
import config
importlib.reload(config)
from config import PATHS, CONFIG, logger

# Configure debug mode if specified
if args.debug:
    logger.info(f"Running in DEBUG mode with sample size: {args.sample_size}")
    CONFIG['debug'] = {
        'enabled': True,
        'sample_size': args.sample_size
    }

#%% Initialize output directory
os.makedirs(PATHS['output_dir'], exist_ok=True)

#%% Run analysis pipeline
# analysis_results = run_analysis_pipeline(
#     force_recompute=args.force_recompute,
#     n_processes=args.processes,
#     cell_type=args.cell_type,
#     chromosome=args.chromosome
# )
# results = analysis_results['methylation_results']


### Cells based execution ########################################################################################################%% Initialize parameters
#%% Initialize parameters
force_recompute = True
n_processes = None
cell_type = 'NSC'  # or None to analyze all cell types
chromosome = None   # or specify like 'chr1'

if force_recompute:
    clear_cache()
    logger.info("Cache cleared due to force recompute flag")

#%% Load MeCP2 binding data
mecp2_binding = pd.read_csv(os.path.join(PATHS['mecp2_dir'], PATHS['mecp2_file']))
print("MeCP2 binding data:")
print(mecp2_binding.head())

#%% Load gene annotations
genes_df, gene_name_to_id = load_gene_annotations(PATHS['gtf_file'])

# Filter by chromosome if specified
if chromosome:
    mecp2_binding = mecp2_binding[mecp2_binding['chr'] == chromosome]
    genes_df = genes_df[genes_df['chr'] == chromosome]
    logger.info(f"Filtered data for chromosome {chromosome}")

#%% Load RNA-seq expression data
if cell_type:
    expression_data = {
        cell_type: load_expression_data(PATHS['rnaseq'][cell_type], gene_name_to_id, single_file=True)
    }
    logger.info(f"Analyzing cell type: {cell_type}")
else:
    expression_data = load_expression_data(PATHS['rnaseq'], gene_name_to_id)


#%% Print mapping statistics
for ct, expr_df in expression_data.items():
    print(f"\n{ct} statistics:")
    print(f"Total genes: {len(expr_df)}")
    print(f"Genes with valid mapping: {expr_df['gene_id'].notna().sum()}")
    print("\nSample of mapped genes:")
    sample_df = expr_df[expr_df['gene_id'].notna()].head()
    print(pd.DataFrame({
        'Gene Symbol': sample_df['gene'],
        'Ensembl ID': sample_df['gene_id'],
        'log2FC': sample_df['log2FoldChange'],
        'padj': sample_df['padj']
    }))

#%% Run analysis for each cell type
results = {}
tss_results = {}

for ct, expr_df in expression_data.items():
    # Merge datasets
    merged_df = validate_and_merge_data(genes_df, expr_df, mecp2_binding)
    
    # Calculate methylation levels
    results[ct] = calculate_contextual_methylation(
        merged_df,
        PATHS['medip_dir'],
        ct,
        PATHS['genome_fasta'],
        n_processes=n_processes
    )
    
    # Analyze TSS binding patterns
    logger.info(f"Performing TSS binding analysis for {ct}...")
    # Creates:
    # - {output_dir}/{cell_type}_tss_binding_patterns.pdf (Plot showing TSS binding patterns)
    # - {output_dir}/{cell_type}_tss_analysis.txt (Text file with binding statistics)
    # - {output_dir}/tss_binding/{cell_type}_exo_enriched.csv (CSV with exo-enriched genes)
    # - {output_dir}/tss_binding/{cell_type}_endo_only.csv (CSV with endo-only genes) 
    # - {output_dir}/tss_binding/{cell_type}_non_enriched.csv (CSV with non-enriched genes)
    tss_results[ct] = analyze_tss_binding_patterns(results[ct], ct, PATHS['output_dir'])

#%%
# Print all columns without wrapping
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

merged_df.head()

#%%
print(tss_results['NSC'].keys())
tss_results['NSC']['exo_enriched']

#%% Collect all results
analysis_results = {
    'methylation_results': results,
    'tss_results': tss_results,
    'genes_df': genes_df,
    'expression_data': expression_data,
    'mecp2_binding': mecp2_binding
}

##########################################################################################################

#%%
print(analysis_results.keys())

#%%
print(results.keys())

#%%
print(results['NSC'].head())

#%% Save results to pickle file for later use
logger.info("Saving full results dictionary to pickle file...")
import pickle
with open(f"{PATHS['output_dir']}/full_results.pkl", 'wb') as f:
    pickle.dump(results, f)
logger.info("Results saved successfully")

#%% Code for loading previously saved results (commented out)
# logger.info("Loading full results dictionary from pickle file...")
# with open(f"{PATHS['output_dir']}/full_results.pkl", 'rb') as f:
#     results = pickle.load(f)
# logger.info("Results loaded successfully")

### Plots ########################################################################################################%% Initialize parameters

#%% Generate visualizations of methylation patterns
# Creates the following files/directories:
# - {output_dir}/methylation_plots/
#   - {cell_type}_methylation_expression.pdf - Scatterplots of methylation vs expression
#   - {cell_type}_methylation_vs_binding.pdf - Scatterplots of methylation vs binding signal
#   - {cell_type}_methylation_heatmap.pdf - Heatmap of methylation patterns
#   - {cell_type}_binding_proportions.pdf - Bar plots of binding type distributions
create_methylation_plots(results, PATHS['output_dir'])

#%% Perform statistical analysis on methylation data
# Creates the following files/directories:
# - No files or directories created directly - this function only performs statistical tests
#   and returns a dictionary of results that are later saved in statistical_analysis.txt
stats_results = perform_statistical_analysis(results)

#%%
print(stats_results.keys())
print(stats_results['NSC'])

#%% Save detailed results to CSV files
for cell_type, df in results.items():
    df.to_csv(f"{PATHS['output_dir']}/{cell_type}_methylation_analysis.csv", 
                index=False)

#%% Save statistical analysis results
with open(f"{PATHS['output_dir']}/statistical_analysis.txt", 'w') as f:
    for cell_type, stats in stats_results.items():
        f.write(f"\n{cell_type} Analysis:\n")
        f.write("="*20 + "\n")
        for region, results in stats.items():
            f.write(f"\n{region.upper()}:\n")
            f.write("-"*20 + "\n")
            
            for test_name, values in results.items():
                f.write(f"{test_name}:\n")
                f.write(f"Statistic: {values['statistic']}\n")
                f.write(f"P-value: {values['pvalue']}\n") 

#%% Load saved results for additional analysis
import pickle
logger.info("Loading full results dictionary from pickle file...")
with open(f"{PATHS['output_dir']}/full_results.pkl", 'rb') as f:
    results = pickle.load(f)
logger.info("Results loaded successfully")

#%% Analyze genes regulated by MeCP2
# Creates the following files/directories:
# - {output_dir}/mecp2_regulated_analysis/
#   - {cell_type}_regulated_stats.txt - Statistical analysis results
#   - {cell_type}_cpg_density.pdf - CpG density plots
#   - {cell_type}_regulated_genes.csv - Lists of regulated genes
create_mecp2_regulated_analysis(results, PATHS['output_dir'])

#%% Debug gene filtering process
# Creates the following files/directories:
# - {output_dir}/debug_logs/
#   - {cell_type}_filtering_steps.txt - Detailed logs of each filtering step
#   - {cell_type}_filtered_genes.csv - Lists of genes at each filtering stage
debug_regulated_genes_filtering(results)

#%% Print improved summary of regulated genes
# Creates the following files/directories:
# - {output_dir}/regulated_summary/
#   - {cell_type}_regulated_summary.txt - Detailed summary of regulated genes by binding type
#   - {cell_type}_regulated_counts.txt - Gene counts by regulation and binding category
print_regulated_summary_v2(results)

#%% Generate focused visualization of significant differences
# Creates the following files/directories:
# - {output_dir}/methylation_plots/
#   - {cell_type}_methylation_differences.pdf - Violin plots showing methylation differences
#   - {cell_type}_methylation_stats.txt - Statistical test results for methylation differences
plot_significant_differences(results, PATHS['output_dir'])

#%% Generate comprehensive visualization of MeCP2 regulatory patterns
# Creates the following files/directories:
# - {output_dir}/regulatory_patterns/
#   - {cell_type}_binding_correlations.pdf - Plots showing correlations between binding strength and methylation
#   - {cell_type}_distributions.pdf - Distribution plots for methylation and binding signals
plot_mecp2_regulatory_patterns(results, PATHS['output_dir'])

#%% Perform detailed methylation analysis
logger.info("Performing detailed methylation analysis...")
try:
    # Creates the following files/directories:
    # - {output_dir}/methylation_patterns/
    #   - {cell_type}_methylation_stats.txt - Statistical analysis of methylation patterns
    #   - {cell_type}_methylation_heatmap.pdf - Heatmap showing methylation patterns
    #   - {cell_type}_correlation_plots.pdf - Correlation plots between methylation and expression
    #   - {cell_type}_distribution_plots.pdf - Distribution plots of methylation levels
    analyze_methylation_patterns_detailed(results, PATHS['output_dir'])
    logger.info("Completed detailed methylation analysis")
except Exception as e:
    logger.error(f"Error in detailed methylation analysis: {str(e)}")

#%% Analyze binding enrichment groups
binding_analysis = {}
for cell_type, df in results.items():
    logger.info(f"Analyzing binding enrichment for {cell_type}")
    cell_output_dir = os.path.join(PATHS['output_dir'], cell_type)
    os.makedirs(cell_output_dir, exist_ok=True)
    # Creates the following files/directories:
    # - {cell_output_dir}/
    #   - {group_name}_analysis.txt - Detailed analysis results for each binding group (exo_enriched, endo_only, non_enriched)
    #   - {group_name}_{status}_genes.csv - Gene lists for each binding group and regulation status combination
    binding_analysis[cell_type] = analyze_binding_enrichment_groups(df, cell_output_dir)

# Print summary of results
for cell_type, cell_results in binding_analysis.items():
    print(f"\nResults for {cell_type}:")
    for group_name, group_data in cell_results.items():
        print(f"\n{group_name}:")
        print("Regulation counts:")
        print(group_data['regulation_counts'])

#%% Special analysis for exo-enriched upregulated genes
logger.info("Performing special analysis of exo-enriched upregulated genes...")
# Creates the following files/directories:
# - {output_dir}/exo_upregulated/
#   - {cell_type}_exo_up_methylation.pdf - Methylation patterns in exo-enriched upregulated genes
#   - {cell_type}_exo_up_stats.txt - Statistical analysis of exo-enriched upregulated genes
#   - {cell_type}_exo_up_genes.csv - List of exo-enriched upregulated genes with details
analyze_exo_upregulated(results, PATHS['output_dir'])

#%% Perform detailed binding enrichment analysis
logger.info("Performing detailed binding enrichment analysis...")
# Creates the following files/directories:
# - {output_dir}/binding_enrichment/
#   - {cell_type}_binding_enrichment_stats.txt - Statistical analysis of binding enrichment patterns
#   - {cell_type}_binding_enrichment_plots.pdf - Plots showing binding enrichment patterns
#   - {cell_type}_binding_enrichment_heatmap.pdf - Heatmap of binding enrichment
#   - {cell_type}_binding_enrichment_genes.csv - List of genes with binding enrichment details
analyze_binding_enrichment_detailed(results, PATHS['output_dir'])

#%% Create additional visualizations
logger.info("Creating additional visualizations...")
for cell_type, df in results.items():
    # Creates the following files/directories:
    # - {output_dir}/additional_visualizations/
    #   - {cell_type}_methylation_patterns.pdf - Violin plots of methylation patterns
    #   - {cell_type}_binding_correlations.pdf - Correlation plots between binding signals and methylation
    #   - {cell_type}_distributions.pdf - KDE plots of methylation and binding signal distributions
    create_additional_visualizations(df, cell_type, PATHS['output_dir'])
    # Creates the following files/directories:
    # - {output_dir}/exo_enriched_analysis/{cell_type}/
    #   - exo_enriched_patterns.pdf - Plots showing methylation patterns and signal comparisons
    #   - exo_enriched_stats.txt - Statistical analysis of exo-enriched genes
    create_exo_enriched_comparisons(df, cell_type, PATHS['output_dir'])

#%% Analyze exo-enriched patterns specifically
logger.info("Analyzing exo-enriched patterns...")
# Creates files in:
# - {output_dir}/regulatory_patterns/{cell_type}_methylation_patterns.pdf
#   - Shows methylation patterns for promoter/gene body regions
#   - Includes boxplots and scatterplots colored by expression status
# - {output_dir}/tss_analysis/{cell_type}/{cell_type}_tss_analysis_stats.txt
#   - Contains statistical metrics for binding categories
# - {output_dir}/tss_analysis/{cell_type}/{binding_category}_{regulation_status}_genes.csv
#   - Lists of genes for each binding category and regulation status
analyze_exo_enriched_patterns(results, PATHS['output_dir'], exo_enrichment_threshold=1.5)  # For stricter filtering





#%% Code for debugging data merging issues (commented out)
# genes_df, expression_data, mecp2_binding = debug_merge_issues()

# Print sample gene IDs from each dataset
# print("\nSample gene IDs from each dataset:")
# print("\nGene annotations:")
# print(genes_df['gene_id'].head())
# print("\nExpression data (NEU):")
# print(expression_data['NEU']['gene_id'].head())
# print("\nExpression data (NSC):")
# print(expression_data['NSC']['gene_id'].head())

# Check for string formatting issues
# print("\nChecking for string formatting issues:")
# print("\nGene annotations gene_id example:", genes_df['gene_id'].iloc[0])
# print("Gene annotations gene_id type:", type(genes_df['gene_id'].iloc[0]))
# print("\nExpression data gene_id example:", expression_data['NEU']['gene_id'].iloc[0])
# print("Expression data gene_id type:", type(expression_data['NEU']['gene_id'].iloc[0]))

# Function to fix gene IDs
# def fix_gene_ids(df, column='gene_id'):
#     """Clean and standardize gene IDs"""
#     if column in df.columns:
#         # Remove version numbers if present
#         df[column] = df[column].str.split('.').str[0]
#         # Remove any whitespace
#         df[column] = df[column].str.strip()
#         # Ensure string type
#         df[column] = df[column].astype(str)
#     return df

# Fix gene IDs in all datasets
# genes_df = fix_gene_ids(genes_df)
# expression_data = {
#     cell_type: fix_gene_ids(df)
#     for cell_type, df in expression_data.items()
# }

# Try merge again
# for cell_type, expr_df in expression_data.items():
#     merged_df = genes_df.merge(expr_df, on='gene_id', how='inner')
#     print(f"\nMerged DataFrame for {cell_type} after fixing:")
#     print("Shape:", merged_df.shape)
#     if merged_df.shape[0] > 0:
#         print("Sample rows:")
#         print(merged_df[['gene_id', 'gene_name', 'log2FoldChange', 'padj']].head())

# %%
