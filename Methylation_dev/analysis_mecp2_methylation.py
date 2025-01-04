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
from typing import Dict, List, Tuple
import pyranges as pr
from functools import partial
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import time
import argparse
from functions import *
from config import CONFIG, PATHS, logger
from cache_utils import clear_cache

#%% Configure debug mode
# Set up command line arguments to enable debug mode and specify sample size
parser = argparse.ArgumentParser(description='Analyze MeCP2 methylation patterns')
parser.add_argument('--debug', action='store_true', help='Run in debug mode with minimal dataset')
parser.add_argument('--sample-size', type=int, default=100, help='Number of genes to use in debug mode')
parser.add_argument('--force-recompute', action='store_true', 
                   help='Force recomputation of all analysis stages')
parser.add_argument('--experiment', type=str, default='align2_005_consistent_peaks',
                   help='Experiment name to analyze')
parser.add_argument('--processes', type=int, default=None,
                   help='Number of processes to use for parallel computation')
parser.add_argument('--cell-type', type=str, choices=['NEU', 'NSC'], 
                   help='Cell type to analyze')
parser.add_argument('--chromosome', type=str, 
                   help='Chromosome to analyze (e.g., chr1)')
args = parser.parse_args()

# Modify configuration based on debug mode
if args.debug:
    logger.info(f"Running in DEBUG mode with sample size: {args.sample_size}")
    CONFIG['debug'] = {
        'enabled': True,
        'sample_size': args.sample_size
    }

#%% Initialize output directory
# Create directory for storing analysis results
os.makedirs(PATHS['output_dir'], exist_ok=True)

#%% Load and process gene annotations and expression data
# Load gene annotations and create gene name mapping
genes_df, gene_name_to_id = load_gene_annotations(PATHS['gtf_file'])

# Load RNA-seq expression data with gene ID mapping
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

#%% Load MeCP2 binding data
# Read MeCP2 ChIP-seq binding data from CSV file
mecp2_binding = pd.read_csv(os.path.join(PATHS['mecp2_dir'], PATHS['mecp2_file']))

#%% Run analysis pipeline
def run_analysis_pipeline(force_recompute: bool = False, n_processes: int = None,
                        cell_type: str = None, chromosome: str = None) -> Dict[str, Any]:
    """Run the analysis pipeline with support for cell type and chromosome filtering"""
    if force_recompute:
        clear_cache()
        logger.info("Cache cleared due to force recompute flag")
    
    # Load all data first
    mecp2_binding = pd.read_csv(os.path.join(PATHS['mecp2_dir'], PATHS['mecp2_file']))
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
    for ct, expr_df in expression_data.items():
        merged_df = validate_and_merge_data(genes_df, expr_df, mecp2_binding)
        results[ct] = calculate_contextual_methylation(
            merged_df,
            PATHS['medip_dir'],
            ct,
            PATHS['genome_fasta'],
            n_processes=n_processes
        )
    
    return {
        'methylation_results': results,
        'genes_df': genes_df,
        'expression_data': expression_data,
        'mecp2_binding': mecp2_binding
    }

# Run pipeline with specified filters
analysis_results = run_analysis_pipeline(
    force_recompute=args.force_recompute,
    n_processes=args.processes,
    cell_type=args.cell_type,
    chromosome=args.chromosome
)
results = analysis_results['methylation_results']

#%% Save results to pickle file
# Store full results dictionary for later analysis
logger.info("Saving full results dictionary to pickle file...")
import pickle
with open(f"{PATHS['output_dir']}/full_results.pkl", 'wb') as f:
    pickle.dump(results, f)
logger.info("Results saved successfully")

#%% Load results from pickle (commented out)
# Code for loading previously saved results
# logger.info("Loading full results dictionary from pickle file...")
# with open(f"{PATHS['output_dir']}/full_results.pkl", 'rb') as f:
#     results = pickle.load(f)
# logger.info("Results loaded successfully")

#%% Generate visualizations
# Create plots showing methylation patterns
create_methylation_plots(results, PATHS['output_dir'])

#%% Perform statistical analysis
# Run statistical tests on methylation data
stats_results = perform_statistical_analysis(results)

#%% Save detailed results
# Export results to CSV files for each cell type
for cell_type, df in results.items():
    df.to_csv(f"{PATHS['output_dir']}/{cell_type}_methylation_analysis.csv", 
                index=False)

#%% Save statistical analysis
# Write statistical test results to text file
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

#%% Load results for additional analysis
# Load previously saved results from pickle file
import pickle
logger.info("Loading full results dictionary from pickle file...")
with open(f"{PATHS['output_dir']}/full_results.pkl", 'rb') as f:
    results = pickle.load(f)
logger.info("Results loaded successfully")

#%% Analyze MeCP2-regulated genes
# Perform analysis of genes regulated by MeCP2
create_mecp2_regulated_analysis(results, PATHS['output_dir'])

#%% Print regulation summary
# Display summary of regulated genes
print_regulated_summary(results)

#%% Debug gene filtering
# Run debugging for regulated genes filtering
debug_regulated_genes_filtering(results)

#%% Print modified summary
# Display improved summary of regulated genes
print_regulated_summary_v2(results)

#%% Generate focused visualization
# Create plots highlighting significant differences
plot_significant_differences(results, PATHS['output_dir'])

#%% Generate comprehensive visualization
# Create detailed plots of MeCP2 regulatory patterns
plot_mecp2_regulatory_patterns(results, PATHS['output_dir'])

#%% Print final summary
# Display final summary of regulated genes
print_regulated_summary_v2(results)

#%% Generate final visualization
# Create final plots of significant differences
plot_significant_differences(results, PATHS['output_dir'])


#%% Debug merge issues (commented out)
# Code for debugging data merging issues and gene ID formatting
# genes_df, expression_data, mecp2_binding = debug_merge_issues()

################ Print sample gene IDs from each dataset
# print("\nSample gene IDs from each dataset:")
# print("\nGene annotations:")
# print(genes_df['gene_id'].head())
# print("\nExpression data (NEU):")
# print(expression_data['NEU']['gene_id'].head())
# print("\nExpression data (NSC):")
# print(expression_data['NSC']['gene_id'].head())

################# Check for any string formatting issues
# print("\nChecking for string formatting issues:")
# print("\nGene annotations gene_id example:", genes_df['gene_id'].iloc[0])
# print("Gene annotations gene_id type:", type(genes_df['gene_id'].iloc[0]))
# print("\nExpression data gene_id example:", expression_data['NEU']['gene_id'].iloc[0])
# print("Expression data gene_id type:", type(expression_data['NEU']['gene_id'].iloc[0]))

################# Try to fix the merge issue
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

################# Fix gene IDs in all datasets
# genes_df = fix_gene_ids(genes_df)
# expression_data = {
#     cell_type: fix_gene_ids(df)
#     for cell_type, df in expression_data.items()
# }

################# Try merge again
# for cell_type, expr_df in expression_data.items():
#     merged_df = genes_df.merge(expr_df, on='gene_id', how='inner')
#     print(f"\nMerged DataFrame for {cell_type} after fixing:")
#     print("Shape:", merged_df.shape)
#     if merged_df.shape[0] > 0:
#         print("Sample rows:")
#         print(merged_df[['gene_id', 'gene_name', 'log2FoldChange', 'padj']].head())

# After loading results
logger.info("Performing detailed methylation analysis...")
try:
    analyze_methylation_patterns_detailed(results, PATHS['output_dir'])
    logger.info("Completed detailed methylation analysis")
except Exception as e:
    logger.error(f"Error in detailed methylation analysis: {str(e)}")

#%% Analyze binding enrichment groups
logger.info("Analyzing binding enrichment groups...")
analyze_binding_enrichment_groups(results, PATHS['output_dir'])

#%% Special analysis for exo-enriched upregulated genes
def analyze_exo_upregulated(results: Dict[str, pd.DataFrame], output_dir: str):
    """Detailed analysis of upregulated genes with exo MeCP2 binding"""
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

logger.info("Performing special analysis of exo-enriched upregulated genes...")
analyze_exo_upregulated(results, PATHS['output_dir'])

#%% Perform detailed binding enrichment analysis
logger.info("Performing detailed binding enrichment analysis...")
analyze_binding_enrichment_detailed(results, PATHS['output_dir'])

#%% Create additional visualizations
logger.info("Creating additional visualizations...")
for cell_type, df in results.items():
    create_additional_visualizations(df, cell_type, PATHS['output_dir'])
    create_exo_enriched_comparisons(df, cell_type, PATHS['output_dir'])

# After loading results
for cell_type, cell_data in results.items():
    binding_analysis = analyze_binding_enrichment_groups(cell_data, PATHS['output_dir'])
    
    if binding_analysis:  # Check if we got results
        # Print summary of key findings
        for group_name, group_results in binding_analysis.items():
            print(f"\n{cell_type} - {group_name}:")
            print(f"Total genes: {group_results['regulation_counts'].sum()}")
            print("Regulation distribution:")
            print(group_results['regulation_counts'])