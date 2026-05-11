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
import func_cpg_islands_methylation
import func_tss_analysis
import func_regions_methylation

importlib.reload(functions)
importlib.reload(config) 
importlib.reload(cache_utils)
importlib.reload(func_cpg_islands_methylation)
importlib.reload(func_tss_analysis)
importlib.reload(func_regions_methylation)

from functions import *
from config import *
from cache_utils import *
from func_cpg_islands_methylation import *
from func_tss_analysis import *
from func_regions_methylation import *

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

    # Identify peaks associated with genes
    peak_regions = identify_peak_regions(merged_df, mecp2_binding)
    print("\nPeak regions shape:", peak_regions.shape)
    print("Peak regions columns:", peak_regions.columns.tolist())
    print("\nFirst few rows of peak_regions:")
    print(peak_regions.head())

    # Debug mecp2_binding data
    print("\nMeCP2 binding data shape:", mecp2_binding.shape)
    print("MeCP2 binding columns:", mecp2_binding.columns.tolist())
    print("\nFirst few rows of MeCP2 binding data:")
    print(mecp2_binding.head())

    # Debug merged_df
    print("\nMerged data shape:", merged_df.shape)
    print("Merged data columns:", merged_df.columns.tolist())
    print("\nFirst few rows of merged data:")
    print(merged_df.head())

    if not peak_regions.empty:
        # Calculate methylation in peaks
        peak_methylation = calculate_peak_methylation(
            peak_regions,
            PATHS['medip_dir'],
            ct,
            PATHS['genome_fasta'],
            n_processes=n_processes
        )
        print("\nPeak methylation shape:", peak_methylation.shape)
        print("Peak methylation columns:", peak_methylation.columns.tolist())
        print("\nFirst few rows of peak_methylation:")
        print(peak_methylation.head())
    else:
        print("\nNo peaks were identified - debugging peak identification:")
        print("1. Check if chromosome formats match:")
        print("Merged data chromosome format:", merged_df['chr'].head())
        print("MeCP2 binding chromosome format:", mecp2_binding['chr'].head())

    print("\nExpression data columns:", expression_data[ct].columns.tolist())
    print("\nFirst few rows of expression data:")
    print(expression_data[ct].head())

    # Analyze peak methylation patterns
    peak_results = analyze_peak_methylation_patterns(
        peak_methylation,
        expression_data[ct],
        ct,
        PATHS['output_dir']
    )

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
    from IPython import get_ipython
    get_ipython()
    in_notebook = True
except (NameError, ImportError):
    in_notebook = False

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
            self.experiment = 'align1_005'
            self.processes = None
            self.cell_type = 'NSC'
            self.chromosome = None
    args = Args()
    debug = False
    sample_size = 100
    force_recompute = True
    experiment = 'align1_005'
    n_processes = None
    cell_type = 'NSC'
    chromosome = None

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

#%% ########################## Load MeCP2 binding data ############################################################
mecp2_binding = pd.read_csv(os.path.join(PATHS['mecp2_dir'], PATHS['mecp2_file']))
print("MeCP2 binding data:")
print(mecp2_binding.head())
"""
MeCP2 binding data:
    chr    start      end  exo_signal  endo_signal  enrichment    pvalue  \
0  chr1  3531624  3531843    0.000000    15.177982    0.000000  1.000000   
1  chr1  3670619  3671074   29.287905     0.000000         inf  1.000000   
2  chr1  3671654  3672156   12.519179     6.646521    1.883569  0.193931   
3  chr1  4491701  4493673   18.136526     0.000000         inf  1.000000   
4  chr1  4571641  4572075   13.948418     0.000000         inf  1.000000   

  binding_type  peak_width_exo  peak_width_endo  significant  
0    endo_only            0.00            416.0        False  
1     exo_only          331.75              0.0        False  
2         both          378.00            336.0        False  
3     exo_only          243.00              0.0        False  
4     exo_only          317.50              0.0        False  
"""


#%% ########################## Load gene annotations ############################################################
genes_df, gene_name_to_id = load_gene_annotations(PATHS['gtf_file'])
print(dict(list(gene_name_to_id.items())[:5]))
print(genes_df.head())

# Filter by chromosome if specified
if chromosome:
    mecp2_binding = mecp2_binding[mecp2_binding['chr'] == chromosome]
    genes_df = genes_df[genes_df['chr'] == chromosome]
    logger.info(f"Filtered data for chromosome {chromosome}")

#%% ########################## Load RNA-seq expression data ############################################################
if cell_type:
    expression_data = {
        cell_type: load_expression_data(PATHS['rnaseq'][cell_type], gene_name_to_id, single_file=True)
    }
    logger.info(f"Analyzing cell type: {cell_type}")
else:
    expression_data = load_expression_data(PATHS['rnaseq'], gene_name_to_id)


#%% ########################## Print mapping statistics ############################################################
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

#%% ########################## 
for ct, expr_df in expression_data.items():
    # Merge datasets
    merged_df = validate_and_merge_data(genes_df, expr_df, mecp2_binding)

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
merged_df.head()

#%% ########################## Analyze peak methylation patterns ############################################################
# Identify peaks associated with genes
peak_regions = identify_peak_regions(merged_df, mecp2_binding)
print("\nPeak regions shape:", peak_regions.shape)
print("Peak regions columns:", peak_regions.columns.tolist())
print("\nFirst few rows of peak_regions:")
print(peak_regions.head())

# Debug mecp2_binding data
print("\nMeCP2 binding data shape:", mecp2_binding.shape)
print("MeCP2 binding columns:", mecp2_binding.columns.tolist())
print("\nFirst few rows of MeCP2 binding data:")
print(mecp2_binding.head())

# Debug merged_df
print("\nMerged data shape:", merged_df.shape)
print("Merged data columns:", merged_df.columns.tolist())
print("\nFirst few rows of merged data:")
print(merged_df.head())

if not peak_regions.empty:
    # Calculate methylation in peaks
    peak_methylation = calculate_peak_methylation(
        peak_regions,
        PATHS['medip_dir'],
        ct,
        PATHS['genome_fasta'],
        n_processes=n_processes
    )
    print("\nPeak methylation shape:", peak_methylation.shape)
    print("Peak methylation columns:", peak_methylation.columns.tolist())
    print("\nFirst few rows of peak_methylation:")
    print(peak_methylation.head())
else:
    print("\nNo peaks were identified - debugging peak identification:")
    print("1. Check if chromosome formats match:")
    print("Merged data chromosome format:", merged_df['chr'].head())
    print("MeCP2 binding chromosome format:", mecp2_binding['chr'].head())

print("\nExpression data columns:", expression_data[ct].columns.tolist())
print("\nFirst few rows of expression data:")
print(expression_data[ct].head())

# Analyze peak methylation patterns
peak_results = analyze_peak_methylation_patterns(
    peak_methylation,
    expression_data[ct],
    ct,
    PATHS['output_dir']
)
# #%%
# Print all columns without wrapping
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

#%% ########################## Load CpG islands ############################################################
cpg_islands = load_cpg_islands(PATHS['cpg_islands_file'])

results = {}

for ct, expr_df in expression_data.items():
    # Merge datasets
    merged_df = validate_and_merge_data(genes_df, expr_df, mecp2_binding)
    
    # Get CpG island regions within promoters and gene bodies
    cpg_regions = identify_cpg_associated_regions(merged_df, cpg_islands)
    
    # Filter for exo/both binding
    filtered_df = merged_df[merged_df['binding_type'].isin(['exo', 'both'])]
    
    # Calculate methylation levels only for CpG island regions
    results[ct] = calculate_contextual_methylation_cpg_only(
        filtered_df,
        cpg_regions,
        PATHS['medip_dir'],
        ct,
        PATHS['genome_fasta'],
        n_processes=n_processes
    )

merged_df.head()

#%% Save results to pickle file for later use
logger.info("Saving full results dictionary to pickle file...")
import pickle
with open(f"{PATHS['output_dir']}/full_results.pkl", 'wb') as f:
    pickle.dump(results, f)
logger.info("Results saved successfully")

#%% Load saved results for additional analysis
import pickle
logger.info("Loading full results dictionary from pickle file...")
with open(f"{PATHS['output_dir']}/full_results.pkl", 'rb') as f:
    results = pickle.load(f)
logger.info("Results loaded successfully")

#%%
print(results.keys())

#%%
print(results['NSC'].head())

#%% Debug the results structure
print("Results keys:", results.keys())
for cell_type, result_dict in results.items():
    print(f"\n{cell_type} data:")
    print("Type:", type(result_dict))
    print("Keys:", result_dict.keys())
    print("First item:", next(iter(result_dict.items())))

#%% Save detailed results to CSV files
for cell_type, result_df in results.items():
    # Save to CSV
    output_file = os.path.join(PATHS['output_dir'], f"{cell_type}_methylation_analysis.csv")
    result_df.to_csv(output_file, index=False)
    logger.info(f"Saved results for {cell_type} to {output_file}")

    # Print summary statistics
    print(f"\nSummary for {cell_type}:")
    print("\nColumns in the dataset:")
    for col in result_df.columns:
        print(f"- {col}")
    
    print("\nFirst few rows:")
    print(result_df.head())
    
    print("\nBasic statistics:")
    numeric_cols = result_df.select_dtypes(include=['float64', 'int64']).columns
    print(result_df[numeric_cols].describe())

#%%
tss_results = {}

for ct, expr_df in expression_data.items():
    # Analyze TSS binding patterns
    logger.info(f"Performing TSS binding analysis for {ct}...")
    tss_results[ct] = analyze_tss_binding_patterns(results[ct], ct, PATHS['output_dir'])


#%%
print(tss_results['NSC'].keys())
tss_results['NSC']['exo_enriched']