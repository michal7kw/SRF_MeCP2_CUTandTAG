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