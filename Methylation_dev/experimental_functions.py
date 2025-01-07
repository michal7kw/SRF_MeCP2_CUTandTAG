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

#Print summary for regulated genes
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

#Create detailed visualizations for MeCP2-bound regulated genes
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

#Visualization functions
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
