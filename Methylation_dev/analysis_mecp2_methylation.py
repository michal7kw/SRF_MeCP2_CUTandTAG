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
from functions import *
from config import CONFIG, PATHS, logger, verify_paths
from cache_utils import clear_cache

def analyze_tss_binding_patterns(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze binding patterns specifically around TSS and their relationship with gene regulation"""
    
    # Check if we have methylation data
    has_methylation = all(col in df.columns for col in ['cpg_island_count', 'cpg_island_mean_methylation'])
    logger.info(f"Methylation data available: {has_methylation}")
    
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
        # Skip empty categories
        if len(category_df) == 0:
            logger.warning(f"No genes found in group: {category} for {cell_type}")
            continue
            
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
                stats = {
                    'count': len(group_df)
                }
                
                # Add methylation stats only if available
                if has_methylation:
                    stats.update({
                        'methylation': {
                            'mean': group_df['cpg_island_mean_methylation'].mean(),
                            'std': group_df['cpg_island_mean_methylation'].std(),
                            'median': group_df['cpg_island_mean_methylation'].median()
                        },
                        'cpg_count': {
                            'mean': group_df['cpg_island_count'].mean(),
                            'std': group_df['cpg_island_count'].std(),
                            'median': group_df['cpg_island_count'].median()
                        }
                    })
                
                group_stats[reg_status] = stats
                
                # Save gene lists
                output_file = os.path.join(tss_dir, f'{category}_{reg_status}_genes.csv')
                group_df.to_csv(output_file, index=False)
        
        results[category] = {
            'regulation_counts': pd.Series({k: len(v) for k, v in regulation_groups.items()}),
            'group_stats': group_stats
        }
    
    return results

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
    tss_results = {}
    
    for ct, expr_df in expression_data.items():
        logger.info(f"\nProcessing {ct}...")
        
        # Merge data
        merged_df = validate_and_merge_data(genes_df, expr_df, mecp2_binding)
        
        # Ensure column names match PyRanges requirements
        if 'chr' in merged_df.columns:
            merged_df = merged_df.rename(columns={
                'chr': 'Chromosome',
                'start': 'Start',
                'end': 'End'
            })
        
        logger.info(f"Merged data columns: {merged_df.columns.tolist()}")
        
        # Calculate methylation levels
        methylation_data = calculate_methylation_levels(
            merged_df,
            PATHS['medip_dir'],
            ct,
            PATHS['genome_fasta']
        )
        
        # Merge methylation data with other data if available
        if not methylation_data.empty and len(methylation_data.columns) > 0:
            results[ct] = pd.concat([merged_df, methylation_data], axis=1)
            logger.info("Methylation data added successfully")
        else:
            logger.warning("No methylation data available - proceeding with basic analysis")
            results[ct] = merged_df
        
        # Save intermediate results
        output_dir = os.path.join(PATHS['output_dir'], 'cpg_analysis', ct)
        os.makedirs(output_dir, exist_ok=True)
        
        logger.info(f"Saving results to {output_dir}")
        results[ct].to_csv(os.path.join(output_dir, f'{ct}_full_results.csv'), index=False)
        
        # Log data shape and columns for debugging
        logger.info(f"Results shape for {ct}: {results[ct].shape}")
        logger.info(f"Available columns: {results[ct].columns.tolist()}")
        
        # Perform TSS binding analysis
        logger.info(f"Performing TSS binding analysis for {ct}...")
        tss_results[ct] = analyze_tss_binding_patterns(results[ct], ct, PATHS['output_dir'])
        
        # Add detailed methylation analysis
        methylation_stats = analyze_methylation_patterns(results[ct], ct, output_dir)
        
        # Add new detailed binding-methylation analysis
        binding_results = analyze_binding_methylation_expression(results[ct], ct, output_dir)
        
        # Add new detailed analyses
        detailed_dir = os.path.join(output_dir, 'detailed_analysis')
        os.makedirs(detailed_dir, exist_ok=True)
        
        # Run new analyses
        binding_patterns = analyze_mecp2_detailed_binding(results[ct], ct, detailed_dir)
        occupancy_patterns = analyze_cpg_occupancy(results[ct], ct, detailed_dir)
        regulatory_features = analyze_regulatory_features(results[ct], ct, detailed_dir)
    
    return {
        'methylation_results': results,
        'tss_results': tss_results,
        'genes_df': genes_df,
        'expression_data': expression_data,
        'mecp2_binding': mecp2_binding
    }

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Analyze MeCP2 binding and methylation patterns')
    parser.add_argument('--experiment', type=str, default='align1_005', help='Experiment name')
    parser.add_argument('--processes', type=int, default=None, help='Number of processes to use')
    parser.add_argument('--cell-type', type=str, choices=['NEU', 'NSC', None], default=None, 
                       help='Cell type to analyze')
    parser.add_argument('--chromosome', type=str, default=None, help='Chromosome to analyze')
    parser.add_argument('--debug', action='store_true', help='Run in debug mode')
    parser.add_argument('--sample-size', type=int, default=100, help='Sample size for debug mode')
    parser.add_argument('--force-recompute', action='store_true', help='Force recomputation')
    
    args = parser.parse_args()
    
    # Update CONFIG with command line arguments
    CONFIG['experiment'] = args.experiment
    CONFIG['debug']['enabled'] = args.debug
    CONFIG['debug']['sample_size'] = args.sample_size
    
    # Run the analysis pipeline
    try:
        logger.info("Starting MeCP2 methylation analysis...")
        
        # Verify paths and input files
        verify_paths()
        verify_input_files()
        
        # Run the main analysis pipeline
        results = run_analysis_pipeline(
            force_recompute=args.force_recompute,
            n_processes=args.processes,
            cell_type=args.cell_type,
            chromosome=args.chromosome
        )
        
        if not results:
            logger.error("Analysis pipeline returned no results")
            return
        
        # Process CpG analysis results
        for ct, df in results['methylation_results'].items():
            if args.cell_type and ct != args.cell_type:
                continue
            
            logger.info(f"\nProcessing CpG analysis for {ct}...")
            
            # Debug print of available columns
            logger.info(f"Available columns for {ct}: {df.columns.tolist()}")
            
            # Verify required columns exist
            required_columns = ['cpg_island_count', 'cpg_island_mean_methylation']
            if not all(col in df.columns for col in required_columns):
                missing = [col for col in required_columns if col not in df.columns]
                logger.error(f"Missing required columns: {missing}")
                continue
            
            # Create output directories
            output_dir = os.path.join(PATHS['output_dir'], 'cpg_analysis', ct)
            os.makedirs(output_dir, exist_ok=True)
            
            # Generate CpG analysis visualizations and reports
            try:
                create_cpg_methylation_plots(df, ct, output_dir)
                generate_cpg_report(df, ct, output_dir)
                
                # Additional CpG-specific analyses
                if 'binding_type' in df.columns:
                    analyze_cpg_binding_correlation(df, ct, output_dir)
                if 'expression_status' in df.columns:
                    analyze_cpg_expression_correlation(df, ct, output_dir)
                
                logger.info(f"Completed CpG analysis for {ct}")
                
            except Exception as e:
                logger.error(f"Error in CpG analysis for {ct}: {str(e)}")
                logger.error(f"DataFrame info:\n{df.info()}")
        
        logger.info("\nAnalysis completed successfully")
        
    except Exception as e:
        logger.error(f"Error in main analysis: {str(e)}")
        raise

def verify_input_files():
    """Verify that all required input files exist"""
    required_files = [
        PATHS['cpg_islands_file'],
        PATHS['gtf_file'],
        PATHS['genome_fasta']
    ]
    
    for file_path in required_files:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Required file not found: {file_path}")
        logger.info(f"Found required file: {file_path}")

def analyze_cpg_binding_correlation(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze correlation between CpG methylation and MeCP2 binding"""
    
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df, x='binding_type', y='cpg_island_mean_methylation')
    plt.title(f'{cell_type}: CpG Methylation by MeCP2 Binding Type')
    plt.xlabel('Binding Type')
    plt.ylabel('CpG Island Methylation (%)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cpg_binding_correlation.pdf'))
    plt.close()
    
    # Save statistics
    stats_file = os.path.join(output_dir, 'cpg_binding_stats.txt')
    with open(stats_file, 'w') as f:
        f.write(f"CpG Methylation by Binding Type - {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        for binding_type in df['binding_type'].unique():
            subset = df[df['binding_type'] == binding_type]
            f.write(f"\n{binding_type}:\n")
            f.write(f"Count: {len(subset)}\n")
            f.write(f"Mean methylation: {subset['cpg_island_mean_methylation'].mean():.2f}%\n")
            f.write(f"Median methylation: {subset['cpg_island_mean_methylation'].median():.2f}%\n")
            f.write(f"Std Dev: {subset['cpg_island_mean_methylation'].std():.2f}%\n")

def analyze_cpg_expression_correlation(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze correlation between CpG methylation and gene expression"""
    
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df, x='cpg_island_mean_methylation', y='log2FoldChange', 
                    hue='binding_type', alpha=0.6)
    plt.title(f'{cell_type}: CpG Methylation vs Expression')
    plt.xlabel('CpG Island Methylation (%)')
    plt.ylabel('log2 Fold Change')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cpg_expression_correlation.pdf'))
    plt.close()

if __name__ == "__main__":
    main()