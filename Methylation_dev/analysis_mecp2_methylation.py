#!/usr/bin/env python3
"""
MeCP2 Methylation Analysis Pipeline

This script performs comprehensive analysis of MeCP2 binding patterns,
DNA methylation, and gene expression relationships. It specifically focuses on:
1. CpG island identification in promoter regions
2. Methylation patterns in different binding contexts (Exo/Endo)
3. Correlation with gene expression
4. Detailed analysis of exo-enriched genes

Author: [Your Name]
Date: [Current Date]
"""

#%% Import required libraries
import argparse
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
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from functions import (
    load_gene_annotations, load_expression_data, identify_promoter_cpg_islands,
    validate_and_merge_data, calculate_contextual_methylation,
    analyze_binding_enrichment_groups, analyze_cpg_methylation_patterns,
    plot_methylation_comparison, analyze_cell_types, analyze_exo_enriched_genes,
    analyze_methylation_patterns_parallel, create_group_visualizations,
    create_analysis_visualizations, save_analysis_results, print_summary_statistics
)
from config import CONFIG, PATHS, logger
import multiprocessing

#%% Configure argument parser
def setup_argument_parser():
    """Configure command line argument parser with detailed help messages"""
    parser = argparse.ArgumentParser(
        description='Analyze MeCP2 methylation patterns',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--debug', action='store_true',
                       help='Run in debug mode with reduced dataset')
    parser.add_argument('--sample-size', type=int, default=100,
                       help='Number of genes to analyze in debug mode')
parser.add_argument('--force-recompute', action='store_true', 
                   help='Force recomputation of all analysis stages')
    parser.add_argument('--experiment', type=str, default='align2_005',
                       help='Experiment identifier')
parser.add_argument('--processes', type=int, default=None,
                       help='Number of parallel processes to use')
parser.add_argument('--cell-type', type=str, choices=['NEU', 'NSC'], 
                       help='Specific cell type to analyze')
parser.add_argument('--chromosome', type=str, 
                       help='Specific chromosome to analyze (e.g., chr1)')
    return parser.parse_args()

#%% Main analysis pipeline
def main():
    """Main analysis pipeline"""
    # Parse arguments and configure
    args = setup_argument_parser()
    
    # Configure debug mode if requested
if args.debug:
    logger.info(f"Running in DEBUG mode with sample size: {args.sample_size}")
        CONFIG['debug'] = {'enabled': True, 'sample_size': args.sample_size}

    # Initialize output directories
os.makedirs(PATHS['output_dir'], exist_ok=True)

    #%% Step 1: Load and prepare input data
    logger.info("Loading gene annotations...")
    genes_df, gene_name_to_id = load_gene_annotations(PATHS['gtf_file'])
    
    # Apply chromosome filter if specified
    if args.chromosome:
        genes_df = genes_df[genes_df['chr'] == args.chromosome]
        logger.info(f"Filtered data for chromosome {args.chromosome}")

    #%% Step 2: Load expression data
    logger.info("Loading differential expression data...")
    expression_data = load_expression_data_for_analysis(args, gene_name_to_id)

    #%% Step 3: Identify CpG islands
    logger.info("Identifying CpG islands in promoter regions...")
    cpg_associations = identify_promoter_cpg_islands(
        genes_df, 
        PATHS['cpg_islands_file']
    )

    #%% Step 4: Analyze methylation patterns
    results = analyze_cell_types(
        genes_df,
        expression_data,
        cpg_associations,
        args
    )

    #%% Step 5: Special analysis of exo-enriched genes
    analyze_exo_enriched_genes(results)

    #%% Step 6: Run parallel methylation analysis
    logger.info("Running parallel methylation analysis...")
    analyze_methylation_patterns_parallel(
        results, 
        PATHS['output_dir'],
        n_processes=args.processes
    )

    #%% Step 7: Create visualizations
    create_all_visualizations(results)

    #%% Step 8: Save detailed results
    save_all_results(results)

    logger.info("Analysis completed successfully")

def load_expression_data_for_analysis(args, gene_name_to_id):
    """Load expression data based on specified parameters"""
    if args.cell_type:
        return {
            args.cell_type: load_expression_data(
                PATHS['rnaseq'][args.cell_type],
                gene_name_to_id,
                single_file=True
            )
        }
    return load_expression_data(PATHS['rnaseq'], gene_name_to_id)

def analyze_cell_types(genes_df, expression_data, cpg_associations, args):
    """Analyze methylation patterns for each cell type"""
    results = {}
    
    if not expression_data:
        logger.error("No expression data provided")
        return results
        
    for cell_type, expr_df in expression_data.items():
        logger.info(f"\nAnalyzing {cell_type}...")
        
        try:
            # Merge and validate data
            merged_df = validate_and_merge_data(genes_df, expr_df, cpg_associations)
            if merged_df is None or len(merged_df) == 0:
                logger.warning(f"No valid data for {cell_type}")
                continue
                
            # Calculate methylation levels
            methylation_results = calculate_contextual_methylation(
            merged_df,
            PATHS['medip_dir'],
                cell_type,
            PATHS['genome_fasta'],
                n_processes=args.processes
            )
            
            if methylation_results is None:
                logger.warning(f"No methylation results for {cell_type}")
                continue
                
            # Store results even if some analyses fail
            results[cell_type] = {
                'methylation': methylation_results,
                'binding': None,
                'cpg': None
            }
            
            # Try to run additional analyses
            try:
                results[cell_type]['binding'] = analyze_binding_enrichment_groups(
                    methylation_results, 
                    PATHS['output_dir']
                )
            except Exception as e:
                logger.error(f"Error in binding analysis for {cell_type}: {str(e)}")
                
            try:
                results[cell_type]['cpg'] = analyze_cpg_methylation_patterns(
                    methylation_results, 
                    cpg_associations
                )
            except Exception as e:
                logger.error(f"Error in CpG analysis for {cell_type}: {str(e)}")
                
        except Exception as e:
            logger.error(f"Error analyzing {cell_type}: {str(e)}")
            continue
    
    return results

def store_and_visualize_results(methylation_results, binding_analysis, 
                              cpg_analysis, cell_type):
    """Store results and create initial visualizations"""
    results = {
        'methylation': methylation_results,
        'binding': binding_analysis,
        'cpg': cpg_analysis
    }
    
    # Create initial visualizations
    plot_methylation_comparison(binding_analysis, PATHS['output_dir'])
    
    # Print summary statistics
    print_summary_statistics(cell_type, binding_analysis, cpg_analysis)
    
    return results

def create_all_visualizations(results):
    """Create comprehensive visualizations for all results"""
    if not results:
        logger.warning("No results available for visualization")
        return
        
    logger.info("Creating comprehensive visualizations...")
    for cell_type, cell_results in results.items():
        try:
            # Check if we have binding results
            binding_results = cell_results.get('binding')
            if binding_results and isinstance(binding_results, dict):
                for category, data in binding_results.items():
                    if data and isinstance(data, dict) and 'gene_list' in data:
                        create_group_visualizations(
                            data['gene_list'],
                            category,
                            cell_type,
                            PATHS['output_dir']
                        )
                    else:
                        logger.warning(f"No valid data for {category} in {cell_type}")
            
            # Check if we have methylation results
            methylation_results = cell_results.get('methylation')
            if methylation_results is not None and isinstance(methylation_results, pd.DataFrame):
                create_analysis_visualizations(
                    methylation_results,
                    cell_type,
                    PATHS['output_dir']
                )
            else:
                logger.warning(f"No valid methylation results for {cell_type}")
except Exception as e:
            logger.error(f"Error creating visualizations for {cell_type}: {str(e)}")
            continue

def save_all_results(results):
    """Save all analysis results"""
    logger.info("Saving detailed analysis results...")
    for cell_type, cell_results in results.items():
        save_analysis_results(
            cell_results,
            cell_type,
            PATHS['output_dir']
        )

if __name__ == "__main__":
    main()