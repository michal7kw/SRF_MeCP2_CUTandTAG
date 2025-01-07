#%% Import required libraries
# Core data analysis and visualization libraries
import pandas as pd  # For data manipulation and analysis
import numpy as np   # For numerical operations
import seaborn as sns  # For statistical data visualization
import matplotlib.pyplot as plt  # For creating plots
from scipy import stats  # For statistical functions

# File and system operations
import os  # For operating system operations
import pyBigWig  # For handling BigWig files
import pysam  # For handling SAM/BAM files
import logging  # For logging messages
from typing import Dict, List, Tuple, Any  # For type hints
import pyranges as pr  # For genomic range operations

# Parallel processing and utilities
from functools import partial  # For partial function application
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor  # For parallel processing
from tqdm import tqdm  # For progress bars
import time  # For timing operations
import argparse  # For command line argument parsing

# Local imports
from functions import *  # Import all functions from functions.py
from config import *    # Import all configurations
from cache_utils import *  # Import caching utilities

def main():
    """
    Main function to run the MeCP2 methylation analysis pipeline.
    Handles command line arguments, sets up configurations, and orchestrates the analysis workflow.
    """
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Analyze MeCP2 binding and methylation patterns')
    parser.add_argument('--experiment', type=str, default='align1_005', 
                       help='Experiment name to analyze')
    parser.add_argument('--processes', type=int, default=None, 
                       help='Number of parallel processes to use')
    parser.add_argument('--cell-type', type=str, choices=['NEU', 'NSC', None], default=None, 
                       help='Cell type to analyze (NEU: neurons, NSC: neural stem cells)')
    parser.add_argument('--chromosome', type=str, default=None, 
                       help='Specific chromosome to analyze (e.g., "chr1")')
    parser.add_argument('--debug', action='store_true', 
                       help='Run in debug mode with additional logging')
    parser.add_argument('--sample-size', type=int, default=100, 
                       help='Sample size to use in debug mode')
    parser.add_argument('--force-recompute', action='store_true', 
                       help='Force recomputation instead of using cached results')
    
    args = parser.parse_args()
    
    # Update global configuration with command line arguments
    CONFIG['experiment'] = args.experiment
    CONFIG['debug']['enabled'] = args.debug
    CONFIG['debug']['sample_size'] = args.sample_size
    
    # Update file paths based on experiment name
    update_paths_for_experiment(args.experiment)
    
    # Main analysis workflow
    try:
        # Initial setup and verification
        logger.info("Starting MeCP2 methylation analysis...")
        logger.info(f"Using experiment: {args.experiment}")
        logger.info(f"MeCP2 directory: {PATHS['mecp2_dir']}")
        
        # Verify all required paths and input files exist
        verify_paths()
        
        # Execute the main analysis pipeline
        # This function:
        # 1. Loads and validates input data (MeCP2 binding data, methylation data, expression data)
        # 2. Maps MeCP2 binding regions to genes and standardizes chromosome names
        # 3. Calculates methylation levels in promoter and gene body regions
        # 4. Analyzes CpG islands and their relationship with MeCP2 binding
        # 5. Generates statistics about binding patterns and methylation levels
        # 6. Returns a dictionary containing:
        #    - methylation_results: Methylation data for each cell type
        #    - tss_results: Transcription start site analysis
        #    - genes_df: Gene information and annotations
        #    - expression_data: Gene expression changes
        #    - mecp2_binding: MeCP2 binding regions and signals
        results = run_analysis_pipeline(
            force_recompute=args.force_recompute,  # Skip cache and recompute results
            n_processes=args.processes,            # Number of CPU processes to use
            cell_type=args.cell_type,             # Specific cell type to analyze (NEU/NSC)
            chromosome=args.chromosome             # Specific chromosome to analyze
        )
        
        # Check if pipeline returned results
        if not results:
            logger.error("Analysis pipeline returned no results")
            return
        
        # Process results for each cell type
        for ct, df in results['methylation_results'].items():
            # Skip if specific cell type requested and this isn't it
            if args.cell_type and ct != args.cell_type:
                continue
            
            logger.info(f"\nProcessing CpG analysis for {ct}...")
            
            # Set up output directory for this cell type
            output_dir = os.path.join(PATHS['base_output_dir'], 'cpg_analysis', ct)
            os.makedirs(output_dir, exist_ok=True)
            
            # Generate visualizations and reports
            try:
                # Create standard methylation plots
                # Uses create_methylation_distribution_plot() to create distribution plots for methylation and binding signals
                create_cpg_methylation_plots(df, ct, output_dir)
                
                # Generate detailed methylation report
                # Uses save_analysis_summary() to save statistical summary of the analysis
                generate_cpg_report(df, ct, output_dir)
                
                # Additional analyses based on available data
                if 'binding_type' in df.columns:
                    # Analyze correlation between CpG methylation and binding
                    # Uses analyze_mecp2_binding_patterns() to analyze binding patterns and relationship with methylation
                    analyze_cpg_binding_correlation(df, ct, output_dir)
                if 'expression_status' in df.columns:
                    # Analyze correlation between CpG methylation and expression
                    # Uses plot_mecp2_regulatory_patterns() to visualize MeCP2's regulatory patterns
                    analyze_cpg_expression_correlation(df, ct, output_dir)
                
                logger.info(f"Completed CpG analysis for {ct}")
                
except Exception as e:
                # Log detailed error information for debugging
                logger.error(f"Error in CpG analysis for {ct}: {str(e)}")
                logger.error(f"DataFrame info:\n{df.info()}")
        
        logger.info("\nAnalysis completed successfully")
        
    except Exception as e:
        logger.error(f"Error in main analysis: {str(e)}")
        raise

if __name__ == "__main__":
    main()