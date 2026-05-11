import pandas as pd
import numpy as np
import logging
from typing import Dict, List, Tuple, Any
import os
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from config import CONFIG, PATHS, logger
import pyranges as pr

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
        # genes_df['gene_body_end'] = genes_df['gene_body_end'] + window_size
        
        # Debug coordinate ranges
        logger.info("\nCoordinate ranges after standardization:")
        for chrom in sorted(set(mecp2_binding['chr'].unique()) | set(genes_df['chr'].unique())):
            binding_regions = mecp2_binding[mecp2_binding['chr'] == chrom]
            genes = genes_df[genes_df['chr'] == chrom]
            if not binding_regions.empty and not genes.empty:
                logger.info(f"\nChromosome {chrom}:")
                logger.info(f"Binding regions: {len(binding_regions)}, range: {binding_regions['peak_start'].min()}-{binding_regions['peak_end'].max()}")
                logger.info(f"Genes: {len(genes)}, range: {genes['promoter_start'].min()}-{genes['gene_body_end'].max()}")
        
        # Create PyRanges objects
        binding_pr = pr.PyRanges(
            chromosomes=mecp2_binding['chr'],
            starts=mecp2_binding['peak_start'].astype(int),
            ends=mecp2_binding['peak_end'].astype(int)
        )
        
        genes_pr = pr.PyRanges(
            chromosomes=genes_df['chr'],
            starts=genes_df['promoter_start'].astype(int),
            ends=genes_df['gene_body_end'].astype(int)
        )
        
        # Add metadata
        for col in ['binding_type', 'exo_signal', 'endo_signal', 'peak_width_exo', 'peak_width_endo']:
            if col in mecp2_binding.columns:
                setattr(binding_pr, col, mecp2_binding[col].values)
        
        genes_pr.gene_id = genes_df['gene_id'].values
        
        # Find overlaps
        logger.info("\nFinding overlaps...")
        overlaps = binding_pr.join(genes_pr)
        
        if overlaps is None or len(overlaps) == 0:
            logger.warning("No overlaps found between binding regions and genes")
            return pd.DataFrame(columns=['gene_id', 'binding_type', 'exo_signal', 'endo_signal',
                                      'peak_start', 'peak_end', 'peak_width_exo', 'peak_width_endo'])
        
        # Convert to DataFrame
        binding_genes = overlaps.as_df()
        logger.info(f"Found {len(binding_genes)} overlaps")
        
        # Update groupby to include peak coordinates
        gene_binding = binding_genes.groupby('gene_id').agg({
            'binding_type': lambda x: x.iloc[0] if len(set(x)) == 1 else 'both',
            'exo_signal': 'max',
            'endo_signal': 'max',
            'Start': 'min',  # Keep the start of the first peak
            'End': 'max',    # Keep the end of the last peak
            'peak_width_exo': 'max',
            'peak_width_endo': 'max'
        }).reset_index()
        
        # Rename columns for clarity
        gene_binding = gene_binding.rename(columns={
            'Start': 'peak_start',
            'End': 'peak_end'
        })
        
        return gene_binding
        
    except Exception as e:
        logger.error(f"Error in map_binding_to_genes: {str(e)}")
        raise

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
            merged_df = merged_df.merge(
                binding_data[['gene_id', 'binding_type', 'exo_signal', 'endo_signal',
                            'peak_start', 'peak_end', 'peak_width_exo', 'peak_width_endo']],
                on='gene_id', 
                how='left'
            )
            merged_df['mecp2_bound'] = ~merged_df['binding_type'].isna()
        else:
            merged_df['mecp2_bound'] = False
            merged_df['binding_type'] = None
            merged_df['exo_signal'] = None
            merged_df['endo_signal'] = None
            merged_df['peak_start'] = None
            merged_df['peak_end'] = None
            merged_df['peak_width_exo'] = None
            merged_df['peak_width_endo'] = None
        
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
