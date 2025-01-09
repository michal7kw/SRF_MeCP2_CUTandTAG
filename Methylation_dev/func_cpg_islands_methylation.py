import pandas as pd
import pyranges as pr # type: ignore
from config import CONFIG, PATHS, logger
import numpy as np
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import os
from tqdm import tqdm
import pysam
from typing import Dict, List, Tuple, Any, Union
import pyBigWig


def load_cpg_islands(cpg_file: str) -> pr.PyRanges:
    """Load CpG islands from BED file into PyRanges object
    
    Format:
    chr1    3531624    3531843    611    CpG:    27
    """
    # Read the BED file
    cpg_df = pd.read_csv(cpg_file, sep='\t', header=None,
                        names=['Chromosome', 'Start', 'End', 'ID', 'Type', 'cpg_count'])
    
    # Convert CpG count to numeric (it's the last column)
    cpg_df['cpg_count'] = pd.to_numeric(cpg_df['cpg_count'])
    
    # Debug info
    logger.debug(f"Loaded CpG islands data shape: {cpg_df.shape}")
    logger.debug(f"Columns: {cpg_df.columns.tolist()}")
    logger.debug(f"First few rows:\n{cpg_df.head()}")
    
    # Convert to PyRanges
    cpg_islands = pr.PyRanges(cpg_df)
    logger.info(f"Loaded {len(cpg_islands)} CpG islands")
    return cpg_islands

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
            merged_df = merged_df.merge(binding_data, on='gene_id', how='left')
            merged_df['mecp2_bound'] = ~merged_df['binding_type'].isna()
        else:
            merged_df['mecp2_bound'] = False
            merged_df['binding_type'] = None
            merged_df['exo_signal'] = None
            merged_df['endo_signal'] = None
        
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

def identify_cpg_associated_regions(genes_df: pd.DataFrame, 
                                  cpg_islands: pr.PyRanges) -> pd.DataFrame:
    """Identify the specific CpG island regions within promoters and gene bodies"""
    
    # Create PyRanges objects for promoters and gene bodies
    promoters = pr.PyRanges(
        pd.DataFrame({
            'Chromosome': genes_df['chr'],
            'Start': genes_df['promoter_start'],
            'End': genes_df['promoter_end'],
            'gene_id': genes_df['gene_id'],
            'region_type': 'promoter'
        })
    )
    
    gene_bodies = pr.PyRanges(
        pd.DataFrame({
            'Chromosome': genes_df['chr'],
            'Start': genes_df['gene_body_start'],
            'End': genes_df['gene_body_end'],
            'gene_id': genes_df['gene_id'],
            'region_type': 'gene_body'
        })
    )
    
    # Find overlaps using intersection method
    promoter_cpg = promoters.join(cpg_islands)  # Changed from intersect to join
    gene_body_cpg = gene_bodies.join(cpg_islands)  # Changed from intersect to join
    
    cpg_regions = pd.DataFrame()
    
    if len(promoter_cpg) > 0:
        promoter_regions = promoter_cpg.df
        promoter_regions = pd.DataFrame({
            'gene_id': promoter_regions['gene_id'],
            'chr': promoter_regions['Chromosome'],
            'start': promoter_regions['Start'],
            'end': promoter_regions['End'],
            'region_type': 'promoter',
            'cpg_count': promoter_regions['cpg_count'].astype(int)  # Convert to int
        })
        cpg_regions = pd.concat([cpg_regions, promoter_regions])
    
    if len(gene_body_cpg) > 0:
        gene_body_regions = gene_body_cpg.df
        gene_body_regions = pd.DataFrame({
            'gene_id': gene_body_regions['gene_id'],
            'chr': gene_body_regions['Chromosome'],
            'start': gene_body_regions['Start'],
            'end': gene_body_regions['End'],
            'region_type': 'gene_body',
            'cpg_count': gene_body_regions['cpg_count'].astype(int)  # Convert to int
        })
        cpg_regions = pd.concat([cpg_regions, gene_body_regions])
    
    if cpg_regions.empty:
        logger.warning("No CpG islands found in promoter or gene body regions")
        return pd.DataFrame(columns=['gene_id', 'chr', 'start', 'end', 
                                   'region_type', 'cpg_count'])
    
    # Reset index and sort
    cpg_regions = cpg_regions.reset_index(drop=True)
    cpg_regions = cpg_regions.sort_values(['chr', 'start'])
    
    # Log summary
    logger.info(f"Found {len(cpg_regions)} CpG island regions "
                f"({len(cpg_regions[cpg_regions['region_type'] == 'promoter'])} in promoters, "
                f"{len(cpg_regions[cpg_regions['region_type'] == 'gene_body'])} in gene bodies)")
    logger.info(f"Average CpG count: {cpg_regions['cpg_count'].mean():.1f}")
    
    return cpg_regions

def calculate_contextual_methylation_cpg_only(merged_df: pd.DataFrame,
                                            cpg_regions: pd.DataFrame,
                                            medip_dir: str,
                                            cell_type: str,
                                            genome_fasta: str,
                                            n_processes: int = None) -> pd.DataFrame:
    """Calculate methylation levels only for CpG island regions with CpG density weighting
    
    Args:
        merged_df: DataFrame with gene information and binding data
        cpg_regions: DataFrame with CpG island regions and their CpG counts
        medip_dir: Directory containing MeDIP bigWig files
        cell_type: Cell type code ('NSC' or 'NEU')
        genome_fasta: Path to genome FASTA file
        n_processes: Number of parallel processes
    
    Returns:
        DataFrame with methylation levels weighted by CpG density
    """
    
    logger.info(f"Calculating CpG-weighted methylation for {cell_type}")
    logger.info(f"Number of input regions: {len(cpg_regions)}")
    
    # Calculate methylation for CpG regions
    methylation_results = calculate_methylation_levels_parallel(
        cpg_regions[['chr', 'start', 'end']],
        medip_dir,
        cell_type,
        genome_fasta,
        n_processes
    )
    
    # Combine results by region type
    results_by_gene = {}
    for gene_id in cpg_regions['gene_id'].unique():
        gene_regions = cpg_regions[cpg_regions['gene_id'] == gene_id]
        gene_results = {
            'gene_id': gene_id,
            'promoter_methylation': None,
            'promoter_cpg_count': 0,
            'promoter_cpg_density': 0,
            'gene_body_methylation': None,
            'gene_body_cpg_count': 0,
            'gene_body_cpg_density': 0
        }
        
        for region_type in ['promoter', 'gene_body']:
            type_regions = gene_regions[gene_regions['region_type'] == region_type]
            if not type_regions.empty:
                weighted_methylation = 0
                total_cpgs = 0
                total_length = 0
                
                for idx, region in type_regions.iterrows():
                    region_methylation = methylation_results.iloc[idx]['methylation']
                    region_length = region['end'] - region['start']
                    region_cpgs = region['cpg_count']
                    
                    if region_methylation is not None and region_cpgs > 0:
                        # Weight methylation by number of CpGs in the region
                        weighted_methylation += region_methylation * region_cpgs
                        total_cpgs += region_cpgs
                        total_length += region_length
                
                if total_cpgs > 0:
                    # Calculate CpG-density-weighted methylation
                    gene_results[f'{region_type}_methylation'] = weighted_methylation / total_cpgs
                    gene_results[f'{region_type}_cpg_count'] = total_cpgs
                    # Calculate CpG density per kb
                    gene_results[f'{region_type}_cpg_density'] = (total_cpgs * 1000) / total_length
                    
                    # Log extreme values for debugging
                    if gene_results[f'{region_type}_methylation'] > 90:
                        logger.debug(f"High methylation detected for {gene_id} {region_type}: "
                                   f"{gene_results[f'{region_type}_methylation']:.1f}%")
        
        results_by_gene[gene_id] = gene_results
    
    # Create final DataFrame
    result_df = pd.DataFrame(results_by_gene.values())
    
    # Add summary statistics
    logger.info("Methylation summary statistics:")
    for region_type in ['promoter', 'gene_body']:
        meth_col = f'{region_type}_methylation'
        density_col = f'{region_type}_cpg_density'
        valid_data = result_df[result_df[meth_col].notna()]
        
        if not valid_data.empty:
            logger.info(f"\n{region_type.title()} statistics:")
            logger.info(f"Mean methylation: {valid_data[meth_col].mean():.1f}%")
            logger.info(f"Median methylation: {valid_data[meth_col].median():.1f}%")
            logger.info(f"Mean CpG density: {valid_data[density_col].mean():.1f} CpGs/kb")
            logger.info(f"Number of regions: {len(valid_data)}")
    
    # Merge with original data
    final_df = merged_df[['gene_id', 'gene_name', 'binding_type', 'mecp2_bound', 
                         'log2FoldChange', 'padj']].merge(result_df, on='gene_id', how='right')
    
    # Add quality control flags
    final_df['low_cpg_promoter'] = final_df['promoter_cpg_density'] < 10  # Flag low CpG density
    final_df['high_meth_promoter'] = final_df['promoter_methylation'] > 90  # Flag suspiciously high methylation
    
    logger.info(f"Final output shape: {final_df.shape}")
    logger.info(f"Genes with low CpG density in promoter: {final_df['low_cpg_promoter'].sum()}")
    logger.info(f"Genes with high methylation in promoter: {final_df['high_meth_promoter'].sum()}")
    
    return final_df

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
        genes_df['gene_body_end'] = genes_df['gene_body_end'] + window_size
        
        # Debug coordinate ranges
        logger.info("\nCoordinate ranges after standardization:")
        for chrom in sorted(set(mecp2_binding['chr'].unique()) | set(genes_df['chr'].unique())):
            binding_regions = mecp2_binding[mecp2_binding['chr'] == chrom]
            genes = genes_df[genes_df['chr'] == chrom]
            if not binding_regions.empty and not genes.empty:
                logger.info(f"\nChromosome {chrom}:")
                logger.info(f"Binding regions: {len(binding_regions)}, range: {binding_regions['start'].min()}-{binding_regions['end'].max()}")
                logger.info(f"Genes: {len(genes)}, range: {genes['promoter_start'].min()}-{genes['gene_body_end'].max()}")
        
        # Create PyRanges objects
        binding_pr = pr.PyRanges(
            chromosomes=mecp2_binding['chr'],
            starts=mecp2_binding['start'].astype(int),
            ends=mecp2_binding['end'].astype(int)
        )
        
        genes_pr = pr.PyRanges(
            chromosomes=genes_df['chr'],
            starts=genes_df['promoter_start'].astype(int),
            ends=genes_df['gene_body_end'].astype(int)
        )
        
        # Add metadata
        for col in ['binding_type', 'exo_signal', 'endo_signal']:
            if col in mecp2_binding.columns:
                setattr(binding_pr, col, mecp2_binding[col].values)
        
        genes_pr.gene_id = genes_df['gene_id'].values
        
        # Find overlaps
        logger.info("\nFinding overlaps...")
        overlaps = binding_pr.join(genes_pr)
        
        if overlaps is None or len(overlaps) == 0:
            logger.warning("No overlaps found between binding regions and genes")
            return pd.DataFrame(columns=['gene_id', 'binding_type', 'exo_signal', 'endo_signal'])
        
        # Convert to DataFrame
        binding_genes = overlaps.as_df()
        logger.info(f"Found {len(binding_genes)} overlaps")
        
        # Group by gene_id and aggregate binding info
        gene_binding = binding_genes.groupby('gene_id').agg({
            'binding_type': lambda x: x.iloc[0] if len(set(x)) == 1 else 'both',
            'exo_signal': 'max',
            'endo_signal': 'max'
        }).reset_index()
        
        return gene_binding
        
    except Exception as e:
        logger.error(f"Error in map_binding_to_genes: {str(e)}")
        raise

def calculate_methylation_levels_parallel(region_df: pd.DataFrame,
                                        medip_dir: str,
                                        cell_type: str,
                                        genome_fasta: str,
                                        n_processes: int = None) -> pd.DataFrame:
    """Calculate methylation levels for genomic regions in parallel."""
    
    if n_processes is None:
        n_processes = max(1, multiprocessing.cpu_count() - 1)

    # Split data into chunks
    chunk_size = max(1, len(region_df) // (n_processes * 4))
    chunks = np.array_split(region_df, len(region_df) // chunk_size + 1)

    # Process chunks in parallel
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = [
            executor.submit(calculate_methylation_for_regions,
                          chunk, medip_dir, cell_type, genome_fasta)
            for chunk in chunks
        ]
        
        results = []
        for future in tqdm(futures, desc="Calculating methylation levels"):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                logger.error(f"Error processing chunk: {str(e)}")
                continue

    if not results:
        raise RuntimeError("No valid results obtained from parallel processing")
    
    return pd.concat(results, ignore_index=True)

def calculate_methylation_for_regions(region_df: pd.DataFrame,
                                    medip_dir: str,
                                    cell_type: str,
                                    genome_fasta: str) -> pd.DataFrame:
    """Calculate methylation for a set of genomic regions."""
    
    # Map cell types to file prefixes
    prefix_map = {'NEU': 'N', 'NSC': 'PP'}
    prefix = prefix_map.get(cell_type)
    if not prefix:
        raise ValueError(f"Unknown cell type: {cell_type}")

    # Get replicate file paths
    ip_files = [os.path.join(medip_dir, f"Medip_{prefix}_output_r{i}.bw") 
                for i in range(1, 4)]
    input_files = [os.path.join(medip_dir, f"Medip_{prefix}_input_r{i}.bw") 
                  for i in range(1, 4)]

    # Check files exist
    for f in ip_files + input_files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Missing required file: {f}")

    results = []
    with pysam.FastaFile(genome_fasta) as fasta:
        for _, region in region_df.iterrows():
            try:
                methylation, cpg_count = calculate_region_methylation(
                    region, ip_files, input_files, fasta
                )
                results.append({
                    'methylation': methylation,
                    'cpg_count': cpg_count
                })
            except Exception as e:
                logger.debug(f"Error processing region {region['chr']}:{region['start']}-{region['end']}: {str(e)}")
                results.append({
                    'methylation': 0,
                    'cpg_count': 0
                })

    return pd.DataFrame(results)

def calculate_region_methylation(region: pd.Series,
                               ip_files: List[str],
                               input_files: List[str],
                               fasta: pysam.FastaFile) -> Tuple[float, int]:
    """Calculate methylation with better handling of missing data"""
    
    replicate_values = []
    
    # Get sequence and validate CpG content
    sequence = fasta.fetch(region['chr'], int(region['start']), int(region['end']))
    cpg_count = sequence.upper().count('CG')
    
    if cpg_count == 0:
        logger.warning(f"No CpGs in region {region['chr']}:{region['start']}-{region['end']}")
        return None, cpg_count  # Mark as missing rather than zero

    valid_replicates = 0
    for ip_file, input_file in zip(ip_files, input_files):
        with pyBigWig.open(ip_file) as bw_ip, pyBigWig.open(input_file) as bw_input:
            try:
                ip_values = bw_ip.values(region['chr'], int(region['start']), int(region['end']))
                input_values = bw_input.values(region['chr'], int(region['start']), int(region['end']))
                
                if ip_values and input_values:
                    valid_pairs = [(ip, inp) for ip, inp in zip(ip_values, input_values)
                                 if ip is not None and inp is not None]
                    
                    if valid_pairs:
                        ip_vals, input_vals = zip(*valid_pairs)
                        methylation = calculate_normalized_methylation(ip_vals, input_vals, sequence)
                        if methylation is not None:
                            replicate_values.append(methylation)
                            valid_replicates += 1
                        
            except Exception as e:
                logger.debug(f"Error processing replicate: {str(e)}")
                continue

    # Require minimum number of valid replicates
    MIN_REPLICATES = 2
    if valid_replicates >= MIN_REPLICATES:
        return np.mean(replicate_values), cpg_count
    else:
        logger.warning(f"Insufficient valid replicates ({valid_replicates}) for region")
        return None, cpg_count

def calculate_normalized_methylation(ip_values: List[float],
                                  input_values: List[float],
                                  sequence: str) -> float:
    """Calculate methylation with better handling of edge cases"""
    
    # 1. Signal validation
    if not ip_values or not input_values:
        logger.warning("Missing IP or input values")
        return None  # Better to mark as missing than assume zero
        
    # 2. Calculate mean signals with minimum threshold
    MIN_SIGNAL = 0.001  # Prevent division by very small numbers
    ip_mean = max(MIN_SIGNAL, np.mean(ip_values))
    input_mean = max(MIN_SIGNAL, np.mean(input_values))
    
    # 3. Calculate enrichment
    enrichment = ip_mean / input_mean
    
    # 4. Normalize by CpG density with minimum threshold
    cpg_count = sequence.upper().count('CG')
    region_length = len(sequence)
    MIN_DENSITY = 0.001  # Minimum density threshold
    cpg_density = max(MIN_DENSITY, cpg_count / region_length)
    
    normalized_enrichment = enrichment / cpg_density
    
    # 5. Convert to methylation percentage with baseline
    # Assuming minimum ~20% methylation in mammalian cells
    BASE_METHYLATION = 20
    methylation = BASE_METHYLATION + min(80, normalized_enrichment * 15)
    
    return methylation
