##########################################################################################################################

def calculate_contextual_methylation(merged_df: pd.DataFrame,
                                   medip_dir: str,
                                   cell_type: str,
                                   genome_fasta: str,
                                   n_processes: int = None) -> pd.DataFrame:
    """Calculate methylation levels for promoter and gene body regions.
    
    Args:
        merged_df: DataFrame with columns from validate_and_merge_data():
            - chr, gene_start, gene_end, strand (from GTF)
            - promoter_start, promoter_end (calculated)
            - gene_body_start, gene_body_end (calculated)
            - gene_id, gene_name (identifiers)
            - binding_type (from mecp2_binding)
            - mecp2_bound (boolean)
            - exo_signal, endo_signal (from binding data)
            - log2FoldChange, padj (from expression data)
        medip_dir: Directory containing MeDIP bigWig files
        cell_type: Cell type code ('NSC' or 'NEU')
        genome_fasta: Path to genome FASTA file
        n_processes: Number of parallel processes
    
    Returns:
        DataFrame with all columns from merged_df plus:
        - promoter_methylation: float (0-100)
        - promoter_cpg_count: int
        - gene_body_methylation: float (0-100)
        - gene_body_cpg_count: int
    """
    logger.info(f"Starting methylation calculation for {cell_type}")
    logger.info(f"Input data shape: {merged_df.shape}")
    
    # Verify required columns
    required_columns = ['chr', 'promoter_start', 'promoter_end', 
                       'gene_body_start', 'gene_body_end']
    missing_columns = [col for col in required_columns if col not in merged_df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    # Prepare region data
    promoter_data = merged_df[['chr', 'promoter_start', 'promoter_end']].rename(
        columns={'promoter_start': 'start', 'promoter_end': 'end'}
    )
    gene_body_data = merged_df[['chr', 'gene_body_start', 'gene_body_end']].rename(
        columns={'gene_body_start': 'start', 'gene_body_end': 'end'}
    )

    logger.info("Calculating methylation for promoters and gene bodies")
    # Calculate methylation for both regions in parallel
    with ThreadPoolExecutor(max_workers=2) as executor:
        promoter_future = executor.submit(
            calculate_methylation_levels_parallel,
            promoter_data, medip_dir, cell_type, genome_fasta, n_processes
        )
        gene_body_future = executor.submit(
            calculate_methylation_levels_parallel,
            gene_body_data, medip_dir, cell_type, genome_fasta, n_processes
        )
        
        try:
            promoter_methylation = promoter_future.result()
            gene_body_methylation = gene_body_future.result()
        except Exception as e:
            logger.error(f"Error in parallel methylation calculation: {str(e)}")
            raise

    # Create result DataFrame preserving all original columns
    result_df = merged_df.copy()
    
    # Add methylation data
    result_df['promoter_methylation'] = promoter_methylation['methylation']
    result_df['promoter_cpg_count'] = promoter_methylation['cpg_count']
    result_df['gene_body_methylation'] = gene_body_methylation['methylation']
    result_df['gene_body_cpg_count'] = gene_body_methylation['cpg_count']

    # Fill NaN values if any
    result_df['promoter_methylation'] = result_df['promoter_methylation'].fillna(0)
    result_df['gene_body_methylation'] = result_df['gene_body_methylation'].fillna(0)
    result_df['promoter_cpg_count'] = result_df['promoter_cpg_count'].fillna(0)
    result_df['gene_body_cpg_count'] = result_df['gene_body_cpg_count'].fillna(0)

    # Ensure binding columns exist
    if 'binding_type' not in result_df.columns:
        result_df['binding_type'] = 'none'
    if 'mecp2_bound' not in result_df.columns:
        result_df['mecp2_bound'] = False

    logger.info(f"Completed methylation calculation. Output shape: {result_df.shape}")
    logger.info("Sample of results:")
    logger.info(result_df[['gene_id', 'promoter_methylation', 
                          'gene_body_methylation', 'binding_type']].head())

    return result_df

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

#########################################################################################

def load_cpg_islands1(cpg_file: str) -> pr.PyRanges:
    """Load CpG islands from BED file into PyRanges object"""
    cpg_df = pd.read_csv(cpg_file, sep='\t', header=None,
                        names=['Chromosome', 'Start', 'End', 'ID', 'Type', 'CpG_count'])
    
    # Convert to PyRanges for efficient genomic operations
    cpg_islands = pr.PyRanges(cpg_df)
    logger.info(f"Loaded {len(cpg_islands)} CpG islands")
    return cpg_islands

def identify_cpg_associated_genes(genes_df: pd.DataFrame, 
                                cpg_islands: pr.PyRanges,
                                promoter_window: int = 2000) -> pd.DataFrame:
    """Identify genes with CpG islands in promoter or gene body"""
    
    # Create PyRanges objects for promoters and gene bodies
    promoters = pr.PyRanges(
        pd.DataFrame({
            'Chromosome': genes_df['chr'],
            'Start': genes_df['promoter_start'],
            'End': genes_df['promoter_end'],
            'gene_id': genes_df['gene_id']
        })
    )
    
    gene_bodies = pr.PyRanges(
        pd.DataFrame({
            'Chromosome': genes_df['chr'],
            'Start': genes_df['gene_body_start'],
            'End': genes_df['gene_body_end'],
            'gene_id': genes_df['gene_id']
        })
    )
    
    # Find overlaps
    promoter_cpg = promoters.join(cpg_islands).df
    gene_body_cpg = gene_bodies.join(cpg_islands).df
    
    # Get unique gene IDs with CpG islands
    cpg_genes = pd.DataFrame({
        'gene_id': pd.concat([promoter_cpg['gene_id'], 
                            gene_body_cpg['gene_id']]).unique()
    })
    
    # Add CpG association information
    cpg_genes['has_promoter_cpg'] = cpg_genes['gene_id'].isin(promoter_cpg['gene_id'])
    cpg_genes['has_gene_body_cpg'] = cpg_genes['gene_id'].isin(gene_body_cpg['gene_id'])
    
    return cpg_genes

def filter_for_analysis(merged_df: pd.DataFrame, 
                       cpg_genes: pd.DataFrame) -> pd.DataFrame:
    """Filter genes for analysis based on CpG islands and MeCP2 binding"""
    
    # Merge with CpG information
    filtered_df = merged_df.merge(cpg_genes, on='gene_id', how='inner')
    
    # Filter for exo binding
    filtered_df = filtered_df[
        (filtered_df['binding_type'].isin(['exo', 'both'])) &
        (filtered_df['has_promoter_cpg'] | filtered_df['has_gene_body_cpg'])
    ]
    
    logger.info(f"Found {len(filtered_df)} genes with CpG islands and exo/both binding")
    return filtered_df

#########################################################################################

def add_methylation_qc_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Enhanced quality control metrics"""
    df = df.copy()
    
    try:
        # Calculate region lengths and CpG density
        if 'start' in df.columns and 'end' in df.columns:
            region_lengths = df['end'] - df['start']
        else:
            logger.warning("No start/end columns found. Using default region length.")
            region_lengths = 1000  # default length
            
        # Calculate improved metrics
        df['region_length'] = region_lengths
        df['cpg_density'] = df['cpg_count'] / df['region_length']
        
        # Add coverage uniformity if available
        try:
            df['coverage_uniformity'] = calculate_coverage_uniformity(df)
        except Exception as e:
            logger.warning(f"Could not calculate coverage uniformity: {str(e)}")
            df['coverage_uniformity'] = 1.0  # default value
            
        # Add signal-to-noise ratio if available
        try:
            df['signal_to_noise'] = calculate_signal_to_noise(df)
        except Exception as e:
            logger.warning(f"Could not calculate signal-to-noise ratio: {str(e)}")
            df['signal_to_noise'] = 2.0  # default value
        
        # Calculate composite quality score
        df['methylation_confidence'] = df.apply(
            lambda x: calculate_enhanced_confidence_score(
                x['methylation'],
                x['cpg_density'],
                x['cpg_count'],
                x.get('coverage_uniformity', 1.0),
                x.get('signal_to_noise', 2.0)
            ),
            axis=1
        )
        
        return df
        
    except Exception as e:
        logger.error(f"Error in add_methylation_qc_metrics: {str(e)}")
        logger.error(f"DataFrame columns: {df.columns.tolist()}")
        logger.error(f"DataFrame head:\n{df.head().to_string()}")
        raise

def calculate_methylation_levels(region_df: pd.DataFrame,
                               medip_dir: str,
                               cell_type_prefix: str,
                               genome_fasta: str,
                               use_cache: bool = True) -> pd.DataFrame:
    """Calculate methylation levels for genomic regions"""
    # Ensure we have the required columns
    required_cols = ['chr', 'start', 'end']
    if not all(col in region_df.columns for col in required_cols):
        raise ValueError(f"Missing required columns. Need {required_cols}, got {region_df.columns.tolist()}")
    
    # Calculate methylation using replicates
    methylation_data = calculate_methylation_levels_with_replicates(
        region_df,
        medip_dir,
        cell_type_prefix,
        genome_fasta
    )
    
    # Add region information to methylation data
    methylation_data = methylation_data.copy()
    methylation_data['start'] = region_df['start'].values
    methylation_data['end'] = region_df['end'].values
    
    # Add quality metrics
    methylation_data = add_methylation_qc_metrics(methylation_data)
    
    # Validate results
    methylation_data = validate_methylation_levels(methylation_data)
    
    # Remove temporary columns if needed
    if 'start' in methylation_data.columns:
        methylation_data = methylation_data.drop(['start', 'end'], axis=1)
    
    return methylation_data

def process_chunk(chunk: pd.DataFrame, medip_dir: str, cell_type_prefix: str, genome_fasta: str) -> pd.DataFrame:
    """Process a chunk of regions for methylation calculation"""
    return calculate_methylation_levels_with_replicates(
        chunk, medip_dir, cell_type_prefix, genome_fasta
    )
