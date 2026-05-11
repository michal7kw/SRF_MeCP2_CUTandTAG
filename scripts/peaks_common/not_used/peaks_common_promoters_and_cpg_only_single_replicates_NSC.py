import sys
import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import argparse
from intervaltree import IntervalTree, Interval
import pybedtools
from collections import defaultdict
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import time

# Constants
PROMOTER_UPSTREAM = 2500  # Increased upstream distance
PROMOTER_DOWNSTREAM = 500  # Reduced downstream distance
CpG_WINDOW = 500  # Distance to consider CpG island association
CpG_MERGE_DISTANCE = 100  # Distance within which to merge CpG islands

def standardize_gene_name(gene_name):
    """Standardize gene names to match between DEA and GTF"""
    if pd.isna(gene_name):
        return None
    
    # Convert to string if not already
    gene_name = str(gene_name)
    
    # Remove version numbers if present (e.g., Gene.1 -> Gene)
    gene_name = gene_name.split('.')[0]
    
    # Remove common prefixes/suffixes that might differ between annotations
    prefixes = ['gene-', 'Gene-', 'GENE-']
    for prefix in prefixes:
        if gene_name.startswith(prefix):
            gene_name = gene_name[len(prefix):]
    
    return gene_name.strip()

def validate_peak_file(df, sample_name):
    """Validate peak file contents"""
    if df.empty:
        print(f"Warning: Empty peak file for {sample_name}")
        return df
    
    required_cols = ['chr', 'start', 'end', 'signalValue']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Warning: Missing columns in {sample_name} peak file: {missing_cols}")
        return pd.DataFrame()
    
    # Check for invalid values
    if (df['end'] <= df['start']).any():
        print(f"Warning: Invalid peak coordinates in {sample_name}")
        df = df[df['end'] > df['start']]
    
    if (df['signalValue'] < 0).any():
        print(f"Warning: Negative signal values in {sample_name}")
        df = df[df['signalValue'] >= 0]
    
    return df

def load_peak_file(filepath, sample_name):
    """Load and validate peak file with enhanced error handling"""
    if not os.path.exists(filepath):
        print(f"Warning: Peak file not found: {filepath}")
        return pd.DataFrame()
    
    try:
        df = pd.read_csv(filepath, sep='\t', header=None,
                        names=['chr', 'start', 'end', 'name', 'score',
                              'strand', 'signalValue', 'pValue',
                              'qValue', 'peak'])
        return validate_peak_file(df, sample_name)
    except Exception as e:
        print(f"Error loading peak file {filepath}: {str(e)}")
        return pd.DataFrame()

def load_gene_annotations():
    """Load gene annotations including alternative promoters"""
    print("Loading gene annotations...")
    gtf_file = "../DATA/gencode.vM10.annotation.gtf"
    
    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file not found: {gtf_file}")
    
    # Load both gene and transcript information
    gene_annotations = pd.read_csv(gtf_file, sep='\t', comment='#',
                                 names=['chr', 'source', 'feature', 'start', 'end',
                                       'score', 'strand', 'frame', 'attributes'])
    
    # Process genes and transcripts separately
    genes = gene_annotations[gene_annotations['feature'] == 'gene'].copy()
    transcripts = gene_annotations[gene_annotations['feature'] == 'transcript'].copy()
    
    # Extract gene information as before
    def extract_gene_info(attr):
        info = {}
        for field in attr.split(';'):
            field = field.strip()
            if field.startswith('gene_name'):
                info['gene_name'] = field.split('"')[1]
            elif field.startswith('gene_id'):
                info['gene_id'] = field.split('"')[1]
            elif field.startswith('gene_type'):
                info['gene_type'] = field.split('"')[1]
            elif field.startswith('transcript_id'):
                info['transcript_id'] = field.split('"')[1]
        return pd.Series(info)
    
    # Process genes
    gene_info = genes['attributes'].apply(extract_gene_info)
    genes = pd.concat([genes, gene_info], axis=1)
    
    # Process transcripts to get alternative TSSs
    transcript_info = transcripts['attributes'].apply(extract_gene_info)
    transcripts = pd.concat([transcripts, transcript_info], axis=1)
    
    # Group transcripts by gene and collect all TSSs
    gene_tss = defaultdict(list)
    for _, transcript in transcripts.iterrows():
        if transcript['strand'] == '+':
            tss = transcript['start']
        else:
            tss = transcript['end']
        gene_tss[transcript['gene_id']].append({
            'tss': tss,
            'chr': transcript['chr'],
            'strand': transcript['strand']
        })
    
    # Add TSS information to gene annotations
    genes['alternative_tss'] = genes['gene_id'].map(gene_tss)
    
    return genes, gene_tss

def load_cpg_islands(cpg_file):
    """Load CpG islands and create interval trees for efficient lookup"""
    print("\nLoading CpG islands...")
    
    cpg_islands = pd.read_csv(cpg_file, sep='\t', header=None,
                             names=['chr', 'start', 'end', 'id', 'cpg_info', 'cpg_count'])
    
    # Create interval trees for each chromosome
    cpg_trees = defaultdict(IntervalTree)
    for _, cpg in cpg_islands.iterrows():
        interval = Interval(cpg['start'], cpg['end'], 
                          {'id': cpg['id'], 'cpg_count': cpg['cpg_count']})
        cpg_trees[cpg['chr']].add(interval)
    
    print(f"Loaded {len(cpg_islands)} CpG islands")
    return cpg_trees

def find_associated_cpg(chrom, start, end, cpg_trees, window=CpG_WINDOW):
    """Find and merge CpG islands associated with a given genomic region"""
    if chrom not in cpg_trees:
        return None
    
    # Expand search region by window size
    overlaps = list(cpg_trees[chrom].overlap(start - window, end + window))
    
    if not overlaps:
        return None
    
    # If only one CpG island, return it
    if len(overlaps) == 1:
        cpg = overlaps[0]
        return {
            'cpg_start': cpg.begin,
            'cpg_end': cpg.end,
            'cpg_id': str(cpg.data['id']),
            'cpg_count': cpg.data['cpg_count']
        }
    
    # Sort overlapping CpG islands by position
    overlaps.sort(key=lambda x: x.begin)
    
    # Merge nearby/overlapping CpG islands
    merged_cpgs = []
    current_cpg = overlaps[0]
    current_start = current_cpg.begin
    current_end = current_cpg.end
    current_count = current_cpg.data['cpg_count']
    current_ids = [str(current_cpg.data['id'])]
    
    for cpg in overlaps[1:]:
        if cpg.begin - current_end <= CpG_MERGE_DISTANCE:
            # Merge this CpG island with current one
            current_end = max(current_end, cpg.end)
            current_count += cpg.data['cpg_count']
            current_ids.append(str(cpg.data['id']))
        else:
            # Save current merged region and start new one
            merged_cpgs.append({
                'cpg_start': current_start,
                'cpg_end': current_end,
                'cpg_id': ','.join(current_ids),
                'cpg_count': current_count
            })
            current_start = cpg.begin
            current_end = cpg.end
            current_count = cpg.data['cpg_count']
            current_ids = [str(cpg.data['id'])]
    
    # Add the last merged region
    merged_cpgs.append({
        'cpg_start': current_start,
        'cpg_end': current_end,
        'cpg_id': ','.join(current_ids),
        'cpg_count': current_count
    })
    
    # Return the merged region that covers the largest portion of the target region
    best_cpg = max(merged_cpgs, key=lambda x: min(x['cpg_end'], end) - max(x['cpg_start'], start))
    return best_cpg

def default_interval_tree():
    """Helper function to replace lambda in defaultdict"""
    return defaultdict(IntervalTree)

def process_chromosome_peaks(data):
    """Process peaks for a single chromosome in parallel"""
    chrom, sample_peaks = data
    trees = defaultdict(IntervalTree)
    
    for sample, peaks_df in sample_peaks.items():
        for _, peak in peaks_df.iterrows():
            interval = Interval(peak['start'], peak['end'],
                              {'signal': peak['signalValue'], 'qvalue': peak['qValue']})
            trees[sample].add(interval)
    
    return chrom, dict(trees)

def create_peak_trees(peaks_dict):
    """Create interval trees for peak data using parallel processing"""
    print("Organizing peaks by chromosome...")
    # Group peaks by chromosome
    chrom_peaks = defaultdict(dict)
    for sample, peaks in peaks_dict.items():
        for chrom in peaks['chr'].unique():
            if chrom not in chrom_peaks:
                chrom_peaks[chrom] = {}
            chrom_peaks[chrom][sample] = peaks[peaks['chr'] == chrom]
    
    # Process chromosomes in parallel
    print(f"Processing {len(chrom_peaks)} chromosomes in parallel...")
    with ProcessPoolExecutor(max_workers=mp.cpu_count()) as executor:
        futures = []
        for chrom, sample_peaks in chrom_peaks.items():
            futures.append(executor.submit(process_chromosome_peaks, (chrom, sample_peaks)))
        
        # Collect results
        peak_trees = {}
        for future in tqdm(futures, desc="Creating interval trees"):
            chrom, trees = future.result()
            peak_trees[chrom] = trees
    
    return peak_trees

def process_gene_chunk(gene_chunk, cpg_trees, exo_trees, endo_trees):
    """Process a chunk of genes considering alternative promoters"""
    chunk_results = []
    
    for _, gene in gene_chunk.iterrows():
        if not gene['alternative_tss']:
            continue
            
        # Process each TSS for the gene
        for tss_info in gene['alternative_tss']:
            chrom = tss_info['chr']
            tss = tss_info['tss']
            strand = tss_info['strand']
            
            if chrom not in exo_trees or chrom not in endo_trees:
                continue
            
            # Define asymmetric promoter region based on strand
            if strand == '+':
                promoter_start = max(0, tss - PROMOTER_UPSTREAM)
                promoter_end = tss + PROMOTER_DOWNSTREAM
            else:
                promoter_start = max(0, tss - PROMOTER_DOWNSTREAM)
                promoter_end = tss + PROMOTER_UPSTREAM
            
            # Find associated CpG island
            cpg = find_associated_cpg(gene['chr'], promoter_start, promoter_end, cpg_trees)
            
            if cpg is None:
                continue
            
            # Find peaks in the region
            exo_peaks = []
            endo_peaks = []
            
            # Search region includes both promoter and CpG island
            search_start = min(promoter_start, cpg['cpg_start'])
            search_end = max(promoter_end, cpg['cpg_end'])
            
            # Collect exo peaks with sample tracking
            for sample, tree in exo_trees.get(gene['chr'], defaultdict(IntervalTree)).items():
                overlaps = list(tree.overlap(search_start, search_end))
                if overlaps:
                    exo_peaks.extend([{
                        'sample': sample,
                        'start': o.begin,
                        'end': o.end,
                        'signal': o.data['signal'],
                        'qvalue': o.data['qvalue']
                    } for o in overlaps])
            
            # Collect endo peaks with sample tracking
            for sample, tree in endo_trees.get(gene['chr'], defaultdict(IntervalTree)).items():
                overlaps = list(tree.overlap(search_start, search_end))
                if overlaps:
                    endo_peaks.extend([{
                        'sample': sample,
                        'start': o.begin,
                        'end': o.end,
                        'signal': o.data['signal'],
                        'qvalue': o.data['qvalue']
                    } for o in overlaps])
            
            # Group peaks by sample to count replicates
            exo_by_sample = defaultdict(list)
            endo_by_sample = defaultdict(list)
            
            for peak in exo_peaks:
                exo_by_sample[peak['sample']].append(peak)
            for peak in endo_peaks:
                endo_by_sample[peak['sample']].append(peak)
            
            # Calculate summary statistics
            exo_stats = {
                'peak_count': len(exo_peaks),
                'replicate_count': len(exo_by_sample),
                'mean_signal': np.mean([p['signal'] for p in exo_peaks]) if exo_peaks else 0,
                'min_qvalue': min([p['qvalue'] for p in exo_peaks]) if exo_peaks else 1,
                'samples': ';'.join(sorted(exo_by_sample.keys()))
            }
            
            endo_stats = {
                'peak_count': len(endo_peaks),
                'replicate_count': len(endo_by_sample),
                'mean_signal': np.mean([p['signal'] for p in endo_peaks]) if endo_peaks else 0,
                'min_qvalue': min([p['qvalue'] for p in endo_peaks]) if endo_peaks else 1,
                'samples': ';'.join(sorted(endo_by_sample.keys()))
            }
            
            chunk_results.append({
                'gene': gene['gene_name'],
                'chr': gene['chr'],
                'tss': tss,
                'strand': gene['strand'],
                'promoter_start': promoter_start,
                'promoter_end': promoter_end,
                'cpg_start': cpg['cpg_start'],
                'cpg_end': cpg['cpg_end'],
                'cpg_id': cpg['cpg_id'],
                'cpg_count': cpg['cpg_count'],
                'exo_peak_count': exo_stats['peak_count'],
                'exo_replicate_count': exo_stats['replicate_count'],
                'exo_mean_signal': exo_stats['mean_signal'],
                'exo_min_qvalue': exo_stats['min_qvalue'],
                'exo_samples': exo_stats['samples'],
                'endo_peak_count': endo_stats['peak_count'],
                'endo_replicate_count': endo_stats['replicate_count'],
                'endo_mean_signal': endo_stats['mean_signal'],
                'endo_min_qvalue': endo_stats['min_qvalue'],
                'endo_samples': endo_stats['samples']
            })
    
    return chunk_results

def analyze_promoter_cpg_peaks(gene_annotations, cpg_trees, exo_trees, endo_trees, n_processes=None):
    """Analyze peaks in promoter regions with CpG island association using parallel processing"""
    if n_processes is None:
        n_processes = max(1, mp.cpu_count() - 1)
    
    # Convert trees to regular dicts for pickling
    cpg_trees = dict(cpg_trees)
    exo_trees = dict(exo_trees)
    endo_trees = dict(endo_trees)
    
    # Split genes into chunks for parallel processing
    chunk_size = max(1, len(gene_annotations) // (n_processes * 4))
    gene_chunks = np.array_split(gene_annotations, len(gene_annotations) // chunk_size + 1)
    
    print(f"\nProcessing {len(gene_annotations)} genes in {len(gene_chunks)} chunks using {n_processes} processes")
    
    # Process chunks in parallel
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        process_chunk_partial = partial(process_gene_chunk,
                                     cpg_trees=cpg_trees,
                                     exo_trees=exo_trees,
                                     endo_trees=endo_trees)
        
        results = []
        futures = list(tqdm(
            executor.map(process_chunk_partial, gene_chunks),
            total=len(gene_chunks),
            desc="Processing gene chunks"
        ))
        
        # Flatten results
        results = [item for sublist in futures for item in sublist]
    
    return pd.DataFrame(results)

def categorize_peaks(df):
    """
    Categorize peaks based on presence in exo/endo samples with refined criteria:
    - Common: Strong peaks in both exo and endo (≥2 replicates each)
    - Exo-only: Strong peaks in exo (≥2 replicates) but weak/no peaks in endo
    - Endo-only: Strong peaks in endo (≥2 replicates) but weak/no peaks in exo
    """
    # Initialize category as none
    df['category'] = 'none'
    
    # Define strong peak criteria using only replicate count and signal strength
    strong_exo = (
        (df['exo_replicate_count'] >= 2) 
        # & (df['exo_mean_signal'] > 5.0)
    )
    
    strong_endo = (
        (df['endo_replicate_count'] >= 2) 
        # & (df['endo_mean_signal'] > 5.0)
    )
    
    # Categorize peaks
    df.loc[strong_exo & strong_endo, 'category'] = 'common'
    df.loc[strong_exo & ~strong_endo, 'category'] = 'exo_only'
    df.loc[~strong_exo & strong_endo, 'category'] = 'endo_only'
    
    return df

def save_categorized_results(results_df, output_dir):
    """Save categorized results to separate files"""
    # Create category-specific dataframes
    exo_only_df = results_df[results_df['category'] == 'exo_only']
    endo_only_df = results_df[results_df['category'] == 'endo_only']
    common_df = results_df[results_df['category'] == 'common']
    
    # Save to separate files
    exo_only_df.to_csv(f'{output_dir}/exo_only.csv', index=False)
    endo_only_df.to_csv(f'{output_dir}/endo_only.csv', index=False)
    common_df.to_csv(f'{output_dir}/common.csv', index=False)
    
    # Print statistics
    print("\nDetailed Peak Statistics:")
    print(f"Exo-only peaks: {len(exo_only_df)} regions")
    print(f"Endo-only peaks: {len(endo_only_df)} regions")
    print(f"Common peaks: {len(common_df)} regions")
    print(f"No significant peaks: {len(results_df[results_df['category'] == 'none'])} regions")

def main():
    parser = argparse.ArgumentParser(description='Analyze promoter peaks with CpG islands')
    parser.add_argument('--working-dir', type=str, required=True)
    parser.add_argument('--results-dir', type=str, required=True)
    parser.add_argument('--processes', type=int, default=None,
                       help='Number of processes to use (default: CPU count - 1)')
    args = parser.parse_args()

    os.chdir(args.working_dir)
    os.makedirs(args.results_dir, exist_ok=True)

    # Load data
    print("\nLoading data...")
    gene_annotations, _ = load_gene_annotations()
    cpg_trees = load_cpg_islands(args.working_dir + "/DATA/cpg_islands.bed")
    
    # Load peaks
    exo_samples = ['NSCv1', 'NSCv2', 'NSCv3']
    endo_samples = ['NSCM1', 'NSCM2', 'NSCM3']
    
    print("\nLoading peak files...")
    peaks_exo = {sample: load_peak_file(f"{args.results_dir}/peaks/{sample}_peaks.narrowPeak", sample)
                 for sample in exo_samples}
    peaks_endo = {sample: load_peak_file(f"{args.results_dir}/peaks/{sample}_peaks.narrowPeak", sample)
                  for sample in endo_samples}
    
    # Add debugging information
    print("\nPeak file statistics:")
    for sample, peaks in peaks_exo.items():
        print(f"{sample}: {len(peaks)} peaks")
    for sample, peaks in peaks_endo.items():
        print(f"{sample}: {len(peaks)} peaks")
    
    # Create interval trees for efficient searching
    print("\nCreating interval trees...")
    start_time = time.time()
    
    print("Processing exogenous peaks...")
    exo_trees = create_peak_trees(peaks_exo)
    
    print("Processing endogenous peaks...")
    endo_trees = create_peak_trees(peaks_endo)
    
    print(f"Interval trees created in {time.time() - start_time:.2f} seconds")
    
    # Print tree statistics
    print("\nInterval tree statistics:")
    total_exo_intervals = 0
    total_endo_intervals = 0
    
    for chr in exo_trees:
        for sample in exo_trees[chr]:
            intervals = len(exo_trees[chr][sample])
            total_exo_intervals += intervals
            print(f"Exo - {chr} - {sample}: {intervals} intervals")
    
    for chr in endo_trees:
        for sample in endo_trees[chr]:
            intervals = len(endo_trees[chr][sample])
            total_endo_intervals += intervals
            print(f"Endo - {chr} - {sample}: {intervals} intervals")
    
    print(f"\nTotal intervals - Exo: {total_exo_intervals}, Endo: {total_endo_intervals}")
    
    # Analyze peaks using parallel processing
    results = analyze_promoter_cpg_peaks(
        gene_annotations, 
        cpg_trees, 
        exo_trees, 
        endo_trees,
        n_processes=args.processes
    )
    
    # Categorize and save results
    results = categorize_peaks(results)
    
    # Save all results and category-specific files
    results.to_csv(f'{args.results_dir}/promoter_cpg_peak_analysis.csv', index=False)
    save_categorized_results(results, args.results_dir)
    
    # Print summary
    print("\nAnalysis Summary:")
    print(results['category'].value_counts())
    
    # Save detailed statistics
    stats = {
        'total_genes': len(gene_annotations),
        'genes_with_cpg': len(results),
        'common_peaks': len(results[results['category'] == 'common']),
        'exo_only_peaks': len(results[results['category'] == 'exo_only']),
        'endo_only_peaks': len(results[results['category'] == 'endo_only']),
        'no_peaks': len(results[results['category'] == 'none'])
    }
    
    with open(f'{args.results_dir}/analysis_statistics.txt', 'w') as f:
        for key, value in stats.items():
            f.write(f'{key}: {value}\n')

if __name__ == "__main__":
    main() 