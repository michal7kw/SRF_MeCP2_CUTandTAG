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
PROMOTER_UPSTREAM = 2500
PROMOTER_DOWNSTREAM = 500
CpG_WINDOW = 500

def standardize_gene_name(gene_name):
    """Standardize gene names to match between DEA and GTF"""
    if pd.isna(gene_name):
        return None
    
    gene_name = str(gene_name)
    gene_name = gene_name.split('.')[0]
    
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
    
    if (df['end'] <= df['start']).any():
        print(f"Warning: Invalid peak coordinates in {sample_name}")
        df = df[df['end'] > df['start']]
    
    if (df['signalValue'] < 0).any():
        print(f"Warning: Negative signal values in {sample_name}")
        df = df[df['signalValue'] >= 0]
    
    return df

def load_peak_file(filepath, sample_name):
    """Load and validate peak file"""
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
    """Load gene annotations using only primary TSS"""
    print("Loading gene annotations...")
    gtf_file = "../DATA/gencode.vM10.annotation.gtf"
    
    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file not found: {gtf_file}")
    
    # Load gene information
    gene_annotations = pd.read_csv(gtf_file, sep='\t', comment='#',
                                 names=['chr', 'source', 'feature', 'start', 'end',
                                       'score', 'strand', 'frame', 'attributes'])
    
    genes = gene_annotations[gene_annotations['feature'] == 'gene'].copy()
    
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
        return pd.Series(info)
    
    gene_info = genes['attributes'].apply(extract_gene_info)
    genes = pd.concat([genes, gene_info], axis=1)
    
    # Add TSS information based on gene start/end
    genes['tss'] = genes.apply(lambda x: x['start'] if x['strand'] == '+' else x['end'], axis=1)
    
    return genes

def load_cpg_islands(cpg_file):
    """Load CpG islands and create interval trees"""
    print("\nLoading CpG islands...")
    
    cpg_islands = pd.read_csv(cpg_file, sep='\t', header=None,
                             names=['chr', 'start', 'end', 'id', 'cpg_info', 'cpg_count'])
    
    cpg_trees = defaultdict(IntervalTree)
    for _, cpg in cpg_islands.iterrows():
        interval = Interval(cpg['start'], cpg['end'], 
                          {'id': cpg['id'], 'cpg_count': cpg['cpg_count']})
        cpg_trees[cpg['chr']].add(interval)
    
    print(f"Loaded {len(cpg_islands)} CpG islands")
    return cpg_trees

def find_associated_cpg(chrom, start, end, cpg_trees, window=CpG_WINDOW):
    """Find closest CpG island to a genomic region"""
    if chrom not in cpg_trees:
        return None
    
    overlaps = list(cpg_trees[chrom].overlap(start - window, end + window))
    
    if not overlaps:
        return None
    
    # Find the CpG island with maximum overlap
    best_cpg = max(overlaps, key=lambda x: min(x.end, end) - max(x.begin, start))
    
    return {
        'cpg_start': best_cpg.begin,
        'cpg_end': best_cpg.end,
        'cpg_id': str(best_cpg.data['id']),
        'cpg_count': best_cpg.data['cpg_count']
    }

def default_tree():
    """Helper function to create default IntervalTree"""
    return defaultdict(IntervalTree)

def create_peak_trees(peaks_dict):
    """Create interval trees for peak data"""
    peak_trees = defaultdict(default_tree)
    
    for sample, peaks in peaks_dict.items():
        for _, peak in peaks.iterrows():
            interval = Interval(peak['start'], peak['end'],
                              {'signal': peak['signalValue'], 'qvalue': peak['qValue']})
            peak_trees[peak['chr']][sample].add(interval)
    
    return dict(peak_trees)  # Convert to regular dict for better pickling

def process_gene_chunk(gene_chunk, cpg_trees, exo_trees, endo_trees):
    """Process a chunk of genes"""
    chunk_results = []
    
    for _, gene in gene_chunk.iterrows():
        chrom = gene['chr']
        tss = gene['tss']
        strand = gene['strand']
        
        if chrom not in exo_trees or chrom not in endo_trees:
            continue
        
        # Define promoter region based on strand
        if strand == '+':
            promoter_start = max(0, tss - PROMOTER_UPSTREAM)
            promoter_end = tss + PROMOTER_DOWNSTREAM
        else:
            promoter_start = max(0, tss - PROMOTER_DOWNSTREAM)
            promoter_end = tss + PROMOTER_UPSTREAM
        
        # Find associated CpG island
        cpg = find_associated_cpg(chrom, promoter_start, promoter_end, cpg_trees)
        
        if cpg is None:
            continue
        
        # Find peaks in the region
        search_start = min(promoter_start, cpg['cpg_start'])
        search_end = max(promoter_end, cpg['cpg_end'])
        
        # Get overlapping peaks
        exo_peaks = list(exo_trees[chrom]['exo'].overlap(search_start, search_end))
        endo_peaks = list(endo_trees[chrom]['endo'].overlap(search_start, search_end))
        
        # Calculate peak statistics
        exo_stats = {
            'peak_count': len(exo_peaks),
            'mean_signal': np.mean([p.data['signal'] for p in exo_peaks]) if exo_peaks else 0,
            'min_qvalue': min([p.data['qvalue'] for p in exo_peaks]) if exo_peaks else 1
        }
        
        endo_stats = {
            'peak_count': len(endo_peaks),
            'mean_signal': np.mean([p.data['signal'] for p in endo_peaks]) if endo_peaks else 0,
            'min_qvalue': min([p.data['qvalue'] for p in endo_peaks]) if endo_peaks else 1
        }
        
        chunk_results.append({
            'gene': gene['gene_name'],
            'chr': chrom,
            'tss': tss,
            'strand': strand,
            'promoter_start': promoter_start,
            'promoter_end': promoter_end,
            'cpg_start': cpg['cpg_start'],
            'cpg_end': cpg['cpg_end'],
            'cpg_id': cpg['cpg_id'],
            'cpg_count': cpg['cpg_count'],
            'exo_peak_count': exo_stats['peak_count'],
            'exo_mean_signal': exo_stats['mean_signal'],
            'exo_min_qvalue': exo_stats['min_qvalue'],
            'endo_peak_count': endo_stats['peak_count'],
            'endo_mean_signal': endo_stats['mean_signal'],
            'endo_min_qvalue': endo_stats['min_qvalue']
        })
    
    return chunk_results

def analyze_promoter_cpg_peaks(gene_annotations, cpg_trees, exo_trees, endo_trees, n_processes=None):
    """Analyze peaks in promoter regions with CpG island association"""
    if n_processes is None:
        n_processes = max(1, mp.cpu_count() - 1)
    
    chunk_size = max(1, len(gene_annotations) // (n_processes * 4))
    gene_chunks = np.array_split(gene_annotations, len(gene_annotations) // chunk_size + 1)
    
    print(f"\nProcessing {len(gene_annotations)} genes in {len(gene_chunks)} chunks...")
    
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        process_chunk_partial = partial(process_gene_chunk,
                                     cpg_trees=cpg_trees,
                                     exo_trees=exo_trees,
                                     endo_trees=endo_trees)
        
        results = []
        for chunk_result in tqdm(executor.map(process_chunk_partial, gene_chunks),
                               total=len(gene_chunks)):
            results.extend(chunk_result)
    
    return pd.DataFrame(results)

def categorize_peaks(df):
    """Categorize peaks based on presence in exo/endo samples"""
    df['category'] = 'none'
    
    has_exo = (df['exo_peak_count'] > 0) & (df['exo_mean_signal'] > 0)
    has_endo = (df['endo_peak_count'] > 0) & (df['endo_mean_signal'] > 0)
    
    df.loc[has_exo & has_endo, 'category'] = 'common'
    df.loc[has_exo & ~has_endo, 'category'] = 'exo_only'
    df.loc[~has_exo & has_endo, 'category'] = 'endo_only'
    
    return df

def save_results(results_df, output_dir):
    """Save analysis results"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Save full results
    results_df.to_csv(f'{output_dir}/all_results.csv', index=False)
    
    # Save category-specific results
    for category in ['exo_only', 'endo_only', 'common']:
        category_df = results_df[results_df['category'] == category]
        category_df.to_csv(f'{output_dir}/{category}.csv', index=False)
    
    # Save statistics
    stats = results_df['category'].value_counts()
    stats.to_csv(f'{output_dir}/statistics.csv')
    
    print("\nResults Summary:")
    print(stats)

def main():
    parser = argparse.ArgumentParser(description='Analyze promoter peaks with CpG islands')
    parser.add_argument('--working-dir', type=str, required=True)
    parser.add_argument('--results-dir', type=str, required=True)
    parser.add_argument('--peaks-dir', type=str, required=True)
    parser.add_argument('--processes', type=int, default=None)
    args = parser.parse_args()

    os.chdir(args.working_dir)
    
    # Load data
    gene_annotations = load_gene_annotations()
    cpg_trees = load_cpg_islands(args.working_dir + "/DATA/cpg_islands.bed")
    
    # Load peak files
    peaks_exo = {'exo': load_peak_file(f"{args.peaks_dir}/NPCs_exo_combined.narrowPeak", "exo")}
    peaks_endo = {'endo': load_peak_file(f"{args.peaks_dir}/NPCs_endo_combined.narrowPeak", "endo")}
    
    # Create interval trees
    print("\nCreating interval trees...")
    exo_trees = create_peak_trees(peaks_exo)
    endo_trees = create_peak_trees(peaks_endo)
    
    # Analyze peaks
    results = analyze_promoter_cpg_peaks(
        gene_annotations, 
        cpg_trees,
        exo_trees,
        endo_trees,
        n_processes=args.processes
    )
    
    # Categorize and save results
    results = categorize_peaks(results)
    save_results(results, args.results_dir)

if __name__ == "__main__":
    main() 