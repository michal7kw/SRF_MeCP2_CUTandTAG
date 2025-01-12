#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy import stats
import os
import argparse
from tqdm import tqdm

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Analyze enrichment with combined replicate data')
    parser.add_argument('--exo-peaks', required=True, help='Path to combined exogenous peaks file')
    parser.add_argument('--endo-peaks', required=True, help='Path to combined endogenous peaks file')
    parser.add_argument('--dea', required=True, help='Path to DEA results file')
    parser.add_argument('--gtf', required=True, help='Path to gene annotations GTF file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--promoter-window', type=int, default=2000, 
                       help='Size of promoter window (default: 2000bp)')
    return parser.parse_args()

def load_data(args):
    """Load all required data files"""
    print("Loading data files...")
    
    # Load peaks with correct column names
    peak_columns = ['chrom', 'start', 'end', 'name', 'score', 
                   'strand', 'signalValue', 'pValue', 'qValue', 'peak']
    
    # Load peaks and rename chromosome column
    exo_peaks = pd.read_csv(args.exo_peaks, sep='\t', names=peak_columns)
    endo_peaks = pd.read_csv(args.endo_peaks, sep='\t', names=peak_columns)
    
    # Rename chromosome column to match GTF
    exo_peaks = exo_peaks.rename(columns={'chrom': 'chr'})
    endo_peaks = endo_peaks.rename(columns={'chrom': 'chr'})
    
    # Load DEA results
    dea = pd.read_csv(args.dea)
    
    # Load and process GTF
    gene_info = process_gtf(args.gtf)
    
    return exo_peaks, endo_peaks, dea, gene_info

def process_gtf(gtf_file):
    """Process GTF file to extract gene information"""
    print("Processing GTF file...")
    
    gene_info = pd.read_csv(gtf_file, sep='\t', comment='#',
                           names=['chr', 'source', 'feature', 'start', 'end',
                                 'score', 'strand', 'frame', 'attributes'])
    
    # Filter for genes only
    gene_info = gene_info[gene_info['feature'] == 'gene']
    
    # Extract gene name from attributes
    gene_info['gene_name'] = gene_info['attributes'].str.extract('gene_name "([^"]*)"')
    
    return gene_info[['chr', 'start', 'end', 'strand', 'gene_name']]

def categorize_genes(dea, log2fc_threshold=0.5, padj_threshold=0.05):
    """Categorize genes based on differential expression"""
    dea['regulation'] = 'unchanged'
    
    # Up-regulated
    dea.loc[(dea['log2FoldChange'] > log2fc_threshold) & 
            (dea['padj'] < padj_threshold), 'regulation'] = 'up'
    
    # Down-regulated
    dea.loc[(dea['log2FoldChange'] < -log2fc_threshold) & 
            (dea['padj'] < padj_threshold), 'regulation'] = 'down'
    
    return dea

def associate_peaks_with_genes(peaks, gene_info, promoter_window):
    """Associate peaks with genes based on promoter proximity"""
    peak_gene_associations = []
    
    for _, gene in tqdm(gene_info.iterrows(), total=len(gene_info), 
                       desc="Associating peaks with genes"):
        # Define promoter region based on strand
        if gene['strand'] == '+':
            promoter_start = max(0, gene['start'] - promoter_window)
            promoter_end = gene['start'] + promoter_window
        else:
            promoter_start = max(0, gene['end'] - promoter_window)
            promoter_end = gene['end'] + promoter_window
        
        # Find overlapping peaks
        overlapping_peaks = peaks[
            (peaks['chr'] == gene['chr']) &
            (peaks['start'] <= promoter_end) &
            (peaks['end'] >= promoter_start)
        ]
        
        # Add associations
        for _, peak in overlapping_peaks.iterrows():
            peak_gene_associations.append({
                'chr': peak['chr'],
                'peak_start': peak['start'],
                'peak_end': peak['end'],
                'peak_signal': peak['signalValue'],
                'peak_qvalue': peak['qValue'],
                'gene_name': gene['gene_name'],
                'gene_strand': gene['strand'],
                'distance_to_tss': min(abs(peak['start'] - gene['start']), 
                                    abs(peak['end'] - gene['start']))
            })
    
    return pd.DataFrame(peak_gene_associations)

def categorize_peaks(exo_associations, endo_associations, 
                    signal_ratio_threshold=1.5):
    """Categorize peaks based on binding patterns"""
    # Combine associations
    exo_associations['source'] = 'exo'
    endo_associations['source'] = 'endo'
    all_associations = pd.concat([exo_associations, endo_associations])
    
    # Group by genomic location and gene
    grouped = all_associations.groupby(
        ['chr', 'peak_start', 'peak_end', 'gene_name']
    ).agg({
        'peak_signal': lambda x: dict(zip(x.index, x)),
        'source': list
    }).reset_index()
    
    # Categorize peaks
    def determine_category(row):
        sources = set(row['source'])
        if len(sources) == 1:
            return f"{row['source'][0]}-only"
        else:
            exo_signal = row['peak_signal'].get('exo', 0)
            endo_signal = row['peak_signal'].get('endo', 0)
            ratio = exo_signal / endo_signal if endo_signal > 0 else float('inf')
            
            if ratio > signal_ratio_threshold:
                return 'both-exo-enriched'
            elif ratio < 1/signal_ratio_threshold:
                return 'both-endo-enriched'
            else:
                return 'both-similar'
    
    grouped['binding_category'] = grouped.apply(determine_category, axis=1)
    return grouped

def main():
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data
    exo_peaks, endo_peaks, dea, gene_info = load_data(args)
    
    # Categorize genes by expression
    dea = categorize_genes(dea)
    
    # Associate peaks with genes
    print("\nProcessing exogenous peaks...")
    exo_associations = associate_peaks_with_genes(
        exo_peaks, gene_info, args.promoter_window)
    
    print("\nProcessing endogenous peaks...")
    endo_associations = associate_peaks_with_genes(
        endo_peaks, gene_info, args.promoter_window)
    
    # Categorize peaks
    print("\nCategorizing peaks...")
    categorized_peaks = categorize_peaks(exo_associations, endo_associations)
    
    # Merge with gene regulation information
    categorized_peaks = categorized_peaks.merge(
        dea[['gene_name', 'regulation', 'log2FoldChange', 'padj']],
        on='gene_name', how='left'
    )
    
    # Save results by category and regulation
    print("\nSaving results...")
    for binding_cat in categorized_peaks['binding_category'].unique():
        for reg_cat in categorized_peaks['regulation'].unique():
            subset = categorized_peaks[
                (categorized_peaks['binding_category'] == binding_cat) &
                (categorized_peaks['regulation'] == reg_cat)
            ]
            
            if not subset.empty:
                filename = f"{binding_cat}_{reg_cat}_peaks.csv"
                subset.to_csv(os.path.join(args.output_dir, filename), index=False)
                print(f"Saved {len(subset)} peaks to {filename}")
    
    # Save summary statistics
    summary = categorized_peaks.groupby(
        ['binding_category', 'regulation']
    ).size().unstack(fill_value=0)
    
    summary.to_csv(os.path.join(args.output_dir, 'summary_statistics.csv'))
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()