import sys
import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import argparse
from typing import Dict, Tuple, List
from dataclasses import dataclass

# Constants for asymmetric promoter region
UPSTREAM_WINDOW = 2500  # Upstream from TSS
DOWNSTREAM_WINDOW = 500  # Downstream from TSS

@dataclass
class PromoterRegion:
    chr: str
    start: int
    end: int
    strand: str
    gene_name: str
    gene_type: str

def load_gene_annotations(gtf_file: str) -> pd.DataFrame:
    """Load gene annotations from GTF file with enhanced processing."""
    print("Loading gene annotations...")
    
    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file not found: {gtf_file}")
    
    # Read GTF file
    gene_annotations = pd.read_csv(gtf_file, sep='\t', comment='#',
                                 names=['chr', 'source', 'feature', 'start', 'end',
                                       'score', 'strand', 'frame', 'attributes'])
    
    # Filter for genes
    gene_annotations = gene_annotations[gene_annotations['feature'] == 'gene']
    
    # Extract gene information from attributes
    def extract_gene_info(attr: str) -> pd.Series:
        info = {}
        for field in attr.split(';'):
            field = field.strip()
            if field.startswith('gene_name'):
                info['gene_name'] = field.split('"')[1]
            elif field.startswith('gene_type'):
                info['gene_type'] = field.split('"')[1]
        return pd.Series(info)
    
    # Process attributes
    gene_info = gene_annotations['attributes'].apply(extract_gene_info)
    gene_annotations = pd.concat([gene_annotations, gene_info], axis=1)
    
    # Define promoter regions based on strand
    gene_annotations['promoter_start'] = np.where(
        gene_annotations['strand'] == '+',
        gene_annotations['start'] - UPSTREAM_WINDOW,
        gene_annotations['end'] - DOWNSTREAM_WINDOW
    )
    
    gene_annotations['promoter_end'] = np.where(
        gene_annotations['strand'] == '+',
        gene_annotations['start'] + DOWNSTREAM_WINDOW,
        gene_annotations['end'] + UPSTREAM_WINDOW
    )
    
    # Ensure promoter starts are not negative
    gene_annotations['promoter_start'] = gene_annotations['promoter_start'].clip(lower=0)
    
    print(f"Loaded {len(gene_annotations)} genes from GTF")
    return gene_annotations

def load_narrowpeak(file_path: str) -> pd.DataFrame:
    """Load narrowPeak format file with validation."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Peak file not found: {file_path}")
    
    columns = ['chr', 'start', 'end', 'name', 'score', 
              'strand', 'signalValue', 'pValue', 'qValue', 'peak']
    
    df = pd.read_csv(file_path, sep='\t', names=columns)
    
    # Basic validation
    df = df[df['end'] > df['start']]
    df = df[df['signalValue'] >= 0]
    
    return df

def find_peaks_in_promoters(peaks_df: pd.DataFrame, 
                          gene_annotations: pd.DataFrame) -> pd.DataFrame:
    """Find peaks that overlap with promoter regions."""
    
    # Create intervals for efficient overlap detection
    promoter_intervals = gene_annotations[['chr', 'promoter_start', 'promoter_end', 
                                         'gene_name', 'gene_type', 'strand']].copy()
    
    # Initialize results list
    promoter_peaks = []
    
    # Group peaks by chromosome for efficiency
    peaks_by_chr = dict(tuple(peaks_df.groupby('chr')))
    
    # Process each chromosome's promoters
    for chr_name, chr_promoters in tqdm(promoter_intervals.groupby('chr')):
        if chr_name not in peaks_by_chr:
            continue
            
        chr_peaks = peaks_by_chr[chr_name]
        
        # Process each promoter
        for _, promoter in chr_promoters.iterrows():
            # Find peaks overlapping with this promoter
            overlapping_peaks = chr_peaks[
                (chr_peaks['start'] <= promoter['promoter_end']) &
                (chr_peaks['end'] >= promoter['promoter_start'])
            ]
            
            if not overlapping_peaks.empty:
                # Aggregate peak information if multiple peaks exist
                peak_info = {
                    'chr': promoter['chr'],
                    'gene_name': promoter['gene_name'],
                    'gene_type': promoter['gene_type'],
                    'strand': promoter['strand'],
                    'promoter_start': promoter['promoter_start'],
                    'promoter_end': promoter['promoter_end'],
                    'peak_count': len(overlapping_peaks),
                    'max_signal': overlapping_peaks['signalValue'].max(),
                    'total_signal': overlapping_peaks['signalValue'].sum(),
                    'peak_starts': ','.join(map(str, overlapping_peaks['start'])),
                    'peak_ends': ','.join(map(str, overlapping_peaks['end']))
                }
                promoter_peaks.append(peak_info)
    
    return pd.DataFrame(promoter_peaks)

def analyze_promoter_binding(args):
    """Main analysis function."""
    # Load gene annotations
    gene_annotations = load_gene_annotations("../DATA/gencode.vM10.annotation.gtf")
    
    # Define peak files
    peak_files = {
        'endo': os.path.join(args.peaks_dir, "NPCs_endo_combined.narrowPeak"),
        'exo': os.path.join(args.peaks_dir, "NPCs_exo_combined.narrowPeak")
    }
    
    # Load peaks
    peaks = {
        sample: load_narrowpeak(file_path)
        for sample, file_path in peak_files.items()
    }
    
    # Find peaks in promoters for each sample
    promoter_peaks = {
        sample: find_peaks_in_promoters(peak_df, gene_annotations)
        for sample, peak_df in peaks.items()
    }
    
    # Analyze overlap
    all_genes = set(promoter_peaks['endo']['gene_name']) | set(promoter_peaks['exo']['gene_name'])
    
    results = {
        'common': [],
        'endo_only': [],
        'exo_only': []
    }
    
    for gene in all_genes:
        endo_data = promoter_peaks['endo'][promoter_peaks['endo']['gene_name'] == gene]
        exo_data = promoter_peaks['exo'][promoter_peaks['exo']['gene_name'] == gene]
        
        if not endo_data.empty and not exo_data.empty:
            # Gene has peaks in both samples
            results['common'].append({
                'gene_name': gene,
                'endo_peaks': endo_data.iloc[0]['peak_count'],
                'exo_peaks': exo_data.iloc[0]['peak_count'],
                'endo_signal': endo_data.iloc[0]['max_signal'],
                'exo_signal': exo_data.iloc[0]['max_signal']
            })
        elif not endo_data.empty:
            results['endo_only'].append({
                'gene_name': gene,
                'peak_count': endo_data.iloc[0]['peak_count'],
                'signal': endo_data.iloc[0]['max_signal']
            })
        else:
            results['exo_only'].append({
                'gene_name': gene,
                'peak_count': exo_data.iloc[0]['peak_count'],
                'signal': exo_data.iloc[0]['max_signal']
            })
    
    # Convert results to DataFrames and save
    for category, data in results.items():
        if data:
            df = pd.DataFrame(data)
            df.to_csv(os.path.join(args.output_dir, f'{category}_promoter_peaks.csv'), index=False)
    
    # Print summary
    print("\nAnalysis Summary:")
    print(f"Common promoter peaks: {len(results['common'])}")
    print(f"Endo-only promoter peaks: {len(results['endo_only'])}")
    print(f"Exo-only promoter peaks: {len(results['exo_only'])}")
def main():
    parser = argparse.ArgumentParser(description='Analyze promoter peak binding')
    parser.add_argument('--peaks-dir', required=True, help='Directory containing peak files')
    parser.add_argument('--output-dir', required=True, help='Directory for output files')
    parser.add_argument('--gtf-path', required=True, help='Path to GTF annotation file')
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    analyze_promoter_binding(args)

if __name__ == "__main__":
    main()