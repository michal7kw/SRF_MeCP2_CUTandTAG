import pandas as pd
import matplotlib.pyplot as plt
from gtfparse import read_gtf
import numpy as np
from typing import Tuple, List, Dict
import seaborn as sns
import os
import argparse
from collections import defaultdict
from intervaltree import IntervalTree
from dataclasses import dataclass
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

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

def load_gene_annotations(gtf_path: str) -> Tuple[pd.DataFrame, Dict]:
    """Load gene annotations including alternative promoters."""
    print("Loading gene annotations...")
    df_gtf = read_gtf(gtf_path)
    if hasattr(df_gtf, 'to_pandas'):
        df_gtf = df_gtf.to_pandas()

    # Process genes
    genes = df_gtf[df_gtf["feature"] == "gene"].copy()
    
    # Process transcripts for alternative TSSs
    transcripts = df_gtf[df_gtf["feature"] == "transcript"].copy()
    
    # Create dictionary to store alternative TSSs
    gene_tss = defaultdict(list)
    
    for _, transcript in transcripts.iterrows():
        if transcript['strand'] == '+':
            tss = transcript['start']
        else:
            tss = transcript['end']
        
        gene_tss[transcript['gene_id']].append({
            'tss': tss,
            'chr': transcript['seqname'],
            'strand': transcript['strand']
        })
    
    return genes, gene_tss

def create_promoter_regions(genes: pd.DataFrame, gene_tss: Dict) -> List[PromoterRegion]:
    """Create promoter regions considering alternative TSSs."""
    promoter_regions = []
    
    for _, gene in genes.iterrows():
        if gene['gene_id'] not in gene_tss:
            continue
        
        for tss_info in gene_tss[gene['gene_id']]:
            if tss_info['strand'] == '+':
                start = max(0, tss_info['tss'] - UPSTREAM_WINDOW)
                end = tss_info['tss'] + DOWNSTREAM_WINDOW
            else:
                start = max(0, tss_info['tss'] - DOWNSTREAM_WINDOW)
                end = tss_info['tss'] + UPSTREAM_WINDOW
            
            promoter_regions.append(PromoterRegion(
                chr=tss_info['chr'],
                start=start,
                end=end,
                strand=tss_info['strand'],
                gene_name=gene['gene_name'],
                gene_type=gene['gene_type']
            ))
    
    return promoter_regions

def create_interval_trees(peaks_df: pd.DataFrame) -> Dict[str, IntervalTree]:
    """Create interval trees for efficient peak searching."""
    trees = defaultdict(IntervalTree)
    
    for _, peak in peaks_df.iterrows():
        trees[peak['chr']].add(
            interval=IntervalTree.Interval(
                begin=peak['start'],
                end=peak['end'],
                data={'signal': peak['signalValue']}
            )
        )
    
    return trees

def find_peaks_in_promoters(promoter_regions: List[PromoterRegion],
                          peak_trees: Dict[str, IntervalTree]) -> Dict[str, List[float]]:
    """Find peaks in promoter regions and return gene-peak associations."""
    promoter_peaks = defaultdict(list)
    
    for promoter in promoter_regions:
        if promoter.chr not in peak_trees:
            continue
        
        overlaps = peak_trees[promoter.chr].overlap(promoter.start, promoter.end)
        if overlaps:
            # Store the maximum signal value if multiple peaks exist
            max_signal = max(overlap.data['signal'] for overlap in overlaps)
            promoter_peaks[promoter.gene_name].append(max_signal)
    
    return promoter_peaks

def analyze_promoter_binding(args):
    """Main analysis function."""
    # Load gene annotations and create promoter regions
    genes, gene_tss = load_gene_annotations(args.gtf_path)
    promoter_regions = create_promoter_regions(genes, gene_tss)
    
    # Load peak files
    peak_files = {
        'endo': os.path.join(args.peaks_dir, "NPCs_endo_combined.narrowPeak"),
        'exo': os.path.join(args.peaks_dir, "NPCs_exo_combined.narrowPeak")
    }
    
    peaks = {
        sample: load_narrowpeak(file_path)
        for sample, file_path in peak_files.items()
    }
    
    # Create interval trees for efficient searching
    peak_trees = {
        sample: create_interval_trees(peak_df)
        for sample, peak_df in peaks.items()
    }
    
    # Find peaks in promoters
    promoter_peaks = {
        sample: find_peaks_in_promoters(promoter_regions, trees)
        for sample, trees in peak_trees.items()
    }
    
    # Analyze overlap
    all_genes = set(promoter_peaks['endo'].keys()) | set(promoter_peaks['exo'].keys())
    
    results = {
        'common': [],
        'endo_only': [],
        'exo_only': []
    }
    
    for gene in all_genes:
        has_endo = gene in promoter_peaks['endo']
        has_exo = gene in promoter_peaks['exo']
        
        if has_endo and has_exo:
            results['common'].append({
                'gene': gene,
                'endo_signal': max(promoter_peaks['endo'][gene]),
                'exo_signal': max(promoter_peaks['exo'][gene])
            })
        elif has_endo:
            results['endo_only'].append({
                'gene': gene,
                'signal': max(promoter_peaks['endo'][gene])
            })
        else:
            results['exo_only'].append({
                'gene': gene,
                'signal': max(promoter_peaks['exo'][gene])
            })
    
    # Save results
    for category, data in results.items():
        if data:
            df = pd.DataFrame(data)
            df.to_csv(os.path.join(args.output_dir, f'{category}_promoter_peaks.csv'), index=False)
    
    # Create visualization
    plot_distributions(results, args.output_dir)
    
    # Print summary
    print("\nAnalysis Summary:")
    print(f"Common promoter peaks: {len(results['common'])}")
    print(f"Endo-only promoter peaks: {len(results['endo_only'])}")
    print(f"Exo-only promoter peaks: {len(results['exo_only'])}")

def plot_distributions(results: Dict, output_dir: str):
    """Create pie chart of peak distributions."""
    categories = {
        'Common': len(results['common']),
        'Endogenous Only': len(results['endo_only']),
        'Exogenous Only': len(results['exo_only'])
    }
    
    colors = ['#e5e5b5', '#b5d7e5', '#3182bd']
    
    plt.figure(figsize=(10, 8))
    plt.pie(categories.values(), labels=categories.keys(), colors=colors,
            autopct='%1.1f%%', startangle=90)
    plt.title('Distribution of Promoter Peak Binding')
    
    plt.savefig(os.path.join(output_dir, "promoter_peak_distribution.png"),
                dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Analyze promoter peak binding')
    parser.add_argument('--gtf-path', required=True, help='Path to GENCODE GTF file')
    parser.add_argument('--peaks-dir', required=True, help='Directory containing peak files')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    analyze_promoter_binding(args)

if __name__ == "__main__":
    main() 