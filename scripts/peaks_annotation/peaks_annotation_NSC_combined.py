import pandas as pd
import matplotlib.pyplot as plt
from gtfparse import read_gtf
import numpy as np
from typing import Tuple, List
import seaborn as sns
import os
import argparse

def load_gencode_annotations(gtf_path: str) -> pd.DataFrame:
    """Load and process GENCODE GTF file."""
    print("Loading GENCODE annotations...")
    df_gtf = read_gtf(gtf_path)
    # Convert to pandas if it's a polars dataframe
    if hasattr(df_gtf, 'to_pandas'):
        df_gtf = df_gtf.to_pandas()
    # Filter for genes and basic annotations
    genes = df_gtf[df_gtf["feature"] == "gene"].copy()
    return genes

def categorize_region(position: Tuple[str, int, int], 
                     genes: pd.DataFrame,
                     promoter_size: int = 2000) -> str:
    """
    Categorize a genomic region based on its position relative to genes.
    """
    chrom, start, end = position
    
    # Filter genes for the relevant chromosome
    chr_genes = genes[genes["seqname"] == chrom]
    
    if chr_genes.empty:
        return "Distal Intergenic"
    
    for _, gene in chr_genes.iterrows():
        gene_start = gene["start"]
        gene_end = gene["end"]
        strand = gene["strand"]
        
        # Define promoter region based on strand
        if strand == "+":
            promoter_start = gene_start - promoter_size
            promoter_end = gene_start
        else:
            promoter_start = gene_end
            promoter_end = gene_end + promoter_size
            
        # Check if peak overlaps with promoter
        if max(start, promoter_start) <= min(end, promoter_end):
            return "Promoter"
        
        # Check if peak overlaps with gene body
        if max(start, gene_start) <= min(end, gene_end):
            # For simplicity, we're considering the first exon/intron as one category
            if (strand == "+" and abs(start - gene_start) <= 1000) or \
               (strand == "-" and abs(end - gene_end) <= 1000):
                return "1st Intron/Exon"
            return "Gene body"
            
    return "Distal Intergenic"

def load_narrowpeak(file_path: str) -> pd.DataFrame:
    """Load narrowPeak format file."""
    columns = ['chr', 'start', 'end', 'name', 'score', 
              'strand', 'signalValue', 'pValue', 'qValue', 'peak']
    return pd.read_csv(file_path, sep='\t', names=columns)

def analyze_peaks(peaks_file: str, genes: pd.DataFrame) -> dict:
    """Analyze peaks and categorize them into genomic regions."""
    print(f"Analyzing peaks from {peaks_file}...")
    peaks = load_narrowpeak(peaks_file)
    
    categories = {
        "Promoter": 0,
        "1st Intron/Exon": 0,
        "Gene body": 0,
        "5' UTR": 0,
        "3' UTR": 0,
        "Distal Intergenic": 0
    }
    
    for _, peak in peaks.iterrows():
        category = categorize_region(
            (peak["chr"], peak["start"], peak["end"]),
            genes
        )
        categories[category] += 1
    
    return categories

def plot_distributions(endo_data: dict, exo_data: dict, output_dir: str):
    """Create side-by-side pie charts for NPCs endogenous vs exogenous."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
    
    colors = {
        "Promoter": "#e5e5b5",
        "1st Intron/Exon": "#b5d7e5",
        "Gene body": "#3182bd",
        "5' UTR": "#e6550d",
        "3' UTR": "#fd8d3c",
        "Distal Intergenic": "#e7969c"
    }
    
    # Plot data
    datasets = [
        (endo_data, "NPCs Endogenous", ax1),
        (exo_data, "NPCs Exogenous", ax2)
    ]
    
    for data, title, ax in datasets:
        values = list(data.values())
        labels = list(data.keys())
        colors_list = [colors[label] for label in labels]
        
        ax.pie(values, labels=labels, colors=colors_list,
               autopct='%1.1f%%', startangle=90)
        ax.set_title(title)
    
    plt.savefig(os.path.join(output_dir, "peak_distributions_NPCs_combined.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Analyze gene annotations in peaks.')
    parser.add_argument('--gtf-path', required=True, help='Path to GENCODE GTF file')
    parser.add_argument('--peaks-dir', required=True, help='Directory containing narrowPeak files')
    parser.add_argument('--output-dir', required=True, help='Output directory for plots')
    args = parser.parse_args()
    
    # Load gene annotations
    genes = load_gencode_annotations(args.gtf_path)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Define peak files
    peak_files = {
        'endo': os.path.join(args.peaks_dir, "NPCs_endo_combined.narrowPeak"),
        'exo': os.path.join(args.peaks_dir, "NPCs_exo_combined.narrowPeak")
    }
    
    # Process files
    results = {}
    for condition, file_path in peak_files.items():
        if os.path.exists(file_path):
            results[condition] = analyze_peaks(file_path, genes)
        else:
            print(f"Warning: File not found - {file_path}")
    
    # Create plot if both files were processed successfully
    if len(results) == 2:
        plot_distributions(results['endo'], results['exo'], args.output_dir)
        
        # Print statistics
        for condition, categories in results.items():
            print(f"\nDistribution for {condition}:")
            total = sum(categories.values())
            for category, count in categories.items():
                percentage = (count / total) * 100
                print(f"{category}: {count} peaks ({percentage:.1f}%)")

if __name__ == "__main__":
    main() 