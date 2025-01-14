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

def load_cpg_islands(cpg_path: str) -> pd.DataFrame:
    """Load CpG islands annotation."""
    columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 
              'thickStart', 'thickEnd', 'rgb']
    return pd.read_csv(cpg_path, sep='\t', names=columns, comment='#')

def categorize_region(position: Tuple[str, int, int], 
                     genes: pd.DataFrame,
                     cpg_islands: pd.DataFrame,
                     promoter_size: int = 2000) -> str:
    """
    Categorize a genomic region based on its position relative to genes and CpG islands.
    """
    chrom, start, end = position
    
    # Filter genes and CpG islands for the relevant chromosome first
    chr_genes = genes[genes["seqname"] == chrom]
    cpg_chr = cpg_islands[cpg_islands["chr"] == chrom]
    
    if chr_genes.empty:
        return "Intergenic"
    
    # First check if the peak is in a promoter region
    in_promoter = False
    current_gene = None
    promoter_region = None
    
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
            in_promoter = True
            promoter_region = (promoter_start, promoter_end)
            current_gene = (gene_start, gene_end)
            break
    
    # If in promoter, check for CpG island overlap
    if in_promoter:
        for _, cpg in cpg_chr.iterrows():
            if max(start, cpg["start"]) <= min(end, cpg["end"]):
                return "Promoter CpG Island"
        return "Promoter non-CpG"
    
    # If not in promoter, check for gene body
    for _, gene in chr_genes.iterrows():
        if max(start, gene["start"]) <= min(end, gene["end"]):
            # Check if it's a CpG island in gene body
            for _, cpg in cpg_chr.iterrows():
                if max(start, cpg["start"]) <= min(end, cpg["end"]):
                    return "Other CpG Island"
            return "Gene body"
    
    # If we get here, it's intergenic
    # Check if it's a CpG island in intergenic region
    for _, cpg in cpg_chr.iterrows():
        if max(start, cpg["start"]) <= min(end, cpg["end"]):
            return "Other CpG Island"
    
    return "Intergenic"

def load_narrowpeak(file_path: str) -> pd.DataFrame:
    """Load narrowPeak format file."""
    columns = ['chr', 'start', 'end', 'name', 'score', 
              'strand', 'signalValue', 'pValue', 'qValue', 'peak']
    return pd.read_csv(file_path, sep='\t', names=columns)

def analyze_peaks(peaks_file: str, genes: pd.DataFrame, cpg_islands: pd.DataFrame) -> dict:
    """Analyze peaks and categorize them into genomic regions."""
    print(f"Analyzing peaks from {peaks_file}...")
    peaks = load_narrowpeak(peaks_file)
    
    categories = {
        "Promoter CpG Island": 0,
        "Other CpG Island": 0,
        "Promoter non-CpG": 0,
        "Gene body": 0,
        "Intergenic": 0
    }
    
    for _, peak in peaks.iterrows():
        category = categorize_region(
            (peak["chr"], peak["start"], peak["end"]),
            genes,
            cpg_islands
        )
        categories[category] += 1
    
    return categories

def plot_distributions(endo_data: dict, exo_data: dict, output_dir: str):
    """Create side-by-side pie charts for NPCs endogenous vs exogenous."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
    
    colors = {
        "Promoter CpG Island": "#e5e5b5",
        "Other CpG Island": "#b5d7e5",
        "Promoter non-CpG": "#3182bd",
        "Gene body": "#fd8d3c",
        "Intergenic": "#e7969c"
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
    
    plt.savefig(os.path.join(output_dir, "peak_distributions_NPCs_combined_cpg_specific.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Analyze gene annotations in peaks.')
    parser.add_argument('--gtf-path', required=True, help='Path to GENCODE GTF file')
    parser.add_argument('--peaks-dir', required=True, help='Directory containing narrowPeak files')
    parser.add_argument('--output-dir', required=True, help='Output directory for plots')
    parser.add_argument('--cpg-path', required=True, help='Path to CpG islands annotation file')
    parser.add_argument('--peak-type', required=True, choices=['narrow', 'broad'], help='Type of peaks to process: narrow or broad')
    args = parser.parse_args()
    
    # Load gene annotations
    genes = load_gencode_annotations(args.gtf_path)
    
    # Load CpG islands
    cpg_islands = load_cpg_islands(args.cpg_path)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Define peak files based on peak type
    peak_files = {
        'endo': os.path.join(args.peaks_dir, f"NPCs_endo_combined.{args.peak_type}Peak"),
        'exo': os.path.join(args.peaks_dir, f"NPCs_exo_combined.{args.peak_type}Peak")
    }
    
    # Process files
    results = {}
    for condition, file_path in peak_files.items():
        if os.path.exists(file_path):
            results[condition] = analyze_peaks(file_path, genes, cpg_islands)
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