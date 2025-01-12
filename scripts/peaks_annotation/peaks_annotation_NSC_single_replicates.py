import pandas as pd
import matplotlib.pyplot as plt
from gtfparse import read_gtf
import numpy as np
from typing import Tuple, List
import seaborn as sns
import os
import argparse
import json

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
    
    Args:
        position: Tuple of (chromosome, start, end)
        genes: DataFrame with gene annotations
        promoter_size: Size of promoter region upstream of TSS
    
    Returns:
        String indicating the region category
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

def plot_distributions(npc_endo_reps: list, npc_exo_reps: list, output_dir: str):
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
    
    # Average the replicates
    def average_dicts(dict_list):
        result = {k: 0 for k in dict_list[0].keys()}
        for d in dict_list:
            for k, v in d.items():
                result[k] += v
        for k in result:
            result[k] /= len(dict_list)
        return result
    
    npc_endo_avg = average_dicts(npc_endo_reps)
    npc_exo_avg = average_dicts(npc_exo_reps)
    
    # Plot data
    datasets = [
        (npc_endo_avg, "NPCs Endogenous", ax1),
        (npc_exo_avg, "NPCs Exogenous", ax2)
    ]
    
    for data, title, ax in datasets:
        values = list(data.values())
        labels = list(data.keys())
        colors_list = [colors[label] for label in labels]
        
        ax.pie(values, labels=labels, colors=colors_list,
               autopct='%1.1f%%', startangle=90)
        ax.set_title(title)
    
    plt.savefig(os.path.join(output_dir, "peak_distributions_NPCs.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Analyze gene annotations in peaks.')
    parser.add_argument('--gtf-path', required=True, help='Path to GENCODE GTF file')
    parser.add_argument('--peaks-dir', required=True, help='Directory containing narrowPeak files')
    parser.add_argument('--output-dir', required=True, help='Output directory for plots')
    parser.add_argument('--task-id', type=int, help='SLURM array task ID')
    parser.add_argument('--total-tasks', type=int, help='Total number of SLURM array tasks')
    args = parser.parse_args()
    
    # Load gene annotations
    genes = load_gencode_annotations(args.gtf_path)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.join(args.output_dir, "temp"), exist_ok=True)
    
    # Define all files
    all_files = {
        'npc_endo_1': os.path.join(args.peaks_dir, "NSCv1_narrow_peaks.filtered.narrowPeak"),
        'npc_endo_2': os.path.join(args.peaks_dir, "NSCv2_narrow_peaks.filtered.narrowPeak"),
        'npc_endo_3': os.path.join(args.peaks_dir, "NSCv3_narrow_peaks.filtered.narrowPeak"),
        'npc_exo_1': os.path.join(args.peaks_dir, "NSCM1_narrow_peaks.filtered.narrowPeak"),
        'npc_exo_2': os.path.join(args.peaks_dir, "NSCM2_narrow_peaks.filtered.narrowPeak"),
        'npc_exo_3': os.path.join(args.peaks_dir, "NSCM3_narrow_peaks.filtered.narrowPeak")
    }
    
    # If running as part of array job, process only subset of files
    if args.task_id is not None and args.total_tasks is not None:
        files_list = list(all_files.items())
        chunk_size = len(files_list) // args.total_tasks
        start_idx = args.task_id * chunk_size
        end_idx = start_idx + chunk_size if args.task_id < args.total_tasks - 1 else len(files_list)
        peak_files = dict(files_list[start_idx:end_idx])
    else:
        peak_files = all_files
    
    # Process assigned files
    results = {}
    for condition, file_path in peak_files.items():
        if os.path.exists(file_path):
            results[condition] = analyze_peaks(file_path, genes)
        else:
            print(f"Warning: File not found - {file_path}")
    
    # Save partial results
    if args.task_id is not None:
        results_file = os.path.join(args.output_dir, "temp", f"results_{args.task_id}.json")
        with open(results_file, 'w') as f:
            json.dump(results, f)
    else:
        # If running in single mode, create the plot directly
        endo_reps = [results[k] for k in ['npc_endo_1', 'npc_endo_2', 'npc_endo_3'] 
                     if k in results]
        exo_reps = [results[k] for k in ['npc_exo_1', 'npc_exo_2', 'npc_exo_3'] 
                    if k in results]
        
        if endo_reps and exo_reps:
            plot_distributions(endo_reps, exo_reps, args.output_dir)
            
        # Print statistics
        for condition, categories in results.items():
            print(f"\nDistribution for {condition}:")
            total = sum(categories.values())
            for category, count in categories.items():
                percentage = (count / total) * 100
                print(f"{category}: {count} peaks ({percentage:.1f}%)")

if __name__ == "__main__":
    main() 