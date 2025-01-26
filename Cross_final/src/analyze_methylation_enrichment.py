#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pyBigWig
import pysam
from typing import List, Tuple, Dict
import os
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

def load_gene_lists(base_dir: str) -> Dict[str, pd.DataFrame]:
    """Load the three gene lists with their peak regions."""
    gene_lists = {}
    for category in ['up', 'down', 'no_deg']:
        file_path = os.path.join(base_dir, f'extended_{category}.csv')
        df = pd.read_csv(file_path)
        # Filter out rows without coordinates
        df = df.dropna(subset=['seqnames', 'start', 'end'])
        # Convert coordinates to integers
        df['start'] = df['start'].astype(float).astype(int)
        df['end'] = df['end'].astype(float).astype(int)
        gene_lists[category] = df
    return gene_lists

def calculate_cpg_density(sequence: str) -> float:
    """Calculate CpG density in a sequence."""
    if len(sequence) < 2:
        return 0.0
    cpg_count = sequence.count('CG')
    return (cpg_count * 2) / len(sequence)

def get_region_signal(bw_file: str, chrom: str, start: int, end: int) -> np.ndarray:
    """Get signal values from bigWig file for a region."""
    with pyBigWig.open(bw_file) as bw:
        try:
            values = np.array(bw.values(chrom, start, end))
            values[np.isnan(values)] = 0
            return values
        except:
            return np.zeros(end - start)

def calculate_region_methylation(region: pd.Series,
                               ip_files: List[str],
                               input_files: List[str],
                               fasta: pysam.FastaFile,
                               window_size: int = 500) -> Tuple[float, float]:
    """Calculate methylation level for a region with surrounding context."""
    
    # Extend region by window_size in both directions
    chrom = region['seqnames']
    center = (region['start'] + region['end']) // 2
    start = max(0, center - window_size)
    end = center + window_size
    
    # Get sequence and calculate CpG density
    try:
        sequence = fasta.fetch(chrom, start, end)
        cpg_density = calculate_cpg_density(sequence)
    except:
        return 0.0, 0.0
    
    # Calculate average signal across replicates
    ip_signals = []
    input_signals = []
    
    for ip_file, input_file in zip(ip_files, input_files):
        ip_signal = get_region_signal(ip_file, chrom, start, end)
        input_signal = get_region_signal(input_file, chrom, start, end)
        
        if len(ip_signal) > 0 and len(input_signal) > 0:
            # Normalize IP by input
            norm_factor = np.mean(input_signal) if np.mean(input_signal) > 0 else 1
            normalized_signal = ip_signal / norm_factor
            ip_signals.append(normalized_signal)
            input_signals.append(input_signal)
    
    if not ip_signals:
        return 0.0, 0.0
    
    # Average across replicates
    mean_ip = np.mean(ip_signals, axis=0)
    mean_input = np.mean(input_signals, axis=0)
    
    # Calculate final methylation score
    methylation_score = np.mean(mean_ip) * cpg_density if cpg_density > 0 else 0
    
    return methylation_score, cpg_density

def calculate_smarcb1_enrichment(region: pd.Series,
                               bm_file_cpm: str,
                               bm_file_rpkm: str,
                               bg_files_cpm: List[str],
                               bg_files_rpkm: List[str],
                               window_size: int = 500) -> Tuple[float, float]:
    """Calculate SMARCB1 enrichment (BM vs BG) for a region using both CPM and RPKM."""
    
    chrom = region['seqnames']
    center = (region['start'] + region['end']) // 2
    start = max(0, center - window_size)
    end = center + window_size
    
    # Get BM signals
    bm_signal_cpm = get_region_signal(bm_file_cpm, chrom, start, end)
    bm_signal_rpkm = get_region_signal(bm_file_rpkm, chrom, start, end)
    
    # Get average BG signals for CPM
    bg_signals_cpm = []
    for bg_file in bg_files_cpm:
        bg_signal = get_region_signal(bg_file, chrom, start, end)
        if len(bg_signal) > 0:
            bg_signals_cpm.append(bg_signal)
    
    # Get average BG signals for RPKM
    bg_signals_rpkm = []
    for bg_file in bg_files_rpkm:
        bg_signal = get_region_signal(bg_file, chrom, start, end)
        if len(bg_signal) > 0:
            bg_signals_rpkm.append(bg_signal)
    
    if not bg_signals_cpm or not bg_signals_rpkm or len(bm_signal_cpm) == 0 or len(bm_signal_rpkm) == 0:
        return 0.0, 0.0
    
    # Calculate enrichment for CPM
    mean_bg_cpm = np.mean(bg_signals_cpm, axis=0)
    bg_mean_cpm = np.mean(mean_bg_cpm)
    enrichment_cpm = np.mean(bm_signal_cpm) / bg_mean_cpm if bg_mean_cpm > 0 else 0.0
    
    # Calculate enrichment for RPKM
    mean_bg_rpkm = np.mean(bg_signals_rpkm, axis=0)
    bg_mean_rpkm = np.mean(mean_bg_rpkm)
    enrichment_rpkm = np.mean(bm_signal_rpkm) / bg_mean_rpkm if bg_mean_rpkm > 0 else 0.0
    
    return enrichment_cpm, enrichment_rpkm

def analyze_methylation_and_smarcb1(gene_lists: Dict[str, pd.DataFrame],
                                  medip_dir: str,
                                  smarcb1_dir: str,
                                  genome_fasta: str,
                                  n_processes: int = None) -> pd.DataFrame:
    """Analyze methylation levels and SMARCB1 enrichment for gene lists."""
    
    if n_processes is None:
        n_processes = max(1, os.cpu_count() - 1)
    
    results = []
    
    # Get file paths
    ip_files = [os.path.join(medip_dir, f"Medip_PP_output_r{i}.bw") for i in range(1, 4)]
    input_files = [os.path.join(medip_dir, f"Medip_PP_input_r{i}.bw") for i in range(1, 4)]
    
    # SMARCB1 files with both normalizations
    bm_file_cpm = os.path.join(smarcb1_dir, "BM3_CPM.bw")
    bm_file_rpkm = os.path.join(smarcb1_dir, "BM3_RPKM.bw")
    bg_files_cpm = [os.path.join(smarcb1_dir, f"BG{i}_CPM.bw") for i in range(1, 3)]
    bg_files_rpkm = [os.path.join(smarcb1_dir, f"BG{i}_RPKM.bw") for i in range(1, 3)]
    
    with pysam.FastaFile(genome_fasta) as fasta:
        for category, df in gene_lists.items():
            print(f"Processing {category} genes...")
            
            for _, region in tqdm(df.iterrows(), total=len(df)):
                try:
                    # Calculate methylation
                    methylation_score, cpg_density = calculate_region_methylation(
                        region, ip_files, input_files, fasta
                    )
                    
                    # Calculate SMARCB1 enrichment with both normalizations
                    smarcb1_enrichment_cpm, smarcb1_enrichment_rpkm = calculate_smarcb1_enrichment(
                        region, bm_file_cpm, bm_file_rpkm, bg_files_cpm, bg_files_rpkm
                    )
                    
                    results.append({
                        'category': category,
                        'gene': region['Gene'],
                        'chromosome': region['seqnames'],
                        'start': region['start'],
                        'end': region['end'],
                        'methylation_score': methylation_score,
                        'cpg_density': cpg_density,
                        'smarcb1_enrichment_cpm': smarcb1_enrichment_cpm,
                        'smarcb1_enrichment_rpkm': smarcb1_enrichment_rpkm
                    })
                    
                except Exception as e:
                    print(f"Error processing region {region['seqnames']}:{region['start']}-{region['end']}: {str(e)}")
                    continue
    
    return pd.DataFrame(results)

def plot_results(results: pd.DataFrame, output_dir: str):
    """Create visualization plots for the analysis results."""
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Methylation distribution by category
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=results, x='category', y='methylation_score')
    plt.title('Methylation Score Distribution by Gene Category')
    plt.savefig(os.path.join(output_dir, 'methylation_distribution.png'))
    plt.close()
    
    # 2. SMARCB1 enrichment by category (CPM)
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=results, x='category', y='smarcb1_enrichment_cpm')
    plt.title('SMARCB1 Enrichment (CPM) by Gene Category')
    plt.savefig(os.path.join(output_dir, 'smarcb1_enrichment_cpm.png'))
    plt.close()
    
    # 3. SMARCB1 enrichment by category (RPKM)
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=results, x='category', y='smarcb1_enrichment_rpkm')
    plt.title('SMARCB1 Enrichment (RPKM) by Gene Category')
    plt.savefig(os.path.join(output_dir, 'smarcb1_enrichment_rpkm.png'))
    plt.close()
    
    # 4. Methylation vs SMARCB1 enrichment scatter plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    sns.scatterplot(data=results, x='methylation_score', y='smarcb1_enrichment_cpm', 
                    hue='category', alpha=0.6, ax=ax1)
    ax1.set_title('Methylation vs SMARCB1 (CPM)')
    
    sns.scatterplot(data=results, x='methylation_score', y='smarcb1_enrichment_rpkm', 
                    hue='category', alpha=0.6, ax=ax2)
    ax2.set_title('Methylation vs SMARCB1 (RPKM)')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'methylation_vs_smarcb1.png'))
    plt.close()
    
    # 5. CpG density vs Methylation
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=results, x='cpg_density', y='methylation_score',
                    hue='category', alpha=0.6)
    plt.title('CpG Density vs Methylation Score')
    plt.savefig(os.path.join(output_dir, 'cpg_density_vs_methylation.png'))
    plt.close()
    
    # 6. Correlation between CPM and RPKM enrichment
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=results, x='smarcb1_enrichment_cpm', y='smarcb1_enrichment_rpkm',
                    hue='category', alpha=0.6)
    plt.title('SMARCB1 Enrichment: CPM vs RPKM')
    plt.savefig(os.path.join(output_dir, 'cpm_vs_rpkm.png'))
    plt.close()

def main():
    # Set paths
    base_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Cross_final/data"
    medip_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/DATA/MECP2/MEDIP/output_done/bigwig"
    smarcb1_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1/results/bigwig"
    genome_fasta = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA/mm10.fa"
    output_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Cross_final/results/methylation_analysis"
    
    # Load gene lists
    gene_lists = load_gene_lists(base_dir)
    
    # Run analysis
    results = analyze_methylation_and_smarcb1(
        gene_lists, medip_dir, smarcb1_dir, genome_fasta
    )
    
    # Save results
    results.to_csv(os.path.join(output_dir, 'methylation_smarcb1_results.csv'), index=False)
    
    # Create plots
    plot_results(results, output_dir)

if __name__ == "__main__":
    main()
