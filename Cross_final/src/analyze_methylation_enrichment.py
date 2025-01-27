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
from scipy import stats

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
                               bm_file: str,
                               bg_files: List[str],
                               window_size: int = 500) -> float:
    """Calculate SMARCB1 enrichment (BM vs BG) for a region."""
    
    chrom = region['seqnames']
    center = (region['start'] + region['end']) // 2
    start = max(0, center - window_size)
    end = center + window_size
    
    # Get BM signal
    bm_signal = get_region_signal(bm_file, chrom, start, end)
    
    # Get average BG signals
    bg_signals = []
    for bg_file in bg_files:
        bg_signal = get_region_signal(bg_file, chrom, start, end)
        if len(bg_signal) > 0:
            bg_signals.append(bg_signal)
    
    if not bg_signals or len(bm_signal) == 0:
        return 0.0
    
    # Calculate enrichment
    mean_bg = np.mean(bg_signals, axis=0)
    bg_mean = np.mean(mean_bg)
    enrichment = np.mean(bm_signal) / bg_mean if bg_mean > 0 else 0.0
    
    return enrichment

def remove_outliers(df: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    """Remove outliers using IQR method for specified columns."""
    df_clean = df.copy()
    for col in columns:
        Q1 = df_clean[col].quantile(0.25)
        Q3 = df_clean[col].quantile(0.75)
        IQR = Q3 - Q1
        df_clean = df_clean[
            (df_clean[col] >= Q1 - 1.5 * IQR) & 
            (df_clean[col] <= Q3 + 1.5 * IQR)
        ]
    return df_clean

def calculate_statistics(data: pd.DataFrame, group_col: str, value_col: str) -> pd.DataFrame:
    """Calculate summary statistics for each group."""
    stats_df = data.groupby(group_col)[value_col].agg([
        'count',
        'mean',
        'std',
        'median',
        lambda x: x.quantile(0.25),
        lambda x: x.quantile(0.75)
    ]).round(3)
    stats_df.columns = ['count', 'mean', 'std', 'median', 'q25', 'q75']
    return stats_df

def run_statistical_tests(data: pd.DataFrame, group_col: str, value_col: str) -> pd.DataFrame:
    """Run statistical tests between groups."""
    categories = data[group_col].unique()
    test_results = []
    
    for i, cat1 in enumerate(categories):
        for cat2 in categories[i+1:]:
            group1 = data[data[group_col] == cat1][value_col]
            group2 = data[data[group_col] == cat2][value_col]
            
            # Mann-Whitney U test
            stat, pval = stats.mannwhitneyu(group1, group2, alternative='two-sided')
            
            test_results.append({
                'group1': cat1,
                'group2': cat2,
                'statistic': stat,
                'pvalue': pval
            })
    
    return pd.DataFrame(test_results)

def analyze_methylation_and_smarcb1(gene_lists: Dict[str, pd.DataFrame],
                                  medip_dir: str,
                                  smarcb1_dir: str,
                                  genome_fasta: str,
                                  output_dir: str,
                                  n_processes: int = None) -> pd.DataFrame:
    """Analyze methylation levels and SMARCB1 enrichment for MeCP2-bound regions."""
    
    if n_processes is None:
        n_processes = max(1, os.cpu_count() - 1)
    
    results = []
    
    # Get file paths
    ip_files = [os.path.join(medip_dir, f"Medip_PP_output_r{i}.bw") for i in range(1, 4)]
    input_files = [os.path.join(medip_dir, f"Medip_PP_input_r{i}.bw") for i in range(1, 4)]
    
    # SMARCB1 files (RPKM only)
    bm_file = os.path.join(smarcb1_dir, "BM3_RPKM.bw")  # Exogenous MeCP2
    bg_files = [os.path.join(smarcb1_dir, f"BG{i}_RPKM.bw") for i in range(1, 4)]  # Background
    
    with pysam.FastaFile(genome_fasta) as fasta:
        for category, df in gene_lists.items():
            print(f"Processing {category} genes...")
            
            for _, region in tqdm(df.iterrows(), total=len(df)):
                try:
                    # Calculate methylation
                    methylation_score, cpg_density = calculate_region_methylation(
                        region, ip_files, input_files, fasta
                    )
                    
                    # Calculate SMARCB1 enrichment
                    smarcb1_enrichment = calculate_smarcb1_enrichment(
                        region, bm_file, bg_files
                    )
                    
                    results.append({
                        'category': category,
                        'gene': region['Gene'],
                        'chromosome': region['seqnames'],
                        'start': region['start'],
                        'end': region['end'],
                        'methylation_score': methylation_score,
                        'cpg_density': cpg_density,
                        'smarcb1_enrichment': smarcb1_enrichment
                    })
                    
                except Exception as e:
                    print(f"Error processing region {region['seqnames']}:{region['start']}-{region['end']}: {str(e)}")
                    continue
    
    results_df = pd.DataFrame(results)
    
    # Save intermediate results
    results_df.to_csv(os.path.join(output_dir, 'all_results.csv'), index=False)
    
    # Calculate and save statistics
    stats_dir = os.path.join(output_dir, 'statistics')
    os.makedirs(stats_dir, exist_ok=True)
    
    # Methylation statistics
    meth_stats = calculate_statistics(results_df, 'category', 'methylation_score')
    meth_stats.to_csv(os.path.join(stats_dir, 'methylation_statistics.csv'))
    
    meth_tests = run_statistical_tests(results_df, 'category', 'methylation_score')
    meth_tests.to_csv(os.path.join(stats_dir, 'methylation_statistical_tests.csv'))
    
    # SMARCB1 statistics
    smarcb1_stats = calculate_statistics(results_df, 'category', 'smarcb1_enrichment')
    smarcb1_stats.to_csv(os.path.join(stats_dir, 'smarcb1_statistics.csv'))
    
    smarcb1_tests = run_statistical_tests(results_df, 'category', 'smarcb1_enrichment')
    smarcb1_tests.to_csv(os.path.join(stats_dir, 'smarcb1_statistical_tests.csv'))
    
    return results_df

def plot_results(results: pd.DataFrame, output_dir: str):
    """Create visualization plots with and without outliers."""
    
    plot_dir = os.path.join(output_dir, 'plots')
    os.makedirs(plot_dir, exist_ok=True)
    
    # Define columns for outlier removal
    outlier_cols = ['methylation_score', 'smarcb1_enrichment']
    results_no_outliers = remove_outliers(results, outlier_cols)
    
    # Function to create both versions of each plot
    def create_dual_plots(plot_func, filename_prefix, **kwargs):
        # With outliers
        plot_func(results, os.path.join(plot_dir, f'{filename_prefix}_with_outliers.png'), 
                 title_suffix='(with outliers)', **kwargs)
        # Without outliers
        plot_func(results_no_outliers, os.path.join(plot_dir, f'{filename_prefix}_no_outliers.png'),
                 title_suffix='(without outliers)', **kwargs)
    
    def plot_boxplot(data, output_file, y_col, title_prefix, title_suffix=''):
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=data, x='category', y=y_col)
        plt.title(f'{title_prefix} {title_suffix}')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()
    
    def plot_scatter(data, output_file, title_suffix=''):
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=data, x='methylation_score', y='smarcb1_enrichment',
                       hue='category', alpha=0.6)
        plt.title(f'Methylation vs SMARCB1 Enrichment\nin MeCP2-bound Regions {title_suffix}')
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()
    
    def plot_violin(data, output_file, y_col, title_prefix, title_suffix=''):
        plt.figure(figsize=(10, 6))
        sns.violinplot(data=data, x='category', y=y_col)
        plt.title(f'{title_prefix} {title_suffix}')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()
    
    # Create all plot versions
    create_dual_plots(
        lambda d, f, **k: plot_boxplot(d, f, 'methylation_score', 'Methylation Score Distribution', **k),
        'methylation_boxplot'
    )
    
    create_dual_plots(
        lambda d, f, **k: plot_boxplot(d, f, 'smarcb1_enrichment', 'SMARCB1 Enrichment Distribution', **k),
        'smarcb1_boxplot'
    )
    
    create_dual_plots(plot_scatter, 'methylation_vs_smarcb1')
    
    create_dual_plots(
        lambda d, f, **k: plot_violin(d, f, 'methylation_score', 'Methylation Score Distribution', **k),
        'methylation_violin'
    )
    
    create_dual_plots(
        lambda d, f, **k: plot_violin(d, f, 'smarcb1_enrichment', 'SMARCB1 Enrichment Distribution', **k),
        'smarcb1_violin'
    )
    
    # Save correlation analysis
    correlations = results.groupby('category').apply(
        lambda x: pd.Series({
            'correlation': stats.spearmanr(x['methylation_score'], x['smarcb1_enrichment'])[0],
            'pvalue': stats.spearmanr(x['methylation_score'], x['smarcb1_enrichment'])[1]
        })
    )
    correlations.to_csv(os.path.join(output_dir, 'statistics', 'methylation_smarcb1_correlations.csv'))

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
        gene_lists, medip_dir, smarcb1_dir, genome_fasta, output_dir
    )
    
    # Create plots
    plot_results(results, output_dir)

if __name__ == "__main__":
    main()
