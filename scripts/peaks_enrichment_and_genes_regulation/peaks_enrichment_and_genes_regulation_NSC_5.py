import sys
import subprocess

def check_numpy_version():
    """Check numpy version and downgrade if needed"""
    try:
        import numpy as np
        numpy_version = np.__version__
        if numpy_version.startswith('2.'):
            print("Detected NumPy 2.x - downgrading to 1.x for compatibility")
            subprocess.check_call([sys.executable, "-m", "pip", "install", "numpy<2"])
            print("NumPy downgraded successfully. Please restart the script.")
            sys.exit(0)
    except ImportError:
        print("NumPy not found. Installing compatible version...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "numpy<2"])
        print("NumPy installed successfully. Please restart the script.")
        sys.exit(0)

# Add this at the start of the script
# check_numpy_version()

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools
import os
import pysam
import time
from pybedtools import BedTool
import multiprocessing as mp
from functools import partial
from itertools import chain
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from functools import partial
import numpy as np
from tqdm.auto import tqdm

PROMOTER_WINDOW = 2000  # Define promoter region as ±2kb from TSS

# Set working directory
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Analyze enrichment for NSC samples')
parser.add_argument('--working-dir', type=str, required=True,
                   help='Path to working directory')
parser.add_argument('--data-dir', type=str, required=True,
                   help='Path to data directory')
parser.add_argument('--results-dir', type=str, required=True,
                   help='Path to results directory')
args = parser.parse_args()

os.chdir(args.working_dir)

DATA_DIR = args.data_dir
RESULTS_DIR = args.results_dir

# Add function to calculate sequencing depth
def calculate_sequencing_depth(bam_file):
    """Calculate total mapped reads from BAM file"""
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        return bam.count()

def normalize_signals(df, depth):
    """Apply multiple normalization methods to signals"""
    if df.empty:
        return df
    
    df = df.copy()
    
    # Basic depth normalization
    df['signalValue'] = df['signalValue'].clip(lower=0) * (1e6 / depth)
    
    # Additional normalization methods
    # 1. Quantile normalization
    df['signalValue'] = stats.rankdata(df['signalValue']) / len(df['signalValue'])
    
    # 2. Z-score normalization
    df['signalValue'] = (df['signalValue'] - df['signalValue'].mean()) / df['signalValue'].std()
    
    # 3. Min-max normalization
    df['signalValue'] = (df['signalValue'] - df['signalValue'].min()) / \
                        (df['signalValue'].max() - df['signalValue'].min())
    
    return df

def load_data():
    try:
        # Load DEA results
        dea_nsc = pd.read_csv("../DATA/DEA_NSC.csv")
        print("\nDEA file columns:", dea_nsc.columns.tolist())

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

        def calculate_sequencing_depth_safe(bam_file):
            """Calculate sequencing depth with error handling"""
            try:
                if not os.path.exists(bam_file):
                    print(f"Warning: BAM file not found: {bam_file}")
                    return 1
                
                with pysam.AlignmentFile(bam_file, "rb") as bam:
                    if not bam.check_index():
                        print(f"Warning: BAM index missing for {bam_file}")
                        return 1
                    return max(bam.count(), 1)  # Ensure non-zero return
            except Exception as e:
                print(f"Error calculating depth for {bam_file}: {str(e)}")
                return 1

        # Load samples with enhanced error handling
        exo_samples = ['NSCv1', 'NSCv2', 'NSCv3']
        endo_samples = ['NSCM1', 'NSCM2', 'NSCM3']
        
        peaks_exo = {}
        peaks_endo = {}
        depths_exo = {}
        depths_endo = {}
        
        # Update BAM file path to use ALIGN_DIR
        ALIGN_DIR = f"{args.working_dir}/results_1b/aligned"
        
        # Load exogenous samples
        for sample in exo_samples:
            peak_file = f"{RESULTS_DIR}/peaks/{sample}_peaks.narrowPeak"
            peaks_exo[sample] = load_peak_file(peak_file, sample)
            # Update BAM file path
            bam_file = f"{ALIGN_DIR}/{sample}.bam"
            depths_exo[sample] = calculate_sequencing_depth_safe(bam_file)
            
            if not peaks_exo[sample].empty:
                peaks_exo[sample] = normalize_signals(peaks_exo[sample], depths_exo[sample])
        
        # Load endogenous samples
        for sample in endo_samples:
            peak_file = f"{RESULTS_DIR}/peaks/{sample}_peaks.narrowPeak"
            peaks_endo[sample] = load_peak_file(peak_file, sample)
            # Update BAM file path
            bam_file = f"{ALIGN_DIR}/{sample}.bam"
            depths_endo[sample] = calculate_sequencing_depth_safe(bam_file)
            
            if not peaks_endo[sample].empty:
                peaks_endo[sample] = normalize_signals(peaks_endo[sample], depths_endo[sample])
        
        # Validate we have usable data
        if all(df.empty for df in peaks_exo.values()) or all(df.empty for df in peaks_endo.values()):
            print("Warning: No valid peak data found for analysis")
            raise ValueError("Insufficient peak data for analysis")

        # Load gene annotations
        gene_annotations, name_to_info = load_gene_annotations()
        
        return dea_nsc, peaks_exo, peaks_endo, gene_annotations, name_to_info
        
    except Exception as e:
        print(f"Error in load_data: {str(e)}")
        raise

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

def load_gene_annotations():
    """Load gene annotations from GTF file and extract gene names"""
    print("Loading gene annotations...")
    gtf_file = "../DATA/gencode.vM10.annotation.gtf"
    
    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file not found: {gtf_file}")
    
    gene_annotations = pd.read_csv(gtf_file, sep='\t', comment='#',
                                 names=['chr', 'source', 'feature', 'start', 'end',
                                       'score', 'strand', 'frame', 'attributes'])
    
    # Filter for genes only
    gene_annotations = gene_annotations[gene_annotations['feature'] == 'gene']
    
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
    
    # Extract gene information
    gene_info = gene_annotations['attributes'].apply(extract_gene_info)
    gene_annotations = pd.concat([gene_annotations, gene_info], axis=1)
    
    # Create standardized versions of gene names and IDs
    gene_annotations['gene_name_std'] = gene_annotations['gene_name'].apply(standardize_gene_name)
    gene_annotations['gene_id_std'] = gene_annotations['gene_id'].apply(standardize_gene_name)
    
    # Create a mapping dictionary for quick lookups
    name_to_info = {}
    for _, row in gene_annotations.iterrows():
        if not pd.isna(row['gene_name']):
            name_to_info[standardize_gene_name(row['gene_name'])] = row
        if not pd.isna(row['gene_id']):
            name_to_info[standardize_gene_name(row['gene_id'])] = row
    
    print(f"Loaded {len(gene_annotations)} genes from GTF")
    print(f"Created mapping for {len(name_to_info)} unique gene identifiers")
    
    return gene_annotations, name_to_info

def get_peaks_near_gene(gene, peaks_dict, gene_annotations, name_to_info, window=PROMOTER_WINDOW):
    """Get peaks near a gene's promoter region"""
    gene_std = standardize_gene_name(gene)
    
    if gene_std not in name_to_info:
        return pd.DataFrame()
    
    gene_info = name_to_info[gene_std]
    
    # Create promoter region (TSS ± window)
    if gene_info['strand'] == '+':
        promoter_start = max(0, gene_info['start'] - window)
        promoter_end = gene_info['start'] + window
    else:
        promoter_start = max(0, gene_info['end'] - window)
        promoter_end = gene_info['end'] + window
    
    # Collect overlapping peaks from all samples
    all_peaks = []
    for sample, peaks in peaks_dict.items():
        # Filter peaks for the same chromosome first
        chr_peaks = peaks[peaks['chr'] == gene_info['chr']]
        
        # Find overlapping peaks
        overlapping = chr_peaks.loc[
            (chr_peaks['start'] <= promoter_end) & 
            (chr_peaks['end'] >= promoter_start)
        ].copy()  # Create explicit copy
        
        if not overlapping.empty:
            overlapping.loc[:, 'sample'] = sample
            all_peaks.append(overlapping)
    
    if not all_peaks:
        return pd.DataFrame()
    
    combined_peaks = pd.concat(all_peaks, ignore_index=True)
    return combined_peaks

def calculate_total_genome_peaks(peaks_exo, peaks_endo):
    """Calculate total number of peaks across all samples"""
    total_exo = sum(len(df) for df in peaks_exo.values())
    total_endo = sum(len(df) for df in peaks_endo.values())
    return total_exo + total_endo

def define_enrichment_methods(total_genome_peaks):
    """Define different methods to calculate enrichment"""
    
    methods = {
        # Method 1: Normalized Signal Value Ratio
        'signal_ratio': lambda exo, endo: (
            np.mean(exo['signalValue']) / np.mean(endo['signalValue'])
            if len(exo) > 0 and len(endo) > 0
            else float('inf') if len(exo) > 0
            else 0.0
        ),
        
        # Method 2: Peak Count Ratio
        'peak_count': lambda exo, endo: (
            len(exo.groupby(['chr', 'start', 'end'])) / 
            max(len(endo.groupby(['chr', 'start', 'end'])), 1)  # Prevent division by zero
        ),
        
        # Method 3: Combined Normalized Score
        'combined_score': lambda exo, endo: (
            (np.mean(exo['signalValue']) * len(exo.groupby(['chr', 'start', 'end']))) / 
            max((np.mean(endo['signalValue']) * len(endo.groupby(['chr', 'start', 'end']))), 1)
            if len(exo) > 0
            else 0.0
        ),
        
        # Method 4: Statistical Enrichment
        'statistical': lambda exo, endo: (
            stats.fisher_exact([
                [len(exo.groupby(['chr', 'start', 'end'])), 
                 len(endo.groupby(['chr', 'start', 'end']))],
                [total_genome_peaks - len(exo.groupby(['chr', 'start', 'end'])), 
                 total_genome_peaks - len(endo.groupby(['chr', 'start', 'end']))]
            ])[1] if len(exo) > 0 or len(endo) > 0 else np.nan
        ),
        
        # Method 5: Width-Weighted Signal Ratio
        'width_weighted': lambda exo, endo: (
            np.sum(exo['signalValue'] * (exo['end'] - exo['start'])) / 
            max(np.sum(endo['signalValue'] * (endo['end'] - endo['start'])), 1)
            if len(exo) > 0
            else 0.0
        ),
        
        # Method 6: Width-Normalized Coverage Score
        'coverage_score': lambda exo, endo: (
            (np.sum(exo['end'] - exo['start']) * np.mean(exo['signalValue'])) /
            max((np.sum(endo['end'] - endo['start']) * np.mean(endo['signalValue'])), 1)
            if len(exo) > 0
            else 0.0
        ),
        
        # Method 7: Peak Area Integration
        'area_integration': lambda exo, endo: (
            np.sum(
                (exo['end'] - exo['start']) * 
                exo['signalValue'] * 
                np.clip(-np.log10(exo['qValue']), 0, 50)  # First take -log10, then clip
            ) / max(np.sum(
                (endo['end'] - endo['start']) * 
                endo['signalValue'] * 
                np.clip(-np.log10(endo['qValue']), 0, 50)
            ), 1)
            if len(exo) > 0
            else 0.0
        )
    }
    
    return methods

def analyze_enrichment(dea, peaks_exo, peaks_endo, gene_annotations, name_to_info):
    # Load data
    # print("Loading data...")
    # dea, peaks_exo, peaks_endo, gene_annotations, name_to_info = load_data()
    
    # Print diagnostic information
    print(f"\nDiagnostic information:")
    print(f"DEA genes: {len(dea)}")
    print(f"Annotation genes: {len(gene_annotations)}")
    
    # Standardize DEA gene names
    dea['gene_std'] = dea['gene'].apply(standardize_gene_name)
    
    # Define up-regulated genes (log2FC > 1 and padj < 0.05)
    upreg_genes = dea[
        (dea['log2FoldChange'] > 0.5) & 
        (dea['padj'] < 0.05) &
        (dea['gene_std'].notna())
    ]['gene_std'].tolist()
    
    print(f"Up-regulated genes: {len(upreg_genes)}")
    
    # Calculate total genome peaks
    total_genome_peaks = calculate_total_genome_peaks(peaks_exo, peaks_endo)
    print(f"Total genome peaks: {total_genome_peaks}")
    
    # Calculate enrichment
    methods = define_enrichment_methods(total_genome_peaks)
    results = {}
    
    for method_name, method_func in methods.items():
        print(f"\nProcessing method: {method_name}")
        enrichment_scores = []
        peaks_found = 0
        
        for i, gene in enumerate(upreg_genes, 1):
            if i % 100 == 0:
                print(f"Processing gene {i}/{len(upreg_genes)}")
            
            try:
                exo_peaks = get_peaks_near_gene(gene, peaks_exo, gene_annotations, name_to_info)
                endo_peaks = get_peaks_near_gene(gene, peaks_endo, gene_annotations, name_to_info)
                
                if not exo_peaks.empty or not endo_peaks.empty:
                    peaks_found += 1
                    score = method_func(exo_peaks, endo_peaks)
                    
                    # Get DEA information for this gene
                    dea_info = dea[dea['gene_std'] == gene].iloc[0]
                    
                    if not np.isnan(score):
                        # Create peak coordinate strings
                        exo_coords = []
                        if not exo_peaks.empty:
                            exo_coords = [f"{row['chr']}:{row['start']}-{row['end']}" 
                                        for _, row in exo_peaks.iterrows()]
                        
                        endo_coords = []
                        if not endo_peaks.empty:
                            endo_coords = [f"{row['chr']}:{row['start']}-{row['end']}" 
                                         for _, row in endo_peaks.iterrows()]
                        
                        enrichment_scores.append({
                            'gene': gene,
                            'enrichment_score': score,
                            'log2FoldChange': dea_info['log2FoldChange'],
                            'padj': dea_info['padj'],
                            'exo_peaks': ';'.join(exo_coords) if exo_coords else 'None',
                            'endo_peaks': ';'.join(endo_coords) if endo_coords else 'None',
                            'num_exo_peaks': len(exo_coords),
                            'num_endo_peaks': len(endo_coords)
                        })
                
            except Exception as e:
                print(f"Error processing gene {gene}: {str(e)}")
                continue
        
        print(f"Found peaks for {peaks_found} genes")
        print(f"Found enrichment scores for {len(enrichment_scores)} genes")
        
        if enrichment_scores:
            results[method_name] = pd.DataFrame(enrichment_scores)
    
    # Save detailed results to CSV
    for method_name, df in results.items():
        # Sort by enrichment score in descending order
        df_sorted = df.sort_values('enrichment_score', ascending=False)
        df_sorted.to_csv(f'{RESULTS_DIR}/enrichment_{method_name}_NSC.csv', index=False)
    
    return results

def plot_enrichment(results):
    """Plot enrichment results using different methods"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 15))
    
    for (method_name, df), ax in zip(results.items(), axes.flat):
        sns.histplot(data=df, x='enrichment_score', ax=ax)
        ax.set_title(f'Enrichment Distribution ({method_name})')
        ax.set_xlabel('Enrichment Score')
        ax.set_ylabel('Count')
        
        # Add median line
        median = df['enrichment_score'].median()
        ax.axvline(median, color='red', linestyle='--', 
                  label=f'Median: {median:.2f}')
        ax.legend()
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/enrichment_analysis_NSC.pdf')
    plt.close()

def summarize_results(results):
    """Create summary statistics for each enrichment method"""
    summary = {}
    for method_name, df in results.items():
        summary[method_name] = {
            'median_enrichment': df['enrichment_score'].median(),
            'mean_enrichment': df['enrichment_score'].mean(),
            'std_enrichment': df['enrichment_score'].std(),
            'num_genes': len(df),
            'significant_genes': len(df[df['enrichment_score'] > 1]),
            'median_log2FC': df['log2FoldChange'].median(),
            'mean_log2FC': df['log2FoldChange'].mean(),
            'median_padj': df['padj'].median(),
            'highly_significant': len(df[
                (df['enrichment_score'] > 1) & 
                (df['log2FoldChange'] > 1) & 
                (df['padj'] < 0.01)
            ])
        }
    
    summary_df = pd.DataFrame(summary).T
    summary_df.to_csv(f'{RESULTS_DIR}/enrichment_summary_NSC.csv')
    return summary_df

def print_gene_name_examples(dea, name_to_info):
    """Print examples of gene name matching"""
    print("\nGene name matching examples:")
    for gene in dea['gene'].head(10):
        std_name = standardize_gene_name(gene)
        found = std_name in name_to_info
        print(f"Original: {gene:20} Standardized: {std_name:20} Found: {found}")

def plot_peak_width_distributions(peaks_exo, peaks_endo):
    """Plot histograms of peak widths for exogenous and endogenous samples"""
    plt.figure(figsize=(12, 6))
    
    # Calculate peak widths
    exo_widths = []
    endo_widths = []
    
    # Collect all peak widths
    for sample, peaks in peaks_exo.items():
        widths = peaks['end'] - peaks['start']
        exo_widths.extend(widths)
    
    for sample, peaks in peaks_endo.items():
        widths = peaks['end'] - peaks['start']
        endo_widths.extend(widths)
    
    # Create histogram
    plt.hist(exo_widths, bins=50, alpha=0.5, label='Exogenous', density=True)
    plt.hist(endo_widths, bins=50, alpha=0.5, label='Endogenous', density=True)
    
    # Add labels and title
    plt.xlabel('Peak Width (bp)')
    plt.ylabel('Density')
    plt.title('Distribution of Peak Widths')
    plt.legend()
    
    # Add summary statistics as text
    exo_stats = f'Exo - Mean: {np.mean(exo_widths):.0f}, Median: {np.median(exo_widths):.0f}'
    endo_stats = f'Endo - Mean: {np.mean(endo_widths):.0f}, Median: {np.median(endo_widths):.0f}'
    plt.text(0.02, 0.98, exo_stats + '\n' + endo_stats,
             transform=plt.gca().transAxes, 
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Save plot
    plt.savefig(f'{RESULTS_DIR}/peak_width_distributions.pdf')
    plt.close()

def plot_detailed_peak_width_distributions(peaks_exo, peaks_endo):
    """Plot detailed histograms of peak widths for each sample"""
    # Calculate number of samples
    n_exo = len(peaks_exo)
    n_endo = len(peaks_endo)
    total_samples = n_exo + n_endo
    
    # Create subplot grid
    fig, axes = plt.subplots(total_samples, 1, figsize=(12, 4*total_samples))
    
    # Plot exogenous samples
    for i, (sample, peaks) in enumerate(peaks_exo.items()):
        widths = peaks['end'] - peaks['start']
        axes[i].hist(widths, bins=50, alpha=0.7)
        axes[i].set_title(f'{sample} Peak Width Distribution')
        axes[i].set_xlabel('Peak Width (bp)')
        axes[i].set_ylabel('Count')
        
        # Add statistics
        stats = f'Mean: {np.mean(widths):.0f}\nMedian: {np.median(widths):.0f}\nCount: {len(widths)}'
        axes[i].text(0.98, 0.98, stats,
                    transform=axes[i].transAxes,
                    verticalalignment='top',
                    horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Plot endogenous samples
    for i, (sample, peaks) in enumerate(peaks_endo.items()):
        widths = peaks['end'] - peaks['start']
        axes[i+n_exo].hist(widths, bins=50, alpha=0.7)
        axes[i+n_exo].set_title(f'{sample} Peak Width Distribution')
        axes[i+n_exo].set_xlabel('Peak Width (bp)')
        axes[i+n_exo].set_ylabel('Count')
        
        # Add statistics
        stats = f'Mean: {np.mean(widths):.0f}\nMedian: {np.median(widths):.0f}\nCount: {len(widths)}'
        axes[i+n_exo].text(0.98, 0.98, stats,
                          transform=axes[i+n_exo].transAxes,
                          verticalalignment='top',
                          horizontalalignment='right',
                          bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/peak_width_distributions_detailed.pdf')
    plt.close()

def plot_width_vs_enrichment(results, peaks_exo, peaks_endo, gene_annotations, name_to_info):
    """Plot relationship between peak widths and enrichment scores"""
    plt.figure(figsize=(15, 5))
    
    # Calculate mean peak widths for each gene
    width_ratios = {}
    for gene in results['width_weighted']['gene']:
        exo_peaks = get_peaks_near_gene(gene, peaks_exo, gene_annotations, name_to_info)
        endo_peaks = get_peaks_near_gene(gene, peaks_endo, gene_annotations, name_to_info)
        
        if not exo_peaks.empty and not endo_peaks.empty:
            mean_width_exo = (exo_peaks['end'] - exo_peaks['start']).mean()
            mean_width_endo = (endo_peaks['end'] - endo_peaks['start']).mean()
            if mean_width_endo > 0:  # Prevent division by zero
                width_ratios[gene] = mean_width_exo / mean_width_endo
    
    # Create scatter plots
    for i, method in enumerate(['width_weighted', 'coverage_score', 'area_integration']):
        plt.subplot(1, 3, i+1)
        
        # Get data points and filter out invalid values
        data_points = [(width_ratios[gene], score) 
                      for gene, score in zip(results[method]['gene'], results[method]['enrichment_score'])
                      if gene in width_ratios and score > 0 and width_ratios[gene] > 0]
        
        if data_points:  # Only plot if we have valid data points
            x, y = zip(*data_points)
            
            plt.scatter(x, y, alpha=0.5)
            plt.xlabel('Peak Width Ratio (Exo/Endo)')
            plt.ylabel(f'{method} Enrichment Score')
            plt.title(f'{method} vs Width Ratio\n(n={len(x)} genes)')
            
            try:
                plt.yscale('log')
                plt.xscale('log')
            except ValueError:
                # Fallback to linear scale if log scale fails
                plt.yscale('linear')
                plt.xscale('linear')
                print(f"Warning: Could not use log scale for {method}, using linear scale instead")
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/width_vs_enrichment_NSC.pdf')
    plt.close()

def summarize_peak_distribution(results):
    """Summarize the distribution of peaks across conditions"""
    summary = {}
    for method, df in results.items():
        summary[method] = {
            'total_genes': len(df),
            'exo_only': len(df[df['enrichment_score'] == float('inf')]),
            'endo_only': len(df[df['enrichment_score'] == 0]),
            'both': len(df[(df['enrichment_score'] != float('inf')) & 
                          (df['enrichment_score'] != 0) & 
                          (df['enrichment_score'].notna())])
        }
    
    summary_df = pd.DataFrame(summary).T
    summary_df.to_csv(f'{RESULTS_DIR}/peak_distribution_summary_NSC.csv')
    return summary_df

def load_cpg_islands():
    """Load CpG islands from bed file"""
    try:
        cpg_file = f"../DATA/cpg_islands.bed"
        print(f"\nLoading CpG islands from {cpg_file}")
        
        cpg_islands = pd.read_csv(cpg_file, sep='\t', 
                                 names=['chr', 'start', 'end', 'name', 'score', 'strand'])
        
        print(f"Loaded {len(cpg_islands)} CpG islands")
        return cpg_islands
    except Exception as e:
        print(f"Error loading CpG islands: {str(e)}")
        raise

def analyze_cpg_enrichment(peaks_exo, peaks_endo, cpg_islands):
    """Analyze enrichment of peaks in CpG islands"""
    print("\nAnalyzing CpG island enrichment...")
    
    # Convert dataframes to BedTool objects
    cpg_bed = BedTool.from_dataframe(cpg_islands[['chr', 'start', 'end']])
    
    results = {}
    for sample, peaks in peaks_exo.items():
        peaks_bed = BedTool.from_dataframe(peaks[['chr', 'start', 'end', 'signalValue']])
        intersect = peaks_bed.intersect(cpg_bed, wa=True, wb=True)
        results[f"{sample}_cpg"] = pd.DataFrame(intersect)
    
    for sample, peaks in peaks_endo.items():
        peaks_bed = BedTool.from_dataframe(peaks[['chr', 'start', 'end', 'signalValue']])
        intersect = peaks_bed.intersect(cpg_bed, wa=True, wb=True)
        results[f"{sample}_cpg"] = pd.DataFrame(intersect)
    
    # Calculate enrichment metrics for each CpG island
    cpg_enrichment = []
    
    for _, cpg in cpg_islands.iterrows():
        exo_signal = 0
        endo_signal = 0
        
        # Sum signals for exo samples
        for sample in peaks_exo.keys():
            if f"{sample}_cpg" in results:
                sample_peaks = results[f"{sample}_cpg"]
                mask = (sample_peaks[7] == cpg['chr']) & \
                       (sample_peaks[8] == cpg['start']) & \
                       (sample_peaks[9] == cpg['end'])
                if not sample_peaks[mask].empty:
                    exo_signal += sample_peaks[mask][3].sum()  # signalValue
        
        # Sum signals for endo samples
        for sample in peaks_endo.keys():
            if f"{sample}_cpg" in results:
                sample_peaks = results[f"{sample}_cpg"]
                mask = (sample_peaks[7] == cpg['chr']) & \
                       (sample_peaks[8] == cpg['start']) & \
                       (sample_peaks[9] == cpg['end'])
                if not sample_peaks[mask].empty:
                    endo_signal += sample_peaks[mask][3].sum()  # signalValue
        
        # Calculate enrichment
        enrichment = exo_signal / max(endo_signal, 1)  # Prevent division by zero
        
        cpg_enrichment.append({
            'chr': cpg['chr'],
            'start': cpg['start'],
            'end': cpg['end'],
            'exo_signal': exo_signal,
            'endo_signal': endo_signal,
            'enrichment': enrichment
        })
    
    enrichment_df = pd.DataFrame(cpg_enrichment)
    enrichment_df.to_csv(f'{RESULTS_DIR}/cpg_enrichment_NSC.csv', index=False)
    
    return enrichment_df

def plot_cpg_enrichment(enrichment_df):
    """Create visualizations for CpG island enrichment"""
    # Distribution of enrichment scores
    plt.figure(figsize=(10, 6))
    sns.histplot(data=enrichment_df, x='enrichment', bins=50)
    plt.title('Distribution of CpG Island Enrichment Scores')
    plt.xlabel('Enrichment Score (Exo/Endo)')
    plt.ylabel('Count')
    plt.savefig(f'{RESULTS_DIR}/cpg_enrichment_distribution.pdf')
    plt.close()
    
    # Scatter plot of exo vs endo signals
    plt.figure(figsize=(10, 6))
    plt.scatter(enrichment_df['endo_signal'], enrichment_df['exo_signal'], alpha=0.5)
    plt.xlabel('Endogenous Signal')
    plt.ylabel('Exogenous Signal')
    plt.title('Exogenous vs Endogenous Signal in CpG Islands')
    plt.savefig(f'{RESULTS_DIR}/cpg_signal_comparison.pdf')
    plt.close()

def integrate_with_rna_seq(enrichment_df, dea_nsc, gene_annotations):
    """Integrate CpG enrichment with RNA-seq data"""
    print("\nIntegrating CpG enrichment with RNA-seq data...")
    
    # Create BedTool objects
    cpg_bed = BedTool.from_dataframe(
        enrichment_df[['chr', 'start', 'end']]
    )
    genes_bed = BedTool.from_dataframe(
        gene_annotations[['chr', 'start', 'end', 'gene_name']]
    )
    
    # Find overlaps between CpG islands and genes
    overlaps = cpg_bed.intersect(genes_bed, wa=True, wb=True)
    
    # Create integrated dataset
    integrated_data = []
    for overlap in overlaps:
        cpg_idx = enrichment_df[
            (enrichment_df['chr'] == overlap[0]) & 
            (enrichment_df['start'] == int(overlap[1])) & 
            (enrichment_df['end'] == int(overlap[2]))
        ].index[0]
        
        gene_name = overlap[7]
        if gene_name in dea_nsc['gene'].values:
            dea_info = dea_nsc[dea_nsc['gene'] == gene_name].iloc[0]
            
            integrated_data.append({
                'chr': overlap[0],
                'cpg_start': int(overlap[1]),
                'cpg_end': int(overlap[2]),
                'gene': gene_name,
                'enrichment': enrichment_df.loc[cpg_idx, 'enrichment'],
                'log2FoldChange': dea_info['log2FoldChange'],
                'padj': dea_info['padj']
            })
    
    integrated_df = pd.DataFrame(integrated_data)
    integrated_df.to_csv(f'{RESULTS_DIR}/cpg_rna_integrated_NSC.csv', index=False)
    
    return integrated_df

def calculate_enrichment_scores(exo_row, endo_row):
    """Calculate multiple enrichment scores"""
    
    # Prevent division by zero
    min_signal = 1e-10
    
    # Basic signal ratio
    basic_enrichment = exo_row['signal_basic'] / max(endo_row['signal_basic'], min_signal)
    
    # Width-weighted signal ratio (Method 5)
    ww_enrichment = exo_row['signal_ww'] / max(endo_row['signal_ww'], min_signal)
    
    # Coverage score ratio (Method 6)
    cs_enrichment = exo_row['signal_cs'] / max(endo_row['signal_cs'], min_signal)
    
    return {
        'exo_signal': exo_row['signal_basic'],
        'endo_signal': endo_row['signal_basic'],
        'enrichment': basic_enrichment,
        'exo_ww': exo_row['signal_ww'],
        'endo_ww': endo_row['signal_ww'],
        'enrichment_ww': ww_enrichment,
        'exo_cs': exo_row['signal_cs'],
        'endo_cs': endo_row['signal_cs'],
        'enrichment_cs': cs_enrichment
    }

def process_overlap_chunk(chunk, exo_combined, endo_combined):
    """Process a chunk of overlaps"""
    chunk_results = []
    for overlap in chunk:
        exo_mask = ((exo_combined['chr'] == overlap[0]) & 
                    (exo_combined['start'] == int(overlap[1])) & 
                    (exo_combined['end'] == int(overlap[2])))
        
        endo_mask = ((endo_combined['chr'] == overlap[3]) & 
                     (endo_combined['start'] == int(overlap[4])) & 
                     (endo_combined['end'] == int(overlap[5])))
        
        if not any(exo_mask) or not any(endo_mask):
            continue
        
        exo_row = exo_combined[exo_mask].iloc[0]
        endo_row = endo_combined[endo_mask].iloc[0]
        
        enrichment_scores = calculate_enrichment_scores(exo_row, endo_row)
        
        chunk_results.append({
            'chr': overlap[0],
            'start': int(overlap[1]),
            'end': int(overlap[2]),
            **enrichment_scores,
            'exo_qValue': exo_row['qValue'],
            'endo_qValue': endo_row['qValue'],
            'exo_replicates': exo_row['num_replicates'],
            'endo_replicates': endo_row['num_replicates']
        })
    return chunk_results

def parallel_process_overlaps(overlaps, exo_combined, endo_combined, chunk_size=1000):
    """Process overlaps in parallel chunks"""
    # Split overlaps into chunks
    chunks = [overlaps[i:i + chunk_size] for i in range(0, len(overlaps), chunk_size)]
    
    # Process chunks in parallel
    with ProcessPoolExecutor() as executor:
        process_chunk_partial = partial(process_overlap_chunk, 
                                     exo_combined=exo_combined, 
                                     endo_combined=endo_combined)
        results = list(tqdm(
            executor.map(process_chunk_partial, chunks),
            total=len(chunks),
            desc="Processing overlaps"
        ))
    
    # Flatten results
    return [item for sublist in results for item in sublist]

def process_gene_chunk(gene_chunk, gene_annotations, peaks_exo, peaks_endo):
    """Process a chunk of genes"""
    chunk_results = []
    for gene in gene_chunk:
        if gene not in gene_annotations['gene_name'].values:
            continue
        
        gene_info = gene_annotations[gene_annotations['gene_name'] == gene].iloc[0]
        gene_start = gene_info['start'] - 2000
        gene_end = gene_info['end'] + 2000
        
        # Vectorized peak filtering
        exo_peaks = pd.concat([
            peaks[
                (peaks['chr'] == gene_info['chr']) &
                (peaks['start'].between(gene_start, gene_end))
            ] for peaks in peaks_exo.values()
        ])
        
        endo_peaks = pd.concat([
            peaks[
                (peaks['chr'] == gene_info['chr']) &
                (peaks['start'].between(gene_start, gene_end))
            ] for peaks in peaks_endo.values()
        ])
        
        if len(exo_peaks) > 0 or len(endo_peaks) > 0:
            chunk_results.append({
                'gene': gene,
                'exo_peaks': len(exo_peaks),
                'endo_peaks': len(endo_peaks),
                'exo_signal': exo_peaks['signalValue'].sum() if len(exo_peaks) > 0 else 0,
                'endo_signal': endo_peaks['signalValue'].sum() if len(endo_peaks) > 0 else 0,
                'enrichment': (exo_peaks['signalValue'].sum() / max(endo_peaks['signalValue'].sum(), 1))
                            if len(exo_peaks) > 0 else 0
            })
    return chunk_results

def parallel_analyze_categories(genes, gene_annotations, peaks_exo, peaks_endo, chunk_size=100):
    """Analyze gene categories in parallel"""
    # Split genes into chunks
    gene_chunks = [genes[i:i + chunk_size] for i in range(0, len(genes), chunk_size)]
    
    # Process chunks in parallel
    with ProcessPoolExecutor() as executor:
        process_chunk_partial = partial(process_gene_chunk,
                                     gene_annotations=gene_annotations,
                                     peaks_exo=peaks_exo,
                                     peaks_endo=peaks_endo)
        results = list(tqdm(
            executor.map(process_chunk_partial, gene_chunks),
            total=len(gene_chunks),
            desc="Processing genes"
        ))
    
    # Flatten results
    return [item for sublist in results for item in sublist]

def analyze_mecp2_enrichment_independent(peaks_exo, peaks_endo, n_jobs=None, results_dir=None):
    """Optimized version with parallel processing and intermediate file saving"""
    print("\nAnalyzing Mecp2 enrichment independently...")
    
    # Check for existing intermediate files
    exo_combined_file = f'{results_dir}/intermediate/exo_combined_peaks.csv'
    endo_combined_file = f'{results_dir}/intermediate/endo_combined_peaks.csv'
    enrichment_file = f'{results_dir}/intermediate/enrichment_results.csv'
    
    # Create intermediate directory if it doesn't exist
    os.makedirs(f'{results_dir}/intermediate', exist_ok=True)
    
    # Process exo peaks
    if os.path.exists(exo_combined_file):
        print("Loading existing combined exo peaks...")
        exo_combined = pd.read_csv(exo_combined_file)
    else:
        print("Combining exo peaks in parallel...")
        exo_combined = parallel_combine_peaks(peaks_exo, n_jobs=n_jobs)
        exo_combined.to_csv(exo_combined_file, index=False)
    
    # Process endo peaks
    if os.path.exists(endo_combined_file):
        print("Loading existing combined endo peaks...")
        endo_combined = pd.read_csv(endo_combined_file)
    else:
        print("Combining endo peaks in parallel...")
        endo_combined = parallel_combine_peaks(peaks_endo, n_jobs=n_jobs)
        endo_combined.to_csv(endo_combined_file, index=False)
    
    # Check for existing enrichment results
    if os.path.exists(enrichment_file):
        print("Loading existing enrichment results...")
        enrichment_df = pd.read_csv(enrichment_file)
    else:
        # Convert to BedTool objects
        print("Finding overlapping regions...")
        exo_bed = BedTool.from_dataframe(exo_combined[['chr', 'start', 'end']])
        endo_bed = BedTool.from_dataframe(endo_combined[['chr', 'start', 'end']])
        
        # Find overlapping regions
        overlaps = list(exo_bed.intersect(endo_bed, wa=True, wb=True))
        
        # Process overlaps in parallel
        print("Processing overlaps in parallel...")
        enrichment_data = parallel_process_overlaps(overlaps, exo_combined, endo_combined)
        
        # Convert to DataFrame
        enrichment_df = pd.DataFrame(enrichment_data)
        
        # Calculate significance
        if not enrichment_df.empty:
            # Calculate significance based on enrichment scores
            enrichment_threshold = np.percentile(enrichment_df['enrichment'], 75)  # Top 25%
            qvalue_threshold = 0.05
            
            enrichment_df['significant'] = (
                (enrichment_df['enrichment'] > enrichment_threshold) & 
                (enrichment_df['exo_qValue'] < qvalue_threshold) &
                (enrichment_df['endo_qValue'] < qvalue_threshold)
            )
            
            # Save results
            enrichment_df.to_csv(enrichment_file, index=False)
            
            print(f"Found {len(enrichment_df)} overlapping regions")
            print(f"Identified {enrichment_df['significant'].sum()} significant regions")
            
            # Create visualizations
            with ThreadPoolExecutor() as executor:
                executor.submit(plot_signal_distributions, enrichment_df, results_dir)
                executor.submit(plot_width_analysis, enrichment_df, results_dir)
    
    return enrichment_df

def analyze_by_expression_categories(dea_nsc, peaks_exo, peaks_endo, gene_annotations):
    """Parallelized version of expression category analysis"""
    print("\nAnalyzing Mecp2 binding by expression categories...")
    
    # Categorize genes (vectorized operation)
    dea_nsc['category'] = 'non-deregulated'
    mask_up = (dea_nsc['log2FoldChange'] > 1) & (dea_nsc['padj'] < 0.05)
    mask_down = (dea_nsc['log2FoldChange'] < -1) & (dea_nsc['padj'] < 0.05)
    dea_nsc.loc[mask_up, 'category'] = 'up-regulated'
    dea_nsc.loc[mask_down, 'category'] = 'down-regulated'
    
    # Process categories in parallel
    categories = {
        cat: dea_nsc[dea_nsc['category'] == cat]['gene'].tolist()
        for cat in ['non-deregulated', 'up-regulated', 'down-regulated']
    }
    
    category_results = {}
    for category, genes in categories.items():
        print(f"\nAnalyzing {category} genes ({len(genes)} genes)")
        category_data = parallel_analyze_categories(
            genes, gene_annotations, peaks_exo, peaks_endo)
        category_results[category] = pd.DataFrame(category_data)
        
        # Save category results
        if category_data:
            category_results[category].to_csv(
                f'{RESULTS_DIR}/mecp2_binding_{category.replace("-", "_")}.csv', 
                index=False
            )
    
    return category_results

def plot_category_results(category_results):
    """Create visualizations for expression category analysis"""
    # Binding frequency by category
    plt.figure(figsize=(10, 6))
    categories = []
    exo_freq = []
    endo_freq = []
    both_freq = []
    
    for category, df in category_results.items():
        total_genes = len(df)
        if total_genes > 0:
            categories.append(category)
            exo_freq.append(len(df[df['exo_peaks'] > 0]) / total_genes * 100)
            endo_freq.append(len(df[df['endo_peaks'] > 0]) / total_genes * 100)
            both_freq.append(len(df[(df['exo_peaks'] > 0) & (df['endo_peaks'] > 0)]) / total_genes * 100)
    
    x = np.arange(len(categories))
    width = 0.25
    
    plt.bar(x - width, exo_freq, width, label='Exo')
    plt.bar(x, endo_freq, width, label='Endo')
    plt.bar(x + width, both_freq, width, label='Both')
    
    plt.xlabel('Gene Category')
    plt.ylabel('Percentage of Genes')
    plt.title('Mecp2 Binding Frequency by Gene Category')
    plt.xticks(x, categories, rotation=45)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/binding_frequency_by_category.pdf')
    plt.close()
    
    # Enrichment distribution by category
    plt.figure(figsize=(10, 6))
    for category, df in category_results.items():
        if not df.empty:
            sns.kdeplot(data=df['enrichment'], label=category)
    
    plt.xlabel('Enrichment (Exo/Endo)')
    plt.ylabel('Density')
    plt.title('Mecp2 Enrichment Distribution by Gene Category')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/enrichment_distribution_by_category.pdf')
    plt.close()

def plot_enrichment_method_comparison(enrichment_df):
    """Create visualizations comparing different enrichment calculation methods"""
    plt.figure(figsize=(15, 10))
    
    # Plot correlations between methods
    methods = ['enrichment', 'enrichment_ww', 'enrichment_cs']
    method_names = ['Basic', 'Width-Weighted', 'Coverage Score']
    
    for i, (m1, n1) in enumerate(zip(methods, method_names)):
        for j, (m2, n2) in enumerate(zip(methods, method_names)):
            if i < j:
                plt.subplot(2, 2, i*2 + j)
                plt.scatter(enrichment_df[m1], enrichment_df[m2], alpha=0.5)
                plt.xlabel(f'{n1} Enrichment')
                plt.ylabel(f'{n2} Enrichment')
                
                # Add correlation coefficient
                corr = np.corrcoef(enrichment_df[m1], enrichment_df[m2])[0,1]
                plt.title(f'Correlation: {corr:.3f}')
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/enrichment_method_comparison.pdf')
    plt.close()

def process_sample(sample_data):
    """Process a single sample of peaks"""
    sample, peaks = sample_data
    if not peaks.empty:
        peaks = peaks.copy()
        peaks['sample'] = sample
        return peaks
    return None

def parallel_combine_peaks(peaks_dict, n_jobs=None):
    """Parallelized version of combine_peaks"""
    if not peaks_dict:
        return pd.DataFrame()
    
    # Set number of workers
    if n_jobs is None or n_jobs <= 0:
        n_jobs = max(1, mp.cpu_count() - 1)  # Use all cores except one, minimum 1
    
    # Process samples in parallel
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        processed_peaks = list(executor.map(process_sample, peaks_dict.items()))
        all_peaks = [peaks for peaks in processed_peaks if peaks is not None]
    
    if not all_peaks:
        return pd.DataFrame()
    
    # Combine results
    combined = pd.concat(all_peaks, ignore_index=True)
    
    # Group and calculate metrics (this part is fast enough to keep sequential)
    grouped = combined.groupby(['chr', 'start', 'end'])
    
    # Use numpy operations for faster calculations
    metrics = {
        'signal_basic': grouped['signalValue'].median(),  # Using median instead of mean
        'qValue': grouped['qValue'].min(),
        'num_replicates': grouped.size()
    }
    
    combined_peaks = pd.DataFrame(metrics).reset_index()
    combined_peaks['width'] = combined_peaks['end'] - combined_peaks['start']
    
    # Vectorized calculations for different signal metrics
    # Width-weighted signal
    combined_peaks['signal_ww'] = combined_peaks['signal_basic'] * combined_peaks['width']
    
    # Coverage score
    total_width = combined_peaks['width'].sum()
    combined_peaks['signal_cs'] = combined_peaks['signal_ww'] / total_width
    
    # Additional normalizations
    for signal_col in ['signal_basic', 'signal_ww', 'signal_cs']:
        # Log transform
        combined_peaks[f'{signal_col}_log'] = np.log1p(combined_peaks[signal_col])
        
        # Z-score normalization
        combined_peaks[f'{signal_col}_zscore'] = stats.zscore(
            combined_peaks[signal_col], nan_policy='omit')
        
        # Quantile normalization
        combined_peaks[f'{signal_col}_quantile'] = stats.rankdata(
            combined_peaks[signal_col]) / len(combined_peaks)
    
    return combined_peaks

def plot_signal_distributions(enrichment_df, results_dir):
    """Plot distributions of signal values and enrichment scores"""
    plt.figure(figsize=(15, 10))
    
    # Plot 1: Distribution of exo vs endo signals
    plt.subplot(2, 2, 1)
    plt.hist(enrichment_df['exo_signal'], bins=50, alpha=0.5, label='Exo', density=True)
    plt.hist(enrichment_df['endo_signal'], bins=50, alpha=0.5, label='Endo', density=True)
    plt.xlabel('Signal Value')
    plt.ylabel('Density')
    plt.title('Distribution of Signal Values')
    plt.legend()
    
    # Plot 2: Distribution of enrichment scores
    plt.subplot(2, 2, 2)
    plt.hist(np.log2(enrichment_df['enrichment'].clip(min=0.01)), bins=50)
    plt.xlabel('log2(Enrichment Score)')
    plt.ylabel('Count')
    plt.title('Distribution of Enrichment Scores')
    
    # Plot 3: Signal correlation
    plt.subplot(2, 2, 3)
    plt.scatter(np.log2(enrichment_df['endo_signal'].clip(min=0.01)), 
               np.log2(enrichment_df['exo_signal'].clip(min=0.01)), 
               alpha=0.5)
    plt.xlabel('log2(Endo Signal)')
    plt.ylabel('log2(Exo Signal)')
    plt.title('Exo vs Endo Signal Correlation')
    
    # Plot 4: QValue distribution
    plt.subplot(2, 2, 4)
    plt.hist(-np.log10(enrichment_df['exo_qValue'].clip(max=1)), 
             bins=50, alpha=0.5, label='Exo')
    plt.hist(-np.log10(enrichment_df['endo_qValue'].clip(max=1)), 
             bins=50, alpha=0.5, label='Endo')
    plt.xlabel('-log10(qValue)')
    plt.ylabel('Count')
    plt.title('Distribution of qValues')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(f'{results_dir}/signal_distributions.pdf')
    plt.close()

def plot_width_analysis(enrichment_df, results_dir):
    """Plot analysis of peak widths and their relationship to enrichment"""
    plt.figure(figsize=(15, 10))
    
    # Calculate peak widths
    enrichment_df['width'] = enrichment_df['end'] - enrichment_df['start']
    
    # Plot 1: Width distribution
    plt.subplot(2, 2, 1)
    plt.hist(np.log10(enrichment_df['width']), bins=50)
    plt.xlabel('log10(Peak Width)')
    plt.ylabel('Count')
    plt.title('Distribution of Peak Widths')
    
    # Plot 2: Width vs Enrichment
    plt.subplot(2, 2, 2)
    plt.scatter(np.log10(enrichment_df['width']), 
               np.log2(enrichment_df['enrichment'].clip(min=0.01)), 
               alpha=0.5)
    plt.xlabel('log10(Peak Width)')
    plt.ylabel('log2(Enrichment)')
    plt.title('Peak Width vs Enrichment')
    
    # Plot 3: Width vs Signal
    plt.subplot(2, 2, 3)
    plt.scatter(np.log10(enrichment_df['width']), 
               np.log2(enrichment_df['exo_signal'].clip(min=0.01)), 
               alpha=0.5, label='Exo')
    plt.scatter(np.log10(enrichment_df['width']), 
               np.log2(enrichment_df['endo_signal'].clip(min=0.01)), 
               alpha=0.5, label='Endo')
    plt.xlabel('log10(Peak Width)')
    plt.ylabel('log2(Signal)')
    plt.title('Peak Width vs Signal')
    plt.legend()
    
    # Plot 4: Width vs Replicate Count
    plt.subplot(2, 2, 4)
    plt.scatter(np.log10(enrichment_df['width']), 
               enrichment_df['exo_replicates'], 
               alpha=0.5, label='Exo')
    plt.scatter(np.log10(enrichment_df['width']), 
               enrichment_df['endo_replicates'], 
               alpha=0.5, label='Endo')
    plt.xlabel('log10(Peak Width)')
    plt.ylabel('Number of Replicates')
    plt.title('Peak Width vs Replicate Count')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(f'{results_dir}/width_analysis.pdf')
    plt.close()

def analyze_peak_overlap_in_promoters(peaks_exo, peaks_endo, gene_annotations, name_to_info, window=PROMOTER_WINDOW):
    """
    Analyze peaks in promoter regions to identify common, exo-only, and endo-only peaks.
    Returns three DataFrames containing the peaks in each category.
    """
    print("\nAnalyzing peak overlap in promoter regions...")
    
    # Initialize lists to store results
    common_peaks = []
    exo_only_peaks = []
    endo_only_peaks = []
    
    # Process each gene in the annotations
    for _, gene_info in tqdm(gene_annotations.iterrows(), desc="Processing genes"):
        # Define promoter region
        if gene_info['strand'] == '+':
            promoter_start = max(0, gene_info['start'] - window)
            promoter_end = gene_info['start'] + window
        else:
            promoter_start = max(0, gene_info['end'] - window)
            promoter_end = gene_info['end'] + window
            
        # Get all exo peaks in this promoter
        exo_promoter_peaks = []
        for sample, peaks in peaks_exo.items():
            sample_peaks = peaks[
                (peaks['chr'] == gene_info['chr']) &
                (peaks['start'] <= promoter_end) &
                (peaks['end'] >= promoter_start)
            ]
            if not sample_peaks.empty:
                exo_promoter_peaks.append(sample_peaks)
        
        # Get all endo peaks in this promoter
        endo_promoter_peaks = []
        for sample, peaks in peaks_endo.items():
            sample_peaks = peaks[
                (peaks['chr'] == gene_info['chr']) &
                (peaks['start'] <= promoter_end) &
                (peaks['end'] >= promoter_start)
            ]
            if not sample_peaks.empty:
                endo_promoter_peaks.append(sample_peaks)
        
        # Convert to DataFrames if we found any peaks
        exo_peaks = pd.concat(exo_promoter_peaks) if exo_promoter_peaks else pd.DataFrame()
        endo_peaks = pd.concat(endo_promoter_peaks) if endo_promoter_peaks else pd.DataFrame()
        
        # Record the peak presence for this promoter
        if not exo_peaks.empty and not endo_peaks.empty:
            # Calculate mean signal values if multiple peaks exist
            mean_exo_signal = exo_peaks['signalValue'].mean()
            mean_endo_signal = endo_peaks['signalValue'].mean()
            
            common_peaks.append({
                'gene': gene_info['gene_name'],
                'chr': gene_info['chr'],
                'promoter_start': promoter_start,
                'promoter_end': promoter_end,
                'strand': gene_info['strand'],
                'exo_peaks_count': len(exo_peaks),
                'endo_peaks_count': len(endo_peaks),
                'mean_exo_signal': mean_exo_signal,
                'mean_endo_signal': mean_endo_signal,
                'exo_coordinates': ';'.join([f"{row['start']}-{row['end']}" for _, row in exo_peaks.iterrows()]),
                'endo_coordinates': ';'.join([f"{row['start']}-{row['end']}" for _, row in endo_peaks.iterrows()])
            })
        
        elif not exo_peaks.empty:
            exo_only_peaks.append({
                'gene': gene_info['gene_name'],
                'chr': gene_info['chr'],
                'promoter_start': promoter_start,
                'promoter_end': promoter_end,
                'strand': gene_info['strand'],
                'peaks_count': len(exo_peaks),
                'mean_signal': exo_peaks['signalValue'].mean(),
                'coordinates': ';'.join([f"{row['start']}-{row['end']}" for _, row in exo_peaks.iterrows()])
            })
            
        elif not endo_peaks.empty:
            endo_only_peaks.append({
                'gene': gene_info['gene_name'],
                'chr': gene_info['chr'],
                'promoter_start': promoter_start,
                'promoter_end': promoter_end,
                'strand': gene_info['strand'],
                'peaks_count': len(endo_peaks),
                'mean_signal': endo_peaks['signalValue'].mean(),
                'coordinates': ';'.join([f"{row['start']}-{row['end']}" for _, row in endo_peaks.iterrows()])
            })
    
    # Convert lists to DataFrames
    common_df = pd.DataFrame(common_peaks)
    exo_only_df = pd.DataFrame(exo_only_peaks)
    endo_only_df = pd.DataFrame(endo_only_peaks)
    
    # Save results
    if not common_df.empty:
        common_df.to_csv(f'{RESULTS_DIR}/common_promoter_peaks.csv', index=False)
    if not exo_only_df.empty:
        exo_only_df.to_csv(f'{RESULTS_DIR}/exo_only_promoter_peaks.csv', index=False)
    if not endo_only_df.empty:
        endo_only_df.to_csv(f'{RESULTS_DIR}/endo_only_promoter_peaks.csv', index=False)
    
    # Print summary
    print(f"\nSummary of promoter peak analysis:")
    print(f"Common peaks found in {len(common_df)} promoters")
    print(f"Exo-only peaks found in {len(exo_only_df)} promoters")
    print(f"Endo-only peaks found in {len(endo_only_df)} promoters")
    
    return common_df, exo_only_df, endo_only_df

# Modify the main execution block
if __name__ == "__main__":
    # Create output directory if it doesn't exist
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Set number of workers (minimum 1)
    N_JOBS = max(1, mp.cpu_count() - 1)  # Leave one core free
    print(f"Using {N_JOBS} workers for parallel processing")

    # Load data
    dea, peaks_exo, peaks_endo, gene_annotations, name_to_info = load_data()

    # Approach 1: Independent Mecp2 enrichment analysis
    enrichment_df = analyze_mecp2_enrichment_independent(
        peaks_exo, peaks_endo, n_jobs=N_JOBS, results_dir=RESULTS_DIR
    )
    
    # Integrate with RNA-seq data
    if 'significant' in enrichment_df.columns:
        significant_regions = enrichment_df[enrichment_df['significant']]
        genes_near_regions = []
        
        # Find genes near significant regions
        for _, region in significant_regions.iterrows():
            nearby_genes = gene_annotations[
                (gene_annotations['chr'] == region['chr']) &
                ((gene_annotations['start'] - 2000 <= region['end']) &
                 (gene_annotations['end'] + 2000 >= region['start']))
            ]['gene_name'].tolist()
            
            genes_near_regions.extend(nearby_genes)
        
        # Get expression data for these genes
        genes_with_expression = dea[dea['gene'].isin(genes_near_regions)].copy()
        genes_with_expression.to_csv(f'{RESULTS_DIR}/genes_near_enriched_regions.csv', index=False)
    else:
        print("Warning: No significance information found in enrichment results")
        genes_near_regions = []

    # Approach 2: Analysis by expression categories
    category_results = analyze_by_expression_categories(dea, peaks_exo, peaks_endo, gene_annotations)
    plot_category_results(category_results)

    # Generate summary statistics
    summary_stats = {
        'total_enriched_regions': len(enrichment_df),
        'significant_regions': len(significant_regions) if 'significant' in enrichment_df.columns else 0,
        'genes_near_enriched_regions': len(set(genes_near_regions)),
        'up_regulated_targets': len(category_results['up-regulated']),
        'down_regulated_targets': len(category_results['down-regulated']),
        'non_deregulated_targets': len(category_results['non-deregulated'])
    }
    
    with open(f'{RESULTS_DIR}/analysis_summary.txt', 'w') as f:
        for key, value in summary_stats.items():
            f.write(f'{key}: {value}\n')

    # Generate enrichment analysis
    enrichment_df = analyze_mecp2_enrichment_independent(
        peaks_exo, peaks_endo, n_jobs=N_JOBS, results_dir=RESULTS_DIR
    )
    
    # Plot method comparisons
    plot_enrichment_method_comparison(enrichment_df)
    
    # Calculate summary statistics for each method
    summary_stats.update({
        'median_basic_enrichment': enrichment_df['enrichment'].median(),
        'median_ww_enrichment': enrichment_df['enrichment_ww'].median(),
        'median_cs_enrichment': enrichment_df['enrichment_cs'].median()
    })

    # Analyze peak overlap in promoters
    # common_peaks, exo_only_peaks, endo_only_peaks = analyze_peak_overlap_in_promoters(
    #     peaks_exo, peaks_endo, gene_annotations, name_to_info
    # )
