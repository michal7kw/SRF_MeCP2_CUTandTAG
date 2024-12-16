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

def load_data():
    try:
        # Load DEA results
        dea_nsc = pd.read_csv("../DATA/DEA_NSC.csv")
        print("\nDEA file columns:", dea_nsc.columns.tolist())
        
        # Add timeout for file operations
        def load_with_timeout(filepath, timeout=30):
            """Load file with timeout"""
            start_time = time.time()
            while not os.path.exists(filepath):
                if time.time() - start_time > timeout:
                    print(f"Timeout waiting for file: {filepath}")
                    return None
                time.sleep(1)
            return filepath

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
        
        # Load exogenous samples
        for sample in exo_samples:
            peak_file = f"{RESULTS_DIR}/peaks/{sample}_peaks.narrowPeak"
            peaks_exo[sample] = load_peak_file(peak_file, sample)
            depths_exo[sample] = calculate_sequencing_depth_safe(f"{DATA_DIR}/aligned/{sample}.bam")
            
            if not peaks_exo[sample].empty:
                peaks_exo[sample]['signalValue'] = peaks_exo[sample]['signalValue'].clip(lower=0) * (1e6 / depths_exo[sample])
        
        # Load endogenous samples
        for sample in endo_samples:
            peak_file = f"{RESULTS_DIR}/peaks/{sample}_peaks.narrowPeak"
            peaks_endo[sample] = load_peak_file(peak_file, sample)
            depths_endo[sample] = calculate_sequencing_depth_safe(f"{DATA_DIR}/aligned/{sample}.bam")
            
            if not peaks_endo[sample].empty:
                peaks_endo[sample]['signalValue'] = peaks_endo[sample]['signalValue'].clip(lower=0) * (1e6 / depths_endo[sample])
        
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
    
    # Check if required fields are present
    required_fields = ['chr', 'start', 'end', 'strand']
    if not all(field in gene_info for field in required_fields):
        # Instead of printing warning, silently return empty DataFrame
        return pd.DataFrame()
    
    # Ensure chromosome information is valid
    if pd.isna(gene_info['chr']) or not gene_info['chr']:
        return pd.DataFrame()
    
    try:
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
            if peaks.empty:
                continue
            
            # Ensure required columns exist
            if not all(col in peaks.columns for col in ['chr', 'start', 'end', 'signalValue']):
                continue
                
            # Filter peaks for the same chromosome first
            chr_peaks = peaks[peaks['chr'] == gene_info['chr']]
            
            if chr_peaks.empty:
                continue
            
            # Find overlapping peaks
            overlapping = chr_peaks.loc[
                (chr_peaks['start'] <= promoter_end) & 
                (chr_peaks['end'] >= promoter_start)
            ].copy()
            
            if not overlapping.empty:
                # Ensure signalValue is numeric and non-negative
                overlapping['signalValue'] = pd.to_numeric(overlapping['signalValue'], errors='coerce')
                overlapping['signalValue'] = overlapping['signalValue'].fillna(0).clip(lower=0)
                overlapping['sample'] = sample
                all_peaks.append(overlapping)
        
        if not all_peaks:
            return pd.DataFrame()
        
        combined_peaks = pd.concat(all_peaks, ignore_index=True)
        return combined_peaks
        
    except Exception as e:
        # Instead of printing error, silently return empty DataFrame
        return pd.DataFrame()

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
            
            # Get peaks near gene
            exo_peaks = get_peaks_near_gene(gene, peaks_exo, gene_annotations, name_to_info)
            endo_peaks = get_peaks_near_gene(gene, peaks_endo, gene_annotations, name_to_info)
            
            if not exo_peaks.empty or not endo_peaks.empty:
                peaks_found += 1
                try:
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
                except Exception:
                    # Silently continue if there's an error processing this gene
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
    # Process exogenous samples
    for sample, peaks in peaks_exo.items():
        peaks_bed = BedTool.from_dataframe(peaks[['chr', 'start', 'end', 'signalValue']])
        intersect = peaks_bed.intersect(cpg_bed, wa=True, wb=True)
        # Convert intersection results to DataFrame with named columns
        # Parse the intersection results into a list of lists
        intersect_data = [line.strip().split('\t') for line in str(intersect).strip().split('\n')]
        if intersect_data and intersect_data[0]:  # Check if we have any data
            results[f"{sample}_cpg"] = pd.DataFrame(
                intersect_data,
                columns=['peak_chr', 'peak_start', 'peak_end', 'signalValue',
                        'cpg_chr', 'cpg_start', 'cpg_end']
            )
            # Convert numeric columns
            results[f"{sample}_cpg"]['signalValue'] = pd.to_numeric(results[f"{sample}_cpg"]['signalValue'])
            results[f"{sample}_cpg"]['peak_start'] = pd.to_numeric(results[f"{sample}_cpg"]['peak_start'])
            results[f"{sample}_cpg"]['peak_end'] = pd.to_numeric(results[f"{sample}_cpg"]['peak_end'])
            results[f"{sample}_cpg"]['cpg_start'] = pd.to_numeric(results[f"{sample}_cpg"]['cpg_start'])
            results[f"{sample}_cpg"]['cpg_end'] = pd.to_numeric(results[f"{sample}_cpg"]['cpg_end'])
        else:
            results[f"{sample}_cpg"] = pd.DataFrame(
                columns=['peak_chr', 'peak_start', 'peak_end', 'signalValue',
                        'cpg_chr', 'cpg_start', 'cpg_end']
            )
    
    # Process endogenous samples
    for sample, peaks in peaks_endo.items():
        peaks_bed = BedTool.from_dataframe(peaks[['chr', 'start', 'end', 'signalValue']])
        intersect = peaks_bed.intersect(cpg_bed, wa=True, wb=True)
        # Convert intersection results to DataFrame with named columns
        intersect_data = [line.strip().split('\t') for line in str(intersect).strip().split('\n')]
        if intersect_data and intersect_data[0]:  # Check if we have any data
            results[f"{sample}_cpg"] = pd.DataFrame(
                intersect_data,
                columns=['peak_chr', 'peak_start', 'peak_end', 'signalValue',
                        'cpg_chr', 'cpg_start', 'cpg_end']
            )
            # Convert numeric columns
            results[f"{sample}_cpg"]['signalValue'] = pd.to_numeric(results[f"{sample}_cpg"]['signalValue'])
            results[f"{sample}_cpg"]['peak_start'] = pd.to_numeric(results[f"{sample}_cpg"]['peak_start'])
            results[f"{sample}_cpg"]['peak_end'] = pd.to_numeric(results[f"{sample}_cpg"]['peak_end'])
            results[f"{sample}_cpg"]['cpg_start'] = pd.to_numeric(results[f"{sample}_cpg"]['cpg_start'])
            results[f"{sample}_cpg"]['cpg_end'] = pd.to_numeric(results[f"{sample}_cpg"]['cpg_end'])
        else:
            results[f"{sample}_cpg"] = pd.DataFrame(
                columns=['peak_chr', 'peak_start', 'peak_end', 'signalValue',
                        'cpg_chr', 'cpg_start', 'cpg_end']
            )
    
    # Calculate enrichment metrics for each CpG island
    cpg_enrichment = []
    
    for _, cpg in cpg_islands.iterrows():
        exo_signal = 0
        endo_signal = 0
        
        # Sum signals for exo samples
        for sample in peaks_exo.keys():
            if f"{sample}_cpg" in results:
                sample_peaks = results[f"{sample}_cpg"]
                mask = (sample_peaks['cpg_chr'] == cpg['chr']) & \
                       (sample_peaks['cpg_start'] == cpg['start']) & \
                       (sample_peaks['cpg_end'] == cpg['end'])
                if not sample_peaks[mask].empty:
                    exo_signal += sample_peaks[mask]['signalValue'].sum()
        
        # Sum signals for endo samples
        for sample in peaks_endo.keys():
            if f"{sample}_cpg" in results:
                sample_peaks = results[f"{sample}_cpg"]
                mask = (sample_peaks['cpg_chr'] == cpg['chr']) & \
                       (sample_peaks['cpg_start'] == cpg['start']) & \
                       (sample_peaks['cpg_end'] == cpg['end'])
                if not sample_peaks[mask].empty:
                    endo_signal += sample_peaks[mask]['signalValue'].sum()
        
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
    
    # Create a BedTool object for genes with additional name field
    genes_df = gene_annotations[['chr', 'start', 'end', 'gene_name']].copy()
    # Ensure all fields are strings as required by BedTool
    genes_df['chr'] = genes_df['chr'].astype(str)
    genes_df['start'] = genes_df['start'].astype(str)
    genes_df['end'] = genes_df['end'].astype(str)
    genes_df['gene_name'] = genes_df['gene_name'].astype(str)
    genes_bed = BedTool.from_dataframe(genes_df)
    
    # Find overlaps between CpG islands and genes
    overlaps = cpg_bed.intersect(genes_bed, wa=True, wb=True)
    
    # Create integrated dataset
    integrated_data = []
    
    # Parse intersection results
    for line in str(overlaps).strip().split('\n'):
        if not line:  # Skip empty lines
            continue
            
        fields = line.split('\t')
        if len(fields) < 7:  # Ensure we have enough fields
            continue
            
        cpg_chr, cpg_start, cpg_end = fields[0:3]
        gene_name = fields[6]  # Gene name should be the last field
        
        # Find matching CpG enrichment data
        cpg_data = enrichment_df[
            (enrichment_df['chr'] == cpg_chr) &
            (enrichment_df['start'] == int(cpg_start)) &
            (enrichment_df['end'] == int(cpg_end))
        ]
        
        if not cpg_data.empty:
            # Find matching DEA data
            dea_match = dea_nsc[dea_nsc['gene'] == gene_name]
            
            if not dea_match.empty:
                dea_info = dea_match.iloc[0]
                
                integrated_data.append({
                    'chr': cpg_chr,
                    'cpg_start': int(cpg_start),
                    'cpg_end': int(cpg_end),
                    'gene': gene_name,
                    'enrichment': float(cpg_data.iloc[0]['enrichment']),
                    'exo_signal': float(cpg_data.iloc[0]['exo_signal']),
                    'endo_signal': float(cpg_data.iloc[0]['endo_signal']),
                    'log2FoldChange': float(dea_info['log2FoldChange']),
                    'padj': float(dea_info['padj'])
                })
    
    integrated_df = pd.DataFrame(integrated_data)
    
    # Save to CSV
    if not integrated_df.empty:
        integrated_df.to_csv(f'{RESULTS_DIR}/cpg_rna_integrated_NSC.csv', index=False)
        print(f"Integrated {len(integrated_df)} CpG islands with gene expression data")
    else:
        print("Warning: No integrated data found")
    
    return integrated_df

# if __name__ == "__main__":
# Create output directory if it doesn't exist
os.makedirs('results', exist_ok=True)

# Load data
dea, peaks_exo, peaks_endo, gene_annotations, name_to_info = load_data()

# Load CpG islands
cpg_islands = load_cpg_islands()

# Analyze CpG enrichment
enrichment_df = analyze_cpg_enrichment(peaks_exo, peaks_endo, cpg_islands)

# Create visualizations
plot_cpg_enrichment(enrichment_df)

# Integrate with RNA-seq data
integrated_df = integrate_with_rna_seq(enrichment_df, dea, gene_annotations)

# Plot peak width distributions
plot_peak_width_distributions(peaks_exo, peaks_endo)
plot_detailed_peak_width_distributions(peaks_exo, peaks_endo)

# Run analysis
results = analyze_enrichment(dea, peaks_exo, peaks_endo, gene_annotations, name_to_info)

# Create visualizations
plot_enrichment(results)

# Generate summary statistics
summarize_results(results)

# Save detailed results to CSV
for method_name, df in results.items():
    df.to_csv(f'{RESULTS_DIR}/enrichment_{method_name}_NSC.csv', index=False) 

# Plot width vs enrichment
plot_width_vs_enrichment(results, peaks_exo, peaks_endo, gene_annotations, name_to_info)

# Summarize peak distribution
peak_distribution = summarize_peak_distribution(results)
print("\nPeak Distribution Summary:")
print(peak_distribution)
