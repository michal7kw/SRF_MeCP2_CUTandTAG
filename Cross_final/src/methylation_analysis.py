"""
Analysis of methylation patterns and SMARCB1 binding in different gene sets.
"""

import os
import pandas as pd
import numpy as np
import pyBigWig
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Tuple
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

from config import *
from genomic_utils import *

class MethylationAnalyzer:
    def __init__(self):
        """Initialize the analyzer with gene lists and data paths."""
        self.gene_lists = self._load_gene_lists()
        print("Loading gene coordinates...")
        self.gene_coords = load_gtf(GTF_FILE)
        print(f"Loaded coordinates for {len(self.gene_coords)} genes")
        
    def _load_gene_lists(self) -> Dict[str, List[str]]:
        """Load gene lists from CSV files."""
        gene_lists = {}
        for category, filepath in GENE_LISTS.items():
            try:
                df = pd.read_csv(filepath)
                if 'GENE' in df.columns:
                    genes = df['GENE'].tolist()
                else:
                    df = pd.read_csv(filepath, header=None)
                    genes = df[0].tolist()
                gene_lists[category] = genes
                print(f"Loaded {len(genes)} genes from {category} list")
            except Exception as e:
                print(f"Error loading {filepath}: {str(e)}")
                gene_lists[category] = []
        return gene_lists
    
    def calculate_methylation_profile(self, bigwig_files: List[str], 
                                    genes: List[str], 
                                    window: Tuple[int, int] = PROMOTER_WINDOW) -> pd.DataFrame:
        """Calculate methylation profile around TSS for a set of genes."""
        profiles = []
        window_size = window[1] - window[0]
        
        genes = [g for g in genes if g in self.gene_coords.index]
        print(f"Analyzing methylation for {len(genes)} genes")
        
        for bw_file in bigwig_files:
            print(f"Processing {os.path.basename(bw_file)}")
            try:
                with pyBigWig.open(bw_file) as bw:
                    for gene in genes:
                        chrom = self.gene_coords.loc[gene, 'chromosome']
                        tss = self.gene_coords.loc[gene, 'tss']
                        strand = self.gene_coords.loc[gene, 'strand']
                        
                        start = max(0, tss + window[0])
                        end = tss + window[1]
                        
                        try:
                            values = bw.values(chrom, start, end)
                            values = normalize_bigwig_values(values)
                            if len(values) == window_size:
                                profiles.append({
                                    'gene': gene,
                                    'sample': os.path.basename(bw_file),
                                    'profile': values,
                                    'mean_signal': np.mean(values)
                                })
                        except Exception as e:
                            continue
            except Exception as e:
                print(f"Error opening bigwig file {bw_file}: {str(e)}")
                continue
                            
        return pd.DataFrame(profiles)

    def analyze_methylation_differences(self):
        """Analyze methylation differences between gene categories."""
        results = {}
        stats_results = {}
        all_profiles = []  # Store all profiles for box plots
        
        for condition in MEDIP_BIGWIG.keys():
            print(f"\nAnalyzing methylation for condition: {condition}")
            methylation_data = {}
            condition_stats = {}
            
            for category, genes in self.gene_lists.items():
                print(f"\nProcessing {category} genes...")
                profiles = self.calculate_methylation_profile(
                    MEDIP_BIGWIG[condition], genes
                )
                methylation_data[category] = profiles
                
                # Store profiles for box plots
                if len(profiles) > 0:
                    profiles['condition'] = condition
                    profiles['category'] = category
                    all_profiles.append(profiles)
                
                # Calculate basic statistics
                if len(profiles) > 0:
                    condition_stats[category] = {
                        'mean_signal': profiles['mean_signal'].mean(),
                        'std_signal': profiles['mean_signal'].std(),
                        'n_genes': len(set(profiles['gene']))
                    }
            
            # Perform statistical tests between categories
            if len(methylation_data.get('up', [])) > 0 and len(methylation_data.get('down', [])) > 0:
                stat, pval = stats.mannwhitneyu(
                    methylation_data['up']['mean_signal'],
                    methylation_data['down']['mean_signal']
                )
                condition_stats['up_vs_down'] = {'statistic': stat, 'pvalue': pval}
            
            results[condition] = methylation_data
            stats_results[condition] = condition_stats
        
        # Combine all profiles for plotting
        if all_profiles:
            combined_profiles = pd.concat(all_profiles, ignore_index=True)
            self.plot_methylation_summary(combined_profiles)
            
        return results, stats_results

    def plot_methylation_summary(self, profiles_df: pd.DataFrame):
        """Create summary plots for methylation data."""
        os.makedirs(RESULTS_DIR, exist_ok=True)
        
        # 1. Box plot of mean methylation levels by condition and category
        plt.figure(figsize=(12, 6))
        sns.boxplot(data=profiles_df, x='condition', y='mean_signal', hue='category')
        plt.title('Distribution of Mean Methylation Signals')
        plt.xlabel('Condition')
        plt.ylabel('Mean Methylation Signal (RPM)')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, 'methylation_boxplot.pdf'))
        plt.close()
        
        # 2. Bar plot of mean methylation with error bars
        summary_stats = profiles_df.groupby(['condition', 'category'])['mean_signal'].agg(['mean', 'sem']).reset_index()
        plt.figure(figsize=(12, 6))
        x = np.arange(len(summary_stats['condition'].unique()))
        width = 0.25
        
        categories = summary_stats['category'].unique()
        for i, category in enumerate(categories):
            category_data = summary_stats[summary_stats['category'] == category]
            plt.bar(x + i*width, category_data['mean'], width, 
                   label=category,
                   yerr=category_data['sem'],
                   capsize=5)
        
        plt.xlabel('Condition')
        plt.ylabel('Mean Methylation Signal (RPM)')
        plt.title('Mean Methylation Signal by Condition and Category')
        plt.xticks(x + width, summary_stats['condition'].unique())
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, 'methylation_barplot.pdf'))
        plt.close()
        
        # 3. Violin plot to show distribution
        plt.figure(figsize=(12, 6))
        sns.violinplot(data=profiles_df, x='condition', y='mean_signal', hue='category')
        plt.title('Distribution of Methylation Signals')
        plt.xlabel('Condition')
        plt.ylabel('Mean Methylation Signal (RPM)')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, 'methylation_violin.pdf'))
        plt.close()

    def analyze_smarcb1_binding(self):
        """Analyze SMARCB1 binding in different gene categories."""
        print("\n=== Starting SMARCB1 binding analysis ===")
        results = {}
        enrichment_results = {}
        all_overlaps = []  # Store all overlaps for combined analysis
        
        for condition in ['BG', 'BM']:
            print(f"\n=== Processing condition: {condition} ===")
            binding_data = {}
            condition_enrichment = {}
            
            for category, genes in self.gene_lists.items():
                print(f"\n--- Processing {category} genes ---")
                print(f"Number of genes in category: {len(genes)}")
                gene_coords = self.gene_coords.loc[
                    self.gene_coords.index.isin(genes)
                ]
                print(f"Found coordinates for {len(gene_coords)} genes")
                
                promoter_coords = get_promoter_coords(
                    gene_coords, PROMOTER_WINDOW
                )
                print(f"Generated promoter coordinates for {len(promoter_coords)} genes")
                
                print(f"Looking for peaks in: {SMARCB1_PEAKS[condition]}")
                overlaps = calculate_peak_overlaps(
                    SMARCB1_PEAKS[condition], promoter_coords
                )
                print(f"Found {len(overlaps)} peak overlaps")
                binding_data[category] = overlaps
                
                # Store overlaps for combined analysis
                if not overlaps.empty:
                    overlaps['condition'] = condition
                    overlaps['category'] = category
                    all_overlaps.append(overlaps)
                    print(f"Added overlaps to combined analysis. Total datasets: {len(all_overlaps)}")
                
                # Calculate enrichment statistics
                if len(overlaps) > 0:
                    n_genes_with_peaks = len(set(overlaps['gene']))
                    total_genes = len(genes)
                    mean_overlap = overlaps['overlap_ratio'].mean()
                    std_overlap = overlaps['overlap_ratio'].std()
                    
                    print(f"Statistics for {category}:")
                    print(f"  - Genes with peaks: {n_genes_with_peaks}")
                    print(f"  - Total genes: {total_genes}")
                    print(f"  - Fraction bound: {n_genes_with_peaks/total_genes:.3f}")
                    print(f"  - Mean overlap: {mean_overlap:.3f}")
                    print(f"  - Std overlap: {std_overlap:.3f}")
                    
                    condition_enrichment[category] = {
                        'n_genes_with_peaks': n_genes_with_peaks,
                        'total_genes': total_genes,
                        'fraction_bound': n_genes_with_peaks / total_genes,
                        'mean_overlap': mean_overlap,
                        'std_overlap': std_overlap
                    }
            
            results[condition] = binding_data
            enrichment_results[condition] = condition_enrichment
            print(f"\nCompleted analysis for condition: {condition}")
            print(f"Number of categories with data: {len(condition_enrichment)}")
        
        # Create combined visualizations if we have data
        if all_overlaps:
            print("\n=== Creating combined visualizations ===")
            combined_overlaps = pd.concat(all_overlaps, ignore_index=True)
            print(f"Combined data shape: {combined_overlaps.shape}")
            self.plot_smarcb1_summary(combined_overlaps)
        else:
            print("\nWARNING: No overlap data available for visualization")
            
        return results, enrichment_results

    def plot_smarcb1_summary(self, overlaps_df: pd.DataFrame):
        """Create summary plots for SMARCB1 binding data."""
        os.makedirs(RESULTS_DIR, exist_ok=True)
        
        # 1. Box plot of overlap ratios
        plt.figure(figsize=(12, 6))
        sns.boxplot(data=overlaps_df, x='category', y='overlap_ratio', hue='condition')
        plt.title('Distribution of SMARCB1 Peak Overlap Ratios')
        plt.xlabel('Gene Category')
        plt.ylabel('Overlap Ratio')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, 'smarcb1_overlap_boxplot.pdf'))
        plt.close()
        
        # 2. Violin plot of overlap ratios
        plt.figure(figsize=(12, 6))
        sns.violinplot(data=overlaps_df, x='category', y='overlap_ratio', hue='condition')
        plt.title('Distribution of SMARCB1 Peak Overlap Ratios')
        plt.xlabel('Gene Category')
        plt.ylabel('Overlap Ratio')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, 'smarcb1_overlap_violin.pdf'))
        plt.close()
        
        # 3. Count plot of peaks per gene
        peak_counts = overlaps_df.groupby(['condition', 'category', 'gene']).size().reset_index(name='peak_count')
        plt.figure(figsize=(12, 6))
        sns.boxplot(data=peak_counts, x='category', y='peak_count', hue='condition')
        plt.title('Number of SMARCB1 Peaks per Gene')
        plt.xlabel('Gene Category')
        plt.ylabel('Number of Peaks')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, 'smarcb1_peaks_per_gene.pdf'))
        plt.close()
        
        # 4. Overlap length distribution
        plt.figure(figsize=(12, 6))
        sns.boxplot(data=overlaps_df, x='category', y='overlap_length', hue='condition')
        plt.title('Distribution of SMARCB1 Peak Overlap Lengths')
        plt.xlabel('Gene Category')
        plt.ylabel('Overlap Length (bp)')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, 'smarcb1_overlap_length.pdf'))
        plt.close()
        
        # 5. Statistical summary
        stats_summary = overlaps_df.groupby(['condition', 'category']).agg({
            'overlap_ratio': ['mean', 'std', 'count'],
            'overlap_length': ['mean', 'std'],
            'gene': 'nunique'
        }).round(3)
        
        stats_summary.columns = ['_'.join(col).strip() for col in stats_summary.columns.values]
        stats_summary = stats_summary.reset_index()
        stats_summary.to_csv(os.path.join(RESULTS_DIR, 'smarcb1_detailed_stats.csv'), index=False)

    def plot_smarcb1_enrichment(self, results: Dict, enrichment: Dict):
        """Plot SMARCB1 binding enrichment for different gene categories."""
        os.makedirs(RESULTS_DIR, exist_ok=True)
        
        # Plot binding fractions
        plot_data = []
        for condition in results.keys():
            for category in enrichment[condition].keys():
                stats = enrichment[condition][category]
                plot_data.append({
                    'condition': condition,
                    'category': category,
                    'fraction_bound': stats['fraction_bound'],
                    'n_genes': stats['total_genes'],
                    'mean_overlap': stats['mean_overlap'],
                    'std_overlap': stats['std_overlap']
                })
        
        if plot_data:
            df = pd.DataFrame(plot_data)
            
            # Plot binding fractions
            plt.figure(figsize=(10, 6))
            sns.barplot(data=df, x='category', y='fraction_bound', hue='condition')
            plt.title('Fraction of Genes with SMARCB1 Binding')
            plt.xlabel('Gene Category')
            plt.ylabel('Fraction of Genes')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(os.path.join(RESULTS_DIR, 'smarcb1_binding_fractions.pdf'))
            plt.close()
            
            # Calculate and plot BM/BG enrichment ratios
            enrichment_ratios = []
            for category in df['category'].unique():
                bg_data = df[(df['condition'] == 'BG') & (df['category'] == category)]
                bm_data = df[(df['condition'] == 'BM') & (df['category'] == category)]
                
                if not bg_data.empty and not bm_data.empty:
                    bg_fraction = bg_data['fraction_bound'].values[0]
                    bm_fraction = bm_data['fraction_bound'].values[0]
                    
                    # Calculate ratio and propagate error
                    ratio = bm_fraction / bg_fraction if bg_fraction > 0 else 0
                    
                    # Calculate standard error of the ratio using error propagation
                    if bg_fraction > 0:
                        bg_n = bg_data['n_genes'].values[0]
                        bm_n = bm_data['n_genes'].values[0]
                        bg_se = np.sqrt(bg_fraction * (1 - bg_fraction) / bg_n)
                        bm_se = np.sqrt(bm_fraction * (1 - bm_fraction) / bm_n)
                        ratio_se = ratio * np.sqrt((bm_se/bm_fraction)**2 + (bg_se/bg_fraction)**2)
                    else:
                        ratio_se = 0
                    
                    enrichment_ratios.append({
                        'category': category,
                        'BM_BG_ratio': ratio,
                        'ratio_se': ratio_se
                    })
            
            ratio_df = pd.DataFrame(enrichment_ratios)
            plt.figure(figsize=(8, 6))
            sns.barplot(data=ratio_df, x='category', y='BM_BG_ratio', 
                       yerr=ratio_df['ratio_se'])
            plt.title('SMARCB1 Binding Enrichment (BM/BG Ratio)')
            plt.xlabel('Gene Category')
            plt.ylabel('BM/BG Ratio')
            plt.axhline(y=1, color='r', linestyle='--', alpha=0.5)
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(os.path.join(RESULTS_DIR, 'smarcb1_bm_bg_ratio.pdf'))
            plt.close()

    def plot_methylation_profiles(self, results: Dict):
        """Plot methylation profiles for different gene categories."""
        os.makedirs(RESULTS_DIR, exist_ok=True)
        
        # Plot average profiles
        for condition in results.keys():
            plt.figure(figsize=(10, 6))
            for category in results[condition].keys():
                profiles = results[condition][category]
                if len(profiles) > 0:
                    mean_profile = np.mean([p for p in profiles['profile'].tolist()], axis=0)
                    std_profile = np.std([p for p in profiles['profile'].tolist()], axis=0)
                    x = np.linspace(PROMOTER_WINDOW[0], PROMOTER_WINDOW[1], 
                                  len(mean_profile))
                    
                    plt.plot(x, mean_profile, label=f"{category} (n={len(set(profiles['gene']))})")
                    plt.fill_between(x, 
                                   mean_profile - std_profile,
                                   mean_profile + std_profile,
                                   alpha=0.2)
            
            plt.title(f'Methylation Profiles - {condition}')
            plt.xlabel('Distance from TSS (bp)')
            plt.ylabel('Methylation Level (RPM)')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(RESULTS_DIR, 
                                    f'methylation_profiles_{condition}.pdf'))
            plt.close()
            
        # Plot heatmaps with improved handling
        for condition in results.keys():
            for category in results[condition].keys():
                profiles = results[condition][category]
                if len(profiles) > 0:
                    # Convert profiles to matrix and handle missing values
                    profile_list = profiles['profile'].tolist()
                    profile_matrix = np.array([
                        [v if v is not None else 0 for v in p] 
                        for p in profile_list
                    ])
                    
                    # Normalize each row to improve visualization
                    row_max = np.max(profile_matrix, axis=1)
                    row_max[row_max == 0] = 1  # Avoid division by zero
                    profile_matrix_norm = profile_matrix / row_max[:, np.newaxis]
                    
                    # Sort rows by mean signal for better visualization
                    row_means = np.mean(profile_matrix_norm, axis=1)
                    sort_idx = np.argsort(row_means)
                    profile_matrix_norm = profile_matrix_norm[sort_idx]
                    
                    # Create figure with proper size and DPI
                    plt.figure(figsize=(12, len(profile_matrix_norm)/20 + 2))
                    sns.heatmap(profile_matrix_norm, 
                              cmap='YlOrRd',
                              xticklabels=False,
                              yticklabels=False,
                              rasterized=True)
                    plt.title(f'Methylation Heatmap - {condition} - {category}\n(n={len(profile_matrix_norm)} genes)')
                    plt.xlabel('Distance from TSS (bp)')
                    plt.ylabel('Genes')
                    
                    # Save with high DPI and compress
                    plt.savefig(os.path.join(RESULTS_DIR,
                                           f'methylation_heatmap_{condition}_{category}.pdf'),
                              dpi=300, bbox_inches='tight')
                    plt.close()

def main():
    """Run the methylation and SMARCB1 binding analysis."""
    print("Starting methylation cross-analysis...")
    
    # Initialize analyzer
    print("Initializing analysis...")
    analyzer = MethylationAnalyzer()
    
    # Analyze methylation patterns
    print("\nAnalyzing methylation patterns...")
    methylation_results, methylation_stats = analyzer.analyze_methylation_differences()
    
    # Analyze SMARCB1 binding
    print("\nAnalyzing SMARCB1 binding...")
    smarcb1_results, smarcb1_enrichment = analyzer.analyze_smarcb1_binding()
    
    # Generate plots
    print("\nGenerating plots...")
    analyzer.plot_methylation_profiles(methylation_results)
    analyzer.plot_smarcb1_enrichment(smarcb1_results, smarcb1_enrichment)
    
    # Save methylation statistics
    print("\nSaving statistics...")
    stats_rows = []
    for condition, condition_stats in methylation_stats.items():
        for category, stats in condition_stats.items():
            if isinstance(stats, dict) and 'mean_signal' in stats:
                stats_rows.append({
                    'condition': condition,
                    'category': category,
                    'mean_signal': stats['mean_signal'],
                    'std_signal': stats['std_signal'],
                    'n_genes': stats['n_genes']
                })
    
    stats_df = pd.DataFrame(stats_rows)
    print(f"Saving methylation statistics ({len(stats_df)} rows)")
    stats_df.to_csv(os.path.join(RESULTS_DIR, 'methylation_statistics.csv'), index=False)
    
    # Save SMARCB1 enrichment statistics
    enrichment_rows = []
    print("\nProcessing SMARCB1 enrichment statistics...")
    for condition, condition_enrichment in smarcb1_enrichment.items():
        print(f"Processing condition: {condition}")
        print(f"Categories: {list(condition_enrichment.keys())}")
        for category, stats in condition_enrichment.items():
            print(f"Adding stats for {condition} - {category}")
            enrichment_rows.append({
                'condition': condition,
                'category': category,
                'n_genes_with_peaks': stats['n_genes_with_peaks'],
                'total_genes': stats['total_genes'],
                'fraction_bound': stats['fraction_bound'],
                'mean_overlap': stats['mean_overlap']
            })
    
    enrichment_df = pd.DataFrame(enrichment_rows)
    print(f"Saving SMARCB1 enrichment statistics ({len(enrichment_df)} rows)")
    if len(enrichment_df) > 0:
        enrichment_df.to_csv(os.path.join(RESULTS_DIR, 'smarcb1_enrichment_statistics.csv'), index=False)
        print("Successfully saved SMARCB1 statistics")
    else:
        print("WARNING: No SMARCB1 enrichment data to save!")
    
    print(f"\nAnalysis complete. Results saved in: {RESULTS_DIR}")

if __name__ == "__main__":
    main()
