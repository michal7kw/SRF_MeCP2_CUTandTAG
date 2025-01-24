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
                # Try reading with header
                df = pd.read_csv(filepath)
                if 'GENE' in df.columns:
                    genes = df['GENE'].tolist()
                else:
                    # If no header or different column name, assume first column contains genes
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
        
        # Filter genes to those we have coordinates for
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
                                    'profile': values
                                })
                        except Exception as e:
                            print(f"Error processing gene {gene}: {str(e)}")
                            continue
            except Exception as e:
                print(f"Error opening bigwig file {bw_file}: {str(e)}")
                continue
                            
        return pd.DataFrame(profiles)

    def analyze_methylation_differences(self):
        """Analyze methylation differences between gene categories."""
        results = {}
        
        for condition in MEDIP_BIGWIG.keys():
            print(f"\nAnalyzing methylation for condition: {condition}")
            methylation_data = {}
            for category, genes in self.gene_lists.items():
                print(f"\nProcessing {category} genes...")
                profiles = self.calculate_methylation_profile(
                    MEDIP_BIGWIG[condition], genes
                )
                methylation_data[category] = profiles
            
            results[condition] = methylation_data
            
        return results

    def analyze_smarcb1_binding(self):
        """Analyze SMARCB1 binding in different gene categories."""
        results = {}
        
        for condition in ['BG', 'BM']:
            print(f"\nAnalyzing SMARCB1 binding for condition: {condition}")
            binding_data = {}
            for category, genes in self.gene_lists.items():
                print(f"\nProcessing {category} genes...")
                gene_coords = self.gene_coords.loc[
                    self.gene_coords.index.isin(genes)
                ]
                promoter_coords = get_promoter_coords(
                    gene_coords, PROMOTER_WINDOW
                )
                overlaps = calculate_peak_overlaps(
                    SMARCB1_PEAKS[condition], promoter_coords
                )
                binding_data[category] = overlaps
            
            results[condition] = binding_data
            
        return results

    def plot_methylation_profiles(self, results: Dict):
        """Plot methylation profiles for different gene categories."""
        os.makedirs(RESULTS_DIR, exist_ok=True)
        
        for condition in results.keys():
            plt.figure(figsize=(10, 6))
            for category in results[condition].keys():
                profiles = results[condition][category]
                if len(profiles) > 0:
                    mean_profile = np.mean([p for p in profiles['profile'].tolist()], axis=0)
                    std_profile = np.std([p for p in profiles['profile'].tolist()], axis=0)
                    x = np.linspace(PROMOTER_WINDOW[0], PROMOTER_WINDOW[1], 
                                  len(mean_profile))
                    
                    plt.plot(x, mean_profile, label=f"{category} (n={len(profiles)})")
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

    def plot_smarcb1_enrichment(self, results: Dict):
        """Plot SMARCB1 binding enrichment for different gene categories."""
        os.makedirs(RESULTS_DIR, exist_ok=True)
        
        # Prepare data for plotting
        plot_data = []
        for condition in results.keys():
            for category in results[condition].keys():
                overlaps = results[condition][category]
                if len(overlaps) > 0:
                    plot_data.append({
                        'condition': condition,
                        'category': category,
                        'n_peaks': len(overlaps),
                        'n_genes': len(set(overlaps['gene'])),
                        'mean_overlap': overlaps['overlap_ratio'].mean()
                    })
        
        if plot_data:
            df = pd.DataFrame(plot_data)
            
            # Plot number of genes with peaks
            plt.figure(figsize=(10, 6))
            sns.barplot(data=df, x='category', y='n_genes', hue='condition')
            plt.title('Number of Genes with SMARCB1 Peaks')
            plt.xlabel('Gene Category')
            plt.ylabel('Number of Genes')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(os.path.join(RESULTS_DIR, 'smarcb1_peak_counts.pdf'))
            plt.close()
            
            # Plot mean overlap ratio
            plt.figure(figsize=(10, 6))
            sns.barplot(data=df, x='category', y='mean_overlap', hue='condition')
            plt.title('Mean SMARCB1 Peak Overlap Ratio')
            plt.xlabel('Gene Category')
            plt.ylabel('Mean Overlap Ratio')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(os.path.join(RESULTS_DIR, 'smarcb1_overlap_ratio.pdf'))
            plt.close()

def main():
    print("Initializing analysis...")
    analyzer = MethylationAnalyzer()
    
    # Analyze methylation
    print("\nAnalyzing methylation patterns...")
    methylation_results = analyzer.analyze_methylation_differences()
    analyzer.plot_methylation_profiles(methylation_results)
    
    # Analyze SMARCB1 binding
    print("\nAnalyzing SMARCB1 binding...")
    smarcb1_results = analyzer.analyze_smarcb1_binding()
    analyzer.plot_smarcb1_enrichment(smarcb1_results)
    
    print("\nAnalysis complete. Results saved in:", RESULTS_DIR)

if __name__ == "__main__":
    main()
