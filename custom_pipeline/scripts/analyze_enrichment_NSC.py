import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools
import os

# Define genome-wide constants
GENOME_SIZE = 2.7e9  # Mouse genome size
PROMOTER_WINDOW = 2000  # Define promoter region as ±2kb from TSS

def load_data():
    # Load DEA results
    dea_nsc = pd.read_csv("DATA/DEA_NSC.csv")
    
    # Load peak data for exogenous and endogenous samples
    exo_samples = ['NSCv1', 'NSCv2', 'NSCv3']
    endo_samples = ['NSCM1', 'NSCM2', 'NSCM3']
    
    peaks_exo = {}
    peaks_endo = {}
    
    for sample in exo_samples:
        peaks_exo[sample] = pd.read_csv(f"results/peaks/{sample}_peaks.narrowPeak", 
                                      sep='\t', header=None,
                                      names=['chr', 'start', 'end', 'name', 'score', 
                                            'strand', 'signalValue', 'pValue', 
                                            'qValue', 'peak'])
    
    for sample in endo_samples:
        peaks_endo[sample] = pd.read_csv(f"results/peaks/{sample}_peaks.narrowPeak",
                                       sep='\t', header=None,
                                       names=['chr', 'start', 'end', 'name', 'score',
                                            'strand', 'signalValue', 'pValue',
                                            'qValue', 'peak'])
    
    # Load gene annotations
    gene_annotations = load_gene_annotations()
    
    return dea_nsc, peaks_exo, peaks_endo, gene_annotations

def load_gene_annotations():
    """Load gene annotations from GTF file"""
    gtf_file = "DATA/gencode.vM10.annotation.gtf" 
    gene_annotations = pd.read_csv(gtf_file, sep='\t', comment='#',
                                 names=['chr', 'source', 'feature', 'start', 'end',
                                       'score', 'strand', 'frame', 'attributes'])
    # Filter for genes only
    gene_annotations = gene_annotations[gene_annotations['feature'] == 'gene']
    return gene_annotations

def get_peaks_near_gene(gene, peaks_dict, gene_annotations, window=PROMOTER_WINDOW):
    """Get peaks near a gene's promoter region"""
    # Get gene coordinates
    gene_info = gene_annotations[gene_annotations['attributes'].str.contains(gene)]
    if len(gene_info) == 0:
        return pd.DataFrame()  # Return empty if gene not found
    
    # Create promoter region (TSS ± window)
    if gene_info.iloc[0]['strand'] == '+':
        promoter_start = gene_info.iloc[0]['start'] - window
        promoter_end = gene_info.iloc[0]['start'] + window
    else:
        promoter_start = gene_info.iloc[0]['end'] - window
        promoter_end = gene_info.iloc[0]['end'] + window
    
    # Collect peaks from all samples that overlap with promoter
    all_peaks = []
    for sample, peaks in peaks_dict.items():
        sample_peaks = peaks[
            (peaks['chr'] == gene_info.iloc[0]['chr']) &
            (peaks['start'] <= promoter_end) &
            (peaks['end'] >= promoter_start)
        ]
        all_peaks.append(sample_peaks)
    
    # Combine peaks from all samples
    if all_peaks:
        return pd.concat(all_peaks)
    return pd.DataFrame()

def calculate_total_genome_peaks(peaks_exo, peaks_endo):
    """Calculate total number of peaks across all samples"""
    total_exo = sum(len(df) for df in peaks_exo.values())
    total_endo = sum(len(df) for df in peaks_endo.values())
    return total_exo + total_endo

def define_enrichment_methods(total_genome_peaks):
    """Define different methods to calculate enrichment"""
    
    methods = {
        # Method 1: Signal Value Ratio
        'signal_ratio': lambda exo, endo: np.mean(exo['signalValue']) / np.mean(endo['signalValue']) if len(exo) > 0 and len(endo) > 0 else np.nan,
        
        # Method 2: Peak Count Ratio
        'peak_count': lambda exo, endo: len(exo) / len(endo) if len(endo) > 0 else np.nan,
        
        # Method 3: Combined Score
        'combined_score': lambda exo, endo: (
            np.mean(exo['signalValue']) * len(exo)) / (np.mean(endo['signalValue']) * len(endo)
        ) if len(exo) > 0 and len(endo) > 0 else np.nan,
        
        # Method 4: Statistical Enrichment (Fisher's Exact Test)
        'statistical': lambda exo, endo: stats.fisher_exact([
            [len(exo), len(endo)],
            [total_genome_peaks - len(exo), total_genome_peaks - len(endo)]
        ])[1] if len(exo) > 0 and len(endo) > 0 else np.nan  # Returns p-value
    }
    
    return methods

def analyze_enrichment():
    # Load data
    dea, peaks_exo, peaks_endo, gene_annotations = load_data()
    
    # Calculate total genome peaks
    total_genome_peaks = calculate_total_genome_peaks(peaks_exo, peaks_endo)
    
    # Define up-regulated genes (log2FC > 1 and padj < 0.05)
    upreg_genes = dea[(dea['log2FoldChange'] > 1) & (dea['padj'] < 0.05)]['gene'].tolist()
    
    # Calculate enrichment using different methods
    methods = define_enrichment_methods(total_genome_peaks)
    results = {}
    
    for method_name, method_func in methods.items():
        enrichment_scores = []
        for gene in upreg_genes:
            # Get peaks near gene in exogenous and endogenous samples
            exo_peaks = get_peaks_near_gene(gene, peaks_exo, gene_annotations)
            endo_peaks = get_peaks_near_gene(gene, peaks_endo, gene_annotations)
            
            if not (exo_peaks.empty and endo_peaks.empty):
                score = method_func(exo_peaks, endo_peaks)
                if not np.isnan(score):
                    enrichment_scores.append((gene, score))
        
        results[method_name] = pd.DataFrame(enrichment_scores, 
                                          columns=['gene', 'enrichment_score'])
    
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
    plt.savefig('results/enrichment_analysis_NSC.pdf')
    plt.close()

def summarize_results(results):
    """Create summary statistics for each enrichment method"""
    summary = {}
    for method_name, df in results.items():
        summary[method_name] = {
            'median': df['enrichment_score'].median(),
            'mean': df['enrichment_score'].mean(),
            'std': df['enrichment_score'].std(),
            'num_genes': len(df),
            'significant_genes': len(df[df['enrichment_score'] > 1])
        }
    
    summary_df = pd.DataFrame(summary).T
    summary_df.to_csv('results/enrichment_summary_NSC.csv')

if __name__ == "__main__":
    # Create output directory if it doesn't exist
    os.makedirs('results', exist_ok=True)
    
    # Run analysis
    results = analyze_enrichment()
    
    # Create visualizations
    plot_enrichment(results)
    
    # Generate summary statistics
    summarize_results(results)
    
    # Save detailed results to CSV
    for method_name, df in results.items():
        df.to_csv(f'results/enrichment_{method_name}_NSC.csv', index=False) 
