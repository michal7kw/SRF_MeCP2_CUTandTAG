import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools
import os
import pysam

# Define genome-wide constants
GENOME_SIZE = 2.7e9  # Mouse genome size
PROMOTER_WINDOW = 2000  # Define promoter region as ±2kb from TSS

# Set working directory
os.chdir('/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline')


# Add function to calculate sequencing depth
def calculate_sequencing_depth(bam_file):
    """Calculate total mapped reads from BAM file"""
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        return bam.count()

def load_data():
    # Load DEA results
    dea_nsc = pd.read_csv("DATA/DEA_NSC.csv")
    print("\nDEA file columns:", dea_nsc.columns.tolist())
    
    # Verify required columns exist
    required_columns = ['gene', 'log2FoldChange', 'padj']
    missing_columns = [col for col in required_columns if col not in dea_nsc.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in DEA file: {missing_columns}")
    
    # Load peak data and calculate sequencing depths
    exo_samples = ['NSCv1', 'NSCv2', 'NSCv3']
    endo_samples = ['NSCM1', 'NSCM2', 'NSCM3']
    
    peaks_exo = {}
    peaks_endo = {}
    depths_exo = {}
    depths_endo = {}
    
    # Load exogenous samples and their depths
    for sample in exo_samples:
        peaks_exo[sample] = pd.read_csv(f"results/peaks/{sample}_peaks.narrowPeak", 
                                      sep='\t', header=None,
                                      names=['chr', 'start', 'end', 'name', 'score', 
                                            'strand', 'signalValue', 'pValue', 
                                            'qValue', 'peak'])
        # Calculate sequencing depth
        depths_exo[sample] = calculate_sequencing_depth(f"results/aligned/{sample}.bam")
        # Normalize signal values by sequencing depth
        peaks_exo[sample]['signalValue'] = peaks_exo[sample]['signalValue'] * (1e6 / depths_exo[sample])
    
    # Load endogenous samples and their depths
    for sample in endo_samples:
        peaks_endo[sample] = pd.read_csv(f"results/peaks/{sample}_peaks.narrowPeak",
                                       sep='\t', header=None,
                                       names=['chr', 'start', 'end', 'name', 'score',
                                            'strand', 'signalValue', 'pValue',
                                            'qValue', 'peak'])
        # Calculate sequencing depth
        depths_endo[sample] = calculate_sequencing_depth(f"results/aligned/{sample}.bam")
        # Normalize signal values by sequencing depth
        peaks_endo[sample]['signalValue'] = peaks_endo[sample]['signalValue'] * (1e6 / depths_endo[sample])
    
    # Load gene annotations
    gene_annotations, name_to_info = load_gene_annotations()
    
    return dea_nsc, peaks_exo, peaks_endo, gene_annotations, name_to_info

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
    gtf_file = "DATA/gencode.vM10.annotation.gtf"
    
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
        overlapping = chr_peaks[
            (chr_peaks['start'] <= promoter_end) & 
            (chr_peaks['end'] >= promoter_start)
        ]
        
        if not overlapping.empty:
            overlapping['sample'] = sample
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
                (-np.log10(exo['qValue'].clip(1e-10)))
            ) / max(np.sum(
                (endo['end'] - endo['start']) * 
                endo['signalValue'] * 
                (-np.log10(endo['qValue'].clip(1e-10)))
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
        df_sorted.to_csv(f'results/enrichment_{method_name}_NSC.csv', index=False)
    
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
    summary_df.to_csv('results/enrichment_summary_NSC.csv')
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
    plt.savefig('results/peak_width_distributions.pdf')
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
    plt.savefig('results/peak_width_distributions_detailed.pdf')
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
            width_ratios[gene] = mean_width_exo / mean_width_endo
    
    # Create scatter plots
    for i, method in enumerate(['width_weighted', 'coverage_score', 'area_integration']):
        plt.subplot(1, 3, i+1)
        
        x = [width_ratios[gene] for gene in results[method]['gene'] if gene in width_ratios]
        y = [score for gene, score in zip(results[method]['gene'], 
                                        results[method]['enrichment_score']) 
             if gene in width_ratios]
        
        plt.scatter(x, y, alpha=0.5)
        plt.xlabel('Peak Width Ratio (Exo/Endo)')
        plt.ylabel(f'{method} Enrichment Score')
        plt.title(f'{method} vs Width Ratio')
        plt.yscale('log')
        plt.xscale('log')
    
    plt.tight_layout()
    plt.savefig('results/width_vs_enrichment_NSC.pdf')
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
    summary_df.to_csv('results/peak_distribution_summary_NSC.csv')
    return summary_df

# if __name__ == "__main__":
# Create output directory if it doesn't exist
os.makedirs('results', exist_ok=True)

# Load data
dea, peaks_exo, peaks_endo, gene_annotations, name_to_info = load_data()

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
    df.to_csv(f'results/enrichment_{method_name}_NSC.csv', index=False) 

# Plot width vs enrichment
plot_width_vs_enrichment(results, peaks_exo, peaks_endo, gene_annotations, name_to_info)

# Summarize peak distribution
peak_distribution = summarize_peak_distribution(results)
print("\nPeak Distribution Summary:")
print(peak_distribution)
