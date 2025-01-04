#%%
# Cell 1: Import required libraries
import pandas as pd
import numpy as np
import pyBigWig
import pysam
import pyranges as pr
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional
import logging
from concurrent.futures import ProcessPoolExecutor
from functools import partial, lru_cache
import os

# set working directory
os.chdir("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Methylation")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)
#%%
# Cell 2: Configuration and file paths
CONFIG = {
    'methylation_thresholds': {
        'hypo': 30,
        'hyper': 70
    },
    'quality_thresholds': {
        'coverage': 0.7,
        'min_cpgs': 5
    },
    'expression_thresholds': {
        'log2fc': 1,
        'padj': 0.05
    },
    'standard_chromosomes': [f'chr{i}' for i in range(1, 20)] + ['chrX', 'chrY']
}

# Update these paths according to your environment
PATHS = {
    'mecp2_dir': "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment_align2_005/NSC/mecp2_cpg_enrichment_parallel",
    'mecp2_file': "mecp2_cpg_enrichment_parallel.csv",
    'medip_dir': "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeDIP/output_old/bigwig",
    'cpg_islands_file': "../DATA/cpg_islands.bed",
    'rnaseq_files': {
        'NEU': '../iterative_alternative/DATA/DEA_NEU.csv',
        'NSC': '../iterative_alternative/DATA/DEA_NSC.csv'
    },
    'genome_fasta': "../DATA/mm10.fa",
    'gtf_file': "../DATA/gencode.vM10.annotation.gtf",
    'output_dir': "second_analysis"
}
#%%
# Cell 3: Utility functions for gene region extraction
def load_gene_annotations(gtf_file: str) -> pd.DataFrame:
    """
    Load gene annotations and extract promoter and gene body regions
    """
    genes = pr.read_gtf(gtf_file)
    genes = genes[genes.Feature == 'gene']
    
    # Create promoter regions (2kb upstream of TSS)
    promoters = genes.copy()
    promoters.Start = genes.Start - 2000
    promoters.End = genes.Start
    
    return {
        'genes': genes.as_df(),
        'promoters': promoters.as_df()
    }

def extract_methylation_regions(gene_regions: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Extract promoter and gene body regions for methylation analysis
    """
    promoters = gene_regions['promoters'].copy()
    genes = gene_regions['genes'].copy()
    
    # Add region type
    promoters['region_type'] = 'promoter'
    genes['region_type'] = 'gene_body'
    
    # Combine regions
    regions = pd.concat([promoters, genes])
    
    # Add necessary columns
    regions = regions[['Chromosome', 'Start', 'End', 'gene_name', 'region_type']]
    regions.columns = ['chr', 'start', 'end', 'gene_name', 'region_type']
    
    return regions

#%%
# Cell 4: Methylation analysis functions
@lru_cache(maxsize=1000)
def count_cpgs_in_sequence(sequence: str) -> int:
    """Count CpG dinucleotides in a DNA sequence"""
    return sequence.upper().count('CG')

def calculate_methylation_levels(bw_file: str, regions: pd.DataFrame, 
                               fasta_file: str) -> pd.DataFrame:
    """
    Calculate methylation levels for given regions
    """
    bw = pyBigWig.open(bw_file)
    fasta = pysam.FastaFile(fasta_file)
    
    results = regions.copy()
    
    # Calculate methylation levels
    methylation_values = []
    coverage_values = []
    cpg_counts = []
    
    for _, row in results.iterrows():
        # Get sequence and count CpGs
        seq = fasta.fetch(row['chr'], row['start'], row['end'])
        cpg_counts.append(count_cpgs_in_sequence(seq))
        
        # Get methylation values
        vals = bw.values(row['chr'], row['start'], row['end'])
        if vals and any(v is not None for v in vals):
            valid_vals = [v for v in vals if v is not None]
            methylation_values.append(np.mean(valid_vals))
            coverage_values.append(len(valid_vals) / len(vals))
        else:
            methylation_values.append(0)
            coverage_values.append(0)
    
    results['methylation'] = methylation_values
    results['coverage'] = coverage_values
    results['cpg_count'] = cpg_counts
    
    return results
#%%
# Cell 5: Expression analysis functions
def analyze_expression_patterns(dea_file: str, methylation_data: pd.DataFrame) -> pd.DataFrame:
    """
    Analyze expression patterns in relation to methylation and MeCP2 binding
    """
    # Load differential expression data
    dea = pd.read_csv(dea_file)
    
    # Classify genes based on expression
    dea['expression_status'] = 'unchanged'
    dea.loc[(dea['log2FoldChange'] > CONFIG['expression_thresholds']['log2fc']) & 
            (dea['padj'] < CONFIG['expression_thresholds']['padj']), 'expression_status'] = 'upregulated'
    dea.loc[(dea['log2FoldChange'] < -CONFIG['expression_thresholds']['log2fc']) & 
            (dea['padj'] < CONFIG['expression_thresholds']['padj']), 'expression_status'] = 'downregulated'
    
    # Merge with methylation data
    merged = methylation_data.merge(dea[['gene_name', 'log2FoldChange', 'padj', 'expression_status']], 
                                  on='gene_name', how='inner')
    
    return merged

#%%
# Cell 6: Visualization functions
def plot_methylation_expression_patterns(data: pd.DataFrame, output_dir: str):
    """
    Create visualizations for methylation and expression patterns
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Methylation levels in different regions by expression status
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=data, x='expression_status', y='methylation', hue='region_type')
    plt.title('Methylation Levels by Expression Status and Region Type')
    plt.savefig(f'{output_dir}/methylation_by_expression.pdf')
    plt.close()
    
    # 2. CpG density vs methylation for different expression categories
    plt.figure(figsize=(12, 6))
    for status in data['expression_status'].unique():
        subset = data[data['expression_status'] == status]
        plt.scatter(subset['cpg_count'], subset['methylation'], 
                   label=status, alpha=0.5)
    plt.xlabel('CpG Count')
    plt.ylabel('Methylation Level')
    plt.title('CpG Density vs Methylation by Expression Status')
    plt.legend()
    plt.savefig(f'{output_dir}/cpg_vs_methylation.pdf')
    plt.close()

#%%
# Cell 7: Statistical analysis functions
def perform_statistical_analysis(data: pd.DataFrame) -> Dict:
    """
    Perform statistical tests to compare patterns between groups
    """
    results = {}
    
    # Compare methylation levels between up/down regulated genes
    for region in ['promoter', 'gene_body']:
        region_data = data[data['region_type'] == region]
        
        up_meth = region_data[region_data['expression_status'] == 'upregulated']['methylation']
        down_meth = region_data[region_data['expression_status'] == 'downregulated']['methylation']
        
        stat, pval = stats.mannwhitneyu(up_meth, down_meth)
        
        results[f'{region}_methylation_comparison'] = {
            'statistic': stat,
            'pvalue': pval,
            'up_mean': up_meth.mean(),
            'down_mean': down_meth.mean()
        }
    
    return results

#%%
##########################################################################################
# 1. Load gene annotations
gene_regions = load_gene_annotations(PATHS['gtf_file'])
regions = extract_methylation_regions(gene_regions)

#%%
# 2. Calculate methylation levels
methylation_data = calculate_methylation_levels(
    os.path.join(PATHS['medip_dir'], 'methylation.bw'),
    regions,
    PATHS['genome_fasta']
)

#%%
# 3. Analyze expression patterns
results = {}
for cell_type, dea_file in PATHS['rnaseq_files'].items():
    integrated_data = analyze_expression_patterns(dea_file, methylation_data)
    
    # 4. Create visualizations
    plot_methylation_expression_patterns(
        integrated_data,
        os.path.join(PATHS['output_dir'], cell_type)
    )
    
    # 5. Perform statistical analysis
    stats_results = perform_statistical_analysis(integrated_data)
    results[cell_type] = stats_results

#%%
results