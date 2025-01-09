import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Any, Union
import logging
from config import CONFIG, PATHS, logger

def create_tss_analysis_plots(results: Dict[str, Dict], cell_type: str, output_dir: str):
    """Create visualization plots for TSS binding analysis with CpG density consideration
    
    Args:
        results: Dictionary containing analysis results with structure:
            {category: {
                reg_status: {
                    'count': int,
                    'high_quality_count': int,
                    'promoter_metrics': {
                        'methylation': {'mean': float, 'std': float, 'median': float},
                        'cpg_density': {'mean': float, 'std': float, 'median': float},
                        'cpg_count': {'total': int, 'mean': float}
                    },
                    'gene_body_metrics': {
                        'methylation': {'mean': float, 'std': float, 'median': float},
                        'cpg_density': {'mean': float, 'std': float, 'median': float},
                        'cpg_count': {'total': int, 'mean': float}
                    }
                }
            }}
        cell_type: Cell type being analyzed
        output_dir: Directory to save plots
    """
    # Create plot directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Prepare data for plotting
    plot_data = {
        'Category': [],
        'Regulation': [],
        'Region': [],
        'Methylation': [],
        'StdDev': [],
        'CpG_Density': [],
        'Sample_Size': []
    }
    
    # Extract methylation data for each category, regulation status and genomic region
    for category, stats in results.items():
        for reg_status, reg_stats in stats.items():
            # Add promoter data
            plot_data['Category'].append(category)
            plot_data['Regulation'].append(reg_status)
            plot_data['Region'].append('Promoter')
            plot_data['Methylation'].append(
                reg_stats['promoter_metrics']['methylation']['mean']
            )
            plot_data['StdDev'].append(
                reg_stats['promoter_metrics']['methylation']['std']
            )
            plot_data['CpG_Density'].append(
                reg_stats['promoter_metrics']['cpg_density']['mean']
            )
            plot_data['Sample_Size'].append(reg_stats['high_quality_count'])
            
            # Add gene body data
            plot_data['Category'].append(category)
            plot_data['Regulation'].append(reg_status)
            plot_data['Region'].append('Gene Body')
            plot_data['Methylation'].append(
                reg_stats['gene_body_metrics']['methylation']['mean']
            )
            plot_data['StdDev'].append(
                reg_stats['gene_body_metrics']['methylation']['std']
            )
            plot_data['CpG_Density'].append(
                reg_stats['gene_body_metrics']['cpg_density']['mean']
            )
            plot_data['Sample_Size'].append(reg_stats['high_quality_count'])
    
    # Convert to DataFrame for easier plotting
    df = pd.DataFrame(plot_data)
    
    # Create figure with two subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # Create promoter methylation plot
    promoter_data = df[df['Region'] == 'Promoter']
    sns.barplot(data=promoter_data, x='Category', y='Methylation',
                hue='Regulation', ax=ax1)
    ax1.set_title(f'{cell_type} - Promoter Methylation')
    ax1.set_xlabel('Binding Category')
    ax1.set_ylabel('Methylation Level (%)')
    ax1.tick_params(axis='x', rotation=45)
    
    # Create gene body methylation plot
    gene_body_data = df[df['Region'] == 'Gene Body']
    sns.barplot(data=gene_body_data, x='Category', y='Methylation',
                hue='Regulation', ax=ax2)
    ax2.set_title(f'{cell_type} - Gene Body Methylation')
    ax2.set_xlabel('Binding Category')
    ax2.set_ylabel('Methylation Level (%)')
    ax2.tick_params(axis='x', rotation=45)
    
    # Add sample size annotations
    for ax, data in [(ax1, promoter_data), (ax2, gene_body_data)]:
        for i, category in enumerate(data['Category'].unique()):
            category_data = data[data['Category'] == category]
            y_max = category_data['Methylation'].max()
            ax.text(i, y_max + 2, f'n={category_data["Sample_Size"].iloc[0]}',
                   ha='center', va='bottom')
    
    plt.suptitle(f'{cell_type} - Methylation Analysis by Binding Category', y=1.05)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{cell_type}_methylation_analysis.pdf'),
                bbox_inches='tight')
    plt.close()
    
    # Create CpG density plot
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=df, x='Category', y='CpG_Density', hue='Region')
    plt.title(f'{cell_type} - CpG Density by Category and Region')
    plt.xlabel('Binding Category')
    plt.ylabel('CpG Density (CpGs/kb)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{cell_type}_cpg_density.pdf'))
    plt.close()

def save_tss_analysis_results(results: Dict[str, Dict], cell_type: str, output_dir: str):
    """Save detailed statistics from TSS binding analysis
    
    Args:
        results: Dictionary containing analysis results with structure:
            {category: {
                reg_status: {
                    'count': int,
                    'high_quality_count': int,
                    'promoter_metrics': {
                        'methylation': {'mean': float, 'std': float, 'median': float},
                        'cpg_density': {'mean': float, 'std': float, 'median': float},
                        'cpg_count': {'total': int, 'mean': float}
                    },
                    'gene_body_metrics': {
                        'methylation': {'mean': float, 'std': float, 'median': float},
                        'cpg_density': {'mean': float, 'std': float, 'median': float},
                        'cpg_count': {'total': int, 'mean': float}
                    }
                }
            }}
        cell_type: Cell type being analyzed
        output_dir: Directory to save results
    """
    # Create output file
    stats_file = os.path.join(output_dir, f'{cell_type}_tss_analysis_stats.txt')
    
    with open(stats_file, 'w') as f:
        f.write(f"TSS Binding Analysis Results - {cell_type}\n")
        f.write("="*50 + "\n\n")
        
        for category, stats in results.items():
            f.write(f"\n{category.upper()}\n")
            f.write("-"*30 + "\n")
            
            total_genes = sum(stats.get(status, {}).get('count', 0) 
                            for status in ['not_deregulated', 'upregulated', 'downregulated'])
            
            f.write(f"Total genes: {total_genes}\n")
            f.write(f"High quality genes: {sum(stats.get(status, {}).get('high_quality_count', 0) for status in stats)}\n\n")
            
            for status, measurements in stats.items():
                f.write(f"{status}:\n")
                f.write(f"  Count: {measurements['count']}\n")
                f.write(f"  High quality count: {measurements['high_quality_count']}\n")
                f.write(f"  Percentage: {(measurements['count']/total_genes)*100:.2f}%\n")
                
                # Promoter metrics
                f.write("\n  Promoter Metrics:\n")
                f.write("    Methylation:\n")
                f.write(f"      Mean ± SD: {measurements['promoter_metrics']['methylation']['mean']:.2f} ± "
                       f"{measurements['promoter_metrics']['methylation']['std']:.2f}\n")
                f.write(f"      Median: {measurements['promoter_metrics']['methylation']['median']:.2f}\n")
                f.write("    CpG Density:\n")
                f.write(f"      Mean: {measurements['promoter_metrics']['cpg_density']['mean']:.2f} CpGs/kb\n")
                f.write(f"      Total CpGs: {measurements['promoter_metrics']['cpg_count']['total']}\n")
                
                # Gene body metrics
                f.write("\n  Gene Body Metrics:\n")
                f.write("    Methylation:\n")
                f.write(f"      Mean ± SD: {measurements['gene_body_metrics']['methylation']['mean']:.2f} ± "
                       f"{measurements['gene_body_metrics']['methylation']['std']:.2f}\n")
                f.write(f"      Median: {measurements['gene_body_metrics']['methylation']['median']:.2f}\n")
                f.write("    CpG Density:\n")
                f.write(f"      Mean: {measurements['gene_body_metrics']['cpg_density']['mean']:.2f} CpGs/kb\n")
                f.write(f"      Total CpGs: {measurements['gene_body_metrics']['cpg_count']['total']}\n")
                f.write("\n")


def analyze_tss_binding_patterns(df: pd.DataFrame, cell_type: str, output_dir: str):
    """Analyze MeCP2 binding patterns around transcription start sites (TSS) and their relationship with gene regulation.
    
    This function performs a comprehensive analysis of how MeCP2 binds around gene transcription start sites
    and how this binding correlates with gene regulation. It:

    1. Categorizes genes into binding groups based on MeCP2 ChIP-seq data:
       - exo_enriched: Genes bound by MeCP2 in both exogenous and endogenous conditions
       - exo_only: Genes bound only in exogenous MeCP2 conditions
       - endo_only: Genes bound only in endogenous MeCP2 conditions  
       - non_enriched: Genes not bound by MeCP2

    2. For each binding category, further subdivides genes by expression status:
       - upregulated: Genes with increased expression
       - downregulated: Genes with decreased expression
       - not_deregulated: Genes with unchanged expression

    3. Calculates detailed statistics for each subgroup including:
       - Gene counts and percentages
       - Promoter methylation levels (mean, median, std)
       - Gene body methylation levels (mean, median, std)

    4. Creates visualization plots showing:
       - Distribution of binding patterns
       - Correlation between binding and regulation
       - Methylation patterns in different groups

    5. Saves results to disk:
       - Gene lists as CSV files
       - Statistical metrics in text files
       - Visualization plots as PDFs

    The results provide insights into how MeCP2 binding patterns around TSS may influence
    gene regulation and how this relates to DNA methylation levels.

    Args:
        df: DataFrame containing gene info, binding data and expression status
        cell_type: Cell type being analyzed (e.g. 'NSC', 'NEU') 
        output_dir: Base output directory for saving results
        
    Returns:
        Dictionary containing statistical results for each binding category and regulation status
    """
    
    # Create output directories
    tss_dir = os.path.join(output_dir, 'tss_analysis', cell_type)
    os.makedirs(tss_dir, exist_ok=True)
    
    # Define binding categories
    binding_categories = {
        'exo_enriched': df[df['binding_type'].isin(['exo', 'both'])],
        'exo_only': df[df['binding_type'] == 'exo'],
        'endo_only': df[df['binding_type'] == 'endo'],
        'non_enriched': df[~df['mecp2_bound']]
    }
    
    # Analyze each category
    results = {}
    for category, category_df in binding_categories.items():
        # Split by regulation status
        regulation_groups = {
            'not_deregulated': category_df[category_df['expression_status'] == 'unchanged'],
            'upregulated': category_df[category_df['expression_status'] == 'upregulated'],
            'downregulated': category_df[category_df['expression_status'] == 'downregulated']
        }
        
        # Calculate statistics for each group
        group_stats = {}
        for reg_status, group_df in regulation_groups.items():
            if len(group_df) > 0:
                group_stats[reg_status] = {
                    'count': len(group_df),
                    'promoter_methylation': {
                        'mean': group_df['promoter_methylation'].mean(),
                        'std': group_df['promoter_methylation'].std(),
                        'median': group_df['promoter_methylation'].median()
                    },
                    'gene_body_methylation': {
                        'mean': group_df['gene_body_methylation'].mean(),
                        'std': group_df['gene_body_methylation'].std(),
                        'median': group_df['gene_body_methylation'].median()
                    }
                }
                
                # Save gene lists
                output_file = os.path.join(tss_dir, f'{category}_{reg_status}_genes.csv')
                group_df.to_csv(output_file, index=False)
        
        results[category] = group_stats
    
    # Create visualization
    create_tss_analysis_plots(results, cell_type, tss_dir)
    
    # Save detailed statistics
    save_tss_analysis_results(results, cell_type, tss_dir)
    
    return results
