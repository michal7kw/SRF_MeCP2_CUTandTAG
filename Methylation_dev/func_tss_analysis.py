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
    """Analyze MeCP2 binding patterns around TSS with CpG density considerations
    
    Args:
        df: DataFrame containing:
            - gene info and binding data
            - methylation levels weighted by CpG density
            - CpG density metrics
            - quality control flags
        cell_type: Cell type being analyzed (e.g. 'NSC', 'NEU') 
        output_dir: Base output directory for saving results
    """
    
    # Create output directories
    tss_dir = os.path.join(output_dir, 'tss_analysis', cell_type)
    os.makedirs(tss_dir, exist_ok=True)
    
    # Ensure required columns exist
    required_columns = ['log2FoldChange', 'padj', 'mecp2_bound']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        logger.error(f"Missing required columns: {missing_columns}")
        return None
    
    # Create a copy to avoid modifying original
    df = df.copy()
    
    # Add expression_status column
    df['expression_status'] = 'unchanged'
    df.loc[(df['log2FoldChange'] > 1) & (df['padj'] < 0.05), 'expression_status'] = 'upregulated'
    df.loc[(df['log2FoldChange'] < -1) & (df['padj'] < 0.05), 'expression_status'] = 'downregulated'
    
    # Log expression status distribution
    logger.info("\nExpression status distribution:")
    logger.info(df['expression_status'].value_counts().to_string())
    
    # Handle mecp2_bound column more robustly
    if 'mecp2_bound' in df.columns:
        # First, check the current values
        logger.debug(f"mecp2_bound unique values before conversion: {df['mecp2_bound'].unique()}")
        logger.debug(f"mecp2_bound dtype before conversion: {df['mecp2_bound'].dtype}")
        
        # Convert to boolean more carefully
        try:
            # Handle NaN values first
            df['mecp2_bound'] = df['mecp2_bound'].fillna(False)
            
            # Convert numeric values to boolean
            if df['mecp2_bound'].dtype in ['float64', 'float32', 'int64', 'int32']:
                df['mecp2_bound'] = df['mecp2_bound'] > 0
            else:
                # Try string conversion for any other type
                df['mecp2_bound'] = df['mecp2_bound'].map({'True': True, 'False': False, '1': True, '0': False})
                df['mecp2_bound'] = df['mecp2_bound'].fillna(False)
        
        except Exception as e:
            logger.error(f"Error converting mecp2_bound to boolean: {str(e)}")
            logger.error("Setting all values to False")
            df['mecp2_bound'] = False
        
        logger.debug(f"mecp2_bound unique values after conversion: {df['mecp2_bound'].unique()}")
        logger.debug(f"mecp2_bound dtype after conversion: {df['mecp2_bound'].dtype}")
    else:
        logger.error("mecp2_bound column not found in DataFrame")
        return None
    
    # Define binding categories
    binding_categories = {
        'exo_enriched': df[df['binding_type'].isin(['exo', 'both'])],
        'exo_only': df[df['binding_type'] == 'exo'],
        'endo_only': df[df['binding_type'] == 'endo'],
        'non_enriched': df[df['mecp2_bound'] == False]  # Changed from ~ operator to == False
    }
    
    # Log category sizes
    for category, subset in binding_categories.items():
        logger.info(f"{category}: {len(subset)} genes")
    
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
                # Filter out low quality data
                high_quality_df = group_df[
                    ~group_df['low_cpg_promoter'] & 
                    ~group_df['high_meth_promoter']
                ]
                
                group_stats[reg_status] = {
                    'count': len(group_df),
                    'high_quality_count': len(high_quality_df),
                    'promoter_metrics': {
                        'methylation': {
                            'mean': high_quality_df['promoter_methylation'].mean(),
                            'std': high_quality_df['promoter_methylation'].std(),
                            'median': high_quality_df['promoter_methylation'].median()
                        },
                        'cpg_density': {
                            'mean': high_quality_df['promoter_cpg_density'].mean(),
                            'std': high_quality_df['promoter_cpg_density'].std(),
                            'median': high_quality_df['promoter_cpg_density'].median()
                        },
                        'cpg_count': {
                            'total': high_quality_df['promoter_cpg_count'].sum(),
                            'mean': high_quality_df['promoter_cpg_count'].mean()
                        }
                    },
                    'gene_body_metrics': {
                        'methylation': {
                            'mean': high_quality_df['gene_body_methylation'].mean(),
                            'std': high_quality_df['gene_body_methylation'].std(),
                            'median': high_quality_df['gene_body_methylation'].median()
                        },
                        'cpg_density': {
                            'mean': high_quality_df['gene_body_cpg_density'].mean(),
                            'std': high_quality_df['gene_body_cpg_density'].std(),
                            'median': high_quality_df['gene_body_cpg_density'].median()
                        },
                        'cpg_count': {
                            'total': high_quality_df['gene_body_cpg_count'].sum(),
                            'mean': high_quality_df['gene_body_cpg_count'].mean()
                        }
                    }
                }
                
                # Save gene lists with quality metrics
                output_file = os.path.join(tss_dir, f'{category}_{reg_status}_genes.csv')
                group_df.to_csv(output_file, index=False)
                
                # Save high-quality subset
                hq_output_file = os.path.join(tss_dir, f'{category}_{reg_status}_genes_high_quality.csv')
                high_quality_df.to_csv(hq_output_file, index=False)
        
        results[category] = group_stats
    
    # Create visualization with CpG density consideration
    create_tss_analysis_plots(results, cell_type, tss_dir)
    
    # Save detailed statistics
    save_tss_analysis_results(results, cell_type, tss_dir)
    
    # Log quality metrics
    logger.info(f"\nQuality metrics for {cell_type}:")
    total_genes = len(df)
    low_cpg_genes = df['low_cpg_promoter'].sum()
    high_meth_genes = df['high_meth_promoter'].sum()
    logger.info(f"Total genes analyzed: {total_genes}")
    logger.info(f"Genes with low CpG density: {low_cpg_genes} ({low_cpg_genes/total_genes*100:.1f}%)")
    logger.info(f"Genes with suspiciously high methylation: {high_meth_genes} ({high_meth_genes/total_genes*100:.1f}%)")
    
    return results
