#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
from typing import Dict, List, Tuple
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def load_profile_data(base_dir: str, category: str) -> Dict[str, np.ndarray]:
    """Load saved profile data for a category."""
    data = {}
    
    # Load methylation data
    meth_types = ['basic', 'gc_norm', 'cpg_norm']
    for mtype in meth_types:
        try:
            path = os.path.join(base_dir, f'{category}_methylation_{mtype}.npy')
            data[f'methylation_{mtype}'] = np.load(path)
        except Exception as e:
            logger.warning(f"Could not load {path}: {str(e)}")
    
    # Load SMARCB1 data
    smarcb1_types = ['basic', 'replicate_norm']
    for stype in smarcb1_types:
        try:
            path = os.path.join(base_dir, f'{category}_smarcb1_{stype}.npy')
            data[f'smarcb1_{stype}'] = np.load(path)
        except Exception as e:
            logger.warning(f"Could not load {path}: {str(e)}")
    
    return data

def calculate_peak_metrics(profile: np.ndarray, window_size: int = 2000, bin_size: int = 50) -> Tuple[float, float, float]:
    """Calculate peak height, width, and area."""
    n_bins = window_size * 2 // bin_size
    center_idx = n_bins // 2
    
    # Find peak height
    peak_height = np.max(profile[center_idx-20:center_idx+20])
    
    # Calculate peak width (full width at half maximum)
    half_max = peak_height / 2
    above_half = profile > half_max
    if np.any(above_half):
        left_idx = np.where(above_half)[0][0]
        right_idx = np.where(above_half)[0][-1]
        peak_width = (right_idx - left_idx) * bin_size
    else:
        peak_width = 0
    
    # Calculate area under the curve
    peak_area = np.trapz(profile)
    
    return peak_height, peak_width, peak_area

def prepare_boxplot_data(profiles_data: Dict[str, Dict[str, np.ndarray]]) -> pd.DataFrame:
    """Prepare data for boxplots."""
    data_rows = []
    
    for category, data in profiles_data.items():
        for data_type, profile in data.items():
            height, width, area = calculate_peak_metrics(profile)
            
            data_rows.append({
                'Category': category,
                'Data_Type': data_type,
                'Peak_Height': height,
                'Peak_Width': width,
                'Peak_Area': area
            })
    
    return pd.DataFrame(data_rows)

def create_barplots(df: pd.DataFrame, output_dir: str):
    """Create and save bar plots with error bars."""
    metrics = ['Peak_Height', 'Peak_Width', 'Peak_Area']
    data_types = {
        'methylation': ['methylation_basic', 'methylation_gc_norm', 'methylation_cpg_norm'],
        'smarcb1': ['smarcb1_basic', 'smarcb1_replicate_norm']
    }
    
    for metric in metrics:
        for data_key, types in data_types.items():
            plt.figure(figsize=(10, 6))
            
            # Create bar plot with error bars
            plot_data = df[df['Data_Type'].isin(types)]
            sns.barplot(data=plot_data, x='Category', y=metric, hue='Data_Type', 
                       capsize=0.1, errwidth=1, ci=95)
            
            # Add statistical annotations
            categories = plot_data['Category'].unique()
            for type_name in types:
                y_max = plot_data[plot_data['Data_Type'] == type_name][metric].max()
                plt.text(0, y_max, f'{type_name}', ha='center', va='bottom')
            
            plt.title(f'{data_key.capitalize()} {metric} Distribution')
            plt.xticks(rotation=45)
            plt.tight_layout()
            
            # Save plot
            plt.savefig(os.path.join(output_dir, f'{data_key}_{metric.lower()}_barplot.png'))
            plt.close()

def run_statistical_tests(df: pd.DataFrame) -> pd.DataFrame:
    """Run statistical tests between categories for each metric and data type."""
    results = []
    metrics = ['Peak_Height', 'Peak_Width', 'Peak_Area']
    categories = df['Category'].unique()
    
    for metric in metrics:
        for data_type in df['Data_Type'].unique():
            data_subset = df[df['Data_Type'] == data_type]
            
            for i, cat1 in enumerate(categories):
                for cat2 in categories[i+1:]:
                    group1 = data_subset[data_subset['Category'] == cat1][metric]
                    group2 = data_subset[data_subset['Category'] == cat2][metric]
                    
                    # Mann-Whitney U test
                    stat, pval = stats.mannwhitneyu(group1, group2, alternative='two-sided')
                    
                    results.append({
                        'Metric': metric,
                        'Data_Type': data_type,
                        'Group1': cat1,
                        'Group2': cat2,
                        'Statistic': stat,
                        'P_value': pval,
                        'Significant': pval < 0.05
                    })
    
    return pd.DataFrame(results)

def main():
    # Set paths
    base_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
    working_dir = os.path.join(base_dir, "Cross_final")
    profiles_dir = os.path.join(working_dir, "results/methylation_analysis/integrated_profiles")
    output_dir = os.path.join(working_dir, "results/methylation_analysis/barplots")
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data for each category
    categories = ['up', 'down', 'no_deg']
    profiles_data = {
        category: load_profile_data(profiles_dir, category)
        for category in categories
    }
    
    # Prepare data for plotting
    plot_data = prepare_boxplot_data(profiles_data)
    
    # Create bar plots
    create_barplots(plot_data, output_dir)
    
    # Run statistical tests
    stats_results = run_statistical_tests(plot_data)
    
    # Save statistical results
    stats_results.to_csv(os.path.join(output_dir, 'statistical_tests.csv'), index=False)
    
    logger.info("Bar plot analysis completed successfully")
    logger.info(f"Results saved to {output_dir}")

if __name__ == "__main__":
    main()
