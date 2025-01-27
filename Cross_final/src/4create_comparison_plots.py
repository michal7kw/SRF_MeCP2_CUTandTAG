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
    try:
        path = os.path.join(base_dir, f'{category}_methylation_basic.npy')
        data['methylation'] = np.load(path)
    except Exception as e:
        logger.warning(f"Could not load {path}: {str(e)}")
    
    # Load SMARCB1 data
    try:
        path = os.path.join(base_dir, f'{category}_smarcb1_basic.npy')
        data['smarcb1'] = np.load(path)
    except Exception as e:
        logger.warning(f"Could not load {path}: {str(e)}")
    
    return data

def calculate_signal_metrics(profile: np.ndarray, window_size: int = 2000, bin_size: int = 50) -> Dict[str, np.ndarray]:
    """
    Calculate signal metrics from binned profile data.
    The profile is a 1D array of 80 bins representing ±2000bp around peak centers.
    """
    logger.info(f"Profile shape: {profile.shape}")
    
    try:
        # Calculate center region indices (±500bp)
        n_bins = len(profile)
        center_idx = n_bins // 2
        region_size = (500 // bin_size)  # 500bp converted to bins
        
        # Extract center region
        center_region = profile[center_idx-region_size:center_idx+region_size]
        
        # Calculate metrics for the center region
        metrics = {
            'mean_signal': np.array([np.mean(center_region)]),
            'max_signal': np.array([np.max(center_region)]),
            'total_signal': np.array([np.sum(center_region)])
        }
        
        # Log the shape of the calculated metrics
        for metric_name, metric_values in metrics.items():
            logger.info(f"{metric_name} shape: {metric_values.shape}")
        
        return metrics
        
    except Exception as e:
        logger.error(f"Error calculating metrics: {str(e)}")
        return {}

def create_comparison_plots(data: Dict[str, Dict[str, Dict[str, np.ndarray]]], output_dir: str):
    """Create separate comparison plots for methylation and SMARCB1 binding."""
    categories = ['up', 'down', 'no_deg']
    signals = ['methylation', 'smarcb1']
    metrics = ['mean_signal', 'max_signal', 'total_signal']
    
    # Create separate plots for each signal type
    for signal in signals:
        # Prepare data for this signal type
        plot_data = []
        for category in categories:
            if signal in data[category]:
                metrics_dict = data[category][signal]
                for metric, values in metrics_dict.items():
                    if values is not None and len(values) > 0:
                        plot_data.append({
                            'Category': category,
                            'Metric': metric,
                            'Value': float(values[0])
                        })
        
        # Create DataFrame for this signal
        df = pd.DataFrame(plot_data)
        if df.empty:
            logger.warning(f"No valid data for {signal}")
            continue
        
        # Create plots for each metric for this signal
        for metric in metrics:
            metric_data = df[df['Metric'] == metric]
            if metric_data.empty:
                logger.warning(f"No data for metric: {metric} in {signal}")
                continue
            
            plt.figure(figsize=(8, 6))
            
            # Create barplot
            sns.barplot(data=metric_data,
                       x='Category', y='Value',
                       capsize=0.1, err_kws={'linewidth': 1})
            
            plt.title(f'{signal.upper()} {metric.replace("_", " ")} Around Peak Centers (±500bp)')
            plt.ylabel('Signal Intensity')
            plt.xlabel('Category')
            plt.xticks(rotation=45)
            
            # Add statistical annotations
            categories = metric_data['Category'].unique()
            y_max = metric_data['Value'].max()
            y_offset = y_max * 0.1
            
            for i, cat1 in enumerate(categories):
                for cat2 in categories[i+1:]:
                    group1 = metric_data[metric_data['Category'] == cat1]['Value']
                    group2 = metric_data[metric_data['Category'] == cat2]['Value']
                    
                    if len(group1) > 0 and len(group2) > 0:
                        stat, pval = stats.mannwhitneyu(group1, group2, alternative='two-sided')
                        if pval < 0.05:  # Only add significant results
                            # Calculate position for annotation
                            x_pos = (i + 1) / 2
                            y_pos = y_max + y_offset
                            
                            # Add asterisks based on significance level
                            sig_symbol = '***' if pval < 0.001 else ('**' if pval < 0.01 else '*')
                            
                            plt.text(x_pos, y_pos,
                                   f'{cat1} vs {cat2}\n{sig_symbol} p={pval:.2e}',
                                   ha='center', va='bottom',
                                   bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
                            
                            y_offset += y_max * 0.15
            
            # Adjust plot limits to accommodate annotations
            plt.ylim(0, y_max * (1.5 if y_offset > y_max * 0.5 else 1.2))
            
            plt.tight_layout()
            output_path = os.path.join(output_dir, f'{signal}_{metric.lower()}.png')
            plt.savefig(output_path, bbox_inches='tight', dpi=300)
            logger.info(f"Saved plot to: {output_path}")
            plt.close()

def run_statistical_analysis(data: Dict[str, Dict[str, Dict[str, np.ndarray]]]) -> pd.DataFrame:
    """Run statistical tests between categories for each signal type."""
    results = []
    categories = ['up', 'down', 'no_deg']
    signals = ['methylation', 'smarcb1']
    metrics = ['mean_signal', 'max_signal', 'total_signal']
    
    for signal in signals:
        for metric in metrics:
            for i, cat1 in enumerate(categories):
                for cat2 in categories[i+1:]:
                    if signal in data[cat1] and signal in data[cat2]:
                        vals1 = data[cat1][signal][metric]
                        vals2 = data[cat2][signal][metric]
                        
                        # Mann-Whitney U test on full distributions
                        stat, pval = stats.mannwhitneyu(vals1, vals2, alternative='two-sided')
                        
                        # Calculate fold change using means
                        fold_change = np.mean(vals1) / np.mean(vals2) if np.mean(vals2) != 0 else float('inf')
                        
                        results.append({
                            'Signal': signal,
                            'Metric': metric,
                            'Group1': cat1,
                            'Group2': cat2,
                            'Group1_Mean': np.mean(vals1),
                            'Group2_Mean': np.mean(vals2),
                            'Group1_SD': np.std(vals1),
                            'Group2_SD': np.std(vals2),
                            'Fold_Change': fold_change,
                            'P_value': pval,
                            'Significant': pval < 0.05
                        })
    
    return pd.DataFrame(results)

def main():
    # Set paths
    base_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
    working_dir = os.path.join(base_dir, "Cross_final")
    profiles_dir = os.path.join(working_dir, "results/3integrated_analysis")
    output_dir = os.path.join(working_dir, "results/4create_comparison_plots")
    os.makedirs(output_dir, exist_ok=True)
    
    # Add logging of directories
    logger.info(f"Input directory: {profiles_dir}")
    logger.info(f"Output directory: {output_dir}")
    
    # Load and process data
    categories = ['up', 'down', 'no_deg']
    data = {}
    
    for category in categories:
        # Load basic profiles (we'll use the basic normalization)
        profiles = load_profile_data(profiles_dir, category)
        
        # Calculate metrics for each signal type
        data[category] = {}
        for signal_type, profile in profiles.items():
            if profile is not None:  # Check if profile was loaded successfully
                data[category][signal_type] = calculate_signal_metrics(profile)
    
    # Create plots
    create_comparison_plots(data, output_dir)
    
    # Run statistical analysis
    # stats_results = run_statistical_analysis(data)
    # stats_results.to_csv(os.path.join(output_dir, 'statistical_analysis.csv'), index=False)
    
    # Create summary table
    summary = pd.DataFrame([
        {
            'Category': cat,
            'Signal': signal,
            **metrics
        }
        for cat, signals in data.items()
        for signal, metrics in signals.items()
    ])
    summary.to_csv(os.path.join(output_dir, 'signal_metrics_summary.csv'), index=False)
    
    # Add file existence checks after saving
    expected_files = [
        'mean_signal_comparison.png',
        'max_signal_comparison.png',
        'total_signal_comparison.png',
        'signal_metrics_summary.csv'
    ]
    
    for filename in expected_files:
        filepath = os.path.join(output_dir, filename)
        if os.path.exists(filepath):
            logger.info(f"Successfully created: {filename}")
        else:
            logger.warning(f"Failed to create: {filename}")
    
    logger.info("Comparison analysis completed successfully")
    logger.info(f"Results saved to {output_dir}")

if __name__ == "__main__":
    main()
