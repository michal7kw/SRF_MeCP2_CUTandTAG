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

def calculate_signal_metrics(profile: np.ndarray, window_size: int = 2000, bin_size: int = 50) -> Dict[str, float]:
    """Calculate signal metrics focusing on the peak region."""
    n_bins = window_size * 2 // bin_size
    center_idx = n_bins // 2
    
    # Focus on central region (±500bp around peak center)
    region_size = 500 // bin_size
    central_region = profile[center_idx-region_size:center_idx+region_size]
    
    metrics = {
        'mean_signal': np.mean(central_region),
        'max_signal': np.max(central_region),
        'total_signal': np.sum(central_region)
    }
    
    return metrics

def create_comparison_plots(data: Dict[str, Dict[str, Dict[str, float]]], output_dir: str):
    """Create comparison plots for methylation and SMARCB1 binding."""
    categories = ['up', 'down', 'no_deg']
    signals = ['methylation', 'smarcb1']
    metrics = ['mean_signal', 'max_signal', 'total_signal']
    
    # Prepare data for plotting
    plot_data = []
    for category in categories:
        for signal in signals:
            if signal in data[category]:
                metrics_dict = data[category][signal]
                for metric, value in metrics_dict.items():
                    plot_data.append({
                        'Category': category,
                        'Signal': signal,
                        'Metric': metric,
                        'Value': value
                    })
    
    df = pd.DataFrame(plot_data)
    
    # Create plots for each metric
    for metric in metrics:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Boxplot
        sns.boxplot(data=df[df['Metric'] == metric], 
                   x='Category', y='Value', hue='Signal',
                   ax=ax1)
        ax1.set_title(f'{metric.replace("_", " ").title()} Distribution')
        ax1.set_ylabel('Signal Intensity')
        ax1.tick_params(axis='x', rotation=45)
        
        # Barplot with error bars
        sns.barplot(data=df[df['Metric'] == metric],
                   x='Category', y='Value', hue='Signal',
                   ax=ax2, errorbar='sd')
        ax2.set_title(f'Mean {metric.replace("_", " ").title()} with SD')
        ax2.set_ylabel('Signal Intensity')
        ax2.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{metric}_comparison.png'))
        plt.close()
        
        # Create separate plots for each signal type
        for signal in signals:
            signal_data = df[(df['Metric'] == metric) & (df['Signal'] == signal)]
            
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # Boxplot
            sns.boxplot(data=signal_data, x='Category', y='Value',
                       ax=ax1, color=sns.color_palette()[0])
            ax1.set_title(f'{signal.title()} {metric.replace("_", " ").title()} Distribution')
            ax1.set_ylabel('Signal Intensity')
            ax1.tick_params(axis='x', rotation=45)
            
            # Add statistical annotations
            cat_pairs = [(categories[i], categories[j]) 
                        for i in range(len(categories)) 
                        for j in range(i+1, len(categories))]
            
            max_y = signal_data['Value'].max()
            y_positions = np.linspace(max_y*1.05, max_y*1.2, len(cat_pairs))
            
            for (cat1, cat2), y_pos in zip(cat_pairs, y_positions):
                stat, pval = stats.mannwhitneyu(
                    signal_data[signal_data['Category'] == cat1]['Value'],
                    signal_data[signal_data['Category'] == cat2]['Value']
                )
                x1 = categories.index(cat1)
                x2 = categories.index(cat2)
                
                # Draw significance bars
                ax1.plot([x1, x2], [y_pos, y_pos], '-k', linewidth=1)
                ax1.plot([x1, x1], [y_pos*0.98, y_pos], '-k', linewidth=1)
                ax1.plot([x2, x2], [y_pos*0.98, y_pos], '-k', linewidth=1)
                
                # Add p-value text
                ax1.text((x1+x2)/2, y_pos*1.01, f'p={pval:.2e}', 
                        ha='center', va='bottom')
            
            # Barplot with error bars
            sns.barplot(data=signal_data, x='Category', y='Value',
                       ax=ax2, color=sns.color_palette()[0], errorbar='sd')
            ax2.set_title(f'Mean {signal.title()} {metric.replace("_", " ").title()} with SD')
            ax2.set_ylabel('Signal Intensity')
            ax2.tick_params(axis='x', rotation=45)
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'{signal}_{metric}_comparison.png'))
            plt.close()

def run_statistical_analysis(data: Dict[str, Dict[str, Dict[str, float]]]) -> pd.DataFrame:
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
                        val1 = data[cat1][signal][metric]
                        val2 = data[cat2][signal][metric]
                        
                        # Mann-Whitney U test
                        stat, pval = stats.mannwhitneyu([val1], [val2], alternative='two-sided')
                        
                        # Calculate fold change
                        fold_change = val1 / val2 if val2 != 0 else float('inf')
                        
                        results.append({
                            'Signal': signal,
                            'Metric': metric,
                            'Group1': cat1,
                            'Group2': cat2,
                            'Group1_Value': val1,
                            'Group2_Value': val2,
                            'Fold_Change': fold_change,
                            'P_value': pval,
                            'Significant': pval < 0.05
                        })
    
    return pd.DataFrame(results)

def main():
    # Set paths
    base_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
    working_dir = os.path.join(base_dir, "Cross_final")
    profiles_dir = os.path.join(working_dir, "results/methylation_analysis/integrated_profiles")
    output_dir = os.path.join(working_dir, "results/methylation_analysis/comparison_plots")
    os.makedirs(output_dir, exist_ok=True)
    
    # Load and process data
    categories = ['up', 'down', 'no_deg']
    data = {}
    
    for category in categories:
        # Load profiles
        profiles = load_profile_data(profiles_dir, category)
        
        # Calculate metrics for each signal type
        data[category] = {}
        for signal_type, profile in profiles.items():
            data[category][signal_type] = calculate_signal_metrics(profile)
    
    # Create plots
    create_comparison_plots(data, output_dir)
    
    # Run statistical analysis
    stats_results = run_statistical_analysis(data)
    stats_results.to_csv(os.path.join(output_dir, 'statistical_analysis.csv'), index=False)
    
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
    
    logger.info("Comparison analysis completed successfully")
    logger.info(f"Results saved to {output_dir}")

if __name__ == "__main__":
    main()
