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
import config  # Add this import

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def load_profile_data(base_dir: str, category: str) -> Dict[str, np.ndarray]:
    """Load saved profile data for a category."""
    data = {}
    logger.info(f"Loading data for category: {category}")
    logger.info(f"Looking in directory: {base_dir}")
    
    # Load methylation data
    try:
        path = os.path.join(base_dir, f'{category}_methylation_basic.npy')
        logger.info(f"Attempting to load: {path}")
        if os.path.exists(path):
            data['methylation'] = np.load(path)
            logger.info(f"Successfully loaded {path}")
        else:
            logger.warning(f"File does not exist: {path}")
    except Exception as e:
        logger.error(f"Error loading {path}: {str(e)}")
    
    # Load SMARCB1 data
    try:
        path = os.path.join(base_dir, f'{category}_smarcb1_basic.npy')
        logger.info(f"Attempting to load: {path}")
        if os.path.exists(path):
            data['smarcb1'] = np.load(path)
            logger.info(f"Successfully loaded {path}")
        else:
            logger.warning(f"File does not exist: {path}")
    except Exception as e:
        logger.error(f"Error loading {path}: {str(e)}")
    
    logger.info(f"Loaded {len(data)} datasets for category {category}")
    return data

def calculate_peak_metrics(profile: np.ndarray, window_size: int = 2000, bin_size: int = 50) -> Tuple[float, float, float]:
    """
    Calculate peak metrics from binned profile data.
    The profile is a 1D array of 80 bins representing ±2000bp around peak centers.
    """
    try:
        # Calculate center region indices (±500bp)
        n_bins = len(profile)
        center_idx = n_bins // 2
        region_size = (500 // bin_size)  # 500bp converted to bins
        
        # Extract center region
        center_region = profile[center_idx-region_size:center_idx+region_size]
        
        # Find peak height
        peak_height = np.max(center_region)
        
        # Calculate peak width (full width at half maximum)
        half_max = peak_height / 2
        above_half = center_region > half_max
        if np.any(above_half):
            left_idx = np.where(above_half)[0][0]
            right_idx = np.where(above_half)[0][-1]
            peak_width = (right_idx - left_idx) * bin_size
        else:
            peak_width = 0
        
        # Calculate area under the curve
        peak_area = np.trapz(center_region) * bin_size
        
        return peak_height, peak_width, peak_area
        
    except Exception as e:
        logger.error(f"Error calculating peak metrics: {str(e)}")
        return 0.0, 0.0, 0.0

def prepare_boxplot_data(profiles_data: Dict[str, Dict[str, np.ndarray]]) -> pd.DataFrame:
    """Prepare data for boxplots."""
    data_rows = []
    
    for category, data in profiles_data.items():
        for data_type, profile in data.items():
            height, width, area = calculate_peak_metrics(profile)
            
            data_rows.append({
                'Category': category,
                'DataType': data_type,
                'Peak_Height': height,
                'Peak_Width': width,
                'Peak_Area': area
            })
    
    df = pd.DataFrame(data_rows)
    logger.info(f"DataFrame columns: {df.columns.tolist()}")
    logger.info(f"First few rows of DataFrame:\n{df.head()}")
    return df

def create_barplots(df: pd.DataFrame, output_dir: str):
    """Create and save separate bar plots for methylation and SMARCB1."""
    logger.info(f"Columns in DataFrame passed to create_barplots: {df.columns.tolist()}")
    
    metrics = ['Peak_Area', 'Peak_Height']  # Only using Area and Height metrics
    signals = ['methylation', 'smarcb1']
    
    # Create separate plots for each signal type
    for signal in signals:
        signal_data = df[df['DataType'] == signal]
        
        for metric in metrics:
            plt.figure(figsize=(8, 6))
            
            # Create barplot for single signal type
            sns.barplot(data=signal_data, x='Category', y=metric,
                       capsize=0.1, err_kws={'linewidth': 1})
            
            plt.title(f'{signal.upper()} {metric.replace("_", " ")} Distribution')
            plt.ylabel(f'{metric.replace("_", " ")}')
            plt.xlabel('Category')
            plt.xticks(rotation=45)
            
            # Statistical annotations
            categories = signal_data['Category'].unique()
            y_max = signal_data[metric].max()
            y_offset = y_max * 0.1
            
            # Perform statistical tests between categories
            for i, cat1 in enumerate(categories):
                for cat2 in categories[i+1:]:
                    group1 = signal_data[signal_data['Category'] == cat1][metric]
                    group2 = signal_data[signal_data['Category'] == cat2][metric]
                    
                    if len(group1) > 0 and len(group2) > 0:
                        stat, pval = stats.mannwhitneyu(group1, group2)
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
            output_path = os.path.join(output_dir, f'{signal}_{metric.lower()}_barplot.png')
            plt.savefig(output_path, bbox_inches='tight', dpi=300)
            logger.info(f"Saved plot to: {output_path}")
            plt.close()

def run_statistical_tests(df: pd.DataFrame) -> pd.DataFrame:
    """Run statistical tests between categories for Peak_Area and Peak_Height."""
    results = []
    metrics = ['Peak_Area', 'Peak_Height']  # Updated to use only Area and Height
    categories = df['Category'].unique()
    
    for metric in metrics:
        for data_type in df['DataType'].unique():
            data_subset = df[df['DataType'] == data_type]
            
            for i, cat1 in enumerate(categories):
                for cat2 in categories[i+1:]:
                    group1 = data_subset[data_subset['Category'] == cat1][metric]
                    group2 = data_subset[data_subset['Category'] == cat2][metric]
                    
                    # Mann-Whitney U test
                    stat, pval = stats.mannwhitneyu(group1, group2, alternative='two-sided')
                    
                    results.append({
                        'Metric': metric,
                        'DataType': data_type,
                        'Group1': cat1,
                        'Group2': cat2,
                        'Statistic': stat,
                        'P_value': pval,
                        'Significant': pval < 0.05
                    })
    
    return pd.DataFrame(results)

def main():
    # Update paths to use config
    base_dir = config.BASE_DIR
    working_dir = config.BASE_DIR  # Since BASE_DIR already points to Cross_final
    profiles_dir = os.path.join(config.RESULTS_DIR, "3integrated_analysis")
    output_dir = os.path.join(config.RESULTS_DIR, "5create_boxplots")
    os.makedirs(output_dir, exist_ok=True)
    
    logger.info(f"Working with profiles directory: {profiles_dir}")
    logger.info(f"Output directory: {output_dir}")
    
    # Load data for each category
    categories = ['up', 'down', 'no_deg']
    plot_data = []
    
    for category in categories:
        # Load profiles
        profiles = load_profile_data(profiles_dir, category)
        
        # Calculate metrics for each signal type
        for signal_type, profile in profiles.items():
            if profile is not None:
                height, width, area = calculate_peak_metrics(profile)
                plot_data.append({
                    'Category': category,
                    'DataType': signal_type,
                    'Peak_Height': height,
                    'Peak_Width': width,
                    'Peak_Area': area
                })
    
    # Create DataFrame
    plot_df = pd.DataFrame(plot_data)
    
    if plot_df.empty:
        logger.error("No data was loaded successfully. Check file paths and data files.")
        return
    
    # Create plots
    create_barplots(plot_df, output_dir)
    
    # Save summary statistics
    summary_stats = plot_df.groupby(['Category', 'DataType']).agg({
        'Peak_Height': ['mean', 'std'],
        'Peak_Width': ['mean', 'std'],
        'Peak_Area': ['mean', 'std']
    }).round(3)
    
    summary_stats.to_csv(os.path.join(output_dir, 'peak_metrics_summary.csv'))
    
    logger.info("Bar plot analysis completed successfully")
    logger.info(f"Results saved to {output_dir}")

if __name__ == "__main__":
    main()
