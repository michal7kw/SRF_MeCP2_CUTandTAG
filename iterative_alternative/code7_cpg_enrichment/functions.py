import pandas as pd
import yaml
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def plot_enrichment_distribution(df, factor=3, title=None):
    """
    Create histogram of enrichment values for regions bound by both Exo and Endo.
    
    Args:
        df: DataFrame containing enrichment values
        factor: Factor to multiply IQR by for determining plot range (default=3)
        title: Optional custom title for the plot
    """
    plt.figure(figsize=(10, 6))

    # Calculate reasonable range based on data distribution
    q1, q3 = df['enrichment'].quantile([0.25, 0.75])
    iqr = q3 - q1
    range_min = max(0, q1 - factor * iqr)  # Don't go below 0 for enrichment
    range_max = q3 + factor * iqr

    # Plot histogram with better binning and transparency
    n, bins, patches = plt.hist(df['enrichment'], bins=50, edgecolor='black', alpha=0.7,
                              color='#2196F3', density=False, range=(range_min, range_max))

    # Add grid for better readability
    plt.grid(True, alpha=0.3)

    # Improve axis labels and title
    plt.xlabel('Enrichment Score', fontsize=12)
    plt.ylabel('Number of Regions', fontsize=12)
    if title is None:
        title = 'Distribution of Enrichment Values\nin Regions Bound by Both Exo and Endo'
    else:
        title = f'Distribution of Enrichment Values\nin Regions Bound by Both Exo and Endo\n{title}'
    plt.title(title, fontsize=14, pad=15)

    # Add vertical line at enrichment = 1
    plt.axvline(x=1, color='red', linestyle='--', linewidth=2,
                label='Enrichment = 1 (No difference)')

    # Set x-axis limits based on the actual data range
    plt.xlim(range_min, range_max)

    # Add legend with better positioning
    plt.legend(loc='upper right')

    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    plt.show()
    
    return range_min, range_max

def plot_enrichment_distribution_neu_vs_nsc(df, factor=3, title=None):
    """
    Create histogram of enrichment values for regions bound by both Neu and NSC.
    
    Args:
        df: DataFrame containing enrichment values
        factor: Factor to multiply IQR by for determining plot range (default=3)
        title: Optional custom title for the plot
    """
    plt.figure(figsize=(10, 6))

    # Calculate reasonable range based on data distribution
    q1, q3 = df['enrichment'].quantile([0.25, 0.75])
    iqr = q3 - q1
    range_min = max(0, q1 - factor * iqr)  # Don't go below 0 for enrichment
    range_max = q3 + factor * iqr

    # Plot histogram with better binning and transparency
    n, bins, patches = plt.hist(df['enrichment'], bins=50, edgecolor='black', alpha=0.7,
                              color='#2196F3', density=False, range=(range_min, range_max))

    # Add grid for better readability
    plt.grid(True, alpha=0.3)

    # Improve axis labels and title
    plt.xlabel('Enrichment Score', fontsize=12)
    plt.ylabel('Number of Regions', fontsize=12)
    if title is None:
        title = 'Distribution of Enrichment Values\nin Regions Bound by Both Neu and NSC'
    else:
        title = f'Distribution of Enrichment Values\nin Regions Bound by Both Neu and NSC\n{title}'
    plt.title(title, fontsize=14, pad=15)

    # Add vertical line at enrichment = 1
    plt.axvline(x=1, color='red', linestyle='--', linewidth=2,
                label='Enrichment = 1 (No difference)')

    # Set x-axis limits based on the actual data range
    plt.xlim(range_min, range_max)

    # Add legend with better positioning
    plt.legend(loc='upper right')

    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    plt.show()
    
    return range_min, range_max

def plot_enrichment_boxplot(df, factor=None, title=None):
    """
    Create a box plot to visualize distribution and outliers of enrichment values.
    
    Parameters:
        df (pd.DataFrame): DataFrame containing enrichment values
        factor (float): Factor controlling whisker length in boxplot
        title (str): Custom title for the plot. If None, uses default title.
    """
    plt.figure(figsize=(10, 6))
    # The whis parameter controls the whisker length in the boxplot
    # whis=factor means whiskers extend to points within factor * IQR of Q1/Q3
    # Points beyond the whiskers are considered outliers and plotted individually
    plt.boxplot(df['enrichment'], whis=factor)
    plt.ylabel('Enrichment Score', fontsize=12)
    plt.title(f'Box Plot of Enrichment Values\nShowing Distribution and Outliers\n{title}' if title else 'Box Plot of Enrichment Values\nShowing Distribution and Outliers',
              fontsize=14, pad=15)
    plt.grid(True, alpha=0.3)

    # Add horizontal line at enrichment = 1
    plt.axhline(y=1, color='red', linestyle='--', linewidth=2,
                label='Enrichment = 1 (No difference)')

    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_by_peaks(outliers_df, column_name='enrichment', peaks_column='endo_replicates_with_peaks', title=None):
    """
    Classify regions based on number of endo peaks and plot overlapping histograms.
    
    Parameters:
        outliers_df (pd.DataFrame): DataFrame containing outlier regions
        column_name (str): Name of the column to plot distribution for
        peaks_column (str): Name of the column containing number of peaks
        title (str): Custom title for the plot. If None, uses default title.
    """

    # Classify regions based on number of endo peaks
    peaks_zero_peaks = outliers_df[outliers_df[peaks_column] == 0]
    peaks_one_peaks = outliers_df[outliers_df[peaks_column] == 1] 
    peaks_two_peaks = outliers_df[outliers_df[peaks_column] == 2]
    peaks_three_peaks = outliers_df[outliers_df[peaks_column] == 3]

    plt.figure(figsize=(10, 6))

    # Plot overlapping histograms with transparency
    # Handle infinite values by filtering them out
    for data, label, color in [
        (peaks_zero_peaks[column_name], '0 Peaks', 'blue'),
        (peaks_one_peaks[column_name], '1 Peaks', 'red'),
        (peaks_two_peaks[column_name], '2 Peaks', 'green'),
        (peaks_three_peaks[column_name], '3 Peaks', 'purple')
    ]:
        # Filter out infinite values
        finite_data = data[np.isfinite(data)]
        if len(finite_data) > 0:
            plt.hist(finite_data, bins=30, alpha=0.5, label=label, color=color)

    plt.xlabel(column_name.replace('_',' ').title())
    plt.ylabel('Frequency')
    plt.title(title if title else f'Distribution of {column_name.replace("_"," ").title()} by Number of Peaks')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    return peaks_zero_peaks, peaks_one_peaks, peaks_two_peaks, peaks_three_peaks

def plot_region_length_distribution(df, title=None):
    plt.figure(figsize=(10, 6))

    # Calculate reasonable range based on data distribution
    q1, q3 = df['region_length'].quantile([0.25, 0.75])
    factor = 3
    iqr = q3 - q1
    range_min = max(0, q1 - factor * iqr)  # Don't go below 0 for region_length
    range_max = q3 + factor * iqr

    n, bins, patches = plt.hist(df['region_length'], bins=50, edgecolor='black', alpha=0.7,
                              color='#2196F3', density=False, range=(range_min, range_max),
                              label='Region Length Distribution')

    plt.grid(True, alpha=0.3)

    plt.xlabel('Region Length (bp)', fontsize=12)
    plt.ylabel('Number of Regions', fontsize=12)
    if title is None:
        title = 'Distribution of Region Lengths'
    else:
        title = f'Distribution of Region Lengths\n{title}'
    plt.title(title, fontsize=14, pad=15)

    plt.xlim(range_min, range_max)

    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.show()

def plot_region_length_comparison(df_by_peaks, df_by_signal, region_length_col='region_length'):
    plt.figure(figsize=(10, 6))

    # Calculate reasonable ranges based on both datasets
    q1_peaks, q3_peaks = df_by_peaks[region_length_col].quantile([0.25, 0.75])
    q1_signal, q3_signal = df_by_signal[region_length_col].quantile([0.25, 0.75])

    factor = 3
    iqr_peaks = q3_peaks - q1_peaks
    iqr_signal = q3_signal - q1_signal

    range_min = max(0, min(q1_peaks - factor * iqr_peaks, q1_signal - factor * iqr_signal))
    range_max = max(q3_peaks + factor * iqr_peaks, q3_signal + factor * iqr_signal)

    # Plot both distributions
    plt.hist(df_by_peaks[region_length_col], bins=50, edgecolor='black', alpha=0.5,
             color='#2196F3', density=False, range=(range_min, range_max),
             label=f'Selected by Peaks (n={len(df_by_peaks)})')

    plt.hist(df_by_signal[region_length_col], bins=50, edgecolor='black', alpha=0.5,
             color='#FF9800', density=False, range=(range_min, range_max),
             label=f'Selected by Signal (n={len(df_by_signal)})')

    plt.grid(True, alpha=0.3)
    plt.xlabel('Region Length (bp)', fontsize=12)
    plt.ylabel('Number of Regions', fontsize=12)
    plt.title('Distribution of Region Lengths\nComparison between Peak and Signal Selection', fontsize=14, pad=15)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()

def plot_exo_only_distributions(exo_only_df, title=None):
    # Distribution plots
    plt.figure(figsize=(15, 5))

    # Plot overlapped signal distributions
    plt.subplot(1, 3, 1)
    sns.histplot(data=exo_only_df, x='exo_signal', bins=50, alpha=0.5, label='Exo Signal', color='blue')
    sns.histplot(data=exo_only_df, x='endo_signal', bins=50, alpha=0.5, label='Endo Signal', color='red')
    if title is None:
        _title = 'Distribution of Signals'
    else:
        _title = f'Distribution of Signals\n{title}'
    plt.title(_title)
    plt.xlabel('Signal Value')
    plt.ylabel('Count')
    plt.legend()

    # Plot overlapped length distributions
    plt.subplot(1, 3, 2)
    sns.histplot(data=exo_only_df, x='region_length', bins=50, alpha=0.5, label='Region Length', color='blue')
    sns.histplot(data=exo_only_df, x='cpg_length', bins=50, alpha=0.5, label='CpG Length', color='red')
    if title is None:
        _title = 'Distribution of Lengths'
    else:
        _title = f'Distribution of Lengths\n{title}'
    plt.title(_title)
    plt.xlabel('Length')
    plt.ylabel('Count')
    plt.legend()

    # Count peaks per region
    plt.subplot(1, 3, 3)
    exo_peak_counts = exo_only_df['exo_replicates_with_peaks'].value_counts()
    endo_peak_counts = exo_only_df['endo_replicates_with_peaks'].value_counts()

    x = np.arange(max(max(exo_peak_counts.index), max(endo_peak_counts.index)) + 1)
    width = 0.35

    # Create bars and store them to add labels
    exo_bars = plt.bar(x - width/2, [exo_peak_counts.get(i, 0) for i in x], width, label='Exo', color='blue', alpha=0.5)
    endo_bars = plt.bar(x + width/2, [endo_peak_counts.get(i, 0) for i in x], width, label='Endo', color='red', alpha=0.5)

    # Add count labels on top of each bar
    for bars in [exo_bars, endo_bars]:
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}',
                    ha='center', va='bottom')
    if title is None:
        _title = 'Number of Replicates with Peaks'
    else:
        _title = f'Number of Replicates with Peaks\n{title}'
    plt.title(_title)
    plt.xlabel('Number of Replicates')
    plt.ylabel('Count')
    plt.legend()

    plt.xticks(x, x.astype(int))
    plt.tight_layout()
    plt.show()

def print_outlier_groups(df, selection_type="signal"):
    n_rows = len(df)
    
    print(f"\nFirst 5 outliers selected by {selection_type}:")
    for _, row in df[:5].iterrows():
        print(f"{row['chr']}:{row['start']:,}-{row['end']:,}".ljust(50) + f"enrichment: {row['enrichment']}")

    if n_rows > 10:
        mid_start = n_rows//2 - 2
        print(f"\nMiddle 5 outliers selected by {selection_type}:")
        for _, row in df[mid_start:mid_start+5].iterrows():
            print(f"{row['chr']}:{row['start']:,}-{row['end']:,}".ljust(50) + f"enrichment: {row['enrichment']}")

    print(f"\nLast 5 outliers selected by {selection_type}:")
    for _, row in df[-5:].iterrows():
        print(f"{row['chr']}:{row['start']:,}-{row['end']:,}".ljust(50) + f"enrichment: {row['enrichment']}")

def plot_endo_only_distributions(endo_only_df, title=None):
    # Distribution plots
    plt.figure(figsize=(15, 5))

    # Plot overlapped signal distributions
    plt.subplot(1, 3, 1)
    sns.histplot(data=endo_only_df, x='endo_signal', bins=50, alpha=0.5, label='endo Signal', color='blue')
    sns.histplot(data=endo_only_df, x='exo_signal', bins=50, alpha=0.5, label='Exo Signal', color='red')
    if title is None:
        _title = 'Distribution of Signals'
    else:
        _title = f'Distribution of Signals\n{title}'
    plt.title(_title)
    plt.xlabel('Signal Value')
    plt.ylabel('Count')
    plt.legend()

    # Plot overlapped length distributions
    plt.subplot(1, 3, 2)
    sns.histplot(data=endo_only_df, x='region_length', bins=50, alpha=0.5, label='Region Length', color='blue')
    sns.histplot(data=endo_only_df, x='cpg_length', bins=50, alpha=0.5, label='CpG Length', color='red')
    if title is None:
        _title = 'Distribution of Lengths'
    else:
        _title = f'Distribution of Lengths\n{title}'
    plt.title(_title)
    plt.xlabel('Length')
    plt.ylabel('Count')
    plt.legend()

    # Count peaks per region
    plt.subplot(1, 3, 3)
    exo_peak_counts = endo_only_df['exo_replicates_with_peaks'].value_counts()
    endo_peak_counts = endo_only_df['endo_replicates_with_peaks'].value_counts()

    x = np.arange(max(max(exo_peak_counts.index), max(endo_peak_counts.index)) + 1)
    width = 0.35

    # Create bars and store them to add labels
    exo_bars = plt.bar(x - width/2, [exo_peak_counts.get(i, 0) for i in x], width, label='Exo', color='blue', alpha=0.5)
    endo_bars = plt.bar(x + width/2, [endo_peak_counts.get(i, 0) for i in x], width, label='Endo', color='red', alpha=0.5)

    # Add count labels on top of each bar
    for bars in [endo_bars, endo_bars]:
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}',
                    ha='center', va='bottom')
    if title is None:
        _title = 'Number of Replicates with Peaks'
    else:
        _title = f'Number of Replicates with Peaks\n{title}'
    plt.title(_title)
    plt.xlabel('Number of Replicates')
    plt.ylabel('Count')
    plt.legend()

    plt.xticks(x, x.astype(int))
    plt.tight_layout()
    plt.show()

def print_across_distribution(df, selection_type="signal"):
    n_rows = len(df)
    
    print(f"\nFirst 5 elements selected by {selection_type}:")
    for _, row in df[:5].iterrows():
        print(f"{row['chr']}:{row['start']:,}-{row['end']:,}".ljust(50) + f"enrichment: {row['enrichment']}")

    if n_rows > 10:
        mid_start = n_rows//2 - 2
        print(f"\nMiddle 5 elements selected by {selection_type}:")
        for _, row in df[mid_start:mid_start+5].iterrows():
            print(f"{row['chr']}:{row['start']:,}-{row['end']:,}".ljust(50) + f"enrichment: {row['enrichment']}")

    print(f"\nLast 5 elements selected by {selection_type}:")
    for _, row in df[-5:].iterrows():
        print(f"{row['chr']}:{row['start']:,}-{row['end']:,}".ljust(50) + f"enrichment: {row['enrichment']}")


def plot_exo_only_distributions_no_length_distribution(exo_only_df, title=None):
    # Distribution plots
    plt.figure(figsize=(10, 5))

    # Plot overlapped signal distributions
    plt.subplot(1, 2, 1)
    sns.histplot(data=exo_only_df, x='exo_signal', bins=50, alpha=0.5, label='Exo Signal', color='blue')
    sns.histplot(data=exo_only_df, x='endo_signal', bins=50, alpha=0.5, label='Endo Signal', color='red')
    if title is None:
        _title = 'Distribution of Signals'
    else:
        _title = f'Distribution of Signals\n{title}'
    plt.title(_title)
    plt.xlabel('Signal Value')
    plt.ylabel('Count')
    plt.legend()

    # Count peaks per region
    plt.subplot(1, 2, 2)
    exo_peak_counts = exo_only_df['exo_replicates_with_peaks'].value_counts()
    endo_peak_counts = exo_only_df['endo_replicates_with_peaks'].value_counts()

    x = np.arange(max(max(exo_peak_counts.index), max(endo_peak_counts.index)) + 1)
    width = 0.35

    # Create bars and store them to add labels
    exo_bars = plt.bar(x - width/2, [exo_peak_counts.get(i, 0) for i in x], width, label='Exo', color='blue', alpha=0.5)
    endo_bars = plt.bar(x + width/2, [endo_peak_counts.get(i, 0) for i in x], width, label='Endo', color='red', alpha=0.5)

    # Add count labels on top of each bar
    for bars in [exo_bars, endo_bars]:
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}',
                    ha='center', va='bottom')
    if title is None:
        _title = 'Number of Replicates with Peaks'
    else:
        _title = f'Number of Replicates with Peaks\n{title}'
    plt.title(_title)
    plt.xlabel('Number of Replicates')
    plt.ylabel('Count')
    plt.legend()

    plt.xticks(x, x.astype(int))
    plt.tight_layout()
    plt.show()


def plot_neu_only_distributions(neu_only_df, title=None):
    # Distribution plots
    plt.figure(figsize=(15, 5))

    # Plot overlapped signal distributions
    plt.subplot(1, 3, 1)
    sns.histplot(data=neu_only_df, x='neu_signal', bins=50, alpha=0.5, label='neu Signal', color='blue')
    sns.histplot(data=neu_only_df, x='nsc_signal', bins=50, alpha=0.5, label='nsc Signal', color='red')
    if title is None:
        _title = 'Distribution of Signals'
    else:
        _title = f'Distribution of Signals\n{title}'
    plt.title(_title)
    plt.xlabel('Signal Value')
    plt.ylabel('Count')
    plt.legend()

    # Plot overlapped length distributions
    plt.subplot(1, 3, 2)
    sns.histplot(data=neu_only_df, x='region_length', bins=50, alpha=0.5, label='Region Length', color='blue')
    sns.histplot(data=neu_only_df, x='cpg_length', bins=50, alpha=0.5, label='CpG Length', color='red')
    if title is None:
        _title = 'Distribution of Lengths'
    else:
        _title = f'Distribution of Lengths\n{title}'
    plt.title(_title)
    plt.xlabel('Length')
    plt.ylabel('Count')
    plt.legend()

    # Count peaks per region
    plt.subplot(1, 3, 3)
    neu_peak_counts = neu_only_df['neu_replicates_with_peaks'].value_counts()
    nsc_peak_counts = neu_only_df['nsc_replicates_with_peaks'].value_counts()

    x = np.arange(max(max(neu_peak_counts.index), max(nsc_peak_counts.index)) + 1)
    width = 0.35

    # Create bars and store them to add labels
    neu_bars = plt.bar(x - width/2, [neu_peak_counts.get(i, 0) for i in x], width, label='neu', color='blue', alpha=0.5)
    nsc_bars = plt.bar(x + width/2, [nsc_peak_counts.get(i, 0) for i in x], width, label='nsc', color='red', alpha=0.5)

    # Add count labels on top of each bar
    for bars in [neu_bars, nsc_bars]:
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}',
                    ha='center', va='bottom')
    if title is None:
        _title = 'Number of Replicates with Peaks'
    else:
        _title = f'Number of Replicates with Peaks\n{title}'
    plt.title(_title)
    plt.xlabel('Number of Replicates')
    plt.ylabel('Count')
    plt.legend()

    plt.xticks(x, x.astype(int))
    plt.tight_layout()
    plt.show()



def plot_nsc_only_distributions(nsc_only_df, title=None):
    # Distribution plots
    plt.figure(figsize=(15, 5))

    # Plot overlapped signal distributions
    plt.subplot(1, 3, 1)
    sns.histplot(data=nsc_only_df, x='nsc_signal', bins=50, alpha=0.5, label='nsc Signal', color='blue')
    sns.histplot(data=nsc_only_df, x='neu_signal', bins=50, alpha=0.5, label='neu Signal', color='red')
    if title is None:
        _title = 'Distribution of Signals'
    else:
        _title = f'Distribution of Signals\n{title}'
    plt.title(_title)
    plt.xlabel('Signal Value')
    plt.ylabel('Count')
    plt.legend()

    # Plot overlapped length distributions
    plt.subplot(1, 3, 2)
    sns.histplot(data=nsc_only_df, x='region_length', bins=50, alpha=0.5, label='Region Length', color='blue')
    sns.histplot(data=nsc_only_df, x='cpg_length', bins=50, alpha=0.5, label='CpG Length', color='red')
    if title is None:
        _title = 'Distribution of Lengths'
    else:
        _title = f'Distribution of Lengths\n{title}'
    plt.title(_title)
    plt.xlabel('Length')
    plt.ylabel('Count')
    plt.legend()

    # Count peaks per region
    plt.subplot(1, 3, 3)
    neu_peak_counts = nsc_only_df['neu_replicates_with_peaks'].value_counts()
    nsc_peak_counts = nsc_only_df['nsc_replicates_with_peaks'].value_counts()

    x = np.arange(max(max(neu_peak_counts.index), max(nsc_peak_counts.index)) + 1)
    width = 0.35

    # Create bars and store them to add labels
    neu_bars = plt.bar(x - width/2, [neu_peak_counts.get(i, 0) for i in x], width, label='neu', color='blue', alpha=0.5)
    nsc_bars = plt.bar(x + width/2, [nsc_peak_counts.get(i, 0) for i in x], width, label='nsc', color='red', alpha=0.5)

    # Add count labels on top of each bar
    for bars in [nsc_bars, nsc_bars]:
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}',
                    ha='center', va='bottom')
    if title is None:
        _title = 'Number of Replicates with Peaks'
    else:
        _title = f'Number of Replicates with Peaks\n{title}'
    plt.title(_title)
    plt.xlabel('Number of Replicates')
    plt.ylabel('Count')
    plt.legend()

    plt.xticks(x, x.astype(int))
    plt.tight_layout()
    plt.show()
