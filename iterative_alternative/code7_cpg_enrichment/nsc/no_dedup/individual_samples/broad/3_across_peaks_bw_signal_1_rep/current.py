# %% [markdown]
# # Environment

# %%
# import pandas as pd
# import yaml
# import os
# import sys
# import importlib
# import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt

# # Load config
# CONFIG_PATH = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/code7_cpg_enrichment/config.yaml"

# # Specify which configuration to use (1-based index)
# # Change this to select different configurations
# CONFIG_NUMBER = 3

# with open(CONFIG_PATH, 'r') as f:
#     # Load all documents from the YAML file
#     configs = list(yaml.safe_load_all(f))

# # Select the specific configuration (subtract 1 because list is 0-based)
# config = configs[CONFIG_NUMBER - 1]

# # Get values from the selected config
# BASE_DIR = config['base_dir']
# RUN_NAME = config['run_name']
# CELL_TYPE = config['cell_type']
# ALIGNMENT_TYPE = config['alignment_type']
# PEAKS_TYPE = config['peaks_type']

# ENRICHMENT_FILE = f"{BASE_DIR}/results/{ALIGNMENT_TYPE}/cpg_enrichment/{CELL_TYPE}/{PEAKS_TYPE}/{RUN_NAME}/cpg_enrichment_parallel.csv"
# print(ENRICHMENT_FILE)

# OUTPUT_LISTS_PATH = f"{BASE_DIR}/results/{ALIGNMENT_TYPE}/cpg_enrichment/{CELL_TYPE}/{PEAKS_TYPE}/{RUN_NAME}/lists"
# print(OUTPUT_LISTS_PATH)
# os.makedirs(OUTPUT_LISTS_PATH, exist_ok=True)

# # Set pandas display options to show all columns without wrapping
# pd.set_option('display.max_columns', None)  # Show all columns
# pd.set_option('display.width', None)        # Don't wrap wide DataFrames
# pd.set_option('display.max_colwidth', None) # Don't truncate column contents

# sys.path.append("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/code7_cpg_enrichment")
# try:
#     import functions
#     importlib.reload(functions)
#     from functions import *
# except ModuleNotFoundError:
#     print("Warning: Could not import functions module. Make sure functions.py exists in the code7_cpg_enrichment directory.")


# %%
import pandas as pd
import yaml
import os
import sys
import importlib
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load config
# CONFIG_PATH = "/home/michal/Github/SRF_MeCP2_CUTandTAG/iterative_alternative/code7_cpg_enrichment/config.yaml"
CONFIG_PATH = "D:/Github/SRF_MeCP2_cut_tag/iterative_alternative/code7_cpg_enrichment/config.yaml"
# Specify which configuration to use (1-based index)
# Change this to select different configurations
CONFIG_NUMBER = 3

with open(CONFIG_PATH, 'r') as f:
    # Load all documents from the YAML file
    configs = list(yaml.safe_load_all(f))

# Select the specific configuration (subtract 1 because list is 0-based)
config = configs[CONFIG_NUMBER - 1]

# Get values from the selected config
BASE_DIR = config['base_dir']
RUN_NAME = config['run_name']
CELL_TYPE = config['cell_type']
ALIGNMENT_TYPE = config['alignment_type']
PEAKS_TYPE = config['peaks_type']

ENRICHMENT_FILE = f"{BASE_DIR}/results/{ALIGNMENT_TYPE}/cpg_enrichment/{CELL_TYPE}/{PEAKS_TYPE}/{RUN_NAME}/cpg_enrichment_parallel.csv"
print(ENRICHMENT_FILE)

OUTPUT_LISTS_PATH = f"{BASE_DIR}/results/{ALIGNMENT_TYPE}/cpg_enrichment/{CELL_TYPE}/{PEAKS_TYPE}/{RUN_NAME}/lists"
print(OUTPUT_LISTS_PATH)
os.makedirs(OUTPUT_LISTS_PATH, exist_ok=True)

# Set pandas display options to show all columns without wrapping
pd.set_option('display.max_columns', None)  # Show all columns
pd.set_option('display.width', None)        # Don't wrap wide DataFrames
pd.set_option('display.max_colwidth', None) # Don't truncate column contents

sys.path.append("/home/michal/Github/SRF_MeCP2_CUTandTAG/iterative_alternative/code7_cpg_enrichment")
try:
    import functions
    importlib.reload(functions)
    from functions import *
except ModuleNotFoundError:
    print("Warning: Could not import functions module. Make sure functions.py exists in the code7_cpg_enrichment directory.")


# %% [markdown]
# # Load Data

# %%
# Load the CpG enrichment results
cpg_enrichment_df = pd.read_csv(ENRICHMENT_FILE)

# %%
# Display first few rows and basic info
print("DataFrame shape:", cpg_enrichment_df.shape)
print("\nFirst few rows:")
display(cpg_enrichment_df.head())
print("\nColumn names:")
print(cpg_enrichment_df.columns.tolist())

# %% [markdown]
# **Columns explanation:**
# 
# 1. Signal vs Peaks columns:
# 
# - `endo_replicates_with_signal`: Counts how many replicates have a non-zero signal value when measuring the bigWig signal in the region.
# - `endo_replicates_with_peaks`: Counts how many replicates have an overlapping peak in the broadPeak files.
# 
# So, e.g. it's possible to have `endo_replicates_with_peaks=0` but `endo_replicates_with_signal=2` because:
# - Peaks represent regions that passed the peak caller's statistical threshold for significance
# - Signal values represent the raw enrichment data before any statistical filtering
# - So you can have detectable signal in a region that wasn't strong/consistent enough to be called as a peak
# 
# 2. Scores vs Signals:
# 
# - `endo_peak_scores_by_rep`: Contains the `signalValue` scores from the broadPeak files for each replicate that has a peak overlapping the region. These scores are peak caller-specific enrichment metrics.
# - `endo_replicate_signals`: Contains the actual signal values extracted from the bigWig files for each replicate in that region. These are the raw signal values.
# 
# `endo_peak_scores_by_rep` is NaN when there are no peaks overlapping the region in any replicate (`endo_replicates_with_peaks=0`). 
# 
# **Handling multiple replicates:**
# 
# 1. For Signals (from bigWig files):
# - Each replicate's signal is obtained independently using `get_signal_from_bigwig()`
# - The signals are stored as individual values in `endo_replicate_signals` as a comma-separated string
# - The final `endo_signal` used for enrichment calculations is the **mean across all replicates**
# 
# 2. For Peak Scores:
# - Peak scores are stored in `endo_peak_scores_by_rep` using a specific format:
#   - Scores from different peaks within the same replicate are **comma-separated**
#   - Different replicates are **semicolon-separated**
#   - If a replicate has no peaks, it's simply not included in the string
# - **There's no averaging of peak scores**
# 
# For example:
# - If replicate 1 has two peaks with scores 5.0 and 6.0, and replicate 2 has one peak with score 4.0, while replicate 3 has no peaks:
#   - `endo_peak_scores_by_rep` would be: `"5.0, 6.0; 4.0"`
#   - `endo_replicates_with_peaks` would be `2`
#   - If all replicates had signals of `2.0`, `3.0`, and `1.0`:
#     - `endo_replicate_signals` would be `"2.0, 3.0, 1.0"`
#     - `endo_signal` would be `2.0` (the mean)
#     - `endo_replicates_with_signal` would be `3`

# %%
# Create single plot
plt.figure(figsize=(10, 6))

# Get data for plotting
x = np.arange(4)  # 0 through 3 replicates
exo_counts = [len(cpg_enrichment_df[cpg_enrichment_df['exo_replicates_with_peaks'] == i]) for i in range(4)]
endo_counts = [len(cpg_enrichment_df[cpg_enrichment_df['endo_replicates_with_peaks'] == i]) for i in range(4)]

# Plot bars side by side
width = 0.35
plt.bar(x - width/2, exo_counts, width, label=f'Exo (n={sum(exo_counts)})', color='red', alpha=0.6, edgecolor='black')
plt.bar(x + width/2, endo_counts, width, label=f'Endo (n={sum(endo_counts)})', color='green', alpha=0.6, edgecolor='black')

# Add value labels on top of each bar
for i in range(len(x)):
    plt.text(x[i] - width/2, exo_counts[i], str(exo_counts[i]), 
             ha='center', va='bottom')
    plt.text(x[i] + width/2, endo_counts[i], str(endo_counts[i]),
             ha='center', va='bottom')

plt.title('Distribution of Replicates with Peaks')
plt.xlabel('Number of Replicates')
plt.ylabel('Count')
plt.xticks(x)
plt.legend()

plt.tight_layout()
plt.show()

# %%
# Sort the DataFrame by enrichment value in descending order
cpg_enrichment_df_sorted = cpg_enrichment_df.sort_values(by='enrichment', ascending=False)

# %% [markdown]
# # Split data based on binding type

# %% [markdown]
# ## By signal value

# %% [markdown]
# Non zero signal in minumum 2 replicates

# %%
# Split data based on binding type
exo_only_df_by_signal = cpg_enrichment_df_sorted[cpg_enrichment_df_sorted['binding_type'] == 'exo_only']
endo_only_df_by_signal = cpg_enrichment_df_sorted[cpg_enrichment_df_sorted['binding_type'] == 'endo_only'] 
both_df_by_signal = cpg_enrichment_df_sorted[cpg_enrichment_df_sorted['binding_type'] == 'both']

# Print sizes of each group
print(f"Number of CpG islands bound by exo only: {len(exo_only_df_by_signal)}")
print(f"Number of CpG islands bound by endo only: {len(endo_only_df_by_signal)}")
print(f"Number of CpG islands bound by both: {len(both_df_by_signal)}")


# %% [markdown]
# ## By peaks number

# %% [markdown]
# Minimum 2 replicates with peaks

# %%
# Split data based on binding type
exo_only_df_by_peaks = cpg_enrichment_df_sorted[cpg_enrichment_df_sorted['binding_type_by_peaks'] == 'exo_only']
endo_only_df_by_peaks = cpg_enrichment_df_sorted[cpg_enrichment_df_sorted['binding_type_by_peaks'] == 'endo_only'] 
both_df_by_peaks = cpg_enrichment_df_sorted[cpg_enrichment_df_sorted['binding_type_by_peaks'] == 'both']

# Print sizes of each group
print(f"Number of CpG islands bound by exo only: {len(exo_only_df_by_peaks)}")
print(f"Number of CpG islands bound by endo only: {len(endo_only_df_by_peaks)}")
print(f"Number of CpG islands bound by both: {len(both_df_by_peaks)}")


# %% [markdown]
# # Data analysis

# %%
NSC_COLOR = "#201cf4"  # Blue
NEU_COLOR = "#f04c3c"  # Orange/Red

# %% [markdown]
# ## Both: enrichment distribution

# %%
print("Summary statistics of enrichment values selected by signal:")
print(both_df_by_signal['enrichment'].describe())

print("\nSummary statistics of enrichment values selected by peaks:")
print(both_df_by_peaks['enrichment'].describe())

# %%
range_min_signal, range_max_signal = plot_enrichment_distribution(both_df_by_signal, title=None, color=NSC_COLOR, alpha=0.7, enrichment_lines=[2])
range_min_peaks, range_max_peaks = plot_enrichment_distribution(both_df_by_peaks, title="selected by peaks", color='#2196F3', alpha=0.7, enrichment_lines=[1])


# %% [markdown]
# ## Compare with NEU

# %%
neu_enrichment_df = pd.read_csv("/home/michal/Github/SRF_MeCP2_CUTandTAG/iterative_alternative/code7_cpg_enrichment/neu/no_dedup/individual_samples/broad/3_across_peaks_bw_signal_1_rep/both_df_by_signal.csv")
neu_enrichment_df.head()

# %%
range_min_signal, range_max_signal = plot_enrichment_distribution(neu_enrichment_df, title=None, color=NEU_COLOR, alpha=0.7, enrichment_lines=[2])

# %%
def plot_enrichment_distribution_mod(df_list, labels, factor=3, title=None, 
                                     colors=None, alpha=0.5, enrichment_lines=[2]):
    """
    Create overlaid histograms of enrichment values for regions bound by both Exo and Endo.
    
    Args:
        df_list: List of DataFrames containing enrichment values
        labels: List of labels for each DataFrame
        factor: Factor to multiply IQR by for determining plot range (default=3)
        title: Optional custom title for the plot
        colors: List of colors for each histogram (default: [NEU_COLOR, NSC_COLOR])
        alpha: Transparency value for histograms (default: 0.5)
        enrichment_lines: List of enrichment values at which to draw vertical lines (default: [1])
    """
    plt.figure(figsize=(10, 6))

    # Determine global range based on all data distributions
    all_enrichments = []
    for df in df_list:
        all_enrichments.extend(df['enrichment'])
    
    q1 = np.quantile(all_enrichments, 0.25)
    q3 = np.quantile(all_enrichments, 0.75)
    iqr = q3 - q1
    range_min = max(0, q1 - factor * iqr)  # Don't go below 0 for enrichment
    range_max = q3 + factor * iqr

    if colors is None:
        colors = [NEU_COLOR, NSC_COLOR]  # Assign default colors


    # Plot overlaid histograms
    for i, df in enumerate(df_list):
        plt.hist(df['enrichment'], bins=50, alpha=alpha,
                 label=labels[i], density=False, range=(range_min, range_max), color=colors[i])

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

    # Add vertical lines at specified enrichment values
    for enrichment_val in enrichment_lines:
        plt.axvline(x=enrichment_val, color='red', linestyle='--', linewidth=2,
                    label=f'Enrichment = {enrichment_val}')

    # Set x-axis limits based on the calculated data range
    plt.xlim(range_min, range_max)

    # Add legend with better positioning
    plt.legend(loc='upper right')

    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    plt.show()
    
    return range_min, range_max


range_min_signal, range_max_signal = plot_enrichment_distribution_mod([neu_enrichment_df, both_df_by_signal], 
                                                                     labels=["NEU", "NSC"], 
                                                                     title=None)


# %%
def plot_signal_distribution_mod(df_list, labels, signal_type, factor=3, title=None, colors=None, alpha=0.5, percentile=0.85):
    """
    Create overlaid histograms of signal values for regions bound by both Exo and Endo.
    
    Args:
        df_list: List of DataFrames containing signal values
        labels: List of labels for each DataFrame
        signal_type: 'exo_signal' or 'endo_signal'
        factor: Factor to multiply IQR by for determining plot range (default=3)
        title: Optional custom title for the plot
        colors: List of colors for each histogram (default=None, uses NSC_COLOR and NEU_COLOR)
        alpha: Transparency value for histograms (default=0.5)
        percentile: Percentile value for vertical lines (default=0.95, i.e., 95th percentile)
    """
    plt.figure(figsize=(10, 6))

    # Determine global range based on all data distributions
    all_signals = []
    for df in df_list:
        all_signals.extend(df[signal_type])
    
    q1 = np.quantile(all_signals, 0.25)
    q3 = np.quantile(all_signals, 0.75)
    iqr = q3 - q1
    range_min = max(0, q1 - factor * iqr)  # Don't go below 0 for signal
    range_max = q3 + factor * iqr

    if colors is None:
        colors = [NEU_COLOR, NSC_COLOR]  # Assign default colors


    # Plot overlaid histograms
    for i, df in enumerate(df_list):
        plt.hist(df[signal_type], bins=50, alpha=alpha,
                 label=labels[i], density=False, range=(range_min, range_max), color=colors[i])

        # Add vertical lines at specified enrichment values
        enrichment_percentile = np.quantile(df[signal_type], percentile)
        plt.axvline(x=enrichment_percentile, color=colors[i], linestyle='--', linewidth=2,
                    label=f'{labels[i]} {percentile*100:.0f}th percentile = {enrichment_percentile:.2f}')

    # Add grid for better readability
    plt.grid(True, alpha=0.3)

    # Improve axis labels and title
    plt.xlabel(f'{signal_type} Score', fontsize=12)
    plt.ylabel('Number of Regions', fontsize=12)
    if title is None:
        title = f'Distribution of {signal_type} Values'
    else:
        title = f'Distribution of {signal_type} Values\n{title}'
    plt.title(title, fontsize=14, pad=15)


    # Set x-axis limits based on the calculated data range
    plt.xlim(range_min, range_max)

    # Add legend with better positioning
    plt.legend(loc='upper right')

    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    plt.show()
    
    return range_min, range_max

range_min_signal, range_max_signal = plot_signal_distribution_mod([neu_enrichment_df, both_df_by_signal], labels=["NEU", "NSC"], signal_type='exo_signal', title="selected by signal")


# %%
range_min_signal, range_max_signal = plot_signal_distribution_mod([neu_enrichment_df, both_df_by_signal], labels=["NEU", "NSC"], signal_type='endo_signal', title="selected by signal")

# %% [markdown]
# ## Both: enrichment outliers

# %%
outliers_df_signal = both_df_by_signal[both_df_by_signal['enrichment'] > range_max_signal].copy()
print(f"\nNumber of outliers selected by signal(enrichment > {range_max_signal}):", len(outliers_df_signal))

outliers_df_peaks = both_df_by_peaks[both_df_by_peaks['enrichment'] > range_max_peaks].copy()
print(f"\nNumber of outliers selected by peaks(enrichment > {range_max_peaks}):", len(outliers_df_peaks))

# %%
outliers_df_signal.head()

# %%
print_outlier_groups(outliers_df_signal, "signal")

# %%
n_rows = len(outliers_df_signal)

print("\nFirst 5 outliers selected by signal:")
for _, row in outliers_df_signal[:5].iterrows():
    print(f"{row['chr']}:{row['start']:,}-{row['end']:,}".ljust(50) + f"enrichment: {row['enrichment']}")

if n_rows > 10:
    mid_start = n_rows//2 - 2
    print("\nMiddle 5 outliers selected by signal:")
    for _, row in outliers_df_signal[mid_start:mid_start+5].iterrows():
        print(f"{row['chr']}:{row['start']:,}-{row['end']:,}".ljust(50) + f"enrichment: {row['enrichment']}")

print("\nLast 5 outliers selected by signal:")
for _, row in outliers_df_signal[-5:].iterrows():
    print(f"{row['chr']}:{row['start']:,}-{row['end']:,}".ljust(50) + f"enrichment: {row['enrichment']}")

# %%
print_outlier_groups(outliers_df_peaks, "peaks")

# %%
plot_enrichment_boxplot(both_df_by_peaks, title="selected by peaks")
plot_enrichment_boxplot(both_df_by_signal, title="selected by signal")

# %%
outliers_df_signal["endo_replicates_with_peaks"].value_counts()

signal_endo_zero_peaks  = outliers_df_signal[outliers_df_signal["endo_replicates_with_peaks"] == 0]
signal_endo_one_peaks   = outliers_df_signal[outliers_df_signal["endo_replicates_with_peaks"] == 1]
signal_endo_two_peaks   = outliers_df_signal[outliers_df_signal["endo_replicates_with_peaks"] == 2]
signal_endo_three_peaks = outliers_df_signal[outliers_df_signal["endo_replicates_with_peaks"] == 3]

outliers_df_peaks["endo_replicates_with_peaks"].value_counts()

peaks_endo_zero_peaks  = outliers_df_peaks[outliers_df_peaks["endo_replicates_with_peaks"] == 0]
peaks_endo_one_peaks   = outliers_df_peaks[outliers_df_peaks["endo_replicates_with_peaks"] == 1]
peaks_endo_two_peaks   = outliers_df_peaks[outliers_df_peaks["endo_replicates_with_peaks"] == 2]
peaks_endo_three_peaks = outliers_df_peaks[outliers_df_peaks["endo_replicates_with_peaks"] == 3]

# %%
# Classify and plot for signal and peaks outliers
peaks_zero_peaks_signal, peaks_one_peaks_signal, peaks_two_peaks_signal, peaks_three_peaks_signal = plot_by_peaks(outliers_df_signal, peaks_column='endo_replicates_with_peaks', title="Endo peaks - selected by signal")
peaks_zero_peaks_peaks, peaks_one_peaks_peaks, peaks_two_peaks_peaks, peaks_three_peaks_peaks = plot_by_peaks(outliers_df_peaks, peaks_column='endo_replicates_with_peaks', title="Endo peaks - selected by peaks")

# %%
print_outlier_groups(peaks_two_peaks_peaks, "peaks")

# %%
print_outlier_groups(peaks_two_peaks_signal, "signal")

# %% [markdown]
# ## Regions length distribution

# %% [markdown]
# ### Region length - defined by the outermost peaks coordinates

# %%
plot_region_length_comparison(both_df_by_peaks, both_df_by_signal, region_length_col='region_length')

# %%
both_df_sorted_by_region_length_signal = both_df_by_signal.sort_values(by='region_length', ascending=False)
both_df_sorted_by_region_length_signal.head()

# %%
both_df_sorted_by_region_length_peaks = both_df_by_peaks.sort_values(by='region_length', ascending=False)
both_df_sorted_by_region_length_peaks.head()

# %%
print_outlier_groups(both_df_sorted_by_region_length_peaks, "peaks")

# %%
print_outlier_groups(both_df_sorted_by_region_length_signal, "signal")

# %% [markdown]
# ### CpG length - defined by the CpG coordinates

# %%
plot_region_length_comparison(both_df_by_peaks, both_df_by_signal, region_length_col='cpg_length')

# %%
both_df_sorted_by_cpg_length_peaks = both_df_by_peaks.sort_values(by='cpg_length', ascending=False)
both_df_sorted_by_cpg_length_peaks.head()
print_outlier_groups(both_df_sorted_by_cpg_length_peaks, "peaks")

# %%
both_df_sorted_by_cpg_length_signal = both_df_by_signal.sort_values(by='cpg_length', ascending=False)
both_df_sorted_by_cpg_length_signal.head()
print_outlier_groups(both_df_sorted_by_cpg_length_signal, "signal")

# %% [markdown]
# # Exo only

# %%
exo_only_df_by_peaks.head()

# %%
# Generate summary statistics for numeric columns, ignoring inf values
# Replace inf values with NaN before calculating statistics
numeric_cols = ['exo_signal', 'endo_signal', 'enrichment', 'region_length', 'cpg_length', 'pvalue']

stats_df = exo_only_df_by_peaks[numeric_cols].replace([np.inf, -np.inf], np.nan).describe()
print("Exo only - by peaks, size:", len(exo_only_df_by_peaks))
stats_df

# %%
stats_df = exo_only_df_by_signal[numeric_cols].replace([np.inf, -np.inf], np.nan).describe()
print("Exo only - by signal, size:", len(exo_only_df_by_signal))
stats_df

# %%
print_across_distribution(exo_only_df_by_peaks, selection_type="peaks")

# %%
print_across_distribution(exo_only_df_by_signal, selection_type="signal")

# %%
plot_exo_only_distributions(exo_only_df_by_peaks, title="by peaks")
plot_exo_only_distributions(exo_only_df_by_signal, title="by signal")

# %% [markdown]
# # Endo only

# %%
endo_only_df_by_peaks.head()

# %%
# Generate summary statistics for numeric columns, ignoring inf values
# Replace inf values with NaN before calculating statistics
numeric_cols = ['exo_signal', 'endo_signal', 'enrichment', 'region_length', 'cpg_length', 'pvalue']

stats_df = endo_only_df_by_peaks[numeric_cols].replace([np.inf, -np.inf], np.nan).describe()
print("endo only - by peaks, size:", len(endo_only_df_by_peaks))
stats_df

# %%
stats_df = endo_only_df_by_signal[numeric_cols].replace([np.inf, -np.inf], np.nan).describe()
print("endo only - by signal, size:", len(endo_only_df_by_signal))
stats_df

# %%
print_across_distribution(endo_only_df_by_peaks, selection_type="peaks")

# %%
print_across_distribution(endo_only_df_by_signal, selection_type="signal")

# %%
plot_endo_only_distributions(endo_only_df_by_peaks, title="by peaks")
plot_endo_only_distributions(endo_only_df_by_signal, title="by signal")

# %% [markdown]
# # Create output files

# %% [markdown]
# ## UP in Exo

# %%
up_enriched_signal_1 = both_df_by_signal[(both_df_by_signal['enrichment'] > 1) & (both_df_by_signal['enrichment'] < range_max_signal)].copy()
up_enriched_signal_1_5 = both_df_by_signal[(both_df_by_signal['enrichment'] > 1.5) & (both_df_by_signal['enrichment'] < range_max_signal)].copy()    
up_enriched_signal_2 = both_df_by_signal[(both_df_by_signal['enrichment'] > 2) & (both_df_by_signal['enrichment'] < range_max_signal)].copy()
up_enriched_peaks_1 = both_df_by_peaks[(both_df_by_peaks['enrichment'] > 1) & (both_df_by_peaks['enrichment'] < range_max_peaks)].copy()
up_enriched_peaks_1_5 = both_df_by_peaks[(both_df_by_peaks['enrichment'] > 1.5) & (both_df_by_peaks['enrichment'] < range_max_peaks)].copy()    
up_enriched_peaks_2 = both_df_by_peaks[(both_df_by_peaks['enrichment'] > 2) & (both_df_by_peaks['enrichment'] < range_max_peaks)].copy()

# %%
print("up_enriched_signal_1.shape:",   up_enriched_signal_1.shape)
print("up_enriched_signal_1_5.shape:", up_enriched_signal_1_5.shape)
print("up_enriched_signal_2.shape:",   up_enriched_signal_2.shape)
print("up_enriched_peaks_1.shape:",    up_enriched_peaks_1.shape)
print("up_enriched_peaks_1_5.shape:",  up_enriched_peaks_1_5.shape)
print("up_enriched_peaks_2.shape:",    up_enriched_peaks_2.shape)

# %% [markdown]
# ## UP in Endo

# %%
down_enriched_signal_1 = both_df_by_signal[(both_df_by_signal['enrichment'] < 1.0 ) & (both_df_by_signal['enrichment'] < range_max_signal)].copy()
down_enriched_signal_08 = both_df_by_signal[(both_df_by_signal['enrichment'] < 0.8 ) & (both_df_by_signal['enrichment'] < range_max_signal)].copy()
down_enriched_signal_05 = both_df_by_signal[(both_df_by_signal['enrichment'] < 0.5 ) & (both_df_by_signal['enrichment'] < range_max_signal)].copy()
down_enriched_peaks_1 = both_df_by_peaks[(both_df_by_peaks['enrichment'] < 1.0) & (both_df_by_peaks['enrichment'] < range_max_peaks)].copy()
down_enriched_peaks_08 = both_df_by_peaks[(both_df_by_peaks['enrichment'] < 0.8) & (both_df_by_peaks['enrichment'] < range_max_peaks)].copy()
down_enriched_peaks_05 = both_df_by_peaks[(both_df_by_peaks['enrichment'] < 0.5) & (both_df_by_peaks['enrichment'] < range_max_peaks)].copy()


# %%
print("down_enriched_signal_1.shape:",  down_enriched_signal_1.shape)
print("down_enriched_signal_08.shape:", down_enriched_signal_08.shape)
print("down_enriched_signal_05.shape:", down_enriched_signal_05.shape)
print("down_enriched_peaks_1.shape:",   down_enriched_peaks_1.shape)
print("down_enriched_peaks_08.shape:",  down_enriched_peaks_08.shape)
print("down_enriched_peaks_05.shape:",  down_enriched_peaks_05.shape)

# %% [markdown]
# ## Save output files

# %%
up_enriched_signal_1.to_csv(f'{OUTPUT_LISTS_PATH}/up_enriched_signal_1.csv', index=False)
up_enriched_signal_1_5.to_csv(f'{OUTPUT_LISTS_PATH}/up_enriched_signal_1_5.csv', index=False)
up_enriched_signal_2.to_csv(f'{OUTPUT_LISTS_PATH}/up_enriched_signal_2.csv', index=False)
up_enriched_peaks_1.to_csv(f'{OUTPUT_LISTS_PATH}/up_enriched_peaks_1.csv', index=False)
up_enriched_peaks_1_5.to_csv(f'{OUTPUT_LISTS_PATH}/up_enriched_peaks_1_5.csv', index=False)
up_enriched_peaks_2.to_csv(f'{OUTPUT_LISTS_PATH}/up_enriched_peaks_2.csv', index=False)

down_enriched_signal_1.to_csv(f'{OUTPUT_LISTS_PATH}/down_enriched_signal_1.csv', index=False)
down_enriched_signal_08.to_csv(f'{OUTPUT_LISTS_PATH}/down_enriched_signal_08.csv', index=False)
down_enriched_signal_05.to_csv(f'{OUTPUT_LISTS_PATH}/down_enriched_signal_05.csv', index=False)
down_enriched_peaks_1.to_csv(f'{OUTPUT_LISTS_PATH}/down_enriched_peaks_1.csv', index=False)
down_enriched_peaks_08.to_csv(f'{OUTPUT_LISTS_PATH}/down_enriched_peaks_08.csv', index=False)
down_enriched_peaks_05.to_csv(f'{OUTPUT_LISTS_PATH}/down_enriched_peaks_05.csv', index=False)

exo_only_df_by_signal.to_csv(f'{OUTPUT_LISTS_PATH}/exo_only_df_by_signal.csv', index=False)
endo_only_df_by_signal.to_csv(f'{OUTPUT_LISTS_PATH}/endo_only_df_by_signal.csv', index=False)
exo_only_df_by_peaks.to_csv(f'{OUTPUT_LISTS_PATH}/exo_only_df_by_peaks.csv', index=False)
endo_only_df_by_peaks.to_csv(f'{OUTPUT_LISTS_PATH}/endo_only_df_by_peaks.csv', index=False)

# %% [markdown]
# # Exo min expression filtering

# %%
up_enriched_signal_1_exo_over_20 = up_enriched_signal_1[up_enriched_signal_1['exo_signal'] >= 20]
up_enriched_signal_1_5_exo_over_20 = up_enriched_signal_1_5[up_enriched_signal_1_5['exo_signal'] >= 20]
up_enriched_signal_2_exo_over_20 = up_enriched_signal_2[up_enriched_signal_2['exo_signal'] >= 20]
up_enriched_peaks_1_exo_over_20 = up_enriched_peaks_1[up_enriched_peaks_1['exo_signal'] >= 20]
up_enriched_peaks_1_5_exo_over_20 = up_enriched_peaks_1_5[up_enriched_peaks_1_5['exo_signal'] >= 20]
up_enriched_peaks_2_exo_over_20 = up_enriched_peaks_2[up_enriched_peaks_2['exo_signal'] >= 20]

down_enriched_signal_1_exo_over_20 = down_enriched_signal_1[down_enriched_signal_1['exo_signal'] >= 20]
down_enriched_signal_08_exo_over_20 = down_enriched_signal_08[down_enriched_signal_08['exo_signal'] >= 20]
down_enriched_signal_05_exo_over_20 = down_enriched_signal_05[down_enriched_signal_05['exo_signal'] >= 20]
down_enriched_peaks_1_exo_over_20 = down_enriched_peaks_1[down_enriched_peaks_1['exo_signal'] >= 20]
down_enriched_peaks_08_exo_over_20 = down_enriched_peaks_08[down_enriched_peaks_08['exo_signal'] >= 20]
down_enriched_peaks_05_exo_over_20 = down_enriched_peaks_05[down_enriched_peaks_05['exo_signal'] >= 20]

exo_only_df_by_signal_exo_over_20 = exo_only_df_by_signal[exo_only_df_by_signal['exo_signal'] >= 20]
endo_only_df_by_signal_exo_over_20 = endo_only_df_by_signal[endo_only_df_by_signal['exo_signal'] >= 20]
exo_only_df_by_peaks_exo_over_20 = exo_only_df_by_peaks[exo_only_df_by_peaks['exo_signal'] >= 20]
endo_only_df_by_peaks_exo_over_20 = endo_only_df_by_peaks[endo_only_df_by_peaks['exo_signal'] >= 20]


up_enriched_signal_1_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/up_enriched_signal_1_exo_over_20.csv', index=False)
up_enriched_signal_1_5_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/up_enriched_signal_1_5_exo_over_20.csv', index=False)
up_enriched_signal_2_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/up_enriched_signal_2_exo_over_20.csv', index=False)
up_enriched_peaks_1_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/up_enriched_peaks_1_exo_over_20.csv', index=False)
up_enriched_peaks_1_5_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/up_enriched_peaks_1_5_exo_over_20.csv', index=False)
up_enriched_peaks_2_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/up_enriched_peaks_2_exo_over_20.csv', index=False)

down_enriched_signal_1_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/down_enriched_signal_1_exo_over_20.csv', index=False)
down_enriched_signal_08_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/down_enriched_signal_08_exo_over_20.csv', index=False)
down_enriched_signal_05_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/down_enriched_signal_05_exo_over_20.csv', index=False)
down_enriched_peaks_1_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/down_enriched_peaks_1_exo_over_20.csv', index=False)
down_enriched_peaks_08_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/down_enriched_peaks_08_exo_over_20.csv', index=False)
down_enriched_peaks_05_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/down_enriched_peaks_05_exo_over_20.csv', index=False)

exo_only_df_by_signal_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/exo_only_df_by_signal_exo_over_20.csv', index=False)
endo_only_df_by_signal_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/endo_only_df_by_signal_exo_over_20.csv', index=False)
exo_only_df_by_peaks_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/exo_only_df_by_peaks_exo_over_20.csv', index=False)
endo_only_df_by_peaks_exo_over_20.to_csv(f'{OUTPUT_LISTS_PATH}/endo_only_df_by_peaks_exo_over_20.csv', index=False)


