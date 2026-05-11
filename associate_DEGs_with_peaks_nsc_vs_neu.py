import pandas as pd
import numpy as np
import os

# Define file paths
peaks_file = r"D:\Github\SRF_MeCP2_cut_tag\iterative_alternative\results\no_dedup\cpg_enrichment\neu_vs_nsc_endo\peaks_annotation\up_enriched_signal_2_tss_peaks.csv"
dea_file = r"D:\Github\SRF_MeCP2_cut_tag\data\DEA_NEU.csv"

# --- Determine output file path ---
peaks_dir = os.path.dirname(peaks_file)
peaks_basename = os.path.basename(peaks_file)
peaks_name_without_ext = os.path.splitext(peaks_basename)[0]
output_file = os.path.join(peaks_dir, f"{peaks_name_without_ext}_with_DEGs.csv")

print(f"Reading peaks file: {peaks_file}")
print(f"Reading DEA file: {dea_file}")
print(f"Output will be saved to: {output_file}")

# --- Read the CSV files ---
try:
    peaks_df = pd.read_csv(peaks_file)
    print(f"Successfully read peaks file. Shape: {peaks_df.shape}")
    print("Peaks file columns:", peaks_df.columns.tolist())
    print("First 5 rows of peaks data:\n", peaks_df.head().to_string())

    dea_df = pd.read_csv(dea_file)
    print(f"\nSuccessfully read DEA file. Shape: {dea_df.shape}")
    print("DEA file columns:", dea_df.columns.tolist())
    print("First 5 rows of DEA data:\n", dea_df.head().to_string())

except FileNotFoundError as e:
    print(f"Error: File not found - {e}")
    exit()
except Exception as e:
    print(f"Error reading CSV files: {e}")
    exit()

# --- Prepare DEA data ---
# Ensure 'gene' column exists and handle potential missing values if necessary
if 'gene' not in dea_df.columns:
    print("Error: 'gene' column not found in DEA file.")
    exit()
if 'log2FoldChange' not in dea_df.columns:
    print("Error: 'log2FoldChange' column not found in DEA file.")
    exit()
if 'padj' not in dea_df.columns:
    print("Error: 'padj' column not found in DEA file.")
    exit()

# Select relevant columns and rename 'gene' for merging
dea_subset = dea_df[['gene', 'log2FoldChange', 'padj']].copy()
print(f"\nSelected relevant columns from DEA data. Shape: {dea_subset.shape}")

# --- Prepare Peaks data ---
if 'SYMBOL' not in peaks_df.columns:
    print("Error: 'SYMBOL' column not found in peaks file.")
    exit()

# --- Merge the dataframes ---
# Use 'SYMBOL' from peaks_df and 'gene' from dea_subset
print(f"\nMerging peaks data (on 'SYMBOL') with DEA data (on 'gene')...")
merged_df = pd.merge(peaks_df, dea_subset, left_on='SYMBOL', right_on='gene', how='left')
print(f"Merge complete. Shape of merged data: {merged_df.shape}")
print("Columns after merge:", merged_df.columns.tolist())
print("First 5 rows of merged data:\n", merged_df.head().to_string())


# --- Define function for regulation status ---
def get_regulation(row):
    # Check if padj is NaN (meaning the gene wasn't in the DEA results or failed merge)
    if pd.isna(row['padj']):
        return 'Not in DEA'
    # Check for significance
    if row['padj'] < 0.05:
        if row['log2FoldChange'] > 0:
            return 'UP'
        elif row['log2FoldChange'] < 0:
            return 'DOWN'
        else:
            return 'Significant (No Change)' # Edge case, unlikely but possible
    else:
        return 'Not Significant'

# --- Apply the function to create the 'Regulation' column ---
print("\nApplying regulation status logic...")
merged_df['Regulation'] = merged_df.apply(get_regulation, axis=1)
print("Added 'Regulation' column.")
print("Value counts for 'Regulation':\n", merged_df['Regulation'].value_counts())

# --- Final cleanup (optional: remove the extra 'gene' column from DEA) ---
if 'gene' in merged_df.columns:
     merged_df = merged_df.drop(columns=['gene'])
     print("\nRemoved redundant 'gene' column.")

# --- Save the result ---
try:
    merged_df.to_csv(output_file, index=False)
    print(f"\nSuccessfully saved the results to: {output_file}")
    print("Final DataFrame shape:", merged_df.shape)
    print("Final DataFrame columns:", merged_df.columns.tolist())
    print("First 5 rows of final data:\n", merged_df.head().to_string())

except Exception as e:
    print(f"Error writing output file: {e}")