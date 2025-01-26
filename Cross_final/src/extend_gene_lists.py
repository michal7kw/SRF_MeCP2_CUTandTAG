import pandas as pd
import os

# Define paths
base_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
data_dir = os.path.join(base_dir, "Cross_final/data")
annotation_file = os.path.join(base_dir, "iterative_alternative/results/no_dedup/cpg_enrichment/NSC/broad/cpg_enrichment_2_rep_in_peaks/peaks_annotation_combined/complete_peak_annotation.csv")

# Read the complete peak annotation file
print("Reading complete peak annotation file...")
peak_annotation = pd.read_csv(annotation_file)
print("\nPeak annotation columns:", peak_annotation.columns.tolist())

# Create a mapping dictionary with only the required columns
gene_info = peak_annotation[['SYMBOL', 'seqnames', 'start', 'end', 'distanceToTSS', 'geneId']].copy()
gene_info = gene_info.rename(columns={'SYMBOL': 'Gene'})  

# Process each input file
input_files = ['down.csv', 'no_deg.csv', 'up.csv']

for input_file in input_files:
    print(f"\nProcessing {input_file}...")
    
    # Read input gene list with explicit column name
    input_path = os.path.join(data_dir, input_file)
    genes = pd.read_csv(input_path, names=['Gene'])  
    print(f"Input file columns: {genes.columns.tolist()}")
    print(f"First few rows of {input_file}:")
    print(genes.head())
    
    # Merge with annotation information
    extended_genes = pd.merge(genes, gene_info, on='Gene', how='left')
    print(f"\nFirst few rows of extended {input_file}:")
    print(extended_genes.head())
    
    # Save extended file
    output_file = f"extended_{input_file}"
    output_path = os.path.join(data_dir, output_file)
    extended_genes.to_csv(output_path, index=False)
    print(f"Created {output_file}")

print("\nProcessing completed!")
