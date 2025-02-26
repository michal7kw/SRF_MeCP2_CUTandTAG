#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os

# Define paths
neu_path = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA/DEA_NEU.csv"
nsc_path = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA/DEA_NSC.csv"
output_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/paper_finals/outputs"
output_file = os.path.join(output_dir, "DEGs_lists.xlsx")

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Parameters
min_reads = 50
fc_threshold = 0.5  # log2FoldChange threshold

def process_deg_file(file_path, sheet_name):
    """
    Process a DEG file and return DataFrames for DEGs and non-DEGs
    """
    print(f"Processing {file_path}...")
    
    # Read the CSV file
    df = pd.read_csv(file_path)
    
    # Convert scientific notation to float if needed
    if 'log2FoldChange' in df.columns:
        df['log2FoldChange'] = df['log2FoldChange'].astype(float)
    
    # Filter for minimum reads
    df_filtered = df[df['baseMean'] >= min_reads].copy()
    
    # Remove RIK and GM genes
    df_valid = df_filtered[~df_filtered['gene'].str.contains('Rik|Gm\d+', regex=True)].copy()
    
    # Separate DEGs and non-DEGs
    degs = df_valid[abs(df_valid['log2FoldChange']) > fc_threshold].copy()
    non_degs = df_valid[abs(df_valid['log2FoldChange']) <= fc_threshold].copy()
    
    # Add a column to indicate if it's a DEG
    degs['is_DEG'] = 'Yes'
    non_degs['is_DEG'] = 'No'
    
    # Add a column to indicate up or down regulation
    degs['regulation'] = np.where(degs['log2FoldChange'] > 0, 'Up', 'Down')
    non_degs['regulation'] = 'None'
    
    # Add a column to indicate the source file
    degs['source'] = sheet_name
    non_degs['source'] = sheet_name
    
    # Sort DEGs by absolute fold change (descending)
    degs = degs.assign(abs_log2FC=abs(degs['log2FoldChange']))
    degs = degs.sort_values('abs_log2FC', ascending=False)
    degs = degs.drop('abs_log2FC', axis=1)
    
    # Sort non-DEGs by gene name
    non_degs = non_degs.sort_values('gene')
    
    print(f"Found {len(degs)} DEGs and {len(non_degs)} non-DEGs in {sheet_name}")
    
    return degs, non_degs

def main():
    # Process NEU data
    neu_degs, neu_non_degs = process_deg_file(neu_path, "NEU")
    
    # Process NSC data
    nsc_degs, nsc_non_degs = process_deg_file(nsc_path, "NSC")
    
    # Create a Pandas Excel writer using XlsxWriter as the engine
    print(f"Writing results to {output_file}...")
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        # Write each DataFrame to a different worksheet
        neu_degs.to_excel(writer, sheet_name='NEU_DEGs', index=False)
        neu_non_degs.to_excel(writer, sheet_name='NEU_non_DEGs', index=False)
        nsc_degs.to_excel(writer, sheet_name='NSC_DEGs', index=False)
        nsc_non_degs.to_excel(writer, sheet_name='NSC_non_DEGs', index=False)
        
        # Create a summary sheet
        summary = pd.DataFrame({
            'Category': ['NEU DEGs', 'NEU non-DEGs', 'NSC DEGs', 'NSC non-DEGs'],
            'Count': [len(neu_degs), len(neu_non_degs), len(nsc_degs), len(nsc_non_degs)]
        })
        summary.to_excel(writer, sheet_name='Summary', index=False)
        
        # Get the xlsxwriter workbook and worksheet objects
        workbook = writer.book
        
        # Add some cell formats
        header_format = workbook.add_format({
            'bold': True,
            'text_wrap': True,
            'valign': 'top',
            'border': 1
        })
        
        # Apply the header format to each sheet
        for sheet_name in writer.sheets:
            worksheet = writer.sheets[sheet_name]
            # Get the dataframe corresponding to this sheet
            if sheet_name == 'NEU_DEGs':
                df = neu_degs
            elif sheet_name == 'NEU_non_DEGs':
                df = neu_non_degs
            elif sheet_name == 'NSC_DEGs':
                df = nsc_degs
            elif sheet_name == 'NSC_non_DEGs':
                df = nsc_non_degs
            elif sheet_name == 'Summary':
                df = summary
            
            # Write the column headers with the format
            for col_num, value in enumerate(df.columns.values):
                worksheet.write(0, col_num, value, header_format)
            
            # Adjust columns' width (alternative to autofit which may not be available in all versions)
            for col_num, value in enumerate(df.columns.values):
                # Set width based on column name length (with some padding)
                col_width = max(len(str(value)) + 2, 10)
                worksheet.set_column(col_num, col_num, col_width)
    
    print(f"Excel file created successfully: {output_file}")

if __name__ == "__main__":
    main()