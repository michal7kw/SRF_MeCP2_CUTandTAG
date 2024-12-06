#!/usr/bin/env python3

import pandas as pd
import glob
import os

def convert_peaks_to_csv():
    # Define input and output directories
    input_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/iterative_processing/results/peaks/seacr"
    output_dir = os.path.join(input_dir, "csv")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Define columns for SEACR output
    columns = [
        'chromosome',
        'start',
        'end',
        'signal_value',
        'area',
        'peak_score'
    ]
    
    # Process both control and no-control peak files
    for pattern in ["*.peaks.stringent.bed", "*_no_control.peaks.stringent.bed"]:
        for bed_file in glob.glob(os.path.join(input_dir, pattern)):
            # Read the BED file
            df = pd.read_csv(bed_file, sep='\t', header=None, names=columns)
            
            # Create output filename
            base_name = os.path.basename(bed_file)
            sample_name = base_name.replace('.peaks.stringent.bed', '')
            csv_file = os.path.join(output_dir, f"{sample_name}.csv")
            
            # Add additional useful columns
            df['peak_length'] = df['end'] - df['start']
            df['peak_center'] = df['start'] + (df['peak_length'] // 2)
            
            # Sort by peak score (descending)
            df = df.sort_values('peak_score', ascending=False)
            
            # Add rank column
            df['rank'] = range(1, len(df) + 1)
            
            # Save to CSV
            df.to_csv(csv_file, index=False)
            print(f"Converted {bed_file} to {csv_file}")

if __name__ == "__main__":
    convert_peaks_to_csv()