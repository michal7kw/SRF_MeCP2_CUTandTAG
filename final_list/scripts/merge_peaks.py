#!/usr/bin/env python3

import pandas as pd
import pybedtools
import os
import sys

def merge_replicate_peaks(peak_files, output_file, min_overlap=0.5):
    """Merge peaks from replicates with minimum overlap requirement."""
    try:
        if not peak_files:
            raise ValueError("No peak files provided")

        # Read peaks from all files
        peaks = [pybedtools.BedTool(f) for f in peak_files]
        
        # Start with first file
        merged = peaks[0]
        
        # Merge with other replicates
        for peak in peaks[1:]:
            merged = merged.intersect(peak, 
                                   wa=True,
                                   f=min_overlap,
                                   r=True)
        
        # Sort and save
        merged = merged.sort()
        merged.saveas(output_file)
        
        print(f"Successfully merged {len(peak_files)} peak files into {output_file}")
        return True
    except Exception as e:
        print(f"Error merging peaks: {str(e)}")
        return False

if __name__ == "__main__":
    try:
        # Get input/output from snakemake
        peak_files = snakemake.input.peaks
        output_file = snakemake.output.merged_peaks
        min_overlap = snakemake.params.min_overlap
        
        if not merge_replicate_peaks(peak_files, output_file, min_overlap):
            sys.exit(1)
                
    except Exception as e:
        print(f"Script failed: {str(e)}")
        sys.exit(1)
