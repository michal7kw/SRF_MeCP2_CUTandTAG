import os
import glob
import pandas as pd
from config import PATHS, logger

def merge_results():
    """Merge results from individual chromosome analyses"""
    results = {}
    
    # For each cell type
    for cell_type in ['NEU', 'NSC']:
        cell_results = []
        
        # Find all chromosome results for this cell type
        pattern = f"{PATHS['output_dir']}/{cell_type}_chr*_methylation_analysis.csv"
        chr_files = glob.glob(pattern)
        
        # Merge chromosome results
        for chr_file in chr_files:
            df = pd.read_csv(chr_file)
            cell_results.append(df)
        
        if cell_results:
            # Combine all chromosomes
            results[cell_type] = pd.concat(cell_results, ignore_index=True)
            
            # Save combined results
            output_file = f"{PATHS['output_dir']}/{cell_type}_methylation_analysis.csv"
            results[cell_type].to_csv(output_file, index=False)
            logger.info(f"Saved merged results for {cell_type} to {output_file}")
            
            # Clean up individual chromosome files
            for f in chr_files:
                os.remove(f)
    
    return results

if __name__ == "__main__":
    merge_results() 