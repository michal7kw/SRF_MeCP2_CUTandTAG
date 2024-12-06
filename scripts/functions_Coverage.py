# %%
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import pickle
import os
from pathlib import Path
import mygene


######################## Per CpG Coverage ########################################################################################################################################################################
# %%
wd_dir = '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline'
os.chdir(wd_dir)

# %%
def filter_standard_chromosomes(input_file, output_file):
    """
    Filter BED file to keep only standard chromosomes (chr1-chr22, chrX, chrY, chrM)
    
    Args:
        input_file: Path to input BED file
        output_file: Path to output filtered BED file
    
    Returns:
        tuple: (filtered_count, total_count) or (0, 0) if file doesn't exist
    """
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist")
        return 0, 0
    
    try:
        # Count total peaks
        total = subprocess.run(f"wc -l {input_file}", shell=True, capture_output=True, text=True)
        if total.returncode != 0:
            print(f"Error: Failed to count lines in '{input_file}'")
            return 0, 0
            
        total = int(total.stdout.split()[0])
        
        # Filter and count remaining peaks
        cmd = f"grep -P '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\t' {input_file} > {output_file}"
        subprocess.run(cmd, shell=True)
        
        filtered = subprocess.run(f"wc -l {output_file}", shell=True, capture_output=True, text=True)
        filtered = int(filtered.stdout.split()[0])
        
        print(f"Filtered {input_file}: kept {filtered}/{total} peaks ({filtered/total*100:.1f}%)")
        return filtered, total
        
    except Exception as e:
        print(f"Error processing '{input_file}': {str(e)}")
        return 0, 0

# %%
def calculate_peak_cpg_coverage_exact(peak_file, cpg_file):
    """
    Calculate what percentage of each peak overlaps with CpG islands
    """
    # Use bedtools intersect to find overlaps
    cmd = f"bedtools intersect -a {peak_file} -b {cpg_file} -wao"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    # Parse the output
    overlaps = []
    for line in result.stdout.strip().split('\n'):
        fields = line.split('\t')
        if len(fields) >= 13:  # Ensure we have enough fields
            # peak_length = int(fields[2]) - int(fields[1])
            cpg_length = abs(int(fields[2]) - int(fields[1]))
            overlap_length = int(fields[-1])
            if cpg_length > 0:
                coverage_percent = (overlap_length / cpg_length) * 100
                overlaps.append(coverage_percent)
    
    return overlaps

# %%
import subprocess
import os

def calculate_peak_cpg_coverage(peak_file: str, cpg_file: str, extend: int = 300) -> list:
    """
    Calculate what percentage of each peak overlaps with CpG islands
    CpG islands are extended by extend bp on each side
    """
    try:
        # First extend CpG islands
        cmd_extend = f"bedtools slop -i {cpg_file} -g DATA/genome.size -b {extend} > temp_extended_cpg.bed"
        result = subprocess.run(cmd_extend, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print("Error with bedtools slop:", result.stderr)
            return []
        
        # Use bedtools intersect with extended CpG islands
        cmd = f"bedtools intersect -a {peak_file} -b temp_extended_cpg.bed -wao"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print("Error with bedtools intersect:", result.stderr)
            return []
            
        # Parse the output
        overlaps = []
        for line in result.stdout.strip().split('\n'):
            fields = line.split('\t')
            
            # Find CpG coordinates by looking for 'CpG:' marker
            cpg_start = None
            cpg_end = None
            
            for i in range(len(fields)):
                if fields[i] == 'CpG:':
                    try:
                        cpg_start = int(fields[i-3])
                        cpg_end = int(fields[i-2])
                        overlap_length = int(fields[-1])
                        break
                    except (ValueError, IndexError):
                        continue
            
            # Calculate coverage only if we found valid CpG coordinates
            if cpg_start is not None and cpg_end is not None:
                cpg_length = abs(cpg_end - cpg_start)
                if cpg_length > 0 and overlap_length > 0:
                    coverage_percent = (overlap_length / cpg_length) * 100
                    overlaps.append(coverage_percent)
        
        return overlaps
        
    finally:
        # Clean up temporary file
        if os.path.exists("temp_extended_cpg.bed"):
            os.remove("temp_extended_cpg.bed")


# %%
def plot_coverage_histograms(exo_coverage, endo_coverage, min_coverage=0, max_coverage=100, n_bins=50):
    """
    Create histograms for both types of coverage in separate subplots
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Filter out zeros and 100% coverage values
    exo_nonzero = [x for x in exo_coverage if x > min_coverage and x < max_coverage]
    endo_nonzero = [x for x in endo_coverage if x > min_coverage and x < max_coverage]
    
    # Exogenous plot
    ax1.hist(exo_nonzero, bins=n_bins, color='blue', density=True)
    ax1.set_xlabel('CpG Island Coverage (%)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Exogenous Peaks')
    
    # Add exogenous summary statistics
    exo_mean = sum(x > 0 for x in exo_coverage) / len(exo_coverage) * 100
    ax1.text(0.02, 0.98,
             f'Peaks with CpG: {exo_mean:.1f}%',
             transform=ax1.transAxes,
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Endogenous plot  
    ax2.hist(endo_nonzero, bins=n_bins, color='red', density=True)
    ax2.set_xlabel('CpG Island Coverage (%)')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Endogenous Peaks')
    
    # Add endogenous summary statistics
    endo_mean = sum(x > 80 for x in endo_coverage) / len(endo_coverage) * 100
    ax2.text(0.02, 0.98,
             f'Peaks with CpG: {endo_mean:.1f}%',
             transform=ax2.transAxes,
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    # plt.savefig('cpg_coverage_histogram.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

# %%
def plot_coverage_histograms_by_count(exo_coverage, endo_coverage, min_coverage=0, max_coverage=100, n_bins=50):
    """
    Create histograms for both types of coverage in separate subplots
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Filter out zeros and 100% coverage values
    exo_nonzero = [x for x in exo_coverage if x > min_coverage and x < max_coverage]
    endo_nonzero = [x for x in endo_coverage if x > min_coverage and x < max_coverage]
    
    # Exogenous plot
    ax1.hist(exo_nonzero, bins=n_bins, color='blue')
    ax1.set_xlabel('CpG Island Coverage (%)')
    ax1.set_ylabel('Number of Peaks')
    ax1.set_title('Exogenous Peaks')
    
    # Add exogenous summary statistics
    exo_mean = sum(x > 0 for x in exo_coverage) / len(exo_coverage) * 100
    ax1.text(0.02, 0.98,
             f'Peaks with CpG: {exo_mean:.1f}%',
             transform=ax1.transAxes,
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Endogenous plot  
    ax2.hist(endo_nonzero, bins=50, color='red')
    ax2.set_xlabel('CpG Island Coverage (%)')
    ax2.set_ylabel('Number of Peaks')
    ax2.set_title('Endogenous Peaks')
    
    # Add endogenous summary statistics
    prc_t = 80
    endo_mean = sum(x > prc_t for x in endo_coverage) / len(endo_coverage) * 100
    ax2.text(0.02, 0.98,
             f'Peaks with CpG (>{prc_t}%): {endo_mean:.1f}%',
             transform=ax2.transAxes,
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    # plt.savefig('cpg_coverage_histogram.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


# %%
def plot_coverage_histograms_overlayed(exo_coverage, endo_coverage, min_coverage=0, max_coverage=100, n_bins=50):
    """
    Create an overlayed histogram for both types of coverage
    """
    fig, ax = plt.subplots(figsize=(16, 6))
    
    # Filter out zeros and 100% coverage values
    exo_nonzero = [x for x in exo_coverage if x > min_coverage and x < max_coverage]
    endo_nonzero = [x for x in endo_coverage if x > min_coverage and x < max_coverage]
    
    # Overlayed plot
    ax.hist(exo_nonzero, bins=n_bins, color='blue', alpha=0.5, label='Exogenous Peaks', density=True)
    ax.hist(endo_nonzero, bins=n_bins, color='red', alpha=0.5, label='Endogenous Peaks', density=True)
    ax.set_xlabel('CpG Island Coverage (%)')
    ax.set_ylabel('Frequency')
    ax.set_title('CpG Island Coverage Histogram')
    
    # Add summary statistics
    exo_mean = sum(x > 0 for x in exo_coverage) / len(exo_coverage) * 100
    prc_t = 80
    endo_mean = sum(x > prc_t for x in endo_coverage) / len(endo_coverage) * 100
    ax.text(0.02, 0.98,
            f'Exogenous Peaks with CpG: {exo_mean:.1f}%\nEndogenous Peaks with CpG (>{prc_t}%): {endo_mean:.1f}%',
            transform=ax.transAxes,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.legend()
    plt.tight_layout()
    # plt.savefig('cpg_coverage_histogram.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


# %%
def convert_symbols_to_ensembl(gene_symbols):
    """Convert gene symbols to Ensembl IDs"""
    mg = mygene.MyGeneInfo()
    results = mg.querymany(gene_symbols, scopes='symbol', fields='ensembl.gene', species='mouse')
    
    # Create mapping dictionary
    id_map = {}
    for r in results:
        if 'ensembl' in r:
            if isinstance(r['ensembl'], list):
                # Take the first ensembl ID if there are multiple
                id_map[r['query']] = r['ensembl'][0]['gene']
            else:
                id_map[r['query']] = r['ensembl']['gene']
    return id_map

def get_common_peaks(peak_file, common_genes, genome):
    """
    Filter peaks to keep only those associated with genes that have both 
    Endogenous and Exogenous promoters
    """
    # Create temporary files
    temp_genes = "temp_common_genes.bed"
    temp_out = "temp_filtered_peaks.bed"
    
    # Convert gene symbols to Ensembl IDs
    ensembl_ids = convert_symbols_to_ensembl(common_genes['gene'].tolist())
    print(f"Converted {len(ensembl_ids)} genes to Ensembl IDs")
    
    # Create pattern for grep
    pattern = '|'.join(ensembl_ids.values())
    
    # Extract relevant genes from genome file
    with open(temp_genes, 'w') as f:
        cmd = f"grep -E '{pattern}' {genome}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        f.write(result.stdout)
    
    # Add 2kb upstream and downstream to include promoter regions
    cmd_extend = f"bedtools slop -i {temp_genes} -g DATA/genome.size -b 2000 > temp_extended_genes.bed"
    subprocess.run(cmd_extend, shell=True)
    
    # Intersect peaks with extended gene regions
    cmd2 = f"bedtools intersect -a {peak_file} -b temp_extended_genes.bed -wa | sort | uniq > {temp_out}"
    subprocess.run(cmd2, shell=True)
    
    # Count peaks for logging
    result_total = subprocess.run(f"wc -l {peak_file}", shell=True, capture_output=True, text=True)
    total_peaks = int(result_total.stdout.split()[0])
    
    result_filtered = subprocess.run(f"wc -l {temp_out}", shell=True, capture_output=True, text=True)
    filtered_peaks = int(result_filtered.stdout.split()[0])
    
    print(f"Filtered {peak_file}: kept {filtered_peaks}/{total_peaks} peaks ({filtered_peaks/total_peaks*100:.1f}%)")
    
    # Clean up intermediate files
    subprocess.run("rm temp_common_genes.bed temp_extended_genes.bed", shell=True)
    
    return temp_out


# %%
def get_expression_level(baseMean, q33, q66):
    """
    Assign expression level based on baseMean value and quantile thresholds
    """
    if baseMean <= q33:
        return 'Low'
    elif baseMean <= q66:
        return 'Medium'
    else:
        return 'High'

# %%
def plot_coverage_histograms_expression(high, medium, low, min_coverage=0, max_coverage=100, p_t = 80, n_bins=30):
    """
    Create histograms for both types of coverage in separate subplots
    """
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 6))
    
    # Filter out zeros and 100% coverage values
    nonzero_high = [x for x in high if x > min_coverage and x < max_coverage]
    nonzero_medium = [x for x in medium if x > min_coverage and x < max_coverage]
    nonzero_low = [x for x in low if x > min_coverage and x < max_coverage]
    
    # High expression
    ax1.hist(nonzero_high, bins=n_bins, color='blue', density=True)
    ax1.set_xlabel('CpG Island Coverage (%)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('High Expression')
    
    # Add exogenous summary statistics
    exo_mean = sum(x > p_t for x in nonzero_high) / len(nonzero_high) * 100
    ax1.text(0.02, 0.98,
             f'Peaks with CpG (>{p_t}%): {exo_mean:.1f}% \ntotal peaks: {len(nonzero_high)}',
             transform=ax1.transAxes,
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Medium expression
    ax2.hist(nonzero_medium, bins=n_bins, color='red', density=True)
    ax2.set_xlabel('CpG Island Coverage (%)')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Medium Expression')
    
    # Add endogenous summary statistics
    endo_mean = sum(x > p_t for x in nonzero_medium) / len(nonzero_medium) * 100
    ax2.text(0.02, 0.98,
             f'Peaks with CpG (>{p_t}%): {endo_mean:.1f}% \ntotal peaks: {len(nonzero_medium)}',
             transform=ax2.transAxes,
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Low expression
    ax3.hist(nonzero_low, bins=n_bins, color='green', density=True)
    ax3.set_xlabel('CpG Island Coverage (%)')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Low Expression')
    
    # Add summary statistics
    endo_mean = sum(x > p_t for x in nonzero_low) / len(nonzero_low) * 100
    ax3.text(0.02, 0.98,
             f'Peaks with CpG (>{p_t}%): {endo_mean:.1f}% \ntotal peaks: {len(nonzero_low)}',
             transform=ax3.transAxes,
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    
    plt.tight_layout()
    # plt.savefig('cpg_coverage_histogram.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


######################## Per Peak Coverage ########################################################################################################################################################################
def calculate_peak_cpg_coverage_per_peak(peak_file: str, cpg_file: str, extend: int = 300) -> dict:
    """
    Calculate what percentage of each peak overlaps with CpG islands
    CpG islands are extended by extend bp on each side
    """
    try:
        # First extend CpG islands
        cmd_extend = f"bedtools slop -i {cpg_file} -g DATA/genome.size -b {extend} > temp_extended_cpg.bed"
        result = subprocess.run(cmd_extend, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print("Error with bedtools slop:", result.stderr)
            subprocess.run(f"cp {cpg_file} temp_extended_cpg.bed", shell=True, check=True)
        
        # Use bedtools intersect with extended CpG islands
        cmd = f"bedtools intersect -a {peak_file} -b temp_extended_cpg.bed -wao"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        
        # Process overlaps and calculate coverage
        peak_coverages = {}
        current_peak = {'id': None, 'overlaps': []}
        
        for line in result.stdout.strip().split('\n'):
            fields = line.split('\t')
            if len(fields) < 13:  # Skip malformed lines
                continue
                
            peak_id = f"{fields[0]}:{fields[1]}-{fields[2]}"
            
            # If we encounter a new peak
            if peak_id != current_peak['id']:
                # Process previous peak if it exists
                if current_peak['id'] is not None and current_peak['overlaps']:
                    max_coverage = max(current_peak['overlaps'])
                    peak_coverages[current_peak['id']] = {
                        'chrom': current_peak['id'].split(':')[0],
                        'start': int(current_peak['id'].split(':')[1].split('-')[0]),
                        'end': int(current_peak['id'].split(':')[1].split('-')[1]),
                        'cpg_overlaps': current_peak['overlaps'],
                        'max_coverage': max_coverage,
                        'total_overlaps': len(current_peak['overlaps'])
                    }
                
                # Start new peak
                current_peak = {'id': peak_id, 'overlaps': []}
            
            # Find CpG information - look for it in the line
            cpg_start = None
            cpg_end = None
            overlap_length = 0
            
            # Try to find CpG information by looking for 'CpG:' marker
            for i in range(len(fields)):
                if fields[i] == 'CpG:':
                    try:
                        cpg_start = int(fields[i-3])  # Three fields before 'CpG:'
                        cpg_end = int(fields[i-2])    # Two fields before 'CpG:'
                        overlap_length = int(fields[-1])  # Last field is overlap length
                        break
                    except (ValueError, IndexError):
                        continue
            
            # Only calculate coverage if we found valid CpG information
            if cpg_start is not None and cpg_end is not None and overlap_length > 0:
                cpg_length = abs(cpg_end - cpg_start)
                if cpg_length > 0:  # Avoid division by zero
                    coverage = (overlap_length / cpg_length) * 100
                    current_peak['overlaps'].append(coverage)
        
        # Process the last peak
        if current_peak['id'] is not None and current_peak['overlaps']:
            max_coverage = max(current_peak['overlaps'])
            peak_coverages[current_peak['id']] = {
                'chrom': current_peak['id'].split(':')[0],
                'start': int(current_peak['id'].split(':')[1].split('-')[0]),
                'end': int(current_peak['id'].split(':')[1].split('-')[1]),
                'cpg_overlaps': current_peak['overlaps'],
                'max_coverage': max_coverage,
                'total_overlaps': len(current_peak['overlaps'])
            }
        
        # Add peaks with no overlaps (zero coverage)
        with open(peak_file) as f:
            for line in f:
                fields = line.strip().split('\t')
                peak_id = f"{fields[0]}:{fields[1]}-{fields[2]}"
                if peak_id not in peak_coverages:
                    peak_coverages[peak_id] = {
                        'chrom': fields[0],
                        'start': int(fields[1]),
                        'end': int(fields[2]),
                        'cpg_overlaps': [],
                        'max_coverage': 0.0,
                        'total_overlaps': 0
                    }
        
        return peak_coverages
    
    finally:
        # Clean up temporary file
        if os.path.exists("temp_extended_cpg.bed"):
            os.remove("temp_extended_cpg.bed")





def analyze_coverage_stats_per_peak(peak_coverages: dict, threshold: float = 0.0) -> dict:
    """
    Analyze coverage statistics for all peaks
    
    Args:
        peak_coverages: Dictionary of peak coverages
        threshold: Minimum coverage threshold to consider
        
    Returns:
        Dictionary containing coverage statistics
    """
    total_peaks = len(peak_coverages)
    peaks_with_overlap = sum(1 for p in peak_coverages.values() if p['cpg_overlaps'])
    peaks_above_threshold = sum(1 for p in peak_coverages.values() if p['max_coverage'] >= threshold)
    
    # Get list of maximum coverages for peaks with any overlap
    max_coverages = [p['max_coverage'] for p in peak_coverages.values() if p['cpg_overlaps']]
    
    stats = {
        'total_peaks': total_peaks,
        'peaks_with_overlap': peaks_with_overlap,
        'peaks_above_threshold': peaks_above_threshold,
        'percent_with_overlap': (peaks_with_overlap / total_peaks * 100) if total_peaks > 0 else 0,
        'percent_above_threshold': (peaks_above_threshold / total_peaks * 100) if total_peaks > 0 else 0,
        'mean_max_coverage': sum(max_coverages) / len(max_coverages) if max_coverages else 0,
        'max_coverage': max(max_coverages) if max_coverages else 0,
        'min_coverage': min(max_coverages) if max_coverages else 0,
        'peaks_by_overlap_count': {}
    }
    
    # Count peaks by number of overlaps
    for peak in peak_coverages.values():
        overlap_count = peak['total_overlaps']
        stats['peaks_by_overlap_count'][overlap_count] = stats['peaks_by_overlap_count'].get(overlap_count, 0) + 1
    
    return stats

def plot_coverage_distribution_per_peak(peak_coverages: dict, title: str = "CpG Coverage Distribution"):
    """
    Plot the distribution of maximum CpG coverage across peaks
    
    Args:
        peak_coverages: Dictionary of peak coverages
        title: Plot title
    """
    import matplotlib.pyplot as plt
    
    # Get maximum coverage values for peaks with any overlap
    max_coverages = [p['max_coverage'] for p in peak_coverages.values() if p['cpg_overlaps']]
    
    plt.figure(figsize=(10, 6))
    plt.hist(max_coverages, bins=50, alpha=0.75, edgecolor='black')
    plt.title(title)
    plt.xlabel("Maximum CpG Coverage (%)")
    plt.ylabel("Number of Peaks")
    
    # Add statistics text
    stats_text = (f"Total peaks: {len(peak_coverages)}\n"
                 f"Peaks with overlap: {len(max_coverages)}\n"
                 f"Mean coverage: {sum(max_coverages)/len(max_coverages):.1f}%")
    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes,
             verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8))
    
    plt.grid(True, alpha=0.3)
    plt.show()