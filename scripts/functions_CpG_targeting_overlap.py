
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import os
from pathlib import Path
from matplotlib_venn import venn2 
from upsetplot import from_contents, UpSet
import numpy as np

# wd_dir = '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline'
# os.chdir(wd_dir)


def get_peaks_with_cpg(peak_file, cpg_file, extend=300, coverage_threshold=20, genome_size_file="DATA/genome.size"):
    """
    Identifies peaks that overlap with CpG islands above a coverage threshold.
    
    Args:
        peak_file: BED file containing peak regions
        cpg_file: BED file containing CpG islands
        extend: Number of base pairs to extend CpG islands (default: 100bp)
        coverage_threshold: Minimum percentage overlap required (default: 20%)
    
    Returns:
        set: Peak coordinates that meet the CpG overlap criteria
    """
    print(f"\nProcessing {peak_file}")
    
    # Validate input files
    if not all(os.path.exists(f) for f in [peak_file, cpg_file]):
        print("Error: Input file(s) missing!")
        return set()
    
    # Filter for standard chromosomes (chr1-22, X, Y)
    standard_chroms = [f"chr{i}" for i in range(1, 23)] + ['chrX', 'chrY']
    
    # Create temporary filtered files
    filtered_peaks = "temp_filtered_peaks.bed"
    filtered_cpg = "temp_filtered_cpg.bed"
    
    # Filter peaks and count
    with open(peak_file) as fin, open(filtered_peaks, 'w') as fout:
        peaks = [line for line in fin if line.split('\t')[0] in standard_chroms]
        fout.writelines(peaks)
        print(f"Peaks after chromosome filtering: {len(peaks)}")
    
    # Filter and extend CpG islands
    with open(cpg_file) as fin, open(filtered_cpg, 'w') as fout:
        fout.writelines(line for line in fin if line.split('\t')[0] in standard_chroms)
    
    # Extend CpG regions and find overlaps using bedtools
    subprocess.run(f"bedtools slop -i {filtered_cpg} -g {genome_size_file} -b {extend} > temp_extended_cpg.bed", shell=True)
    result = subprocess.run(f"bedtools intersect -a {filtered_peaks} -b temp_extended_cpg.bed -wao", 
                          shell=True, capture_output=True, text=True)
    
    # Process overlaps and calculate coverage
    peaks_with_cpg = set()
    coverage_stats = []
    current_peak = {'id': None, 'overlaps': []}
    qualified_peaks = 0
    total_peaks_with_overlap = 0
    
    for line in result.stdout.strip().split('\n'):
            fields = line.split('\t')
            peak_id = f"{fields[0]}:{fields[1]}-{fields[2]}"
            
            # If we encounter a new peak
            if peak_id != current_peak['id']:
                # Process previous peak if it exists
                if current_peak['overlaps']:
                    max_coverage = max(current_peak['overlaps'])
                    total_peaks_with_overlap += 1
                    if max_coverage >= coverage_threshold:
                        coverage_stats.append(max_coverage)
                        peaks_with_cpg.add(current_peak['id'])
                        qualified_peaks += 1
                
                current_peak = {'id': peak_id, 'overlaps': []}
                
            # Calculate overlap coverage if exists
            if fields[4] != "." and int(fields[-1]) > 0:
                # Calculate CpG length instead of peak length
                cpg_length = int(fields[6]) - int(fields[5])
                if cpg_length > 0:  # Avoid division by zero
                    coverage = (int(fields[-1]) / cpg_length) * 100
                    current_peak['overlaps'].append(coverage)
    
    # Process the last peak
    if current_peak['overlaps']:
        max_coverage = max(current_peak['overlaps'])
        total_peaks_with_overlap += 1
        if max_coverage >= coverage_threshold:
            coverage_stats.append(max_coverage)
            peaks_with_cpg.add(current_peak['id'])
            qualified_peaks += 1
    
    # Print coverage statistics
    if coverage_stats:
        print(f"\nPeaks with any CpG overlap: {total_peaks_with_overlap}")
        print(f"Peaks meeting {coverage_threshold}% coverage threshold: {qualified_peaks}")
        print(f"Coverage range: {min(coverage_stats):.2f}% - {max(coverage_stats):.2f}%")
        print(f"Mean coverage of qualified peaks: {sum(coverage_stats)/len(coverage_stats):.2f}%")
    
    # Cleanup temporary files
    for f in [filtered_peaks, filtered_cpg, "temp_extended_cpg.bed"]:
        if os.path.exists(f):
            os.remove(f)
    
    return peaks_with_cpg

def analyze_cpg_overlap(exo_peaks, endo_peaks, cpg_file, output_dir, extend=300, coverage_threshold=20, genome_size_file="DATA/genome.size"):
    """
    Analyzes overlap between exogenous and endogenous peaks at CpG islands.
    
    Args:
        exo_peaks: BED file with exogenous peaks
        endo_peaks: BED file with endogenous peaks
        cpg_file: BED file with CpG islands
        output_dir: Output directory for results
    """
    # Get peaks overlapping CpG islands
    exo_cpg_peaks = get_peaks_with_cpg(exo_peaks, cpg_file, extend, coverage_threshold, genome_size_file)
    endo_cpg_peaks = get_peaks_with_cpg(endo_peaks, cpg_file, extend, coverage_threshold, genome_size_file)
    
    # Find overlapping peaks between Exo and Endo
    common_peaks = find_overlapping_peaks(exo_cpg_peaks, endo_cpg_peaks)
    
    # Get specific peaks
    exo_specific = exo_cpg_peaks - common_peaks['exo']
    endo_specific = endo_cpg_peaks - common_peaks['endo']
    
    # Create visualization
    create_venn_diagram(
        exo_specific=len(exo_specific),
        endo_specific=len(endo_specific),
        common=len(common_peaks['exo']),
        total_exo=len(exo_cpg_peaks),
        total_endo=len(endo_cpg_peaks),
        output_dir=output_dir
    )
    
    return {
        'common_exo': common_peaks['exo'],
        'common_endo': common_peaks['endo'],
        'exo_specific': exo_specific,
        'endo_specific': endo_specific
    }

def find_overlapping_peaks(exo_peaks, endo_peaks):
    """Helper function to find overlapping peaks using bedtools"""
    # Write peaks to temporary BED files
    for peaks, filename in [(exo_peaks, 'temp_exo.bed'), (endo_peaks, 'temp_endo.bed')]:
        with open(filename, 'w') as f:
            for peak in peaks:
                chrom, pos = peak.split(':')
                start, end = pos.split('-')
                f.write(f"{chrom}\t{start}\t{end}\n")
    
    # Find overlaps using bedtools
    subprocess.run("bedtools intersect -a temp_exo.bed -b temp_endo.bed -wa -wb > temp_overlaps.bed", shell=True)
    
    # Parse overlaps
    common_exo = set()
    common_endo = set()
    with open('temp_overlaps.bed') as f:
        for line in f:
            fields = line.strip().split('\t')
            common_exo.add(f"{fields[0]}:{fields[1]}-{fields[2]}")
            common_endo.add(f"{fields[3]}:{fields[4]}-{fields[5]}")
    
    # Cleanup
    for f in ['temp_exo.bed', 'temp_endo.bed', 'temp_overlaps.bed']:
        if os.path.exists(f):
            os.remove(f)
    
    return {'exo': common_exo, 'endo': common_endo}

def create_venn_diagram(exo_specific, endo_specific, common, total_exo, total_endo, output_dir):
    """Creates a Venn diagram showing peak overlaps"""
    plt.figure(figsize=(10, 10))
    venn2(subsets=(exo_specific, endo_specific, common),
          set_labels=('Exogenous', 'Endogenous'),
          set_colors=('lightblue', 'lightgreen'))
    
    plt.title(f'CpG Peak Overlap - {os.path.basename(output_dir)}')
    plt.figtext(0.1, 0.02,
                f'Exo shared: {common/total_exo*100:.1f}%\n'
                f'Endo shared: {common/total_endo*100:.1f}%',
                fontsize=10)
    
    plt.savefig(os.path.join(output_dir, "cpg_overlap_venn.pdf"),
                bbox_inches='tight', dpi=300)
    plt.close()

def analyze_coverage_distribution(peak_file, cpg_file, extend=300, genome_size_file="DATA/genome.size"):
    """
    Analyzes and returns the coverage distribution of peaks overlapping with CpG islands.
    
    Args:
        peak_file: BED file containing peak regions
        cpg_file: BED file containing CpG islands
        extend: Number of base pairs to extend CpG islands
    
    Returns:
        list: Coverage percentages for all peaks with any CpG overlap
    """
    # Filter for standard chromosomes
    standard_chroms = [f"chr{i}" for i in range(1, 23)] + ['chrX', 'chrY']
    
    # Create temporary filtered files
    filtered_peaks = "temp_filtered_peaks.bed"
    filtered_cpg = "temp_filtered_cpg.bed"
    
    # Filter files
    with open(peak_file) as fin, open(filtered_peaks, 'w') as fout:
        fout.writelines(line for line in fin if line.split('\t')[0] in standard_chroms)
    
    with open(cpg_file) as fin, open(filtered_cpg, 'w') as fout:
        fout.writelines(line for line in fin if line.split('\t')[0] in standard_chroms)
    
    # Extend CpG regions and find overlaps
    subprocess.run(f"bedtools slop -i {filtered_cpg} -g {genome_size_file} -b {extend} > temp_extended_cpg.bed", shell=True)
    result = subprocess.run(f"bedtools intersect -a {filtered_peaks} -b temp_extended_cpg.bed -wao", 
                          shell=True, capture_output=True, text=True)
    
    # Process overlaps
    coverage_values = []
    current_peak = {'id': None, 'overlaps': []}
    
    for line in result.stdout.strip().split('\n'):
        fields = line.split('\t')
        peak_id = f"{fields[0]}:{fields[1]}-{fields[2]}"
        
        if peak_id != current_peak['id']:
            if current_peak['overlaps']:
                coverage_values.append(max(current_peak['overlaps']))
            current_peak = {'id': peak_id, 'overlaps': []}
        
        if fields[4] != "." and int(fields[-1]) > 0:
            cpg_length = int(fields[6]) - int(fields[5])
            if cpg_length > 0:
                coverage = (int(fields[-1]) / cpg_length) * 100
                current_peak['overlaps'].append(coverage)
    
    # Add last peak
    if current_peak['overlaps']:
        coverage_values.append(max(current_peak['overlaps']))
    
    # Cleanup
    for f in [filtered_peaks, filtered_cpg, "temp_extended_cpg.bed"]:
        if os.path.exists(f):
            os.remove(f)
    
    return coverage_values

def plot_coverage_distributions(genome_size_file="DATA/genome.size", cpg_file="DATA/cpg_islands.bed"):
    """
    Creates overlaid histogram plots comparing Endo vs Exo coverage distributions.
    """
    conditions = {
        'NSC': {
            'Exo': "results/consensus_peaks/NSC_Exo_consensus.bed",
            'Endo': "results/consensus_peaks/NSC_Endo_consensus.bed"
        },
        'Neuron': {
            'Exo': "results/consensus_peaks/Neuron_Exo_consensus.bed",
            'Endo': "results/consensus_peaks/Neuron_Endo_consensus.bed"
        }
    }
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    axes = {'NSC': ax1, 'Neuron': ax2}
    colors = {'Exo': 'lightblue', 'Endo': 'lightgreen'}
    alpha = 0.5
    
    # Plot overlaid histograms
    for cell_type, peaks in conditions.items():
        ax = axes[cell_type]
        
        for condition, peak_file in peaks.items():
            coverage_values = analyze_coverage_distribution(peak_file, cpg_file, genome_size_file=genome_size_file)
            color = colors[condition]
            
            # Calculate normalized histogram values
            hist, bins = np.histogram(coverage_values, bins=50, density=True)
            bin_centers = (bins[:-1] + bins[1:]) / 2
            
            ax.hist(coverage_values, bins=50, color=color, edgecolor='black', 
                   alpha=alpha, label=f'{condition} (n={len(coverage_values)})',
                   density=True)
            
            # Calculate statistics
            mean_cov = np.mean(coverage_values)
            median_cov = np.median(coverage_values)
            stats_text = f'{condition}:\nMean: {mean_cov:.1f}%\nMedian: {median_cov:.1f}%'
            
            # Add statistics text (Exo on top left, Endo on top right)
            if condition == 'Exo':
                ax.text(0.02, 0.95, stats_text, transform=ax.transAxes, 
                       bbox=dict(facecolor=color, alpha=0.8))
            else:
                ax.text(0.75, 0.95, stats_text, transform=ax.transAxes,
                       bbox=dict(facecolor=color, alpha=0.8))
        
        ax.set_title(f'{cell_type} Coverage Distribution')
        ax.set_xlabel('Coverage (%)')
        ax.set_ylabel('Frequency (density)')
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    plt.tight_layout()
    plt.show()