#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool
from pybedtools import BedTool
import argparse
from tqdm import tqdm
import time
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Color palette for different sample types
COLORS = {
    "NSC_ENDO": "#83CBEB",
    "NSC_EXO": "#0070C0",
    "NEU_ENDO": "#FF9999",
    "NEU_EXO": "#FF3300"
}

def extract_signal(args):
    """Extract signal from a bigWig file for a given region"""
    bw_file, chrom, start, end, bin_size = args
    
    try:
        bw = pyBigWig.open(bw_file)
        
        # Calculate bins
        bins = np.linspace(start, end, int((end - start) / bin_size) + 1, dtype=int)
        signal = np.zeros(len(bins) - 1)
        
        # Extract signal for each bin
        for i in range(len(bins) - 1):
            bin_start = bins[i]
            bin_end = bins[i + 1]
            
            try:
                # Get mean signal for the bin
                signal[i] = bw.stats(chrom, bin_start, bin_end, type="mean")[0] or 0
            except:
                signal[i] = 0
                
        bw.close()
        return signal
    
    except Exception as e:
        logger.error(f"Error processing {bw_file}: {str(e)}")
        # Return zeros if there's an error
        return np.zeros(int((end - start) / bin_size))

def compute_matrix(bw_files, regions, upstream=3000, downstream=3000, bin_size=50, n_processes=8):
    """Compute signal matrix for multiple bigWig files and regions"""
    start_time = time.time()
    logger.info(f"Computing matrix for {len(bw_files)} bigWig files and {len(regions)} regions")
    
    # Prepare arguments for parallel processing
    args_list = []
    for region in regions:
        chrom = region.chrom
        # Calculate center position
        center = int((region.start + region.end) / 2)
        start = center - upstream
        end = center + downstream
        
        for bw_file in bw_files:
            args_list.append((bw_file, chrom, start, end, bin_size))
    
    # Process in parallel
    with Pool(processes=n_processes) as pool:
        results = list(tqdm(pool.imap(extract_signal, args_list), total=len(args_list)))
    
    # Reshape results into a matrix [n_regions × n_files, n_bins]
    n_regions = len(regions)
    n_files = len(bw_files)
    n_bins = int((upstream + downstream) / bin_size)
    
    # Reshape and reorder to get [n_files, n_regions, n_bins]
    matrix = np.zeros((n_files, n_regions, n_bins))
    
    for i in range(n_regions):
        for j in range(n_files):
            idx = i * n_files + j
            if idx < len(results):
                signal = results[idx]
                if len(signal) == n_bins:
                    matrix[j, i, :] = signal
    
    # Average over regions to get [n_files, n_bins]
    avg_matrix = np.mean(matrix, axis=1)
    
    logger.info(f"Matrix computation completed in {time.time() - start_time:.2f} seconds")
    
    return avg_matrix

def plot_profiles(matrices, labels, colors, output_file, title="MeCP2 at CpG Islands",
                 xlabel="Distance from CpG Island Center (bp)", ylabel="Signal"):
    """Create a profile plot from computed matrices"""
    plt.figure(figsize=(10, 6))
    
    x = np.linspace(-3000, 3000, matrices[0].shape[0])
    
    for i, (matrix, label) in enumerate(zip(matrices, labels)):
        plt.plot(x, matrix, label=label, color=colors[i], linewidth=2)
    
    plt.title(title, fontsize=14)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.legend(loc='upper right')
    plt.grid(alpha=0.3)
    
    # Add vertical line at center
    plt.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    logger.info(f"Plot saved to {output_file}")
    
    return plt.gcf()

def plot_side_by_side(matrices, labels, colors, output_file, title_prefix="MeCP2",
                    xlabel="Distance from CpG Island Center (bp)", ylabel="Signal"):
    """Create side-by-side profile plots similar to the provided example image"""
    fig, axs = plt.subplots(1, len(matrices), figsize=(15, 5), sharey=True)
    
    x = np.linspace(-3000, 3000, matrices[0].shape[0])
    
    for i, (matrix, label, color) in enumerate(zip(matrices, labels, colors)):
        axs[i].plot(x, matrix, color=color, linewidth=2)
        axs[i].set_title(f"{title_prefix}{label}", fontsize=14)
        axs[i].grid(alpha=0.3)
        axs[i].axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        
        # Add x-axis label
        if i == 0:
            axs[i].set_ylabel(ylabel, fontsize=12)
        
        # Format x-axis
        axs[i].set_xlabel(xlabel, fontsize=12)
        axs[i].set_xticks([-3000, 0, 3000])
        axs[i].set_xticklabels(["-3.0Kb", "CpG Island Center", "3.0Kb"])
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    logger.info(f"Side-by-side plot saved to {output_file}")
    
    return fig

def main():
    # Define sample groups
    sample_groups = {
        "neuron_endo": ["NeuM2.bw", "NeuM3.bw"],
        "neuron_exo": ["NeuV1.bw", "NeuV2.bw", "NeuV3.bw"],
        "nsc_endo": ["NSCM1.bw", "NSCM2.bw", "NSCM3.bw"],
        "nsc_exo": ["NSCv1.bw", "NSCv2.bw", "NSCv3.bw"]
    }
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="MeCP2 CpG island profiling")
    parser.add_argument("--bigwig-dir", required=True, help="Directory containing bigWig files")
    parser.add_argument("--cpg-file", required=True, help="BED file with CpG islands")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--upstream", type=int, default=3000, help="Upstream distance (bp)")
    parser.add_argument("--downstream", type=int, default=3000, help="Downstream distance (bp)")
    parser.add_argument("--bin-size", type=int, default=50, help="Bin size (bp)")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load CpG islands
    logger.info(f"Loading CpG islands from {args.cpg_file}")
    cpg_islands = BedTool(args.cpg_file)
    
    # Process each sample group
    matrices = {}
    for group_name, sample_files in sample_groups.items():
        # Get full paths to bigWig files
        bw_files = [os.path.join(args.bigwig_dir, f) for f in sample_files]
        
        # Compute matrix
        logger.info(f"Computing matrix for {group_name}")
        matrix = compute_matrix(
            bw_files, 
            cpg_islands, 
            upstream=args.upstream, 
            downstream=args.downstream,
            bin_size=args.bin_size,
            n_processes=args.threads
        )
        
        matrices[group_name] = np.mean(matrix, axis=0)
    
    # Create individual plots
    for cell_type in ["neuron", "nsc"]:
        endo_matrix = matrices[f"{cell_type}_endo"]
        exo_matrix = matrices[f"{cell_type}_exo"]
        
        output_file = os.path.join(args.output_dir, f"{cell_type}_endo_exo_cpg_profile.png")
        plot_profiles(
            [endo_matrix, exo_matrix],
            ["Endogenous", "Exogenous"],
            [COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_ENDO"], 
             COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_EXO"]],
            output_file,
            title=f"{cell_type.title()} MeCP2 at CpG Islands"
        )

        # NEW: Create side-by-side plots
        side_by_side_output = os.path.join(args.output_dir, f"{cell_type}_endo_exo_side_by_side.png")
        plot_side_by_side(
            [endo_matrix, exo_matrix],
            ["M", "V"],  # Using M for endogenous and V for exogenous as in example
            [COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_ENDO"], 
             COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_EXO"]],
            side_by_side_output,
            title_prefix=f"{cell_type.title()}"
        )
    
    # Create combined plot
    all_matrices = [
        matrices["neuron_endo"],
        matrices["neuron_exo"],
        matrices["nsc_endo"],
        matrices["nsc_exo"]
    ]
    
    all_labels = ["Neuron Endo", "Neuron Exo", "NSC Endo", "NSC Exo"]
    all_colors = [COLORS["NEU_ENDO"], COLORS["NEU_EXO"], COLORS["NSC_ENDO"], COLORS["NSC_EXO"]]
    
    output_file = os.path.join(args.output_dir, "all_samples_cpg_profile.png")
    plot_profiles(
        all_matrices,
        all_labels,
        all_colors,
        output_file,
        title="MeCP2 at CpG Islands"
    )
    
    # NEW: Create side-by-side plot for all samples
    all_side_by_side_output = os.path.join(args.output_dir, "all_samples_side_by_side.png")
    plot_side_by_side(
        all_matrices,
        ["NEU_M", "NEU_V", "NSC_M", "NSC_V"],
        all_colors,
        all_side_by_side_output,
        title_prefix=""  # Each subplot will have just the label as title
    )

    logger.info("Analysis completed successfully")

if __name__ == "__main__":
    main()