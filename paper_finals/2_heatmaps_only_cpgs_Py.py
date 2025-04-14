#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')  # Must be set before importing pyplot

import os
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib.pyplot as plt
import seaborn as sns
import concurrent.futures
try:
    import numba
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
try:
    import cupy as cp
    HAS_CUDA = True
except ImportError:
    HAS_CUDA = False
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
    n_bins = int((end - start) / bin_size)
    signal = np.zeros(n_bins)
    
    try:
        bw = pyBigWig.open(bw_file)
        chroms = bw.chroms()
        
        # Check if chromosome exists and bounds are valid
        if chrom not in chroms:
            logger.warning(f"Chromosome {chrom} not found in {bw_file}")
            bw.close()
            return signal
            
        chrom_length = chroms[chrom]
        start = max(0, start)
        end = min(chrom_length, end)
        
        if start >= end:
            logger.warning(f"Invalid bounds {start}-{end} for {chrom} in {bw_file}")
            bw.close()
            return signal
            
        try:
            # First try vectorized approach
            means = bw.stats(chrom, start, end, type="mean", nBins=n_bins)
            if means is not None:
                signal = np.nan_to_num(means, nan=0.0, posinf=0.0, neginf=0.0)
            else:
                # Fall back to per-bin if vectorized fails
                bins = np.linspace(start, end, n_bins + 1, dtype=int)
                for i in range(n_bins):
                    try:
                        signal[i] = bw.stats(chrom, bins[i], bins[i+1], type="mean")[0] or 0
                    except:
                        signal[i] = 0
                        
        except Exception as e:
            logger.warning(f"Vectorized extraction failed for {bw_file}: {str(e)}")
            # Fall back to per-bin if vectorized fails
            bins = np.linspace(start, end, n_bins + 1, dtype=int)
            for i in range(n_bins):
                try:
                    signal[i] = bw.stats(chrom, bins[i], bins[i+1], type="mean")[0] or 0
                except Exception as e:
                    if "invalid" in str(e).lower():
                        logger.debug(f"Bin {i} error in {bw_file}: {str(e)}")
                    signal[i] = 0
                    
        bw.close()
        return signal
        
    except Exception as e:
        logger.error(f"Error opening {bw_file}: {str(e)}")
        return signal

def compute_matrix(bw_files, regions, upstream=3000, downstream=3000, bin_size=50, n_processes=8):
    """Compute signal matrix for multiple bigWig files and regions"""
    start_time = time.time()
    logger.info(f"Computing matrix for {len(bw_files)} bigWig files and {len(regions)} regions")
    
    n_bins = int((upstream + downstream) / bin_size)
    n_regions = len(regions)
    n_files = len(bw_files)
    
    # Pre-allocate result matrix
    matrix = np.zeros((n_files, n_regions, n_bins))
    
    # Process regions in chunks to reduce memory usage
    chunk_size = min(1000, n_regions)
    
    for chunk_start in range(0, n_regions, chunk_size):
        chunk_end = min(chunk_start + chunk_size, n_regions)
        chunk_regions = regions[chunk_start:chunk_end]
        
        # Prepare arguments for parallel processing
        args_list = []
        for region in chunk_regions:
            chrom = region.chrom
            center = int((region.start + region.end) / 2)
            start = center - upstream
            end = center + downstream
            
            for bw_file in bw_files:
                args_list.append((bw_file, chrom, start, end, bin_size))
        
        # Process in parallel with concurrent.futures for better performance
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_processes) as executor:
            results = list(tqdm(executor.map(extract_signal, args_list), total=len(args_list)))
        
        # Fill matrix chunk
        for i, region_idx in enumerate(range(chunk_start, chunk_end)):
            for j in range(n_files):
                idx = i * n_files + j
                if idx < len(results):
                    signal = results[idx]
                    if len(signal) == n_bins:
                        matrix[j, region_idx, :] = signal
    
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
    fig, axs = plt.subplots(1, len(matrices), figsize=(15, 5))
    
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
    
    # Process and plot each sample group incrementally
    for cell_type in ["neuron", "nsc"]:
        # Process endogenous first
        endo_files = sample_groups[f"{cell_type}_endo"]
        bw_files = [os.path.join(args.bigwig_dir, f) for f in endo_files]
        
        logger.info(f"Computing matrix for {cell_type} endogenous")
        endo_matrix = compute_matrix(
            bw_files,
            cpg_islands,
            upstream=args.upstream,
            downstream=args.downstream,
            bin_size=args.bin_size,
            n_processes=args.threads
        )
        endo_avg = np.mean(endo_matrix, axis=0)
        
        # Process exogenous next
        exo_files = sample_groups[f"{cell_type}_exo"]
        bw_files = [os.path.join(args.bigwig_dir, f) for f in exo_files]
        
        logger.info(f"Computing matrix for {cell_type} exogenous")
        exo_matrix = compute_matrix(
            bw_files,
            cpg_islands,
            upstream=args.upstream,
            downstream=args.downstream,
            bin_size=args.bin_size,
            n_processes=args.threads
        )
        exo_avg = np.mean(exo_matrix, axis=0)
        
        # Save individual plots and profiles
        output_file = os.path.join(args.output_dir, f"{cell_type}_endo_exo_cpg_profile.png")
        plot_profiles(
            [endo_avg, exo_avg],
            ["Endogenous", "Exogenous"],
            [COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_ENDO"],
             COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_EXO"]],
            output_file,
            title=f"{cell_type.title()} MeCP2 at CpG Islands"
        )
        
        side_by_side_output = os.path.join(args.output_dir, f"{cell_type}_endo_exo_side_by_side.png")
        plot_side_by_side(
            [endo_avg, exo_avg],
            ["Endogenous", "Exogenous"],
            [COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_ENDO"],
             COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_EXO"]],
            side_by_side_output,
            title_prefix=f"{cell_type.title()}"
        )
        
        # Save profile data for potential combined plot
        np.savetxt(os.path.join(args.output_dir, f"{cell_type}_endo_profile.txt"), endo_avg)
        np.savetxt(os.path.join(args.output_dir, f"{cell_type}_exo_profile.txt"), exo_avg)
    
    # Create combined plot if all components exist
    try:
        neuron_endo = np.loadtxt(os.path.join(args.output_dir, "neuron_endo_profile.txt"))
        neuron_exo = np.loadtxt(os.path.join(args.output_dir, "neuron_exo_profile.txt"))
        nsc_endo = np.loadtxt(os.path.join(args.output_dir, "nsc_endo_profile.txt"))
        nsc_exo = np.loadtxt(os.path.join(args.output_dir, "nsc_exo_profile.txt"))
        
        output_file = os.path.join(args.output_dir, "all_samples_cpg_profile.png")
        plot_profiles(
            [neuron_endo, neuron_exo, nsc_endo, nsc_exo],
            ["Neuron Endo", "Neuron Exo", "NSC Endo", "NSC Exo"],
            [COLORS["NEU_ENDO"], COLORS["NEU_EXO"], COLORS["NSC_ENDO"], COLORS["NSC_EXO"]],
            output_file,
            title="MeCP2 at CpG Islands"
        )
    except FileNotFoundError as e:
        logger.warning(f"Skipping combined plot: {str(e)}")

    logger.info("Analysis completed successfully")

if __name__ == "__main__":
    main()