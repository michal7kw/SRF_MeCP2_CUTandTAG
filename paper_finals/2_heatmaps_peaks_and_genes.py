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
import glob
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Color palette for different sample types - same as in R script
COLORS = {
    "NSC_ENDO": "#83CBEB",
    "NSC_EXO": "#0070C0",
    "NEU_ENDO": "#FF9999",
    "NEU_EXO": "#FF3300"
}

def extract_signal_reference_point(args):
    """Extract signal from a bigWig file for a given region around a reference point (TSS)"""
    bw_file, chrom, center, upstream, downstream, bin_size = args
    
    try:
        bw = pyBigWig.open(bw_file)
        
        # Calculate start and end positions
        start = max(0, center - upstream)
        end = center + downstream
        
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
        logger.error(f"Error processing {bw_file} at {chrom}:{center}: {str(e)}")
        # Return zeros if there's an error
        return np.zeros(int((upstream + downstream) / bin_size))

def extract_signal_scale_regions(args):
    """Extract signal from a bigWig file for a gene body (TSS to TES) with scaled body length"""
    bw_file, chrom, start, end, upstream, downstream, body_bins, bin_size = args
    
    try:
        bw = pyBigWig.open(bw_file)
        
        # Calculate actual regions
        upstream_start = max(0, start - upstream)
        upstream_end = start
        body_start = start
        body_end = end
        downstream_start = end
        downstream_end = end + downstream
        
        # Initialize signal arrays for each region
        upstream_signal = np.zeros(int(upstream / bin_size))
        body_signal = np.zeros(body_bins)
        downstream_signal = np.zeros(int(downstream / bin_size))
        
        # Extract upstream signal
        if upstream > 0:
            upstream_bins = np.linspace(upstream_start, upstream_end, int(upstream / bin_size) + 1, dtype=int)
            for i in range(len(upstream_bins) - 1):
                bin_start = upstream_bins[i]
                bin_end = upstream_bins[i + 1]
                try:
                    upstream_signal[i] = bw.stats(chrom, bin_start, bin_end, type="mean")[0] or 0
                except:
                    upstream_signal[i] = 0
        
        # Extract body signal (scaled)
        if body_end > body_start:
            # Create evenly spaced bins across the gene body
            body_bin_size = (body_end - body_start) / body_bins
            for i in range(body_bins):
                bin_start = body_start + int(i * body_bin_size)
                bin_end = body_start + int((i + 1) * body_bin_size)
                if bin_end > bin_start:
                    try:
                        body_signal[i] = bw.stats(chrom, bin_start, bin_end, type="mean")[0] or 0
                    except:
                        body_signal[i] = 0
        
        # Extract downstream signal
        if downstream > 0:
            downstream_bins = np.linspace(downstream_start, downstream_end, int(downstream / bin_size) + 1, dtype=int)
            for i in range(len(downstream_bins) - 1):
                bin_start = downstream_bins[i]
                bin_end = downstream_bins[i + 1]
                try:
                    downstream_signal[i] = bw.stats(chrom, bin_start, bin_end, type="mean")[0] or 0
                except:
                    downstream_signal[i] = 0
        
        # Combine all signals
        signal = np.concatenate([upstream_signal, body_signal, downstream_signal])
        
        bw.close()
        return signal
    
    except Exception as e:
        logger.error(f"Error processing {bw_file} at {chrom}:{start}-{end}: {str(e)}")
        # Return zeros if there's an error
        total_bins = int(upstream / bin_size) + body_bins + int(downstream / bin_size)
        return np.zeros(total_bins)

def compute_matrix_reference_point(bw_files, regions, upstream=5000, downstream=5000, bin_size=50, n_processes=8):
    """Compute signal matrix for multiple bigWig files and regions around a reference point (TSS)"""
    start_time = time.time()
    logger.info(f"Computing reference-point matrix for {len(bw_files)} bigWig files and {len(regions)} regions")
    
    # Prepare arguments for parallel processing
    args_list = []
    for region in regions:
        chrom = region.chrom
        # For TSS, use the start position as the reference point
        center = region.start
        
        for bw_file in bw_files:
            args_list.append((bw_file, chrom, center, upstream, downstream, bin_size))
    
    # Process in parallel
    with Pool(processes=n_processes) as pool:
        results = list(tqdm(pool.imap(extract_signal_reference_point, args_list), total=len(args_list)))
    
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
    
    logger.info(f"Matrix computation completed in {time.time() - start_time:.2f} seconds")
    
    return matrix

def compute_matrix_scale_regions(bw_files, regions, upstream=5000, downstream=5000, body_bins=100, bin_size=50, n_processes=8):
    """Compute signal matrix for multiple bigWig files and gene bodies with scaled body length"""
    start_time = time.time()
    logger.info(f"Computing scale-regions matrix for {len(bw_files)} bigWig files and {len(regions)} regions")
    
    # Prepare arguments for parallel processing
    args_list = []
    for region in regions:
        chrom = region.chrom
        start = region.start
        end = region.end
        
        for bw_file in bw_files:
            args_list.append((bw_file, chrom, start, end, upstream, downstream, body_bins, bin_size))
    
    # Process in parallel
    with Pool(processes=n_processes) as pool:
        results = list(tqdm(pool.imap(extract_signal_scale_regions, args_list), total=len(args_list)))
    
    # Reshape results into a matrix [n_regions × n_files, n_bins]
    n_regions = len(regions)
    n_files = len(bw_files)
    n_bins = int(upstream / bin_size) + body_bins + int(downstream / bin_size)
    
    # Reshape and reorder to get [n_files, n_regions, n_bins]
    matrix = np.zeros((n_files, n_regions, n_bins))
    
    for i in range(n_regions):
        for j in range(n_files):
            idx = i * n_files + j
            if idx < len(results):
                signal = results[idx]
                if len(signal) == n_bins:
                    matrix[j, i, :] = signal
    
    logger.info(f"Matrix computation completed in {time.time() - start_time:.2f} seconds")
    
    return matrix

def plot_tss_profile(matrices, labels, colors, output_file, title="MeCP2 around TSS",
                    xlabel="Distance from TSS (bp)", ylabel="Signal"):
    """Create a profile plot for TSS-centered data"""
    plt.figure(figsize=(10, 6))
    
    # Create x-axis values representing distance from TSS
    x = np.linspace(-5000, 5000, matrices[0].shape[0])
    
    for i, (matrix, label) in enumerate(zip(matrices, labels)):
        plt.plot(x, matrix, label=label, color=colors[i], linewidth=2)
    
    plt.title(title, fontsize=14)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.legend(loc='upper right')
    plt.grid(alpha=0.3)
    
    # Add vertical line at TSS
    plt.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    logger.info(f"TSS profile plot saved to {output_file}")
    
    return plt.gcf()

def plot_gene_body_profile(matrices, labels, colors, output_file, title="MeCP2 gene body coverage",
                         xlabel="", ylabel="Signal"):
    """Create a profile plot for gene body data (TSS to TES)"""
    plt.figure(figsize=(10, 6))
    
    # For gene body plots, we need to create a custom x-axis
    n_bins = matrices[0].shape[0]
    
    # Assuming the structure is [upstream bins, body bins, downstream bins]
    upstream_bins = 100  # 5000bp / 50bp
    body_bins = 100      # Fixed body bins
    downstream_bins = 100  # 5000bp / 50bp
    
    # Create x-axis labels
    x_labels = []
    x_ticks = []
    
    # Add upstream region
    x_ticks.append(0)
    x_labels.append("-5kb")
    
    # Add TSS
    x_ticks.append(upstream_bins)
    x_labels.append("TSS")
    
    # Add TES
    x_ticks.append(upstream_bins + body_bins)
    x_labels.append("TES")
    
    # Add downstream end
    x_ticks.append(n_bins - 1)
    x_labels.append("+5kb")
    
    # Create continuous x-axis
    x = np.arange(n_bins)
    
    for i, (matrix, label) in enumerate(zip(matrices, labels)):
        plt.plot(x, matrix, label=label, color=colors[i], linewidth=2)
    
    plt.title(title, fontsize=14)
    plt.ylabel(ylabel, fontsize=12)
    plt.legend(loc='upper right')
    plt.grid(alpha=0.3)
    
    # Set custom x-axis
    plt.xticks(x_ticks, x_labels)
    
    # Add vertical lines at TSS and TES
    plt.axvline(x=upstream_bins, color='gray', linestyle='--', alpha=0.5)
    plt.axvline(x=upstream_bins + body_bins, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    logger.info(f"Gene body profile plot saved to {output_file}")
    
    return plt.gcf()

def plot_side_by_side_tss(matrices, labels, colors, output_file, title_prefix="MeCP2"):
    """Create side-by-side profile plots for TSS data"""
    fig, axs = plt.subplots(1, len(matrices), figsize=(15, 5), sharey=True)
    
    # Create x-axis values representing distance from TSS
    x = np.linspace(-5000, 5000, matrices[0].shape[0])
    
    for i, (matrix, label, color) in enumerate(zip(matrices, labels, colors)):
        axs[i].plot(x, matrix, color=color, linewidth=2)
        axs[i].set_title(f"{title_prefix} {label}", fontsize=14)
        axs[i].grid(alpha=0.3)
        axs[i].axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        
        # Add x-axis label
        if i == 0:
            axs[i].set_ylabel("Signal", fontsize=12)
        
        # Format x-axis
        axs[i].set_xlabel("Distance from TSS (bp)", fontsize=12)
        axs[i].set_xticks([-5000, 0, 5000])
        axs[i].set_xticklabels(["-5.0Kb", "TSS", "5.0Kb"])
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    logger.info(f"Side-by-side TSS plot saved to {output_file}")
    
    return fig

def plot_side_by_side_gene_body(matrices, labels, colors, output_file, title_prefix="MeCP2"):
    """Create side-by-side profile plots for gene body data"""
    fig, axs = plt.subplots(1, len(matrices), figsize=(15, 5), sharey=True)
    
    # For gene body plots, we need to create a custom x-axis
    n_bins = matrices[0].shape[0]
    
    # Assuming the structure is [upstream bins, body bins, downstream bins]
    upstream_bins = 100  # 5000bp / 50bp
    body_bins = 100      # Fixed body bins
    downstream_bins = 100  # 5000bp / 50bp
    
    # Create x-axis labels
    x_labels = ["-5kb", "TSS", "TES", "+5kb"]
    x_ticks = [0, upstream_bins, upstream_bins + body_bins, n_bins - 1]
    
    # Create continuous x-axis
    x = np.arange(n_bins)
    
    for i, (matrix, label, color) in enumerate(zip(matrices, labels, colors)):
        axs[i].plot(x, matrix, color=color, linewidth=2)
        axs[i].set_title(f"{title_prefix} {label}", fontsize=14)
        axs[i].grid(alpha=0.3)
        
        # Add vertical lines at TSS and TES
        axs[i].axvline(x=upstream_bins, color='gray', linestyle='--', alpha=0.5)
        axs[i].axvline(x=upstream_bins + body_bins, color='gray', linestyle='--', alpha=0.5)
        
        # Add y-axis label for first plot only
        if i == 0:
            axs[i].set_ylabel("Signal", fontsize=12)
        
        # Set custom x-axis
        axs[i].set_xticks(x_ticks)
        axs[i].set_xticklabels(x_labels)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    logger.info(f"Side-by-side gene body plot saved to {output_file}")
    
    return fig

def load_peaks(peaks_dir, sample_name):
    """Load peak files for a given sample"""
    peak_file = None
    for ext in ["narrowPeak", "broadPeak"]:
        temp_peak_file = os.path.join(peaks_dir, f"{sample_name}_{ext}s.filtered.{ext}")
        if os.path.exists(temp_peak_file):
            peak_file = temp_peak_file
            break
    
    if not os.path.exists(peak_file):
        logger.error(f"Peak file not found: {peak_file}")
        return None
    
    logger.info(f"Loading peaks from {peak_file}")
    peaks = BedTool(peak_file)
    logger.info(f"Loaded {len(peaks)} peaks for {sample_name}")
    
    return peaks

def main():
    # Define sample groups - same as in R script
    sample_groups = {
        "neuron_endo": ["NeuM2", "NeuM3"],
        "neuron_exo": ["NeuV1", "NeuV2", "NeuV3"],
        "nsc_endo": ["NSCM1", "NSCM2", "NSCM3"],
        "nsc_exo": ["NSCv1", "NSCv2", "NSCv3"]
    }
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="MeCP2 peaks and gene body profiling")
    parser.add_argument("--bigwig-dir", required=True, help="Directory containing bigWig files")
    parser.add_argument("--peaks-dir", required=True, help="Directory containing peak files")
    parser.add_argument("--gene-regions-bed", required=True, help="BED file with gene regions (TSS to TES)")
    parser.add_argument("--tss-bed", required=True, help="BED file with TSS regions")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--upstream", type=int, default=5000, help="Upstream distance (bp)")
    parser.add_argument("--downstream", type=int, default=5000, help="Downstream distance (bp)")
    parser.add_argument("--bin-size", type=int, default=50, help="Bin size (bp)")
    parser.add_argument("--body-bins", type=int, default=100, help="Number of bins for gene body")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load gene regions and TSS
    logger.info(f"Loading gene regions from {args.gene_regions_bed}")
    gene_regions = BedTool(args.gene_regions_bed)
    logger.info(f"Loaded {len(gene_regions)} gene regions")
    
    logger.info(f"Loading TSS regions from {args.tss_bed}")
    tss_regions = BedTool(args.tss_bed)
    logger.info(f"Loaded {len(tss_regions)} TSS regions")
    
    # Process each sample group for TSS profiles
    tss_matrices = {}
    for group_name, sample_files in sample_groups.items():
        # Get full paths to bigWig files
        bw_files = [os.path.join(args.bigwig_dir, f"{sample}.bw") for sample in sample_files]
        
        # Compute matrix for TSS
        logger.info(f"Computing TSS matrix for {group_name}")
        matrix = compute_matrix_reference_point(
            bw_files, 
            tss_regions, 
            upstream=args.upstream, 
            downstream=args.downstream,
            bin_size=args.bin_size,
            n_processes=args.threads
        )
        
        # Average over samples and regions
        tss_matrices[group_name] = np.mean(matrix, axis=(0, 1))
    
    # Process each sample group for gene body profiles
    gene_body_matrices = {}
    for group_name, sample_files in sample_groups.items():
        # Get full paths to bigWig files
        bw_files = [os.path.join(args.bigwig_dir, f"{sample}.bw") for sample in sample_files]
        
        # Compute matrix for gene bodies
        logger.info(f"Computing gene body matrix for {group_name}")
        matrix = compute_matrix_scale_regions(
            bw_files, 
            gene_regions, 
            upstream=args.upstream, 
            downstream=args.downstream,
            body_bins=args.body_bins,
            bin_size=args.bin_size,
            n_processes=args.threads
        )
        
        # Average over samples and regions
        gene_body_matrices[group_name] = np.mean(matrix, axis=(0, 1))
    
    # Create individual TSS plots
    for cell_type in ["neuron", "nsc"]:
        endo_matrix = tss_matrices[f"{cell_type}_endo"]
        exo_matrix = tss_matrices[f"{cell_type}_exo"]
        
        output_file = os.path.join(args.output_dir, f"{cell_type}_endo_exo_tss_profile.png")
        plot_tss_profile(
            [endo_matrix, exo_matrix],
            ["Endogenous", "Exogenous"],
            [COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_ENDO"], 
             COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_EXO"]],
            output_file,
            title=f"{cell_type.title()} MeCP2 around TSS"
        )
        
        # Create side-by-side TSS plots
        side_by_side_output = os.path.join(args.output_dir, f"{cell_type}_endo_exo_tss_side_by_side.png")
        plot_side_by_side_tss(
            [endo_matrix, exo_matrix],
            ["Endogenous", "Exogenous"],
            [COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_ENDO"], 
             COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_EXO"]],
            side_by_side_output,
            title_prefix=f"{cell_type.title()}"
        )
    
    # Create individual gene body plots
    for cell_type in ["neuron", "nsc"]:
        endo_matrix = gene_body_matrices[f"{cell_type}_endo"]
        exo_matrix = gene_body_matrices[f"{cell_type}_exo"]
        
        output_file = os.path.join(args.output_dir, f"{cell_type}_endo_exo_gene_body_profile.png")
        plot_gene_body_profile(
            [endo_matrix, exo_matrix],
            ["Endogenous", "Exogenous"],
            [COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_ENDO"], 
             COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_EXO"]],
            output_file,
            title=f"{cell_type.title()} MeCP2 gene body coverage"
        )
        
        # Create side-by-side gene body plots
        side_by_side_output = os.path.join(args.output_dir, f"{cell_type}_endo_exo_gene_body_side_by_side.png")
        plot_side_by_side_gene_body(
            [endo_matrix, exo_matrix],
            ["Endogenous", "Exogenous"],
            [COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_ENDO"], 
             COLORS[f"{'NEU' if cell_type == 'neuron' else 'NSC'}_EXO"]],
            side_by_side_output,
            title_prefix=f"{cell_type.title()}"
        )
    
    # Create combined TSS plot
    all_tss_matrices = [
        tss_matrices["neuron_endo"],
        tss_matrices["neuron_exo"],
        tss_matrices["nsc_endo"],
        tss_matrices["nsc_exo"]
    ]
    
    all_labels = ["Neuron Endo", "Neuron Exo", "NSC Endo", "NSC Exo"]
    all_colors = [COLORS["NEU_ENDO"], COLORS["NEU_EXO"], COLORS["NSC_ENDO"], COLORS["NSC_EXO"]]
    
    output_file = os.path.join(args.output_dir, "all_samples_tss_profile.png")
    plot_tss_profile(
        all_tss_matrices,
        all_labels,
        all_colors,
        output_file,
        title="MeCP2 around TSS"
    )
    
    # Create combined gene body plot
    all_gene_body_matrices = [
        gene_body_matrices["neuron_endo"],
        gene_body_matrices["neuron_exo"],
        gene_body_matrices["nsc_endo"],
        gene_body_matrices["nsc_exo"]
    ]
    
    output_file = os.path.join(args.output_dir, "all_samples_gene_body_profile.png")
    plot_gene_body_profile(
        all_gene_body_matrices,
        all_labels,
        all_colors,
        output_file,
        title="MeCP2 gene body coverage"
    )
    
    # Process peaks for each sample
    for group_name, sample_files in sample_groups.items():
        for sample in sample_files:
            peaks = load_peaks(args.peaks_dir, sample)
            if peaks is None:
                continue
            
            # Get bigWig file for this sample
            bw_file = os.path.join(args.bigwig_dir, f"{sample}.bw")
            
            # Compute matrix for peaks
            logger.info(f"Computing peak matrix for {sample}")
            matrix = compute_matrix_reference_point(
                [bw_file], 
                peaks, 
                upstream=args.upstream, 
                downstream=args.downstream,
                bin_size=args.bin_size,
                n_processes=args.threads
            )
            
            # Average over regions
            peak_matrix = np.mean(matrix, axis=(0, 1))
            
            # Plot peak profile
            output_file = os.path.join(args.output_dir, f"{sample}_peak_profile.png")
            
            # Determine color based on sample name
            if sample.startswith("Neu"):
                if "M" in sample:
                    color = COLORS["NEU_ENDO"]
                    sample_type = "Neuron Endogenous"
                else:
                    color = COLORS["NEU_EXO"]
                    sample_type = "Neuron Exogenous"
            else:  # NSC
                if "M" in sample:
                    color = COLORS["NSC_ENDO"]
                    sample_type = "NSC Endogenous"
                else:
                    color = COLORS["NSC_EXO"]
                    sample_type = "NSC Exogenous"
            
            plot_tss_profile(
                [peak_matrix],
                [sample],
                [color],
                output_file,
                title=f"{sample_type} MeCP2 around peaks",
                xlabel="Distance from peak center (bp)"
            )
    
    logger.info("Analysis completed successfully")

if __name__ == "__main__":
    main() 