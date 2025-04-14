#!/usr/bin/env python3

"""
python 0_2_plot_comparison_targeted_cpg.py --input-dir ../results_local/3cpg_metaprofiles_medip_targeted --output-dir ../results_local/4cpg_metaprofiles_medip_targeted_comparison
    # Optionally add --bw-files Medip_DP_output_merged Medip_N_output_merged ... if needed
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- Parameters from the previous script --- 
# These are needed to reconstruct the x-axis (bins)
# TODO: Consider passing these as arguments for flexibility
WINDOW_SIZE = 5000
# BIN_SIZE = 300 # Bin size might not be strictly needed if we use linspace
# --- End Parameters ---

def plot_comparison(nsc_data_path, neuron_data_path, output_plot_path, bw_label):
    """
    Loads metaprofile data (from .npy) for NSC and Neuron targets and plots a comparison.

    Args:
        nsc_data_path (str): Path to the NSC target metaprofile .npy file.
        neuron_data_path (str): Path to the Neuron target metaprofile .npy file.
        output_plot_path (str): Path to save the output comparison plot (.png).
        bw_label (str): Label for the BigWig file (used in the title).
    """
    try:
        # Load profile data directly from .npy files
        nsc_profile = np.load(nsc_data_path)
        neuron_profile = np.load(neuron_data_path)
        logging.info(f"Loaded data for {bw_label}: NSC ({len(nsc_profile)} points) from {os.path.basename(nsc_data_path)}, Neuron ({len(neuron_profile)} points) from {os.path.basename(neuron_data_path)}")

        # Check if profiles have the same length (implies same binning)
        if len(nsc_profile) != len(neuron_profile):
            logging.error(f"Error: Profile lengths differ for {bw_label} (NSC: {len(nsc_profile)}, Neuron: {len(neuron_profile)}). Cannot compare directly. Skipping.")
            return # Skip plotting this comparison

        num_points = len(nsc_profile)
        if num_points == 0:
            logging.warning(f"Warning: Profiles for {bw_label} have zero length. Skipping plot.")
            return

        # Generate the x-axis bins based on window size and number of points
        # Assumes points are evenly spaced across the window centered at 0
        bins = np.linspace(-WINDOW_SIZE / 2, WINDOW_SIZE / 2, num_points)

        plt.figure(figsize=(10, 6))
        plt.plot(bins, nsc_profile, color='lightblue', label='NSC Targets')
        plt.plot(bins, neuron_profile, color='lightcoral', label='Neuron Targets') # lightcoral is close to light red

        # --- Adjust x-axis ticks --- 
        # Keep it simple: label start, center, and end
        tick_positions = [bins[0], 0, bins[-1]]
        tick_labels = [f'{-WINDOW_SIZE/2000:.1f} kb', 'Center', f'{WINDOW_SIZE/2000:.1f} kb']
        # Use plt.xticks to set positions and labels
        # We need the indices corresponding to the tick positions for plt.xticks
        # Find the index closest to 0 for the center tick
        center_index = np.argmin(np.abs(bins))
        # Ensure start, center, end indices are used
        actual_tick_indices = [0, center_index, len(bins) - 1]
        actual_tick_positions = bins[actual_tick_indices]

        plt.xticks(actual_tick_positions, tick_labels)
        # --- End x-axis tick adjustment ---

        plt.xlabel("Position relative to CpG Target Center")
        plt.ylabel("Average Signal")
        plt.title(f"Targeted CpG Metaprofile Comparison ({bw_label})")
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        plt.savefig(output_plot_path)
        plt.close()
        logging.info(f"Saved comparison plot to {output_plot_path}")

    except FileNotFoundError as e:
        logging.error(f"Error: Input file not found - {e}. Skipping comparison for {bw_label}.")
    # No KeyError needed for .npy as we load the whole array
    except Exception as e:
        logging.error(f"An unexpected error occurred while plotting for {bw_label}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Plot comparison metaprofiles for targeted CpGs in NSC vs Neurons.")
    parser.add_argument("--input-dir", required=True, help="Directory containing the _cpg_metaprofile.npy files from the previous step.")
    parser.add_argument("--output-dir", required=True, help="Directory to save the comparison plots.")
    parser.add_argument("--bw-files", nargs='+', default=["Medip_DP_output_merged", "Medip_N_output_merged", "Medip_PP_output_merged"], help="List of BigWig file basenames to process (without .bw extension).")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Define the patterns for the target files using the correct suffix
    nsc_target_suffix = "_vs_nsc_cpg_targets_cpg_metaprofile.npy"
    neuron_target_suffix = "_vs_neu_cpg_targets_cpg_metaprofile.npy"

    processed_count = 0
    error_count = 0

    for bw_base in args.bw_files:
        logging.info(f"Processing BigWig: {bw_base}")
        nsc_file = os.path.join(args.input_dir, f"{bw_base}{nsc_target_suffix}")
        neuron_file = os.path.join(args.input_dir, f"{bw_base}{neuron_target_suffix}")
        output_file = os.path.join(args.output_dir, f"{bw_base}_NSC_vs_Neuron_target_comparison.png")

        if os.path.exists(nsc_file) and os.path.exists(neuron_file):
            plot_comparison(nsc_file, neuron_file, output_file, bw_base)
            processed_count += 1
        else:
            logging.warning(f"Skipping {bw_base}: Required input file(s) not found.")
            if not os.path.exists(nsc_file):
                logging.warning(f"  Missing: {nsc_file}")
            if not os.path.exists(neuron_file):
                logging.warning(f"  Missing: {neuron_file}")
            error_count += 1

    logging.info("="*30)
    logging.info(f"Processing complete. Successfully generated {processed_count} comparison plots.")
    if error_count > 0:
        logging.warning(f"{error_count} comparisons were skipped due to missing input files.")
    logging.info(f"Output plots saved in: {args.output_dir}")
    logging.info("="*30)


if __name__ == "__main__":
    main() 