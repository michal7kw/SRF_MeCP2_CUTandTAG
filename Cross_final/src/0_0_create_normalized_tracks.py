#!/usr/bin/env python3

# Example command:
# python 0_0_create_normalized_tracks.py --auto --dedup-dir ../../data/dedup --output-dir ../results_local/normalized_bigwig

import os
import sys
import argparse
import subprocess
import glob
import re

def create_normalized_bigwig(ip_bw_path, input_bw_path, output_bw, threads=8):
    """Create normalized bigwig file using bigwigCompare"""
    print(f"Comparing BigWig files: {os.path.basename(ip_bw_path)} vs {os.path.basename(input_bw_path)}")

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_bw), exist_ok=True)

    # Run bigwigCompare
    # Using log2 ratio, similar to the original bamCompare approach.
    # Added a pseudocount of 1 to avoid log2(0) issues.
    cmd = [
        'bigwigCompare',
        '--bigwig1', ip_bw_path,
        '--bigwig2', input_bw_path,
        '--operation', 'log2',
        '--pseudocount', '1', # Add pseudocount
        '--binSize', '10',
        '--numberOfProcessors', str(threads),
        '--outFileName', output_bw
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"Created normalized bigwig: {output_bw}")
    except subprocess.CalledProcessError as e:
        print(f"Error running bigwigCompare for {output_bw}:")
        print(f"Command: {' '.join(e.cmd)}")
        print(f"Stderr: {e.stderr}")
        print(f"Stdout: {e.stdout}")
        # Decide if you want to raise the error or just continue
        # raise e # Uncomment to stop execution on error

def auto_match_samples(dedup_dir, output_dir):
    """Automatically match IP and input BigWig samples from the dedup directory"""
    # Get all IP output BigWig samples
    ip_samples = glob.glob(os.path.join(dedup_dir, "*output*_markdup.bw")) # Changed to .bw

    pairs = []
    for ip_bw in ip_samples:
        basename = os.path.basename(ip_bw)
        # Extract condition and replicate from output BigWig filename
        match = re.search(r"Medip_(\w+)_output_(r\d+)_markdup.bw", basename) # Changed to .bw
        if match:
            condition, replicate = match.groups()
            # Find matching input BigWig
            input_bw = os.path.join(dedup_dir, f"Medip_{condition}_input_{replicate}_markdup.bw") # Changed to .bw
            if os.path.exists(input_bw):
                # Keep the same output naming convention
                output_bw_name = os.path.join(output_dir, f"{condition}_{replicate}_methylation_signal.bw")
                pairs.append((ip_bw, input_bw, output_bw_name))
            else:
                print(f"Warning: No matching input BigWig found for {basename}")
                print(f"Expected: {input_bw}")

    return pairs

def main():
    parser = argparse.ArgumentParser(description='Create normalized bigwig files by comparing IP and Input BigWig files using bigwigCompare.')
    # Updated help descriptions
    parser.add_argument('--ip', help='IP (output) BigWig file')
    parser.add_argument('--input', help='Input control BigWig file')
    parser.add_argument('--output', help='Output normalized bigwig file')
    parser.add_argument('--dedup-dir', help='Directory containing input BigWig files (generated from deduplicated BAMs)')
    parser.add_argument('--output-dir', help='Directory to save output normalized bigwig files')
    parser.add_argument('--threads', type=int, default=8, help='Number of threads to use (passed to bigwigCompare)')
    parser.add_argument('--auto', action='store_true', help='Automatically process all sample pairs found in dedup-dir')

    args = parser.parse_args()

    if args.auto and args.dedup_dir and args.output_dir:
        # Auto-match mode
        pairs = auto_match_samples(args.dedup_dir, args.output_dir)
        if not pairs:
            print(f"No matching sample pairs found in {args.dedup_dir}")
            return

        print(f"Found {len(pairs)} sample pairs to process")
        for ip_bw_path, input_bw_path, output_bw_file in pairs:
            create_normalized_bigwig(ip_bw_path, input_bw_path, output_bw_file, args.threads)

    elif args.ip and args.input and args.output:
        # Single sample mode
        create_normalized_bigwig(args.ip, args.input, args.output, args.threads)

    else:
        parser.print_help()
        print("\nEither provide --ip, --input, and --output for a single sample pair,")
        print("or use --auto with --dedup-dir and --output-dir to process all samples.")
        return

if __name__ == "__main__":
    main()