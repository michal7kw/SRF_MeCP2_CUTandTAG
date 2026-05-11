#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pyBigWig
import pysam
from typing import List, Tuple, Dict
import os
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def load_and_clean_regions(file_path: str) -> pd.DataFrame:
    """Load and clean region data, ensuring proper coordinate formatting."""
    df = pd.read_csv(file_path)
    
    # Drop rows with NaN values in essential columns
    df = df.dropna(subset=['seqnames', 'start', 'end'])
    
    # Convert coordinates to integers
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    
    # Ensure chromosome names are strings without decimal points
    df['seqnames'] = df['seqnames'].astype(str).str.replace('.0', '')
    
    return df

def calculate_meta_profile(regions: pd.DataFrame,
                         bigwig_file: str,
                         window_size: int = 2000,
                         bin_size: int = 50) -> np.ndarray:
    """Calculate meta-profile around region centers."""
    n_bins = window_size * 2 // bin_size
    profile_matrix = np.zeros((len(regions), n_bins))
    
    with pyBigWig.open(bigwig_file) as bw:
        for i, (_, region) in enumerate(regions.iterrows()):
            try:
                # Calculate window coordinates
                center = (region['start'] + region['end']) // 2
                start = max(0, center - window_size)
                end = center + window_size
                
                # Get chromosome name without any decimal points
                chrom = str(region['seqnames']).replace('.0', '')
                
                # Verify coordinates
                if not chrom or chrom == 'nan' or start >= end:
                    continue
                
                # Get raw signal using stats
                values = bw.stats(chrom, start, end, nBins=n_bins, type="mean")
                if values and not all(v is None for v in values):
                    # Replace None values with 0 and store
                    profile_matrix[i] = [v if v is not None else 0 for v in values]
            except Exception as e:
                print(f"Error processing region {region['seqnames']}:{start}-{end}: {str(e)}")
                continue
    
    # Average across all regions, excluding zeros
    mean_profile = np.zeros(n_bins)
    for i in range(n_bins):
        col_values = profile_matrix[:, i]
        valid_values = col_values[col_values != 0]
        mean_profile[i] = np.mean(valid_values) if len(valid_values) > 0 else 0
    
    return mean_profile

def analyze_methylation_profiles(gene_lists: Dict[str, str],
                               medip_dir: str,
                               smarcb1_dir: str,
                               output_dir: str):
    """Analyze methylation and SMARCB1 profiles around peaks."""
    
    profiles_dir = os.path.join(output_dir, 'meta_profiles')
    os.makedirs(profiles_dir, exist_ok=True)
    
    # Get file paths
    ip_files = [os.path.join(medip_dir, f"Medip_PP_output_r{i}.bw") for i in range(1, 4)]
    input_files = [os.path.join(medip_dir, f"Medip_PP_input_r{i}.bw") for i in range(1, 4)]
    bm_file = os.path.join(smarcb1_dir, "BM3_RPKM.bw")
    bg_files = [os.path.join(smarcb1_dir, f"BG{i}_RPKM.bw") for i in range(1, 4)]
    
    # Calculate profiles for each category
    for category, file_path in gene_lists.items():
        print(f"\nProcessing {category} regions...")
        
        # Load and clean regions
        regions = load_and_clean_regions(file_path)
        print(f"Number of valid regions after cleaning: {len(regions)}")
        
        if len(regions) == 0:
            print(f"Warning: No valid regions found for {category}")
            continue
        
        # Calculate average profiles
        print("Calculating IP profiles...")
        ip_profiles = []
        for ip_file in ip_files:
            if not os.path.exists(ip_file):
                print(f"Warning: File not found - {ip_file}")
                continue
                
            print(f"Processing {os.path.basename(ip_file)}")
            profile = calculate_meta_profile(regions, ip_file)
            if np.any(profile != 0):  # Only add non-zero profiles
                ip_profiles.append(profile)
        
        print("Calculating input profiles...")
        input_profiles = []
        for input_file in input_files:
            if not os.path.exists(input_file):
                print(f"Warning: File not found - {input_file}")
                continue
                
            print(f"Processing {os.path.basename(input_file)}")
            profile = calculate_meta_profile(regions, input_file)
            if np.any(profile != 0):  # Only add non-zero profiles
                input_profiles.append(profile)
        
        if not ip_profiles or not input_profiles:
            print(f"Warning: No valid profiles found for {category}")
            continue
        
        # Average replicates
        mean_ip = np.mean(ip_profiles, axis=0)
        mean_input = np.mean(input_profiles, axis=0)
        
        # Calculate normalized methylation profile
        methylation_profile = np.zeros_like(mean_ip)
        valid_idx = mean_input != 0
        methylation_profile[valid_idx] = mean_ip[valid_idx] / mean_input[valid_idx]
        
        # Calculate SMARCB1 profiles
        print("Calculating SMARCB1 profiles...")
        if not os.path.exists(bm_file):
            print(f"Warning: File not found - {bm_file}")
            continue
            
        print(f"Processing {os.path.basename(bm_file)}")
        bm_profile = calculate_meta_profile(regions, bm_file)
        
        bg_profiles = []
        for bg_file in bg_files:
            if not os.path.exists(bg_file):
                print(f"Warning: File not found - {bg_file}")
                continue
                
            print(f"Processing {os.path.basename(bg_file)}")
            profile = calculate_meta_profile(regions, bg_file)
            if np.any(profile != 0):  # Only add non-zero profiles
                bg_profiles.append(profile)
        
        if not bg_profiles:
            print(f"Warning: No valid BG profiles found for {category}")
            continue
        
        mean_bg = np.mean(bg_profiles, axis=0)
        
        # Calculate enrichment profile
        smarcb1_profile = np.zeros_like(bm_profile)
        valid_idx = mean_bg != 0
        smarcb1_profile[valid_idx] = bm_profile[valid_idx] / mean_bg[valid_idx]
        
        # Plot profiles
        plot_meta_profiles(methylation_profile, smarcb1_profile, category, profiles_dir)
        
        # Save profile data
        np.save(os.path.join(profiles_dir, f'{category}_methylation_profile.npy'), methylation_profile)
        np.save(os.path.join(profiles_dir, f'{category}_smarcb1_profile.npy'), smarcb1_profile)
        
        # Print some statistics
        print(f"\nProfile statistics for {category}:")
        print(f"Methylation profile - Mean: {np.mean(methylation_profile):.4f}, Max: {np.max(methylation_profile):.4f}")
        print(f"SMARCB1 profile - Mean: {np.mean(smarcb1_profile):.4f}, Max: {np.max(smarcb1_profile):.4f}")

def plot_meta_profiles(methylation_profile: np.ndarray,
                      smarcb1_profile: np.ndarray,
                      category: str,
                      output_dir: str):
    """Create meta-profile plots."""
    
    # Create x-axis coordinates (distance from center)
    x = np.linspace(-2000, 2000, len(methylation_profile))
    
    # Plot methylation profile
    plt.figure(figsize=(10, 6))
    plt.plot(x, methylation_profile, label='Methylation', color='blue')
    plt.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    plt.title(f'Methylation Profile Around {category} Peaks')
    plt.xlabel('Distance from Peak Center (bp)')
    plt.ylabel('Normalized Methylation Signal')
    plt.legend()
    plt.savefig(os.path.join(output_dir, f'{category}_methylation_profile.png'))
    plt.close()
    
    # Plot SMARCB1 profile
    plt.figure(figsize=(10, 6))
    plt.plot(x, smarcb1_profile, label='SMARCB1', color='red')
    plt.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    plt.title(f'SMARCB1 Enrichment Around {category} Peaks')
    plt.xlabel('Distance from Peak Center (bp)')
    plt.ylabel('SMARCB1 Enrichment (RPKM)')
    plt.legend()
    plt.savefig(os.path.join(output_dir, f'{category}_smarcb1_profile.png'))
    plt.close()
    
    # Combined plot
    plt.figure(figsize=(12, 6))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Methylation subplot
    ax1.plot(x, methylation_profile, color='blue')
    ax1.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax1.set_title('Methylation Profile')
    ax1.set_xlabel('Distance from Peak Center (bp)')
    ax1.set_ylabel('Normalized Methylation Signal')
    
    # SMARCB1 subplot
    ax2.plot(x, smarcb1_profile, color='red')
    ax2.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax2.set_title('SMARCB1 Enrichment')
    ax2.set_xlabel('Distance from Peak Center (bp)')
    ax2.set_ylabel('SMARCB1 Enrichment (RPKM)')
    
    plt.suptitle(f'Signal Profiles Around {category} Peaks')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{category}_combined_profile.png'))
    plt.close()

def analyze_peak_characteristics(gene_lists: Dict[str, pd.DataFrame],
                               output_dir: str):
    """Analyze peak characteristics across categories."""
    
    results = []
    for category, df in gene_lists.items():
        # Calculate peak lengths
        df['peak_length'] = df['end'] - df['start']
        
        # Calculate basic statistics
        stats = {
            'category': category,
            'mean_length': df['peak_length'].mean(),
            'median_length': df['peak_length'].median(),
            'std_length': df['peak_length'].std(),
            'total_peaks': len(df)
        }
        results.append(stats)
    
    # Save statistics
    stats_df = pd.DataFrame(results)
    stats_df.to_csv(os.path.join(output_dir, 'peak_characteristics.csv'), index=False)
    
    # Plot peak length distribution
    plt.figure(figsize=(10, 6))
    for category, df in gene_lists.items():
        sns.kdeplot(data=df['peak_length'], label=category)
    plt.title('Peak Length Distribution by Category')
    plt.xlabel('Peak Length (bp)')
    plt.ylabel('Density')
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'peak_length_distribution.png'))
    plt.close()

def main():
    # Set paths
    base_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Cross_final/data"
    medip_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/DATA/MECP2/MEDIP/output_done/bigwig"
    smarcb1_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1/results/bigwig"
    output_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Cross_final/results/6analyze_methylation_profiles"
    
    # Define gene lists with file paths
    gene_lists = {
        'up': os.path.join(base_dir, 'extended_up.csv'),
        'down': os.path.join(base_dir, 'extended_down.csv'),
        'no_deg': os.path.join(base_dir, 'extended_no_deg.csv')
    }
    
    # Run meta-profile analysis
    analyze_methylation_profiles(gene_lists, medip_dir, smarcb1_dir, output_dir)
    
    # Analyze peak characteristics
    gene_lists = {
        category: load_and_clean_regions(file_path)
        for category, file_path in gene_lists.items()
    }
    analyze_peak_characteristics(gene_lists, output_dir)

if __name__ == "__main__":
    main()
