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
from Bio import SeqIO
from Bio.Seq import Seq

class EnhancedProfileAnalyzer:
    def __init__(self, 
                 genome_fasta: str,
                 window_size: int = 2000,
                 bin_size: int = 50):
        """
        Initialize the analyzer with genome sequence for GC/CpG calculations.
        
        Args:
            genome_fasta: Path to genome FASTA file
            window_size: Size of the window around peak center (default: 2000bp)
            bin_size: Size of bins for signal calculation (default: 50bp)
        """
        self.genome_fasta = genome_fasta
        self.window_size = window_size
        self.bin_size = bin_size
        self.n_bins = window_size * 2 // bin_size
        self.genome = pysam.FastaFile(genome_fasta)
    
    def calculate_gc_content(self, chrom: str, start: int, end: int) -> float:
        """Calculate GC content for a genomic region."""
        try:
            seq = self.genome.fetch(chrom, start, end).upper()
            gc_count = seq.count('G') + seq.count('C')
            total = len(seq)
            return gc_count / total if total > 0 else 0
        except Exception as e:
            print(f"Error calculating GC content for {chrom}:{start}-{end}: {str(e)}")
            return 0

    def calculate_cpg_density(self, chrom: str, start: int, end: int) -> float:
        """Calculate CpG density (observed/expected ratio) for a genomic region."""
        try:
            seq = self.genome.fetch(chrom, start, end).upper()
            length = len(seq)
            if length < 2:
                return 0
                
            # Count actual CpG dinucleotides
            cpg_count = seq.count('CG')
            
            # Calculate expected CpG based on C and G frequencies
            c_freq = (seq.count('C') / length)
            g_freq = (seq.count('G') / length)
            expected_cpg = length * c_freq * g_freq
            
            # Calculate observed/expected ratio
            if expected_cpg > 0:
                return (cpg_count / expected_cpg)
            return 0
            
        except Exception as e:
            print(f"Error calculating CpG density for {chrom}:{start}-{end}: {str(e)}")
            return 0

    def calculate_meta_profile(self,
                             regions: pd.DataFrame,
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
                    chrom = str(region['seqnames'])
                    
                    # Verify coordinates
                    if not chrom or chrom == 'nan' or start >= end:
                        continue
                    
                    # Get raw signal using stats
                    values = bw.stats(chrom, start, end, nBins=n_bins, type="mean")
                    if values and not all(v is None for v in values):
                        # Replace None values with 0 and store
                        profile_matrix[i] = [v if v is not None else 0 for v in values]
                except Exception as e:
                    print(f"Error processing region {chrom}:{start}-{end}: {str(e)}")
                    continue
        
        # Average across all regions, excluding zeros
        mean_profile = np.zeros(n_bins)
        for i in range(n_bins):
            col_values = profile_matrix[:, i]
            valid_values = col_values[col_values != 0]
            mean_profile[i] = np.mean(valid_values) if len(valid_values) > 0 else 0
        
        return mean_profile

    def calculate_methylation_profile(self,
                                   regions: pd.DataFrame,
                                   ip_file: str,
                                   input_file: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Calculate methylation profile with GC and CpG normalization.
        
        Returns:
            Tuple of (basic_profile, gc_normalized_profile, cpg_normalized_profile)
        """
        profile_matrix = np.zeros((len(regions), self.n_bins))
        gc_contents = np.zeros((len(regions), self.n_bins))
        cpg_densities = np.zeros((len(regions), self.n_bins))
        
        with pyBigWig.open(ip_file) as ip_bw, pyBigWig.open(input_file) as input_bw:
            for i, (_, region) in enumerate(regions.iterrows()):
                try:
                    # Calculate window coordinates
                    center = (region['start'] + region['end']) // 2
                    start = max(0, center - self.window_size)
                    end = center + self.window_size
                    
                    # Get chromosome name
                    chrom = str(region['seqnames'])
                    
                    # Verify coordinates
                    if not chrom or chrom == 'nan' or start >= end:
                        continue
                    
                    # Calculate signals
                    ip_values = ip_bw.stats(chrom, start, end, nBins=self.n_bins, type="mean")
                    input_values = input_bw.stats(chrom, start, end, nBins=self.n_bins, type="mean")
                    
                    if not ip_values or not input_values:
                        continue
                        
                    # Replace None values with 0
                    ip_values = [v if v is not None else 0 for v in ip_values]
                    input_values = [v if v is not None else 0 for v in input_values]
                    
                    # Calculate basic methylation profile
                    methylation = np.array(ip_values) / (np.array(input_values) + 1e-10)
                    profile_matrix[i] = methylation
                    
                    # Calculate GC content and CpG density for each bin
                    bin_size = (end - start) // self.n_bins
                    for j in range(self.n_bins):
                        bin_start = start + j * bin_size
                        bin_end = bin_start + bin_size
                        gc_contents[i, j] = self.calculate_gc_content(chrom, bin_start, bin_end)
                        cpg_densities[i, j] = self.calculate_cpg_density(chrom, bin_start, bin_end)
                    
                except Exception as e:
                    print(f"Error processing region {chrom}:{start}-{end}: {str(e)}")
                    continue
        
        # Calculate average profiles
        basic_profile = np.nanmean(profile_matrix, axis=0)
        
        # GC-normalized profile
        gc_norm_profile = np.zeros_like(basic_profile)
        for i in range(self.n_bins):
            valid_idx = gc_contents[:, i] > 0
            if np.any(valid_idx):
                gc_norm_profile[i] = np.mean(profile_matrix[valid_idx, i] / gc_contents[valid_idx, i])
        
        # CpG-normalized profile
        cpg_norm_profile = np.zeros_like(basic_profile)
        for i in range(self.n_bins):
            valid_idx = cpg_densities[:, i] > 0
            if np.any(valid_idx):
                cpg_norm_profile[i] = np.mean(profile_matrix[valid_idx, i] / cpg_densities[valid_idx, i])
        
        return basic_profile, gc_norm_profile, cpg_norm_profile

    def calculate_smarcb1_profile(self,
                                regions: pd.DataFrame,
                                bm_file: str,
                                bg_files: List[str]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate SMARCB1 profile using BG samples as input controls.
        
        Args:
            regions: DataFrame with genomic regions
            bm_file: Path to BM (MeCP2-expressing) bigWig file
            bg_files: List of paths to BG (input control) bigWig files
        
        Returns:
            Tuple of:
            - basic_profile: Simple BM/mean(BG) ratio
            - replicate_norm_profile: Average of BM/BG ratios across replicates
        """
        bm_matrix = np.zeros((len(regions), self.n_bins))
        bg_matrices = []
        
        # Process BM signal
        print("Processing BM signal...")
        with pyBigWig.open(bm_file) as bm_bw:
            for i, (_, region) in enumerate(regions.iterrows()):
                try:
                    center = (region['start'] + region['end']) // 2
                    start = max(0, center - self.window_size)
                    end = center + self.window_size
                    
                    chrom = str(region['seqnames'])
                    if not chrom or chrom == 'nan' or start >= end:
                        continue
                    
                    values = bm_bw.stats(chrom, start, end, nBins=self.n_bins, type="mean")
                    if values:
                        bm_matrix[i] = [v if v is not None else 0 for v in values]
                except Exception as e:
                    print(f"Error processing BM region {chrom}:{start}-{end}: {str(e)}")
                    continue
        
        # Process each BG replicate
        print("Processing BG signals...")
        for bg_file in bg_files:
            print(f"Processing {os.path.basename(bg_file)}")
            bg_matrix = np.zeros((len(regions), self.n_bins))
            
            with pyBigWig.open(bg_file) as bg_bw:
                for i, (_, region) in enumerate(regions.iterrows()):
                    try:
                        center = (region['start'] + region['end']) // 2
                        start = max(0, center - self.window_size)
                        end = center + self.window_size
                        
                        chrom = str(region['seqnames'])
                        if not chrom or chrom == 'nan' or start >= end:
                            continue
                        
                        values = bg_bw.stats(chrom, start, end, nBins=self.n_bins, type="mean")
                        if values:
                            bg_matrix[i] = [v if v is not None else 0 for v in values]
                    except Exception as e:
                        print(f"Error processing BG region {chrom}:{start}-{end}: {str(e)}")
                        continue
            
            bg_matrices.append(bg_matrix)
        
        # Calculate profiles
        basic_profile = np.zeros(self.n_bins)
        replicate_norm_profiles = []
        
        # Calculate normalization for each BG replicate
        for bg_matrix in bg_matrices:
            # Create a normalized profile for this replicate
            norm_profile = np.zeros(self.n_bins)
            for i in range(self.n_bins):
                valid_idx = (bg_matrix[:, i] > 0) & (bm_matrix[:, i] > 0)
                if np.any(valid_idx):
                    ratios = bm_matrix[valid_idx, i] / bg_matrix[valid_idx, i]
                    # Remove extreme outliers (ratios > 99th percentile)
                    ratios = ratios[ratios <= np.percentile(ratios, 99)]
                    if len(ratios) > 0:
                        norm_profile[i] = np.mean(ratios)
            replicate_norm_profiles.append(norm_profile)
        
        # Calculate basic profile (BM/mean(BG))
        for i in range(self.n_bins):
            mean_bg = np.mean([bg[:, i] for bg in bg_matrices], axis=0)
            valid_idx = (mean_bg > 0) & (bm_matrix[:, i] > 0)
            if np.any(valid_idx):
                ratios = bm_matrix[valid_idx, i] / mean_bg[valid_idx]
                # Remove extreme outliers
                ratios = ratios[ratios <= np.percentile(ratios, 99)]
                if len(ratios) > 0:
                    basic_profile[i] = np.mean(ratios)
        
        # Average the replicate-normalized profiles
        replicate_norm_profile = np.mean(replicate_norm_profiles, axis=0)
        
        return basic_profile, replicate_norm_profile

    def plot_methylation_profiles(self,
                                basic_profile: np.ndarray,
                                gc_norm_profile: np.ndarray,
                                cpg_norm_profile: np.ndarray,
                                category: str,
                                output_dir: str):
        """Plot methylation profiles with different normalizations."""
        x = np.linspace(-self.window_size, self.window_size, self.n_bins)
        
        plt.figure(figsize=(12, 6))
        plt.plot(x, basic_profile, label='Basic (IP/Input)', color='blue')
        plt.plot(x, gc_norm_profile, label='GC-normalized', color='red')
        plt.plot(x, cpg_norm_profile, label='CpG-normalized', color='green')
        plt.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        plt.title(f'Methylation Profiles Around {category} Peaks')
        plt.xlabel('Distance from Peak Center (bp)')
        plt.ylabel('Normalized Methylation Signal')
        plt.legend()
        plt.savefig(os.path.join(output_dir, f'{category}_methylation_profiles_comparison.png'))
        plt.close()
        
        # Save profile data
        np.save(os.path.join(output_dir, f'{category}_methylation_basic.npy'), basic_profile)
        np.save(os.path.join(output_dir, f'{category}_methylation_gc_norm.npy'), gc_norm_profile)
        np.save(os.path.join(output_dir, f'{category}_methylation_cpg_norm.npy'), cpg_norm_profile)

    def plot_smarcb1_profiles(self,
                            basic_profile: np.ndarray,
                            replicate_norm_profile: np.ndarray,
                            category: str,
                            output_dir: str):
        """Plot SMARCB1 profiles with different normalizations."""
        x = np.linspace(-self.window_size, self.window_size, self.n_bins)
        
        plt.figure(figsize=(12, 6))
        plt.plot(x, basic_profile, label='Basic (BM/mean(BG))', color='blue')
        plt.plot(x, replicate_norm_profile, label='Replicate-normalized', color='red')
        plt.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        plt.title(f'SMARCB1 Profiles Around {category} Peaks')
        plt.xlabel('Distance from Peak Center (bp)')
        plt.ylabel('SMARCB1 Enrichment')
        plt.legend()
        plt.savefig(os.path.join(output_dir, f'{category}_smarcb1_profiles_comparison.png'))
        plt.close()
        
        # Save profile data
        np.save(os.path.join(output_dir, f'{category}_smarcb1_basic.npy'), basic_profile)
        np.save(os.path.join(output_dir, f'{category}_smarcb1_replicate_norm.npy'), replicate_norm_profile)

def load_and_clean_regions(file_path: str) -> pd.DataFrame:
    """Load and clean region data, ensuring proper coordinate formatting."""
    df = pd.read_csv(file_path)
    
    # Drop rows with NaN values in essential columns
    df = df.dropna(subset=['seqnames', 'start', 'end'])
    
    # Clean chromosome names
    df['seqnames'] = df['seqnames'].astype(str).str.replace('.0', '')
    
    # Convert coordinates to integers (handling floating point values)
    df['start'] = df['start'].astype(float).astype(int)
    df['end'] = df['end'].astype(float).astype(int)
    
    # Validate coordinates
    df = df[df['start'] >= 0]
    df = df[df['end'] > df['start']]
    
    # Remove any remaining invalid entries
    df = df[~df['seqnames'].str.contains('nan')]
    
    print(f"Loaded {len(df)} valid regions after cleaning")
    return df

def main():
    # Set paths
    base_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Cross_final"
    genome_fasta = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA/mm10.fa"
    medip_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/DATA/MECP2/MEDIP/output_done/bigwig"
    smarcb1_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1/results/bigwig"
    output_dir = os.path.join(base_dir, "results/methylation_analysis/enhanced_profiles")
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize analyzer
    analyzer = EnhancedProfileAnalyzer(genome_fasta)
    
    # Process each category
    for category in ['up', 'down', 'no_deg']:
        print(f"\nProcessing {category} regions...")
        
        # Load and clean regions
        regions_file = os.path.join(base_dir, "data", f'extended_{category}.csv')
        if not os.path.exists(regions_file):
            print(f"Warning: File not found - {regions_file}")
            continue
            
        regions = load_and_clean_regions(regions_file)
        if len(regions) == 0:
            print(f"Warning: No valid regions found in {regions_file}")
            continue
        
        # Calculate methylation profiles
        print("Calculating methylation profiles...")
        basic_meth, gc_norm_meth, cpg_norm_meth = analyzer.calculate_methylation_profile(
            regions,
            os.path.join(medip_dir, "Medip_PP_output_r1.bw"),
            os.path.join(medip_dir, "Medip_PP_input_r1.bw")
        )
        
        # Plot methylation profiles
        analyzer.plot_methylation_profiles(
            basic_meth, gc_norm_meth, cpg_norm_meth,
            category, output_dir
        )
        
        # Calculate SMARCB1 profiles
        print("Calculating SMARCB1 profiles...")
        basic_smarcb1, replicate_norm_smarcb1 = analyzer.calculate_smarcb1_profile(
            regions,
            os.path.join(smarcb1_dir, "BM3_RPKM.bw"),
            [os.path.join(smarcb1_dir, f"BG{i}_RPKM.bw") for i in range(1, 4)]
        )
        
        # Plot SMARCB1 profiles
        analyzer.plot_smarcb1_profiles(
            basic_smarcb1, replicate_norm_smarcb1,
            category, output_dir
        )
        
        print(f"Completed analysis for {category}")

if __name__ == "__main__":
    main()
