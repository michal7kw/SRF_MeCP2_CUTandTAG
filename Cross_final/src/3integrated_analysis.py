#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pyBigWig
import pysam
from typing import List, Tuple, Dict, Optional
import os
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from Bio import SeqIO
from Bio.Seq import Seq
from dataclasses import dataclass
import logging
from pathlib import Path
import config  # Add this import

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@dataclass
class AnalysisConfig:
    """Configuration for the analysis."""
    window_size: int = 2000
    bin_size: int = 50
    n_processes: int = 16
    categories: List[str] = ('up', 'down', 'no_deg')

class IntegratedAnalyzer:
    def __init__(self, config: AnalysisConfig, genome_fasta: str):
        """
        Initialize the integrated analyzer.
        
        Args:
            config: Analysis configuration
            genome_fasta: Path to genome FASTA file
        """
        self.config = config
        self.n_bins = config.window_size * 2 // config.bin_size
        self.genome = pysam.FastaFile(genome_fasta)
        
    def calculate_gc_content(self, chrom: str, start: int, end: int) -> float:
        """Calculate GC content for a genomic region."""
        try:
            seq = self.genome.fetch(chrom, start, end).upper()
            gc_count = seq.count('G') + seq.count('C')
            total = len(seq)
            return gc_count / total if total > 0 else 0
        except Exception as e:
            logger.warning(f"Error calculating GC content for {chrom}:{start}-{end}: {str(e)}")
            return 0

    def calculate_cpg_density(self, chrom: str, start: int, end: int) -> float:
        """Calculate CpG density (observed/expected ratio) for a genomic region."""
        try:
            seq = self.genome.fetch(chrom, start, end).upper()
            length = len(seq)
            if length < 2:
                return 0.0
            
            # Count CpG, C, and G occurrences
            cpg_count = seq.count('CG')
            c_count = seq.count('C')
            g_count = seq.count('G')
            
            # Calculate observed/expected ratio
            observed = cpg_count
            expected = (c_count * g_count) / length if length > 0 else 0
            return (observed / expected) if expected > 0 else 0
            
        except Exception as e:
            logger.warning(f"Error calculating CpG density for {chrom}:{start}-{end}: {str(e)}")
            return 0

    def calculate_methylation_profile(self,
                                   regions: pd.DataFrame,
                                   ip_file: str,
                                   input_file: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Calculate methylation profiles with different normalizations.
        
        Returns:
            Tuple of (basic_profile, gc_norm_profile, cpg_norm_profile)
        """
        ip_matrix = np.zeros((len(regions), self.n_bins))
        input_matrix = np.zeros((len(regions), self.n_bins))
        gc_content = np.zeros((len(regions), self.n_bins))
        cpg_density = np.zeros((len(regions), self.n_bins))
        
        # Process IP and input signals
        for i, (_, region) in enumerate(tqdm(regions.iterrows(), desc="Processing regions")):
            try:
                center = (region['start'] + region['end']) // 2
                start = max(0, center - self.config.window_size)
                end = center + self.config.window_size
                
                chrom = str(region['seqnames'])
                if not chrom or chrom == 'nan' or start >= end:
                    continue
                
                # Calculate bin coordinates
                bin_edges = np.linspace(start, end, self.n_bins + 1).astype(int)
                
                # Get signals
                with pyBigWig.open(ip_file) as ip_bw, pyBigWig.open(input_file) as input_bw:
                    for j in range(self.n_bins):
                        bin_start, bin_end = bin_edges[j], bin_edges[j+1]
                        
                        # Get signals
                        ip_values = ip_bw.stats(chrom, bin_start, bin_end, type="mean")
                        input_values = input_bw.stats(chrom, bin_start, bin_end, type="mean")
                        
                        ip_matrix[i, j] = ip_values[0] if ip_values and ip_values[0] is not None else 0
                        input_matrix[i, j] = input_values[0] if input_values and input_values[0] is not None else 0
                        
                        # Calculate sequence features
                        gc_content[i, j] = self.calculate_gc_content(chrom, bin_start, bin_end)
                        cpg_density[i, j] = self.calculate_cpg_density(chrom, bin_start, bin_end)
                
            except Exception as e:
                logger.warning(f"Error processing region {chrom}:{start}-{end}: {str(e)}")
                continue
        
        # Calculate profiles
        basic_profile = np.zeros(self.n_bins)
        gc_norm_profile = np.zeros(self.n_bins)
        cpg_norm_profile = np.zeros(self.n_bins)
        
        # Process each bin
        for i in range(self.n_bins):
            # Get valid indices where all required values are present
            valid_idx = (input_matrix[:, i] > 0) & (gc_content[:, i] > 0) & (cpg_density[:, i] > 0)
            
            if np.any(valid_idx):
                # Calculate initial ratios
                initial_ratios = ip_matrix[valid_idx, i] / input_matrix[valid_idx, i]
                
                # Remove outliers from initial ratios
                valid_ratios = initial_ratios[initial_ratios <= np.percentile(initial_ratios, 99)]
                
                if len(valid_ratios) > 0:
                    # Basic profile
                    basic_profile[i] = np.mean(valid_ratios)
                    
                    # Get corresponding GC and CpG values for valid ratios
                    valid_gc = gc_content[valid_idx, i][:len(valid_ratios)]
                    valid_cpg = cpg_density[valid_idx, i][:len(valid_ratios)]
                    
                    # GC normalization
                    gc_ratios = valid_ratios / valid_gc
                    gc_ratios = gc_ratios[~np.isinf(gc_ratios) & ~np.isnan(gc_ratios)]
                    if len(gc_ratios) > 0:
                        gc_norm_profile[i] = np.mean(gc_ratios)
                    
                    # CpG normalization
                    cpg_ratios = valid_ratios / valid_cpg
                    cpg_ratios = cpg_ratios[~np.isinf(cpg_ratios) & ~np.isnan(cpg_ratios)]
                    if len(cpg_ratios) > 0:
                        cpg_norm_profile[i] = np.mean(cpg_ratios)
        
        return basic_profile, gc_norm_profile, cpg_norm_profile

    def calculate_smarcb1_profile(self,
                                regions: pd.DataFrame,
                                bm_file: str,
                                bg_files: List[str]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate SMARCB1 profiles with replicate handling.
        
        Returns:
            Tuple of (basic_profile, replicate_norm_profile)
        """
        bm_matrix = np.zeros((len(regions), self.n_bins))
        bg_matrices = []
        
        # Process BM signal
        logger.info("Processing BM signal...")
        with pyBigWig.open(bm_file) as bm_bw:
            for i, (_, region) in enumerate(tqdm(regions.iterrows())):
                try:
                    center = (region['start'] + region['end']) // 2
                    start = max(0, center - self.config.window_size)
                    end = center + self.config.window_size
                    
                    chrom = str(region['seqnames'])
                    if not chrom or chrom == 'nan' or start >= end:
                        continue
                    
                    values = bm_bw.stats(chrom, start, end, nBins=self.n_bins, type="mean")
                    if values:
                        bm_matrix[i] = [v if v is not None else 0 for v in values]
                except Exception as e:
                    logger.warning(f"Error processing BM region {chrom}:{start}-{end}: {str(e)}")
                    continue
        
        # Process each BG replicate
        logger.info("Processing BG signals...")
        for bg_file in bg_files:
            logger.info(f"Processing {os.path.basename(bg_file)}")
            bg_matrix = np.zeros((len(regions), self.n_bins))
            
            with pyBigWig.open(bg_file) as bg_bw:
                for i, (_, region) in enumerate(tqdm(regions.iterrows())):
                    try:
                        center = (region['start'] + region['end']) // 2
                        start = max(0, center - self.config.window_size)
                        end = center + self.config.window_size
                        
                        chrom = str(region['seqnames'])
                        if not chrom or chrom == 'nan' or start >= end:
                            continue
                        
                        values = bg_bw.stats(chrom, start, end, nBins=self.n_bins, type="mean")
                        if values:
                            bg_matrix[i] = [v if v is not None else 0 for v in values]
                    except Exception as e:
                        logger.warning(f"Error processing BG region {chrom}:{start}-{end}: {str(e)}")
                        continue
            
            bg_matrices.append(bg_matrix)
        
        # Calculate profiles
        basic_profile = np.zeros(self.n_bins)
        replicate_norm_profiles = []
        
        # Calculate normalization for each BG replicate
        for bg_matrix in bg_matrices:
            norm_profile = np.zeros(self.n_bins)
            for i in range(self.n_bins):
                valid_idx = (bg_matrix[:, i] > 0) & (bm_matrix[:, i] > 0)
                if np.any(valid_idx):
                    ratios = bm_matrix[valid_idx, i] / bg_matrix[valid_idx, i]
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
        x = np.linspace(-self.config.window_size, self.config.window_size, self.n_bins)
        
        plt.figure(figsize=(12, 6))
        plt.plot(x, basic_profile, label='Basic (IP/Input)', color='blue')
        plt.plot(x, gc_norm_profile, label='GC-normalized', color='red')
        plt.plot(x, cpg_norm_profile, label='CpG-normalized', color='green')
        plt.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        plt.title(f'Methylation Profiles Around {category} Peaks')
        plt.xlabel('Distance from Peak Center (bp)')
        plt.ylabel('Methylation Enrichment')
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
        x = np.linspace(-self.config.window_size, self.config.window_size, self.n_bins)
        
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

    def calculate_statistics(self, 
                           data: pd.DataFrame,
                           group_col: str,
                           value_col: str) -> pd.DataFrame:
        """Calculate summary statistics for each group."""
        stats_df = data.groupby(group_col)[value_col].agg([
            'count',
            'mean',
            'std',
            'min',
            lambda x: x.quantile(0.25),
            'median',
            lambda x: x.quantile(0.75),
            'max'
        ]).round(3)
        stats_df.columns = ['Count', 'Mean', 'Std', 'Min', 'Q1', 'Median', 'Q3', 'Max']
        return stats_df

    def run_statistical_tests(self,
                            data: pd.DataFrame,
                            group_col: str,
                            value_col: str) -> pd.DataFrame:
        """Run statistical tests between groups."""
        categories = data[group_col].unique()
        results = []
        
        for i, cat1 in enumerate(categories):
            for cat2 in categories[i+1:]:
                group1 = data[data[group_col] == cat1][value_col]
                group2 = data[data[group_col] == cat2][value_col]
                
                # Mann-Whitney U test
                stat, pval = stats.mannwhitneyu(group1, group2, alternative='two-sided')
                
                results.append({
                    'Group1': cat1,
                    'Group2': cat2,
                    'Mann-Whitney_statistic': stat,
                    'p-value': pval,
                    'Significant': pval < 0.05
                })
        
        return pd.DataFrame(results)

def load_and_clean_regions(file_path: str) -> pd.DataFrame:
    """Load and clean region data, ensuring proper coordinate formatting."""
    logger.info(f"Loading regions from {file_path}")
    df = pd.read_csv(file_path)
    
    # Drop rows with NaN values in essential columns
    df = df.dropna(subset=['seqnames', 'start', 'end'])
    
    # Convert coordinates to integers
    df['start'] = df['start'].astype(float).astype(int)
    df['end'] = df['end'].astype(float).astype(int)
    
    # Clean chromosome names
    df['seqnames'] = df['seqnames'].astype(str)
    df = df[~df['seqnames'].str.contains('nan')]
    
    logger.info(f"Loaded {len(df)} valid regions")
    return df

def main():
    # Update paths to use config
    base_dir = config.BASE_DIR
    working_dir = config.BASE_DIR  # Since BASE_DIR already points to Cross_final
    genome_fasta = config.GENOME_FASTA
    medip_dir = config.MEDIP_DIR
    smarcb1_dir = config.SMARCB1_DIR
    output_dir = os.path.join(config.RESULTS_DIR, "3integrated_analysis")
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize config and analyzer
    config = AnalysisConfig()
    analyzer = IntegratedAnalyzer(config, genome_fasta)
    
    # Process each category
    for category in config.categories:
        logger.info(f"\nProcessing {category} regions...")
        
        # Load and clean regions
        regions_file = os.path.join(working_dir, "data", f'extended_{category}.csv')
        if not os.path.exists(regions_file):
            logger.warning(f"File not found - {regions_file}")
            continue
            
        regions = load_and_clean_regions(regions_file)
        if len(regions) == 0:
            logger.warning(f"No valid regions found in {regions_file}")
            continue
        
        # Calculate methylation profiles
        logger.info("Calculating methylation profiles...")
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
        logger.info("Calculating SMARCB1 profiles...")
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
        
        logger.info(f"Completed analysis for {category}")

if __name__ == "__main__":
    main()
