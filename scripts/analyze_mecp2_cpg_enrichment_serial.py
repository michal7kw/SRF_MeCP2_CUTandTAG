import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools
import os
import pysam
from pybedtools import BedTool
import logging
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path
import subprocess

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class EnrichmentConfig:
    """Configuration for enrichment analysis"""
    min_signal_threshold: float = 1.0
    min_fold_change: float = 2.0
    max_qvalue: float = 0.05
    min_peak_width: int = 50

class MeCP2CpGAnalyzer:
    def __init__(self, config: EnrichmentConfig):
        self.config = config
        
    def load_peak_data(self, peak_file: str) -> pd.DataFrame:
        """Load and validate peak file"""
        try:
            peaks = pd.read_csv(peak_file, sep='\t', header=None,
                              names=['chr', 'start', 'end', 'name', 'score',
                                    'strand', 'signalValue', 'pValue',
                                    'qValue', 'peak'])
            
            # Validate required columns
            required_cols = ['chr', 'start', 'end', 'signalValue', 'qValue']
            if not all(col in peaks.columns for col in required_cols):
                raise ValueError(f"Missing required columns in {peak_file}")
            
            # Basic quality filters
            peaks = peaks[
                (peaks['end'] > peaks['start']) &
                (peaks['signalValue'] > 0) &
                (peaks['end'] - peaks['start'] >= self.config.min_peak_width)
            ]
            
            return peaks
            
        except Exception as e:
            logger.error(f"Error loading peak file {peak_file}: {str(e)}")
            raise

    def load_cpg_islands(self, cpg_file: str) -> pd.DataFrame:
        """Load CpG islands data"""
        try:
            cpg = pd.read_csv(cpg_file, sep='\t',
                            names=['chr', 'start', 'end', 'name', 'score', 'strand'])
            logger.info(f"Loaded {len(cpg)} CpG islands")
            return cpg
        except Exception as e:
            logger.error(f"Error loading CpG islands: {str(e)}")
            raise

    def normalize_signal(self, peaks: pd.DataFrame, depth: int) -> pd.DataFrame:
        """Normalize peak signals by sequencing depth"""
        peaks = peaks.copy()
        peaks['signalValue'] = peaks['signalValue'] * (1e6 / depth)  # RPM normalization
        return peaks

    def get_peak_signal(self, peaks_dict: Dict[str, pd.DataFrame], cpg: pd.Series) -> float:
        """
        Calculate average signal across replicates for a CpG region
        
        Parameters:
        - peaks_dict: Dictionary of sample names to peak DataFrames
        - cpg: Series containing CpG island information (chr, start, end)
        
        Returns:
        - Average normalized signal across replicates
        """
        signals = []
        for sample, peaks_df in peaks_dict.items():
            # Find peaks overlapping with CpG island
            overlaps = peaks_df[
                (peaks_df['chr'] == cpg['chr']) &
                (peaks_df['start'] <= cpg['end']) &
                (peaks_df['end'] >= cpg['start'])
            ]
            
            # Sum signal values for overlapping peaks
            signal = overlaps['signalValue'].sum() if not overlaps.empty else 0
            signals.append(signal)
        
        # Return average signal across replicates
        return np.mean(signals) if signals else 0

    def calculate_enrichment(self, 
                           exo_peaks: Dict[str, pd.DataFrame],
                           endo_peaks: Dict[str, pd.DataFrame],
                           cpg_islands: pd.DataFrame,
                           exo_depths: Dict[str, int],
                           endo_depths: Dict[str, int]) -> pd.DataFrame:
        """Calculate enrichment for all CpG islands"""
        # Normalize signals by sequencing depth
        normalized_exo = {
            sample: self.normalize_signal(peaks, exo_depths[sample])
            for sample, peaks in exo_peaks.items()
        }
        
        normalized_endo = {
            sample: self.normalize_signal(peaks, endo_depths[sample])
            for sample, peaks in endo_peaks.items()
        }
        
        enrichment_data = []
        total_cpgs = len(cpg_islands)
        
        logger.info(f"Processing {total_cpgs} CpG islands...")
        
        for idx, cpg in cpg_islands.iterrows():
            if idx % 1000 == 0:
                logger.info(f"Processed {idx}/{total_cpgs} CpG islands...")
            
            exo_signal = self.get_peak_signal(normalized_exo, cpg)
            endo_signal = self.get_peak_signal(normalized_endo, cpg)
            
            enrichment = 0
            if endo_signal > 0:
                enrichment = exo_signal / endo_signal
            elif exo_signal > 0:
                enrichment = float('inf')
            
            enrichment_data.append({
                'chr': cpg['chr'],
                'start': cpg['start'],
                'end': cpg['end'],
                'exo_signal': exo_signal,
                'endo_signal': endo_signal,
                'enrichment': enrichment,
                'significant': (
                    enrichment >= self.config.min_fold_change and 
                    exo_signal > self.config.min_signal_threshold
                )
            })
        
        logger.info("Enrichment calculation complete")
        return pd.DataFrame(enrichment_data)

    def plot_enrichment_results(self, 
                              enrichment_df: pd.DataFrame, 
                              output_dir: str):
        """Create visualization of enrichment results"""
        os.makedirs(output_dir, exist_ok=True)
        
        # Plot 1: Enrichment distribution
        plt.figure(figsize=(10, 6))
        sns.histplot(data=enrichment_df, x='enrichment', bins=50)
        plt.title('Distribution of MeCP2 Enrichment at CpG Islands')
        plt.xlabel('Enrichment (Exo/Endo)')
        plt.ylabel('Count')
        plt.savefig(os.path.join(output_dir, 'enrichment_distribution.pdf'))
        plt.close()
        
        # Plot 2: Signal comparison
        plt.figure(figsize=(10, 6))
        plt.scatter(enrichment_df['endo_signal'], 
                   enrichment_df['exo_signal'], 
                   alpha=0.5)
        plt.xlabel('Endogenous Signal')
        plt.ylabel('Exogenous Signal')
        plt.title('Exogenous vs Endogenous MeCP2 Signal at CpG Islands')
        plt.savefig(os.path.join(output_dir, 'signal_comparison.pdf'))
        plt.close()

def get_sequencing_depths(bam_dir: str) -> Dict[str, int]:
    """
    Calculate sequencing depths from BAM files
    
    Parameters:
    - bam_dir: Directory containing BAM files
    
    Returns:
    - Dictionary mapping sample names to their sequencing depths
    """
    depths = {}
    bam_dir = Path(bam_dir)
    
    # Check if directory exists
    if not bam_dir.exists():
        logger.error(f"BAM directory does not exist: {bam_dir}")
        return depths
        
    # List all BAM files
    bam_files = list(bam_dir.glob('*.bam'))
    if not bam_files:
        logger.error(f"No BAM files found in {bam_dir}")
        return depths
        
    logger.info(f"Found {len(bam_files)} BAM files")
    
    for bam_file in bam_files:
        try:
            # Check if BAM index exists
            if not Path(f"{bam_file}.bai").exists():
                logger.warning(f"BAM index not found for {bam_file}, creating...")
                subprocess.run(['samtools', 'index', str(bam_file)], check=True)
            
            # Get mapped reads using samtools
            cmd = ['samtools', 'view', '-c', '-F', '0x904', str(bam_file)]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Convert output to integer
            depth = int(result.stdout.strip())
            
            # Store both with and without _narrow suffix to match all peak filenames
            sample_name = bam_file.stem
            depths[f"{sample_name}_peaks"] = depth
            depths[f"{sample_name}_narrow_peaks"] = depth
            
            logger.info(f"Calculated depth for {sample_name}: {depth:,} reads")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"samtools error for {bam_file}: {e.stderr}")
        except Exception as e:
            logger.error(f"Error processing {bam_file}: {str(e)}")
            
    return depths

def main():
    # Parse arguments
    import argparse
    parser = argparse.ArgumentParser(description='Analyze MeCP2 enrichment at CpG islands')
    parser.add_argument('--exo-peaks-dir', type=str, required=True,
                       help='Directory containing exogenous peak files')
    parser.add_argument('--endo-peaks-dir', type=str, required=True,
                       help='Directory containing endogenous peak files')
    parser.add_argument('--cpg-islands', type=str, required=True,
                       help='Path to CpG islands file')
    parser.add_argument('--output-dir', type=str, required=True,
                       help='Output directory')
    parser.add_argument('--bam-dir', type=str, required=True,
                       help='Directory containing BAM files')
    args = parser.parse_args()

    # Initialize analyzer
    config = EnrichmentConfig()
    analyzer = MeCP2CpGAnalyzer(config)

    # Debug directories
    logger.info(f"Checking exo peaks directory: {args.exo_peaks_dir}")
    logger.info(f"Directory exists: {os.path.exists(args.exo_peaks_dir)}")
    logger.info(f"Directory contents: {list(Path(args.exo_peaks_dir).glob('*'))}")
    
    logger.info(f"Checking endo peaks directory: {args.endo_peaks_dir}")
    logger.info(f"Directory exists: {os.path.exists(args.endo_peaks_dir)}")
    logger.info(f"Directory contents: {list(Path(args.endo_peaks_dir).glob('*'))}")

    # Load data
    logger.info("Loading peak files...")
    exo_peaks = {}
    endo_peaks = {}
    
    # Load exogenous peaks
    exo_files = list(Path(args.exo_peaks_dir).glob('*.narrowPeak'))
    logger.info(f"Found {len(exo_files)} exo peak files: {[f.name for f in exo_files]}")
    
    for peak_file in exo_files:
        sample_name = peak_file.stem  # e.g., 'NeuV1_peaks'
        try:
            peaks_df = analyzer.load_peak_data(str(peak_file))
            logger.info(f"Successfully loaded {len(peaks_df)} peaks from {peak_file.name}")
            exo_peaks[sample_name] = peaks_df
        except Exception as e:
            logger.error(f"Failed to load {peak_file}: {str(e)}")
    
    # Load endogenous peaks
    endo_files = list(Path(args.endo_peaks_dir).glob('*.narrowPeak'))
    logger.info(f"Found {len(endo_files)} endo peak files: {[f.name for f in endo_files]}")
    
    for peak_file in endo_files:
        sample_name = peak_file.stem
        try:
            peaks_df = analyzer.load_peak_data(str(peak_file))
            logger.info(f"Successfully loaded {len(peaks_df)} peaks from {peak_file.name}")
            endo_peaks[sample_name] = peaks_df
        except Exception as e:
            logger.error(f"Failed to load {peak_file}: {str(e)}")

    # Load CpG islands
    logger.info(f"Loading CpG islands from: {args.cpg_islands}")
    logger.info(f"CpG file exists: {os.path.exists(args.cpg_islands)}")
    cpg_islands = analyzer.load_cpg_islands(args.cpg_islands)

    # Get sequencing depths with better error handling
    logger.info(f"Checking BAM directory: {args.bam_dir}")
    logger.info(f"BAM directory exists: {os.path.exists(args.bam_dir)}")
    if os.path.exists(args.bam_dir):
        logger.info(f"BAM directory contents: {list(Path(args.bam_dir).glob('*.bam'))}")
    
    all_depths = get_sequencing_depths(args.bam_dir)
    
    if not all_depths:
        logger.error("No sequencing depths could be calculated. Check BAM files and permissions.")
        logger.error("Expected BAM files:")
        for peak_name in exo_peaks.keys():
            bam_name = peak_name.replace('_peaks', '')
            logger.error(f"  {bam_name}.bam")
        for peak_name in endo_peaks.keys():
            bam_name = peak_name.replace('_peaks', '')
            logger.error(f"  {bam_name}.bam")
        raise ValueError("Failed to calculate sequencing depths")
    
    # Debug output
    logger.info("\nPeak files found:")
    logger.info(f"Exo peaks: {list(exo_peaks.keys())}")
    logger.info(f"Endo peaks: {list(endo_peaks.keys())}")
    logger.info("\nBAM depths found:")
    logger.info(f"All depths: {list(all_depths.keys())}")
    
    # Split depths into exo and endo, matching peak file names
    exo_depths = {
        sample: depth 
        for sample, depth in all_depths.items() 
        if any(sample.startswith(prefix) for prefix in ['NSCv', 'NeuV'])
    }
    
    endo_depths = {
        sample: depth 
        for sample, depth in all_depths.items() 
        if any(sample.startswith(prefix) for prefix in ['NSCM', 'NeuM'])
    }
    
    # Verify all peak files have corresponding depth values
    missing_exo = set(exo_peaks.keys()) - set(exo_depths.keys())
    missing_endo = set(endo_peaks.keys()) - set(endo_depths.keys())
    
    if missing_exo or missing_endo:
        logger.error("Missing depth values for some peak files:")
        if missing_exo:
            logger.error(f"Exo samples: {missing_exo}")
        if missing_endo:
            logger.error(f"Endo samples: {missing_endo}")
        raise ValueError("Missing depth values for some samples")
    
    logger.info("\nExogenous sample depths:")
    for sample, depth in exo_depths.items():
        logger.info(f"{sample}: {depth:,} reads")
        
    logger.info("\nEndogenous sample depths:")
    for sample, depth in endo_depths.items():
        logger.info(f"{sample}: {depth:,} reads")
    
    # Add to main() before enrichment calculation
    if not all_depths:
        raise ValueError(f"No BAM files found in {args.bam_dir}")

    logger.info("\nFile name mapping:")
    for peak_name in exo_peaks.keys():
        bam_name = peak_name.replace('_peaks', '')
        logger.info(f"Peak file: {peak_name} -> BAM file: {bam_name}.bam")
    
    # Calculate enrichment with depth normalization
    enrichment_results = analyzer.calculate_enrichment(
        exo_peaks, 
        endo_peaks, 
        cpg_islands,
        exo_depths,
        endo_depths
    )

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Save results
    enrichment_results.to_csv(
        os.path.join(args.output_dir, 'mecp2_cpg_enrichment.csv'),
        index=False
    )

    # Create visualizations
    logger.info("Creating visualizations...")
    analyzer.plot_enrichment_results(enrichment_results, args.output_dir)

    # Print summary
    logger.info("\nAnalysis Summary:")
    logger.info(f"Total CpG islands analyzed: {len(enrichment_results)}")
    logger.info(f"Significantly enriched regions: {enrichment_results['significant'].sum()}")
    logger.info(f"Mean enrichment: {enrichment_results['enrichment'].mean():.2f}")

if __name__ == "__main__":
    main() 