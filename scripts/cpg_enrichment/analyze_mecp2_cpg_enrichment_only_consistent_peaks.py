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

    def check_peak_consistency(self, overlaps_dict: Dict[str, pd.DataFrame]) -> bool:
        """Check if a region has consistent peaks across replicates"""
        replicates_with_peaks = sum(1 for df in overlaps_dict.values() if not df.empty)
        
        # For M2 endo group (2 replicates), require both to have peaks
        if len(overlaps_dict) == 2:
            return replicates_with_peaks == 2
        
        # For groups with 3 replicates, require at least 2
        return replicates_with_peaks >= 2

    def calculate_enrichment_chunk(self, 
                                   chunk_id: int,
                                   total_chunks: int,
                                   exo_peaks: Dict[str, pd.DataFrame],
                                   endo_peaks: Dict[str, pd.DataFrame],
                                   cpg_islands: pd.DataFrame,
                                   exo_depths: Dict[str, int],
                                   endo_depths: Dict[str, int]) -> pd.DataFrame:
        """Calculate enrichment for a subset of CpG islands"""
        # Calculate chunk boundaries
        chunk_size = len(cpg_islands) // total_chunks
        start_idx = chunk_id * chunk_size
        end_idx = start_idx + chunk_size if chunk_id < total_chunks - 1 else len(cpg_islands)
        
        cpg_chunk = cpg_islands.iloc[start_idx:end_idx]
        
        # Normalize signals
        normalized_exo = {
            sample: self.normalize_signal(peaks, exo_depths[sample])
            for sample, peaks in exo_peaks.items()
        }
        
        normalized_endo = {
            sample: self.normalize_signal(peaks, endo_depths[sample])
            for sample, peaks in endo_peaks.items()
        }
        
        enrichment_data = []
        for _, cpg in cpg_chunk.iterrows():
            # Get overlapping peaks
            exo_overlaps = {
                sample: peaks[
                    (peaks['chr'] == cpg['chr']) &
                    (peaks['start'] <= cpg['end']) &
                    (peaks['end'] >= cpg['start'])
                ]
                for sample, peaks in normalized_exo.items()
            }
            
            endo_overlaps = {
                sample: peaks[
                    (peaks['chr'] == cpg['chr']) &
                    (peaks['start'] <= cpg['end']) &
                    (peaks['end'] >= cpg['start'])
                ]
                for sample, peaks in normalized_endo.items()
            }
            
            # Check consistency of binding
            has_consistent_exo = self.check_peak_consistency(exo_overlaps)
            has_consistent_endo = self.check_peak_consistency(endo_overlaps)
            
            # Skip regions with no consistent binding
            if not (has_consistent_exo or has_consistent_endo):
                continue
            
            # Calculate signals only for consistent peaks
            def calculate_weighted_signal(overlaps_dict):
                signals = []
                for df in overlaps_dict.values():
                    if df.empty:
                        signals.append(0.0)
                    else:
                        weighted_signal = np.sum(
                            df['signalValue'] * (df['end'] - df['start'])
                        )
                        signals.append(weighted_signal)
                return np.mean(signals) if signals else 0.0
            
            exo_signal = calculate_weighted_signal(exo_overlaps)
            endo_signal = calculate_weighted_signal(endo_overlaps)
            
            # Determine binding type based on consistent peaks
            if has_consistent_exo and has_consistent_endo:
                binding_type = 'both'
            elif has_consistent_exo:
                binding_type = 'exo_only'
            elif has_consistent_endo:
                binding_type = 'endo_only'
            else:
                continue  # Skip inconsistent regions
            
            # Collect all width-weighted signals for statistical testing
            all_exo_signals = []
            all_endo_signals = []
            
            for df in exo_overlaps.values():
                if not df.empty:
                    weighted_signals = df['signalValue'] * (df['end'] - df['start'])
                    all_exo_signals.extend(weighted_signals)
            
            for df in endo_overlaps.values():
                if not df.empty:
                    weighted_signals = df['signalValue'] * (df['end'] - df['start'])
                    all_endo_signals.extend(weighted_signals)
            
            if all_exo_signals and all_endo_signals:
                _, pvalue = stats.mannwhitneyu(all_exo_signals, all_endo_signals)
            else:
                pvalue = 1.0
            
            # Calculate enrichment with safety check
            if endo_signal > 0:
                enrichment = exo_signal / endo_signal
            else:
                enrichment = float('inf') if exo_signal > 0 else 0.0
            
            # Update significance check to avoid division by zero
            is_significant = (
                (enrichment >= self.config.min_fold_change if endo_signal > 0 else exo_signal > 0) and
                exo_signal > self.config.min_signal_threshold and
                pvalue < 0.05
            )
            
            enrichment_data.append({
                'chr': cpg['chr'],
                'start': cpg['start'],
                'end': cpg['end'],
                'exo_signal': exo_signal,
                'endo_signal': endo_signal,
                'enrichment': enrichment,
                'pvalue': pvalue,
                'binding_type': binding_type,
                'peak_width_exo': np.mean([
                    df['end'].mean() - df['start'].mean() 
                    for df in exo_overlaps.values() if not df.empty
                ]) if has_consistent_exo else 0,
                'peak_width_endo': np.mean([
                    df['end'].mean() - df['start'].mean() 
                    for df in endo_overlaps.values() if not df.empty
                ]) if has_consistent_endo else 0,
                'significant': is_significant
            })
        
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
    Calculate sequencing depths from BAM files and map to all possible peak file names
    """
    depths = {}
    bam_dir = Path(bam_dir)
    
    if not bam_dir.exists():
        logger.error(f"BAM directory does not exist: {bam_dir}")
        return depths
        
    bam_files = list(bam_dir.glob('*.bam'))
    if not bam_files:
        logger.error(f"No BAM files found in {bam_dir}")
        return depths
        
    logger.info(f"Found {len(bam_files)} BAM files")
    
    for bam_file in bam_files:
        try:
            if not Path(f"{bam_file}.bai").exists():
                logger.warning(f"BAM index not found for {bam_file}, creating...")
                subprocess.run(['samtools', 'index', str(bam_file)], check=True)
            
            cmd = ['samtools', 'view', '-c', '-F', '0x904', str(bam_file)]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            depth = int(result.stdout.strip())
            
            # Map BAM file depth to all possible peak file names
            base_name = bam_file.stem
            depths[f"{base_name}_peaks"] = depth
            depths[f"{base_name}_narrow_peaks"] = depth
            
            logger.info(f"Calculated depth for {base_name}: {depth:,} reads")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"samtools error for {bam_file}: {e.stderr}")
        except Exception as e:
            logger.error(f"Error processing {bam_file}: {str(e)}")
            
    return depths

def load_peak_files(peak_dir: Path, pattern: str = "*_peaks.narrowPeak") -> Dict[str, BedTool]:
    """Load peak files from directory"""
    logger.info(f"Checking peaks directory: {peak_dir}")
    logger.info(f"Directory exists: {peak_dir.exists()}")
    
    if not peak_dir.exists():
        raise ValueError(f"Peak directory does not exist: {peak_dir}")
        
    # Use glob to find peak files
    peak_files = list(peak_dir.glob(pattern))
    logger.info(f"Directory contents: {peak_files}")
    
    # Create dictionary of BedTool objects
    peaks = {}
    for peak_file in peak_files:
        try:
            # Extract sample name from filename
            sample_name = peak_file.stem  # Remove .narrowPeak extension
            peaks[sample_name] = BedTool(str(peak_file))
            logger.info(f"Loaded peak file: {peak_file}")
        except Exception as e:
            logger.error(f"Error loading {peak_file}: {e}")
            
    logger.info(f"Found {len(peaks)} peak files: {list(peaks.keys())}")
    return peaks

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
    parser.add_argument('--chunk-id', type=int, required=True,
                       help='Chunk ID for parallel processing')
    parser.add_argument('--total-chunks', type=int, required=True,
                       help='Total number of chunks')
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

    # Load peak files with more detailed logging
    logger.info("Loading peak files...")
    exo_peaks = load_peak_files(Path(args.exo_peaks_dir))
    endo_peaks = load_peak_files(Path(args.endo_peaks_dir))
    
    if not exo_peaks or not endo_peaks:
        raise ValueError("No peak files found in one or both directories")
        
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
    
    # Calculate enrichment for this chunk
    chunk_results = analyzer.calculate_enrichment_chunk(
        args.chunk_id,
        args.total_chunks,
        exo_peaks,
        endo_peaks,
        cpg_islands,
        exo_depths,
        endo_depths
    )
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Save chunk results
    chunk_file = os.path.join(args.output_dir, f'chunk_{args.chunk_id}.csv')
    chunk_results.to_csv(chunk_file, index=False)
    
    logger.info(f"Chunk {args.chunk_id} complete. Results saved to {chunk_file}")

if __name__ == "__main__":
    main() 