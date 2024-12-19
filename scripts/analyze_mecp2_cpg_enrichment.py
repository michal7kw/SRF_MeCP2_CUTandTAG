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
        peaks['signalValue'] = peaks['signalValue'] * (1e6 / depth)
        return peaks

    def calculate_enrichment(self, 
                           exo_peaks: Dict[str, pd.DataFrame],
                           endo_peaks: Dict[str, pd.DataFrame],
                           cpg_islands: pd.DataFrame) -> pd.DataFrame:
        """Calculate MeCP2 enrichment at CpG islands"""
        
        # Convert CpG islands to BedTool
        cpg_bed = BedTool.from_dataframe(cpg_islands[['chr', 'start', 'end']])
        
        # Process and combine replicates
        def combine_replicate_peaks(peaks_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
            combined = []
            for sample, peaks in peaks_dict.items():
                if not peaks.empty:
                    peaks_bed = BedTool.from_dataframe(peaks[['chr', 'start', 'end', 'signalValue']])
                    intersect = peaks_bed.intersect(cpg_bed, wa=True, wb=True)
                    
                    # Convert intersection results to DataFrame
                    if str(intersect):
                        df = pd.DataFrame([
                            str(interval).strip().split('\t') 
                            for interval in str(intersect).strip().split('\n')
                        ])
                        df.columns = ['peak_chr', 'peak_start', 'peak_end', 'signalValue',
                                    'cpg_chr', 'cpg_start', 'cpg_end']
                        df['sample'] = sample
                        combined.append(df)
            
            return pd.concat(combined) if combined else pd.DataFrame()

        # Process exo and endo peaks
        exo_combined = combine_replicate_peaks(exo_peaks)
        endo_combined = combine_replicate_peaks(endo_peaks)

        # Calculate enrichment for each CpG island
        enrichment_data = []
        for _, cpg in cpg_islands.iterrows():
            # Get overlapping peaks
            exo_overlaps = exo_combined[
                (exo_combined['cpg_chr'] == cpg['chr']) &
                (exo_combined['cpg_start'].astype(int) == cpg['start']) &
                (exo_combined['cpg_end'].astype(int) == cpg['end'])
            ]
            
            endo_overlaps = endo_combined[
                (endo_combined['cpg_chr'] == cpg['chr']) &
                (endo_combined['cpg_start'].astype(int) == cpg['start']) &
                (endo_combined['cpg_end'].astype(int) == cpg['end'])
            ]
            
            # Calculate signals
            exo_signal = exo_overlaps['signalValue'].astype(float).sum()
            endo_signal = endo_overlaps['signalValue'].astype(float).sum()
            
            # Calculate enrichment
            enrichment = exo_signal / max(endo_signal, self.config.min_signal_threshold)
            
            enrichment_data.append({
                'chr': cpg['chr'],
                'start': cpg['start'],
                'end': cpg['end'],
                'exo_signal': exo_signal,
                'endo_signal': endo_signal,
                'enrichment': enrichment,
                'n_exo_peaks': len(exo_overlaps),
                'n_endo_peaks': len(endo_overlaps)
            })
        
        enrichment_df = pd.DataFrame(enrichment_data)
        
        # Add significance flags
        enrichment_df['significant'] = (
            (enrichment_df['enrichment'] >= self.config.min_fold_change) &
            (enrichment_df['n_exo_peaks'] > 0)
        )
        
        return enrichment_df

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
    args = parser.parse_args()

    # Initialize analyzer
    config = EnrichmentConfig()
    analyzer = MeCP2CpGAnalyzer(config)

    # Load data
    logger.info("Loading peak files...")
    exo_peaks = {}
    endo_peaks = {}
    
    # Load exogenous peaks
    for peak_file in Path(args.exo_peaks_dir).glob('*.narrowPeak'):
        sample_name = peak_file.stem
        exo_peaks[sample_name] = analyzer.load_peak_data(str(peak_file))
    
    # Load endogenous peaks
    for peak_file in Path(args.endo_peaks_dir).glob('*.narrowPeak'):
        sample_name = peak_file.stem
        endo_peaks[sample_name] = analyzer.load_peak_data(str(peak_file))

    # Load CpG islands
    cpg_islands = analyzer.load_cpg_islands(args.cpg_islands)

    # Calculate enrichment
    logger.info("Calculating enrichment...")
    enrichment_results = analyzer.calculate_enrichment(exo_peaks, endo_peaks, cpg_islands)

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