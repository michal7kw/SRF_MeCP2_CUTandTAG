#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import pyBigWig
import os
import logging
import argparse
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path
from glob import glob

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

CPG_WINDOW = 0
NUM_REPS = 2

@dataclass
class EnrichmentConfig:
    """Configuration for enrichment analysis"""
    min_signal_threshold: float = 0.0
    min_fold_change: float = 1.0
    max_qvalue: float = 0.05

class CpGAnalyzer:
    def __init__(self, config: EnrichmentConfig):
        self.config = config
        
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

    def get_signal_from_bigwig(
        self,
        bigwig_files: Dict[str, str],
        chrom: str,
        start: int,
        end: int,
        method: str = "mean"
    ) -> Tuple[float, List[float]]:
        """Calculate signal from multiple bigWig replicates"""
        signals = []
        
        for sample, bw_path in bigwig_files.items():
            with pyBigWig.open(bw_path) as bw:
                try:
                    signal = bw.stats(
                        chrom,
                        start,
                        end,
                        type=method
                    )[0]
                    signals.append(signal if signal is not None else 0.0)
                except RuntimeError:
                    signals.append(0.0)
        
        return np.mean(signals) if signals else 0.0, signals

    def load_peak_files(self, peak_dir: str) -> List[pd.DataFrame]:
        """Load individual peak files from directory"""
        peak_files = glob(f"{peak_dir}/*_peaks.broadPeak")
        peaks = []
        for peak_file in peak_files:
            df = pd.read_csv(peak_file, sep='\t',
                            names=['chr', 'start', 'end', 'name', 'score',
                                  'strand', 'signalValue', 'pValue',
                                  'qValue', 'peak'])
            peaks.append(df)
        return peaks

    def calculate_enrichment_chunk(self, 
                                 chunk_id: int,
                                 total_chunks: int,
                                 neu_bigwigs: Dict[str, str],
                                 nsc_bigwigs: Dict[str, str],
                                 cpg_islands: pd.DataFrame,
                                 neu_peaks_list: List[pd.DataFrame],
                                 nsc_peaks_list: List[pd.DataFrame]) -> pd.DataFrame:
        """Calculate enrichment comparing Neuron vs NSC exo samples using peak-defined regions"""
        
        chunk_size = len(cpg_islands) // total_chunks
        start_idx = chunk_id * chunk_size
        end_idx = start_idx + chunk_size if chunk_id < total_chunks - 1 else len(cpg_islands)
        
        cpg_chunk = cpg_islands.iloc[start_idx:end_idx]
        
        enrichment_data = []
        for _, cpg in cpg_chunk.iterrows():
            # Find overlapping peaks in each replicate
            neu_overlaps_by_rep = []
            for neu_peaks in neu_peaks_list:
                overlaps = neu_peaks[
                    (neu_peaks['chr'] == cpg['chr']) &
                    (neu_peaks['start'] <= cpg['end'] + CPG_WINDOW) &
                    (neu_peaks['end'] >= cpg['start'] - CPG_WINDOW)
                ]
                neu_overlaps_by_rep.append(overlaps)
                
            nsc_overlaps_by_rep = []
            for nsc_peaks in nsc_peaks_list:
                overlaps = nsc_peaks[
                    (nsc_peaks['chr'] == cpg['chr']) &
                    (nsc_peaks['start'] <= cpg['end'] + CPG_WINDOW) &
                    (nsc_peaks['end'] >= cpg['start'] - CPG_WINDOW)
                ]
                nsc_overlaps_by_rep.append(overlaps)
            
            # Check if we have peaks in at least NUM_REPS replicates
            neu_rep_with_peaks = sum(1 for overlaps in neu_overlaps_by_rep if not overlaps.empty)
            nsc_rep_with_peaks = sum(1 for overlaps in nsc_overlaps_by_rep if not overlaps.empty)
            
            # Skip if neither condition has at least NUM_REPS replicates with peaks
            if neu_rep_with_peaks < NUM_REPS and nsc_rep_with_peaks < NUM_REPS:
                continue
            
            # Get the outer bounds of all overlapping peaks
            all_starts = []
            all_ends = []
            
            for overlaps in neu_overlaps_by_rep + nsc_overlaps_by_rep:
                if not overlaps.empty:
                    all_starts.extend(overlaps['start'].tolist())
                    all_ends.extend(overlaps['end'].tolist())
            
            # Use the outermost boundaries of overlapping peaks
            region_start = min(all_starts)
            region_end = max(all_ends)
            
            # Get signals using the peak-defined boundaries
            neu_signal, neu_replicates = self.get_signal_from_bigwig(
                neu_bigwigs,
                cpg['chr'],
                region_start,
                region_end
            )
            
            nsc_signal, nsc_replicates = self.get_signal_from_bigwig(
                nsc_bigwigs,
                cpg['chr'],
                region_start,
                region_end
            )
            
            # Check for signal in at least NUM_REPS replicate
            neu_has_signal = sum(1 for s in neu_replicates if s > 0) >= NUM_REPS
            nsc_has_signal = sum(1 for s in nsc_replicates if s > 0) >= NUM_REPS
            
            if not (neu_has_signal or nsc_has_signal):
                continue
            
            # Statistical test between replicates
            if neu_replicates and nsc_replicates:
                nonzero_neu = [s for s in neu_replicates if s > 0]
                nonzero_nsc = [s for s in nsc_replicates if s > 0]
                
                if len(nonzero_neu) >= 1 and len(nonzero_nsc) >= 1:
                    _, pvalue = stats.mannwhitneyu(
                        nonzero_neu,
                        nonzero_nsc,
                        alternative='two-sided'
                    )
                else:
                    pvalue = 1.0
            else:
                pvalue = 1.0
            
            # Calculate fold change (Neuron/NSC)
            if nsc_signal > 0:
                fold_change = neu_signal / nsc_signal
            else:
                fold_change = float('inf') if neu_signal > 0 else 0.0
            
            # Determine binding type
            if neu_has_signal and nsc_has_signal:
                binding_type = 'both'
            elif neu_has_signal:
                binding_type = 'neu_only'
            elif nsc_has_signal:
                binding_type = 'nsc_only'
            else:
                continue
            
            # Determine binding type by peaks
            if neu_rep_with_peaks >= NUM_REPS and nsc_rep_with_peaks >= NUM_REPS:
                binding_type_by_peaks = 'both'
            elif neu_rep_with_peaks >= NUM_REPS:
                binding_type_by_peaks = 'neu_only'
            elif nsc_rep_with_peaks >= NUM_REPS:
                binding_type_by_peaks = 'nsc_only'
            else:
                binding_type_by_peaks = 'none'
            
            # Store results with replicate information
            enrichment_data.append({
                'chr': cpg['chr'],
                'start': cpg['start'],
                'end': cpg['end'],
                'neu_signal': neu_signal,
                'nsc_signal': nsc_signal,
                'fold_change': fold_change,
                'pvalue': pvalue,
                'binding_type': binding_type,
                'binding_type_by_peaks': binding_type_by_peaks,
                'significant': (
                    (fold_change >= self.config.min_fold_change if nsc_signal > 0 else neu_signal > 0) and
                    neu_signal > self.config.min_signal_threshold and
                    pvalue < self.config.max_qvalue
                ),
                'neu_replicates_with_signal': sum(1 for s in neu_replicates if s > 0),
                'nsc_replicates_with_signal': sum(1 for s in nsc_replicates if s > 0),
                'neu_replicate_signals': ','.join(map(str, neu_replicates)),
                'nsc_replicate_signals': ','.join(map(str, nsc_replicates)),
                'region_length': region_end - region_start,
                'cpg_length': cpg['end'] - cpg['start'],
                'cpg_score': cpg['score'],
                'cpg_name': cpg['name'],
                'neu_replicates_with_peaks': neu_rep_with_peaks,
                'nsc_replicates_with_peaks': nsc_rep_with_peaks,
                'neu_peak_scores_by_rep': ';'.join(','.join(map(str, rep['signalValue'].tolist())) 
                                                  for rep in neu_overlaps_by_rep if not rep.empty),
                'nsc_peak_scores_by_rep': ';'.join(','.join(map(str, rep['signalValue'].tolist()))
                                                  for rep in nsc_overlaps_by_rep if not rep.empty),
                'region_start': region_start,
                'region_end': region_end
            })
        
        return pd.DataFrame(enrichment_data)

def main():
    parser = argparse.ArgumentParser(description='Compare CpG enrichment between Neuron and NSC exo samples')
    parser.add_argument('--bigwig-dir', required=True, help='Directory containing bigWig files')
    parser.add_argument('--cpg-islands', required=True, help='Path to CpG islands bed file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--chunk-id', type=int, required=True, help='Chunk ID for parallel processing')
    parser.add_argument('--total-chunks', type=int, required=True, help='Total number of chunks')
    parser.add_argument('--neu-peaks', required=True, help='Directory containing Neuron peak files')
    parser.add_argument('--nsc-peaks', required=True, help='Directory containing NSC peak files')
    parser.add_argument('--min-signal', type=float, default=0.1, help='Minimum signal threshold')
    parser.add_argument('--min-fold-change', type=float, default=1.5, help='Minimum fold change threshold')
    parser.add_argument('--max-qvalue', type=float, default=0.05, help='Maximum q-value threshold')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Initialize analyzer with custom thresholds
    config = EnrichmentConfig(
        min_signal_threshold=args.min_signal,
        min_fold_change=args.min_fold_change,
        max_qvalue=args.max_qvalue
    )
    analyzer = CpGAnalyzer(config)
    
    # Load CpG islands
    cpg_islands = analyzer.load_cpg_islands(args.cpg_islands)
    
    # Define sample mappings
    neu_samples = ["NeuV1", "NeuV2", "NeuV3"]
    nsc_samples = ["NSCv1", "NSCv2", "NSCv3"]
    
    # Create bigWig file mappings
    neu_bigwigs = {
        sample: os.path.join(args.bigwig_dir, f"{sample}.bw")
        for sample in neu_samples
    }
    nsc_bigwigs = {
        sample: os.path.join(args.bigwig_dir, f"{sample}.bw")
        for sample in nsc_samples
    }
    
    # Load peak files
    neu_peaks = analyzer.load_peak_files(args.neu_peaks)
    nsc_peaks = analyzer.load_peak_files(args.nsc_peaks)
    
    # Calculate enrichment for this chunk
    results = analyzer.calculate_enrichment_chunk(
        args.chunk_id,
        args.total_chunks,
        neu_bigwigs,
        nsc_bigwigs,
        cpg_islands,
        neu_peaks,
        nsc_peaks
    )
    
    # Save results
    output_file = os.path.join(args.output_dir, f"enrichment_chunk_{args.chunk_id}.csv")
    results.to_csv(output_file, index=False)
    
    logger.info(f"Completed chunk {args.chunk_id}. Found {len(results)} enriched regions.")

if __name__ == "__main__":
    main()
