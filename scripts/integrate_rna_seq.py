import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import logging
from typing import Dict, List, Tuple
from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class RNASeqIntegrator:
    def __init__(self, enrichment_file: str, rna_seq_file: str):
        self.enrichment_df = pd.read_csv(enrichment_file)
        self.rna_seq_df = pd.read_csv(rna_seq_file)
        
    def categorize_genes(self) -> pd.DataFrame:
        """Categorize genes based on expression changes"""
        self.rna_seq_df['category'] = 'non-deregulated'
        self.rna_seq_df.loc[
            (self.rna_seq_df['log2FoldChange'] > 1) & 
            (self.rna_seq_df['padj'] < 0.05),
            'category'
        ] = 'up-regulated'
        self.rna_seq_df.loc[
            (self.rna_seq_df['log2FoldChange'] < -1) & 
            (self.rna_seq_df['padj'] < 0.05),
            'category'
        ] = 'down-regulated'
        
        return self.rna_seq_df
    
    def integrate_data(self, gene_annotations: pd.DataFrame) -> pd.DataFrame:
        """Integrate enrichment data with RNA-seq results"""
        # Find genes near enriched regions
        integrated_data = []
        
        for _, region in self.enrichment_df[self.enrichment_df['significant']].iterrows():
            nearby_genes = gene_annotations[
                (gene_annotations['chr'] == region['chr']) &
                ((gene_annotations['start'] - 2000 <= region['end']) &
                 (gene_annotations['end'] + 2000 >= region['start']))
            ]
            
            for _, gene in nearby_genes.iterrows():
                if gene['gene_name'] in self.rna_seq_df['gene'].values:
                    expression_data = self.rna_seq_df[
                        self.rna_seq_df['gene'] == gene['gene_name']
                    ].iloc[0]
                    
                    integrated_data.append({
                        'gene': gene['gene_name'],
                        'chr': region['chr'],
                        'start': region['start'],
                        'end': region['end'],
                        'enrichment': region['enrichment'],
                        'log2FoldChange': expression_data['log2FoldChange'],
                        'padj': expression_data['padj'],
                        'category': expression_data['category']
                    })
        
        return pd.DataFrame(integrated_data)

def main():
    # Parse arguments
    import argparse
    parser = argparse.ArgumentParser(description='Integrate MeCP2 enrichment with RNA-seq data')
    parser.add_argument('--enrichment-file', type=str, required=True,
                       help='Path to enrichment results file')
    parser.add_argument('--rna-seq-file', type=str, required=True,
                       help='Path to RNA-seq results file')
    parser.add_argument('--gene-annotations', type=str, required=True,
                       help='Path to gene annotations file')
    parser.add_argument('--output-dir', type=str, required=True,
                       help='Output directory')
    args = parser.parse_args()

    # Initialize integrator
    integrator = RNASeqIntegrator(args.enrichment_file, args.rna_seq_file)

    # Load gene annotations
    gene_annotations = pd.read_csv(args.gene_annotations)

    # Categorize genes
    categorized_genes = integrator.categorize_genes()

    # Integrate data
    integrated_results = integrator.integrate_data(gene_annotations)

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Save results
    integrated_results.to_csv(
        os.path.join(args.output_dir, 'integrated_results.csv'),
        index=False
    )

    # Print summary
    logger.info("\nIntegration Summary:")
    logger.info(f"Total genes near enriched regions: {len(integrated_results)}")
    category_counts = integrated_results['category'].value_counts()
    for category, count in category_counts.items():
        logger.info(f"{category}: {count} genes")

if __name__ == "__main__":
    main() 