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
        """Initialize with enrichment and RNA-seq data files"""
        self.enrichment_df = pd.read_csv(enrichment_file)
        self.rna_seq_df = pd.read_csv(rna_seq_file)
        
        # Verify RNA-seq data format
        required_cols = ['gene', 'baseMean', 'log2FoldChange', 'padj']
        missing_cols = [col for col in required_cols if col not in self.rna_seq_df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in RNA-seq data: {missing_cols}")
        
        logger.info(f"Loaded {len(self.rna_seq_df)} genes from RNA-seq data")
        
    def categorize_genes(self, log2fc_threshold: float = 1.0, 
                        padj_threshold: float = 0.05) -> pd.DataFrame:
        """
        Categorize genes based on expression changes
        
        Parameters:
        - log2fc_threshold: Minimum absolute log2 fold change to consider (default: 1.0)
        - padj_threshold: Maximum adjusted p-value to consider significant (default: 0.05)
        """
        self.rna_seq_df['category'] = 'non-deregulated'
        
        # Handle NA values in padj
        self.rna_seq_df['padj'] = self.rna_seq_df['padj'].fillna(1.0)
        
        # Up-regulated genes
        self.rna_seq_df.loc[
            (self.rna_seq_df['log2FoldChange'] >= log2fc_threshold) & 
            (self.rna_seq_df['padj'] < padj_threshold),
            'category'
        ] = 'up-regulated'
        
        # Down-regulated genes
        self.rna_seq_df.loc[
            (self.rna_seq_df['log2FoldChange'] <= -log2fc_threshold) & 
            (self.rna_seq_df['padj'] < padj_threshold),
            'category'
        ] = 'down-regulated'
        
        # Log category counts
        category_counts = self.rna_seq_df['category'].value_counts()
        logger.info("\nGene categorization summary:")
        for category, count in category_counts.items():
            logger.info(f"{category}: {count} genes")
        
        return self.rna_seq_df
    
    def integrate_data(self, gene_annotations: pd.DataFrame) -> pd.DataFrame:
        """
        Integrate enrichment data with RNA-seq results
        
        Parameters:
        - gene_annotations: DataFrame with gene coordinates (chr, start, end, gene_name)
        """
        if not {'chr', 'start', 'end', 'gene_name'}.issubset(gene_annotations.columns):
            raise ValueError("Gene annotations must contain: chr, start, end, gene_name")
        
        integrated_data = []
        significant_regions = self.enrichment_df[self.enrichment_df['significant']]
        
        logger.info(f"Processing {len(significant_regions)} significant regions...")
        
        for _, region in significant_regions.iterrows():
            nearby_genes = gene_annotations[
                (gene_annotations['chr'] == region['chr']) &
                ((gene_annotations['start'] - 2000 <= region['end']) &
                 (gene_annotations['end'] + 2000 >= region['start']))
            ]
            
            for _, gene in nearby_genes.iterrows():
                expression_data = self.rna_seq_df[
                    self.rna_seq_df['gene'] == gene['gene_name']
                ]
                
                if not expression_data.empty:
                    expr = expression_data.iloc[0]
                    integrated_data.append({
                        'gene': gene['gene_name'],
                        'chr': region['chr'],
                        'start': region['start'],
                        'end': region['end'],
                        'enrichment': region['enrichment'],
                        'baseMean': expr['baseMean'],
                        'log2FoldChange': expr['log2FoldChange'],
                        'padj': expr['padj'],
                        'category': expr['category']
                    })
        
        return pd.DataFrame(integrated_data)

    def plot_integration_summary(self, integrated_results: pd.DataFrame, 
                               output_dir: str):
        """Create summary plots of the integrated analysis"""
        os.makedirs(output_dir, exist_ok=True)
        
        # Plot 1: Enrichment vs Expression Change
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=integrated_results, 
                       x='log2FoldChange', 
                       y='enrichment',
                       hue='category',
                       alpha=0.6)
        plt.title('MeCP2 Enrichment vs Expression Change')
        plt.xlabel('log2 Fold Change')
        plt.ylabel('MeCP2 Enrichment')
        plt.savefig(os.path.join(output_dir, 'enrichment_vs_expression.pdf'))
        plt.close()
        
        # Plot 2: Category Distribution
        plt.figure(figsize=(8, 6))
        sns.countplot(data=integrated_results, x='category')
        plt.title('Distribution of Gene Categories')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'category_distribution.pdf'))
        plt.close()

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
    parser.add_argument('--log2fc-threshold', type=float, default=1.0,
                       help='log2 fold change threshold (default: 1.0)')
    parser.add_argument('--padj-threshold', type=float, default=0.05,
                       help='Adjusted p-value threshold (default: 0.05)')
    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    try:
        # Initialize integrator
        integrator = RNASeqIntegrator(args.enrichment_file, args.rna_seq_file)

        # Load gene annotations
        gene_annotations = pd.read_csv(args.gene_annotations)

        # Categorize genes
        categorized_genes = integrator.categorize_genes(
            log2fc_threshold=args.log2fc_threshold,
            padj_threshold=args.padj_threshold
        )

        # Integrate data
        integrated_results = integrator.integrate_data(gene_annotations)

        # Save results
        integrated_results.to_csv(
            os.path.join(args.output_dir, 'integrated_results.csv'),
            index=False
        )

        # Create visualizations
        integrator.plot_integration_summary(integrated_results, args.output_dir)

        # Generate summary statistics
        summary_stats = {
            'total_genes': len(integrated_results),
            'up_regulated': len(integrated_results[integrated_results['category'] == 'up-regulated']),
            'down_regulated': len(integrated_results[integrated_results['category'] == 'down-regulated']),
            'non_deregulated': len(integrated_results[integrated_results['category'] == 'non-deregulated']),
            'median_enrichment': integrated_results['enrichment'].median()
        }

        # Save summary statistics
        with open(os.path.join(args.output_dir, 'integration_summary.txt'), 'w') as f:
            for key, value in summary_stats.items():
                f.write(f"{key}: {value}\n")

        logger.info("\nIntegration Summary:")
        logger.info(f"Total genes near enriched regions: {len(integrated_results)}")
        category_counts = integrated_results['category'].value_counts()
        for category, count in category_counts.items():
            logger.info(f"{category}: {count} genes")

    except Exception as e:
        logger.error(f"Error during integration: {str(e)}")
        raise

if __name__ == "__main__":
    main() 