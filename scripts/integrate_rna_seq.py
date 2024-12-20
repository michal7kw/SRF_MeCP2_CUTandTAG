import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import logging
from typing import Dict, List, Tuple
from pathlib import Path

# Set up more detailed logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class RNASeqIntegrator:
    def __init__(self, enrichment_file: str, rna_seq_file: str):
        """Initialize with enrichment and RNA-seq data files"""
        logger.info(f"Loading enrichment data from: {enrichment_file}")
        self.enrichment_df = pd.read_csv(enrichment_file)
        logger.info(f"Enrichment data shape: {self.enrichment_df.shape}")
        logger.info(f"Enrichment columns: {self.enrichment_df.columns.tolist()}")
        
        logger.info(f"Loading RNA-seq data from: {rna_seq_file}")
        self.rna_seq_df = pd.read_csv(rna_seq_file)
        logger.info(f"RNA-seq data shape: {self.rna_seq_df.shape}")
        logger.info(f"RNA-seq columns: {self.rna_seq_df.columns.tolist()}")
        
        # Filter for significant MeCP2 enrichment
        self.significant_regions = self.enrichment_df[self.enrichment_df['significant']]
        logger.info(f"Found {len(self.significant_regions)} significantly enriched CpG regions")
        logger.info("Sample of significant regions:")
        logger.info(self.significant_regions.head().to_string())
        
        # Verify RNA-seq data format
        required_cols = ['gene', 'baseMean', 'log2FoldChange', 'padj']
        missing_cols = [col for col in required_cols if col not in self.rna_seq_df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in RNA-seq data: {missing_cols}")
        
    def integrate_data(self, gene_annotations: pd.DataFrame) -> pd.DataFrame:
        """Map significantly enriched CpG regions to nearby genes and integrate with RNA-seq data"""
        logger.info(f"Gene annotations shape: {gene_annotations.shape}")
        logger.info(f"Gene annotations columns: {gene_annotations.columns.tolist()}")
        logger.info("Sample of gene annotations:")
        logger.info(gene_annotations.head().to_string())
        
        if not {'chr', 'start', 'end', 'gene_name'}.issubset(gene_annotations.columns):
            raise ValueError("Gene annotations must contain: chr, start, end, gene_name")
        
        integrated_data = []
        total_regions = len(self.significant_regions)
        
        logger.info(f"Processing {total_regions} significant regions...")
        for idx, region in self.significant_regions.iterrows():
            if idx % 100 == 0:
                logger.info(f"Processing region {idx}/{total_regions}")
            
            nearby_genes = gene_annotations[
                (gene_annotations['chr'] == region['chr']) &
                ((gene_annotations['start'] - 2000 <= region['end']) &
                 (gene_annotations['end'] + 2000 >= region['start']))
            ]
            
            if len(nearby_genes) > 0:
                logger.debug(f"Found {len(nearby_genes)} genes near region {region['chr']}:{region['start']}-{region['end']}")
            
            for _, gene in nearby_genes.iterrows():
                gene_data = {
                    'gene': gene['gene_name'],
                    'chr': region['chr'],
                    'cpg_start': region['start'],
                    'cpg_end': region['end'],
                    'mecp2_enrichment': region['enrichment'],
                    'exo_signal': region['exo_signal'],
                    'endo_signal': region['endo_signal'],
                    'distance_to_gene': min(
                        abs(region['start'] - gene['start']),
                        abs(region['end'] - gene['end'])
                    )
                }
                
                # Add RNA-seq data if available
                rna_data = self.rna_seq_df[self.rna_seq_df['gene'] == gene['gene_name']]
                if not rna_data.empty:
                    expr = rna_data.iloc[0]
                    gene_data.update({
                        'baseMean': expr['baseMean'],
                        'log2FoldChange': expr['log2FoldChange'],
                        'padj': expr['padj']
                    })
                else:
                    gene_data.update({
                        'baseMean': np.nan,
                        'log2FoldChange': np.nan,
                        'padj': np.nan
                    })
                
                integrated_data.append(gene_data)
        
        result_df = pd.DataFrame(integrated_data)
        logger.info(f"Integration complete. Found {len(result_df)} gene-CpG associations")
        logger.info("Sample of integrated results:")
        logger.info(result_df.head().to_string())
        return result_df
    
    def plot_integration_summary(self, integrated_results: pd.DataFrame, output_dir: str):
        """Create summary plots of the integrated analysis"""
        os.makedirs(output_dir, exist_ok=True)
        
        # Plot 1: MeCP2 Enrichment Distribution
        plt.figure(figsize=(10, 6))
        sns.histplot(data=integrated_results, x='mecp2_enrichment', bins=50)
        plt.title('Distribution of MeCP2 Enrichment at CpG Islands')
        plt.xlabel('MeCP2 Enrichment (Exo/Endo)')
        plt.ylabel('Count')
        plt.savefig(os.path.join(output_dir, 'mecp2_enrichment_distribution.pdf'))
        plt.close()
        
        # Plot 2: MeCP2 Enrichment vs Expression Change
        plt.figure(figsize=(10, 6))
        mask = ~integrated_results['log2FoldChange'].isna()
        sns.scatterplot(
            data=integrated_results[mask],
            x='log2FoldChange',
            y='mecp2_enrichment',
            alpha=0.6
        )
        plt.title('MeCP2 Enrichment vs Gene Expression Change')
        plt.xlabel('log2 Fold Change (RNA-seq)')
        plt.ylabel('MeCP2 Enrichment (Exo/Endo)')
        plt.savefig(os.path.join(output_dir, 'mecp2_vs_expression.pdf'))
        plt.close()

def load_gtf(gtf_file: str) -> pd.DataFrame:
    """Load and parse GTF file, extracting only gene records"""
    logger.info(f"Loading GTF file: {gtf_file}")
    
    # Read GTF file with correct format
    df = pd.read_csv(gtf_file, sep='\t', comment='#',
                     names=['chr', 'source', 'feature', 'start', 'end', 
                           'score', 'strand', 'frame', 'attributes'])
    
    # Filter for gene entries only
    genes = df[df['feature'] == 'gene'].copy()
    logger.info(f"Found {len(genes)} gene records")
    
    # Extract gene_name from attributes
    def extract_gene_name(attr_str):
        for attr in attr_str.split('; '):
            if attr.startswith('gene_name'):
                return attr.split('"')[1]
        return None
    
    genes['gene_name'] = genes['attributes'].apply(extract_gene_name)
    
    # Select required columns
    result = genes[['chr', 'start', 'end', 'gene_name']].copy()
    logger.info(f"Processed gene annotations shape: {result.shape}")
    logger.info("Sample of processed annotations:")
    logger.info(result.head().to_string())
    
    return result

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Integrate MeCP2 enrichment with RNA-seq data')
    parser.add_argument('--enrichment-file', type=str, required=True,
                       help='Path to MeCP2 enrichment results')
    parser.add_argument('--rna-seq-file', type=str, required=True,
                       help='Path to RNA-seq results')
    parser.add_argument('--gene-annotations', type=str, required=True,
                       help='Path to gene annotations file')
    parser.add_argument('--output-dir', type=str, required=True,
                       help='Output directory')
    args = parser.parse_args()
    
    try:
        # Initialize integrator
        integrator = RNASeqIntegrator(args.enrichment_file, args.rna_seq_file)
        
        # Load gene annotations using custom GTF parser
        logger.info("Loading gene annotations...")
        gene_annotations = load_gtf(args.gene_annotations)
        
        # Integrate data
        logger.info("Integrating MeCP2 enrichment with gene annotations and RNA-seq data...")
        integrated_results = integrator.integrate_data(gene_annotations)
        
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Save full results
        integrated_results.to_csv(
            os.path.join(args.output_dir, 'mecp2_enriched_genes.csv'),
            index=False
        )
        
        # Save MeCP2-enriched genes only (without RNA-seq data)
        mecp2_enriched = integrated_results[['gene', 'chr', 'cpg_start', 'cpg_end', 
                                           'mecp2_enrichment', 'exo_signal', 'endo_signal']]
        mecp2_enriched = mecp2_enriched.drop_duplicates()
        mecp2_enriched.to_csv(
            os.path.join(args.output_dir, 'mecp2_enriched_genes_only.csv'),
            index=False
        )
        
        # Create visualizations
        integrator.plot_integration_summary(integrated_results, args.output_dir)
        
        # Generate summary statistics
        logger.info("\nAnalysis Summary:")
        logger.info(f"Total enriched CpG regions: {len(integrator.significant_regions)}")
        logger.info(f"Total genes near enriched regions: {len(mecp2_enriched)}")
        logger.info(f"Genes with RNA-seq data: {integrated_results['log2FoldChange'].notna().sum()}")
        
    except Exception as e:
        logger.error(f"Error during integration: {str(e)}")
        raise

if __name__ == "__main__":
    main() 