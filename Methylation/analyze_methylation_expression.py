import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
from typing import Dict, Tuple, List
import logging
from pathlib import Path
import pyBigWig
import pysam
from concurrent.futures import ProcessPoolExecutor
from functools import partial, lru_cache

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class MethylationExpressionAnalyzer:
    def __init__(self, output_dir: str = "methylation_expression_analysis", debug: bool = False):
        """
        Initialize the analyzer with output directory.
        
        Parameters:
        -----------
        output_dir : str
            Directory for output files
        debug : bool
            If True, runs analysis on a subset of data for faster debugging
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.debug = debug
        
        # Set plotting style
        plt.style.use('seaborn-v0_8')
        self.colors = {
            'UP': '#ff7f0e',
            'DOWN': '#1f77b4',
            'UNCHANGED': '#2ca02c'
        }
        
        # Define conditions
        self.conditions = ['exogenous', 'endogenous']
        
        # Configure debug logging
        if self.debug:
            logging.getLogger().setLevel(logging.DEBUG)
            logger.info("Running in DEBUG mode")

    def load_data(self, 
                 expression_files: Dict[str, str],
                 medip_dir: str,
                 cpg_islands_file: str,
                 gtf_file: str,
                 genome_fasta: str) -> None:
        """
        Load and prepare the data for analysis.
        
        Parameters:
        -----------
        expression_files : Dict[str, str]
            Dictionary with cell types as keys and paths to expression data
        medip_dir : str
            Path to directory containing MeDIP bigWig files
        cpg_islands_file : str
            Path to CpG islands BED file
        gtf_file : str
            Path to GTF annotation file
        genome_fasta : str
            Path to genome FASTA file
        """
        logger.info("Loading data...")
        
        # In debug mode, process only a subset of chromosomes
        self.debug_chroms = ['chr1', 'chr2'] if self.debug else None
        
        # Load GTF file for gene annotations and create gene ID mapping
        logger.info("Reading GTF file...")
        self.gene_info = self._parse_gtf(gtf_file)
        self.gene_id_map = self._create_gene_id_map(self.gene_info)
        
        # Load expression data
        self.expression_data = {}
        for cell_type, file_path in expression_files.items():
            df = pd.read_csv(file_path)
            
            # Add Ensembl IDs to expression data
            df['ensembl_id'] = df['gene'].map(self.gene_id_map.get)
            
            # Classify genes based on expression
            df['regulation'] = 'UNCHANGED'
            df.loc[(df['log2FoldChange'] > 1) & (df['padj'] < 0.05), 'regulation'] = 'UP'
            df.loc[(df['log2FoldChange'] < -1) & (df['padj'] < 0.05), 'regulation'] = 'DOWN'
            
            self.expression_data[cell_type] = df
        
        # Load and classify CpG islands
        self.cpg_df = self.classify_cpg_islands(gtf_file, cpg_islands_file)
        
        if self.debug:
            # Subset data for debugging
            self.cpg_df = self.cpg_df.head(1000)
            logger.debug(f"Debug mode: Using {len(self.cpg_df)} CpG islands")
        
        # Define cell types and conditions mapping
        cell_types_mapping = {
            'NEU': 'N',  # Neurons
            'NSC': 'PP'  # Neural Stem Cells
        }
        
        # Load methylation data for each cell type and condition
        self.methylation_data = {}
        for cell_type, medip_prefix in cell_types_mapping.items():
            self.methylation_data[cell_type] = {}
            
            # Process exogenous condition
            logger.info(f"Processing {cell_type} exogenous condition...")
            exo_replicates = [
                os.path.join(medip_dir, f"Medip_{medip_prefix}_output_r{i}.bw")
                for i in range(1, 4)  # r1, r2, r3
            ]
            
            # Calculate methylation using parallel processing for exogenous
            exo_methylation = process_cell_type_parallel(
                cell_type=cell_type,
                medip_prefix=medip_prefix,
                mecp2_df=self.cpg_df,
                medip_dir=medip_dir,
                genome_fasta=genome_fasta
            )
            self.methylation_data[cell_type]['exogenous'] = exo_methylation
            
            # Process endogenous condition similarly
            logger.info(f"Processing {cell_type} endogenous condition...")
            endo_replicates = [
                os.path.join(medip_dir, f"Medip_{medip_prefix}_output_r{i}.bw")
                for i in range(1, 4)
            ]
            
            # Calculate methylation for endogenous
            endo_methylation = process_cell_type_parallel(
                cell_type=cell_type,
                medip_prefix=medip_prefix,
                mecp2_df=self.cpg_df,
                medip_dir=medip_dir,
                genome_fasta=genome_fasta
            )
            self.methylation_data[cell_type]['endogenous'] = endo_methylation

    def _parse_gtf(self, gtf_file: str) -> pd.DataFrame:
        """
        Parse GTF file to extract gene information.
        """
        genes = []
        with open(gtf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if fields[2] == 'gene':
                    info = dict(item.strip().split(' ') for item in fields[8].rstrip(';').split('; '))
                    genes.append({
                        'gene_id': info['gene_id'].strip('"'),
                        'gene_name': info['gene_name'].strip('"'),
                        'chr': fields[0],
                        'start': int(fields[3]),
                        'end': int(fields[4]),
                        'strand': fields[6]
                    })
        return pd.DataFrame(genes)

    def _create_gene_id_map(self, gene_info: pd.DataFrame) -> Dict[str, str]:
        """
        Create a mapping between gene symbols and Ensembl IDs.
        """
        return dict(zip(gene_info['gene_name'], gene_info['gene_id']))

    def analyze_methylation_differences(self) -> Dict:
        """
        Analyze methylation differences between conditions and gene expression.
        """
        results = {}
        
        # Create plots directory
        plots_dir = self.output_dir / 'plots'
        plots_dir.mkdir(exist_ok=True)
        
        # Process each cell type
        for cell_type in ['NEU', 'NSC']:
            results[cell_type] = {}
            
            # Get expression data for this cell type
            expr_data = self.expression_data[cell_type]
            logger.info(f"\nExpression data for {cell_type}:")
            logger.info(f"Total genes: {len(expr_data)}")
            logger.info(f"Sample gene IDs: {expr_data['gene'].head().tolist()}")
            
            # Get methylation data for this cell type
            try:
                if cell_type not in self.methylation_data:
                    logger.error(f"Could not find methylation data for {cell_type}")
                    logger.debug(f"Available keys: {list(self.methylation_data.keys())}")
                    continue
                    
                cell_type_data = self.methylation_data[cell_type]
                
                # Process each condition
                for condition in ['exogenous', 'endogenous']:
                    logger.info(f"Processing {cell_type} - {condition}")
                    
                    if condition not in cell_type_data:
                        logger.error(f"Could not find {condition} data for {cell_type}")
                        logger.debug(f"Available conditions: {list(cell_type_data.keys())}")
                        continue
                    
                    # Get methylation data for this condition
                    condition_data = cell_type_data[condition]
                    
                    # Add location information if not present
                    if 'location' not in condition_data.columns:
                        condition_data = self.classify_cpg_islands(condition_data)
                    
                    # Split into promoter and gene body regions
                    promoter_data = condition_data[condition_data['location'] == 'promoter'].copy()
                    gene_body_data = condition_data[condition_data['location'] == 'gene_body'].copy()
                    
                    logger.info(f"\nMethylation data for {cell_type} - {condition}:")
                    logger.info(f"Total CpG islands: {len(condition_data)}")
                    logger.info(f"Promoter CpGs: {len(promoter_data)}")
                    logger.info(f"Gene body CpGs: {len(gene_body_data)}")
                    if 'associated_gene' in condition_data.columns:
                        logger.info(f"Sample associated genes: {condition_data['associated_gene'].head().tolist()}")

                    # Calculate average methylation per gene
                    if len(promoter_data) > 0:
                        promoter_grouped = promoter_data.groupby('associated_gene')['methylation'].agg(['mean', 'std']).reset_index()
                        promoter_grouped.columns = ['gene_id', 'promoter_mean', 'promoter_std']
                    else:
                        promoter_grouped = pd.DataFrame(columns=['gene_id', 'promoter_mean', 'promoter_std'])

                    if len(gene_body_data) > 0:
                        gene_body_grouped = gene_body_data.groupby('associated_gene')['methylation'].agg(['mean', 'std']).reset_index()
                        gene_body_grouped.columns = ['gene_id', 'gene_body_mean', 'gene_body_std']
                    else:
                        gene_body_grouped = pd.DataFrame(columns=['gene_id', 'gene_body_mean', 'gene_body_std'])

                    # Match with expression data using Ensembl IDs
                    expr_with_meth = expr_data.merge(
                        promoter_grouped, 
                        left_on='ensembl_id',  # Changed from 'gene' to 'ensembl_id'
                        right_on='gene_id', 
                        how='left'
                    )
                    expr_with_meth = expr_with_meth.merge(
                        gene_body_grouped, 
                        left_on='ensembl_id',  # Changed from 'gene' to 'ensembl_id'
                        right_on='gene_id', 
                        how='left'
                    )

                    logger.info("\nGene matching info:")
                    logger.info(f"Genes with promoter data: {expr_with_meth['promoter_mean'].notna().sum()}")
                    logger.info(f"Genes with gene body data: {expr_with_meth['gene_body_mean'].notna().sum()}")

                    # Calculate statistics for each regulation group
                    results[cell_type][condition] = self._calculate_group_statistics(expr_with_meth)
                    
            except Exception as e:
                logger.error(f"Error processing {cell_type}: {str(e)}")
                if self.debug:
                    import traceback
                    logger.error(traceback.format_exc())
                continue

        return results

    def _analyze_region_methylation(self, data: pd.DataFrame, region: str) -> Dict:
        """
        Analyze methylation data for a specific region.
        """
        results = {}
        methylation_col = f'{region}_methylation'
        
        # Calculate statistics for each regulation group
        for reg in ['UP', 'DOWN', 'UNCHANGED']:
            group_data = data[data['regulation'] == reg][methylation_col]
            if len(group_data) > 0:
                results[reg] = {
                    'n': len(group_data),
                    'mean': group_data.mean(),
                    'median': group_data.median(),
                    'std': group_data.std()
                }
        
        # Calculate statistical tests between UP and DOWN regulated genes
        if 'UP' in results and 'DOWN' in results:
            up_data = data[data['regulation'] == 'UP'][methylation_col]
            down_data = data[data['regulation'] == 'DOWN'][methylation_col]
            
            stat, pval = stats.mannwhitneyu(up_data, down_data, alternative='two-sided')
            cohens_d = (up_data.mean() - down_data.mean()) / np.sqrt((up_data.var() + down_data.var()) / 2)
            
            results['statistical_tests'] = {
                'mannwhitney_pval': pval,
                'cohens_d': cohens_d
            }
        
        return results

    def _create_methylation_plots(self, data: pd.DataFrame, cell_type: str,
                                 condition: str):
        """
        Create visualization plots for methylation analysis.
        """
        # Create condition-specific directory
        condition_dir = self.output_dir / cell_type / condition
        condition_dir.mkdir(parents=True, exist_ok=True)
        
        # Create plots
        self._plot_methylation_distribution(data, cell_type, condition, condition_dir)
        self._plot_methylation_boxplots(data, cell_type, condition, condition_dir)
        self._plot_methylation_expression_scatter(data, cell_type, condition, condition_dir)
        self._plot_methylation_violin(data, cell_type, condition, condition_dir)
        
        # Create combined plots comparing promoter vs gene body
        self._plot_combined_methylation(data, cell_type, condition, condition_dir)

    def _plot_combined_methylation(self, data: pd.DataFrame, cell_type: str, 
                                 condition: str, output_dir: Path):
        """Create combined plots comparing promoter and gene body methylation."""
        # Prepare data in long format
        plot_data = pd.melt(
            data,
            id_vars=['regulation', 'log2FoldChange'],
            value_vars=['promoter_methylation', 'gene_body_methylation'],
            var_name='region',
            value_name='methylation'
        )
        plot_data['region'] = plot_data['region'].map({
            'promoter_methylation': 'Promoter',
            'gene_body_methylation': 'Gene Body'
        })
        
        # Create violin plot
        plt.figure(figsize=(12, 6))
        sns.violinplot(data=plot_data, x='regulation', y='methylation',
                      hue='region', split=True,
                      order=['UP', 'DOWN', 'UNCHANGED'])
        plt.title(f'{cell_type} - {condition}: Methylation Comparison')
        plt.xlabel('Gene Regulation')
        plt.ylabel('Methylation Level (%)')
        plt.tight_layout()
        plt.savefig(output_dir / 'combined_methylation_violin.pdf')
        plt.close()

    def _plot_methylation_distribution(self, data: pd.DataFrame, region_type: str):
        """
        Create violin and box plots for methylation distribution.
        
        Parameters:
        -----------
        data : pd.DataFrame
            DataFrame containing methylation data and regulation groups
        region_type : str
            Type of region ('promoter' or 'gene_body')
        """
        methylation_col = f'{region_type}_methylation'
        
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        
        # Violin plot
        sns.violinplot(data=data, x='regulation', y=methylation_col, 
                      order=['UP', 'DOWN', 'UNCHANGED'],
                      palette=self.colors, ax=ax1)
        ax1.set_title(f'{region_type.title()} Methylation Distribution')
        ax1.set_ylabel('Methylation Level')
        
        # Box plot
        sns.boxplot(data=data, x='regulation', y=methylation_col,
                    order=['UP', 'DOWN', 'UNCHANGED'],
                    palette=self.colors, ax=ax2)
        ax2.set_title(f'{region_type.title()} Methylation Box Plot')
        ax2.set_ylabel('Methylation Level')
        
        # Add statistical annotation if applicable
        if len(data[data['regulation'] == 'UP']) > 0 and len(data[data['regulation'] == 'DOWN']) > 0:
            stat, pval = stats.mannwhitneyu(
                data[data['regulation'] == 'UP'][methylation_col],
                data[data['regulation'] == 'DOWN'][methylation_col],
                alternative='two-sided'
            )
            ax2.text(0.5, 1.1, f'Mann-Whitney p = {pval:.2e}',
                    horizontalalignment='center', transform=ax2.transAxes)
        
        # Save the plot
        plot_dir = self.output_dir / 'plots'
        plot_dir.mkdir(exist_ok=True)
        plt.tight_layout()
        plt.savefig(plot_dir / f'{region_type}_methylation_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()

    def _plot_methylation_boxplots(self, data: pd.DataFrame, cell_type: str,
                                 condition: str, output_dir: Path):
        """Create box plots comparing methylation levels."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Promoter methylation
        sns.boxplot(data=data, x='regulation', y='promoter_methylation',
                   order=['UP', 'DOWN', 'UNCHANGED'],
                   palette=self.colors, ax=ax1)
        ax1.set_title(f'{cell_type} - {condition}: Promoter Methylation')
        
        # Gene body methylation
        sns.boxplot(data=data, x='regulation', y='gene_body_methylation',
                   order=['UP', 'DOWN', 'UNCHANGED'],
                   palette=self.colors, ax=ax2)
        ax2.set_title(f'{cell_type} - {condition}: Gene Body Methylation')
        
        for ax in [ax1, ax2]:
            ax.set_xlabel('Gene Regulation')
            ax.set_ylabel('Methylation Level (%)')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'methylation_boxplots.pdf')
        plt.close()

    def _plot_methylation_expression_scatter(self, data: pd.DataFrame, cell_type: str,
                                          condition: str, output_dir: Path):
        """Create scatter plots of methylation vs expression."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Promoter methylation vs expression
        sns.scatterplot(data=data, x='promoter_methylation', y='log2FoldChange',
                       hue='regulation', palette=self.colors, alpha=0.6, ax=ax1)
        ax1.set_title(f'{cell_type} - {condition}: Promoter Methylation vs Expression')
        ax1.set_xlabel('Promoter Methylation Level (%)')
        ax1.set_ylabel('log2 Fold Change')
        
        # Gene body methylation vs expression
        sns.scatterplot(data=data, x='gene_body_methylation', y='log2FoldChange',
                       hue='regulation', palette=self.colors, alpha=0.6, ax=ax2)
        ax2.set_title(f'{cell_type} - {condition}: Gene Body Methylation vs Expression')
        ax2.set_xlabel('Gene Body Methylation Level (%)')
        ax2.set_ylabel('log2 Fold Change')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'methylation_expression_scatter.pdf')
        plt.close()

    def _plot_methylation_violin(self, data: pd.DataFrame, cell_type: str,
                               condition: str, output_dir: Path):
        """Create violin plots for methylation distribution."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Promoter methylation
        sns.violinplot(data=data, x='regulation', y='promoter_methylation',
                      order=['UP', 'DOWN', 'UNCHANGED'],
                      palette=self.colors, ax=ax1)
        ax1.set_title(f'{cell_type} - {condition}: Promoter Methylation')
        
        # Gene body methylation
        sns.violinplot(data=data, x='regulation', y='gene_body_methylation',
                      order=['UP', 'DOWN', 'UNCHANGED'],
                      palette=self.colors, ax=ax2)
        ax2.set_title(f'{cell_type} - {condition}: Gene Body Methylation')
        
        for ax in [ax1, ax2]:
            ax.set_xlabel('Gene Regulation')
            ax.set_ylabel('Methylation Level (%)')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'methylation_violin.pdf')
        plt.close()

    def classify_cpg_islands(self, gtf_file: str, cpg_file: str) -> pd.DataFrame:
        """
        Classify CpG islands based on their genomic location relative to genes.
        """
        logger.info("Reading GTF file...")
        
        # Read GTF file
        genes_df = pd.read_csv(gtf_file, sep='\t', comment='#',
                              names=['chr', 'source', 'feature', 'start', 'end',
                                    'score', 'strand', 'frame', 'attributes'])
        
        # Filter for genes and debug chromosomes if in debug mode
        genes_df = genes_df[genes_df['feature'] == 'gene']
        if self.debug_chroms:
            genes_df = genes_df[genes_df['chr'].isin(self.debug_chroms)]
            logger.debug(f"Debug mode: Using chromosomes {self.debug_chroms}")
        
        # Extract and clean gene_id
        genes_df['gene_id'] = genes_df['attributes'].str.extract('gene_id "([^"]*)"')
        genes_df['gene_id'] = genes_df['gene_id'].str.split('.').str[0]
        
        logger.debug(f"Total genes from GTF: {len(genes_df)}")
        logger.debug(f"Sample gene IDs: {genes_df['gene_id'].head().tolist()}")
        
        # Read CpG islands
        logger.info("Reading CpG islands...")
        cpg_df = pd.read_csv(cpg_file, sep='\t',
                            names=['chr', 'start', 'end', 'id', 'cpg_info', 'cpg_count'])
        
        if self.debug_chroms:
            cpg_df = cpg_df[cpg_df['chr'].isin(self.debug_chroms)]
        
        logger.debug(f"Total CpG islands: {len(cpg_df)}")
        
        # Add progress bar for classification
        from tqdm import tqdm
        
        # Create empty columns for classification
        cpg_df['location'] = 'intergenic'
        cpg_df['associated_gene'] = None
        cpg_df['strand'] = '+'
        
        # Define promoter distance
        promoter_distance = 2000
        
        # Classify each CpG island with progress bar
        logger.info("Classifying CpG islands...")
        for idx, cpg in tqdm(cpg_df.iterrows(), total=len(cpg_df), 
                            desc="Classifying CpG islands"):
            # Get genes on the same chromosome
            chr_genes = genes_df[genes_df['chr'] == cpg['chr']]
            
            if idx % 1000 == 0:
                logger.debug(f"Processing CpG island {idx}/{len(cpg_df)}")
                logger.debug(f"Current CpG: {cpg['chr']}:{cpg['start']}-{cpg['end']}")
            
            for _, gene in chr_genes.iterrows():
                # Check if CpG island overlaps with promoter region
                promoter_start = gene['start'] - promoter_distance if gene['strand'] == '+' else gene['end'] - promoter_distance
                promoter_end = gene['start'] + promoter_distance if gene['strand'] == '+' else gene['end'] + promoter_distance
                
                if (cpg['start'] <= promoter_end and cpg['end'] >= promoter_start):
                    cpg_df.at[idx, 'location'] = 'promoter'
                    cpg_df.at[idx, 'associated_gene'] = gene['gene_id']
                    cpg_df.at[idx, 'strand'] = gene['strand']
                    break
                
                # Check if CpG island overlaps with gene body
                elif (cpg['start'] <= gene['end'] and cpg['end'] >= gene['start']):
                    cpg_df.at[idx, 'location'] = 'gene_body'
                    cpg_df.at[idx, 'associated_gene'] = gene['gene_id']
                    cpg_df.at[idx, 'strand'] = gene['strand']
                    break
        
        # Add summary statistics
        location_stats = cpg_df['location'].value_counts()
        logger.info("\nCpG Island Classification Summary:")
        logger.info(f"Total CpG islands: {len(cpg_df)}")
        logger.info(f"Promoter CpGs: {location_stats.get('promoter', 0)}")
        logger.info(f"Gene body CpGs: {location_stats.get('gene_body', 0)}")
        logger.info(f"Intergenic CpGs: {location_stats.get('intergenic', 0)}")
        
        return cpg_df

    def _calculate_group_statistics(self, data: pd.DataFrame) -> Dict:
        """
        Calculate methylation statistics and create plots for each region type.
        """
        results = {}
        
        # Ensure we have the regulation column
        if 'regulation' not in data.columns:
            logger.error("No regulation information found in data")
            return results
        
        # Analyze promoter methylation
        if 'promoter_mean' in data.columns:
            promoter_data = data[['regulation', 'promoter_mean']].copy()
            promoter_data.rename(columns={'promoter_mean': 'promoter_methylation'}, inplace=True)
            promoter_data = promoter_data.dropna()
            if len(promoter_data) > 0:
                results['promoter'] = self._analyze_region_methylation(promoter_data, 'promoter')
                self._plot_methylation_distribution(promoter_data, 'promoter')

        # Analyze gene body methylation
        if 'gene_body_mean' in data.columns:
            gene_body_data = data[['regulation', 'gene_body_mean']].copy()
            gene_body_data.rename(columns={'gene_body_mean': 'gene_body_methylation'}, inplace=True)
            gene_body_data = gene_body_data.dropna()
            if len(gene_body_data) > 0:
                results['gene_body'] = self._analyze_region_methylation(gene_body_data, 'gene_body')
                self._plot_methylation_distribution(gene_body_data, 'gene_body')
        
        return results

def process_cell_type_parallel(cell_type: str, medip_prefix: str, 
                             mecp2_df: pd.DataFrame, medip_dir: str, 
                             genome_fasta: str) -> pd.DataFrame:
    """
    Process methylation data for a cell type in parallel
    """
    replicates = ['r1', 'r2', 'r3']
    
    # Create partial function for parallel processing
    process_replicate = partial(
        calculate_methylation_batch,
        regions_df=mecp2_df,
        fasta_file=genome_fasta
    )
    
    # Process replicates in parallel
    bw_files = [
        os.path.join(medip_dir, f"Medip_{medip_prefix}_output_{rep}.bw")
        for rep in replicates
    ]
    
    with ProcessPoolExecutor() as executor:
        replicate_results = list(executor.map(process_replicate, bw_files))
    
    # Combine results
    result_df = mecp2_df.copy()
    
    # Calculate mean methylation and coverage across replicates
    result_df['methylation'] = np.mean([
        rep['methylation'].values for rep in replicate_results
    ], axis=0)
    
    result_df['coverage'] = np.mean([
        rep['coverage'].values for rep in replicate_results
    ], axis=0)
    
    # Copy other metrics from first replicate (they should be the same across replicates)
    result_df['cpg_count'] = replicate_results[0]['cpg_count']
    result_df['cpg_density'] = replicate_results[0]['cpg_density']
    
    return result_df

def calculate_methylation_batch(bw_file: str, regions_df: pd.DataFrame, 
                              fasta_file: str) -> pd.DataFrame:
    """
    Calculate methylation levels and coverage for multiple regions efficiently
    with proper signal normalization
    """
    bw = pyBigWig.open(bw_file)
    fasta = pysam.FastaFile(fasta_file)
    
    # Create a copy of the input DataFrame to preserve original index
    result_df = regions_df.copy()
    
    # Pre-filter valid chromosomes but keep track of invalid ones
    valid_chroms = set(fasta.references) & set(regions_df['chr'].unique())
    valid_mask = regions_df['chr'].isin(valid_chroms)
    
    # Initialize all values with zeros/NaN
    result_df['methylation'] = 0.0
    result_df['coverage'] = 0.0
    result_df['cpg_density'] = 0.0
    
    # Process only valid chromosomes
    valid_regions = regions_df[valid_mask].copy()
    
    if len(valid_regions) > 0:
        # Calculate CpG density using the count from BED file
        valid_regions['region_length'] = valid_regions['end'] - valid_regions['start']
        valid_regions['cpg_density'] = (valid_regions['cpg_count'].astype(float) * 100) / valid_regions['region_length']
        
        # Calculate methylation and coverage in batches
        batch_size = 1000
        methylation_values = []
        coverage_values = []
        
        # Get global statistics for normalization
        all_values = []
        for _, row in valid_regions.iterrows():
            vals = bw.values(row['chr'], row['start'], row['end'])
            if vals and any(v is not None for v in vals):
                all_values.extend([v for v in vals if v is not None])
        
        if all_values:  # Only proceed if we have values
            # Calculate normalization factors
            global_median = np.median(all_values)
            global_99th = np.percentile(all_values, 99)
            
            for i in range(0, len(valid_regions), batch_size):
                batch = valid_regions.iloc[i:i+batch_size]
                values = [
                    bw.values(row['chr'], row['start'], row['end'])
                    for _, row in batch.iterrows()
                ]
                
                # Calculate methylation and coverage for each region
                for vals in values:
                    if vals is None or len(vals) == 0:
                        methylation_values.append(0)
                        coverage_values.append(0)
                    else:
                        valid_vals = [v for v in vals if v is not None]
                        total_positions = len(vals)
                        covered_positions = len(valid_vals)
                        
                        if covered_positions > 0:
                            # Normalize values to 0-100 range
                            normalized_vals = np.clip(
                                [100 * (v / global_99th) for v in valid_vals],
                                0, 100
                            )
                            methylation = np.mean(normalized_vals)
                        else:
                            methylation = 0
                        
                        coverage = covered_positions / total_positions if total_positions > 0 else 0
                        
                        methylation_values.append(methylation)
                        coverage_values.append(coverage)
            
            # Update values for valid regions
            valid_regions['methylation'] = methylation_values
            valid_regions['coverage'] = coverage_values
            
            # Update the result DataFrame with calculated values
            result_df.loc[valid_regions.index, 'methylation'] = valid_regions['methylation']
            result_df.loc[valid_regions.index, 'coverage'] = valid_regions['coverage']
            result_df.loc[valid_regions.index, 'cpg_density'] = valid_regions['cpg_density']
    
    # Cleanup
    bw.close()
    fasta.close()
    
    return result_df

def count_cpgs_in_sequence(sequence: str) -> int:
    """
    Count CpG dinucleotides in a DNA sequence with caching
    """
    return sequence.upper().count('CG')

os.chdir("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Methylation")

def main():
    # Add argument parsing for debug mode
    import argparse
    parser = argparse.ArgumentParser(description='Analyze methylation and expression data')
    parser.add_argument('--debug', action='store_true', help='Run in debug mode with subset of data')
    args = parser.parse_args()
    
    os.chdir("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Methylation")

    # Define input files
    expression_files = {
        'NEU': '../iterative_alternative/DATA/DEA_NEU.csv',
        'NSC': '../iterative_alternative/DATA/DEA_NSC.csv'
    }
    
    medip_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeDIP/output_old/bigwig"
    cpg_islands_file = "../DATA/cpg_islands.bed"
    gtf_file = "../DATA/gencode.vM10.annotation.gtf"
    genome_fasta = "../DATA/mm10.fa"
    output_dir = "methylation_expression_analysis"
    
    try:
        # Initialize analyzer with debug mode if specified
        analyzer = MethylationExpressionAnalyzer(output_dir, debug=args.debug)
        
        # Load data using MeDIP bigWig files
        analyzer.load_data(
            expression_files=expression_files,
            medip_dir=medip_dir,
            cpg_islands_file=cpg_islands_file,
            gtf_file=gtf_file,
            genome_fasta=genome_fasta
        )
        
        # Perform analysis
        results = analyzer.analyze_methylation_differences()
        
        # Save results
        results_file = Path(output_dir) / 'methylation_analysis_results.txt'
        with open(results_file, 'w') as f:
            for cell_type, conditions in results.items():
                f.write(f"\n{cell_type} Results:\n")
                f.write("================\n")
                
                for condition, regions in conditions.items():
                    f.write(f"\n{condition.upper()}:\n")
                    f.write("-" * len(condition) + "\n")
                    
                    for region_type, stats in regions.items():
                        f.write(f"\n{region_type.upper()}:\n")
                        f.write("-" * len(region_type) + "\n")
                        
                        # Write basic statistics
                        for reg in ['UP', 'DOWN', 'UNCHANGED']:
                            if reg in stats:
                                f.write(f"\n{reg} regulated genes:\n")
                                f.write(f"Number of genes: {stats[reg]['n']}\n")
                                f.write(f"Mean methylation: {stats[reg]['mean']:.3f}\n")
                                f.write(f"Median methylation: {stats[reg]['median']:.3f}\n")
                                f.write(f"Standard deviation: {stats[reg]['std']:.3f}\n")
                        
                        # Write statistical tests if they exist
                        if 'statistical_tests' in stats:
                            f.write("\nStatistical Tests:\n")
                            f.write(f"Mann-Whitney U p-value: {stats['statistical_tests']['mannwhitney_pval']:.2e}\n")
                            f.write(f"Effect size (Cohen's d): {stats['statistical_tests']['cohens_d']:.3f}\n")
                    f.write("\n" + "="*50 + "\n")
        
        logger.info("Analysis completed successfully")
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        if args.debug:
            import traceback
            logger.error(traceback.format_exc())
        raise

if __name__ == "__main__":
    main() 