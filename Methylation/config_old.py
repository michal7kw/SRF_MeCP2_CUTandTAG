import os
import logging
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Set default plotting style
sns.set_theme(style="whitegrid")
plt.rcParams['figure.figsize'] = [10, 6]
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300

# Set working directory
WORKING_DIR = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Methylation"
os.chdir(WORKING_DIR)

# Configuration settings
CONFIG = {
    'methylation_thresholds': {
        'hypo': 30,
        'hyper': 70
    },
    'quality_thresholds': {
        'coverage': 0.7,
        'min_cpgs': 5
    },
    'expression_thresholds': {
        'log2fc': 0.5,
        'padj': 0.05
    },
    'promoter_region': (-2000, 500),  # relative to TSS
    'gene_body_start': 500,  # after TSS
    'standard_chromosomes': [f'chr{i}' for i in range(1, 20)] + ['chrX', 'chrY'],
    'debug': {
        'enabled': False,
        'sample_size': None,
        'random_seed': 42,
        'cache_suffix': '_debug'
    },
    'cache': {
        'enabled': True,
        'directory': 'analysis_cache',
        'methylation_file': 'methylation_analysis.pkl',
        'binding_file': 'binding_analysis.pkl',
        'expression_file': 'expression_analysis.pkl'
    }
}

# File paths
PATHS = {
    'mecp2_dir': "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment_align1_005/NSC/mecp2_cpg_enrichment_parallel",
    'mecp2_file': "mecp2_cpg_enrichment_parallel.csv",
    'medip_dir': "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeDIP/output_old/bigwig",
    'gtf_file': "../DATA/gencode.vM10.annotation.gtf",
    'genome_fasta': "../DATA/mm10.fa",
    'output_dir': "analyze_mecp2_cpg_enrichment_align1_005",
    'cpg_islands_file': "../DATA/cpg_islands.bed",
    'rnaseq': {
        'NEU': '../iterative_alternative/DATA/DEA_NEU.csv',
        'NSC': '../iterative_alternative/DATA/DEA_NSC.csv'
    },
    'cache_dir': os.path.join(WORKING_DIR, "analysis_cache")
} 