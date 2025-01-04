import os
import logging
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List, Dict
import argparse

# Get experiment name from command line arguments if running as main script
def get_experiment_name():
    try:
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument('--experiment', type=str, default='align2_005_consistent_peaks')
        args, _ = parser.parse_known_args()
        return args.experiment
    except:
        return 'align2_005_consistent_peaks'

EXPERIMENT = get_experiment_name()

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
WORKING_DIR = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Methylation_dev"
os.chdir(WORKING_DIR)

# Configuration settings
CONFIG = {
    'methylation_thresholds': {
        'hypo': 30,
        'hyper': 70,
        'intermediate': (30, 70),
        'high_confidence': 0.8  # Minimum confidence score for high-quality measurements
    },
    'quality_thresholds': {
        'coverage': 0.7,
        'min_cpgs': 5,
        'min_signal_to_noise': 2.0,
        'min_coverage_uniformity': 0.6,
        'max_local_bias': 0.3
    },
    'expression_thresholds': {
        'log2fc': 0.5,
        'padj': 0.05,
        'min_expression': 1.0  # Minimum expression level to consider
    },
    'genomic_regions': {
        'promoter': (-2000, 500),  # relative to TSS
        'gene_body_start': 500,  # after TSS
        'cpg_island_shore': 2000,  # bp from CpG island
        'cpg_island_shelf': 4000   # bp from shore
    },
    'methylation_analysis': {
        'window_size': 50,  # bp for sliding window analysis
        'smoothing_window': 3,  # windows for signal smoothing
        'min_region_size': 100,  # minimum size for methylated regions
        'background_percentile': 25  # percentile for local background
    },
    'replicates': {
        'min_replicates': 2,  # minimum number of valid replicates
        'max_replicate_cv': 0.3,  # maximum coefficient of variation
        'correlation_threshold': 0.7  # minimum correlation between replicates
    },
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
        'expression_file': 'expression_analysis.pkl',
        'validation_file': 'validation_results.pkl'
    },
    'validation': {
        'known_methylated_regions': {
            'imprinted_genes': ['Snrpn', 'H19', 'Igf2r'],
            'housekeeping_promoters': ['Actb', 'Gapdh', 'Hprt'],
            'tissue_specific': ['Neurod1', 'Gfap', 'Syp']
        },
        'expected_patterns': {
            'cpg_island_promoters': 'low',
            'gene_bodies': 'intermediate',
            'intergenic': 'high'
        }
    }
}

# File paths
PATHS = {
    'mecp2_dir': f"/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment_{EXPERIMENT}/NSC/mecp2_cpg_enrichment_parallel",
    'mecp2_file': "mecp2_cpg_enrichment_parallel.csv",
    'medip_dir': "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_medip/output_done/bigwig",
    'gtf_file': "../DATA/gencode.vM10.annotation.gtf",
    'genome_fasta': "../DATA/mm10.fa",
    'output_dir': f"analyze_mecp2_cpg_enrichment_{EXPERIMENT}",
    'cpg_islands_file': "../DATA/cpg_islands.bed",
    'rnaseq': {
        'NEU': '../iterative_alternative/DATA/DEA_NEU.csv',
        'NSC': '../iterative_alternative/DATA/DEA_NSC.csv'
    },
    'cache_dir': os.path.join(WORKING_DIR, "analysis_cache"),
    'validation_data': {
        'bisulfite_data': '../DATA/validation/bisulfite_data.csv',
        'cpg_islands': '../DATA/validation/cpg_islands.bed',
        'chromatin_state': '../DATA/validation/chromatin_states.bed'
    },
    'genome_features': {
        'cpg_islands': '../DATA/features/cpg_islands.bed',
        'regulatory_elements': '../DATA/features/regulatory_elements.bed',
        'chromatin_states': '../DATA/features/chromatin_states.bed'
    }
}

# Quality control parameters
QC_PARAMS = {
    'coverage': {
        'minimum': 5,  # Minimum read coverage
        'optimal': 30  # Optimal read coverage
    },
    'replicates': {
        'correlation_threshold': 0.8,
        'cv_threshold': 0.2
    },
    'signal': {
        'min_snr': 2.0,  # Minimum signal-to-noise ratio
        'background_method': 'local_median',
        'smoothing_window': 3
    }
}

# Biological validation parameters
BIO_VALIDATION = {
    'expected_patterns': {
        'active_promoters': {
            'methylation_range': (0, 30),
            'expression_correlation': 'negative'
        },
        'gene_bodies': {
            'methylation_range': (20, 80),
            'expression_correlation': 'positive'
        }
    },
    'control_regions': {
        'constitutive_unmethylated': ['housekeeping_promoters'],
        'constitutive_methylated': ['intergenic_repeats'],
        'tissue_specific': ['neural_enhancers']
    }
} 