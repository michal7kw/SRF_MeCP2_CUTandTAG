#%% Import required libraries
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import os
import pyBigWig
import pysam
import logging
from typing import Dict, List, Tuple, Any
import pyranges as pr
from functools import partial
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from tqdm import tqdm
import time
import argparse

# Set working directory
os.chdir('/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Methylation_dev')

import importlib
import functions
import config
import cache_utils
import func_cpg_islands_methylation
import func_tss_analysis
import func_regions_methylation
import functions_general

importlib.reload(functions)
importlib.reload(config) 
importlib.reload(cache_utils)
importlib.reload(func_cpg_islands_methylation)
importlib.reload(func_tss_analysis)
importlib.reload(func_regions_methylation)
importlib.reload(functions_general)

from functions import *
from config import *
from cache_utils import *
from func_cpg_islands_methylation import *
# from func_tss_analysis import *
from func_regions_methylation import *
from functions_general import *

#%% Configure debug mode and parse command line arguments
# Default arguments for notebook
class Args:
    def __init__(self):
        self.debug = False
        self.sample_size = 100
        self.force_recompute = True
        self.experiment = 'align1_005'
        self.processes = None
        self.cell_type = 'NSC'
        self.chromosome = None
args = Args()
debug = False
sample_size = 100
force_recompute = True
experiment = 'align1_005'
n_processes = None
cell_type = 'NSC'
chromosome = None

# Set experiment name in environment variable for config.py to use
os.environ['EXPERIMENT_NAME'] = args.experiment

# Reload config to update paths with new experiment name
import importlib
import config
importlib.reload(config)
from config import PATHS, CONFIG, logger

# Configure debug mode if specified
if args.debug:
    logger.info(f"Running in DEBUG mode with sample size: {args.sample_size}")
    CONFIG['debug'] = {
        'enabled': True,
        'sample_size': args.sample_size
    }

#%% Initialize output directory
os.makedirs(PATHS['output_dir'], exist_ok=True)

#%% ########################## Load MeCP2 binding data ############################################################
mecp2_binding = pd.read_csv(os.path.join(PATHS['mecp2_dir'], PATHS['mecp2_file']))
print("MeCP2 binding data:")
print(mecp2_binding.head())
# Filter by chromosome if specified
if chromosome:
    mecp2_binding = mecp2_binding[mecp2_binding['chr'] == chromosome]
    logger.info(f"Filtered data for chromosome {chromosome}")

"""mecp2_binding
MeCP2 binding data:
    chr    start      end  exo_signal  endo_signal  enrichment    pvalue  \
0  chr1  3531624  3531843    0.000000    15.177982    0.000000  1.000000   
1  chr1  3670619  3671074   29.287905     0.000000         inf  1.000000   
2  chr1  3671654  3672156   12.519179     6.646521    1.883569  0.193931   
3  chr1  4491701  4493673   18.136526     0.000000         inf  1.000000   
4  chr1  4571641  4572075   13.948418     0.000000         inf  1.000000   

  binding_type  peak_width_exo  peak_width_endo  significant  
0    endo_only            0.00            416.0        False  
1     exo_only          331.75              0.0        False  
2         both          378.00            336.0        False  
3     exo_only          243.00              0.0        False  
4     exo_only          317.50              0.0        False  
"""
mecp2_binding = mecp2_binding.rename(columns={'start': 'peak_start', 'end': 'peak_end'})

#%% ########################## Load gene annotations ############################################################
genes_df, gene_name_to_id = load_gene_annotations(PATHS['gtf_file'])
print(dict(list(gene_name_to_id.items())[1000:1003]))
print(genes_df.head())
# Filter by chromosome if specified
if chromosome:
    genes_df = genes_df[genes_df['chr'] == chromosome]
    logger.info(f"Filtered data for chromosome {chromosome}")

"""genes_df
{'Celrr': 'ENSMUSG00000097881', 'Gm29359': 'ENSMUSG00000101700', 'Gm37174': 'ENSMUSG00000104253'}
               gene_id      gene_name   chr    start      end strand  \
0   ENSMUSG00000102693  4933401J01Rik  chr1  3073252  3074322      +   
3   ENSMUSG00000064842        Gm26206  chr1  3102015  3102125      +   
24  ENSMUSG00000102851        Gm18956  chr1  3252756  3253236      +   
36  ENSMUSG00000089699         Gm1992  chr1  3466586  3513553      +   
43  ENSMUSG00000103147         Gm7341  chr1  3531794  3532720      +   

    promoter_start  promoter_end  gene_body_start  gene_body_end  
0          3071252       3073752          3073752        3074322  
3          3100015       3102515          3102515        3102125  
24         3250756       3253256          3253256        3253236  
36         3464586       3467086          3467086        3513553  
43         3529794       3532294          3532294        3532720  
"""

#%% ########################## Load RNA-seq expression data ############################################################
if cell_type:
    expression_data = {
        cell_type: load_expression_data(PATHS['rnaseq'][cell_type], gene_name_to_id, single_file=True)
    }
    logger.info(f"Analyzing cell type: {cell_type}")
else:
    expression_data = load_expression_data(PATHS['rnaseq'], gene_name_to_id)

print(expression_data.keys())
print(expression_data['NSC'].head())
expr_df = expression_data['NSC']

"""expr_df
dict_keys(['NSC'])
     gene      baseMean  log2FoldChange     lfcSE      stat    pvalue  \
0   Top2b  13465.464930       -0.129776  0.055509 -2.337916  0.019392   
1  Zfp113    974.648226       -0.204025  0.087200 -2.339752  0.019297   
3    Pfkm   2126.364649        0.232159  0.099092  2.342858  0.019137   
4   Syde2    306.400724       -0.346175  0.147748 -2.343007  0.019129   
5  Fam8a1   3567.752059        0.170360  0.072699  2.343368  0.019111   

       padj             gene_id expression_status  
0  0.049904  ENSMUSG00000017485         unchanged  
1  0.049670  ENSMUSG00000037007         unchanged  
3  0.049290  ENSMUSG00000033065         unchanged  
4  0.049280  ENSMUSG00000036863         unchanged  
5  0.049237  ENSMUSG00000069237         unchanged  
"""

#%% ########################## Print mapping statistics ############################################################
print(f"Total genes: {len(expr_df)}")
print(f"Genes with valid mapping: {expr_df['gene_id'].notna().sum()}")
print("\nSample of mapped genes:")
sample_df = expr_df[expr_df['gene_id'].notna()].head()
print(pd.DataFrame({
    'Gene Symbol': sample_df['gene'],
    'Ensembl ID': sample_df['gene_id'],
    'log2FC': sample_df['log2FoldChange'],
    'padj': sample_df['padj']
}))

"""
Total genes: 8025
Genes with valid mapping: 8025

Sample of mapped genes:
  Gene Symbol          Ensembl ID    log2FC      padj
0       Top2b  ENSMUSG00000017485 -0.129776  0.049904
1      Zfp113  ENSMUSG00000037007 -0.204025  0.049670
3        Pfkm  ENSMUSG00000033065  0.232159  0.049290
4       Syde2  ENSMUSG00000036863 -0.346175  0.049280
5      Fam8a1  ENSMUSG00000069237  0.170360  0.049237
"""

#%% ########################## Merge datasets ############################################################
merged_df = validate_and_merge_data(genes_df, expr_df, mecp2_binding)

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

"""
	gene_id	gene_name	chr	start	end	strand	promoter_start	promoter_end	gene_body_start	gene_body_end	tss	gene	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	expression_status	binding_type	exo_signal	endo_signal	peak_start	peak_end	peak_width_exo	peak_width_endo	binding_place	mecp2_bound
0	ENSMUSG00000025903	Lypla1	chr1	4807787	4848410	+	4805787	4808287	4808287	4848410	4807787	Lypla1	1606.272899	-0.605672	0.076884	-7.877686	3.335000e-15	3.664540e-14	downregulated	both	51.533421	37.104364	4807559.0	4808103.0	519.0	509.333333	promoter	True
1	ENSMUSG00000033813	Tcea1	chr1	4857813	4897909	+	4855813	4858313	4858313	4897909	4857813	Tcea1	3420.565606	-0.536336	0.060458	-8.871212	7.235410e-19	1.024950e-17	downregulated	both	509.548336	348.137860	4857465.0	4858372.0	1023.5	845.500000	promoter	True
2	ENSMUSG00000033793	Atp6v1h	chr1	5070017	5162529	+	5068017	5070517	5070517	5162529	5070017	Atp6v1h	2461.465043	-0.417803	0.070693	-5.910114	3.418710e-09	2.307570e-08	unchanged	both	49.335831	39.023866	5083039.0	5083536.0	676.5	531.000000	gene_body	True
3	ENSMUSG00000051285	Pcmtd1	chr1	7088919	7173628	+	7086919	7089419	7089419	7173628	7088919	Pcmtd1	2843.258046	-0.447877	0.069469	-6.447139	1.139820e-10	8.749560e-10	unchanged	both	1178.923170	985.072503	7088491.0	7089271.0	1326.0	1147.000000	promoter	True
4	ENSMUSG00000103509	Gm38372	chr1	7148109	7152137	+	7146109	7148609	7148609	7152137	7148109	Gm38372	477.814110	-0.297845	0.117813	-2.528109	1.146786e-02	3.126582e-02	unchanged	NaN	NaN	NaN	NaN	NaN	NaN	NaN	NaN	False
"""

merged_df.head()

#%% ########################## analyze_peak_methylation ############################################################
methylation_results = analyze_peak_methylation(
    merged_df=merged_df,
    medip_dir=PATHS['medip_dir'],
    cell_type=cell_type,
    genome_fasta=PATHS['genome_fasta']
)

"""
	chr	start	end	gene_id	binding_type	expression_status	log2FoldChange	methylation	cpg_count
0	chr1	4807559.0	4808103.0	ENSMUSG00000025903	both	downregulated	-0.605672	32.307881	85
1	chr1	4857465.0	4858372.0	ENSMUSG00000033813	both	downregulated	-0.536336	24.982538	95
2	chr1	5083039.0	5083536.0	ENSMUSG00000033793	both	unchanged	-0.417803	20.052000	53
3	chr1	7088491.0	7089271.0	ENSMUSG00000051285	both	unchanged	-0.447877	37.984170	85
5	chr1	9545269.0	9546268.0	ENSMUSG00000025911	both	upregulated	0.697750	100.000000	543
"""
methylation_results.head()

#%% ########################## plot_peak_methylation_by_category ############################################################
peak_stats = plot_peak_methylation_by_category(methylation_results, cell_type)

#%% ########################## analyze_cpg_methylation ############################################################
cpg_islands = load_cpg_islands(PATHS['cpg_islands_file'])

cpg_methylation_results = analyze_cpg_methylation(
    merged_df=merged_df,
    cpg_islands=cpg_islands,
    medip_dir=PATHS['medip_dir'],
    cell_type=cell_type,
    genome_fasta=PATHS['genome_fasta']
)
cpg_methylation_results.head()
#%% ########################## plot_cpg_methylation_by_category ############################################################
cpg_stats = plot_cpg_methylation_by_category(cpg_methylation_results, cell_type)


#%% ########################## Documentation ############################################################

"""
This script analyzes methylation patterns in MeCP2 binding regions and CpG islands, comparing between different cell types (NSC vs NEU).

Input files:
- MeDIP-seq bigWig files for IP and input samples (specified in medip_dir)
- CpG islands BED file (specified in PATHS['cpg_islands_file']) 
- Genome FASTA file (specified in PATHS['genome_fasta'])
- Gene expression and MeCP2 binding data (in merged_df DataFrame)

The script performs two main analyses:

1. Peak-level methylation analysis:
- Calculates methylation levels in MeCP2 binding peaks using MeDIP-seq data
- Compares methylation between different binding categories (exo_only, both, etc.)
- Analyzes relationship between peak methylation and gene expression changes
- Generates plots showing methylation levels in promoters vs gene bodies

2. CpG island methylation analysis:  
- Identifies CpG islands overlapping with gene promoters and bodies
- Calculates methylation levels in these CpG islands
- Compares methylation patterns between binding categories
- Analyzes how CpG island methylation relates to gene expression
- Creates visualization plots

Output files:
- PDF plots in plots/peak_methylation/<cell_type>/ showing peak methylation analysis
- PDF plots in plots/cpg_methylation/<cell_type>/ showing CpG island methylation analysis
- DataFrames with detailed methylation results for peaks and CpG islands

The analysis helps understand how MeCP2 binding and DNA methylation patterns
are related to gene regulation in different cellular contexts.
"""
