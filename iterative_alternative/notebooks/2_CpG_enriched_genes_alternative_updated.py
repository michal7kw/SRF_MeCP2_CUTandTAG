# %% [markdown]
# # Documentation

# %% [markdown]
# ## CpG-Gene Enrichment Analysis Pipeline
# ## Overview
# Analyze the relationship between genomic *peaks*, *CpG islands*, and gene *transcription start sites* (TSS) in different cell types and conditions. 
# 
# >Identify genes with CpG-overlapping peaks near their TSS regions and perform comprehensive comparative analysis.
# 
# Purposes:
# 1. Identify genomic peaks that overlap with CpG islands
# 2. Analyze the relationship between CpG-associated peaks and gene TSS regions
# 3. Compare differences in CpG-gene associations between:
#    - Different cell types (NSC and Neurons)
#    - Different conditions (Exogenous and Endogenous)
# 
# Features:
# - Filters for standard chromosomes (chr1-22, X, Y)
# - Extends CpG islands and TSS regions by customizable distances
# - Calculates coverage statistics for peak-CpG overlaps
# - Identifies genes with CpG-overlapping peaks near their TSS
# - Generates comparative analysis between different conditions
# 
# Input:
# - Peak files in BED format (`results/consensus_peaks/{cell_type}_{condition}_consensus.bed`)
# - CpG islands file in BED format (`DATA/cpg_islands.bed`)
# - Gene annotation file in GTF format (`DATA/gencode.vM10.annotation.gtf`)
# - Genome size file (`DATA/genome.size`)
# 
# Key Parameters:
# - `extend_cpg`: Base pairs to extend CpG islands (default: 300bp)
# - `extend_tss`: Base pairs to extend TSS regions (default: 2000bp)
# - `coverage_threshold`: Minimum percentage overlap required (adjustable, examples use 20% and 80%)
# 
# Output Files:
# The pipeline generates several output files across multiple directories:
# 
# ### Main Analysis Output (`results/cpg_gene_analysis/`)
# - `{cell_type}_{condition}_cpg_genes.tsv`: Detailed gene-level statistics
# - `{cell_type}_common_genes.txt`: Genes common between conditions
# - `comparison_summary.txt`: Overall comparison statistics
# 
# ### Visualization Output (`results/visualizations/`)
# - `peak_distribution.png`: Violin plots of peaks per gene
# - `distance_distribution.png`: KDE plots of TSS distances
# - `{cell_type}_venn.png`: Venn diagrams of condition overlaps
# - `{cell_type}_{condition}_top_genes_heatmap.png`: Heatmaps of top genes
# - `total_genes_summary.png`: Bar plots of total genes
# - `peak_size_distribution.png`: Box plots of peak sizes
# 
# ### Comprehensive Analysis Output (`results/comprehensive_analysis/`)
# - `endo_cpg_coverage_comparison.png`: Comparison of Endo coverage between cell types
# - `cpg_target_overlap.png`: Four-way Venn diagram of all conditions
# - `endo_coverage_stats.txt`: Statistical analysis of coverage patterns
# - `cpg_overlap_stats.txt`: Detailed overlap statistics
# - `{cell_type}_exo_enriched_genes.tsv`: Exo vs Endo enrichment analysis
# - `common_enriched_genes.txt`: Genes enriched in both cell types
# 
# ## Main Functions
# 
# ### Data Processing Functions
# 1. **get_peaks_with_cpg()**
#    - *Identifies peaks overlapping with CpG islands*
#    - Calculates coverage statistics
#    - Filters based on coverage threshold
# 
# 2. **get_genes_with_cpg_enrichment()**
#    - *Maps CpG-overlapping peaks to nearby genes*
#    - Calculates various metrics including:
#      - Number of peaks per gene
#      - Total coverage
#      - Average peak size
#      - Minimum distance to TSS
# 
# 3. **extract_tss_regions()**
#    - Extracts TSS locations from GTF file
#    - Filters for protein-coding genes
#    - Handles both forward and reverse strands
# 
# ### Analysis Functions
# 4. **compare_endo_cpg_coverage()**
#    - *Compares CpG coverage between NPCs and Neurons*
#    - Performs statistical testing
#    - Generates coverage distribution plots
# 
# 5. **analyze_common_cpg_targets()**
#    - *Identifies targets common across conditions*
#    - Creates four-way Venn diagrams
#    - Calculates pairwise overlaps
# 
# 6. **analyze_exo_vs_endo_enrichment()**
#    - Calculates enrichment ratios
#    - Identifies significantly enriched genes
#    - Finds common enriched genes across cell types
# 
# ### Visualization Functions
# 7. **plot_peak_distribution()**
#    - Creates violin plots of peak distributions
#    - Compares across all conditions
# 
# 8. **plot_distance_distribution()**
#    - Generates KDE plots of TSS distances
#    - Shows distribution patterns across conditions
# 
# 9. **create_venn_diagrams()**
#    - Visualizes overlap between conditions
#    - Separate diagrams for each cell type
# 
# 10. **plot_top_genes_heatmap()**
#     - Creates heatmaps of top genes
#     - Shows normalized coverage scores
# 
# ### Utility Functions
# 11. **load_and_process_data()**
#     - Loads analysis results
#     - Calculates additional metrics
#     - Prepares data for visualization
# 
# 12. **create_comparison_summary()**
#     - Generates comparative analysis
#     - Identifies condition-specific genes
#     - Calculates overlap statistics
# 
# ## Usage
# The pipeline can be run with different coverage thresholds:
# ```python
# # For 20% coverage threshold
# get_genes_with_cpg_enrichment(..., coverage_threshold=20)
# 
# # For 80% coverage threshold
# get_genes_with_cpg_enrichment(..., coverage_threshold=80)
# ```
# 
# Visualizations can be generated using:
# ```python
# generate_all_visualizations(results_dir, visualization_dir)
# generate_comprehensive_analysis(results_dir, analysis_dir)
# ```

# %% [markdown]
# ![image.png](attachment:image.png)

# %% [markdown]
# # Environment

# %%
# Standard library imports
import os
import subprocess
from collections import defaultdict
from pathlib import Path

# Third party imports
import numpy as np
import pandas as pd
from scipy import stats

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn3
import glob
from IPython.display import Image, display
# from venn import venn

wd_dir = '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/notebooks'
os.chdir(wd_dir)

root_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"

current_dir = os.getcwd()

# Import functions module using full path
import sys
sys.path.append("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/scripts/Archive")
import functions_CpG_enriched_genes
import importlib
importlib.reload(functions_CpG_enriched_genes)
from functions_CpG_enriched_genes import *


# %% [markdown]
# # Run analysis

# %% [markdown]
# ## 20% coverage threshold

# %%
# Setup paths
cpg_file = f"{root_dir}/DATA/cpg_islands.bed"
gtf_file = f"{root_dir}/DATA/gencode.vM10.annotation.gtf"
output_dir = f"{root_dir}/results/cpg_gene_analysis"
os.makedirs(output_dir, exist_ok=True)

results = {}

for cell_type in ['NSC', 'Neuron']:
    results[cell_type] = {}
    for condition in ['Exo', 'Endo']:
        peak_file = f"{root_dir}/results/consensus_peaks/{cell_type}_{condition}_consensus.bed"
        
        result = get_genes_with_cpg_enrichment(
            peak_file=peak_file,
            cpg_file=cpg_file,
            gtf_file=gtf_file,
            output_dir=output_dir,
            cell_type=cell_type,
            condition=condition,
            extend_cpg=0,
            extend_tss=2000,
            coverage_threshold=20,
            genome_size_file=f"{root_dir}/DATA/genome.size"
        )
        
        if result is not None:
            results[cell_type][condition] = result
        else:
            print(f"Skipping comparison for {cell_type} {condition} due to processing error")
            results[cell_type][condition] = {'genes': pd.DataFrame(), 'total_peaks': 0}

# %%
create_comparison_summary(results, output_dir)

# %% [markdown]
# # Visualization 1

# %%
results_dir = f"{root_dir}/results/cpg_gene_analysis"
data = load_and_process_data(results_dir)

# %%
# visualization_dir = "results/visualizations"
# generate_all_visualizations(results_dir)

# %%
# plot_peak_distribution(data)

# %%
# plot_distance_distribution(data)

# %%
import importlib
import functions_CpG_enriched_genes
importlib.reload(functions_CpG_enriched_genes)
from functions_CpG_enriched_genes import *

# %%
create_venn_diagrams_percentages(data)

# %%
import importlib
import functions_CpG_enriched_genes
importlib.reload(functions_CpG_enriched_genes)
from functions_CpG_enriched_genes import *

create_venn_diagrams_endo_comparison(data)

# %%
# plot_top_genes_heatmap(data)


# %%
create_summary_statistics(data)


# %%
# plot_peak_size_distribution(data)

# %%
# print("Displaying generated visualizations:")
# for img_path in glob.glob(f"{visualization_dir}/*.png"):
#     print(f"\n{os.path.basename(img_path)}:")
#     display(Image(filename=img_path, width=500))


# %%
# Create results directory if it doesn't exist
os.makedirs("results", exist_ok=True)
os.makedirs(f"{root_dir}/results/cpg_gene_analysis", exist_ok=True)

results_dir = f"{root_dir}/results/cpg_gene_analysis"
generate_comprehensive_analysis(results_dir)

# %% [markdown]
# ## NSCs vs Neurons

# %%
# compare_endo_cpg_coverage(data)

# %%
compare_endo_cpg_coverage_common(data)

# %% [markdown]
# ## Exo vs Endo

# %%
compare_exo_endo_coverage_common(data)


