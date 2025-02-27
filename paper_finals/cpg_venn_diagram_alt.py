#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import numpy as np
import os

# Define file paths
BASE_DIR = '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/cpg_enrichment'

NSC_EXPERIMENTS_DIRS = [
    f'{BASE_DIR}/NSC/broad/cpg_enrichment_1_rep_in_cpg/cpg_enrichment_parallel.csv',
    f'{BASE_DIR}/NSC/broad/cpg_enrichment_2_rep_in_cpg/cpg_enrichment_parallel.csv',
    f'{BASE_DIR}/NSC/broad/cpg_enrichment_1_rep_in_peaks/cpg_enrichment_parallel.csv',
    f'{BASE_DIR}/NSC/broad/cpg_enrichment_2_rep_in_peaks/cpg_enrichment_parallel.csv'
]

NEU_EXPERIMENTS_DIRS = [
    f'{BASE_DIR}/Neu/broad/cpg_enrichment_1_rep_in_cpg/cpg_enrichment_parallel.csv',
    f'{BASE_DIR}/Neu/broad/cpg_enrichment_2_rep_in_cpg/cpg_enrichment_parallel.csv',
    f'{BASE_DIR}/Neu/broad/cpg_enrichment_1_rep_in_peaks/cpg_enrichment_parallel.csv',
    f'{BASE_DIR}/Neu/broad/cpg_enrichment_2_rep_in_peaks/cpg_enrichment_parallel.csv'
]

ROOT_DIR = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
OUTPUT_DIR = [
    f'{ROOT_DIR}/paper_finals/outputs/cpg_venn_diagrams/cpg_enrichment_1_rep_in_peaks',
    f'{ROOT_DIR}/paper_finals/outputs/cpg_venn_diagrams/cpg_enrichment_2_rep_in_peaks',
    f'{ROOT_DIR}/paper_finals/outputs/cpg_venn_diagrams/cpg_enrichment_1_rep_in_cpg',
    f'{ROOT_DIR}/paper_finals/outputs/cpg_venn_diagrams/cpg_enrichment_2_rep_in_cpg'
]

def categorize_binding(row):
    binding_type = row['binding_type']
    if binding_type == 'both':
        return 'both'
    elif binding_type == 'exo_only':
        return 'exo_only'
    elif binding_type == 'endo_only':
        return 'endo_only'
    else:
        return 'none'

def process_experiment(nsc_file, neu_file, output_dir):
    # Read the CSV files
    print(f"Reading Neural Stem Cells (NSC) CpG data from {nsc_file}...")
    nsc_df = pd.read_csv(nsc_file)
    print(f"Reading Neurons (Neu) CpG data from {neu_file}...")
    neu_df = pd.read_csv(neu_file)

    # Print total counts
    print(f"Total CpG islands in Neural Stem Cells: {len(nsc_df)}")
    print(f"Total CpG islands in Neurons: {len(neu_df)}")

    # Apply categorization
    nsc_df['binding_category'] = nsc_df.apply(categorize_binding, axis=1)
    neu_df['binding_category'] = neu_df.apply(categorize_binding, axis=1)

    # Create unique identifiers for each CpG island
    nsc_df['cpg_id'] = nsc_df['chr'] + '_' + nsc_df['start'].astype(str) + '_' + nsc_df['end'].astype(str)
    neu_df['cpg_id'] = neu_df['chr'] + '_' + neu_df['start'].astype(str) + '_' + neu_df['end'].astype(str)

    # Create sets for each condition
    nsc_endo = set(nsc_df[nsc_df['binding_category'].isin(['both', 'endo_only'])]['cpg_id'])
    nsc_exo = set(nsc_df[nsc_df['binding_category'].isin(['both', 'exo_only'])]['cpg_id'])
    neu_endo = set(neu_df[neu_df['binding_category'].isin(['both', 'endo_only'])]['cpg_id'])
    neu_exo = set(neu_df[neu_df['binding_category'].isin(['both', 'exo_only'])]['cpg_id'])

    # Print counts for each condition
    print(f"NSC Endo: {len(nsc_endo)}")
    print(f"NSC Exo: {len(nsc_exo)}")
    print(f"Neurons Endo: {len(neu_endo)}")
    print(f"Neurons Exo: {len(neu_exo)}")

    # Create plots
    plt.figure(figsize=(20, 15))

    # 1. Compare NSC Endo vs NSC Exo
    plt.subplot(2, 2, 1)
    venn2([nsc_endo, nsc_exo], ('NSC Endo', 'NSC Exo'))
    plt.title('NSCs: Endo vs Exo', fontsize=14)

    # 2. Compare Neu Endo vs Neu Exo
    plt.subplot(2, 2, 2)
    venn2([neu_endo, neu_exo], ('Neu Endo', 'Neu Exo'))
    plt.title('Neurons: Endo vs Exo', fontsize=14)

    # 3. Compare NSC (Endo+Exo) vs Neu (Endo+Exo)
    plt.subplot(2, 2, 3)
    nsc_all = nsc_endo.union(nsc_exo)
    neu_all = neu_endo.union(neu_exo)
    venn2([nsc_all, neu_all], ('NSCs', 'Neurons'))
    plt.title('All NSCs vs All Neurons', fontsize=14)

    # 4. Compare Endo (NSC+Neu) vs Exo (NSC+Neu)
    plt.subplot(2, 2, 4)
    all_endo = nsc_endo.union(neu_endo)
    all_exo = nsc_exo.union(neu_exo)
    venn2([all_endo, all_exo], ('Endo', 'Exo'))
    plt.title('All Endo vs All Exo', fontsize=14)

    plt.suptitle('CpG Islands Distribution Across Conditions', fontsize=18)
    plt.tight_layout()
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    plt.savefig(f'{output_dir}/cpg_venn_diagrams.png', dpi=300, bbox_inches='tight')
    print(f"Venn diagrams saved to: {output_dir}/cpg_venn_diagrams.png")

    # Create a 3-way Venn diagram for Cell Type and Condition
    plt.figure(figsize=(20, 8))

    # 1. NSC Endo, NSC Exo, Neu Endo
    plt.subplot(1, 2, 1)
    venn3([nsc_endo, nsc_exo, neu_endo], ('NSC Endo', 'NSC Exo', 'Neu Endo'))
    plt.title('NSC Endo vs NSC Exo vs Neu Endo', fontsize=14)

    # 2. NSC Endo, NSC Exo, Neu Exo
    plt.subplot(1, 2, 2)
    venn3([nsc_endo, nsc_exo, neu_exo], ('NSC Endo', 'NSC Exo', 'Neu Exo'))
    plt.title('NSC Endo vs NSC Exo vs Neu Exo', fontsize=14)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/cpg_venn_diagrams_3way.png', dpi=300, bbox_inches='tight')
    print(f"3-way Venn diagrams saved to: {output_dir}/cpg_venn_diagrams_3way.png")

    # Calculate the number of CpGs common to all four conditions
    common_to_all = len(nsc_endo.intersection(nsc_exo, neu_endo, neu_exo))
    print(f"\nCpG Islands Common to All Four Conditions: {common_to_all}")

    # Save the results to a text file
    with open(f'{output_dir}/cpg_overlap_summary.txt', 'w') as f:
        f.write(f"Total CpG islands in Neural Stem Cells: {len(nsc_df)}\n")
        f.write(f"Total CpG islands in Neurons: {len(neu_df)}\n\n")
        f.write(f"NSC Endo: {len(nsc_endo)}\n")
        f.write(f"NSC Exo: {len(nsc_exo)}\n")
        f.write(f"Neurons Endo: {len(neu_endo)}\n")
        f.write(f"Neurons Exo: {len(neu_exo)}\n\n")
        f.write(f"CpG Islands Common to All Four Conditions: {common_to_all}\n\n")
        
        # Additional pairwise overlaps
        f.write("Pairwise Overlaps:\n")
        f.write(f"NSC Endo ∩ NSC Exo: {len(nsc_endo.intersection(nsc_exo))}\n")
        f.write(f"NSC Endo ∩ Neu Endo: {len(nsc_endo.intersection(neu_endo))}\n")
        f.write(f"NSC Endo ∩ Neu Exo: {len(nsc_endo.intersection(neu_exo))}\n")
        f.write(f"NSC Exo ∩ Neu Endo: {len(nsc_exo.intersection(neu_endo))}\n")
        f.write(f"NSC Exo ∩ Neu Exo: {len(nsc_exo.intersection(neu_exo))}\n")
        f.write(f"Neu Endo ∩ Neu Exo: {len(neu_endo.intersection(neu_exo))}\n\n")
        
        # Triple overlaps
        f.write("Triple Overlaps:\n")
        f.write(f"NSC Endo ∩ NSC Exo ∩ Neu Endo: {len(nsc_endo.intersection(nsc_exo, neu_endo))}\n")
        f.write(f"NSC Endo ∩ NSC Exo ∩ Neu Exo: {len(nsc_endo.intersection(nsc_exo, neu_exo))}\n")
        f.write(f"NSC Endo ∩ Neu Endo ∩ Neu Exo: {len(nsc_endo.intersection(neu_endo, neu_exo))}\n")
        f.write(f"NSC Exo ∩ Neu Endo ∩ Neu Exo: {len(nsc_exo.intersection(neu_endo, neu_exo))}\n")

    print(f"Detailed summary saved to: {output_dir}/cpg_overlap_summary.txt")

# Process each experiment
for nsc_file, neu_file, output_dir in zip(NSC_EXPERIMENTS_DIRS, NEU_EXPERIMENTS_DIRS, OUTPUT_DIR):
    print(f"\nProcessing experiment:")
    print(f"NSC file: {nsc_file}")
    print(f"Neu file: {neu_file}")
    print(f"Output directory: {output_dir}")
    process_experiment(nsc_file, neu_file, output_dir)
