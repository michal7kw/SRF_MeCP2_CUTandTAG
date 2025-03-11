#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn4
import numpy as np
from upsetplot import UpSet
from upsetplot import from_memberships
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
OUTPUT_DIRS = [
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

def process_experiment(neu_file, nsc_file, output_dir):
    # Read the CSV files
    print(f"Reading Neurons (Neu) CpG data from {neu_file}...")
    neu_df = pd.read_csv(neu_file)
    print(f"Reading Neural Stem Cells (NSC) CpG data from {nsc_file}...")
    nsc_df = pd.read_csv(nsc_file)

    # Print total counts
    print(f"Total CpG islands in Neurons: {len(neu_df)}")
    print(f"Total CpG islands in Neural Stem Cells: {len(nsc_df)}")

    # Apply categorization
    neu_df['binding_category'] = neu_df.apply(categorize_binding, axis=1)
    nsc_df['binding_category'] = nsc_df.apply(categorize_binding, axis=1)

    # Create unique identifiers for each CpG island
    neu_df['cpg_id'] = neu_df['chr'] + '_' + neu_df['start'].astype(str) + '_' + neu_df['end'].astype(str)
    nsc_df['cpg_id'] = nsc_df['chr'] + '_' + nsc_df['start'].astype(str) + '_' + nsc_df['end'].astype(str)

    # Create sets for each condition
    neu_endo = set(neu_df[neu_df['binding_category'].isin(['both', 'endo_only'])]['cpg_id'])
    neu_exo = set(neu_df[neu_df['binding_category'].isin(['both', 'exo_only'])]['cpg_id'])
    nsc_endo = set(nsc_df[nsc_df['binding_category'].isin(['both', 'endo_only'])]['cpg_id'])
    nsc_exo = set(nsc_df[nsc_df['binding_category'].isin(['both', 'exo_only'])]['cpg_id'])

    # Print counts for each condition
    print(f"Neurons Endo: {len(neu_endo)}")
    print(f"Neurons Exo: {len(neu_exo)}")
    print(f"NSC Endo: {len(nsc_endo)}")
    print(f"NSC Exo: {len(nsc_exo)}")

    # Create a dictionary to store the membership of each CpG island
    all_cpgs = neu_endo.union(neu_exo, nsc_endo, nsc_exo)
    memberships = []

    for cpg in all_cpgs:
        membership = []
        if cpg in neu_endo:
            membership.append("Neu Endo")
        if cpg in neu_exo:
            membership.append("Neu Exo")
        if cpg in nsc_endo:
            membership.append("NSC Endo")
        if cpg in nsc_exo:
            membership.append("NSC Exo")
        memberships.append(membership)

    # Create UpSet plot
    data = from_memberships(memberships, data=np.ones(len(memberships)))

    fig = plt.figure(figsize=(12, 8))
    upset = UpSet(data, sort_by='cardinality', show_counts=True)
    upset.plot()
    plt.title('CpG Islands Overlap Across Four Conditions', fontsize=16)
    upset_plot_path = os.path.join(output_dir, 'cpg_upset_plot.png')
    plt.savefig(upset_plot_path, dpi=300, bbox_inches='tight')
    print(f"UpSet plot saved to: {upset_plot_path}")
    plt.close()

    # Create a 4-way membership matrix for all CpGs
    membership_matrix = pd.DataFrame({
        'Neu_Endo': [1 if cpg in neu_endo else 0 for cpg in all_cpgs],
        'Neu_Exo': [1 if cpg in neu_exo else 0 for cpg in all_cpgs],
        'NSC_Endo': [1 if cpg in nsc_endo else 0 for cpg in all_cpgs],
        'NSC_Exo': [1 if cpg in nsc_exo else 0 for cpg in all_cpgs]
    })

    # Count CpGs in each combination
    combinations = {}
    for i in range(1, 16):  # 2^4 - 1 possible combinations (excluding empty set)
        binary = format(i, '04b')
        mask = (membership_matrix['Neu_Endo'] == int(binary[0])) & \
               (membership_matrix['Neu_Exo'] == int(binary[1])) & \
               (membership_matrix['NSC_Endo'] == int(binary[2])) & \
               (membership_matrix['NSC_Exo'] == int(binary[3]))
        count = mask.sum()
        if count > 0:
            label = []
            if binary[0] == '1': label.append('Neu_Endo')
            if binary[1] == '1': label.append('Neu_Exo')
            if binary[2] == '1': label.append('NSC_Endo')
            if binary[3] == '1': label.append('NSC_Exo')
            combinations['+'.join(label)] = count

    # Calculate the number of CpGs common to all four conditions
    common_to_all = len(neu_endo.intersection(neu_exo, nsc_endo, nsc_exo))

    # Save the results to a text file
    summary_path = os.path.join(output_dir, 'cpg_overlap_summary.txt')
    with open(summary_path, 'w') as f:
        f.write(f"Total CpG islands in Neurons: {len(neu_df)}\n")
        f.write(f"Total CpG islands in Neural Stem Cells: {len(nsc_df)}\n\n")
        f.write(f"Neurons Endo: {len(neu_endo)}\n")
        f.write(f"Neurons Exo: {len(neu_exo)}\n")
        f.write(f"NSC Endo: {len(nsc_endo)}\n")
        f.write(f"NSC Exo: {len(nsc_exo)}\n\n")
        f.write("CpG Island Counts by Condition Combination:\n")
        for combo, count in sorted(combinations.items(), key=lambda x: x[1], reverse=True):
            f.write(f"{combo}: {count}\n")
        f.write(f"\nCpG Islands Common to All Four Conditions: {common_to_all}\n")

    print(f"Summary saved to: {summary_path}")

    # Create Venn diagram
    fig, ax = plt.subplots(figsize=(10, 10))
    venn4([neu_endo, neu_exo, nsc_endo, nsc_exo], 
          set_labels=('Neu Endo', 'Neu Exo', 'NSC Endo', 'NSC Exo'))
    plt.title('CpG Islands Overlap Across Four Conditions', fontsize=16)
    venn_plot_path = os.path.join(output_dir, 'cpg_venn_diagram.png')
    plt.savefig(venn_plot_path, dpi=300, bbox_inches='tight')
    print(f"Venn diagram saved to: {venn_plot_path}")
    plt.close()

# Process each experiment
for neu_file, nsc_file, output_dir in zip(NEU_EXPERIMENTS_DIRS, NSC_EXPERIMENTS_DIRS, OUTPUT_DIRS):
    print(f"\nProcessing experiment:")
    print(f"Neu file: {neu_file}")
    print(f"NSC file: {nsc_file}")
    print(f"Output directory: {output_dir}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    process_experiment(neu_file, nsc_file, output_dir)

print("All experiments processed.")
