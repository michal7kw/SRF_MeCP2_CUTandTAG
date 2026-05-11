#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import numpy as np
import os
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from matplotlib.patches import Patch

# Define file paths
BASE_DIR = '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/cpg_enrichment'

NSC_EXPERIMENTS_DIRS = [
    # f'{BASE_DIR}/NSC/broad/cpg_enrichment_1_rep_in_cpg/cpg_enrichment_parallel.csv',
    # f'{BASE_DIR}/NSC/broad/cpg_enrichment_2_rep_in_cpg/cpg_enrichment_parallel.csv',
    # f'{BASE_DIR}/NSC/broad/cpg_enrichment_2_rep_in_peaks/cpg_enrichment_parallel.csv',
    f'{BASE_DIR}/NSC/broad/cpg_enrichment_1_rep_in_peaks/cpg_enrichment_parallel.csv'
]

NEU_EXPERIMENTS_DIRS = [
    # f'{BASE_DIR}/Neu/broad/cpg_enrichment_1_rep_in_cpg/cpg_enrichment_parallel.csv',
    # f'{BASE_DIR}/Neu/broad/cpg_enrichment_2_rep_in_cpg/cpg_enrichment_parallel.csv',
    # f'{BASE_DIR}/Neu/broad/cpg_enrichment_2_rep_in_peaks/cpg_enrichment_parallel.csv',
    f'{BASE_DIR}/Neu/broad/cpg_enrichment_1_rep_in_peaks/cpg_enrichment_parallel.csv'
]

ROOT_DIR = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
OUTPUT_DIR = [
    # f'{ROOT_DIR}/paper_finals/outputs/cpg_venn_diagrams/cpg_enrichment_1_rep_in_cpg',
    # f'{ROOT_DIR}/paper_finals/outputs/cpg_venn_diagrams/cpg_enrichment_2_rep_in_cpg',
    # f'{ROOT_DIR}/paper_finals/outputs/cpg_venn_diagrams/cpg_enrichment_2_rep_in_peaks',
    f'{ROOT_DIR}/paper_finals/outputs/cpg_venn_diagrams/cpg_enrichment_1_rep_in_peaks'
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

    # Custom function to customize venn2 diagrams
    def customize_venn2(venn_obj):
        # First, set all patches to have black borders
        for patch in venn_obj.patches:
            patch.set_edgecolor('black')
            patch.set_linewidth(2.0)
            patch.set_alpha(1.0)
        
        # Set unique regions to red
        if venn_obj.get_patch_by_id('10'):
            venn_obj.get_patch_by_id('10').set_facecolor('#FF0000')  # Red for unique to first set
            venn_obj.get_patch_by_id('10').set_alpha(0.7)
        
        if venn_obj.get_patch_by_id('01'):
            venn_obj.get_patch_by_id('01').set_facecolor('#FF0000')  # Red for unique to second set
            venn_obj.get_patch_by_id('01').set_alpha(0.7)
        
        # Set overlap region to transparent
        if venn_obj.get_patch_by_id('11'):
            venn_obj.get_patch_by_id('11').set_facecolor('none')
            venn_obj.get_patch_by_id('11').set_alpha(0.0)
            
        # Remove the default labels
        for label in venn_obj.set_labels:
            if label is not None:
                label.set_text('')
        for label in venn_obj.subset_labels:
            if label is not None:
                label.set_text('')
    
    # Custom function to customize venn3 diagrams
    def customize_venn3(venn_obj):
        # First, set all patches to have black borders
        for patch in venn_obj.patches:
            patch.set_edgecolor('black')
            patch.set_linewidth(2.0)
            patch.set_alpha(1.0)
        
        # Set unique regions to red
        for region in ['100', '010', '001']:
            if venn_obj.get_patch_by_id(region):
                venn_obj.get_patch_by_id(region).set_facecolor('#FF0000')  # Red for unique regions
                venn_obj.get_patch_by_id(region).set_alpha(0.7)
        
        # Set double overlap regions to lighter red (gradient)
        for i, region in enumerate(['110', '101', '011']):
            if venn_obj.get_patch_by_id(region):
                # Create a gradient from red to white (75% of the way to white)
                color = to_rgba('#FF8080')
                venn_obj.get_patch_by_id(region).set_facecolor(color)
                venn_obj.get_patch_by_id(region).set_alpha(0.7)
        
        # Set triple overlap region to transparent
        if venn_obj.get_patch_by_id('111'):
            venn_obj.get_patch_by_id('111').set_facecolor('none')
            venn_obj.get_patch_by_id('111').set_alpha(0.0)
            
        # Remove the default labels
        for label in venn_obj.set_labels:
            if label is not None:
                label.set_text('')
        for label in venn_obj.subset_labels:
            if label is not None:
                label.set_text('')

    # 1. Compare NSC Endo vs NSC Exo
    plt.subplot(2, 2, 1)
    v1 = venn2([nsc_endo, nsc_exo], ('NSC Endo', 'NSC Exo'))
    customize_venn2(v1)
    
    # Create custom legend with counts
    nsc_endo_only = len(nsc_endo - nsc_exo)
    nsc_exo_only = len(nsc_exo - nsc_endo)
    nsc_both = len(nsc_endo.intersection(nsc_exo))
    
    legend_elements = [
        Patch(facecolor='#FF0000', edgecolor='black', alpha=0.7, label=f'NSC Endo only: {nsc_endo_only}'),
        Patch(facecolor='#FF0000', edgecolor='black', alpha=0.7, label=f'NSC Exo only: {nsc_exo_only}'),
        Patch(facecolor='none', edgecolor='black', label=f'Both: {nsc_both}')
    ]
    plt.legend(handles=legend_elements, loc='best', fontsize=16)
    plt.title('NSCs: Endo vs Exo', fontsize=18)

    # 2. Compare Neu Endo vs Neu Exo
    plt.subplot(2, 2, 2)
    v2 = venn2([neu_endo, neu_exo], ('Neu Endo', 'Neu Exo'))
    customize_venn2(v2)
    
    # Create custom legend with counts
    neu_endo_only = len(neu_endo - neu_exo)
    neu_exo_only = len(neu_exo - neu_endo)
    neu_both = len(neu_endo.intersection(neu_exo))
    
    legend_elements = [
        Patch(facecolor='#FF0000', edgecolor='black', alpha=0.7, label=f'Neu Endo only: {neu_endo_only}'),
        Patch(facecolor='#FF0000', edgecolor='black', alpha=0.7, label=f'Neu Exo only: {neu_exo_only}'),
        Patch(facecolor='none', edgecolor='black', label=f'Both: {neu_both}')
    ]
    plt.legend(handles=legend_elements, loc='best', fontsize=16)
    plt.title('Neurons: Endo vs Exo', fontsize=18)

    # 3. Compare NSC (Endo+Exo) vs Neu (Endo+Exo)
    plt.subplot(2, 2, 3)
    nsc_all = nsc_endo.union(nsc_exo)
    neu_all = neu_endo.union(neu_exo)
    v3 = venn2([nsc_all, neu_all], ('NSCs', 'Neurons'))
    customize_venn2(v3)
    
    # Create custom legend with counts
    nsc_only = len(nsc_all - neu_all)
    neu_only = len(neu_all - nsc_all)
    both_cells = len(nsc_all.intersection(neu_all))
    
    legend_elements = [
        Patch(facecolor='#FF0000', edgecolor='black', alpha=0.7, label=f'NSCs only: {nsc_only}'),
        Patch(facecolor='#FF0000', edgecolor='black', alpha=0.7, label=f'Neurons only: {neu_only}'),
        Patch(facecolor='none', edgecolor='black', label=f'Both: {both_cells}')
    ]
    plt.legend(handles=legend_elements, loc='best', fontsize=16)
    plt.title('All NSCs vs All Neurons', fontsize=18)

    # 4. Compare Endo (NSC+Neu) vs Exo (NSC+Neu)
    plt.subplot(2, 2, 4)
    all_endo = nsc_endo.union(neu_endo)
    all_exo = nsc_exo.union(neu_exo)
    v4 = venn2([all_endo, all_exo], ('Endo', 'Exo'))
    customize_venn2(v4)
    
    # Create custom legend with counts
    endo_only = len(all_endo - all_exo)
    exo_only = len(all_exo - all_endo)
    both_conditions = len(all_endo.intersection(all_exo))
    
    legend_elements = [
        Patch(facecolor='#FF0000', edgecolor='black', alpha=0.7, label=f'Endo only: {endo_only}'),
        Patch(facecolor='#FF0000', edgecolor='black', alpha=0.7, label=f'Exo only: {exo_only}'),
        Patch(facecolor='none', edgecolor='black', label=f'Both: {both_conditions}')
    ]
    plt.legend(handles=legend_elements, loc='best', fontsize=16)
    plt.title('All Endo vs All Exo', fontsize=18)

    plt.suptitle('CpG Islands Distribution Across Conditions', fontsize=20)
    plt.tight_layout()
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    plt.savefig(f'{output_dir}/cpg_venn_diagrams.png', dpi=300, bbox_inches='tight')
    print(f"Venn diagrams saved to: {output_dir}/cpg_venn_diagrams.png")

    # Create a 3-way Venn diagram for Cell Type and Condition
    plt.figure(figsize=(20, 8))

    # 1. NSC Endo, NSC Exo, Neu Endo
    plt.subplot(1, 2, 1)
    v5 = venn3([nsc_endo, nsc_exo, neu_endo], ('NSC Endo', 'NSC Exo', 'Neu Endo'))
    customize_venn3(v5)
    
    # Create custom legend with counts for 3-way Venn
    a = nsc_endo
    b = nsc_exo
    c = neu_endo
    
    # Calculate all regions
    only_a = len(a - b - c)
    only_b = len(b - a - c)
    only_c = len(c - a - b)
    a_and_b = len((a & b) - c)
    a_and_c = len((a & c) - b)
    b_and_c = len((b & c) - a)
    a_and_b_and_c = len(a & b & c)
    
    legend_elements = [
        Patch(facecolor='#FF0000', edgecolor='black', alpha=0.7, label=f'NSC Endo only: {only_a}'),
        Patch(facecolor='#FF0000', edgecolor='black', alpha=0.7, label=f'NSC Exo only: {only_b}'),
        Patch(facecolor='#FF0000', edgecolor='black', alpha=0.7, label=f'Neu Endo only: {only_c}'),
        Patch(facecolor='#FF8080', edgecolor='black', alpha=0.7, label=f'NSC Endo & NSC Exo: {a_and_b}'),
        Patch(facecolor='#FF8080', edgecolor='black', alpha=0.7, label=f'NSC Endo & Neu Endo: {a_and_c}'),
        Patch(facecolor='#FF8080', edgecolor='black', alpha=0.7, label=f'NSC Exo & Neu Endo: {b_and_c}'),
        Patch(facecolor='none', edgecolor='black', label=f'All three: {a_and_b_and_c}')
    ]
    plt.legend(handles=legend_elements, loc='best', fontsize=14)
    plt.title('NSC Endo vs NSC Exo vs Neu Endo', fontsize=18)

    # 2. NSC Endo, NSC Exo, Neu Exo
    plt.subplot(1, 2, 2)
    v6 = venn3([nsc_endo, nsc_exo, neu_exo], ('NSC Endo', 'NSC Exo', 'Neu Exo'))
    customize_venn3(v6)
    
    # Create custom legend with counts for 3-way Venn
    a = nsc_endo
    b = nsc_exo
    c = neu_exo
    
    # Calculate all regions
    only_a = len(a - b - c)
    only_b = len(b - a - c)
    only_c = len(c - a - b)
    a_and_b = len((a & b) - c)
    a_and_c = len((a & c) - b)
    b_and_c = len((b & c) - a)
    a_and_b_and_c = len(a & b & c)
    
    legend_elements = [
        Patch(facecolor='#FF0000', edgecolor='black', alpha=0.7, label=f'NSC Endo only: {only_a}'),
        Patch(facecolor='#FF0000', edgecolor='black', alpha=0.7, label=f'NSC Exo only: {only_b}'),
        Patch(facecolor='#FF0000', edgecolor='black', alpha=0.7, label=f'Neu Exo only: {only_c}'),
        Patch(facecolor='#FF8080', edgecolor='black', alpha=0.7, label=f'NSC Endo & NSC Exo: {a_and_b}'),
        Patch(facecolor='#FF8080', edgecolor='black', alpha=0.7, label=f'NSC Endo & Neu Exo: {a_and_c}'),
        Patch(facecolor='#FF8080', edgecolor='black', alpha=0.7, label=f'NSC Exo & Neu Exo: {b_and_c}'),
        Patch(facecolor='none', edgecolor='black', label=f'All three: {a_and_b_and_c}')
    ]
    plt.legend(handles=legend_elements, loc='best', fontsize=14)
    plt.title('NSC Endo vs NSC Exo vs Neu Exo', fontsize=18)

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
