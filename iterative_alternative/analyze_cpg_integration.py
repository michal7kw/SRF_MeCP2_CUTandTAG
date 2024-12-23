# %%
import pandas as pd
import glob
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Get required command line arguments
if len(sys.argv) <= 2:
    raise ValueError("Both EXPERIMENT and CELL_LINE must be provided as command line arguments")

EXPERIMENT = sys.argv[1]
CELL_LINE = sys.argv[2]

# Set working directory
os.chdir(f'/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment_{EXPERIMENT}/{CELL_LINE}')

# %%
# Load the integrated results
results_dir = "integrated/"
df = pd.read_csv(f"{results_dir}/mecp2_enriched_genes.csv")

# %%
df = df[df['baseMean'].isna() == False]

# %%
df = df[df['binding_type'] == 'both']

# %%
# Define gene categories based on RNA-seq data
def categorize_gene(row, fc_threshold=0.5, padj_threshold=0.05):
    if pd.isna(row['log2FoldChange']) or pd.isna(row['padj']):
        return 'non-deregulated'
    elif row['padj'] > padj_threshold:
        return 'non-deregulated'
    elif row['log2FoldChange'] >= fc_threshold:
        return 'up-regulated'
    elif row['log2FoldChange'] <= -fc_threshold:
        return 'down-regulated'
    else:
        return 'non-deregulated'

# %%
# Add category column
df['category'] = df.apply(categorize_gene, axis=1)

# %%
# Plot density for each category
for category, color in zip(['non-deregulated', 'up-regulated', 'down-regulated'], 
                         ['blue', 'orange', 'green']):
    subset = df[df['category'] == category]
    if len(subset) > 0:
        sns.kdeplot(data=subset['mecp2_enrichment'], 
                   label=category,
                   color=color)

plt.title('Mecp2 Enrichment Distribution by Gene Category')
plt.xlabel('Enrichment (Exo/Endo)')
plt.ylabel('Density')
plt.xlim(0, 8)  
plt.ylim(0, 1.2)
plt.legend()

# Add some statistics
for category in ['non-deregulated', 'up-regulated', 'down-regulated']:
    subset = df[df['category'] == category]
    print(f"\n{category}:")
    print(f"Number of genes: {len(subset)}")
    print(f"Mean enrichment: {subset['mecp2_enrichment'].mean():.2f}")
    print(f"Median enrichment: {subset['mecp2_enrichment'].median():.2f}")

# Save the plot
plt.savefig(f"{results_dir}/mecp2_enrichment_by_expression.pdf", 
            bbox_inches='tight', dpi=300)
plt.close()

# Export gene lists by category
for category in ['non-deregulated', 'up-regulated', 'down-regulated']:
    subset = df[df['category'] == category]
    # Save each category to a separate CSV file
    subset.to_csv(f"{results_dir}/{category.replace('-', '_')}_genes.csv", index=False)


