# %%
import pandas as pd
import glob
import numpy as np
import os
import matplotlib.pyplot as plt
import sys

# Replace the hard-coded EXPERIMENT with command line argument
EXPERIMENT = sys.argv[1] if len(sys.argv) > 1 else "align2_005"

# Set working directory
os.chdir(f'/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment_{EXPERIMENT}')

# %% [markdown]
# # Parallel
# 

# %%
# Load differential expression analysis results
dea_nsc = pd.read_csv('../../DATA/DEA_NSC.csv')
print(dea_nsc.shape)
dea_nsc.head()

# %%
dea_nsc = dea_nsc[dea_nsc['padj'] < 0.05]
dea_nsc.shape

# %%
expression_threshold = dea_nsc['baseMean'].quantile(0.02)
print(expression_threshold)

# %%
dea_nsc = dea_nsc[dea_nsc['baseMean'] > expression_threshold]
dea_nsc.shape

# %%
dea_nsc.to_csv(f'/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/DATA/DEA_NSC_filtered.csv', index=False)

# %%
# Read CpG islands bed file with tab separator and proper column names
cpg_islands = pd.read_csv('../../DATA/cpg_islands.bed', sep='\t', 
                         names=['chr', 'start', 'end', 'id', 'cpg_label', 'cpg_count'])

# Remove the "CpG:" prefix from cpg_label column
cpg_islands['cpg_label'] = cpg_islands['cpg_label'].str.replace('CpG:', '')

print(cpg_islands.shape)
cpg_islands.head()

# %%
# Get list of all chunk files
chunk_files = glob.glob('mecp2_cpg_enrichment_parallel/chunk_*.csv')
chunk_files

# %%
# Read and concatenate all chunks
df_parallel = pd.concat([pd.read_csv(f) for f in chunk_files], ignore_index=True)

# Sort by chromosome and start position
df_parallel = df_parallel.sort_values(['chr', 'start'])

print(f"Total regions analyzed: {len(df_parallel)}")

# %%
df_parallel.head()

# %%
df_parallel = df_parallel[df_parallel['chr'].isin([f'chr{i}' for i in range(1,20)] + ['chrX', 'chrY'])]
df_parallel.shape

# %%
df_parallel['significant'] = True

# %%
df_parallel.to_csv('mecp2_cpg_enrichment_parallel/mecp2_cpg_enrichment_parallel.csv', index=False)

# %%
# df_parallel = df_parallel[(df_parallel['exo_signal'] > 4.0) | (df_parallel['endo_signal'] > 4.0)]

# %%
df_both = df_parallel[df_parallel['binding_type'] == "both"]
df_both.shape

# %%
df_both = df_both.sort_values('enrichment', ascending=False)
df_both.head()

# %%
df_both.head()

# %%
df_both.to_csv('mecp2_cpg_enrichment_parallel/mecp2_cpg_enrichment_both.csv', index=False)

# %% [markdown]
# # Integrate

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# %%
# Load the integrated results
results_dir = "integrated/"
df = pd.read_csv(f"{results_dir}/mecp2_enriched_genes.csv")

# %%
print(df.shape)
df.head()

# %%
df = df[df['baseMean'].isna() == False]
print(df.shape)
df.head()

# %%
df.binding_type.value_counts()

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
df['category'].value_counts()

# %%
# Create the density plot
plt.figure(figsize=(12, 8))

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


