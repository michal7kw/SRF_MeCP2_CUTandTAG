# %%
import pandas as pd
import glob
import numpy as np
import os
import matplotlib.pyplot as plt
import sys

# Get required command line arguments
if len(sys.argv) <= 2:
    raise ValueError("Both EXPERIMENT and CELL_LINE must be provided as command line arguments")

EXPERIMENT = sys.argv[1]
CELL_LINE = sys.argv[2]

DATA_DIR = f'/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/DATA'

# Set working directory
os.chdir(f'/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment_{EXPERIMENT}/{CELL_LINE}')

# %%
# Get list of all chunk files
chunk_files = glob.glob('mecp2_cpg_enrichment_parallel/chunk_*.csv')

# %%
# Read and concatenate all chunks
df_parallel = pd.concat([pd.read_csv(f) for f in chunk_files], ignore_index=True)

# Sort by chromosome and start position
df_parallel = df_parallel.sort_values(['chr', 'start'])

# %%
df_parallel = df_parallel[df_parallel['chr'].isin([f'chr{i}' for i in range(1,20)] + ['chrX', 'chrY'])]

# %%
df_parallel['significant'] = True

# %%
df_parallel.to_csv('mecp2_cpg_enrichment_parallel/mecp2_cpg_enrichment_parallel.csv', index=False)

# %%
df_both = df_parallel[df_parallel['binding_type'] == "both"]

# %%
df_both = df_both.sort_values('enrichment', ascending=False)

# %%
df_both.to_csv('mecp2_cpg_enrichment_parallel/mecp2_cpg_enrichment_both.csv', index=False)
