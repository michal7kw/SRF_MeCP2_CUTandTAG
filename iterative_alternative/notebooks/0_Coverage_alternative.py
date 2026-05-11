# %% [markdown]
# # Environment

# %%
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import pickle
import os
from pathlib import Path
import mygene

wd_dir = '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/notebooks'
os.chdir(wd_dir)
root_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative"
data_dir = f"/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA"
aligned_dir = f"{root_dir}/results/aligned"
peaks_dir = f"{root_dir}/results_alternative"
broad_peaks_dir = f"{root_dir}/results/no_dedup/peaks/broad"

current_dir = os.getcwd()

# Always reload functions
import importlib
import sys
sys.path.append('/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/scripts/Archive')
import functions_Coverage
importlib.reload(functions_Coverage)
from functions_Coverage import *

# %%
cpg_file = f"{data_dir}/cpg_islands.bed"

# %%
# Read cpg_file
cpg_df = pd.read_csv(cpg_file, sep='\t', header=None)
cpg_df.head()

# %%
# Directly assign endogenous samples from data
endo_samples_nsc = pd.DataFrame({
    'SampleID': ['NSCM1', 'NSCM2', 'NSCM3'],
    'Tissue': ['NSC', 'NSC', 'NSC'],
    'Factor': ['Endo', 'Endo', 'Endo'],
    'Condition': ['Endogenous', 'Endogenous', 'Endogenous'],
    'Replicate': [1, 2, 3],
    'bamReads': [
        f'{aligned_dir}/NSCM1.bed',
        f'{aligned_dir}/NSCM2.bed',
        f'{aligned_dir}/NSCM3.bed'
    ],
    'Peaks': [
        f'{broad_peaks_dir}/NSCM1_broad_peaks.broadPeak',
        f'{broad_peaks_dir}/NSCM2_broad_peaks.broadPeak',
        f'{broad_peaks_dir}/NSCM3_broad_peaks.broadPeak'
    ],
    'PeakCaller': ['broad', 'broad', 'broad']
})

endo_samples_neurons = pd.DataFrame({
    'SampleID': ['NeuM2', 'NeuM3'],
    'Tissue': ['Neuron', 'Neuron'],
    'Factor': ['Endo', 'Endo'],
    'Condition': ['Endogenous', 'Endogenous'],
    'Replicate': [2, 3],
    'bamReads': [
        f'{aligned_dir}/NeuM2.bed',
        f'{aligned_dir}/NeuM3.bed'
    ],
    'Peaks': [
        f'{broad_peaks_dir}/NeuM2_broad_peaks.broadPeak',
        f'{broad_peaks_dir}/NeuM3_broad_peaks.broadPeak'
    ],
    'PeakCaller': ['broad', 'broad']
})


# %%
# Create DataFrames for exogenous samples
exo_samples_nsc = pd.DataFrame({
    'SampleID': ['NSCv1', 'NSCv2', 'NSCv3'],
    'Tissue': ['NSC', 'NSC', 'NSC'],
    'Factor': ['Virus', 'Virus', 'Virus'],
    'Condition': ['Exogenous', 'Exogenous', 'Exogenous'],
    'Replicate': [1, 2, 3],
    'bamReads': [
        f'{aligned_dir}/NSCv1.bam',
        f'{aligned_dir}/NSCv2.bam',
        f'{aligned_dir}/NSCv3.bam'
    ],
    'Peaks': [
        f'{broad_peaks_dir}/NSCv1_broad_peaks.broadPeak',
        f'{broad_peaks_dir}/NSCv2_broad_peaks.broadPeak',
        f'{broad_peaks_dir}/NSCv3_broad_peaks.broadPeak'
    ],
    'PeakCaller': ['broad', 'broad', 'broad']
})

exo_samples_neurons = pd.DataFrame({
    'SampleID': ['NeuV1', 'NeuV2', 'NeuV3'],
    'Tissue': ['Neuron', 'Neuron', 'Neuron'],
    'Factor': ['Virus', 'Virus', 'Virus'],
    'Condition': ['Exogenous', 'Exogenous', 'Exogenous'],
    'Replicate': [1, 2, 3],
    'bamReads': [
        f'{aligned_dir}/NeuV1.bam',
        f'{aligned_dir}/NeuV2.bam',
        f'{aligned_dir}/NeuV3.bam'
    ],
    'Peaks': [
        f'{broad_peaks_dir}/NeuV1_broad_peaks.broadPeak',
        f'{broad_peaks_dir}/NeuV2_broad_peaks.broadPeak',
        f'{broad_peaks_dir}/NeuV3_broad_peaks.broadPeak'
    ],
    'PeakCaller': ['broad', 'broad', 'broad']
})

# %%
print(exo_samples_nsc.shape)
print(endo_samples_nsc.shape)

print(exo_samples_neurons.shape)
print(endo_samples_neurons.shape)

# %% [markdown]
# # All cpg islands

# %% [markdown]
# ## Key differences in the two approaches:
# 
# ### 1. First Implementation (`per CpG Coverage`):
# - Calculates coverage as: overlap_length / cpg_length * 100
# - The coverage is based on how much of the CpG island is covered by the peak 
# - If multiple CpG islands overlap with a peak, it keeps all coverage values
# 
# ### 2. Second Implementation (`per Peak Coverage`):
# - Also calculates coverage as: overlap_length / cpg_length * 100
# - However, it maintains a per-peak tracking system where it:
#   - Groups all overlaps by peak ID
#   - Takes the maximum coverage value for each peak
#   - Only keeps peaks that meet the coverage threshold
#   - Reports the coverage statistics only for qualified peaks
# 
# ---
# 
# ### Peak vs CpG Island Perspective
# - First implementation keeps all overlap values, even if multiple CpG islands overlap with the same peak
# - Second implementation tracks coverage per peak and only keeps the maximum coverage value for each peak
# 
# ---
# 
# ### Key Structural Differences 
# 1. Return Value Structure
#    - First implementation: Returns a simple list of coverage percentages
#    - Second implementation: Returns a detailed dictionary with peak-specific information including:
#      - Chromosome
#      - Start/end positions
#      - List of all CpG overlaps for each peak
#      - Maximum coverage
#      - Total number of overlaps
# 
# ### Peak Identity Tracking
#    - First implementation: Doesn't track individual peaks or their locations
#    - Second implementation: Creates unique peak IDs using format "chromosome:start-end"
#    - Second implementation maintains peak context throughout processing
# 
# --- 
# 
# ### Errata
# 
# 1. Coverage Calculation
#    - First (Old) implementation: Calculates coverage using peak dimensions
#    - First (New) implementation: Identifies CpG regions by looking for 'CpG:' marker

# %% [markdown]
# ## NSC

# %% [markdown]
# ### Per CpG Coverage

# %%
print(exo_samples_nsc['Peaks'].to_string())


# %%
# %%capture
exo_coverage_nsc = []

# First filter the CpG file once (no need to do it for each peak file)
temp_cpg = "temp_filtered_cpg.bed"
filter_standard_chromosomes(input_file=cpg_file, output_file=temp_cpg)

try:
    # Process each peak file
    for peak_file in exo_samples_nsc['Peaks']:
        # Create temporary filtered peak file
        temp_peak = "temp_filtered_peaks.bed"
        filter_standard_chromosomes(peak_file, temp_peak)
        
        # Calculate coverage using filtered peaks and filtered CpG file
        coverage = calculate_peak_cpg_coverage(temp_peak, temp_cpg, extend=100, genome_size_file=f"{data_dir}/genome.size")
        exo_coverage_nsc.extend(coverage)
        
        # Clean up peak file
        if os.path.exists(temp_peak):
            os.remove(temp_peak)
finally:
    # Clean up CpG file
    if os.path.exists(temp_cpg):
        os.remove(temp_cpg)

# %%
exo_coverage_nsc[:10]

# %%
# %%capture
endo_coverage_nsc = []

# First filter the CpG file once (no need to do it for each peak file)
temp_cpg = "temp_filtered_cpg.bed"
filter_standard_chromosomes(input_file=cpg_file, output_file=temp_cpg)

try:
    # Process each peak file
    for peak_file in endo_samples_nsc['Peaks']:
        # Create temporary filtered peak file
        temp_peak = "temp_filtered_peaks.bed"
        filter_standard_chromosomes(peak_file, temp_peak)
        
        # Calculate coverage using filtered peaks and filtered CpG file
        coverage = calculate_peak_cpg_coverage(temp_peak, temp_cpg, extend=100, genome_size_file=f"{data_dir}/genome.size")
        endo_coverage_nsc.extend(coverage)
        
        # Clean up peak file
        if os.path.exists(temp_peak):
            os.remove(temp_peak)
finally:
    # Clean up CpG file
    if os.path.exists(temp_cpg):
        os.remove(temp_cpg)

# %%
endo_coverage_nsc[:10]


# %%
# Create the visualization
plot_coverage_histograms(exo_coverage_nsc, endo_coverage_nsc, min_coverage=0, n_bins=30)

# %%
plot_coverage_histograms_overlayed(exo_coverage_nsc, endo_coverage_nsc, min_coverage=0, n_bins=30)

# %%
os.makedirs(f'{peaks_dir}/coverage', exist_ok=True)
coverage_dir = f'{peaks_dir}/coverage'

# %%
# Save the lists to files
with open(f'{coverage_dir}/exo_coverage_nsc.pkl', 'wb') as f:
    pickle.dump(exo_coverage_nsc, f)

with open(f'{coverage_dir}/endo_coverage_nsc.pkl', 'wb') as f:
    pickle.dump(endo_coverage_nsc, f)


# %%
# Load the lists back
with open(f'{coverage_dir}/exo_coverage_nsc.pkl', 'rb') as f:
    exo_coverage_nsc = pickle.load(f)

with open(f'{coverage_dir}/endo_coverage_nsc.pkl', 'rb') as f:
    endo_coverage_nsc = pickle.load(f)

# %% [markdown]
# ### Per Peak Coverage

# %%
%%capture
exo_coverage_nsc = []
for peak_file in exo_samples_nsc['Peaks']:
    try:
        # Create temporary filtered peak file
        temp_peak = "temp_filtered_peaks.bed"
        filter_standard_chromosomes(peak_file, temp_peak)
        
        # Calculate coverage using filtered peaks
        peak_coverages = calculate_peak_cpg_coverage_per_peak(temp_peak, cpg_file, extend=100)
        
        # Extract coverage values, using 0.0 for peaks with no CpG overlap
        coverage_values = [info['max_coverage'] for info in peak_coverages.values()]
        exo_coverage_nsc.extend(coverage_values)
        
        # Clean up
        subprocess.run("rm temp_filtered_peaks.bed", shell=True)
        
    except Exception as e:
        print(f"Error processing {peak_file}: {str(e)}")
        continue

# %%
exo_coverage_nsc[:10]

# %%
%%capture
endo_coverage_nsc = []
for peak_file in endo_samples_nsc['Peaks']:
    try:
        # Create temporary filtered peak file
        temp_peak = "temp_filtered_peaks.bed"
        filter_standard_chromosomes(peak_file, temp_peak)
        
        # Calculate coverage using filtered peaks
        peak_coverages = calculate_peak_cpg_coverage_per_peak(temp_peak, cpg_file, extend=100)
        
        # Extract coverage values, using 0.0 for peaks with no CpG overlap
        coverage_values = [info['max_coverage'] for info in peak_coverages.values()]
        endo_coverage_nsc.extend(coverage_values)
        
        # Clean up
        subprocess.run("rm temp_filtered_peaks.bed", shell=True)
    except Exception as e:
        print(f"Error processing {peak_file}: {str(e)}")
        continue

# %%
endo_coverage_nsc[:10]

# %%
# Create the visualization
plot_coverage_histograms(exo_coverage_nsc, endo_coverage_nsc, min_coverage=0, n_bins=30)

# %%
plot_coverage_histograms_overlayed(exo_coverage_nsc, endo_coverage_nsc, min_coverage=0, n_bins=30)

# %% [markdown]
# ## Neurons
# 

# %% [markdown]
# ### Per CpG Coverage

# %%
%%capture
exo_coverage_neurons = []

# First filter the CpG file once (no need to do it for each peak file)
temp_cpg = "temp_filtered_cpg.bed"
filter_standard_chromosomes(input_file=cpg_file, output_file=temp_cpg)

try:
    # Process each peak file
    for peak_file in exo_samples_neurons['Peaks']:
        # Create temporary filtered peak file
        temp_peak = "temp_filtered_peaks.bed"
        filter_standard_chromosomes(peak_file, temp_peak)
        
        # Calculate coverage using filtered peaks and filtered CpG file
        coverage = calculate_peak_cpg_coverage(temp_peak, temp_cpg, extend=100, genome_size_file=f"{data_dir}/genome.size")
        exo_coverage_neurons.extend(coverage)
        
        # Clean up peak file
        if os.path.exists(temp_peak):
            os.remove(temp_peak)
finally:
    # Clean up CpG file
    if os.path.exists(temp_cpg):
        os.remove(temp_cpg)

# %%
exo_coverage_neurons[:10]

# %%
%%capture
endo_coverage_neurons = []

# First filter the CpG file once (no need to do it for each peak file)
temp_cpg = "temp_filtered_cpg.bed"
filter_standard_chromosomes(cpg_file, temp_cpg)

try:
    # Process each peak file
    for peak_file in endo_samples_nsc['Peaks']:
        # Create temporary filtered peak file
        temp_peak = "temp_filtered_peaks.bed"
        filter_standard_chromosomes(peak_file, temp_peak)
        
        # Calculate coverage using filtered peaks and filtered CpG file
        coverage = calculate_peak_cpg_coverage(temp_peak, temp_cpg, extend=100, genome_size_file=f"{data_dir}/genome.size")
        endo_coverage_neurons.extend(coverage)
        
        # Clean up peak file
        if os.path.exists(temp_peak):
            os.remove(temp_peak)
finally:
    # Clean up CpG file
    if os.path.exists(temp_cpg):
        os.remove(temp_cpg)

# %%
endo_coverage_neurons[:10]


# %%
# Create the visualization
plot_coverage_histograms(exo_coverage_neurons, endo_coverage_neurons, min_coverage=50, n_bins=30)

# %%
plot_coverage_histograms_overlayed(exo_coverage_neurons, endo_coverage_neurons, min_coverage=0, n_bins=30)

# %%
# Save the lists to files
with open(f'{coverage_dir}/exo_coverage_neurons.pkl', 'wb') as f:
    pickle.dump(exo_coverage_neurons, f)

with open(f'{coverage_dir}/endo_coverage_neurons.pkl', 'wb') as f:
    pickle.dump(endo_coverage_neurons, f)

# %%

# Load the lists back
with open(f'{coverage_dir}/exo_coverage_neurons.pkl', 'rb') as f:
    exo_coverage_neurons = pickle.load(f)

with open(f'{coverage_dir}/endo_coverage_neurons.pkl', 'rb') as f:
    endo_coverage_nsc = pickle.load(f)

# %% [markdown]
# # Only with the common peaks

# %% [markdown]
# ### NSC
# 

# %%
# Read the gene list
genes_df_nsc = pd.read_csv(f"{data_dir}/allgenes_NSC_total.csv")

# Filter for genes that have both types of promoters
common_genes_nsc = genes_df_nsc[
    (genes_df_nsc['Endogenous_Promoter'] == True) & 
    (genes_df_nsc['Exogenous_Promoter'] == True)
]
print(f"Number of genes with both promoter types: {len(common_genes_nsc)}")
print("Example genes:", common_genes_nsc['gene'].head().tolist())

# %%
common_genes_nsc.head()


# %%
%%capture
exo_coverage_nsc = []

# First filter the CpG file once (no need to do it for each peak file)
temp_cpg = "temp_filtered_cpg.bed"
filter_standard_chromosomes(cpg_file, temp_cpg)

try:
    # Process each peak file
    for peak_file in exo_samples_nsc['Peaks']:
        # Create temporary filtered peak file
        temp_peak = "temp_filtered_peaks.bed"
        filter_standard_chromosomes(peak_file, temp_peak)
        
        # Calculate coverage using filtered peaks and filtered CpG file
        coverage = calculate_peak_cpg_coverage(temp_peak, temp_cpg, extend=100, genome_size_file=f"{data_dir}/genome.size")
        exo_coverage_nsc.extend(coverage)
        
        # Clean up peak file
        if os.path.exists(temp_peak):
            os.remove(temp_peak)
finally:
    # Clean up CpG file
    if os.path.exists(temp_cpg):
        os.remove(temp_cpg)

# %%
exo_coverage_nsc[:10]


# %%
%%capture
endo_coverage_nsc = []

# First filter the CpG file once (no need to do it for each peak file)
temp_cpg = "temp_filtered_cpg.bed"
filter_standard_chromosomes(cpg_file, temp_cpg)

try:
    # Process each peak file
    for peak_file in endo_samples_nsc['Peaks']:
        # Create temporary filtered peak file
        temp_peak = "temp_filtered_peaks.bed"
        filter_standard_chromosomes(peak_file, temp_peak)
        
        # Calculate coverage using filtered peaks and filtered CpG file
        coverage = calculate_peak_cpg_coverage(temp_peak, temp_cpg, extend=100, genome_size_file=f"{data_dir}/genome.size")
        endo_coverage_nsc.extend(coverage)
        
        # Clean up peak file
        if os.path.exists(temp_peak):
            os.remove(temp_peak)
finally:
    # Clean up CpG file
    if os.path.exists(temp_cpg):
        os.remove(temp_cpg)

# %%
# Create the visualization
plot_coverage_histograms(exo_coverage_nsc, endo_coverage_nsc, min_coverage=50, n_bins=20)

# %%
plot_coverage_histograms_overlayed(exo_coverage_nsc, endo_coverage_nsc, min_coverage=0, n_bins=20)

# %%
# Save the lists to files
with open(f'{coverage_dir}/selected_exo_coverage_nsc.pkl', 'wb') as f:
    pickle.dump(exo_coverage_nsc, f)

with open(f'{coverage_dir}/selected_endo_coverage_nsc.pkl', 'wb') as f:
    pickle.dump(endo_coverage_nsc, f)

# %%
# Load the lists back
with open(f'{coverage_dir}/selected_exo_coverage_nsc.pkl', 'rb') as f:
    exo_coverage_nsc = pickle.load(f)

with open(f'{coverage_dir}/selected_endo_coverage_nsc.pkl', 'rb') as f:
    endo_coverage_nsc = pickle.load(f)

# %% [markdown]
# ### Neurons
# 

# %%
# Read the gene list
genes_df_neurons = pd.read_csv(f"{data_dir}/allgenes_NEU_total.csv")

# Filter for genes that have both types of promoters
common_genes_neurons = genes_df_neurons[
    (genes_df_neurons['Endogenous_Promoter'] == True) & 
    (genes_df_neurons['Exogenous_Promoter'] == True)
]
print(f"Number of genes with both promoter types: {len(common_genes_neurons)}")
print("Example genes:", common_genes_neurons['gene'].head().tolist())

# %%
data_dir

# %%
%%capture
exo_coverage_neurons = []
for peak_file in exo_samples_neurons['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, common_genes_neurons, f"{data_dir}/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file, genome_size_file=f"{data_dir}/genome.size") # Corrected path
    exo_coverage_neurons.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# %%
exo_coverage_neurons[:10]

# %%
%%capture
endo_coverage_neurons = []
for peak_file in endo_samples_neurons['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, common_genes_neurons, f"{data_dir}/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file, genome_size_file=f"{data_dir}/genome.size")
    endo_coverage_neurons.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# %%
endo_coverage_neurons[:10]

# %%
if len(exo_coverage_neurons) > 0 and len(endo_coverage_neurons) > 0:
    plot_coverage_histograms(exo_coverage_neurons, endo_coverage_neurons, min_coverage=50, n_bins=20)
else:
    print("Warning: Coverage lists are empty, cannot plot histograms")

# %%
if len(exo_coverage_neurons) > 0 and len(endo_coverage_neurons) > 0:
    plot_coverage_histograms_overlayed(exo_coverage_neurons, endo_coverage_neurons, min_coverage=0, n_bins=20)
else:
    print("Warning: Coverage lists are empty, cannot plot histograms")

# %%
# Save the lists to files
with open(f'{coverage_dir}/selected_exo_coverage_neurons.pkl', 'wb') as f:
    pickle.dump(exo_coverage_neurons, f)    

with open(f'{coverage_dir}/selected_endo_coverage_neurons.pkl', 'wb') as f:
    pickle.dump(endo_coverage_neurons, f)

# %%
# Load the lists back
with open(f'{coverage_dir}/selected_exo_coverage_neurons.pkl', 'rb') as f:
    exo_coverage_neurons = pickle.load(f)

with open(f'{coverage_dir}/selected_endo_coverage_neurons.pkl', 'rb') as f:
    endo_coverage_neurons = pickle.load(f)

# %% [markdown]
# # DE expression bins

# %%
DEA_NEU = pd.read_csv(f"{root_dir}/DATA/DEA_NEU.csv", header=0)
DEA_NSC = pd.read_csv(f"{root_dir}/DATA/DEA_NSC.csv", header=0)

# %%
DEA_NSC.head()

# %%
# Calculate quantiles and assign expression levels
q33_, q66_ = DEA_NSC['baseMean'].quantile([0.33, 0.66])
DEA_NSC['expression_level'] = DEA_NSC['baseMean'].apply(lambda x: get_expression_level(x, q33_, q66_))

q33, q66 = DEA_NEU['baseMean'].quantile([0.33, 0.66])
DEA_NEU['expression_level'] = DEA_NEU['baseMean'].apply(lambda x: get_expression_level(x, q33, q66))

# %% [markdown]
# ## UP

# %%
DEA_NEU_up = DEA_NEU[(DEA_NEU['padj'] < 0.05) & (DEA_NEU['log2FoldChange'] > 0)]
DEA_NSC_up = DEA_NSC[(DEA_NSC['padj'] < 0.05) & (DEA_NSC['log2FoldChange'] > 0)]

# %% [markdown]
# ### NSC

# %%
# Calculate quantiles for expression binning
# q33, q66 = DEA_NSC_up['baseMean'].quantile([0.33, 0.66])

# %%
# # Let's see the distribution
# print("Expression level boundaries:")
# print(f"Low: baseMean <= {q33:.2f}")
# print(f"Medium: {q33:.2f} < baseMean <= {q66:.2f}")
# print(f"High: baseMean > {q66:.2f}")
# print("\nNumber of genes in each category:")
# print(DEA_NSC_up['expression_level'].value_counts())

# %%
DEA_NSC_up_high = DEA_NSC_up[DEA_NSC_up.expression_level == 'High']
DEA_NSC_up_medium = DEA_NSC_up[DEA_NSC_up.expression_level == 'Medium']
DEA_NSC_up_low = DEA_NSC_up[DEA_NSC_up.expression_level == 'Low']

# %%
%%capture
# Calculate coverage for filtered peaks
exo_coverage_up_high = []
for peak_file in exo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_up_high, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_up_high.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
exo_coverage_up_medium = []
for peak_file in exo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_up_medium, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_up_medium.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)
    
# Calculate coverage for filtered peaks
exo_coverage_up_low = []
for peak_file in exo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_up_low, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_up_low.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# %%
%%capture
# Calculate coverage for filtered peaks
endo_coverage_up_high = []
for peak_file in endo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_up_high, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_up_high.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
endo_coverage_up_medium = []
for peak_file in endo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_up_medium, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_up_medium.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
endo_coverage_up_low = []
for peak_file in endo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_up_low, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_up_low.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)


# %%
# Exogenous
if len(exo_coverage_up_high) > 0 and len(exo_coverage_up_medium) > 0 and len(exo_coverage_up_low) > 0:
    fig = plot_coverage_histograms_expression(exo_coverage_up_high, exo_coverage_up_medium, exo_coverage_up_low, min_coverage=0, p_t = 80, n_bins=15)
    plt.show()


# %%
# Endogenous
if len(endo_coverage_up_high) > 0 and len(endo_coverage_up_medium) > 0 and len(endo_coverage_up_low) > 0:
    fig = plot_coverage_histograms_expression(endo_coverage_up_high, endo_coverage_up_medium, endo_coverage_up_low, min_coverage=0, p_t = 80, n_bins=25)
    plt.show()


# %%
# Save the lists to files
with open(f'{coverage_dir}/expression_coverage_up_high_nsc.pkl', 'wb') as f:
    pickle.dump(exo_coverage_up_high, f)

with open(f'{coverage_dir}/expression_coverage_up_medium_nsc.pkl', 'wb') as f:
    pickle.dump(exo_coverage_up_medium, f)

with open(f'{coverage_dir}/expression_coverage_up_low_nsc.pkl', 'wb') as f:
    pickle.dump(exo_coverage_up_low, f)

# %% [markdown]
# ### Neurons
# 

# %%
# Calculate quantiles for expression binning
# q33, q66 = DEA_NEU_up['baseMean'].quantile([0.33, 0.66])


# %%
# Let's see the distribution
# print("Expression level boundaries:")
# print(f"Low: baseMean <= {q33:.2f}")
# print(f"Medium: {q33:.2f} < baseMean <= {q66:.2f}")
# print(f"High: baseMean > {q66:.2f}")
# print("\nNumber of genes in each category:")
# print(DEA_NSC_up['expression_level'].value_counts())

# %%
DEA_NEU_up_high = DEA_NEU_up[DEA_NEU_up.expression_level == 'High']
DEA_NEU_up_medium = DEA_NEU_up[DEA_NEU_up.expression_level == 'Medium']
DEA_NEU_up_low = DEA_NEU_up[DEA_NEU_up.expression_level == 'Low']

# %%
%%capture
# Calculate coverage for filtered peaks
exo_coverage_up_high = []
for peak_file in exo_samples_neurons['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NEU_up_high, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_up_high.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
exo_coverage_up_medium = []
for peak_file in exo_samples_neurons['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NEU_up_medium, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_up_medium.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)
    
# Calculate coverage for filtered peaks
exo_coverage_up_low = []
for peak_file in exo_samples_neurons['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NEU_up_low, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_up_low.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# %%
%%capture
# Calculate coverage for filtered peaks
endo_coverage_up_high = []
for peak_file in endo_samples_neurons['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NEU_up_high, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_up_high.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
endo_coverage_up_medium = []
for peak_file in endo_samples_neurons['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NEU_up_medium, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_up_medium.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
endo_coverage_up_low = []
for peak_file in endo_samples_neurons['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NEU_up_low, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_up_low.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)


# %%
# Exogenous
if len(exo_coverage_up_high) > 0 and len(exo_coverage_up_medium) > 0 and len(exo_coverage_up_low) > 0:
    fig = plot_coverage_histograms_expression(exo_coverage_up_high, exo_coverage_up_medium, exo_coverage_up_low, min_coverage=0, p_t = 80, n_bins=15)
    plt.show()


# %%
# Endogenous
if len(endo_coverage_up_high) > 0 and len(endo_coverage_up_medium) > 0 and len(endo_coverage_up_low) > 0:
    fig = plot_coverage_histograms_expression(endo_coverage_up_high, endo_coverage_up_medium, endo_coverage_up_low, min_coverage=0, p_t = 80, n_bins=25)
    plt.show()


# %%
# Save the lists to files
with open(f'{coverage_dir}/expression_coverage_up_high_neurons.pkl', 'wb') as f:
    pickle.dump(exo_coverage_up_high, f)

with open(f'{coverage_dir}/expression_coverage_up_medium_neurons.pkl', 'wb') as f:
    pickle.dump(exo_coverage_up_medium, f)

with open(f'{coverage_dir}/expression_coverage_up_low_neurons.pkl', 'wb') as f:
    pickle.dump(exo_coverage_up_low, f)

# %% [markdown]
# ## Down

# %%
DEA_NEU_down = DEA_NEU[(DEA_NEU['padj'] < 0.05) & (DEA_NEU['log2FoldChange'] < 0)]
DEA_NSC_down = DEA_NSC[(DEA_NSC['padj'] < 0.05) & (DEA_NSC['log2FoldChange'] < 0)]

# %% [markdown]
# ### NSC
# 

# %%
# # Calculate quantiles for expression binning
# q33, q66 = DEA_NSC_down['baseMean'].quantile([0.33, 0.66])

# %%
# # Let's see the distribution
# print("Expression level boundaries:")
# print(f"Low: baseMean <= {q33:.2f}")
# print(f"Medium: {q33:.2f} < baseMean <= {q66:.2f}")
# print(f"High: baseMean > {q66:.2f}")
# print("\nNumber of genes in each category:")
# print(DEA_NSC_down['expression_level'].value_counts())

# %%
DEA_NSC_down_high = DEA_NSC_down[DEA_NSC_down.expression_level == 'High']
DEA_NSC_down_medium = DEA_NSC_down[DEA_NSC_down.expression_level == 'Medium']
DEA_NSC_down_low = DEA_NSC_down[DEA_NSC_down.expression_level == 'Low']

# %%
%%capture
exo_coverage_down_high = []
for peak_file in exo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_down_high, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_down_high.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
exo_coverage_down_medium = []
for peak_file in exo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_down_medium, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_down_medium.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
exo_coverage_down_low = []
for peak_file in exo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_down_low, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_down_low.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# %%
%%capture
# Calculate coverage for filtered peaks
endo_coverage_down_high = []
for peak_file in endo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_down_high, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_down_high.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
endo_coverage_down_medium = []
for peak_file in endo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_down_medium, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_down_medium.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
endo_coverage_down_low = []
for peak_file in endo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_down_low, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_down_low.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)


# %%
# Exogenous
if len(exo_coverage_down_high) > 0 and len(exo_coverage_down_medium) > 0 and len(exo_coverage_down_low) > 0:
    fig = plot_coverage_histograms_expression(exo_coverage_down_high, exo_coverage_down_medium, exo_coverage_down_low, min_coverage=0, p_t = 80, n_bins=15)
    plt.show()


# %%
# Endogenous
if len(endo_coverage_down_high) > 0 and len(endo_coverage_down_medium) > 0 and len(endo_coverage_down_low) > 0:
    fig = plot_coverage_histograms_expression(endo_coverage_down_high, endo_coverage_down_medium, endo_coverage_down_low, min_coverage=0, p_t = 80, n_bins=15)
    plt.show()


# %% [markdown]
# ### Neurons

# %%
# # Calculate quantiles for expression binning
# q33, q66 = DEA_NEU_down['baseMean'].quantile([0.33, 0.66])


# %%
# Let's see the distribution
# print("Expression level boundaries:")
# print(f"Low: baseMean <= {q33:.2f}")
# print(f"Medium: {q33:.2f} < baseMean <= {q66:.2f}")
# print(f"High: baseMean > {q66:.2f}")
# print("\nNumber of genes in each category:")
# print(DEA_NEU_down['expression_level'].value_counts())

DEA_NEU_down_high = DEA_NEU_down[DEA_NEU_down.expression_level == 'High']
DEA_NEU_down_medium = DEA_NEU_down[DEA_NEU_down.expression_level == 'Medium']
DEA_NEU_down_low = DEA_NEU_down[DEA_NEU_down.expression_level == 'Low']

# %%
%%capture
exo_coverage_down_high = []
for peak_file in exo_samples_neurons['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NEU_down_high, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_down_high.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
exo_coverage_down_medium = []
for peak_file in exo_samples_neurons['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NEU_down_medium, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_down_medium.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
exo_coverage_down_low = []
for peak_file in exo_samples_neurons['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NEU_down_low, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_down_low.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# %%
%%capture
# Calculate coverage for filtered peaks
endo_coverage_down_high = []
for peak_file in endo_samples_neurons['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NEU_down_high, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_down_high.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
endo_coverage_down_medium = []
for peak_file in endo_samples_neurons['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NEU_down_medium, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_down_medium.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
endo_coverage_down_low = []
for peak_file in endo_samples_neurons['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NEU_down_low, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_down_low.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)


# %%
# Exogenous
if len(exo_coverage_down_high) > 0 and len(exo_coverage_down_medium) > 0 and len(exo_coverage_down_low) > 0:
    fig = plot_coverage_histograms_expression(exo_coverage_down_high, exo_coverage_down_medium, exo_coverage_down_low, min_coverage=0, p_t = 80, n_bins=15)
    plt.show()


# %%
# Endogenous
if len(endo_coverage_down_high) > 0 and len(endo_coverage_down_medium) > 0 and len(endo_coverage_down_low) > 0:
    fig = plot_coverage_histograms_expression(endo_coverage_down_high, endo_coverage_down_medium, endo_coverage_down_low, min_coverage=0, p_t = 80, n_bins=15)
    plt.show()


# %% [markdown]
# ## Not DE

# %% [markdown]
# ### NSC

# %%
DEA_NEU_not_de = DEA_NEU[DEA_NEU['padj'] > 0.05]
DEA_NSC_not_de = DEA_NSC[DEA_NSC['padj'] > 0.05]

# %%
# Calculate quantiles for expression binning
# q33, q66 = DEA_NSC_not_de['baseMean'].quantile([0.33, 0.66])

# %%
# # Let's see the distribution
# print("Expression level boundaries:")
# print(f"Low: baseMean <= {q33:.2f}")
# print(f"Medium: {q33:.2f} < baseMean <= {q66:.2f}")
# print(f"High: baseMean > {q66:.2f}")
# print("\nNumber of genes in each category:")
# print(DEA_NSC_not_de['expression_level'].value_counts())

# %%
DEA_NSC_not_de_high = DEA_NSC_not_de[DEA_NSC_not_de.expression_level == 'High']
DEA_NSC_not_de_medium = DEA_NSC_not_de[DEA_NSC_not_de.expression_level == 'Medium']
DEA_NSC_not_de_low = DEA_NSC_not_de[DEA_NSC_not_de.expression_level == 'Low']

# %%
%%capture
exo_coverage_not_de_high = []
for peak_file in exo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_not_de_high, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_not_de_high.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
exo_coverage_not_de_medium = []
for peak_file in exo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_not_de_medium, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_not_de_medium.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
exo_coverage_not_de_low = []
for peak_file in exo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_not_de_low, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    exo_coverage_not_de_low.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# %%
%%capture
# Calculate coverage for filtered peaks
endo_coverage_not_de_high = []
for peak_file in endo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_not_de_high, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_not_de_high.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
endo_coverage_not_de_medium = []
for peak_file in endo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_not_de_medium, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_not_de_medium.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)

# Calculate coverage for filtered peaks
endo_coverage_not_de_low = []
for peak_file in endo_samples_nsc['Peaks']:
    filtered_peaks = get_common_peaks(peak_file, DEA_NSC_not_de_low, f"{root_dir}/DATA/genes.bed")
    coverage = calculate_peak_cpg_coverage(filtered_peaks, cpg_file,)
    endo_coverage_not_de_low.extend(coverage)
    subprocess.run("rm temp_filtered_peaks.bed", shell=True)


# %%
# Exogenous
if len(exo_coverage_not_de_high) > 0 and len(exo_coverage_not_de_medium) > 0 and len(exo_coverage_not_de_low) > 0:
    fig = plot_coverage_histograms_expression(exo_coverage_not_de_high, exo_coverage_not_de_medium, exo_coverage_not_de_low, min_coverage=0, p_t = 80, n_bins=15)
    plt.show()


# %%
# Endogenous
if len(endo_coverage_not_de_high) > 0 and len(endo_coverage_not_de_medium) > 0 and len(endo_coverage_not_de_low) > 0:
    fig = plot_coverage_histograms_expression(endo_coverage_not_de_high, endo_coverage_not_de_medium, endo_coverage_not_de_low, min_coverage=0, p_t = 80, n_bins=15)
    plt.show()


# %% [markdown]
# # Plot peaks signals

# %%
import gzip
import numpy as np

def plot_peak_signals(peak_file, bigwig_files, window_size=1000, num_peaks=5):
    """
    Plot signal profiles for selected peaks across multiple bigWig files
    """
    import gzip  # Add this import at the top
    
    # Select random peaks
    cmd = f"shuf -n {num_peaks} {peak_file} > temp_random_peaks.bed"
    subprocess.run(cmd, shell=True)
    
    # Center peaks and extend by window_size
    half_window = window_size // 2
    signals = {}
    
    for bw_file in bigwig_files:
        # Use deepTools computeMatrix to get signal values
        matrix_file = "temp_matrix.gz"
        cmd = (f"computeMatrix reference-point"
               f" --referencePoint center"
               f" -b {half_window} -a {half_window}"
               f" -R temp_random_peaks.bed"
               f" -S {bw_file}"
               f" --skipZeros"
               f" -o {matrix_file}")
        
        # Run command and check for errors
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error running computeMatrix for {bw_file}:")
            print(result.stderr)
            continue
            
        try:
            # Read the gzipped matrix
            data = []
            with gzip.open(matrix_file, 'rt') as f:
                for line in f:
                    if not line.startswith(('@', '#')):
                        values = line.strip().split('\t')[6:]  # Skip region info columns
                        data.append([float(x) for x in values])
            
            if data:  # Only store if we got some data
                signals[os.path.basename(bw_file)] = np.array(data)
            
        except Exception as e:
            print(f"Error processing matrix file for {bw_file}:")
            print(str(e))
            continue
        finally:
            # Clean up matrix file
            if os.path.exists(matrix_file):
                os.remove(matrix_file)
    
    if not signals:
        print("No valid signal data was obtained.")
        return
        
    # Create figure
    fig, axes = plt.subplots(num_peaks, 1, figsize=(10, 3*num_peaks))
    if num_peaks == 1:
        axes = [axes]
    
    # Create x-axis values (distance from center)
    bins = signals[list(signals.keys())[0]].shape[1]
    x = np.linspace(-half_window, half_window, bins)
    
    # Plot each peak
    for i in range(num_peaks):
        for bw_file, signal in signals.items():
            if i < len(signal):  # Make sure we have enough peaks
                axes[i].plot(x, signal[i], label=os.path.basename(bw_file))
                
        axes[i].axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        axes[i].set_xlabel('Distance from peak center (bp)')
        axes[i].set_ylabel('Signal')
        axes[i].legend()
    
    plt.tight_layout()
    plt.show()
    
    # # Final cleanup
    # if os.path.exists("temp_random_peaks.bed"):
    #     os.remove("temp_random_peaks.bed")



# %%
# Get peaks for genes with both promoter types
filtered_peaks = get_common_peaks(exo_samples['Peaks'].iloc[0], common_genes_nsc, f"{root_dir}/DATA/genes.bed")

# Save peaks to file - handle string output from get_common_peaks
with open(f"{root_dir}/results/filtered_peaks.bed", 'w') as f:
    f.write(filtered_peaks)

# Load peaks back
filtered_peaks = pd.read_csv(f"{root_dir}/results/filtered_peaks.bed", sep='\t', names=['chrom', 'start', 'end', 'name', 'score', 'strand'])

# %%
# Use the function
bigwig_dir = f"{root_dir}/results/bigwig"
bigwig_files = [
    f"{bigwig_dir}/Exogenous_Neuron.bw",
    f"{bigwig_dir}/Endogenous_Neuron.bw"
]

# %%
# Plot the signals
plot_peak_signals(filtered_peaks, bigwig_files, window_size=2000, num_peaks=20)

# %%



