#!/bin/bash
#SBATCH --job-name=peaks_annotation_NSC_single_replicates
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_annotation_NSC_single_replicates.err"
#SBATCH --output="logs/peaks_annotation_NSC_single_replicates.out"
#SBATCH --array=0-5  # One job per input file

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

mkdir -p logs

GTF_PATH="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA/gencode.vM10.annotation.gtf"
PEAKS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results_2_align2_005/peaks/narrow"
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results_5_align2_005/peaks_annotation_NSC_single_replicates"

# Run the parallel analysis
python -u ../scripts/peaks_annotation/peaks_annotation_NSC_single_replicates.py \
    --gtf-path "$GTF_PATH" \
    --peaks-dir "$PEAKS_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --task-id "$SLURM_ARRAY_TASK_ID" \
    --total-tasks "$SLURM_ARRAY_TASK_COUNT" \
    2>&1 | tee "logs/peaks_annotation_NSC_single_replicates${SLURM_ARRAY_TASK_ID}.out"

# If this is the last task, combine results and create the final plot
if [ "$SLURM_ARRAY_TASK_ID" -eq "$((SLURM_ARRAY_TASK_COUNT-1))" ]; then
    # Wait for all other tasks to finish
    sleep 30  # Give time for other tasks to finish writing their files
    
    # Combine results and create final plot
    python - <<EOF
import json
import os
from scripts.gene_annotation_analysis_NSC import plot_distributions

# Load all partial results
results = {}
temp_dir = os.path.join("$OUTPUT_DIR", "temp")
for i in range($SLURM_ARRAY_TASK_COUNT):
    with open(os.path.join(temp_dir, f"results_{i}.json")) as f:
        results.update(json.load(f))

# Separate results into endogenous and exogenous replicates
endo_reps = [results[k] for k in ['npc_endo_1', 'npc_endo_2', 'npc_endo_3'] 
             if k in results]
exo_reps = [results[k] for k in ['npc_exo_1', 'npc_exo_2', 'npc_exo_3'] 
            if k in results]

# Create final plot
if endo_reps and exo_reps:
    plot_distributions(endo_reps, exo_reps, "$OUTPUT_DIR")

# Print final statistics
for condition, categories in results.items():
    print(f"\nDistribution for {condition}:")
    total = sum(categories.values())
    for category, count in categories.items():
        percentage = (count / total) * 100
        print(f"{category}: {count} peaks ({percentage:.1f}%)")

# Clean up temporary files
import shutil
shutil.rmtree(temp_dir)
EOF
fi
