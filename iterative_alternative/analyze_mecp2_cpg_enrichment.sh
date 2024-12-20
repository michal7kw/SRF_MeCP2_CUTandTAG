#!/bin/bash
#SBATCH --job-name=mecp2_cpg_parallel
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/mecp2_cpg_parallel.err"
#SBATCH --output="logs/mecp2_cpg_parallel.out"
#SBATCH --array=0-9  # Process in 10 chunks

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create logs directory if it doesn't exist
mkdir -p logs
mkdir -p ${RESULTS_DIR}/mecp2_cpg_enrichment_parallel

WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative"
DATA_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/DATA"
RESULTS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment"

# Process chunk
python ../scripts/analyze_mecp2_cpg_enrichment.py \
    --exo-peaks-dir ${RESULTS_DIR}/exo \
    --endo-peaks-dir ${RESULTS_DIR}/endo \
    --cpg-islands ${DATA_DIR}/cpg_islands.bed \
    --output-dir ${RESULTS_DIR}/mecp2_cpg_enrichment_parallel \
    --bam-dir ${WORKING_DIR}/results_1/aligned \
    --chunk-id $SLURM_ARRAY_TASK_ID \
    --total-chunks 10

# If this is the last chunk, wait for all other chunks and combine results
if [ $SLURM_ARRAY_TASK_ID -eq 9 ]; then
    echo "Waiting for all chunks to complete..."
    # Wait for all chunk files to be created
    while [ $(ls ${RESULTS_DIR}/mecp2_cpg_enrichment_parallel/chunk_*.csv 2>/dev/null | wc -l) -lt 12 ]; do
        sleep 30
    done
    
    echo "Combining results..."
    python - <<EOF
import pandas as pd
import glob
import os

results_dir = "${RESULTS_DIR}/mecp2_cpg_enrichment_parallel"
chunks = glob.glob(os.path.join(results_dir, "chunk_*.csv"))
combined = pd.concat([pd.read_csv(f) for f in chunks])
output_file = os.path.join(results_dir, "mecp2_cpg_enrichment_parallel.csv")
combined.to_csv(output_file, index=False)
print(f"Combined {len(chunks)} chunks into {output_file}")
print(f"Total rows in combined file: {len(combined)}")
EOF
fi
