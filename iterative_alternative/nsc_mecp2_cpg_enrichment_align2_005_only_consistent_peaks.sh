#!/bin/bash
#SBATCH --job-name=nsc_cpg_enrichment_align2_005_consistent_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/nsc_cpg_enrichment_align2_005_consistent_peaks.err"
#SBATCH --output="logs/nsc_cpg_enrichment_align2_005_consistent_peaks.out"
#SBATCH --array=0-9  # Process in 10 chunks

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

EXPERIMENT_IN="align2_005"
EXPERIMENT_OUT="align2_005_consistent_peaks"

PEAKS_EXPERIMENT="results_2_${EXPERIMENT_IN}"
ALIGNMENT_EXPERIMENT="results_1b"
CELL_LINE="NSC"

WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative"
DATA_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA"
RESULTS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment_${EXPERIMENT_OUT}/${CELL_LINE}"

# Create directories
mkdir -p "${RESULTS_DIR}/exo"
mkdir -p "${RESULTS_DIR}/endo"
mkdir -p logs
mkdir -p ${RESULTS_DIR}/mecp2_cpg_enrichment_parallel

# Clean up any existing peak files to avoid duplicates
rm -f "${RESULTS_DIR}/exo"/*.narrowPeak
rm -f "${RESULTS_DIR}/endo"/*.narrowPeak

echo "Copying and organizing peaks from ${PEAKS_EXPERIMENT}..."

# Copy exogenous peaks (virus samples) - removed symlinks
for sample in NSCv1 NSCv2 NSCv3; do
    cp "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/${PEAKS_EXPERIMENT}/peaks/narrow/${sample}_narrow_peaks.narrowPeak" "${RESULTS_DIR}/exo/${sample}_peaks.narrowPeak"
done

# Copy endogenous peaks (M samples) - removed symlinks
for sample in NSCM1 NSCM2 NSCM3; do
    cp "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/${PEAKS_EXPERIMENT}/peaks/narrow/${sample}_narrow_peaks.narrowPeak" "${RESULTS_DIR}/endo/${sample}_peaks.narrowPeak"
done

# Process chunk
python ../scripts/analyze_mecp2_cpg_enrichment_only_consistent_peaks.py \
    --exo-peaks-dir ${RESULTS_DIR}/exo \
    --endo-peaks-dir ${RESULTS_DIR}/endo \
    --cpg-islands ${DATA_DIR}/cpg_islands.bed \
    --output-dir ${RESULTS_DIR}/mecp2_cpg_enrichment_parallel \
    --bam-dir ${WORKING_DIR}/${ALIGNMENT_EXPERIMENT}/aligned \
    --chunk-id $SLURM_ARRAY_TASK_ID \
    --total-chunks 10

# If this is the last chunk, wait for all other chunks and combine results
if [ $SLURM_ARRAY_TASK_ID -eq 9 ]; then
    echo "Waiting for all chunks to complete..."
    
    # Wait for all chunk files
    expected_files=10
    while true; do
        actual_files=$(ls ${RESULTS_DIR}/mecp2_cpg_enrichment_parallel/chunk_*.csv 2>/dev/null | wc -l)
        if [ "$actual_files" -eq "$expected_files" ]; then
            break
        fi
        echo "Found $actual_files/$expected_files chunk files. Waiting..."
        sleep 30
    done
    
    echo "All chunks found. Combining results..."
    python - <<EOF
import pandas as pd
import glob
import os

results_dir = "${RESULTS_DIR}/mecp2_cpg_enrichment_parallel"
chunks = sorted(glob.glob(os.path.join(results_dir, "chunk_*.csv")))
if len(chunks) != 10:
    raise ValueError(f"Expected 10 chunk files, found {len(chunks)}")

combined = pd.concat([pd.read_csv(f) for f in chunks])
output_file = os.path.join(results_dir, "mecp2_cpg_enrichment_parallel.csv")
combined.to_csv(output_file, index=False)
print(f"Combined {len(chunks)} chunks into {output_file}")
print(f"Total rows in combined file: {len(combined)}")
EOF
fi
