#!/bin/bash
#SBATCH --job-name=pa_no_dedup_NSC_combined
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_annotation/no_dedup/peaks_annotation_NSC_combined.err"
#SBATCH --output="logs/peaks_annotation/no_dedup/peaks_annotation_NSC_combined.out"

set -e  # Exit on error
# set -x  # Print commands as they're executed

# Define paths
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
SCRIPT_DIR="${BASE_DIR}/scripts/peaks_annotation"

INPUT_DIR="results/no_dedup/peaks"
RESULTS_DIR="results/no_dedup/peaks_annotation"

cd $WORKING_DIR || exit 1

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Process both narrow and broad peaks
for PEAK_TYPE in narrow broad; do
    echo "Processing ${PEAK_TYPE} peaks..."
    
    # Input paths
    PEAKS_DIR="${INPUT_DIR}/${PEAK_TYPE}"
    
    # Output paths
    OUTPUT_DIR="${RESULTS_DIR}/NSC_combined_${PEAK_TYPE}"

    # Remove OUTPUT_DIR if it exists and recreate it
    mkdir -p "${OUTPUT_DIR}"

    # Run the R script with quoted arguments
    Rscript "${SCRIPT_DIR}/peaks_annotation_NSC_combined.R" \
        --peaks-dir "${PEAKS_DIR}" \
        --output-dir "${OUTPUT_DIR}" \
        --peak-type "${PEAK_TYPE}" \
        2>&1 | tee "logs/peaks_annotation/no_dedup/peaks_annotation_NSC_combined_${PEAK_TYPE}.out"
done
