#!/bin/bash
#SBATCH --job-name=pa_dedup_NSC_combined_cpg_specific
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_annotation/dedup/peaks_annotation_NSC_combined_cpg_specific.err"
#SBATCH --output="logs/peaks_annotation/dedup/peaks_annotation_NSC_combined_cpg_specific.out"

set -e  # Exit on error
# set -x  # Print commands as they're executed

# Define paths
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
SCRIPT_DIR="${BASE_DIR}/scripts/peaks_annotation"

INPUT_DIR="results/dedup/peaks"
RESULTS_DIR="results/dedup/peaks_annotation"

cd $WORKING_DIR || exit 1

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

mkdir -p logs

# Process both narrow and broad peaks
for PEAK_TYPE in narrow broad; do
    echo "Processing ${PEAK_TYPE} peaks..."
    
    # Input paths
    PEAKS_DIR="${INPUT_DIR}/${PEAK_TYPE}"
    
    # Output paths
    OUTPUT_DIR="${RESULTS_DIR}/NSC_combined_cpg_specific_${PEAK_TYPE}"

    # Remove OUTPUT_DIR if it exists and recreate it
    mkdir -p "${OUTPUT_DIR}"

    # Run R script with explicit argument names
    Rscript "${SCRIPT_DIR}/peaks_annotation_NSC_combined_cpg_specific.R" \
        --peaks-dir "${PEAKS_DIR}" \
        --output-dir "${OUTPUT_DIR}" \
        --peak-type "${PEAK_TYPE}" \
        2>&1 | tee "logs/peaks_annotation/dedup/peaks_annotation_NSC_combined_cpg_specific_${PEAK_TYPE}.out"
done
