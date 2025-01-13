#!/bin/bash
#SBATCH --job-name=peaks_enrichment_and_genes_regulation_NSC_5.
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_enrichment_and_genes_regulation_NSC_5.err"
#SBATCH --output="logs/peaks_enrichment_and_genes_regulation_NSC_5.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

mkdir -p logs

# Define base directory
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"

# Define directory structure
WORKING_DIR="${BASE_DIR}/iterative_alternative"
DATA_DIR="${BASE_DIR}/DATA"

# Define experiment-specific directories
ALIGN_DIR="${WORKING_DIR}/results_1b"
PEAKS_DIR="${WORKING_DIR}/results_2_align2_005/peaks/narrow"
RESULTS_DIR="${WORKING_DIR}/results_5_align2_005/peaks_enrichment_and_genes_regulation_NSC_5"

# Create results directory
mkdir -p "${RESULTS_DIR}/peaks"

echo "Copying and renaming peaks from ${PEAKS_DIR} to ${RESULTS_DIR}..."

# Copy and rename peak files with correct naming schema
for sample in NSCv1 NSCv2 NSCv3 NSCM1 NSCM2 NSCM3; do
    peak_file="${PEAKS_DIR}/${sample}_narrow_peaks.filtered.narrowPeak"
    if [ -f "${peak_file}" ]; then
        cp "${peak_file}" "${RESULTS_DIR}/peaks/${sample}_peaks.narrowPeak"
        echo "Copied ${peak_file} to ${RESULTS_DIR}/peaks/${sample}_peaks.narrowPeak"
    else
        echo "Warning: Source file not found for sample ${sample}: ${peak_file}"
    fi
done

# Run analysis script
python -u ../scripts/peaks_enrichment_and_genes_regulation/peaks_enrichment_and_genes_regulation_NSC_5.py \
    --working-dir "${WORKING_DIR}" \
    --data-dir "${DATA_DIR}" \
    --results-dir "${RESULTS_DIR}" \
    2>&1 | tee "logs/peaks_enrichment_and_genes_regulation_NSC_5.out"


#####################################################################
# This script analyzes enrichment of peaks in NSC samples
#
# Input files:
# - Peak files in ${RESULTS_DIR}/peaks/
#   - NSC samples (NSCv1-3, NSCM1-3) narrow peak files
# - Original aligned BAM files in ${DATA_DIR}/aligned/
#
# The script:
# 1. Copies and renames peak files from results_2_align2_new_005 to results_5_align2_new_005
# 2. Calls Python script analyze_enrichment_NSC_5.py to:
#    - Analyze peak overlaps between NSC samples
#    - Calculate enrichment statistics
#    - Generate enrichment plots and reports
#
# Output files in ${RESULTS_DIR}:
# - Renamed peak files
# - Enrichment analysis results
# - Overlap statistics
# - Visualization plots
#####################################################################
