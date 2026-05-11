#!/bin/bash
#SBATCH --job-name=copy_peaks_NEU
#SBATCH --account=kubacki.michal
#SBATCH --mem=4GB
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment/neu/no_dedup/individual_samples/broad/copy_peaks.err"
#SBATCH --output="logs/cpg_enrichment/neu/no_dedup/individual_samples/broad/copy_peaks.out"

PEAKS="broad"
CELL_LINE="Neu"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"

cd $WORKING_DIR || exit 1

PEAKS_DIR="results/no_dedup/peaks/${PEAKS}"
RESULTS_DIR="results/no_dedup/cpg_enrichment/${CELL_LINE}/${PEAKS}"

# Create directories
mkdir -p "${RESULTS_DIR}/exo"
mkdir -p "${RESULTS_DIR}/endo"

# Clean up any existing peak files to avoid duplicates
rm -f "${RESULTS_DIR}/exo"/*.${PEAKS}Peak
rm -f "${RESULTS_DIR}/endo"/*.${PEAKS}Peak

echo "Copying and organizing peaks from ${PEAKS_DIR}..."

# Copy peaks
for sample in NeuV1 NeuV2 NeuV3; do
    cp "${PEAKS_DIR}/${sample}_${PEAKS}_peaks.filtered.${PEAKS}Peak" "${RESULTS_DIR}/exo/${sample}_peaks.${PEAKS}Peak"
done

for sample in NeuM2 NeuM3; do
    cp "${PEAKS_DIR}/${sample}_${PEAKS}_peaks.filtered.${PEAKS}Peak" "${RESULTS_DIR}/endo/${sample}_peaks.${PEAKS}Peak"
done

echo "Peak files copying completed." 