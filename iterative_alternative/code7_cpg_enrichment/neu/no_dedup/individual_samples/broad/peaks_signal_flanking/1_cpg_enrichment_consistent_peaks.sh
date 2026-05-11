#!/bin/bash
#SBATCH --job-name=cpg_enrichment_nsc_no_dedup
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment/nsc/no_dedup/individual_samples/broad/cpg_enrichment_consistent_peaks_with_flanking.err"
#SBATCH --output="logs/cpg_enrichment/nsc/no_dedup/individual_samples/broad/cpg_enrichment_consistent_peaks_with_flanking.out"
#SBATCH --array=0-9  # Process in 10 chunks

ALIGNMENT="results_1b"
PEAKS="broad"
CELL_LINE="NSC"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
DATA_DIR="${BASE_DIR}/DATA"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
SCRIPT_DIR="${BASE_DIR}/scripts/cpg_enrichment"

cd $WORKING_DIR || exit 1

PEAKS_DIR="results/no_dedup/peaks/${PEAKS}"
RESULTS_DIR="results/no_dedup/cpg_enrichment/${CELL_LINE}/${PEAKS}"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create directories
mkdir -p "${RESULTS_DIR}/exo"
mkdir -p "${RESULTS_DIR}/endo"
mkdir -p "${RESULTS_DIR}/mecp2_cpg_enrichment_parallel_with_flanking"

# Clean up any existing peak files to avoid duplicates
rm -f "${RESULTS_DIR}/exo"/*.${PEAKS}Peak
rm -f "${RESULTS_DIR}/endo"/*.${PEAKS}Peak

echo "Copying and organizing peaks from ${PEAKS_DIR}..."

# Copy peaks
for sample in NSCv1 NSCv2 NSCv3; do
    cp "${PEAKS_DIR}/${sample}_${PEAKS}_peaks.filtered.${PEAKS}Peak" "${RESULTS_DIR}/exo/${sample}_peaks.${PEAKS}Peak"
done

for sample in NSCM1 NSCM2 NSCM3; do
    cp "${PEAKS_DIR}/${sample}_${PEAKS}_peaks.filtered.${PEAKS}Peak" "${RESULTS_DIR}/endo/${sample}_peaks.${PEAKS}Peak"
done

# Process chunk
python ${SCRIPT_DIR}/cpg_enrichment_consistent_peaks_with_flanking.py \
    --exo-peaks-dir ${RESULTS_DIR}/exo \
    --endo-peaks-dir ${RESULTS_DIR}/endo \
    --cpg-islands ${DATA_DIR}/cpg_islands.bed \
    --output-dir ${RESULTS_DIR}/cpg_enrichment_consistent_peaks_with_flanking \
    --bam-dir ${WORKING_DIR}/${ALIGNMENT}/aligned \
    --chunk-id $SLURM_ARRAY_TASK_ID \
    --peaks-type ${PEAKS} \
    --total-chunks 10