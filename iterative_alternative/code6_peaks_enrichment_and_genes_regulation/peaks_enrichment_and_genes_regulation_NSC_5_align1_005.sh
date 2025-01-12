#!/bin/bash
#SBATCH --job-name=peaks_enrichment_and_genes_regulation_NSC_5_align1_005
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_enrichment_and_genes_regulation_NSC_5_align1_005.err"
#SBATCH --output="logs/peaks_enrichment_and_genes_regulation_NSC_5_align1_005.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

mkdir -p logs

WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative"
DATA_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/results_1"
INPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/results_2_align1_005"
RESULTS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/results_5_align1_005/peaks_enrichment_and_genes_regulation_NSC_5_align1_005"

echo "Copying and renaming peaks from results_2_align1_005 to ${RESULTS_DIR}..."
rm -rf "${RESULTS_DIR}/peaks"
mkdir -p "${RESULTS_DIR}/peaks"

# Copy and rename files
for sample in NSCv1 NSCv2 NSCv3 NSCM1 NSCM2 NSCM3; do
    if [ -f "${INPUT_DIR}/peaks/narrow/${sample}_narrow_peaks.narrowPeak" ]; then
        cp "${INPUT_DIR}/peaks/narrow/${sample}_narrow_peaks.narrowPeak" \
           "${RESULTS_DIR}/peaks/${sample}_peaks.narrowPeak"
    else
        echo "Warning: Source file not found for sample ${sample}"
    fi
done

python -u ../scripts/analyze_enrichment_NSC_5.py \
    --working-dir $WORKING_DIR \
    --data-dir $DATA_DIR \
    --results-dir $RESULTS_DIR \
    2>&1 | tee "logs/enrichment_NSC_5_align1_005.out"


#####################################################################
# This script analyzes enrichment of peaks in NSC samples
#
# Input files:
# - Peak files in ${RESULTS_DIR}/peaks/
#   - NSC samples (NSCv1-3, NSCM1-3) narrow peak files
# - Original aligned BAM files in ${DATA_DIR}/aligned/
#
# The script:
# 1. Copies and renames peak files from results_2_align1_005 to results_5_align1_005
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
