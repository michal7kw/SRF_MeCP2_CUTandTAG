#!/bin/bash
#SBATCH --job-name=peaks_annotation_NSC_combined_cpg_specific_R
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_annotation_NSC_combined_cpg_specific_R.err"
#SBATCH --output="logs/peaks_annotation_NSC_combined_cpg_specific_R.out"

set -e  # Exit on error
set -x  # Print commands as they're executed

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

mkdir -p logs

SCRIPT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/scripts/peaks_annotation"
PEAKS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results_2_align2_005/peaks/narrow"
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results_5_align2_005/peaks_annotation_NSC_combined_cpg_specific_R"

# Remove OUTPUT_DIR if it exists and recreate it
rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Run R script with explicit argument names
Rscript "${SCRIPT_DIR}/peaks_annotation_NSC_combined_cpg_specific_R.R" \
    --peaks-dir "${PEAKS_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    2>&1 | tee "logs/peaks_annotation_NSC_combined_cpg_specific_R.out" 
