#!/bin/bash
#SBATCH --job-name=peaks_annotation_NSC_combined_R
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_annotation_NSC_combined_R.err"
#SBATCH --output="logs/peaks_annotation_NSC_combined_R.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create logs directory
mkdir -p logs

# Base directories
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
PROJ_DIR="${BASE_DIR}/iterative_alternative"

# Analysis directories
PEAKS_DIR="${PROJ_DIR}/results_2_align2_005/peaks/narrow"
OUTPUT_DIR="${PROJ_DIR}/results_5_align2_005/peaks_annotation_NSC_combined_R"

# Remove OUTPUT_DIR if it exists and recreate it
rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Run the R script with quoted arguments
Rscript ../scripts/peaks_annotation/peaks_annotation_NSC_combined_R.R \
    --peaks-dir "${PEAKS_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    2>&1 | tee "logs/peaks_annotation_NSC_combined_R.out" 
