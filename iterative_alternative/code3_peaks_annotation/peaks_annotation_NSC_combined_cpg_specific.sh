#!/bin/bash
#SBATCH --job-name=peaks_annotation_NSC_combined_cpg_specific
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_annotation_NSC_combined_cpg_specific.err"
#SBATCH --output="logs/peaks_annotation_NSC_combined_cpg_specific.out"

# Base directories
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
PROJ_DIR="${BASE_DIR}/iterative_alternative"

cd ${PROJ_DIR}

# Environment setup
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

mkdir -p logs

# Input data paths
DATA_DIR="${BASE_DIR}/DATA"
GTF_PATH="${DATA_DIR}/gencode.vM10.annotation.gtf"
CPG_PATH="${DATA_DIR}/cpgIslandExt.txt"

# Analysis directories
PEAKS_DIR="${PROJ_DIR}/results_2_align2_005/peaks/narrow"
OUTPUT_DIR="${PROJ_DIR}/results_5_align2_005/peaks_annotation_NSC_combined_cpg_specific"

# Remove OUTPUT_DIR if it exists and recreate it
rm -rf "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

python -u ../scripts/peaks_annotation/peaks_annotation_NSC_combined_cpg_specific.py \
    --gtf-path "$GTF_PATH" \
    --peaks-dir "$PEAKS_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --cpg-path "$CPG_PATH" \
    2>&1 | tee "logs/peaks_annotation_NSC_combined_cpg_specific.out"
