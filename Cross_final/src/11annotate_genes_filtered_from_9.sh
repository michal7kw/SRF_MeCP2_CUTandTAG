#!/bin/bash
#SBATCH --job-name=11annotate_genes_filtered_from_9
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cross_analysis/11annotate_genes_filtered_from_9.err"
#SBATCH --output="logs/cross_analysis/11annotate_genes_filtered_from_9.out"

# Base directories
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/Cross_final"
SCRIPT_DIR="${WORKING_DIR}/src"

# Create logs directory
mkdir -p "${BASE_DIR}/logs/cross_analysis"
mkdir -p "${WORKING_DIR}/results/11annotate_genes_filtered_from_9"
cd $WORKING_DIR || exit 1

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

echo "Starting meta-profile analysis..."

# Run the Python analysis script
Rscript ${SCRIPT_DIR}/11annotate_genes_filtered_from_9.R
