#!/bin/bash
#SBATCH --job-name=4create_comparison_plots
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cross_analysis/4create_comparison_plots.err"
#SBATCH --output="logs/cross_analysis/4create_comparison_plots.out"

# Base directories
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/Cross_final"
SCRIPT_DIR="${WORKING_DIR}/src"

# Create necessary directories
mkdir -p "${BASE_DIR}/logs/cross_analysis"
mkdir -p "${WORKING_DIR}/results/4create_comparison_plots"

# Change to working directory
cd $WORKING_DIR || exit 1

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

echo "Starting comparison plot analysis..."
date

# Run the Python analysis script
python ${SCRIPT_DIR}/4create_comparison_plots.py \
    2>&1 | tee -a "${BASE_DIR}/logs/cross_analysis/comparison_plots.log"

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "Analysis completed successfully"
    date
    echo "Results are available in ${WORKING_DIR}/results/4create_comparison_plots/"
    echo "Generated files:"
    ls -l "${WORKING_DIR}/results/4create_comparison_plots/"
else
    echo "Analysis failed with error code $?"
    date
    exit 1
fi
