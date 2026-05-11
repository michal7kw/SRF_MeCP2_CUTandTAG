#!/bin/bash
#SBATCH --job-name=7analyze_smarcb1_comparison_1
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cross_analysis/7analyze_smarcb1_comparison_1.err"
#SBATCH --output="logs/cross_analysis/7analyze_smarcb1_comparison_1.out"

# Base directories
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/Cross_final"
SCRIPT_DIR="${WORKING_DIR}/src"

# Create logs directory
mkdir -p "${BASE_DIR}/logs/cross_analysis"
mkdir -p "${WORKING_DIR}/results/7analyze_smarcb1_comparison_1"

cd $WORKING_DIR || exit 1

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

echo "Starting meta-profile analysis..."

# Run the Python analysis script
python ${SCRIPT_DIR}/7analyze_smarcb1_comparison_1.py \
    2>&1 | tee -a "${BASE_DIR}/logs/cross_analysis/profile_analysis.log"

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "Analysis completed successfully"
    echo "Results are available in ${WORKING_DIR}/results/7analyze_smarcb1_comparison_1/"
    echo "Generated files:"
    ls -l "${WORKING_DIR}/results/7analyze_smarcb1_comparison_1/"
else
    echo "Analysis failed with error code $?"
    exit 1
fi
