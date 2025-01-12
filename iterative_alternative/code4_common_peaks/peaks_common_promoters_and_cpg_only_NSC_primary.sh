#!/bin/bash
#SBATCH --job-name=peaks_common_promoters_and_cpg_only_NSC_primary
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_common_promoters_and_cpg_only_NSC_primary.err"
#SBATCH --output="logs/peaks_common_promoters_and_cpg_only_NSC_primary.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

mkdir -p logs

# Base paths
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"

# Input paths
PEAKS_DIR="${WORKING_DIR}/results_2_align2_005/peaks/narrow"

# Output paths 
RESULTS_DIR="${WORKING_DIR}/results_5_align2_005/peaks_common_promoters_and_cpg_only_NSC_primary"

# Create output directories
echo "Creating directory structure..."
rm -rf "${RESULTS_DIR}"
mkdir -p "${RESULTS_DIR}"

# Run analysis script
python -u ../scripts/peaks_common/peaks_common_promoters_and_cpg_only_NSC_primary.py \
    --working-dir "${WORKING_DIR}" \
    --results-dir "${RESULTS_DIR}" \
    --peaks-dir "${PEAKS_DIR}" \
    2>&1 | tee "logs/peaks_common_promoters_and_cpg_only_NSC_primary.out"
