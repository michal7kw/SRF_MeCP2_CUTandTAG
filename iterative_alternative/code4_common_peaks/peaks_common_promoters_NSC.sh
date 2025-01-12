#!/bin/bash
#SBATCH --job-name=peaks_common_promoters_NSC
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_common_promoters_NSC.err"
#SBATCH --output="logs/peaks_common_promoters_NSC.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

mkdir -p logs

# Define base directory
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"

# Define directories relative to base
WORKING_DIR="${BASE_DIR}/iterative_alternative"
DATA_DIR="${WORKING_DIR}/results_1b"
PEAKS_DIR="${WORKING_DIR}/results_2_align2_005/peaks/narrow"
OUTPUT_DIR="${WORKING_DIR}/results_5_align2_005/peaks_common_promoters_NSC"
GTF_PATH="${BASE_DIR}/DATA/gencode.vM10.annotation.gtf"

echo "Creating output directory..."
mkdir -p "${OUTPUT_DIR}"

echo "Running gene annotation analysis..."
python -u ../scripts/peaks_common/peaks_common_promoters_NSC.py \
    --gtf-path "${GTF_PATH}" \
    --peaks-dir "${PEAKS_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    2>&1 | tee "logs/peaks_common_promoters_NSC.out"

echo "Analysis complete. Results are in ${OUTPUT_DIR}"
