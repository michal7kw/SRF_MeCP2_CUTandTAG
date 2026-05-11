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

# Base paths
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
GTF_PATH="${BASE_DIR}/DATA/gencode.vM10.annotation.gtf"

# Process both narrow and broad peaks
for PEAK_TYPE in narrow broad; do
    echo "Processing ${PEAK_TYPE} peaks..."
    
    # Input paths
    PEAKS_DIR="${WORKING_DIR}/results_2_align2_005/peaks/${PEAK_TYPE}"

    # Output paths 
    OUTPUT_DIR="${WORKING_DIR}/results_5_align2_005/peaks_common_promoters_NSC_${PEAK_TYPE}"

    # Create output directories
    echo "Creating directory structure for ${PEAK_TYPE} peaks..."
    rm -rf "${OUTPUT_DIR}"
    mkdir -p "${OUTPUT_DIR}"

    # Run analysis script
    python -u ../scripts/peaks_common/peaks_common_promoters_NSC.py \
        --gtf-path "${GTF_PATH}" \
        --peaks-dir "${PEAKS_DIR}" \
        --output-dir "${OUTPUT_DIR}" \
        --peak-type "${PEAK_TYPE}" \
        2>&1 | tee "logs/peaks_common_promoters_NSC_${PEAK_TYPE}.out"

    echo "Analysis complete for ${PEAK_TYPE} peaks. Results are in ${OUTPUT_DIR}"
done
