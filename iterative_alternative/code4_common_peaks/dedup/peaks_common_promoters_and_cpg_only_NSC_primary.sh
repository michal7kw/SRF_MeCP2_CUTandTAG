#!/bin/bash
#SBATCH --job-name=pcp_dedup
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_common/dedup/peaks_common_promoters_and_cpg_only_NSC_primary.err"
#SBATCH --output="logs/peaks_common/dedup/peaks_common_promoters_and_cpg_only_NSC_primary.out"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
SCRIPT_DIR="${BASE_DIR}/scripts/peaks_common"
cd $WORKING_DIR || exit 1

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Base paths
INPUT_DIR="results/dedup/peaks"
RESULTS_DIR="results/dedup/peaks_common/primary_tss"

mkdir -p "${RESULTS_DIR}"

# Process both narrow and broad peaks
for PEAK_TYPE in narrow broad; do
    echo "Processing ${PEAK_TYPE} peaks..."
    
    # Input paths
    PEAKS_DIR="${INPUT_DIR}/${PEAK_TYPE}"

    # Run analysis script
    python -u "${SCRIPT_DIR}/peaks_common_promoters_and_cpg_only_NSC_primary.py" \
        --working-dir "${WORKING_DIR}" \
        --results-dir "${RESULTS_DIR}" \
        --peaks-dir "${PEAKS_DIR}" \
        --peak-type "${PEAK_TYPE}" \
        2>&1 | tee "logs/peaks_common/dedup/peaks_common_promoters_and_cpg_only_NSC_primary_${PEAK_TYPE}.out"
done
