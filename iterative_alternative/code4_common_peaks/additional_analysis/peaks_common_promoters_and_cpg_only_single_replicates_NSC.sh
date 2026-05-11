#!/bin/bash
#SBATCH --job-name=peaks_common_promoters_and_cpg_only_single_replicates_NSC
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_common_promoters_and_cpg_only_single_replicates_NSC.err"
#SBATCH --output="logs/peaks_common_promoters_and_cpg_only_single_replicates_NSC.out"

# Base directory structure
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"

# Input directories
INPUT_DIR="${WORKING_DIR}/results_2_align2_005"
PEAKS_INPUT_DIR="${INPUT_DIR}/peaks/narrow"

# Output directories
RESULTS_DIR="${WORKING_DIR}/results_5_align2_005/peaks_common_promoters_and_cpg_only_single_replicates_NSC"
PEAKS_OUTPUT_DIR="${RESULTS_DIR}/peaks"

# Create required directories
cd ${WORKING_DIR}
mkdir -p logs
rm -rf "${RESULTS_DIR}"
mkdir -p "${PEAKS_OUTPUT_DIR}"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

echo "Creating directory structure..."

echo "Copying and renaming peaks from ${PEAKS_INPUT_DIR} to ${PEAKS_OUTPUT_DIR}..."

# Copy and rename peak files
for sample in NSCv1 NSCv2 NSCv3 NSCM1 NSCM2 NSCM3; do
    source_file="${PEAKS_INPUT_DIR}/${sample}_narrow_peaks.filtered.narrowPeak"
    target_file="${PEAKS_OUTPUT_DIR}/${sample}_peaks.narrowPeak"
    
    if [ -f "${source_file}" ]; then
        echo "Copying ${source_file} to ${target_file}"
        cp "${source_file}" "${target_file}"
    else
        echo "Warning: Source file not found: ${source_file}"
    fi
done

# Verify directory contents
echo "Contents of source directory:"
ls -l "${PEAKS_INPUT_DIR}"

echo "Contents of target directory:"
ls -l "${PEAKS_OUTPUT_DIR}"

python -u ../scripts/peaks_common/peaks_common_promoters_and_cpg_only_single_replicates_NSC.py \
    --working-dir ${WORKING_DIR} \
    --results-dir ${RESULTS_DIR} \
    2>&1 | tee "logs/peaks_common_promoters_and_cpg_only_single_replicates_NSC.out"
