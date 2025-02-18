#!/bin/bash
#SBATCH --job-name=peaks_loc_comb
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment/nsc/no_dedup/individual_samples/broad/log_combined.err"
#SBATCH --output="logs/cpg_enrichment/nsc/no_dedup/individual_samples/broad/log_combined.out"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
CONFIGURATION_NUMBER=4

# Load parameters from config.yaml
CONFIG_FILE="${BASE_DIR}/iterative_alternative/code7_cpg_enrichment/config.yaml"

# Using awk to get specific fields based on document number
DOC_NUM=$CONFIGURATION_NUMBER
CELL_LINE=$(awk -v doc="$DOC_NUM" 'BEGIN{doc_count=0} /^---/{doc_count++} doc_count==doc && /cell_type:/{print $2}' "$CONFIG_FILE")
ALIGNMENT_TYPE=$(awk -v doc="$DOC_NUM" 'BEGIN{doc_count=0} /^---/{doc_count++} doc_count==doc && /alignment_type:/{print $2}' "$CONFIG_FILE")
PEAKS=$(awk -v doc="$DOC_NUM" 'BEGIN{doc_count=0} /^---/{doc_count++} doc_count==doc && /peaks_type:/{print $2}' "$CONFIG_FILE")
RUN_NAME=$(awk -v doc="$DOC_NUM" 'BEGIN{doc_count=0} /^---/{doc_count++} doc_count==doc && /run_name:/{print $2}' "$CONFIG_FILE")

# Remove any quotes from the values
CELL_LINE=$(echo "$CELL_LINE" | tr -d '"')
ALIGNMENT_TYPE=$(echo "$ALIGNMENT_TYPE" | tr -d '"')
PEAKS=$(echo "$PEAKS" | tr -d '"')
RUN_NAME=$(echo "$RUN_NAME" | tr -d '"')

# Print parameters for verification
echo "Parameters loaded from config file:"
echo "BASE_DIR: ${BASE_DIR}"
echo "CELL_LINE: ${CELL_LINE}"
echo "ALIGNMENT_TYPE: ${ALIGNMENT_TYPE}"
echo "PEAKS: ${PEAKS}"
echo "RUN_NAME: ${RUN_NAME}"

WORKING_DIR="${BASE_DIR}/iterative_alternative"
SCRIPT_DIR="${BASE_DIR}/scripts/cpg_enrichment"

echo "Derived parameters:"
echo "WORKING_DIR: ${WORKING_DIR}"
echo "SCRIPT_DIR: ${SCRIPT_DIR}"

cd $WORKING_DIR || exit 1

RESULTS_DIR="results/${ALIGNMENT_TYPE}/cpg_enrichment/${CELL_LINE}/${PEAKS}/${RUN_NAME}"
echo "RESULTS_DIR: ${RESULTS_DIR}"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

Rscript "${SCRIPT_DIR}/peaks_localization_combined.R" \
    --work-dir "${RESULTS_DIR}" \
    --data-file "cpg_enrichment_parallel.csv" 