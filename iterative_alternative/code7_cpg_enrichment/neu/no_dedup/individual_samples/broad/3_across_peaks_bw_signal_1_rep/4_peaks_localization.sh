#!/bin/bash
#SBATCH --job-name=peaks_localization
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-9
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment/neu/no_dedup/individual_samples/broad/log_%a.err"
#SBATCH --output="logs/cpg_enrichment/neu/no_dedup/individual_samples/broad/log_%a.out"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
CONFIGURATION_NUMBER=7

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

# Convert the multi-line string to an array
readarray -t DATA_ARRAY <<EOF
up_enriched_signal_1.csv
up_enriched_signal_1_5.csv
up_enriched_signal_2.csv
down_enriched_signal_1.csv
down_enriched_signal_08.csv
down_enriched_signal_05.csv
exo_only_df_by_signal.csv
endo_only_df_by_signal.csv
EOF

# Get the current file to process based on array task ID
CURRENT_FILE="${DATA_ARRAY[$SLURM_ARRAY_TASK_ID]}"

echo "Processing file: ${CURRENT_FILE}"

Rscript "${SCRIPT_DIR}/peaks_localization.R" \
    --work-dir "${RESULTS_DIR}" \
    --data-to-analyze "${CURRENT_FILE}"

peaks_annotation/


# up_enriched_signal_1.csv
# up_enriched_signal_1_5.csv
# up_enriched_signal_2.csv
# down_enriched_signal_1.csv
# down_enriched_signal_08.csv
# down_enriched_signal_05.csv
# exo_only_df_by_signal.csv
# endo_only_df_by_signal.csv

# up_enriched_peaks_1.csv
# up_enriched_peaks_1_5.csv
# up_enriched_peaks_2.csv
# down_enriched_peaks_1.csv
# down_enriched_peaks_08.csv
# down_enriched_peaks_05.csv
# exo_only_df_by_peaks.csv
# endo_only_df_by_peaks.csv