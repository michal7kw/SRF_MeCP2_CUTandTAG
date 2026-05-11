#!/bin/bash
#SBATCH --job-name=analyze_cpg_genes
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-7
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment/nsc/no_dedup/individual_samples/broad/cpg_enrichment_%a.err"
#SBATCH --output="logs/cpg_enrichment/nsc/no_dedup/individual_samples/broad/cpg_enrichment_%a.out"

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

cd $WORKING_DIR || exit 1

RESULTS_DIR="results/${ALIGNMENT_TYPE}/cpg_enrichment/${CELL_LINE}/${PEAKS}/${RUN_NAME}"

# Convert the multi-line string to an array of analysis directories
readarray -t DATA_ARRAY <<EOF
up_enriched_signal_1
up_enriched_signal_1_5
up_enriched_signal_2
down_enriched_signal_1
down_enriched_signal_08
down_enriched_signal_05
exo_only_df_by_signal
endo_only_df_by_signal
EOF

# Get the current directory to process based on array task ID
CURRENT_DIR="${DATA_ARRAY[$SLURM_ARRAY_TASK_ID]}"

echo "Processing directory: ${CURRENT_DIR}"

# Create the full results directory path
ANALYSIS_DIR="${RESULTS_DIR}/peaks_annotation/${CURRENT_DIR}"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

Rscript "${SCRIPT_DIR}/analyze_cpg_genes.R" \
    --work-dir "${ANALYSIS_DIR}"
