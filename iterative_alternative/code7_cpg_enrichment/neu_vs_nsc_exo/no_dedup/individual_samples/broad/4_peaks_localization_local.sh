#!/bin/sh
# This script is adapted for local execution.

# Set BASE_DIR relative to the script's location or assume it's run from the project root.
# Assuming the script is run from d:/Github/SRF_MeCP2_cut_tag
BASE_DIR="d:/Github/SRF_MeCP2_cut_tag" # Use correct project root path and case
CONFIGURATION_NUMBER=10

# Load parameters from config.yaml
# Ensure config.yaml path is correct relative to BASE_DIR
CONFIG_FILE="${BASE_DIR}/iterative_alternative/code7_cpg_enrichment/config.yaml"

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file not found at ${CONFIG_FILE}"
    exit 1
fi

# Using awk to get specific fields based on document number
# Note: Requires awk to be available in your local environment (e.g., Git Bash, WSL)
DOC_NUM=$CONFIGURATION_NUMBER
ALIGNMENT_TYPE=$(awk -v doc="$DOC_NUM" 'BEGIN{doc_count=0} /^---/{doc_count++} doc_count==doc && /alignment_type:/{print $2}' "$CONFIG_FILE")
PEAKS=$(awk -v doc="$DOC_NUM" 'BEGIN{doc_count=0} /^---/{doc_count++} doc_count==doc && /peaks_type:/{print $2}' "$CONFIG_FILE")
RUN_NAME=$(awk -v doc="$DOC_NUM" 'BEGIN{doc_count=0} /^---/{doc_count++} doc_count==doc && /run_name:/{print $2}' "$CONFIG_FILE")

# Remove any quotes from the values
ALIGNMENT_TYPE=$(echo "$ALIGNMENT_TYPE" | tr -d '"')
PEAKS=$(echo "$PEAKS" | tr -d '"')
RUN_NAME=$(echo "$RUN_NAME" | tr -d '"')

# Print parameters for verification
echo "Parameters loaded from config file:"
echo "BASE_DIR: ${BASE_DIR}"
echo "ALIGNMENT_TYPE: ${ALIGNMENT_TYPE}"
echo "PEAKS: ${PEAKS}"
echo "RUN_NAME: ${RUN_NAME}"

# Define directories relative to BASE_DIR
WORKING_DIR="${BASE_DIR}/iterative_alternative"
SCRIPT_DIR="${BASE_DIR}/scripts/cpg_enrichment" # Adjusted path for R script
RESULTS_DIR="${WORKING_DIR}/results/${ALIGNMENT_TYPE}/cpg_enrichment/${RUN_NAME}" # Adjusted path relative to WORKING_DIR

echo "Derived parameters:"
echo "WORKING_DIR: ${WORKING_DIR}"
echo "SCRIPT_DIR: ${SCRIPT_DIR}"
echo "RESULTS_DIR: ${RESULTS_DIR}"

# Verify R script exists
if [ ! -f "${SCRIPT_DIR}/peaks_localization_neu_vs_nsc_endo_tune.R" ]; then
    echo "Error: R script not found at ${SCRIPT_DIR}/peaks_localization_neu_vs_nsc_endo_tune.R"
    exit 1
fi

# Navigate to the working directory relative to the script's execution location
cd "$WORKING_DIR" || { echo "Error: Failed to change directory to ${WORKING_DIR}"; exit 1; }

# Activate your local conda environment
# Replace 'snakemake' if your environment has a different name
# This command might need adjustment based on your shell (e.g., use 'conda activate snakemake' directly in bash/zsh)
# echo "Attempting to activate conda environment 'snakemake'..."
# If using standard cmd or PowerShell, conda init might be needed first.
# For Git Bash or WSL, 'conda activate' should work if conda is in PATH.
# conda activate snakemake || { echo "Error: Failed to activate conda environment 'snakemake'. Make sure conda is initialized and the environment exists."; exit 1; }


# Define the array of files to process
DATA_ARRAY=(
"up_enriched_signal_1.csv"
"up_enriched_signal_1_5.csv"
"up_enriched_signal_2.csv"
"down_enriched_signal_1.csv"
"down_enriched_signal_08.csv"
"down_enriched_signal_05.csv"
"neu_only_df_by_signal.csv"
"nsc_only_df_by_signal.csv"
"endo_list_neu.csv"
"endo_list_nsc.csv"
)

# Loop through the files and process each one
for CURRENT_FILE in "${DATA_ARRAY[@]}"; do
    echo "----------------------------------------"
    echo "Processing file: ${CURRENT_FILE}"
    echo "----------------------------------------"

    # Verify input file exists before processing
    INPUT_FILE_PATH="${RESULTS_DIR}/lists/${CURRENT_FILE}"
    if [ ! -f "${INPUT_FILE_PATH}" ]; then
        echo "Warning: Input file not found at ${INPUT_FILE_PATH}. Skipping..."
        continue # Skip to the next file
    fi

    # Ensure Rscript is in PATH or provide full path
    # Use paths relative to the current directory (which is now WORKING_DIR)
    Rscript "${SCRIPT_DIR}/peaks_localization_neu_vs_nsc_endo_tune.R" \
        --work-dir "${RESULTS_DIR}" \
        --data-to-analyze "${CURRENT_FILE}"

    # Check Rscript exit status
    if [ $? -ne 0 ]; then
        echo "Error: R script failed for file ${CURRENT_FILE}"
        # Decide if you want to stop or continue with the next file
        # exit 1 # Uncomment to stop on first error
    fi
done

echo "----------------------------------------"
echo "Script finished."
echo "----------------------------------------"

# Deactivate conda environment (optional)
# conda deactivate
