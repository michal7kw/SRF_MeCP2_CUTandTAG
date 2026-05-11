#!/bin/bash

# --- !!! ADJUST THESE PATHS FOR YOUR LOCAL SETUP !!! ---
# Base directories
BASE_DIR="/path/to/your/local/SRF_MeCP2_CUTandTAG" # Example: Update with your base project directory
WORKING_DIR="${BASE_DIR}/Cross_final"
SCRIPT_DIR="${WORKING_DIR}/src"
CONDA_ACTIVATION_CMD="source /path/to/your/miniconda3/bin/activate your_env_name" # Example: Update with your conda path and environment name
# --- !!! END OF PATH ADJUSTMENTS !!! ---

# Define log and results directories locally
LOG_DIR="${WORKING_DIR}/logs_local/cross_analysis"
RESULTS_DIR="${WORKING_DIR}/results_local/2analyze_methylation_enhanced"

# Create log and results directories
mkdir -p "${LOG_DIR}"
mkdir -p "${RESULTS_DIR}" # Python script might recreate subdirs, but ensure base exists

# Define log file
LOG_FILE="${LOG_DIR}/2analyze_methylation_enhanced_local.log"

# Optional: Change to working directory if needed
# echo "Changing to working directory: ${WORKING_DIR}"
# cd "${WORKING_DIR}" || { echo "Failed to change directory to ${WORKING_DIR}"; exit 1; }

# Activate conda environment
echo "Activating conda environment..."
eval $CONDA_ACTIVATION_CMD
if [ $? -ne 0 ]; then
    echo "Failed to activate conda environment. Check CONDA_ACTIVATION_CMD in the script."
    exit 1
fi
echo "Conda environment activated."

# --- IMPORTANT --- 
# Ensure that a 'config.py' file exists and is accessible by the Python script below.
# This config.py should contain the CORRECT LOCAL PATHS for:
# BASE_DIR, GENOME_FASTA, MEDIP_DIR, SMARCB1_DIR, RESULTS_DIR (pointing to ${RESULTS_DIR} or similar)
# The Python script uses these paths from config.py to find data and save results.
# -----------------

echo "Starting enhanced meta-profile analysis... (Logging to ${LOG_FILE})"

# Run the Python analysis script - assumes it's run from WORKING_DIR or SCRIPT_DIR is correctly pathed
# Using python -u for unbuffered output might help in seeing logs in real-time
python -u ${SCRIPT_DIR}/2analyze_methylation_enhanced.py \
    2>&1 | tee -a "${LOG_FILE}"

# Check if the script completed successfully
SCRIPT_EXIT_CODE=${PIPESTATUS[0]}
if [ ${SCRIPT_EXIT_CODE} -eq 0 ]; then
    echo "Analysis completed successfully"
    echo "Results should be available in the directory specified within your config.py (likely related to ${RESULTS_DIR})"
    echo "Check the log file for details: ${LOG_FILE}"
else
    echo "Analysis failed with error code ${SCRIPT_EXIT_CODE}"
    echo "Check the log file for errors: ${LOG_FILE}"
    exit 1
fi

echo "Local script finished." 