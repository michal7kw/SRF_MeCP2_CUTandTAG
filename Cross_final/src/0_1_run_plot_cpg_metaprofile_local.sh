#!/bin/bash

# Check for --force argument
FORCE_FLAG=""
if [[ " $@ " =~ " --force " ]]; then
  FORCE_FLAG="--force"
  echo "Force recalculation requested."
fi

# --- !!! ADJUST THESE PATHS FOR YOUR LOCAL SETUP !!! ---
# Base project directory
BASE_DIR="/mnt/d/Github/SRF_MeCP2_cut_tag" # Example: Update with your base project directory
WORKING_DIR="${BASE_DIR}/Cross_final"
SCRIPT_DIR="${WORKING_DIR}/src"
# CONDA_ACTIVATION_CMD="source /path/to/your/miniconda3/bin/activate your_env_name" # Example: Update with your conda path and environment name

# Input files
CPG_ISLANDS_BED="${BASE_DIR}/data/cpg_islands.bed" # Path to your CpG islands BED file
# --- Directory containing BigWig files to process --- 
# BIGWIG_DIR="${BASE_DIR}/data/bigwig" 
BIGWIG_DIR="${WORKING_DIR}/results_local/normalized_bigwig"

# Output directory
OUTPUT_DIR="${WORKING_DIR}/results_local/3cpg_metaprofiles"

# Other settings
NUM_THREADS=8 # Adjust number of threads as needed
OUTPUT_PREFIX="$(basename ${BIGWIG_FILE} .bw)" # Optional: Prefix for output files based on BW name
# --- !!! END OF PATH ADJUSTMENTS !!! ---

# Define log directory and file
LOG_DIR="${WORKING_DIR}/logs_local/cpg_metaprofiles"
# Log file will be generated per job below

# Create log and output directories
mkdir -p "${LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Activate conda environment
# echo "Activating conda environment..."
# eval ${CONDA_ACTIVATION_CMD}
# if [ $? -ne 0 ]; then
#     echo "Failed to activate conda environment. Check CONDA_ACTIVATION_CMD in the script."
#     exit 1
# fi
# echo "Conda environment activated."

# Check if input files/directories exist
if [ ! -f "${CPG_ISLANDS_BED}" ]; then
    echo "Error: CpG islands BED file not found at ${CPG_ISLANDS_BED}"
    exit 1
fi
if [ ! -d "${BIGWIG_DIR}" ]; then
    echo "Error: BigWig directory not found at ${BIGWIG_DIR}"
    exit 1
fi

echo "Starting parallel CpG meta-profile generation for *.bw files in: ${BIGWIG_DIR}"
echo "Maximum parallel jobs: ${NUM_THREADS}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Log directory: ${LOG_DIR}"

# Export variables needed by the subshell command in xargs
export SCRIPT_DIR CPG_ISLANDS_BED OUTPUT_DIR FORCE_FLAG LOG_DIR

# Function to process a single BigWig file
process_bw_file() {
    local bigwig_file="$1"
    local bw_basename=$(basename "${bigwig_file}" .bw)
    local output_prefix="${bw_basename}" 
    local log_file="${LOG_DIR}/3run_plot_${bw_basename}.log"

    echo "Processing: ${bw_basename} -> Log: ${log_file}"

    python -u "${SCRIPT_DIR}/0_1_plot_cpg_metaprofile.py" \
        --cpg-bed "${CPG_ISLANDS_BED}" \
        --bigwig "${bigwig_file}" \
        --output-dir "${OUTPUT_DIR}" \
        --threads 4 \
        --window-size 5000 \
        --bin-size 300 \
        --output-prefix "${output_prefix}" \
        ${FORCE_FLAG} \
        >"${log_file}" 2>&1 # Redirect stdout and stderr to the specific log file

    # Check exit code and report status
    local exit_code=$?
    if [ ${exit_code} -eq 0 ]; then
        echo "Successfully processed: ${bw_basename}"
    else
        echo "Error processing ${bw_basename} (Exit Code: ${exit_code}). Check log: ${log_file}"
    fi
    return ${exit_code}
}

# Export the function so xargs can use it
export -f process_bw_file

# Find all .bw files and process them in parallel using xargs
# -print0 and -0 handle filenames with spaces or special characters
find "${BIGWIG_DIR}" -maxdepth 1 -name '*.bw' -print0 | xargs -0 -n 1 -P "${NUM_THREADS}" -I {} bash -c 'process_bw_file "$@"' _ {}

# Note: Error checking for xargs itself is complex. Individual job errors are reported by process_bw_file.
# You might want to check the logs directory for any non-zero exit codes reported above.

echo "Parallel processing script finished. Check individual logs in ${LOG_DIR} for details and potential errors."

# Old single file execution removed
# echo "Starting CpG meta-profile generation for: $(basename ${BIGWIG_FILE})"
# echo "Logging to: ${LOG_FILE}"
# 
# # Run the Python script
# python -u ${SCRIPT_DIR}/3plot_cpg_metaprofile.py \
#     --cpg-bed "${CPG_ISLANDS_BED}" \
#     --bigwig "${BIGWIG_FILE}" \
#     --output-dir "${OUTPUT_DIR}" \
#     --threads 16 \
#     --window-size 5000 \
#     --bin-size 300 \
#     --output-prefix "${OUTPUT_PREFIX}" \
#     ${FORCE_FLAG} \
#     2>&1 | tee -a "${LOG_FILE}"
# 
# # Check if the script completed successfully
# SCRIPT_EXIT_CODE=${PIPESTATUS[0]}
# if [ ${SCRIPT_EXIT_CODE} -eq 0 ]; then
#     echo "Meta-profile generation completed successfully for $(basename ${BIGWIG_FILE})"
#     echo "Output plot/data should be in: ${OUTPUT_DIR}"
#     echo "Check the log file for details: ${LOG_FILE}"
# else
#     echo "Meta-profile generation failed with error code ${SCRIPT_EXIT_CODE} for $(basename ${BIGWIG_FILE})"
#     echo "Check the log file for errors: ${LOG_FILE}"
#     exit 1
# fi
# 
# echo "Local script finished." 