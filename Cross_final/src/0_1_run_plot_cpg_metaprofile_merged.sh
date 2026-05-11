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
# --- Directory containing BigWig files ---
BIGWIG_DIR="${BASE_DIR}/data/bigwig/merged"

# --- Specific BigWig files to process ---
SPECIFIC_BW_FILES=(
    "Medip_DP_output_merged.bw"
    "Medip_N_output_merged.bw"
    "Medip_PP_output_merged.bw"
)

# Output directory (modified for specific files)
OUTPUT_DIR="${WORKING_DIR}/results_local/3cpg_metaprofiles_medip"

# Other settings
NUM_THREADS_PER_JOB=5 # Threads for each python script instance
# --- !!! END OF PATH ADJUSTMENTS !!! ---

# Define log directory and file (modified for specific files)
LOG_DIR="${WORKING_DIR}/logs_local/cpg_metaprofiles_medip"
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

echo "Starting CpG meta-profile generation for specific Medip files:"
printf " - %s\\n" "${SPECIFIC_BW_FILES[@]}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Log directory: ${LOG_DIR}"

# Export variables needed by the subshell command in process_bw_file
export SCRIPT_DIR CPG_ISLANDS_BED OUTPUT_DIR FORCE_FLAG LOG_DIR NUM_THREADS_PER_JOB

# Function to process a single BigWig file
process_bw_file() {
    local bw_filename="$1"
    local bigwig_file="${BIGWIG_DIR}/${bw_filename}" # Construct full path
    local bw_basename=$(basename "${bigwig_file}" .bw)
    local output_prefix="${bw_basename}"
    local log_file="${LOG_DIR}/3run_plot_${bw_basename}.log"

    # Check if the specific BigWig file exists
    if [ ! -f "${bigwig_file}" ]; then
        echo "Error: BigWig file not found at ${bigwig_file}"
        return 1 # Return error code
    fi

    echo "Processing: ${bw_basename} -> Log: ${log_file}"

    python -u "${SCRIPT_DIR}/0_1_plot_cpg_metaprofile.py" \
        --cpg-bed "${CPG_ISLANDS_BED}" \
        --bigwig "${bigwig_file}" \
        --output-dir "${OUTPUT_DIR}" \
        --threads ${NUM_THREADS_PER_JOB} \
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

# Export the function so background processes can use it
export -f process_bw_file

# Process the specific files in parallel
PIDS=() # Array to store background process IDs
EXIT_CODES=() # Array to store exit codes

for bw_file in "${SPECIFIC_BW_FILES[@]}"; do
    # Run processing in the background
    process_bw_file "${bw_file}" &
    PIDS+=($!) # Store the PID of the background process
done

# Wait for all background jobs to finish and collect exit codes
all_success=true
for pid in "${PIDS[@]}"; do
    wait $pid
    exit_code=$?
    EXIT_CODES+=(${exit_code})
    if [ ${exit_code} -ne 0 ]; then
        all_success=false
    fi
done

echo "----------------------------------------"
echo "Summary of processing:"
for i in ${!SPECIFIC_BW_FILES[@]}; do
    bw_file="${SPECIFIC_BW_FILES[$i]}"
    exit_code="${EXIT_CODES[$i]}"
    if [ ${exit_code} -eq 0 ]; then
        echo " - ${bw_file}: Success"
    else
        echo " - ${bw_file}: Failed (Exit Code: ${exit_code})"
    fi
done
echo "----------------------------------------"

if $all_success; then
    echo "All specific Medip files processed successfully."
    echo "Output plots/data should be in: ${OUTPUT_DIR}"
    echo "Check individual logs in ${LOG_DIR} for details."
else
    echo "One or more files failed during processing."
    echo "Check the logs in ${LOG_DIR} for errors."
    # Optionally exit with an error code if any job failed
    # exit 1
fi

echo "Specific Medip processing script finished." 