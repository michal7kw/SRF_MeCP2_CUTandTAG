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

# --- Directory containing CpG BED files ---
# ASSUMPTION: Bed files are in a subdirectory. Adjust if needed.
BED_DIR="${BASE_DIR}/data"
# BED_DIR="${WORKING_DIR}/results_local/normalized_bigwig"

# --- CpG BED file paths (must match expected names) ---
NSC_TARGETS_BED="${BED_DIR}/nsc_cpg_targets.bed"
NEU_TARGETS_BED="${BED_DIR}/neu_cpg_targets.bed"
NSC_NOT_TARGETS_BED="${BED_DIR}/nsc_cpg_not_targets.bed"
NEU_NOT_TARGETS_BED="${BED_DIR}/neu_cpg_not_targets.bed"

CPG_BED_FILES=(
    "${NEU_NOT_TARGETS_BED}"
    "${NEU_TARGETS_BED}"
    "${NSC_NOT_TARGETS_BED}"
    "${NSC_TARGETS_BED}"
)

# --- Directory containing BigWig files ---
BIGWIG_DIR="${BASE_DIR}/data/bigwig/merged"

# --- Specific BigWig files to process (Basenames only) ---
SPECIFIC_BW_BASENAMES=(
    "Medip_DP_output_merged"
    "Medip_N_output_merged"
    "Medip_PP_output_merged"
    # "Medip_DP_input_merged"
    # "Medip_N_input_merged"
    # "Medip_PP_input_merged"
)

# Output directory for individual profiles/data
OUTPUT_DIR_INDIVIDUAL="${WORKING_DIR}/results_local/3cpg_metaprofiles_medip_targeted"
# Output directory for comparison plots
OUTPUT_DIR_COMPARISON="${WORKING_DIR}/results_local/4cpg_metaprofiles_medip_targeted_comparison"

# Other settings
NUM_THREADS_PER_JOB=5 # Threads for each python script instance
WINDOW_SIZE=10000 # Set full window size to 10000 bp for +/- 5kb
BIN_SIZE=300 # Bin size for metaprofile calculation

# --- Define Suffixes for Comparison Plotting (Must match python script expectations) ---
# NSC_TARGET_SMOOTHED_SUFFIX="_vs_nsc_cpg_targets_cpg_metaprofile_smoothed.npy"
# NEURON_TARGET_SMOOTHED_SUFFIX="_vs_neu_cpg_targets_cpg_metaprofile_smoothed.npy"

# --- Define Specific Files for Final Comparison --- 
BW_NEURON="Medip_N_output_merged"
BW_NSC="Medip_PP_output_merged"
BED_NEURON_TARGET_BASENAME="neu_cpg_targets"
BED_NSC_TARGET_BASENAME="nsc_cpg_targets"

# Construct the specific smoothed .npy file paths needed for the final comparison
NEURON_SMOOTHED_NPY="${OUTPUT_DIR_INDIVIDUAL}/${BW_NEURON}_vs_${BED_NEURON_TARGET_BASENAME}_cpg_metaprofile_smoothed.npy"
NSC_SMOOTHED_NPY="${OUTPUT_DIR_INDIVIDUAL}/${BW_NSC}_vs_${BED_NSC_TARGET_BASENAME}_cpg_metaprofile_smoothed.npy"
COMPARISON_PLOT_OUTPUT_FILE="${OUTPUT_DIR_COMPARISON}/Specific_MedipN_vs_MedipPP_Target_Comparison.png"

# --- !!! END OF PATH ADJUSTMENTS !!! ---

# Define log directory and file
LOG_DIR_INDIVIDUAL="${WORKING_DIR}/logs_local/3cpg_metaprofiles_medip_targeted_individual"
LOG_DIR_COMPARISON="${WORKING_DIR}/logs_local/4cpg_metaprofiles_medip_targeted_comparison"
# Log file will be generated per job below

# Create log and output directories
mkdir -p "${LOG_DIR_INDIVIDUAL}"
mkdir -p "${LOG_DIR_COMPARISON}"
mkdir -p "${OUTPUT_DIR_INDIVIDUAL}"
mkdir -p "${OUTPUT_DIR_COMPARISON}"

# Activate conda environment
# echo "Activating conda environment..."
# eval ${CONDA_ACTIVATION_CMD}
# if [ $? -ne 0 ]; then
#     echo "Failed to activate conda environment. Check CONDA_ACTIVATION_CMD in the script."
#     exit 1
# fi
# echo "Conda environment activated."

# Check if input files/directories exist
if [ ! -d "${BIGWIG_DIR}" ]; then
    echo "Error: BigWig directory not found at ${BIGWIG_DIR}"
    exit 1
fi
for cpg_bed_file in "${CPG_BED_FILES[@]}"; do
    if [ ! -f "${cpg_bed_file}" ]; then
        echo "Error: CpG BED file not found at ${cpg_bed_file}"
        exit 1
    fi
done
for bw_base in "${SPECIFIC_BW_BASENAMES[@]}"; do
    bw_file="${BIGWIG_DIR}/${bw_base}.bw"
    if [ ! -f "${bw_file}" ]; then
        echo "Error: BigWig file not found at ${bw_file}"
        exit 1
    fi
done

echo "Starting CpG meta-profile generation for:"
echo " BigWig files (basenames):"
printf "  - %s\n" "${SPECIFIC_BW_BASENAMES[@]}"
echo " CpG BED files:"
printf "  - %s\n" "${CPG_BED_FILES[@]}"
echo "Output directory (Individual): ${OUTPUT_DIR_INDIVIDUAL}"
echo "Output directory (Comparison): ${OUTPUT_DIR_COMPARISON}"
echo "Log directory (Individual): ${LOG_DIR_INDIVIDUAL}"
echo "Log directory (Comparison): ${LOG_DIR_COMPARISON}"

# Export variables needed by the subshell command in process_bw_file
export SCRIPT_DIR OUTPUT_DIR_INDIVIDUAL FORCE_FLAG LOG_DIR_INDIVIDUAL NUM_THREADS_PER_JOB BIGWIG_DIR WINDOW_SIZE BIN_SIZE

# Function to process a single BigWig file against a single CpG BED file (Individual Calculation)
process_bw_cpg_pair_individual() {
    local bw_basename="$1"
    local cpg_bed_path="$2"
    local bigwig_file="${BIGWIG_DIR}/${bw_basename}.bw" # Construct full path
    local cpg_basename=$(basename "${cpg_bed_path}" .bed)
    local output_prefix="${bw_basename}_vs_${cpg_basename}"
    local log_file="${LOG_DIR_INDIVIDUAL}/3run_plot_${output_prefix}.log"

    # Check if the specific BigWig file exists (already checked globally, but good practice)
    if [ ! -f "${bigwig_file}" ]; then
        echo "Error (Worker): BigWig file not found at ${bigwig_file}. Skipping job."
        return 1 # Return error code
    fi
     # Check if the specific CpG BED file exists (already checked globally, but good practice)
    if [ ! -f "${cpg_bed_path}" ]; then
        echo "Error (Worker): CpG BED file not found at ${cpg_bed_path}. Skipping job."
        return 1 # Return error code
    fi

    echo "Processing Individual: ${bw_basename} vs ${cpg_basename} -> Log: ${log_file}"

    # --- Call the calculation/plotting mode of the python script --- 
    python -u "${SCRIPT_DIR}/0_1_plot_cpg_metaprofile.py" \
        --cpg-bed "${cpg_bed_path}" \
        --bigwig "${bigwig_file}" \
        --output-dir "${OUTPUT_DIR_INDIVIDUAL}" \
        --threads ${NUM_THREADS_PER_JOB} \
        --window-size $((WINDOW_SIZE / 2)) \
        --bin-size ${BIN_SIZE} \
        --output-prefix "${output_prefix}" \
        ${FORCE_FLAG} \
        >"${log_file}" 2>&1 # Redirect stdout and stderr to the specific log file

    # Check exit code and report status
    local exit_code=$?
    if [ ${exit_code} -eq 0 ]; then
        echo "Successfully processed individual: ${output_prefix}"
    else
        echo "Error processing individual ${output_prefix} (Exit Code: ${exit_code}). Check log: ${log_file}"
    fi
    return ${exit_code}
}

# Export the function so background processes can use it
export -f process_bw_cpg_pair_individual

# --- Run Individual Calculations First --- 
PIDS_INDIVIDUAL=() # Array to store background process IDs
JOB_INFO_INDIVIDUAL=() # Array to store job descriptions (BW vs CpG)
EXIT_CODES_INDIVIDUAL=() # Array to store exit codes

for bw_base in "${SPECIFIC_BW_BASENAMES[@]}"; do
    for cpg_bed_file_path in "${CPG_BED_FILES[@]}"; do
        # Run processing in the background
        process_bw_cpg_pair_individual "${bw_base}" "${cpg_bed_file_path}" & 
        pid=$!
        PIDS_INDIVIDUAL+=(${pid}) # Store the PID of the background process
        JOB_INFO_INDIVIDUAL+=("${bw_base} vs $(basename ${cpg_bed_file_path})") # Store job description
    done
done

# Wait for all individual background jobs to finish and collect exit codes
all_individual_success=true
failed_individual_jobs=()
echo "Waiting for individual calculation jobs to complete..."
for i in ${!PIDS_INDIVIDUAL[@]}; do
    pid=${PIDS_INDIVIDUAL[$i]}
    job_desc=${JOB_INFO_INDIVIDUAL[$i]}
    wait ${pid}
    exit_code=$?
    EXIT_CODES_INDIVIDUAL+=(${exit_code})
    if [ ${exit_code} -ne 0 ]; then
        all_individual_success=false
        failed_individual_jobs+=("  - ${job_desc} (Exit Code: ${exit_code})")
    fi
done

echo "----------------------------------------"
echo "Summary of Individual Processing:"
if $all_individual_success; then
    echo " All individual jobs completed successfully."
else
    echo " Errors occurred in individual processing:"
    printf '%s\n' "${failed_individual_jobs[@]}"
    echo " Individual logs in: ${LOG_DIR_INDIVIDUAL}"
fi
echo "----------------------------------------"

# --- Run Comparison Plotting --- 
if ! $all_individual_success; then
    echo "Skipping comparison plotting due to errors in individual processing."
else
    echo "Proceeding with specific comparison plotting..."
    # PIDS_COMPARISON=() # Array to store background process IDs - Not needed for single job
    # JOB_INFO_COMPARISON=() # Array to store job descriptions (BW Comparison) - Not needed
    # EXIT_CODES_COMPARISON=() # Array to store exit codes - Just need one

    # Export comparison-specific variables - Simpler now
    export SCRIPT_DIR OUTPUT_DIR_INDIVIDUAL OUTPUT_DIR_COMPARISON FORCE_FLAG LOG_DIR_COMPARISON WINDOW_SIZE
    # export NSC_TARGET_SMOOTHED_SUFFIX NEURON_TARGET_SMOOTHED_SUFFIX WINDOW_SIZE # Not needed

    # Function to generate comparison plot for a single BigWig - REMOVED
    # export -f generate_comparison_plot - REMOVED

    # --- Run the specific comparison --- 
    log_file="${LOG_DIR_COMPARISON}/4run_compare_Specific_MedipN_vs_MedipPP.log"
    echo "Generating Specific Comparison Plot: ${COMPARISON_PLOT_OUTPUT_FILE} -> Log: ${log_file}"

    # Check if the required input smoothed .npy files exist
    if [[ -f "${NEURON_SMOOTHED_NPY}" && -f "${NSC_SMOOTHED_NPY}" ]]; then
        python -u "${SCRIPT_DIR}/0_1_plot_cpg_metaprofile.py" \
            --plot-comparison \
            --input-npy-file1 "${NEURON_SMOOTHED_NPY}" \
            --label1 "Medip_N (Neuron Targets)" \
            --input-npy-file2 "${NSC_SMOOTHED_NPY}" \
            --label2 "Medip_PP (NSC Targets)" \
            --comparison-output-file "${COMPARISON_PLOT_OUTPUT_FILE}" \
            --window-size $((WINDOW_SIZE / 2)) \
            >"${log_file}" 2>&1
        
        comparison_exit_code=$?
    else
        echo "Error: Skipping specific comparison plot. Missing required smoothed input file(s):" 
        if [ ! -f "${NEURON_SMOOTHED_NPY}" ]; then echo "  - ${NEURON_SMOOTHED_NPY}"; fi
        if [ ! -f "${NSC_SMOOTHED_NPY}" ]; then echo "  - ${NSC_SMOOTHED_NPY}"; fi
        comparison_exit_code=1 # Mark as failed
    fi

    # Wait for comparison jobs - Not needed for single job

    # Report comparison status
    echo "----------------------------------------"
    echo "Summary of Specific Comparison Plotting:"
    if [ ${comparison_exit_code} -eq 0 ]; then
        echo " Specific comparison plot generated successfully: ${COMPARISON_PLOT_OUTPUT_FILE}"
    else
        echo " Error occurred during specific comparison plotting (Exit Code: ${comparison_exit_code})."
        echo " Check log: ${log_file}"
    fi
    echo "----------------------------------------"

fi # End check for individual success

# Final overall status
# Need to check both individual success and the specific comparison success
final_success=false
if $all_individual_success && [ ${comparison_exit_code:-1} -eq 0 ]; then # Default exit code to 1 if not set
    final_success=true
fi

if $final_success; then
     echo "All processing (Individual and Specific Comparison) completed successfully."
     echo "Check outputs in ${OUTPUT_DIR_INDIVIDUAL} and ${OUTPUT_DIR_COMPARISON}."
else
     echo "One or more steps failed during processing. Check logs."
     # Optionally exit with an error code if any job failed
     # exit 1
fi

echo "Targeted CpG metaprofile processing script finished." 