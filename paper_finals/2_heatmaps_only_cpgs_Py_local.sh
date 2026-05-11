#!/bin/bash

# Set up logging functions
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $1"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $1" >&2
}

log_warning() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARNING] $1" >&2
}

log_success() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [SUCCESS] $1"
}

# Log system information (adjust if commands are not available locally)
log_info "=== SYSTEM INFORMATION ==="
log_info "Hostname: $(hostname)"
# log_info "CPU: $(lscpu | grep 'Model name' | sed 's/Model name:[[:space:]]*//')" # May need adjustment or removal depending on OS
log_info "Memory: $(free -h | grep Mem | awk '{print $2}') total" # May need adjustment or removal depending on OS / command availability
NUM_THREADS=16 # Default number of threads for local execution, adjust as needed
log_info "CPU Threads: $NUM_THREADS (Set manually in script)"

# --- !!! ADJUST THESE PATHS FOR YOUR LOCAL SETUP !!! ---
# Define working directory and data paths
WORKING_DIR="." # Assumes you run this script from the 'paper_finals' directory
BIGWIG_DIR="../data/bigwig_mecp2" # Example: Relative path to your bigwig files
CPG_ISLANDS_BED="../data/cpg_islands.bed" # Example: Relative path to your CpG islands bed file
OUTPUT_DIR="${WORKING_DIR}/outputs_python_local"
# CONDA_ACTIVATION_CMD="source /path/to/your/miniconda3/bin/activate snakemake" # Example: Update with your conda path and environment name
# --- !!! END OF PATH ADJUSTMENTS !!! ---

# Create output and log directories if they don't exist
mkdir -p logs_local # Use a different log directory for local runs
mkdir -p ${OUTPUT_DIR}

# Redirect logs to local files
LOG_FILE="logs_local/2_heatmaps_only_cpgs_Py_local.log"
exec > >(tee -a "${LOG_FILE}") 2> >(tee -a "${LOG_FILE}" >&2)

# Log start of execution
log_info "Starting Python local job execution"
log_info "Working directory: $WORKING_DIR"
log_info "BigWig directory: $BIGWIG_DIR"
log_info "CpG islands file: $CPG_ISLANDS_BED"
log_info "Output directory: $OUTPUT_DIR"

# Change to working directory (if needed, script assumes it runs from paper_finals)
# log_info "Changing to working directory"
# if cd $WORKING_DIR; then
#     log_success "Successfully changed to working directory"
# else
#     log_error "Failed to change to working directory"
#     exit 1
# fi

# Activate conda environment
# log_info "Activating conda environment"
# eval $CONDA_ACTIVATION_CMD # Use eval to handle potential complex activation commands
# if [ $? -eq 0 ]; then
#     log_success "Successfully activated specified conda environment"
# else
#     log_error "Failed to activate specified conda environment. Check CONDA_ACTIVATION_CMD."
#     exit 1
# fi

# Record start time
START_TIME=$(date +%s)
log_info "Starting Python script execution at $(date '+%Y-%m-%d %H:%M:%S')"

# Run Python script with timing
log_info "Executing: python 2_heatmaps_only_cpgs_Py.py --bigwig-dir ${BIGWIG_DIR} --cpg-file ${CPG_ISLANDS_BED} --output-dir ${OUTPUT_DIR} --threads $NUM_THREADS"

if python 2_heatmaps_only_cpgs_Py.py \
    --bigwig-dir ${BIGWIG_DIR} \
    --cpg-file ${CPG_ISLANDS_BED} \
    --output-dir ${OUTPUT_DIR} \
    --threads $NUM_THREADS; then
    
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log_success "Python script execution completed successfully"
    log_info "Execution time: $((DURATION / 3600)) hours, $(((DURATION % 3600) / 60)) minutes and $((DURATION % 60)) seconds"
    
    # Check if output files were created
    log_info "Checking output files..."
    # Check for expected output files (adjust if names change)
    if [ -f "${OUTPUT_DIR}/all_samples_cpg_profile.png" ] && \
       [ -f "${OUTPUT_DIR}/neuron_endo_exo_side_by_side.png" ] && \
       [ -f "${OUTPUT_DIR}/nsc_endo_exo_side_by_side.png" ]; then
        log_success "Expected output files generated successfully"
        ls -lh ${OUTPUT_DIR}/*.png | awk '{print $9, "("$5")"}' | while read line; do
            log_info "Generated: $line"
        done
    else
        log_warning "One or more expected output files may not have been generated correctly"
        log_info "Listing contents of ${OUTPUT_DIR}:"
        ls -lh ${OUTPUT_DIR}
    fi
else
    EXIT_CODE=$?
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log_error "Python script execution failed with exit code: $EXIT_CODE"
    log_info "Execution time before failure: $((DURATION / 3600)) hours, $(((DURATION % 3600) / 60)) minutes and $((DURATION % 60)) seconds"
    exit $EXIT_CODE
fi

# Report memory usage at the end (adjust if command not available)
# log_info "Final memory usage: $(free -h | grep Mem | awk '{print $3}') used of $(free -h | grep Mem | awk '{print $2}') total"

log_info "Local job finished at $(date '+%Y-%m-%d %H:%M:%S')" 