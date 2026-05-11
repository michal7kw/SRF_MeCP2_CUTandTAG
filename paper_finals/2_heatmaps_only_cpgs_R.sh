#!/bin/bash
#SBATCH --job-name=2_heatmaps_only_cpgs_R
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/2_heatmaps_only_cpgs_R.err"
#SBATCH --output="logs/2_heatmaps_only_cpgs_R.out"

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

# Log system information
log_info "=== SYSTEM INFORMATION ==="
log_info "Hostname: $(hostname)"
log_info "CPU: $(lscpu | grep 'Model name' | sed 's/Model name:[[:space:]]*//')"
log_info "Memory: $(free -h | grep Mem | awk '{print $2}') total"
log_info "Slurm Job ID: $SLURM_JOB_ID"
log_info "CPU Threads: $SLURM_NTASKS"

# Define working directory and data paths
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/paper_finals"
BIGWIG_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results_1b/bigwig"
CPG_ISLANDS_BED="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA/cpg_islands.bed"
OUTPUT_DIR="${WORKING_DIR}/outputs_r"

# Create output and log directories if they don't exist
mkdir -p logs
mkdir -p ${OUTPUT_DIR}

# Log start of execution
log_info "Starting R implementation job execution"
log_info "Working directory: $WORKING_DIR"
log_info "BigWig directory: $BIGWIG_DIR"
log_info "CpG islands file: $CPG_ISLANDS_BED"
log_info "Output directory: $OUTPUT_DIR"

# Change to working directory
log_info "Changing to working directory"
if cd $WORKING_DIR; then
    log_success "Successfully changed to working directory"
else
    log_error "Failed to change to working directory"
    exit 1
fi

# Activate conda environment
log_info "Activating conda environment"
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
if conda activate snakemake; then
    log_success "Successfully activated snakemake conda environment"
else
    log_error "Failed to activate snakemake conda environment"
    exit 1
fi

# Set R-specific environment variables
export R_MAX_VSIZE=64Gb  # Maximum memory size for R vector heap
export R_MAX_NUM_THREADS=$SLURM_NTASKS  # Set max number of threads

# Record start time
START_TIME=$(date +%s)
log_info "Starting R script execution at $(date '+%Y-%m-%d %H:%M:%S')"

# Run R script with timing and parameters
log_info "Executing: Rscript 2_heatmaps_only_cpgs_R.R"

if Rscript 2_heatmaps_only_cpgs_R.R \
    --bw-dir=${BIGWIG_DIR} \
    --cpg-file=${CPG_ISLANDS_BED} \
    --output-dir=${OUTPUT_DIR} \
    --cores=$SLURM_NTASKS; then
    
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log_success "R script execution completed successfully"
    log_info "Execution time: $((DURATION / 3600)) hours, $(((DURATION % 3600) / 60)) minutes and $((DURATION % 60)) seconds"
    
    # Check if output files were created
    log_info "Checking output files..."
    if [ -d "${OUTPUT_DIR}" ] && [ "$(ls -A ${OUTPUT_DIR})" ]; then
        log_success "Output files generated successfully"
        find ${OUTPUT_DIR} -name "*.png" -type f -exec ls -lh {} \; | awk '{print $9, "("$5")"}' | while read line; do
            log_info "Generated: $line"
        done
    else
        log_warning "Output files may not have been generated correctly"
    fi
else
    EXIT_CODE=$?
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log_error "R script execution failed with exit code: $EXIT_CODE"
    log_info "Execution time before failure: $((DURATION / 3600)) hours, $(((DURATION % 3600) / 60)) minutes and $((DURATION % 60)) seconds"
    exit $EXIT_CODE
fi

# Report memory usage at the end
log_info "Final memory usage: $(free -h | grep Mem | awk '{print $3}') used of $(free -h | grep Mem | awk '{print $2}') total"

# Report R session info
log_info "R session information:"
Rscript -e "writeLines(capture.output(sessionInfo()))"

log_info "Job completed at $(date '+%Y-%m-%d %H:%M:%S')"