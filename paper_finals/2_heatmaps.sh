#!/bin/bash
#SBATCH --job-name=2_heatmaps
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/2_heatmaps.err"
#SBATCH --output="logs/2_heatmaps.out"

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

# Define working directory
WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/paper_finals"

# Log start of execution
log_info "Starting job execution"
log_info "Working directory: $WORKING_DIR"

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

# Record start time
START_TIME=$(date +%s)
log_info "Starting R script execution at $(date '+%Y-%m-%d %H:%M:%S')"

# Run R script with timing
log_info "Executing: Rscript 2_heatmaps.R"
if Rscript 2_heatmaps.R; then
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log_success "R script execution completed successfully"
    log_info "Execution time: $((DURATION / 60)) minutes and $((DURATION % 60)) seconds"
else
    EXIT_CODE=$?
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log_error "R script execution failed with exit code: $EXIT_CODE"
    log_info "Execution time before failure: $((DURATION / 60)) minutes and $((DURATION % 60)) seconds"
    exit $EXIT_CODE
fi

log_info "Job completed at $(date '+%Y-%m-%d %H:%M:%S')"