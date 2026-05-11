#!/bin/bash
#SBATCH --job-name=2_heatmaps_only_cpgs_R_clean
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/2_heatmaps_only_cpgs_R_clean.err"
#SBATCH --output="logs/2_heatmaps_only_cpgs_R_clean.out"

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
OUTPUT_DIR="${WORKING_DIR}/outputs_r_clean"
CACHE_DIR="${WORKING_DIR}/cache_r"
SCRIPT_FILE="2_heatmaps_only_cpgs_R_clean.R"

# Create necessary directories
mkdir -p logs
mkdir -p ${OUTPUT_DIR}
mkdir -p ${CACHE_DIR}

# Log start of execution
log_info "Starting MeCP2 CpG island analysis job"
log_info "Working directory: $WORKING_DIR"
log_info "BigWig directory: $BIGWIG_DIR"
log_info "CpG islands file: $CPG_ISLANDS_BED"
log_info "Output directory: $OUTPUT_DIR"
log_info "Cache directory: $CACHE_DIR"
log_info "Script file: $SCRIPT_FILE"

# Change to working directory
log_info "Changing to working directory"
if cd $WORKING_DIR; then
    log_success "Successfully changed to working directory"
else
    log_error "Failed to change to working directory"
    exit 1
fi

# Check if script exists
if [ ! -f "$SCRIPT_FILE" ]; then
    log_error "R script file not found: $SCRIPT_FILE"
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
export R_MAX_VSIZE=64Gb
export R_MAX_NUM_THREADS=$SLURM_NTASKS
export OMP_NUM_THREADS=$SLURM_NTASKS
export OMP_THREAD_LIMIT=$SLURM_NTASKS
export BIOCPARALLEL_WORKER_NUMBER=$SLURM_NTASKS

# Record memory usage before execution
log_info "Initial memory usage: $(free -h | grep Mem | awk '{print $3}') used of $(free -h | grep Mem | awk '{print $2}') total"

# Record start time
START_TIME=$(date +%s)
log_info "Starting R script execution at $(date '+%Y-%m-%d %H:%M:%S')"

# Run R script
log_info "Executing: Rscript $SCRIPT_FILE"

if Rscript $SCRIPT_FILE \
    --bigwig-dir="${BIGWIG_DIR}" \
    --cpg-file="${CPG_ISLANDS_BED}" \
    --output-dir="${OUTPUT_DIR}" \
    --cache-dir="${CACHE_DIR}" \
    --cores=$SLURM_NTASKS \
    --force-recompute=FALSE; then
    
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log_success "R script execution completed successfully"
    log_info "Execution time: $((DURATION / 3600)) hours, $(((DURATION % 3600) / 60)) minutes and $((DURATION % 60)) seconds"
    
    # Check output files
    log_info "Checking output files..."
    if [ -d "${OUTPUT_DIR}" ] && [ "$(ls -A ${OUTPUT_DIR})" ]; then
        OUTPUT_COUNT=$(find ${OUTPUT_DIR} -type f | wc -l)
        log_success "Generated $OUTPUT_COUNT output files"
        
        # List PDF and PNG files
        log_info "PDF outputs:"
        find ${OUTPUT_DIR} -name "*.pdf" -type f -exec ls -lh {} \; | awk '{print $9, "("$5")"}' | while read line; do
            log_info "  $line"
        done
        
        log_info "PNG outputs:"
        find ${OUTPUT_DIR} -name "*.png" -type f -exec ls -lh {} \; | awk '{print $9, "("$5")"}' | while read line; do
            log_info "  $line"
        done
    else
        log_warning "No output files were generated in ${OUTPUT_DIR}"
    fi
else
    EXIT_CODE=$?
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log_error "R script execution failed with exit code: $EXIT_CODE"
    log_info "Execution time before failure: $((DURATION / 3600)) hours, $(((DURATION % 3600) / 60)) minutes and $((DURATION % 60)) seconds"
    
    # Check log files for errors
    R_LOG="${WORKING_DIR}/r_script_execution.log"
    if [ -f "$R_LOG" ]; then
        log_info "Last 20 lines of R log file:"
        tail -n 20 "$R_LOG" | while read line; do
            log_info "  $line"
        done
    fi
    
    exit $EXIT_CODE
fi

# Report final memory usage
log_info "Final memory usage: $(free -h | grep Mem | awk '{print $3}') used of $(free -h | grep Mem | awk '{print $2}') total"

# Report R session info
log_info "R session information:"
Rscript -e "writeLines(capture.output(sessionInfo()))"

# Report cache statistics
if [ -d "${CACHE_DIR}" ]; then
    CACHE_SIZE=$(du -sh "${CACHE_DIR}" | cut -f1)
    CACHE_FILES=$(find "${CACHE_DIR}" -type f | wc -l)
    log_info "Cache statistics: ${CACHE_SIZE} used by ${CACHE_FILES} files"
fi

log_info "Job completed at $(date '+%Y-%m-%d %H:%M:%S')"