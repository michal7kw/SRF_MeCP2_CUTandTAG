#!/bin/bash
#SBATCH --job-name=2_heatmaps_only_cpgs_Py
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/2_heatmaps_only_cpgs_Py.err"
#SBATCH --output="logs/2_heatmaps_only_cpgs_Py.out"

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
OUTPUT_DIR="${WORKING_DIR}/outputs_python"

# Create output and log directories if they don't exist
mkdir -p logs
mkdir -p ${OUTPUT_DIR}

# Log start of execution
log_info "Starting Python implementation job execution"
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

# Record start time
START_TIME=$(date +%s)
log_info "Starting Python script execution at $(date '+%Y-%m-%d %H:%M:%S')"

# Run Python script with timing
log_info "Executing: python 2_heatmaps_only_cpgs_Py.py --bigwig-dir ${BIGWIG_DIR} --cpg-file ${CPG_ISLANDS_BED} --output-dir ${OUTPUT_DIR} --threads $SLURM_NTASKS"

if python 2_heatmaps_only_cpgs_Py.py \
    --bigwig-dir ${BIGWIG_DIR} \
    --cpg-file ${CPG_ISLANDS_BED} \
    --output-dir ${OUTPUT_DIR} \
    --threads $SLURM_NTASKS; then
    
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log_success "Python script execution completed successfully"
    log_info "Execution time: $((DURATION / 3600)) hours, $(((DURATION % 3600) / 60)) minutes and $((DURATION % 60)) seconds"
    
    # Check if output files were created
    log_info "Checking output files..."
    expected_files=(
        "${OUTPUT_DIR}/neuron_endo_exo_cpg_profile.png"
        "${OUTPUT_DIR}/neuron_endo_exo_side_by_side.png"
        "${OUTPUT_DIR}/nsc_endo_exo_cpg_profile.png"
        "${OUTPUT_DIR}/nsc_endo_exo_side_by_side.png"
    )
    
    all_exist=true
    for file in "${expected_files[@]}"; do
        if [ ! -f "$file" ]; then
            log_warning "Missing expected output file: $file"
            all_exist=false
        fi
    done
    
    if $all_exist; then
        log_success "All expected output files generated successfully"
    else
        log_warning "Some expected output files may be missing"
    fi
    
    # List all generated PNG files
    log_info "Generated PNG files:"
    ls -lh ${OUTPUT_DIR}/*.png | awk '{print $9, "("$5")"}' | while read line; do
        log_info "$line"
    done
else
    EXIT_CODE=$?
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log_error "Python script execution failed with exit code: $EXIT_CODE"
    log_info "Execution time before failure: $((DURATION / 3600)) hours, $(((DURATION % 3600) / 60)) minutes and $((DURATION % 60)) seconds"
    exit $EXIT_CODE
fi

# Report memory usage at the end
log_info "Final memory usage: $(free -h | grep Mem | awk '{print $3}') used of $(free -h | grep Mem | awk '{print $2}') total"

log_info "Job completed at $(date '+%Y-%m-%d %H:%M:%S')"