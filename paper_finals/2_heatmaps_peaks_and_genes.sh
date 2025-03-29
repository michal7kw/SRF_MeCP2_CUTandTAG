#!/bin/bash
#SBATCH --job-name=2_heatmaps_peaks_and_genes
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/2_heatmaps_peaks_and_genes.err"
#SBATCH --output="logs/2_heatmaps_peaks_and_genes.out"

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
PEAKS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/peaks/broad"
OUTPUT_DIR="${WORKING_DIR}/outputs_peaks_and_genes"

# Create temporary directory for gene regions and TSS files
TEMP_DIR="${OUTPUT_DIR}/temp"
mkdir -p ${TEMP_DIR}

# Define paths for gene regions and TSS files
GENE_REGIONS_BED="${TEMP_DIR}/gene_regions.bed"
TSS_BED="${TEMP_DIR}/tss_regions.bed"

# Create output and log directories if they don't exist
mkdir -p logs
mkdir -p ${OUTPUT_DIR}

# Log start of execution
log_info "Starting Python implementation job execution"
log_info "Working directory: $WORKING_DIR"
log_info "BigWig directory: $BIGWIG_DIR"
log_info "Peaks directory: $PEAKS_DIR"
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

# Generate gene regions and TSS files using R script
log_info "Generating gene regions and TSS files"
cat > ${TEMP_DIR}/generate_regions.R << 'EOF'
#!/usr/bin/env Rscript

library(GenomicRanges)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
gene_regions_file <- args[1]
tss_file <- args[2]

# Load TxDb for mouse genome
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Extract gene regions
genes <- genes(txdb)
cat(paste("Extracted", length(genes), "gene regions\n"))

# Export gene regions
export(genes, gene_regions_file)
cat(paste("Gene regions exported to:", gene_regions_file, "\n"))

# Extract TSS regions
tss_regions <- promoters(genes, upstream=1, downstream=1)
cat(paste("Created", length(tss_regions), "TSS regions\n"))

# Export TSS regions
export(tss_regions, tss_file)
cat(paste("TSS regions exported to:", tss_file, "\n"))
EOF

# Make the R script executable
chmod +x ${TEMP_DIR}/generate_regions.R

# Run the R script to generate regions
log_info "Running R script to generate gene regions and TSS files"
if Rscript ${TEMP_DIR}/generate_regions.R ${GENE_REGIONS_BED} ${TSS_BED}; then
    log_success "Successfully generated gene regions and TSS files"
else
    log_error "Failed to generate gene regions and TSS files"
    exit 1
fi

# Record start time
START_TIME=$(date +%s)
log_info "Starting Python script execution at $(date '+%Y-%m-%d %H:%M:%S')"

# Run Python script with timing
log_info "Executing: python 2_heatmaps_peaks_and_genes.py"

if python 2_heatmaps_peaks_and_genes.py \
    --bigwig-dir ${BIGWIG_DIR} \
    --peaks-dir ${PEAKS_DIR} \
    --gene-regions-bed ${GENE_REGIONS_BED} \
    --tss-bed ${TSS_BED} \
    --output-dir ${OUTPUT_DIR} \
    --upstream 5000 \
    --downstream 5000 \
    --bin-size 50 \
    --body-bins 100 \
    --threads $SLURM_NTASKS; then
    
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log_success "Python script execution completed successfully"
    log_info "Execution time: $((DURATION / 3600)) hours, $(((DURATION % 3600) / 60)) minutes and $((DURATION % 60)) seconds"
    
    # Check debug files
    log_info "Checking debug files..."
    if [ -f "${OUTPUT_DIR}/debug_raw_signals.png" ] && [ -f "${OUTPUT_DIR}/debug_tss_tes_signals.csv" ]; then
        log_info "Debug data generated successfully:"
        log_info "Debug Plot: ${OUTPUT_DIR}/debug_raw_signals.png"
        log_info "Debug Data: ${OUTPUT_DIR}/debug_tss_tes_signals.csv"
        
        # Print TSS/TES ratio information
        if [ -f "${OUTPUT_DIR}/debug_tss_tes_signals.csv" ]; then
            log_info "TSS/TES Signal Ratios:"
            tail -n +2 "${OUTPUT_DIR}/debug_tss_tes_signals.csv" | while IFS=, read -r sample tss_signal tes_signal ratio; do
                log_info "$sample: TSS=$tss_signal, TES=$tes_signal, Ratio=$ratio"
            done
        fi
    fi
    
    # Check if output files were created
    log_info "Checking output files..."
    if [ -f "${OUTPUT_DIR}/all_samples_tss_profile.png" ] && [ -f "${OUTPUT_DIR}/all_samples_gene_body_profile.png" ]; then
        log_success "Output files generated successfully"
        ls -lh ${OUTPUT_DIR}/*.png | awk '{print $9, "("$5")"}' | while read line; do
            log_info "Generated: $line"
        done
    else
        log_warning "Output files may not have been generated correctly"
    fi
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