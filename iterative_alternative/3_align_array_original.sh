#!/bin/bash
#SBATCH --job-name=align
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/align_%a.err"
#SBATCH --output="logs/align_%a.out"
#SBATCH --array=0-11

# Function for logging with timestamps
log_progress() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to validate input files
validate_inputs() {
    local sample=$1
    local r1="results_1/trimmed/${sample}_R1_001_val_1.fq.gz"
    local r2="results_1/trimmed/${sample}_R2_001_val_2.fq.gz"
    [[ -f "$r1" && -f "$r2" ]] || return 1
}

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative || {
    log_progress "Error: Failed to change to working directory"
    exit 1
}

# Activate conda environment
if ! source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake; then
    log_progress "Error: Failed to activate conda environment"
    exit 1
fi

# Create required directories
mkdir -p results_1/aligned logs

# Get sample names from trimmed files (more reliable than fastqc results)
ALL_SAMPLES=($(find results_1/trimmed -name "*_R1_001_val_1.fq.gz" -type f | sed 's|results_1/trimmed/||;s|_R1_001_val_1.fq.gz||' | sort))

# Verify we have samples before continuing
if [ ${#ALL_SAMPLES[@]} -eq 0 ]; then
    log_progress "Error: No samples found to process"
    exit 1
fi

# Check if array task ID is valid
if [ $SLURM_ARRAY_TASK_ID -ge ${#ALL_SAMPLES[@]} ]; then
    log_progress "Error: Array task ID ($SLURM_ARRAY_TASK_ID) exceeds number of samples (${#ALL_SAMPLES[@]})"
    exit 1
fi

# Get current sample
SAMPLE=${ALL_SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Set parameters
GENOME_INDEX="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/mm10_bowtie2_index/mm10"
MAX_FRAGMENT=1000
SORT_MEMORY="32G"
THREADS=32
TMP_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/iterative_processing/tmp"

# Create temporary directory if it doesn't exist
if [ ! -d "$TMP_DIR" ]; then
    log_progress "Creating temporary directory..."
    mkdir -p "$TMP_DIR" || {
        log_progress "Error: Failed to create temporary directory at $TMP_DIR"
        exit 1
    }
fi

# Create fixed temporary directory for this sample
TEMP_DIR="$TMP_DIR/${SAMPLE}" 
mkdir -p "$TEMP_DIR" || {
    log_progress "Error: Failed to create job-specific temporary directory"
    exit 1
}

# Ensure temporary directory is cleaned up on exit
# trap 'log_progress "Cleaning up temporary files..."; rm -rf "$TEMP_DIR"' EXIT

# Validate input files
if ! validate_inputs "$SAMPLE"; then
    log_progress "Error: Input files not found for ${SAMPLE}"
    exit 1
fi

# Check if genome index exists
if [ ! -f "${GENOME_INDEX}.1.bt2" ]; then
    log_progress "Error: Genome index not found at ${GENOME_INDEX}"
    exit 1
fi

log_progress "Checking for existing SAM file..."
if [ -s "$TEMP_DIR/${SAMPLE}.sam" ]; then
    log_progress "Found existing SAM file, skipping Bowtie2 alignment..."
else
    log_progress "Starting alignment for ${SAMPLE}..."
    # Run Bowtie2 alignment
    if ! bowtie2 \
        -p 32 \
        -x $GENOME_INDEX \
        -1 results_1/trimmed/${SAMPLE}_R1_001_val_1.fq.gz \
        -2 results_1/trimmed/${SAMPLE}_R2_001_val_2.fq.gz \
        --local --very-sensitive-local \
        --no-mixed --no-discordant \
        --maxins $MAX_FRAGMENT \
        --mm \
        -S "$TEMP_DIR/${SAMPLE}.sam" \
        2> logs/align_${SAMPLE}.log; then
        echo "Error: Bowtie2 alignment failed for ${SAMPLE}"
        echo "Check logs/align_${SAMPLE}.log for details"
        exit 1
    fi
fi

# Check if SAM file exists and has size
if [ ! -s "$TEMP_DIR/${SAMPLE}.sam" ]; then
    echo "Error: SAM file not created or empty for ${SAMPLE}"
    exit 1
fi

log_progress "Converting SAM to BAM and filtering..."

# First convert SAM to BAM with basic filtering
if ! samtools view -@ "$THREADS" -b -h -q 20 "$TEMP_DIR/${SAMPLE}.sam" > "$TEMP_DIR/${SAMPLE}.temp1.bam"; then
    log_progress "Error: Initial SAM to BAM conversion failed for ${SAMPLE}"
    exit 1
fi

# Then apply additional filters in separate steps
log_progress "Applying additional filters..."
if ! samtools view -@ "$THREADS" -b -f 2 -F 1804 "$TEMP_DIR/${SAMPLE}.temp1.bam" > "$TEMP_DIR/${SAMPLE}.temp2.bam"; then
    log_progress "Error: Filtering step 1 failed for ${SAMPLE}"
    exit 1
fi

# Final filtering step for mapping quality and read length
if ! samtools view -@ "$THREADS" -h -b -q 20 "$TEMP_DIR/${SAMPLE}.temp2.bam" > "$TEMP_DIR/${SAMPLE}.filtered.bam"; then
    log_progress "Error: Final filtering step failed for ${SAMPLE}"
    exit 1
fi

# Clean up temporary files
rm -f "$TEMP_DIR/${SAMPLE}.temp1.bam" "$TEMP_DIR/${SAMPLE}.temp2.bam"

log_progress "Sorting BAM file..."

# Sort BAM file
if ! samtools sort \
    -@ $THREADS \
    -m $SORT_MEMORY \
    -T "$TEMP_DIR/${SAMPLE}" \
    "$TEMP_DIR/${SAMPLE}.filtered.bam" \
    -o results_1/aligned/${SAMPLE}.bam; then
    log_progress "Error: BAM sorting failed for ${SAMPLE}"
    exit 1
fi

log_progress "Indexing BAM file..."
samtools index -@ $THREADS results_1/aligned/${SAMPLE}.bam

log_progress "Removing PCR duplicates..."
# Remove duplicates and save metrics
if ! samtools markdup -@ $THREADS -r \
    results_1/aligned/${SAMPLE}.bam \
    results_1/aligned/${SAMPLE}.dedup.bam \
    2> results_1/aligned/${SAMPLE}.markdup_metrics; then
    log_progress "Error: Duplicate marking failed for ${SAMPLE}"
    exit 1
fi

# Index the deduplicated BAM
samtools index -@ $THREADS results_1/aligned/${SAMPLE}.dedup.bam

log_progress "Generating QC metrics..."
# Generate comprehensive QC metrics
samtools flagstat results_1/aligned/${SAMPLE}.dedup.bam > results_1/aligned/${SAMPLE}.flagstat
samtools idxstats results_1/aligned/${SAMPLE}.dedup.bam > results_1/aligned/${SAMPLE}.idxstats
samtools stats results_1/aligned/${SAMPLE}.dedup.bam > results_1/aligned/${SAMPLE}.stats

# Calculate alignment rate from bowtie2 logs
alignment_rate=$(grep "overall alignment rate" logs/align_${SAMPLE}.log | tail -n1 | awk '{print $1}')
echo "Overall alignment rate: ${alignment_rate}" >> results_1/aligned/${SAMPLE}.alignment_summary

log_progress "Calculating library complexity metrics..."
preseq lc_extrap -B results_1/aligned/${SAMPLE}.dedup.bam \
    -o results_1/aligned/${SAMPLE}.complexity_estimates.txt

# log_progress "Cleaning up intermediate files..."
# rm -f "$TEMP_DIR/${SAMPLE}.filtered.bam"
# rm -f "$TEMP_DIR/${SAMPLE}.sam"

log_progress "Processing completed successfully for ${SAMPLE}"