#!/bin/bash
#SBATCH --job-name=align_2
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/align_2.err"
#SBATCH --output="logs/align_2.out"
#SBATCH --array=0-11

# Function for logging with timestamps
log_progress() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

INPUT_DIR="results_1"
RESULTS_DIR="results_1b"

# Function to validate input files
validate_inputs() {
    local sample=$1
    local r1="${INPUT_DIR}/trimmed/${sample}_R1_001_val_1.fq.gz"
    local r2="${INPUT_DIR}/trimmed/${sample}_R2_001_val_2.fq.gz"
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
mkdir -p ${RESULTS_DIR}/aligned logs

# Get sample names from trimmed files (more reliable than fastqc results)
ALL_SAMPLES=($(find ${INPUT_DIR}/trimmed -name "*_R1_001_val_1.fq.gz" -type f | sed 's|${INPUT_DIR}/trimmed/||;s|_R1_001_val_1.fq.gz||' | sort))

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
MAX_FRAGMENT=2000
MIN_FRAGMENT=150
SORT_MEMORY="2G"
THREADS=32
TMP_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/tmpb"
TOTAL_MEMORY="30G"

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
        -p $THREADS \
        -x $GENOME_INDEX \
        -1 ${INPUT_DIR}/trimmed/${SAMPLE}_R1_001_val_1.fq.gz \
        -2 ${INPUT_DIR}/trimmed/${SAMPLE}_R2_001_val_2.fq.gz \
        --local --very-sensitive-local \
        --no-mixed --no-discordant \
        --maxins $MAX_FRAGMENT \
        --minins $MIN_FRAGMENT \
        --dovetail \
        --mm \
        -S "$TEMP_DIR/${SAMPLE}.sam" \
        2> ${RESULTS_DIR}/logs/align_${SAMPLE}.log; then
        echo "Error: Bowtie2 alignment failed for ${SAMPLE}"
        echo "Check ${RESULTS_DIR}/logs/align_${SAMPLE}.log for details"
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

# Modified filtering steps with fragment size selection
log_progress "Applying filters and size selection..."
if ! samtools view -@ "$THREADS" -h "$TEMP_DIR/${SAMPLE}.temp1.bam" | \
    awk -v min=$MIN_FRAGMENT -v max=$MAX_FRAGMENT \
    'BEGIN {OFS="\t"} 
     /^@/ {print; next}
     ($9 >= min && $9 <= max) || ($9 <= -min && $9 >= -max)' | \
    samtools view -@ "$THREADS" -b -f 2 -F 1804 -q 20 > "$TEMP_DIR/${SAMPLE}.temp2.bam"; then
    log_progress "Error: Filtering and size selection failed for ${SAMPLE}"
    exit 1
fi

# Add fragment size distribution analysis
log_progress "Analyzing fragment size distribution..."
samtools view -@ "$THREADS" "$TEMP_DIR/${SAMPLE}.temp2.bam" | \
    awk '{print sqrt($9^2)}' | \
    sort -n | \
    uniq -c > ${RESULTS_DIR}/aligned/${SAMPLE}.fragment_sizes.txt

# Add more detailed QC metrics
log_progress "Generating detailed QC metrics..."
# Calculate percentage of reads in peaks (FRIP)
if [ -f "${INPUT_DIR}/peaks/${SAMPLE}_peaks.narrowPeak" ]; then
    bedtools intersect -a "$TEMP_DIR/${SAMPLE}.temp2.bam" \
        -b "${INPUT_DIR}/peaks/${SAMPLE}_peaks.narrowPeak" \
        -bed -c > ${RESULTS_DIR}/aligned/${SAMPLE}.frip.txt
fi

# Generate insert size metrics
picard CollectInsertSizeMetrics \
    I="$TEMP_DIR/${SAMPLE}.temp2.bam" \
    O=${RESULTS_DIR}/aligned/${SAMPLE}.insert_metrics.txt \
    H=${RESULTS_DIR}/aligned/${SAMPLE}.insert_histogram.pdf

# Remove or comment out the phantompeakqualtools section since it's not installed
# If you need this tool, you should install it first through conda or other means
# phantompeakqualtools run \
#     -c="$TEMP_DIR/${SAMPLE}.temp2.bam" \
#     -p=$THREADS \
#     -out=${RESULTS_DIR}/aligned/${SAMPLE}.cc.qc

# Create comprehensive QC report with error checking
{
    echo "Sample: ${SAMPLE}"
    echo "Date: $(date)"
    
    # Alignment rate
    if [ -f "${RESULTS_DIR}/logs/align_${SAMPLE}.log" ]; then
        echo "Alignment Rate: $(grep "overall alignment rate" ${RESULTS_DIR}/logs/align_${SAMPLE}.log | tail -n1)"
    else
        echo "Alignment Rate: Not available"
    fi
    
    # Fragment Size Statistics
    echo "Fragment Size Statistics:"
    if [ -f "${RESULTS_DIR}/aligned/${SAMPLE}.insert_metrics.txt" ]; then
        awk 'BEGIN{OFS="\t"} NR==1,NR==5{print $1, $2}' "${RESULTS_DIR}/aligned/${SAMPLE}.insert_metrics.txt"
    else
        echo "Not available"
    fi
    
    # Library Complexity
    echo "Library Complexity:"
    if [ -f "${RESULTS_DIR}/aligned/${SAMPLE}.complexity_estimates.txt" ]; then
        head -n 5 "${RESULTS_DIR}/aligned/${SAMPLE}.complexity_estimates.txt"
    else
        echo "Not available"
    fi
} > "${RESULTS_DIR}/aligned/${SAMPLE}.qc_report.txt"

log_progress "Sorting BAM file..."

# Fix the path to use the correct temporary directory structure
if ! samtools sort \
    -@ $THREADS \
    -m $SORT_MEMORY \
    -T "$TMP_DIR/${SAMPLE}" \
    --write-index \
    "$TMP_DIR/${SAMPLE}/${SAMPLE}.temp2.bam" \
    -o ${RESULTS_DIR}/aligned/${SAMPLE}.bam; then
    log_progress "Error: BAM sorting failed for ${SAMPLE}"
    exit 1
fi

log_progress "Indexing BAM file..."
samtools index -@ $THREADS ${RESULTS_DIR}/aligned/${SAMPLE}.bam

log_progress "Removing PCR duplicates..."
# Remove duplicates and save metrics
if ! samtools markdup -@ $THREADS -r \
    ${RESULTS_DIR}/aligned/${SAMPLE}.bam \
    ${RESULTS_DIR}/aligned/${SAMPLE}.dedup.bam \
    2> ${RESULTS_DIR}/aligned/${SAMPLE}.markdup_metrics; then
    log_progress "Error: Duplicate marking failed for ${SAMPLE}"
    exit 1
fi

# Index the deduplicated BAM
samtools index -@ $THREADS ${RESULTS_DIR}/aligned/${SAMPLE}.dedup.bam

log_progress "Generating QC metrics..."
# Generate comprehensive QC metrics
samtools flagstat ${RESULTS_DIR}/aligned/${SAMPLE}.dedup.bam > ${RESULTS_DIR}/aligned/${SAMPLE}.flagstat
samtools idxstats ${RESULTS_DIR}/aligned/${SAMPLE}.dedup.bam > ${RESULTS_DIR}/aligned/${SAMPLE}.idxstats
samtools stats ${RESULTS_DIR}/aligned/${SAMPLE}.dedup.bam > ${RESULTS_DIR}/aligned/${SAMPLE}.stats

# Calculate alignment rate from bowtie2 logs
alignment_rate=$(grep "overall alignment rate" ${RESULTS_DIR}/logs/align_${SAMPLE}.log | tail -n1 | awk '{print $1}')
echo "Overall alignment rate: ${alignment_rate}" >> ${RESULTS_DIR}/aligned/${SAMPLE}.alignment_summary

log_progress "Calculating library complexity metrics..."
preseq lc_extrap -B ${RESULTS_DIR}/aligned/${SAMPLE}.dedup.bam \
    -o ${RESULTS_DIR}/aligned/${SAMPLE}.complexity_estimates.txt

# log_progress "Cleaning up intermediate files..."
# rm -f "$TEMP_DIR/${SAMPLE}.filtered.bam"
# rm -f "$TEMP_DIR/${SAMPLE}.sam"

log_progress "Processing completed successfully for ${SAMPLE}"