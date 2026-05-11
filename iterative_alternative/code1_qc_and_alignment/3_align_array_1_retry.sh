#!/bin/bash
#SBATCH --job-name=align_1_retry
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/alignment/align_1_retry.err"
#SBATCH --output="logs/alignment/align_1_retry.out"
#SBATCH --array=0-11

# Function for logging with timestamps
log_progress() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

RESULTS_DIR="results_1"
TMP_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/tmp"
SORT_MEMORY="8G"
THREADS=32

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative || exit 1

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake || exit 1

# Get failed samples (all samples from the original run)
FAILED_SAMPLES=(
    "NeuV2" "NeuM2" "NeuV3" "NeuM3" "NeuV1" 
    "NSCv1" "NSCM1" "NSCM3" "NSCM2" "IgM" 
    "NSCv3" "NSCv2"
)

# Check if array task ID is valid
if [ $SLURM_ARRAY_TASK_ID -ge ${#FAILED_SAMPLES[@]} ]; then
    log_progress "Error: Array task ID ($SLURM_ARRAY_TASK_ID) exceeds number of samples (${#FAILED_SAMPLES[@]})"
    exit 1
fi

# Get current sample
SAMPLE=${FAILED_SAMPLES[$SLURM_ARRAY_TASK_ID]}
TEMP_DIR="$TMP_DIR/${SAMPLE}"

# Create temporary directories
MARKDUP_TMP="$TEMP_DIR/markdup_tmp"
mkdir -p "$MARKDUP_TMP"

# Ensure the filtered BAM exists before proceeding
if [ ! -f "$TEMP_DIR/${SAMPLE}.filtered.bam" ]; then
    log_progress "Error: Filtered BAM file not found for ${SAMPLE}"
    exit 1
fi

log_progress "Processing recovery for ${SAMPLE}..."

log_progress "Sorting BAM file by name for fixmate..."
# Sort by name first (required for fixmate)
if ! samtools sort -n \
    -@ $THREADS \
    -m $SORT_MEMORY \
    -T "$TEMP_DIR/${SAMPLE}_namesort" \
    "$TEMP_DIR/${SAMPLE}.filtered.bam" \
    -o "$TEMP_DIR/${SAMPLE}.namesorted.bam"; then
    log_progress "Error: Name sorting failed for ${SAMPLE}"
    exit 1
fi

log_progress "Running fixmate..."
# Run fixmate to add MS and MC tags
if ! samtools fixmate -m \
    -@ $THREADS \
    "$TEMP_DIR/${SAMPLE}.namesorted.bam" \
    "$TEMP_DIR/${SAMPLE}.fixmate.bam"; then
    log_progress "Error: Fixmate failed for ${SAMPLE}"
    exit 1
fi

log_progress "Sorting BAM file by position..."
# Sort by position (required for markdup)
if ! samtools sort \
    -@ $THREADS \
    -m $SORT_MEMORY \
    -T "$TEMP_DIR/${SAMPLE}_sort" \
    "$TEMP_DIR/${SAMPLE}.fixmate.bam" \
    -o ${RESULTS_DIR}/aligned/${SAMPLE}.bam; then
    log_progress "Error: Position sorting failed for ${SAMPLE}"
    exit 1
fi

log_progress "Indexing BAM file..."
if ! samtools index -@ $THREADS ${RESULTS_DIR}/aligned/${SAMPLE}.bam; then
    log_progress "Error: BAM indexing failed for ${SAMPLE}"
    exit 1
fi

log_progress "Removing PCR duplicates..."
# Add -T flag for temporary directory and --verbosity 3 for detailed output
if ! samtools markdup \
    -@ $THREADS \
    -r \
    -T "$MARKDUP_TMP" \
    --verbosity 3 \
    ${RESULTS_DIR}/aligned/${SAMPLE}.bam \
    ${RESULTS_DIR}/aligned/${SAMPLE}.dedup.bam \
    2> ${RESULTS_DIR}/aligned/${SAMPLE}.markdup_metrics; then
    log_progress "Error: Duplicate marking failed for ${SAMPLE}"
    # Print the last few lines of the markdup metrics for debugging
    tail -n 20 ${RESULTS_DIR}/aligned/${SAMPLE}.markdup_metrics
    exit 1
fi

# Index the deduplicated BAM
log_progress "Indexing deduplicated BAM..."
if ! samtools index -@ $THREADS ${RESULTS_DIR}/aligned/${SAMPLE}.dedup.bam; then
    log_progress "Error: Dedup BAM indexing failed for ${SAMPLE}"
    exit 1
fi

log_progress "Generating QC metrics..."
samtools flagstat ${RESULTS_DIR}/aligned/${SAMPLE}.dedup.bam > ${RESULTS_DIR}/aligned/${SAMPLE}.flagstat
samtools idxstats ${RESULTS_DIR}/aligned/${SAMPLE}.dedup.bam > ${RESULTS_DIR}/aligned/${SAMPLE}.idxstats
samtools stats ${RESULTS_DIR}/aligned/${SAMPLE}.dedup.bam > ${RESULTS_DIR}/aligned/${SAMPLE}.stats

# Calculate alignment rate from bowtie2 logs if they exist
if [ -f "${RESULTS_DIR}/logs/align_${SAMPLE}.log" ]; then
    alignment_rate=$(grep "overall alignment rate" ${RESULTS_DIR}/logs/align_${SAMPLE}.log | tail -n1 | awk '{print $1}')
    echo "Overall alignment rate: ${alignment_rate}" >> ${RESULTS_DIR}/aligned/${SAMPLE}.alignment_summary
fi

log_progress "Calculating library complexity metrics..."
if ! preseq lc_extrap -B ${RESULTS_DIR}/aligned/${SAMPLE}.dedup.bam \
    -o ${RESULTS_DIR}/aligned/${SAMPLE}.complexity_estimates.txt; then
    log_progress "Warning: Library complexity estimation failed for ${SAMPLE}"
fi

# Clean up temporary files
rm -f "$TEMP_DIR/${SAMPLE}.namesorted.bam"
rm -f "$TEMP_DIR/${SAMPLE}.fixmate.bam"
rm -rf "$MARKDUP_TMP"

log_progress "Processing completed successfully for ${SAMPLE}" 