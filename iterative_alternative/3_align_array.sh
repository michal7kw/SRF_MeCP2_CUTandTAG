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

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Source modules
# source /etc/profile.d/modules.sh
# module load bowtie2/2.4.2
# module load samtools/1.13

# Get sample names directly from fastqc results
ALL_SAMPLES=($(ls results/fastqc/*_R1_001_fastqc.html | xargs -n 1 basename | sed 's/_R1_001_fastqc.html//'))

# Verify we have samples before continuing
if [ ${#ALL_SAMPLES[@]} -eq 0 ]; then
    echo "Error: No samples found to process"
    exit 1
fi

# Check if array task ID is valid
if [ $SLURM_ARRAY_TASK_ID -ge ${#ALL_SAMPLES[@]} ]; then
    echo "Error: Array task ID ($SLURM_ARRAY_TASK_ID) exceeds number of samples (${#ALL_SAMPLES[@]})"
    exit 1
fi

# Get current sample
SAMPLE=${ALL_SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Set parameters
GENOME_INDEX="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/mm10_bowtie2_index/mm10"
MAX_FRAGMENT=1000
SORT_MEMORY="2G"
TMP_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/iterative_processing/tmp"

# Create temporary directory for this job
TEMP_DIR=$(mktemp -d -p $TMP_DIR)
if [ ! -d "$TEMP_DIR" ]; then
    echo "Error: Failed to create temporary directory"
    exit 1
fi

# Ensure temporary directory is cleaned up on exit
trap "rm -rf $TEMP_DIR" EXIT

# Create output directories
mkdir -p results/aligned
mkdir -p logs

# Check if input files exist
if [ ! -f "results/trimmed/${SAMPLE}_R1_001_val_1.fq.gz" ] || [ ! -f "results/trimmed/${SAMPLE}_R2_001_val_2.fq.gz" ]; then
    echo "Error: Input files not found for ${SAMPLE}"
    echo "Looking for:"
    echo "  results/trimmed/${SAMPLE}_R1_001_val_1.fq.gz"
    echo "  results/trimmed/${SAMPLE}_R2_001_val_2.fq.gz"
    exit 1
fi

# Check if genome index exists
if [ ! -f "${GENOME_INDEX}.1.bt2" ]; then
    echo "Error: Genome index not found at ${GENOME_INDEX}"
    echo "Expected files like: ${GENOME_INDEX}.1.bt2, ${GENOME_INDEX}.2.bt2, etc."
    exit 1
fi

echo "Starting alignment for ${SAMPLE}..."
echo "Input files:"
echo "  R1: results/trimmed/${SAMPLE}_R1_001_val_1.fq.gz"
echo "  R2: results/trimmed/${SAMPLE}_R2_001_val_2.fq.gz"
echo "Genome index: ${GENOME_INDEX}"
echo "Output SAM: $TEMP_DIR/${SAMPLE}.sam"

# Run Bowtie2 alignment
if ! bowtie2 \
    -p 32 \
    -x $GENOME_INDEX \
    -1 results/trimmed/${SAMPLE}_R1_001_val_1.fq.gz \
    -2 results/trimmed/${SAMPLE}_R2_001_val_2.fq.gz \
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

# Check if SAM file exists and has size
if [ ! -s "$TEMP_DIR/${SAMPLE}.sam" ]; then
    echo "Error: SAM file not created or empty for ${SAMPLE}"
    exit 1
fi

# Convert SAM to BAM and filter
if ! samtools view -@ 32 -bS -h -q 30 \
    "$TEMP_DIR/${SAMPLE}.sam" \
    > results/aligned/${SAMPLE}.unsorted.bam; then
    echo "Error: SAM to BAM conversion failed for ${SAMPLE}"
    exit 1
fi

# Remove intermediate SAM file
rm -f "$TEMP_DIR/${SAMPLE}.sam"

# Check if the unsorted BAM was created successfully
if [ ! -s results/aligned/${SAMPLE}.unsorted.bam ]; then
    echo "Error: Alignment failed for ${SAMPLE}"
    exit 1
fi

# Sort BAM file
samtools sort \
    -@ 32 \
    -m $SORT_MEMORY \
    -T $TEMP_DIR/${SAMPLE} \
    results/aligned/${SAMPLE}.unsorted.bam \
    -o results/aligned/${SAMPLE}.bam

# Check if the sorted BAM was created successfully
if [ ! -s results/aligned/${SAMPLE}.bam ]; then
    echo "Error: Sorting failed for ${SAMPLE}"
    exit 1
fi

# Index BAM file
samtools index -@ 32 results/aligned/${SAMPLE}.bam

# Cleanup
rm -rf $TEMP_DIR
rm results/aligned/${SAMPLE}.unsorted.bam 