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
#SBATCH --error="logs/align_%A_%a.err"
#SBATCH --output="logs/align_%A_%a.out"
#SBATCH --array=0-11%4

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline

# Source modules
source /etc/profile.d/modules.sh
module load bowtie2/2.4.2
module load samtools/1.13

# Get sample names
EXOGENOUS_SAMPLES=($(ls DATA/EXOGENOUS/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))
ENDOGENOUS_SAMPLES=($(ls DATA/ENDOGENOUS/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))
ALL_SAMPLES=("${EXOGENOUS_SAMPLES[@]}" "${ENDOGENOUS_SAMPLES[@]}")

# Get current sample
SAMPLE=${ALL_SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Set parameters
GENOME_INDEX="/path/to/bowtie2/index/hg38"
MAX_FRAGMENT=1000
SORT_MEMORY="2G"
TMP_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/tmp"

# Create output directories
mkdir -p results/aligned
mkdir -p $TMP_DIR

# Create temporary directory for this job
TEMP_DIR=$(mktemp -d -p $TMP_DIR)

# Run Bowtie2 alignment
bowtie2 \
    -p 32 \
    -x $GENOME_INDEX \
    -1 results/trimmed/${SAMPLE}_R1_001_val_1.fq.gz \
    -2 results/trimmed/${SAMPLE}_R2_001_val_2.fq.gz \
    --local --very-sensitive-local \
    --no-mixed --no-discordant \
    --maxins $MAX_FRAGMENT \
    --mm \
    2> logs/align_${SAMPLE}.log | \
samtools view -@ 32 -bS -q 30 - > results/aligned/${SAMPLE}.unsorted.bam

# Sort BAM file
samtools sort \
    -@ 32 \
    -m $SORT_MEMORY \
    -T $TEMP_DIR/${SAMPLE} \
    -o results/aligned/${SAMPLE}.bam \
    results/aligned/${SAMPLE}.unsorted.bam

# Index BAM file
samtools index -@ 32 results/aligned/${SAMPLE}.bam

# Cleanup
rm -rf $TEMP_DIR
rm results/aligned/${SAMPLE}.unsorted.bam 