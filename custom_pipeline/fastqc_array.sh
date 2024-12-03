#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/fastqc_%A_%a.err"
#SBATCH --output="logs/fastqc_%A_%a.out"
#SBATCH --array=0-11%4

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline

# Source modules
source /etc/profile.d/modules.sh
module load fastqc/0.11.9

# Get sample names
EXOGENOUS_SAMPLES=($(ls DATA/EXOGENOUS/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))
ENDOGENOUS_SAMPLES=($(ls DATA/ENDOGENOUS/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))
ALL_SAMPLES=("${EXOGENOUS_SAMPLES[@]}" "${ENDOGENOUS_SAMPLES[@]}")

# Get current sample
SAMPLE=${ALL_SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Determine input directory
if [[ -f "DATA/EXOGENOUS/${SAMPLE}_R1_001.fastq.gz" ]]; then
    INPUT_DIR="DATA/EXOGENOUS"
else
    INPUT_DIR="DATA/ENDOGENOUS"
fi

# Create output directory
mkdir -p results/fastqc

# Run FastQC
fastqc \
    ${INPUT_DIR}/${SAMPLE}_R1_001.fastq.gz \
    ${INPUT_DIR}/${SAMPLE}_R2_001.fastq.gz \
    --outdir=results/fastqc \
    --threads=4 