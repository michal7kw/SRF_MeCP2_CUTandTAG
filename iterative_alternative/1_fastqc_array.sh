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
#SBATCH --error="logs/fastqc_%a.err"
#SBATCH --output="logs/fastqc_%a.out"
#SBATCH --array=0-11

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative

RESULTS_DIR="results_1"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Get sample names
EXOGENOUS_SAMPLES=($(ls ../DATA/EXOGENOUS/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))
ENDOGENOUS_SAMPLES=($(ls ../DATA/ENDOGENOUS/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))
ALL_SAMPLES=("${EXOGENOUS_SAMPLES[@]}" "${ENDOGENOUS_SAMPLES[@]}")

# Get current sample
SAMPLE=${ALL_SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Determine input directory
if [[ -f "../DATA/EXOGENOUS/${SAMPLE}_R1_001.fastq.gz" ]]; then
    INPUT_DIR="../DATA/EXOGENOUS"
else
    INPUT_DIR="../DATA/ENDOGENOUS"
fi

# Create output directory
mkdir -p ${RESULTS_DIR}/fastqc

# Run FastQC
fastqc \
    ${INPUT_DIR}/${SAMPLE}_R1_001.fastq.gz \
    ${INPUT_DIR}/${SAMPLE}_R2_001.fastq.gz \
    --outdir=${RESULTS_DIR}/fastqc \
    --threads=4