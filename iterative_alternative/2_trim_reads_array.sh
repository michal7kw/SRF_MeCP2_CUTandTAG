#!/bin/bash
#SBATCH --job-name=trim_reads
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/trim_%a.err"
#SBATCH --output="logs/trim_%a.out"
#SBATCH --array=0-11

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

RESULTS_DIR="results_1"

# Get sample names (same as in fastqc script)
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
mkdir -p ${RESULTS_DIR}/trimmed

# Run Trim Galore
trim_galore \
    --paired \
    --gzip \
    --fastqc \
    --cores 16 \
    --output_dir ${RESULTS_DIR}/trimmed \
    ${INPUT_DIR}/${SAMPLE}_R1_001.fastq.gz \
    ${INPUT_DIR}/${SAMPLE}_R2_001.fastq.gz 