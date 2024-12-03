#!/bin/bash
#SBATCH --job-name=peaks_seacr
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_%A_%a.err"
#SBATCH --output="logs/peaks_%A_%a.out"
#SBATCH --array=0-9%4  # Excluding IgM controls

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline

# Source modules
source /etc/profile.d/modules.sh
module load bedtools/2.29.1
module load R/4.1.0

# Path to SEACR script
SEACR="/path/to/SEACR/SEACR_1.3.sh"

# Get sample names (excluding IgM)
EXOGENOUS_SAMPLES=($(ls DATA/EXOGENOUS/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//' | grep -v "IgM"))
ENDOGENOUS_SAMPLES=($(ls DATA/ENDOGENOUS/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//' | grep -v "IgM"))
ALL_SAMPLES=("${EXOGENOUS_SAMPLES[@]}" "${ENDOGENOUS_SAMPLES[@]}")

# Get current sample
SAMPLE=${ALL_SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Create output directories
mkdir -p results/peaks/seacr
mkdir -p results/bedgraph

# Convert BAM to bedgraph
bedtools genomecov -bg -ibam results/aligned/${SAMPLE}.bam \
    > results/bedgraph/${SAMPLE}.bedgraph

# Get IgM control from appropriate directory
if [[ -f "DATA/EXOGENOUS/${SAMPLE}_R1_001.fastq.gz" ]]; then
    CONTROL_BG="results/bedgraph/IgM_exo.bedgraph"
else
    CONTROL_BG="results/bedgraph/IgM_endo.bedgraph"
fi

# Run SEACR peak calling
bash $SEACR \
    results/bedgraph/${SAMPLE}.bedgraph \
    $CONTROL_BG \
    non \
    stringent \
    results/peaks/seacr/${SAMPLE}

# Also run without control (top 1% of peaks)
bash $SEACR \
    results/bedgraph/${SAMPLE}.bedgraph \
    0.01 \
    non \
    stringent \
    results/peaks/seacr/${SAMPLE}_no_control 