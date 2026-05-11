#!/bin/bash
#SBATCH --job-name=seq_depth
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/seq_depth_%a.err"
#SBATCH --output="logs/seq_depth_%a.out"
#SBATCH --array=0-23

# ls -1 ./DATA/{EXOGENOUS,ENDOGENOUS}/*.fastq.gz | wc -l

# Create logs directory if it doesn't exist
mkdir -p logs

# Set base directory and create required directories
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG"
cd $BASE_DIR
# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Define input directories
EXO_DIR="./DATA/EXOGENOUS"
ENDO_DIR="./DATA/ENDOGENOUS"

# Create arrays of input files from both directories
EXO_FILES=(${EXO_DIR}/*.fastq.gz)
ENDO_FILES=(${ENDO_DIR}/*.fastq.gz)

# Combine the arrays
FILES=("${EXO_FILES[@]}" "${ENDO_FILES[@]}")

# Get the current file based on SLURM array task ID
CURRENT_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

echo "Processing file: $CURRENT_FILE"

# Create output directory for FastQC
mkdir -p fastqc_results

# Run both seqkit stats and FastQC
echo "Running seqkit stats..."
seqkit stats "$CURRENT_FILE" > "${CURRENT_FILE%.fastq.gz}_seqkit_stats.txt"

echo "Running FastQC..."
fastqc -o fastqc_results -t 4 "$CURRENT_FILE"

# Also calculate basic stats using zcat
echo "Calculating basic read count..."
READ_COUNT=$(zcat "$CURRENT_FILE" | wc -l | awk '{print $1/4}')
echo "Number of reads in $CURRENT_FILE: $READ_COUNT" > "${CURRENT_FILE%.fastq.gz}_basic_stats.txt"

echo "Done processing $CURRENT_FILE"