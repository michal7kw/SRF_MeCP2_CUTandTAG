#!/bin/bash
#SBATCH --job-name=peak_annotation
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peak_annotation.err"
#SBATCH --output="logs/peak_annotation.out"

# Set up logging
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $1"
}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Run the R script
log_info "Starting peak annotation analysis..."
Rscript peak_annotation.R

log_info "Analysis completed"
