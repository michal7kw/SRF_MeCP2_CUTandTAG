#!/bin/bash
#SBATCH --job-name=meta_1
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/meta_1.err"
#SBATCH --output="logs/meta_1.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Visualization

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create logs directory if it doesn't exist
mkdir -p logs

Rscript metaprofile_analysis.R
