#!/bin/bash
#SBATCH --job-name=metaprofile_old
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/metaprofile_old.err"
#SBATCH --output="logs/metaprofile_old.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Visualization

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create logs directory if it doesn't exist
mkdir -p logs

# Run both analyses
Rscript metaprofile_analysis_old.R