#!/bin/bash
#SBATCH --job-name=1_DEGs_lists
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/1_DEGs_lists.err"
#SBATCH --output="logs/1_DEGs_lists.out"

WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/paper_finals"

cd $WORKING_DIR || exit 1

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

python 1_DEGs_lists.py