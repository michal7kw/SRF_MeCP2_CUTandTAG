#!/bin/bash
#SBATCH --job-name=2R_heatmaps_peaks_and_genes
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/2R_heatmaps_peaks_and_genes.err"
#SBATCH --output="logs/2R_heatmaps_peaks_and_genes.out"

WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/paper_finals"

cd $WORKING_DIR || exit 1

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake


Rscript 2_heatmaps_peaks_and_genes_with_profiles.R

