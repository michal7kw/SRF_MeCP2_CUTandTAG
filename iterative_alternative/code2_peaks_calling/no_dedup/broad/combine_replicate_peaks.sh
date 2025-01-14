#!/bin/bash
#SBATCH --job-name=crp_no_dedup_broad
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_calling/no_dedup/combine_replicate_peaks_broad.err"
#SBATCH --output="logs/peaks_calling/no_dedup/combine_replicate_peaks_broad.out"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"

cd $WORKING_DIR || exit 1

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

mkdir -p logs

PEAKS_DIR="results/no_dedup/peaks/broad"

# Process exogenous peaks
python -u ../scripts/combine_peaks/combine_replicate_peaks_broad.py \
    --peaks-dir "$PEAKS_DIR" \
    --condition exo \
    2>&1 | tee "logs/peaks_calling/no_dedup/combine_replicate_peaks_broad_exo.out" 

# Process exogenous peaks
python -u ../scripts/combine_peaks/combine_replicate_peaks_broad.py \
    --peaks-dir "$PEAKS_DIR" \
    --condition endo \
    2>&1 | tee "logs/peaks_calling/no_dedup/combine_replicate_peaks_broad_endo.out" 
