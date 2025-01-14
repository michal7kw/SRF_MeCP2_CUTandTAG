#!/bin/bash
#SBATCH --job-name=combine_replicate_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/combine_replicate_peaks.err"
#SBATCH --output="logs/combine_replicate_peaks.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

mkdir -p logs

PEAKS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results_2_align2_005/peaks/broad"

# Process endogenous peaks
python -u ../scripts/combine_replicate_peaks.py \
    --peaks-dir "$PEAKS_DIR" \
    --condition endo \
    --peak-type broad \
    2>&1 | tee "logs/combine_replicate_peaks_endo.out"

# Process exogenous peaks
python -u ../scripts/combine_replicate_peaks.py \
    --peaks-dir "$PEAKS_DIR" \
    --condition exo \
    --peak-type broad \
    2>&1 | tee "logs/combine_replicate_peaks_exo.out" 