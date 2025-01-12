#!/bin/bash
#SBATCH --job-name=peaks_annotation_NSC_combined
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_annotation_NSC_combined.err"
#SBATCH --output="logs/peaks_annotation_NSC_combined.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

mkdir -p logs

GTF_PATH="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA/gencode.vM10.annotation.gtf"
PEAKS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results_2_align2_005/peaks/narrow"
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results_5_align2_005/peaks_annotation_NSC_combined"

python -u ../scripts/peaks_annotation/peaks_annotation_NSC_combined.py \
    --gtf-path "$GTF_PATH" \
    --peaks-dir "$PEAKS_DIR" \
    --output-dir "$OUTPUT_DIR" \
    2>&1 | tee "logs/peaks_annotation_NSC_combined.out" 
