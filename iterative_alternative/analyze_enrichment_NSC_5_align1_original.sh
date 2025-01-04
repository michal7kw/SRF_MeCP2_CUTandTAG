#!/bin/bash
#SBATCH --job-name=enrichment_NSC_5_align1_original
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/enrichment_NSC_5_align1_original.err"
#SBATCH --output="logs/enrichment_NSC_5_align1_original.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

mkdir -p logs

WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative"
DATA_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/results_1"
INPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/results_2_align1_original"
RESULTS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/results_5_align1_original"

echo "Copying peaks folder from results_2_align1_original to ${RESULTS_DIR}..."
rm -rf "${RESULTS_DIR}/peaks"
cp -r "${INPUT_DIR}/peaks" "${RESULTS_DIR}/"

python -u ../scripts/analyze_enrichment_NSC_5.py \
    --working-dir $WORKING_DIR \
    --data-dir $DATA_DIR \
    --results-dir $RESULTS_DIR \
    2>&1 | tee "logs/enrichment_NSC_5_align1_original.out"