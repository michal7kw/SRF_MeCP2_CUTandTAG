#!/bin/bash
#SBATCH --job-name=cpg_enrichment_2_rep_in_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment/neu/no_dedup/individual_samples/broad/cpg_enrichment_2_rep_in_peaks.err"
#SBATCH --output="logs/cpg_enrichment/neu/no_dedup/individual_samples/broad/cpg_enrichment_2_rep_in_peaks.out"

CELL_LINE="Neu"
PEAKS="broad"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
SCRIPT_DIR="${BASE_DIR}/scripts/cpg_enrichment"

cd $WORKING_DIR || exit 1

RESULTS_DIR="results/no_dedup/cpg_enrichment/${CELL_LINE}/${PEAKS}/cpg_enrichment_2_rep_in_peaks"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

python ${SCRIPT_DIR}/combine_chunks.py \
    --work-dir "${RESULTS_DIR}"

echo "Chunk combination complete" 