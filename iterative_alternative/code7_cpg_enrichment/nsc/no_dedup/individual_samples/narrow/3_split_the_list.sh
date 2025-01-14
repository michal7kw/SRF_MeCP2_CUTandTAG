#!/bin/bash
#SBATCH --job-name=split_the_list
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment/nsc/no_dedup/individual_samples/narrow/split_the_list.err"
#SBATCH --output="logs/cpg_enrichment/nsc/no_dedup/individual_samples/narrow/split_the_list.out"

CELL_LINE="NSC"
PEAKS="narrow"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
SCRIPT_DIR="${BASE_DIR}/scripts/cpg_enrichment"

cd $WORKING_DIR || exit 1

RESULTS_DIR="results/no_dedup/cpg_enrichment/${CELL_LINE}/${PEAKS}/mecp2_cpg_enrichment_parallel"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

python ${SCRIPT_DIR}/split_the_list.py \
    --work-dir "${RESULTS_DIR}"