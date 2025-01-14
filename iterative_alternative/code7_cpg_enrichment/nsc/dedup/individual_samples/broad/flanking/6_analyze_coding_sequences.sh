#!/bin/bash
#SBATCH --job-name=analyze_coding_sequences
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment/nsc/dedup/individual_samples/broad/analyze_coding_sequences.err"
#SBATCH --output="logs/cpg_enrichment/nsc/dedup/individual_samples/broad/analyze_coding_sequences.out"

CELL_LINE="NSC"
PEAKS="broad"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
SCRIPT_DIR="${BASE_DIR}/scripts/cpg_enrichment"

cd $WORKING_DIR || exit 1

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

RESULTS_DIR="results/dedup/cpg_enrichment/${CELL_LINE}/${PEAKS}/mecp2_cpg_enrichment_parallel_with_flanking"

Rscript "${SCRIPT_DIR}/analyze_coding_sequences.R" \
    --work-dir "${RESULTS_DIR}" \