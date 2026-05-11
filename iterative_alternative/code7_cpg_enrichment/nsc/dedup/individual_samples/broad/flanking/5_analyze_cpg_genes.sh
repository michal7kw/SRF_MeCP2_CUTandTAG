#!/bin/bash
#SBATCH --job-name=analyze_cpg_genes
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment/nsc/dedup/individual_samples/broad/analyze_cpg_genes_with_flanking.err"
#SBATCH --output="logs/cpg_enrichment/nsc/dedup/individual_samples/broad/analyze_cpg_genes_with_flanking.out"

CELL_LINE="NSC"
PEAKS="broad"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
SCRIPT_DIR="${BASE_DIR}/scripts/cpg_enrichment"

cd $WORKING_DIR || exit 1

RESULTS_DIR="results/dedup/cpg_enrichment/${CELL_LINE}/${PEAKS}/mecp2_cpg_enrichment_parallel_with_flanking"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

Rscript "${SCRIPT_DIR}/analyze_cpg_genes.R" \
    --work-dir "${RESULTS_DIR}" 