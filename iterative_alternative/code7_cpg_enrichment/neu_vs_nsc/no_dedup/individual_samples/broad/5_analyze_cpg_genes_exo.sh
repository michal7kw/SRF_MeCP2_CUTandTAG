#!/bin/bash
#SBATCH --job-name=analyze_cpg_genes_neu_vs_nsc_exo
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment/neu_vs_nsc_exo/analyze_cpg_genes_neu_vs_nsc_exo.err"
#SBATCH --output="logs/cpg_enrichment/neu_vs_nsc_exo/analyze_cpg_genes_neu_vs_nsc_exo.out"

CELL_LINE="NSC"
PEAKS="broad"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
SCRIPT_DIR="${BASE_DIR}/scripts/cpg_enrichment"

cd $WORKING_DIR || exit 1

RESULTS_DIR="results/no_dedup/cpg_enrichment/neu_vs_nsc_exo"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

Rscript "${SCRIPT_DIR}/analyze_cpg_genes_neu_vs_nsc_exo.R" \
    --work-dir "${RESULTS_DIR}" 