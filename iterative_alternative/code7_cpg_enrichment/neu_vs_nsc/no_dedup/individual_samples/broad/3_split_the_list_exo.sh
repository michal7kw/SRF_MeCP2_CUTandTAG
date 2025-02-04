#!/bin/bash
#SBATCH --job-name=split_the_list_neu_vs_nsc_exo
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment/neu_vs_nsc_exo/split_the_list_neu_vs_nsc_exo.err"
#SBATCH --output="logs/cpg_enrichment/neu_vs_nsc_exo/split_the_list_neu_vs_nsc_exo.out"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
SCRIPT_DIR="${BASE_DIR}/scripts/cpg_enrichment"

cd $WORKING_DIR || exit 1

# Set paths for peak files and results
PEAKS_BASE_DIR="${WORKING_DIR}/results/no_dedup/cpg_enrichment"
RESULTS_DIR="${PEAKS_BASE_DIR}/neu_vs_nsc_exo"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

python ${SCRIPT_DIR}/split_the_list_neu_vs_nsc_exo.py \
    --work-dir "${RESULTS_DIR}"