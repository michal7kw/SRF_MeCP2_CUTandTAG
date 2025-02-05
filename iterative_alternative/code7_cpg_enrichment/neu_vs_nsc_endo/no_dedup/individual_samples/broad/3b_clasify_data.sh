#!/bin/bash
#SBATCH --job-name=Enrichment_analysis_intra_cell
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment/neu_vs_nsc_endo/Enrichment_analysis_intra_cell.err"
#SBATCH --output="logs/cpg_enrichment/neu_vs_nsc_endo/Enrichment_analysis_intra_cell.out"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
EXPERIMENT="neu_vs_nsc_endo/no_dedup/individual_samples/broad"
SCRIPT_DIR="${BASE_DIR}/iterative_alternative/code7_cpg_enrichment/${EXPERIMENT}"
cd $WORKING_DIR || exit 1

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

python ${SCRIPT_DIR}/Enrichment_analysis_intra_cell.py

