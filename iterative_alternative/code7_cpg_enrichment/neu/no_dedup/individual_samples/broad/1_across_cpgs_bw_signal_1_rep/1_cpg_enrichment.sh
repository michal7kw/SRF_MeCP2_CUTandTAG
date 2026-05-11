#!/bin/bash
#SBATCH --job-name=cpg_enrichment_1_rep_in_cpg
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment/neu/no_dedup/individual_samples/broad/cpg_enrichment_1_rep_in_cpg.err"
#SBATCH --output="logs/cpg_enrichment/neu/no_dedup/individual_samples/broad/cpg_enrichment_1_rep_in_cpg.out"
#SBATCH --array=0-9  # Process in 10 chunks

ALIGNMENT="results_1b"
PEAKS="broad"
CELL_LINE="Neu"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
DATA_DIR="${BASE_DIR}/DATA"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
SCRIPT_DIR="${BASE_DIR}/scripts/cpg_enrichment/bw_signal/across_cpg_islands"

PEAKS_BASE_DIR="${WORKING_DIR}/results/no_dedup/cpg_enrichment/${CELL_LINE}/${PEAKS}"
RESULTS_DIR="${PEAKS_BASE_DIR}/cpg_enrichment_1_rep_in_cpg"

BIGWIG_DIR="${WORKING_DIR}/${ALIGNMENT}/bigwig"

cd $WORKING_DIR || exit 1

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create output directory
mkdir -p "${RESULTS_DIR}"

echo "Running CpG-focused enrichment analysis..."

# Run the Python script with peak directories
python ${SCRIPT_DIR}/cpg_enrichment_1_rep_in_cpg.py \
    --bigwig-dir ${BIGWIG_DIR} \
    --cpg-islands ${DATA_DIR}/cpg_islands.bed \
    --output-dir ${RESULTS_DIR} \
    --chunk-id $SLURM_ARRAY_TASK_ID \
    --total-chunks 10 \
    --exo-peaks-dir ${PEAKS_BASE_DIR}/exo \
    --endo-peaks-dir ${PEAKS_BASE_DIR}/endo \
    --cell-type ${CELL_LINE}