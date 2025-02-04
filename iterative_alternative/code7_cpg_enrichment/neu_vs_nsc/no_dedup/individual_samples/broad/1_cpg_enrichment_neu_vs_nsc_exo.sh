#!/bin/bash
#SBATCH --job-name=cpg_enrichment_neu_vs_nsc_exo
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment/neu_vs_nsc_exo/cpg_enrichment_neu_vs_nsc_exo.err"
#SBATCH --output="logs/cpg_enrichment/neu_vs_nsc_exo/cpg_enrichment_neu_vs_nsc_exo.out"
#SBATCH --array=0-9  # Process in 10 chunks

ALIGNMENT="results_1b"
PEAKS="broad"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
DATA_DIR="${BASE_DIR}/DATA"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
SCRIPT_DIR="${BASE_DIR}/scripts/cpg_enrichment/bw_signal/across_peaks"

# Set paths for peak files and results
PEAKS_BASE_DIR="${WORKING_DIR}/results/no_dedup/cpg_enrichment"
RESULTS_DIR="${PEAKS_BASE_DIR}/neu_vs_nsc_exo_comparison"

BIGWIG_DIR="${WORKING_DIR}/${ALIGNMENT}/bigwig"

cd $WORKING_DIR || exit 1

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create output directories
mkdir -p "${RESULTS_DIR}"
mkdir -p "logs/cpg_enrichment/neu_vs_nsc_exo"

echo "Running Neuron vs NSC exo comparison enrichment analysis..."

# Run the Python script with peak directories
python ${SCRIPT_DIR}/cpg_enrichment_neu_vs_nsc_exo.py \
    --bigwig-dir ${BIGWIG_DIR} \
    --cpg-islands ${DATA_DIR}/cpg_islands.bed \
    --output-dir ${RESULTS_DIR} \
    --chunk-id $SLURM_ARRAY_TASK_ID \
    --total-chunks 10 \
    --neu-peaks ${PEAKS_BASE_DIR}/Neu/${PEAKS}/exo \
    --nsc-peaks ${PEAKS_BASE_DIR}/NSC/${PEAKS}/exo \
    --min-signal 0.05 \
    --min-fold-change 1.2 \
    --max-qvalue 0.1

echo "Analysis complete for chunk $SLURM_ARRAY_TASK_ID"
