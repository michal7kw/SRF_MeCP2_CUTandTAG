#!/bin/bash
#SBATCH --job-name=alternative_gene_enrichment_NSC_combined
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/alternative_gene_enrichment_NSC_combined.err"
#SBATCH --output="logs/alternative_gene_enrichment_NSC_combined.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create required directories
mkdir -p logs

# Define directory structure
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
DATA_DIR="${BASE_DIR}/DATA"

# Define experiment-specific directories
ALIGN_DIR="${WORKING_DIR}/results_1b"
PEAKS_DIR="${WORKING_DIR}/results_2_align2_005"
RESULTS_DIR="${WORKING_DIR}/results_5_align2_005/alternative_gene_enrichment_NSC_combined"

# Create results directory
mkdir -p "${RESULTS_DIR}"

# Define input files (fix the array declaration)
EXO_PEAKS="${PEAKS_DIR}/peaks/narrow/NPCs_exo_combined.narrowPeak"
ENDO_PEAKS="${PEAKS_DIR}/peaks/narrow/NPCs_endo_combined.narrowPeak"
DEA_FILE="${DATA_DIR}/DEA_NSC.csv"
GTF_FILE="${DATA_DIR}/gencode.vM10.annotation.gtf"

# Run enrichment analysis (use the variables directly)
python -u ../scripts/alternative_gene_enrichment/alternative_gene_enrichment_NSC_combined.py \
    --exo-peaks "${EXO_PEAKS}" \
    --endo-peaks "${ENDO_PEAKS}" \
    --dea "${DEA_FILE}" \
    --gtf "${GTF_FILE}" \
    --output-dir "${RESULTS_DIR}" \
    --promoter-window 2000 \
    2>&1 | tee "logs/alternative_gene_enrichment_NSC_combined.out"

# Check exit status
if [ $? -eq 0 ]; then
    echo "Enrichment analysis completed successfully"
else
    echo "Error: Enrichment analysis failed"
    exit 1
fi
