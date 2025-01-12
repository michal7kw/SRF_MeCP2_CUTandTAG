#!/bin/bash
#SBATCH --job-name=peaks_enrichment_and_genes_regulation_NSC_5_align2_005_combined
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_enrichment_and_genes_regulation_NSC_5_align2_005_combined.err"
#SBATCH --output="logs/peaks_enrichment_and_genes_regulation_NSC_5_align2_005_combined.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

mkdir -p logs

# Define base directory
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"

# Define directories relative to base
WORKING_DIR="${BASE_DIR}/iterative_alternative"
DATA_DIR="${WORKING_DIR}/results_1b"
INPUT_DIR="${WORKING_DIR}/results_2_align2_005"
RESULTS_DIR="${WORKING_DIR}/results_5_align2_005"

# Create results directory if it doesn't exist
mkdir -p "${RESULTS_DIR}/enrichment_analysis"

# Define input files
EXO_PEAKS="${INPUT_DIR}/peaks/narrow/NPCs_exo_combined.narrowPeak"
ENDO_PEAKS="${INPUT_DIR}/peaks/narrow/NPCs_endo_combined.narrowPeak"
DEA_FILE="${WORKING_DIR}/DATA/DEA_NSC.csv"
GTF_FILE="${BASE_DIR}/DATA/gencode.vM10.annotation.gtf"

# Run the enrichment analysis
python -u ../scripts/peaks_enrichment_and_genes_regulation/peaks_enrichment_and_genes_regulation_NSC_5_combined.py \
    --exo-peaks "${EXO_PEAKS}" \
    --endo-peaks "${ENDO_PEAKS}" \
    --dea "${DEA_FILE}" \
    --gtf "${GTF_FILE}" \
    --output-dir "${RESULTS_DIR}/enrichment_analysis" \
    --promoter-window 2000 \
    2>&1 | tee "logs/peaks_enrichment_and_genes_regulation_NSC_5_combined.out"

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "Enrichment analysis completed successfully"
else
    echo "Error: Enrichment analysis failed"
    exit 1
fi
