#!/bin/bash
#SBATCH --job-name=process_cpg_enrichment_NSC_5
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/process_cpg_enrichment_NSC_5_align2_005.err"
#SBATCH --output="logs/process_cpg_enrichment_NSC_5_align2_005.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

mkdir -p logs

# Define base directory
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"

# Define directory structure
WORKING_DIR="${BASE_DIR}/iterative_alternative"
DATA_DIR="${BASE_DIR}/DATA"
RESULTS_DIR="${WORKING_DIR}/results_5_align2_005/peaks_enrichment_and_genes_regulation_NSC_5"

# Check if the input file exists
if [ ! -f "${RESULTS_DIR}/cpg_enrichment_NSC.csv" ]; then
    echo "Error: Input file ${RESULTS_DIR}/cpg_enrichment_NSC.csv not found"
    echo "Please run peaks_enrichment_and_genes_regulation_NSC_5.py first"
    exit 1
fi

# Check if GTF file exists
if [ ! -f "${DATA_DIR}/gencode.vM10.annotation.gtf" ]; then
    echo "Error: GTF file ${DATA_DIR}/gencode.vM10.annotation.gtf not found"
    exit 1
fi

# Run analysis script
python -u ../scripts/peaks_enrichment_and_genes_regulation/process_cpg_enrichment.py \
    --working-dir "${WORKING_DIR}" \
    --data-dir "${DATA_DIR}" \
    --results-dir "${RESULTS_DIR}" \
    2>&1 | tee "logs/process_cpg_enrichment_NSC_5_align2_005.out"

# Check if the output file was created
if [ ! -f "${RESULTS_DIR}/cpg_enrichment_annotated.csv" ]; then
    echo "Error: Output file was not created"
    exit 1
else
    echo "Successfully created ${RESULTS_DIR}/cpg_enrichment_annotated.csv"
fi

