#!/bin/bash
#SBATCH --job-name=mecp2_integrate
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/mecp2_integrate.err"
#SBATCH --output="logs/mecp2_integrate.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create logs directory if it doesn't exist
mkdir -p logs

WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative"
DATA_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/DATA"
RESULTS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment"

# Ensure enrichment analysis is complete
ENRICHMENT_FILE="${RESULTS_DIR}/mecp2_cpg_enrichment_parallel/mecp2_cpg_enrichment_parallel.csv"
if [ ! -f "$ENRICHMENT_FILE" ]; then
    echo "Error: Enrichment analysis results not found at $ENRICHMENT_FILE"
    exit 1
fi

# Run the integration script
python ../scripts/integrate_rna_seq.py \
    --enrichment-file ${ENRICHMENT_FILE} \
    --rna-seq-file ${DATA_DIR}/DEA_NSC.csv \
    --gene-annotations ${DATA_DIR}/gencode.vM10.annotation.gtf \
    --output-dir ${RESULTS_DIR}/integrated/

echo "Integration complete. Check results in ${RESULTS_DIR}/integrated/"
