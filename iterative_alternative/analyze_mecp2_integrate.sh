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

# Set base directories
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
DATA_DIR="${BASE_DIR}/DATA"
RESULTS_DIR="${WORKING_DIR}/analyze_mecp2_cpg_enrichment"

cd ${WORKING_DIR}

# Create necessary directories
mkdir -p logs
mkdir -p ${RESULTS_DIR}/integrated

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define input files
ENRICHMENT_FILE="${RESULTS_DIR}/mecp2_cpg_enrichment_parallel/mecp2_cpg_enrichment_parallel.csv"
RNA_SEQ_FILE="${WORKING_DIR}/analyze_mecp2_cpg_enrichment/integrated/DEA_NEU_filtered.csv"
GTF_FILE="${DATA_DIR}/gencode.vM10.annotation.gtf"

# Validate input files
echo "Validating input files..."

if [ ! -f "$ENRICHMENT_FILE" ]; then
    echo "Error: Enrichment analysis results not found at $ENRICHMENT_FILE"
    exit 1
fi

if [ ! -f "$RNA_SEQ_FILE" ]; then
    echo "Error: RNA-seq results not found at $RNA_SEQ_FILE"
    exit 1
fi

if [ ! -f "$GTF_FILE" ]; then
    echo "Error: Gene annotations not found at $GTF_FILE"
    exit 1
fi

echo "Found all required input files:"
echo "- Enrichment data: $ENRICHMENT_FILE"
echo "- RNA-seq data: $RNA_SEQ_FILE"
echo "- Gene annotations: $GTF_FILE"

# Create output directory
OUTPUT_DIR="${RESULTS_DIR}/integrated"
mkdir -p ${OUTPUT_DIR}

echo "Starting integration analysis..."
echo "Output will be saved to: ${OUTPUT_DIR}"

# Run the integration script
python ../scripts/integrate_rna_seq.py \
    --enrichment-file ${ENRICHMENT_FILE} \
    --rna-seq-file ${RNA_SEQ_FILE} \
    --gene-annotations ${GTF_FILE} \
    --output-dir ${OUTPUT_DIR}

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "Integration completed successfully"
    echo "Results are available in: ${OUTPUT_DIR}"
    
    # List generated files
    echo -e "\nGenerated files:"
    ls -lh ${OUTPUT_DIR}
else
    echo "Error: Integration script failed"
    exit 1
fi
