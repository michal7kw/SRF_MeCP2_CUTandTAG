#!/bin/bash
#SBATCH --job-name=nsc_integrate_align2_005_consistent_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/nsc_integrate_align2_005_consistent_peaks.err"
#SBATCH --output="logs/nsc_integrate_align2_005_consistent_peaks.out"

# Set base directories
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
WORKING_DIR="${BASE_DIR}/iterative_alternative"
DATA_DIR="${BASE_DIR}/DATA"
EXPERIMENT="align2_005_consistent_peaks"
CELL_LINE="NSC"
RESULTS_DIR="${WORKING_DIR}/analyze_mecp2_cpg_enrichment_${EXPERIMENT}/${CELL_LINE}"

# Define input files
ENRICHMENT_FILE="${RESULTS_DIR}/mecp2_cpg_enrichment_parallel/mecp2_cpg_enrichment_parallel.csv"
RNA_SEQ_FILE="${WORKING_DIR}/DATA/DEA_${CELL_LINE}_filtered.csv"
GTF_FILE="${DATA_DIR}/gencode.vM10.annotation.gtf"

cd ${WORKING_DIR}

# Create necessary directories
mkdir -p logs
mkdir -p ${RESULTS_DIR}/integrated

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

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
python ../scripts/cpg_enrichment_integration/integrate_rna_seq_only_consistent_peaks.py \
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

#####################################################################
# This script integrates MeCP2 enrichment data with RNA-seq results
# to analyze relationships between MeCP2 binding and gene expression
#
# Input files:
# - Enrichment analysis results from analyze_mecp2_cpg_enrichment.py:
#   Contains CpG regions with MeCP2 binding data including:
#   - Enrichment scores
#   - Signal values for exo/endo samples
#   - Binding type classifications
#   - Statistical significance
#
# - RNA-seq differential expression results containing:
#   - Gene names
#   - Base mean expression
#   - Log2 fold changes
#   - Adjusted p-values
#
# - Gene annotations from GTF file with:
#   - Gene coordinates
#   - Gene names and IDs
#
# The script:
# 1. Validates input files
# 2. Creates output directory structure
# 3. Runs integrate_rna_seq_only_consistent_peaks.py which:
#    - Maps enriched CpG regions to nearby genes
#    - Integrates with RNA-seq expression data
#    - Categorizes genes by expression changes
#    - Generates enrichment distribution plots
#    - Performs statistical analyses
#
# Output files in ${OUTPUT_DIR}:
# - mecp2_enriched_genes.csv: Full integration results with:
#   - Gene information
#   - CpG coordinates
#   - Enrichment scores
#   - Expression changes
# - mecp2_enriched_genes_only.csv: Simplified list of enriched genes
# - Visualization plots:
#   - Enrichment distributions
#   - Expression correlation plots
#   - Category-based analyses
#####################################################################

