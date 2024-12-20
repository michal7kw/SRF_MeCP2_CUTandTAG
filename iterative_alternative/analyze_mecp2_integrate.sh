#!/bin/bash
#SBATCH --job-name=mecp2_cpg_enrichment
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/mecp2_cpg_enrichment.err"
#SBATCH --output="logs/mecp2_cpg_enrichment.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create logs directory if it doesn't exist
mkdir -p logs

WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative"
DATA_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/DATA"
RESULTS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment"

# Run the script with working directory argument and full error traceback
python integrate_rna_seq.py \
    --enrichment-file ${RESULTS_DIR}/mecp2_cpg_enrichment.csv \
    --rna-seq-file ${RESULTS_DIR}/rna_seq_results.csv \
    --gene-annotations ${DATA_DIR}/gencode.vM10.annotation.gtf \
    --output-dir ${RESULTS_DIR}/integrated/

# peaks were copied to results_5 from results_2_new_001
