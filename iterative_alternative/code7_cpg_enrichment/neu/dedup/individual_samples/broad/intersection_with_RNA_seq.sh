#!/bin/bash
#SBATCH --job-name=intersection_with_RNA_seq
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/intersection_with_RNA_seq.err"
#SBATCH --output="logs/intersection_with_RNA_seq.out"

# Create logs directory if it doesn't exist
mkdir -p logs

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

SCRIPT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/scripts/cpg_enrichment"
WORK_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment_align2_005_consistent_peaks/NSC/mecp2_cpg_enrichment_parallel"

# Ensure the work directory exists
mkdir -p "${WORK_DIR}"

# Run R script with explicit argument names (space between option and value)
Rscript "${SCRIPT_DIR}/intersection_with_RNA_seq.R" \
    --work-dir "${WORK_DIR}" \
    2>&1 | tee "logs/intersection_with_RNA_seq.out" 
