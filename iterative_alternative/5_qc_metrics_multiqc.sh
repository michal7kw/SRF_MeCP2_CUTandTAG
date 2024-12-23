#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/multiqc.err"
#SBATCH --output="logs/multiqc.out"

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

ALIGNED_DIR="results_1b/aligned"
PEAKS_DIR="results_2_new_005_align2/peaks"
QC_DIR="results_1/qc"
LOGS_DIR="results_1/logs"
RESULTS_DIR="results_1"

# Create output directory for MultiQC
mkdir -p ${RESULTS_DIR}/multiqc

# Run MultiQC
multiqc \
    --force \
    --outdir ${RESULTS_DIR}/multiqc \
    --filename multiqc_report \
    --title "CUT&Tag QC Report" \
    --comment "Quality control metrics for CUT&Tag data" \
    ${ALIGNED_DIR} \
    ${PEAKS_DIR} \
    ${QC_DIR} \
    ${LOGS_DIR}
