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

# Create output directory for MultiQC
mkdir -p results_1/multiqc

# Run MultiQC
multiqc \
    --force \
    --outdir results_1/multiqc \
    --filename multiqc_report \
    --title "CUT&Tag QC Report" \
    --comment "Quality control metrics for CUT&Tag data" \
    results_1/aligned \
    results_1/peaks \
    results_1/qc \
    logs