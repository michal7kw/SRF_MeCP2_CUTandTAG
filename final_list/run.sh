#!/bin/bash
#SBATCH --job-name=CutTag_Pipeline
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/final_list/logs/Snake_pipeline.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/final_list/logs/Snake_pipeline.out"

# Set pipeline directory path
PIPELINE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/final_list"
export PIPELINE_DIR

# Activate the conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Install required executor plugin if not already installed
# pip install snakemake-executor-plugin-slurm

# Set temporary directories
export TMPDIR="${PIPELINE_DIR}/tmp"
export XDG_CACHE_HOME="/beegfs/scratch/ric.broccoli/kubacki.michal/cache"
export SNAKEMAKE_PERSISTENT_CACHE="/beegfs/scratch/ric.broccoli/kubacki.michal/snakemake_cache"

# # Create necessary directories
# mkdir -p "$TMPDIR" && chmod 755 "$TMPDIR"
# mkdir -p "$XDG_CACHE_HOME" && chmod 755 "$XDG_CACHE_HOME"
# mkdir -p "$SNAKEMAKE_PERSISTENT_CACHE" && chmod 755 "$SNAKEMAKE_PERSISTENT_CACHE"
# mkdir -p ${PIPELINE_DIR}/logs/slurm && chmod 755 ${PIPELINE_DIR}/logs/slurm
# mkdir -p ${PIPELINE_DIR}/results/{fastqc,trimmed,aligned,peaks/{macs2,seacr,consensus,final},qc,multiqc} \
#     && chmod -R 755 ${PIPELINE_DIR}/results

# Clean up any existing temporary files
# find ${PIPELINE_DIR}/tmp -type f -delete
# find ${PIPELINE_DIR}/results/aligned -name "*.unsorted.bam" -delete

# Unlock the working directory if needed
snakemake --unlock

# Run snakemake with SLURM executor
snakemake \
    --executor slurm \
    --jobs 32 \
    --default-resources \
        "slurm_time='4:00:00'" \
        "mem_mb=16000" \
        "threads=4" \
    --resources "mem_mb=500000" \
    --set-resources "align:mem_mb=96000" \
    --set-resources "macs2_peaks:mem_mb=32000" \
    --set-resources "seacr_peaks:mem_mb=32000" \
    --set-threads "align=32" \
    --set-threads "macs2_peaks=16" \
    --set-threads "seacr_peaks=4" \
    --rerun-incomplete \
    --latency-wait 60 \
    --conda-frontend conda \
    --snakefile ${PIPELINE_DIR}/Snakefile \
    2>&1 | tee ${PIPELINE_DIR}/logs/snakemake.log

# Create summary of run
snakemake --summary > "${PIPELINE_DIR}/logs/workflow_summary.txt"


    # --keep-going \
    # --restart-times 3 \
