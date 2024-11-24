#!/bin/bash
#SBATCH --job-name=Snake_pipeline
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/Snake_pipeline_%j.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/Snake_pipeline_%j.out"

# Set up conda cache and remove any stale locks
export CONDA_PKGS_DIRS="/beegfs/scratch/ric.broccoli/kubacki.michal/conda_pkgs"
mkdir -p "$CONDA_PKGS_DIRS"
rm -rf "$CONDA_PKGS_DIRS/cache.lock"
rm -rf "$CONDA_PKGS_DIRS/locks"
mkdir -p "$CONDA_PKGS_DIRS/locks"

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Set temporary directories
export TMPDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/tmp"
export XDG_CACHE_HOME="/beegfs/scratch/ric.broccoli/kubacki.michal/cache"
export SNAKEMAKE_PERSISTENT_CACHE="/beegfs/scratch/ric.broccoli/kubacki.michal/snakemake_cache"

# Create necessary directories
mkdir -p "$TMPDIR"
mkdir -p "$XDG_CACHE_HOME"
mkdir -p "$SNAKEMAKE_PERSISTENT_CACHE"
mkdir -p logs/slurm

# Clean up temporary files
rm -rf "$TMPDIR"/*
rm -rf ~/.cache/snakemake/*

# Unlock the working directory if needed
snakemake --unlock

# Run Snakemake
snakemake \
    --snakefile Snakefile \
    --cluster "sbatch --parsable --account={cluster.account} --mem={cluster.mem} --time={cluster.time} --cpus-per-task={threads}" \
    --cluster-config cluster.json \
    --cores all \
    --jobs 64 \
    --latency-wait 120 \
    --restart-times 3 \
    --keep-going \
    --rerun-incomplete \
    --use-conda \
    --default-resources \
        mem_mb=8000 \
        runtime=240 \
    2>&1 | tee "logs/snakemake_$(date +%Y%m%d_%H%M%S).log"

# Create summary of run
snakemake --summary > "logs/workflow_summary_$(date +%Y%m%d_%H%M%S).txt"

# Cleanup on exit
trap 'rm -rf "$TMPDIR"/*' EXIT