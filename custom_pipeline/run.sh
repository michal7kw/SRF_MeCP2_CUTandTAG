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

# Load required modules
module load trimgalore/0.6.6
module load fastqc/0.11.9
module load bowtie2/2.4.2
module load samtools/1.13
module load macs2/2.2.7.1
module load R/4.1.0

# Set temporary directories
export TMPDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/tmp"
export XDG_CACHE_HOME="/beegfs/scratch/ric.broccoli/kubacki.michal/cache"
export SNAKEMAKE_PERSISTENT_CACHE="/beegfs/scratch/ric.broccoli/kubacki.michal/snakemake_cache"

# Create necessary directories with proper permissions
mkdir -p "$TMPDIR" && chmod 755 "$TMPDIR"
mkdir -p "$XDG_CACHE_HOME" && chmod 755 "$XDG_CACHE_HOME"
mkdir -p "$SNAKEMAKE_PERSISTENT_CACHE" && chmod 755 "$SNAKEMAKE_PERSISTENT_CACHE"
mkdir -p logs/slurm && chmod 755 logs/slurm
mkdir -p results/{fastqc,trimmed,aligned,peaks,multiqc,fragment_sizes,extended_analysis} && chmod -R 755 results

# Clean up temporary files
rm -rf "$TMPDIR"/*
rm -rf ~/.cache/snakemake/*

# Activate the jupyter_nb environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/jupyter_nb

# Unlock the working directory if needed
snakemake --unlock

# Ensure modules are available
module purge
module load trimgalore/0.6.6
module load fastqc/0.11.9
module load bowtie2/2.4.2
module load samtools/1.13
module load macs2/2.2.7.1
module load R/4.1.0

# Verify module loading
command -v trim_galore >/dev/null 2>&1 || { echo "trim_galore not found"; exit 1; }
command -v fastqc >/dev/null 2>&1 || { echo "fastqc not found"; exit 1; }

# Run Snakemake
snakemake \
    --snakefile Snakefile \
    --cluster "sbatch --parsable \
        --account={cluster.account} \
        --mem={cluster.mem} \
        --time={cluster.time} \
        --cpus-per-task={threads} \
        --output=logs/slurm/%j.out \
        --error=logs/slurm/%j.err" \
    --cluster-config cluster.json \
    --cores all \
    --jobs 64 \
    --latency-wait 120 \
    --restart-times 3 \
    --keep-going \
    --rerun-incomplete \
    --use-envmodules \
    --default-resources \
        mem_mb=8000 \
        runtime=240 \
    2>&1 | tee "logs/snakemake_$(date +%Y%m%d_%H%M%S).log"

# Create summary of run
snakemake --summary > "logs/workflow_summary_$(date +%Y%m%d_%H%M%S).txt"

# Cleanup on exit
trap 'rm -rf "$TMPDIR"/*' EXIT