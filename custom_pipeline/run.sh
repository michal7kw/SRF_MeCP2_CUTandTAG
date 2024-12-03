#!/bin/bash
#SBATCH --job-name=Snake_pipeline
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/Snake_pipeline.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/Snake_pipeline.out"

# Activate the jupyter_nb environment first
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Set temporary directories
export TMPDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/tmp"
export XDG_CACHE_HOME="/beegfs/scratch/ric.broccoli/kubacki.michal/cache"
export SNAKEMAKE_PERSISTENT_CACHE="/beegfs/scratch/ric.broccoli/kubacki.michal/snakemake_cache"

# Create necessary directories with proper permissions
mkdir -p "$TMPDIR" && chmod 755 "$TMPDIR"
mkdir -p "$XDG_CACHE_HOME" && chmod 755 "$XDG_CACHE_HOME"
mkdir -p "$SNAKEMAKE_PERSISTENT_CACHE" && chmod 755 "$SNAKEMAKE_PERSISTENT_CACHE"
mkdir -p logs/slurm && chmod 755 logs/slurm
mkdir -p results/{fastqc,trimmed,aligned,peaks,multiqc,fragment_sizes,extended_analysis,qc,plots,peak_analysis} && chmod -R 755 results

# # Clean up temporary files
# rm -rf "$TMPDIR"/*
# rm -rf ~/.cache/snakemake/*

# Unlock the working directory if needed
snakemake --unlock

# Set the ALL_SAMPLES variable
ALL_SAMPLES=($(ls DATA/EXOGENOUS DATA/ENDOGENOUS | grep '_R1_001.fastq.gz' | sed 's/_R1_001.fastq.gz//'))

# Run snakemake with forcerun for peak calling only
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --cluster-config cluster.json \
    --cluster "sbatch --parsable \
        --account={cluster.account} \
        --mem={cluster.mem} \
        --time={cluster.time} \
        --cpus-per-task={threads} \
        --output=logs/slurm/%j.out \
        --error=logs/slurm/%j.err" \
    --cores all \
    --jobs 32 \
    --latency-wait 120 \
    --restart-times 3 \
    --keep-going \
    --rerun-incomplete \
    --use-envmodules \
    --resources skip_fastqc=1 skip_trim=1 skip_align=1 \
    --until call_peaks \
    --forcerun call_peaks \
    all \
    2>&1 | tee "logs/snakemake_peaks_rerun.log"

# Create summary of run
snakemake --summary > "logs/workflow_summary.txt"

# Cleanup on exit
trap 'rm -rf "$TMPDIR"/*' EXIT