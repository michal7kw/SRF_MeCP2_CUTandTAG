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
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/jupyter_nb

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

# Clean up temporary files
rm -rf "$TMPDIR"/*
rm -rf ~/.cache/snakemake/*

# Activate the jupyter_nb environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/jupyter_nb

# Unlock the working directory if needed
snakemake --unlock

# Set the ALL_SAMPLES variable
ALL_SAMPLES=($(ls DATA/EXOGENOUS DATA/ENDOGENOUS | grep '_R1_001.fastq.gz' | sed 's/_R1_001.fastq.gz//'))

# Run Snakemake with environment activation for each job
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
    $(for sample in ${ALL_SAMPLES[@]}; do echo "results/peaks/${sample}_peaks.narrowPeak"; done) \
    2>&1 | tee "logs/snakemake_$(date +%Y%m%d_%H%M%S).log"

# Create summary of run
snakemake --summary > "logs/workflow_summary_$(date +%Y%m%d_%H%M%S).txt"

# Cleanup on exit
trap 'rm -rf "$TMPDIR"/*' EXIT