#!/bin/bash
#SBATCH --job-name=Snake_test
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/Snake_test.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/Snake_test.out"

# Activate environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/jupyter_nb

# Set directories
export TMPDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/tmp"
export XDG_CACHE_HOME="/beegfs/scratch/ric.broccoli/kubacki.michal/cache"
export SNAKEMAKE_PERSISTENT_CACHE="/beegfs/scratch/ric.broccoli/kubacki.michal/snakemake_cache"

# Create directories
mkdir -p "$TMPDIR" && chmod 755 "$TMPDIR"
mkdir -p "$XDG_CACHE_HOME" && chmod 755 "$XDG_CACHE_HOME"
mkdir -p "$SNAKEMAKE_PERSISTENT_CACHE" && chmod 755 "$SNAKEMAKE_PERSISTENT_CACHE"
mkdir -p logs/slurm && chmod 755 logs/slurm
mkdir -p results/{fastqc,trimmed,aligned,peaks,multiqc,fragment_sizes,extended_analysis,qc,plots,peak_analysis} && chmod -R 755 results

# Clean up
rm -rf "$TMPDIR"/*
rm -rf ~/.cache/snakemake/*


snakemake --unlock

# Run Snakemake for a single sample
SAMPLE="IgM"

snakemake \
    --snakefile Snakefile \
    --configfile config_test.yaml \
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
    --config sample=IgM \
    2>&1 | tee "logs/snakemake_test_$(date +%Y%m%d_%H%M%S).log"