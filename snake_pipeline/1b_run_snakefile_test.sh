#!/bin/bash
#SBATCH --job-name=Snake_test
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/snake_pipeline/logs/Snake_test.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/snake_pipeline/logs/Snake_test.out"

# Set pipeline directory path
PIPELINE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/snake_pipeline"
DATA_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/DATA"

# Activate environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Set temporary directories
export TMPDIR="${PIPELINE_DIR}/tmp"
export XDG_CACHE_HOME="/beegfs/scratch/ric.broccoli/kubacki.michal/cache"
export SNAKEMAKE_PERSISTENT_CACHE="/beegfs/scratch/ric.broccoli/kubacki.michal/snakemake_cache"

# Create necessary directories with proper permissions
mkdir -p "$TMPDIR" && chmod 755 "$TMPDIR"
mkdir -p "$XDG_CACHE_HOME" && chmod 755 "$XDG_CACHE_HOME"
mkdir -p "$SNAKEMAKE_PERSISTENT_CACHE" && chmod 755 "$SNAKEMAKE_PERSISTENT_CACHE"
mkdir -p ${PIPELINE_DIR}/logs/slurm && chmod 755 ${PIPELINE_DIR}/logs/slurm
mkdir -p ${PIPELINE_DIR}/results/{fastqc,trimmed,aligned,peaks,multiqc,fragment_sizes,extended_analysis,qc,plots,peak_analysis} && chmod -R 755 ${PIPELINE_DIR}/results

# Clean up temporary files
rm -rf "$TMPDIR"/*
rm -rf ~/.cache/snakemake/*

# Unlock the working directory if needed
snakemake --unlock

# Set test sample
TEST_SAMPLE="IgM"  # You can change this to test different samples

# Run Snakemake for a single test sample
snakemake \
    --snakefile ${PIPELINE_DIR}/Snakefile \
    --configfile ${PIPELINE_DIR}/configs/config_test.yaml \
    --cluster "sbatch --parsable \
        --account=kubacki.michal \
        --mem={resources.mem_mb}MB \
        --time={resources.time} \
        --cpus-per-task={threads} \
        --output=${PIPELINE_DIR}/logs/slurm/%j.out \
        --error=${PIPELINE_DIR}/logs/slurm/%j.err" \
    --cores all \
    --jobs 32 \
    --latency-wait 120 \
    --restart-times 3 \
    --keep-going \
    --rerun-incomplete \
    --config sample=$TEST_SAMPLE \
    2>&1 | tee "${PIPELINE_DIR}/logs/snakemake_test_${TEST_SAMPLE}_$(date +%Y%m%d_%H%M%S).log"

# Create summary of test run
snakemake --summary > "${PIPELINE_DIR}/logs/workflow_summary_test.txt"

# Cleanup on exit
trap 'rm -rf "$TMPDIR"/*' EXIT