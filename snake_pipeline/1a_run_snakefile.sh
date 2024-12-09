#!/bin/bash
#SBATCH --job-name=Snake_pipeline
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/snake_pipeline/logs/Snake_pipeline.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/snake_pipeline/logs/Snake_pipeline.out"

# Set pipeline directory path
PIPELINE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/snake_pipeline"

# Activate the conda environment
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
mkdir -p ${PIPELINE_DIR}/results/{fastqc,trimmed,aligned,peaks,multiqc} && chmod -R 755 ${PIPELINE_DIR}/results

# Clean up any existing temporary files
find ${PIPELINE_DIR}/tmp -type f -delete
find ${PIPELINE_DIR}/results/aligned -name "*.unsorted.bam" -delete

# Unlock the working directory if needed
snakemake --unlock

# # Clean up incomplete files
# rm -f /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/snake_pipeline/results/aligned/*.unsorted.bam
# rm -f /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/snake_pipeline/results/aligned/*.bam
# rm -f /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/snake_pipeline/results/aligned/*.bam.bai

# Set the ALL_SAMPLES variable
# ALL_SAMPLES=($(ls ${DATA_DIR}/EXOGENOUS ${DATA_DIR}/ENDOGENOUS | grep '_R1_001.fastq.gz' | sed 's/_R1_001.fastq.gz//'))

# Run snakemake
snakemake \
    --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -c {threads}" \
    --cluster-config cluster.yaml \
    --jobs 32 \
    --rerun-incomplete \
    --latency-wait 60 \
    --keep-going \
    --restart-times 3 \
    --conda-frontend conda \
    --snakefile ${PIPELINE_DIR}/Snakefile \
    --configfile ${PIPELINE_DIR}/configs/config.yaml \
    2>&1 | tee ${PIPELINE_DIR}/logs/snakemake.log

# Create summary of run
snakemake --summary > "${PIPELINE_DIR}/logs/workflow_summary.txt"

# Cleanup on exit
# trap 'rm -rf "$TMPDIR"/*' EXIT

# --keep-going \
# --forcerun all

# # Run snakemake with forcerun for peak calling only
# snakemake \
#     --snakefile ${PIPELINE_DIR}/Snakefile \
#     --configfile ${PIPELINE_DIR}/configs/config.yaml \
#     --cluster "sbatch --parsable \
#         --account=kubacki.michal \
#         --mem={resources.mem_mb}MB \
#         --time={resources.time} \
#         --cpus-per-task={threads} \
#         --output=${PIPELINE_DIR}/logs/slurm/%j.out \
#         --error=${PIPELINE_DIR}/logs/slurm/%j.err" \
#     --cores all \
#     --jobs 32 \
#     --latency-wait 120 \
#     --restart-times 3 \
#     --keep-going \
#     --rerun-incomplete \
#     --until call_peaks \
#     --forcerun call_peaks \
#     all \
#     2>&1 | tee "${PIPELINE_DIR}/logs/snakemake_peaks_rerun.log"
