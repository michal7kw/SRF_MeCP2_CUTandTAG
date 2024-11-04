#!/bin/bash
#SBATCH --job-name=Snake_pipeline
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_ChipSeq/custom_pipeline/logs/Snake_pipeline.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_ChipSeq/custom_pipeline/logs/Snake_pipeline.out"

# Load the appropriate conda environment (if needed)
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Define the cluster submission command

# Run Snakemake
# snakemake --unlock
snakemake --snakefile Snakefile --latency-wait 60 --configfile config.yaml --use-conda -p --printshellcmds --keep-going --rerun-incomplete
# snakemake -n --snakefile Snakefile -p --verbose --configfile config.yaml





# --ntasks=256
# --nodes=8
# --ntasks-per-node=32

# snakemake --unlock --snakefile Snakefile --configfile config.yaml --use-conda --cores all -p
# snakemake --unlock --snakefile Snakefile --configfile config.yaml --use-conda --cores all -p --printshellcmds --keep-going --rerun-incomplete
# snakemake --snakefile Snakefile --configfile config.yaml --use-conda --cores all -p --printshellcmds --keep-going --rerun-incomplete

# snakemake --unlock --snakefile Snakefile --configfile config.yaml --use-conda \
#     --cluster "srun --nodes=1 --ntasks=1 --cpus-per-task={threads}" \
#     --jobs 128 -p

# snakemake --unlock