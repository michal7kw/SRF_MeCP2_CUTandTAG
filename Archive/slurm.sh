
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