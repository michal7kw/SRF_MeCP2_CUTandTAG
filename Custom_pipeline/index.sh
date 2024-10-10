#!/bin/bash
#SBATCH --job-name=chipseq
#SBATCH --account kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --ntasks=32
#SBATCH --nodes=2
#SBATCH --mail-type=ALL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=kubacki.michal@hst.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/logs/chipseq.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/logs/chipseq.out"

# Load the appropriate conda environment (if needed)
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

bowtie2-build mm10.fa mm10_bowtie2_index/mm10
