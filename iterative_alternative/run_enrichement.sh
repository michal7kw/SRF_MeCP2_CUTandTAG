#!/bin/bash
#SBATCH --job-name=analyze_cpg_enrichment
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/analyze_cpg_enrichment.err"
#SBATCH --output="logs/analyze_cpg_enrichment.out"

EXPERIMENT="align1_005"
CELL_LINE="NSC"

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment_${EXPERIMENT}/${CELL_LINE}

source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Run the Python script with the experiment name
python analyze_cpg_enrichment.py "$EXPERIMENT" "$CELL_LINE"