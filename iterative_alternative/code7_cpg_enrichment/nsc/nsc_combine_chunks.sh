#!/bin/bash
#SBATCH --job-name=nsc_combine_chunks
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/nsc_combine_chunks.err"
#SBATCH --output="logs/nsc_combine_chunks.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

EXPERIMENT_OUT="align2_005_consistent_peaks"
CELL_LINE="NSC"
RESULTS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment_${EXPERIMENT_OUT}/${CELL_LINE}"

# Create logs directory if it doesn't exist
mkdir -p logs

# Run the combine_chunks.py script
python ../scripts/cpg_enrichment/combine_chunks.py \
    --results-dir "${RESULTS_DIR}/mecp2_cpg_enrichment_parallel"

echo "Chunk combination complete" 