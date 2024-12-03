#!/bin/bash
#SBATCH --job-name=analyze_enrichment_NSC
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/analyze_enrichment_NSC.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/analyze_enrichment_NSC.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda deactivate
conda activate jupyter_nb

# Create logs directory if it doesn't exist
mkdir -p logs

# Run the script with full error traceback
python -u scripts/analyze_enrichment_NSC.py 2>&1 | tee "logs/analyze_enrichment_NSC.out"