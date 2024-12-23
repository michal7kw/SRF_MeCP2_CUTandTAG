#!/bin/bash
#SBATCH --job-name=cpg_enrichment_serial
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/cpg_enrichment_serial.err"
#SBATCH --output="logs/cpg_enrichment_serial.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative"
DATA_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/DATA"
RESULTS_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/analyze_mecp2_cpg_enrichment"

# Create exo and endo directories
mkdir -p "${RESULTS_DIR}/exo"
mkdir -p "${RESULTS_DIR}/endo"

# Create logs directory if it doesn't exist
mkdir -p logs
mkdir -p ${RESULTS_DIR}/mecp2_cpg_enrichment

# Copy peaks from results_2_new_005 and organize into exo/endo folders
echo "Copying and organizing peaks from results_2_new_005..."

# Copy exogenous peaks (virus samples)
for sample in NeuV1 NeuV2 NeuV3 NSCv1 NSCv2 NSCv3; do
    cp "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/results_2_new_005/peaks/narrow/${sample}_narrow_peaks.narrowPeak" "${RESULTS_DIR}/exo/"
done

# Copy endogenous peaks (M samples) 
for sample in NeuM2 NeuM3 NSCM1 NSCM2 NSCM3; do
    cp "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/results_2_new_005/peaks/narrow/${sample}_narrow_peaks.narrowPeak" "${RESULTS_DIR}/endo/"
done

# Run the script with working directory argument and full error traceback
python ../scripts/analyze_mecp2_cpg_enrichment_serial.py \
    --exo-peaks-dir ${RESULTS_DIR}/exo \
    --endo-peaks-dir ${RESULTS_DIR}/endo \
    --cpg-islands ${DATA_DIR}/cpg_islands.bed \
    --output-dir ${RESULTS_DIR}/mecp2_cpg_enrichment \
    --bam-dir ${WORKING_DIR}/results_1/aligned \
    2>&1 | tee "logs/mecp2_cpg_enrichment.out"