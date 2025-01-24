#!/bin/bash
#SBATCH --job-name=mecp2_analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --array=0-41  # 2 cell types * 21 chromosomes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/mecp2_analysis_%a.err"
#SBATCH --output="logs/mecp2_analysis_%a.out"

WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Methylation_dev"
cd $WORKING_DIR || exit 1

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate snakemake

# Parse command line arguments
DEBUG_MODE=false
FORCE_RECOMPUTE=false
SAMPLE_SIZE=100
EXPERIMENT="align2_005"  # Default value
CELL_TYPE_FILTER="NSC"   # Default value

CHROMOSOMES=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" \
            "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" \
            "chrX" "chrY")

# Adjust array task count based on cell type filter
if [ "$CELL_TYPE_FILTER" != "ALL" ]; then
    #SBATCH --array=0-20  # 1 cell type * 21 chromosomes
    SLURM_ARRAY_TASK_COUNT=21
else
    #SBATCH --array=0-41  # 2 cell types * 21 chromosomes
    SLURM_ARRAY_TASK_COUNT=42
fi

# Calculate cell type and chromosome from array task ID
CELL_TYPE_IDX=$((SLURM_ARRAY_TASK_ID / ${#CHROMOSOMES[@]}))
CHR_IDX=$((SLURM_ARRAY_TASK_ID % ${#CHROMOSOMES[@]}))

CELL_TYPE=${CELL_TYPES[$CELL_TYPE_IDX]}
CHROMOSOME=${CHROMOSOMES[$CHR_IDX]}

echo "Processing $CELL_TYPE - $CHROMOSOME"

# Build python command
CMD="python analysis_mecp2_methylation.py \
     --experiment $EXPERIMENT \
     --processes $SLURM_NTASKS \
     --cell-type $CELL_TYPE \
     --chromosome $CHROMOSOME"

if [ "$DEBUG_MODE" = true ]; then
    CMD="$CMD --debug --sample-size $SAMPLE_SIZE"
fi

if [ "$FORCE_RECOMPUTE" = true ]; then
    CMD="$CMD --force-recompute"
fi

# Print configuration
echo "Running with configuration:"
echo "Working directory: $WORKING_DIR"
echo "Experiment: $EXPERIMENT"
echo "Cell type filter: $CELL_TYPE_FILTER"
echo "Cell type: $CELL_TYPE"
echo "Chromosome: $CHROMOSOME"
echo "Debug mode: $DEBUG_MODE"
echo "Sample size: $SAMPLE_SIZE"
echo "Force recompute: $FORCE_RECOMPUTE"
echo "Processes per job: $SLURM_NTASKS"
echo "Command: $CMD"
echo "-----------------------------------"

# Execute the command
eval $CMD

# After completion, check if this is the last job in the array
if [ $SLURM_ARRAY_TASK_ID -eq $((SLURM_ARRAY_TASK_COUNT - 1)) ]; then
    echo "All jobs completed, running final merge..."
    python merge_results.py
fi