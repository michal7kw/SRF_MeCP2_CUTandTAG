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
#SBATCH --error="logs/align1_005_mecp2_analysis_%a.err"
#SBATCH --output="logs/align1_005_mecp2_analysis_%a.out"

WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Methylation_dev"
cd $WORKING_DIR || exit 1

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate snakemake

# Function to print usage
print_usage() {
    echo "Usage: sbatch run_mecp2.sh [OPTIONS]"
    echo "Options:"
    echo "  --experiment VALUE     Specify experiment name (e.g., align1_005)"
    echo "  --debug SAMPLE_SIZE    Run in debug mode with specified sample size (default: 100)"
    echo "  --force-recompute      Force recomputation of all analysis stages"
    echo "  --cell-type VALUE      Specify cell type to analyze (NEU, NSC, or ALL) (default: ALL)"
    echo "  --help                 Display this help message"
}

# Parse command line arguments
DEBUG_MODE=false
FORCE_RECOMPUTE=false
SAMPLE_SIZE=100
EXPERIMENT="align1_005"  # Default value
CELL_TYPE_FILTER="NSC"   # Default value

while [[ $# -gt 0 ]]; do
    case $1 in
        --experiment)
            if [[ -n $2 ]]; then
                EXPERIMENT=$2
                shift
            fi
            ;;
        --cell-type)
            if [[ -n $2 ]]; then
                case "${2^^}" in  # Convert to uppercase
                    "NEU"|"NSC"|"ALL")
                        CELL_TYPE_FILTER="${2^^}"
                        shift
                        ;;
                    *)
                        echo "Error: Invalid cell type. Must be NEU, NSC, or ALL"
                        exit 1
                        ;;
                esac
            fi
            ;;
        --debug)
            DEBUG_MODE=true
            if [[ -n $2 && $2 =~ ^[0-9]+$ ]]; then
                SAMPLE_SIZE=$2
                shift
            fi
            ;;
        --force-recompute)
            FORCE_RECOMPUTE=true
            ;;
        --help)
            print_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            print_usage
            exit 1
            ;;
    esac
    shift
done

# Define cell types based on filter
if [ "$CELL_TYPE_FILTER" = "ALL" ]; then
    CELL_TYPES=("NEU" "NSC")
else
    CELL_TYPES=("$CELL_TYPE_FILTER")
fi

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