#!/bin/bash
#SBATCH --job-name=3_mecp2_meth
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/analyze_mecp2_cpg_enrichment_align1_005.err"
#SBATCH --output="logs/analyze_mecp2_cpg_enrichment_align1_005.out"

WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Methylation"
cd $WORKING_DIR || exit 1

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate snakemake

# Function to print usage
print_usage() {
    echo "Usage: sbatch run_mecp2.sh [OPTIONS]"
    echo "Options:"
    echo "  --debug SAMPLE_SIZE    Run in debug mode with specified sample size (default: 100)"
    echo "  --force-recompute      Force recomputation of all analysis stages"
    echo "  --help                 Display this help message"
}

# Parse command line arguments
DEBUG_MODE=false
FORCE_RECOMPUTE=false
SAMPLE_SIZE=100

while [[ $# -gt 0 ]]; do
    case $1 in
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

# Build python command based on options
CMD="python analysis_mecp2_methylation.py"

if [ "$DEBUG_MODE" = true ]; then
    CMD="$CMD --debug --sample-size $SAMPLE_SIZE"
fi

if [ "$FORCE_RECOMPUTE" = true ]; then
    CMD="$CMD --force-recompute"
fi

# Print configuration
echo "Running with configuration:"
echo "Working directory: $WORKING_DIR"
echo "Debug mode: $DEBUG_MODE"
echo "Sample size: $SAMPLE_SIZE"
echo "Force recompute: $FORCE_RECOMPUTE"
echo "Command: $CMD"
echo "-----------------------------------"

# Execute the command
eval $CMD