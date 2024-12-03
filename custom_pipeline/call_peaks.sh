#!/bin/bash
#SBATCH --job-name=call_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/peaks_%A_%a.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/peaks_%A_%a.out"
#SBATCH --array=0-11%4  # Adjust numbers based on your sample count

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Source the module system
if [ -f "/etc/profile.d/modules.sh" ]; then
    source /etc/profile.d/modules.sh
elif [ -f "/usr/share/Modules/init/bash" ]; then
    source /usr/share/Modules/init/bash
fi

# Check if MACS2 is available
if ! command -v macs2 &> /dev/null; then
    echo "Error: macs2 command not found. Loading module..."
    module load macs2/2.2.7.1
    if ! command -v macs2 &> /dev/null; then
        echo "Error: Failed to load macs2. Exiting."
        exit 1
    fi
fi

# Check if samtools is available
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found. Loading module..."
    module load samtools
    if ! command -v samtools &> /dev/null; then
        echo "Error: Failed to load samtools. Control quality checks will be skipped."
        SKIP_CONTROL_QC=true
    fi
fi

# Print environment info
echo "Using MACS2 from: $(which macs2)"
echo "MACS2 version: $(macs2 --version)"

# Set default paths
BAM_DIR="results/aligned"
OUTPUT_DIR="results/peaks"
CONTROL_BAM="results/aligned/IgM.bam"

# Verify input files exist
if [ ! -d "$BAM_DIR" ]; then
    echo "Error: BAM directory not found: $BAM_DIR"
    exit 1
fi

if [ ! -f "$CONTROL_BAM" ]; then
    echo "Error: Control BAM file not found: $CONTROL_BAM"
    exit 1
fi

# Default values
GENOME_SIZE="hs"
QVALUE="0.01"  # More stringent q-value for higher confidence peaks
FORMAT="BAMPE"
THREADS=16
BROAD="--broad"  # Use broad peak calling for MeCP2

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Get list of BAM files (excluding control)
mapfile -t BAM_FILES < <(find "$BAM_DIR" -name "*.bam" ! -name "$(basename $CONTROL_BAM)")

# Get the current sample based on SLURM array task ID
if [ ${#BAM_FILES[@]} -le $SLURM_ARRAY_TASK_ID ]; then
    echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) exceeds number of samples (${#BAM_FILES[@]})"
    exit 1
fi

CURRENT_BAM="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"
SAMPLE=$(basename "$CURRENT_BAM" .bam)

echo "Processing sample $SAMPLE (array task $SLURM_ARRAY_TASK_ID)"

# Prepare control parameter
CONTROL_PARAM=""
if [ -n "$CONTROL_BAM" ]; then
    CONTROL_PARAM="-c $CONTROL_BAM"
fi

# Run MACS2
echo "Running MACS2 for $SAMPLE..."

# Set peak calling parameters based on broad flag
SUMMIT_PARAM=""
if [ -z "$BROAD" ]; then
    SUMMIT_PARAM="--call-summits"
fi

echo "Command: macs2 callpeak \\"
echo "    -t $CURRENT_BAM \\"
echo "    $CONTROL_PARAM \\"
echo "    -f $FORMAT \\"
echo "    -g $GENOME_SIZE \\"
echo "    -n $SAMPLE \\"
echo "    --outdir $OUTPUT_DIR \\"
echo "    -q $QVALUE \\"
echo "    --nomodel \\"
echo "    $SUMMIT_PARAM \\"
echo "    $BROAD"

macs2 callpeak \
    -t "$CURRENT_BAM" \
    $CONTROL_PARAM \
    -f "$FORMAT" \
    -g "$GENOME_SIZE" \
    -n "$SAMPLE" \
    --outdir "$OUTPUT_DIR" \
    -q "$QVALUE" \
    --nomodel \
    --broad \
    --broad-cutoff 0.05 \
    --keep-dup 1 \
    --buffer-size 10000 \
    2>&1 | tee "$OUTPUT_DIR/${SAMPLE}_macs2.log"

if [ $? -ne 0 ]; then
    echo "Error: MACS2 failed for sample $SAMPLE"
    exit 1
fi

echo "Finished processing $SAMPLE"

# Verify output was created
if [ ! -f "$OUTPUT_DIR/${SAMPLE}_peaks.broadPeak" ]; then
    echo "Error: Output file not created for sample $SAMPLE"
    exit 1
fi

echo "Peak calling completed for sample $SAMPLE"

# Function to check control sample quality
check_control() {
    if [ "$SKIP_CONTROL_QC" = true ]; then
        echo "WARNING: Skipping control quality checks (samtools not available)"
        return 0
    fi

    local CONTROL_BAM="$1"
    local MIN_MAPQ=30
    
    echo "Checking control sample quality..."
    
    # Check mapping quality
    local MAPPED_READS=$(samtools view -F 0x4 -q ${MIN_MAPQ} -c "$CONTROL_BAM" || echo "ERROR")
    if [ "$MAPPED_READS" = "ERROR" ]; then
        echo "WARNING: Failed to check control sample quality"
        return 1
    fi
    
    local TOTAL_READS=$(samtools view -c "$CONTROL_BAM")
    
    # Only calculate mapping rate if we got valid numbers
    if [[ "$MAPPED_READS" =~ ^[0-9]+$ ]] && [[ "$TOTAL_READS" =~ ^[0-9]+$ ]]; then
        local MAPPING_RATE=$(awk "BEGIN {print $MAPPED_READS/$TOTAL_READS}")
        
        echo "Control sample QC:"
        echo "Total reads: $TOTAL_READS"
        echo "Mapped reads (MAPQ>=$MIN_MAPQ): $MAPPED_READS"
        echo "Mapping rate: $MAPPING_RATE"
        
        # Add warning if mapping rate is low
        if (( $(echo "$MAPPING_RATE < 0.5" | bc -l) )); then
            echo "WARNING: Low mapping rate in control sample"
        fi
    else
        echo "WARNING: Failed to calculate mapping statistics"
    fi

    if [ "$TOTAL_READS" -lt 10000000 ]; then
        echo "WARNING: Low read count in control sample"
    fi
}