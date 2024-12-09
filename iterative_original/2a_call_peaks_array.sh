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
#SBATCH --error="logs/peaks.err"
#SBATCH --output="logs/peaks.out"
#SBATCH --array=0-9

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_original

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
INPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/snake_pipeline/results"
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_original/results"

BAM_DIR="${INPUT_DIR}/aligned"
CONTROL_BAM="${INPUT_DIR}/aligned/IgM.bam"
OUTPUT_DIR="${OUTPUT_DIR}/peaks"

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
GENOME_SIZE="mm"
QVALUE="0.05"
FORMAT="BAMPE"

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

# First run: Narrow peaks
echo "Calling narrow peaks..."
echo "Command: macs2 callpeak \\"
echo "    -t $CURRENT_BAM \\"
echo "    $CONTROL_PARAM \\"
echo "    -f $FORMAT \\"
echo "    -g $GENOME_SIZE \\"
echo "    -n ${SAMPLE} \\"
echo "    --outdir $OUTPUT_DIR \\"
echo "    -q $QVALUE \\"
echo "    --nomodel \\"
echo "    --keep-dup auto \\"
echo "    --nolambda \\"
echo "    --bdg"

macs2 callpeak \
    -t "$CURRENT_BAM" \
    $CONTROL_PARAM \
    -f "$FORMAT" \
    -g "$GENOME_SIZE" \
    -n "${SAMPLE}" \
    --outdir "$OUTPUT_DIR" \
    -q "$QVALUE" \
    --nomodel \
    --keep-dup auto \
    --nolambda \
    --bdg \
    2>&1 | tee "$OUTPUT_DIR/${SAMPLE}_macs2_narrow.log"

if [ $? -ne 0 ]; then
    echo "Error: MACS2 narrow peak calling failed for sample $SAMPLE"
    exit 1
fi

# Second run: Broad peaks
echo "Calling broad peaks..."
echo "Command: macs2 callpeak \\"
echo "    -t $CURRENT_BAM \\"
echo "    $CONTROL_PARAM \\"
echo "    -f $FORMAT \\"
echo "    -g $GENOME_SIZE \\"
echo "    -n ${SAMPLE}_broad \\"
echo "    --outdir $OUTPUT_DIR \\"
echo "    -q $QVALUE \\"
echo "    --nomodel \\"
echo "    --broad \\"
echo "    --broad-cutoff 0.05 \\"
echo "    --keep-dup auto \\"
echo "    --nolambda \\"
echo "    --bdg"

macs2 callpeak \
    -t "$CURRENT_BAM" \
    $CONTROL_PARAM \
    -f "$FORMAT" \
    -g "$GENOME_SIZE" \
    -n "${SAMPLE}_broad" \
    --outdir "$OUTPUT_DIR" \
    -q "$QVALUE" \
    --nomodel \
    --broad \
    --broad-cutoff 0.05 \
    --keep-dup auto \
    --nolambda \
    --bdg \
    2>&1 | tee "$OUTPUT_DIR/${SAMPLE}_macs2_broad.log"

if [ $? -ne 0 ]; then
    echo "Error: MACS2 broad peak calling failed for sample $SAMPLE"
    exit 1
fi

echo "Finished processing $SAMPLE"

# Verify outputs were created
if [ ! -f "$OUTPUT_DIR/${SAMPLE}_peaks.narrowPeak" ]; then
    echo "Error: Narrow peak file not created for sample $SAMPLE"
    exit 1
fi

if [ ! -f "$OUTPUT_DIR/${SAMPLE}_broad_peaks.broadPeak" ]; then
    echo "Error: Broad peak file not created for sample $SAMPLE"
    exit 1
fi

# Rename broad peak files to match original naming
mv "$OUTPUT_DIR/${SAMPLE}_broad_peaks.broadPeak" "$OUTPUT_DIR/${SAMPLE}_peaks.broadPeak"
mv "$OUTPUT_DIR/${SAMPLE}_broad_peaks.gappedPeak" "$OUTPUT_DIR/${SAMPLE}_peaks.gappedPeak"
mv "$OUTPUT_DIR/${SAMPLE}_broad_peaks.xls" "$OUTPUT_DIR/${SAMPLE}_peaks.xls"
mv "$OUTPUT_DIR/${SAMPLE}_broad_treat_pileup.bdg" "$OUTPUT_DIR/${SAMPLE}_treat_pileup.bdg"
mv "$OUTPUT_DIR/${SAMPLE}_broad_control_lambda.bdg" "$OUTPUT_DIR/${SAMPLE}_control_lambda.bdg"

echo "Peak calling completed for sample $SAMPLE"