#!/bin/bash
#SBATCH --job-name=peaks_original
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_original.err"
#SBATCH --output="logs/peaks_original.out"
#SBATCH --array=0-10  # Excluding IgM controls

# Set working directory
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative"
cd $BASE_DIR || exit 1

# Set results directory
RESULTS_DIR="results_2_original"

# Create required directories
mkdir -p ${RESULTS_DIR}/peaks
mkdir -p ${RESULTS_DIR}/peaks/broad

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

# Define sample names
SAMPLES=(
    "NeuM2"
    "NeuM3"
    "NeuV1"
    "NeuV2"
    "NeuV3"
    "NSCM1"
    "NSCM2"
    "NSCM3"
    "NSCv1"
    "NSCv2"
    "NSCv3"
)

# Get the current sample
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Define input and control files
INPUT_BAM="results_1/aligned/${SAMPLE}.bam"
CONTROL_BAM="results_1/aligned/IgM.bam"

# Run MACS2 for broad peaks
echo "Running MACS2 broad peak calling for ${SAMPLE}..."
macs2 callpeak \
    -t "$INPUT_BAM" \
    -c "$CONTROL_BAM" \
    -n "${SAMPLE}_broad" \
    --outdir "${RESULTS_DIR}/peaks/broad" \
    -g mm \
    -f BAMPE \
    --broad \
    --broad-cutoff 0.1 \
    -q 0.05 \
    --keep-dup all

# Run MACS2 for narrow peaks
echo "Running MACS2 narrow peak calling for ${SAMPLE}..."
macs2 callpeak \
    -t "$INPUT_BAM" \
    -c "$CONTROL_BAM" \
    -n "${SAMPLE}_narrow" \
    --outdir "${RESULTS_DIR}/peaks" \
    -g mm \
    -f BAMPE \
    -q 0.05 \
    --keep-dup all

if [ $? -ne 0 ]; then
    echo "Error: MACS2 failed for sample $SAMPLE"
    exit 1
fi

echo "Finished processing $SAMPLE"

# Verify output was created
if [ ! -f "${RESULTS_DIR}/peaks/broad/${SAMPLE}_broad_peaks.broadPeak" ]; then
    echo "Error: Output file not created for sample $SAMPLE (broad)"
    exit 1
fi

if [ ! -f "${RESULTS_DIR}/peaks/${SAMPLE}_narrow_peaks.narrowPeak" ]; then
    echo "Error: Output file not created for sample $SAMPLE (narrow)"
    exit 1
fi

echo "Peak calling completed for sample $SAMPLE"