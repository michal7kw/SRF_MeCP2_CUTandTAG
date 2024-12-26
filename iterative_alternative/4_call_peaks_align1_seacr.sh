#!/bin/bash
#SBATCH --job-name=peaks_seacr
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_seacr.err"
#SBATCH --output="logs/peaks_seacr.out"
#SBATCH --array=0-10  # Excluding IgM controls

# Set base directory and create required directories
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative"
cd $BASE_DIR || exit 1

INPUT_DIR="results_1"
RESULTS_DIR="results_2_align1_seacr"
SEACR_DIR="SEACR-master"

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Create directories
mkdir -p ${RESULTS_DIR}/peaks/seacr/{auc,no_control} ${RESULTS_DIR}/bedgraph ${RESULTS_DIR}/peaks || exit 1

# Configuration
SEACR="$BASE_DIR/$SEACR_DIR/SEACR_1.3.sh"
TIMEOUT=3600  # 1 hour timeout

# Convert SEACR output to narrowPeak format
seacr_to_narrowpeak() {
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "peak_"NR, $4, ".", $4, -1, -1, -1}' "$1" > "$2"
}

# Process IgM control - add error checking and verbose output
echo "Processing IgM control..."
if [ ! -f "${RESULTS_DIR}/bedgraph/IgM.bedgraph" ]; then
    echo "Creating IgM bedgraph at: ${RESULTS_DIR}/bedgraph/IgM.bedgraph"
    
    # Ensure the bedgraph directory exists
    mkdir -p "${RESULTS_DIR}/bedgraph"
    
    # Check if IgM BAM file exists
    if [ ! -f "${INPUT_DIR}/aligned/IgM.bam" ]; then
        echo "Error: IgM BAM file not found at ${INPUT_DIR}/aligned/IgM.bam"
        exit 1
    fi
    
    # Create bedgraph with error checking
    if ! bedtools genomecov -bg -ibam "${INPUT_DIR}/aligned/IgM.bam" > "${RESULTS_DIR}/bedgraph/IgM.bedgraph"; then
        echo "Error: Failed to create IgM bedgraph"
        exit 1
    fi
    
    # Verify the bedgraph was created and is not empty
    if [ ! -s "${RESULTS_DIR}/bedgraph/IgM.bedgraph" ]; then
        echo "Error: IgM bedgraph was created but is empty"
        exit 1
    fi
fi

echo "Verifying IgM bedgraph exists at: ${RESULTS_DIR}/bedgraph/IgM.bedgraph"
ls -l "${RESULTS_DIR}/bedgraph/IgM.bedgraph"

# Get sample names and current sample - fix path and add error checking
echo "Getting sample names..."
DATA_DIR="${BASE_DIR}/../DATA"
if [ ! -d "$DATA_DIR" ]; then
    echo "Error: DATA directory not found at $DATA_DIR"
    exit 1
fi

# Replace the sample name detection with direct array definition
SAMPLES=(
    "NeuM2" "NeuM3" "NeuV1" "NeuV2" "NeuV3"
    "NSCM1" "NSCM2" "NSCM3" "NSCv1" "NSCv2" "NSCv3"
)

# Get current sample
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
if [ -z "$SAMPLE" ]; then
    echo "Error: No sample found for array task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Processing sample: $SAMPLE"

# Define absolute paths
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative"
INPUT_DIR="${BASE_DIR}/results_1"
RESULTS_DIR="${BASE_DIR}/results_2_align1_seacr"
SEACR_DIR="${BASE_DIR}/SEACR-master"

# Define input paths using absolute paths
SAMPLE_BG=$(readlink -f "${INPUT_DIR}/bedgraph/${SAMPLE}.bedgraph")
IgM_BG=$(readlink -f "${INPUT_DIR}/bedgraph/IgM.bedgraph")

# Debug info
echo "Debug information:"
echo "Sample bedgraph: $SAMPLE_BG"
echo "IgM bedgraph: $IgM_BG"
echo "File exists check:"
ls -l "$SAMPLE_BG"
ls -l "$IgM_BG"

# Verify files exist before running SEACR
if [ ! -f "$SAMPLE_BG" ] || [ ! -f "$IgM_BG" ]; then
    echo "Error: Required bedgraph files not found:"
    echo "Sample bedgraph: $SAMPLE_BG"
    echo "IgM bedgraph: $IgM_BG"
    exit 1
fi

# Run SEACR with control
echo "Running SEACR with IgM control for ${SAMPLE}..."
OUTPUT_PREFIX=$(readlink -f "${RESULTS_DIR}/peaks/seacr/${SAMPLE}")

if ! bash "$SEACR" \
    "$SAMPLE_BG" \
    "$IgM_BG" \
    non \
    stringent \
    "$OUTPUT_PREFIX"; then
    echo "Error: SEACR failed for ${SAMPLE} with IgM control"
    exit 1
fi

echo "Successfully completed SEACR analysis for ${SAMPLE}"

# Run SEACR without control (no-control mode)
echo "Running SEACR without control for ${SAMPLE}..."
OUTPUT_PREFIX_NC=$(readlink -f "${RESULTS_DIR}/peaks/seacr/no_control/${SAMPLE}")

# Create new temporary directory for no-control run
TEMP_DIR_NC="${RESULTS_DIR}/tmp/${SAMPLE}_nc"
mkdir -p "$TEMP_DIR_NC"
cd "$TEMP_DIR_NC" || exit 1

if ! bash "$SEACR" \
    "$SAMPLE_BG" \
    "0.01" \
    non \
    stringent \
    "$OUTPUT_PREFIX_NC"; then
    echo "Error: SEACR failed for ${SAMPLE} without control"
    exit 1
fi

# Clean up temporary directories
rm -rf "$TEMP_DIR" "$TEMP_DIR_NC"

# Convert outputs to narrowPeak format
for mode in "" "_no_control"; do
    input_bed="${RESULTS_DIR}/peaks/seacr/${SAMPLE}${mode}.stringent.bed"
    output_peak="${RESULTS_DIR}/peaks/seacr/${SAMPLE}${mode}.stringent.narrowPeak"
    
    if [ -f "$input_bed" ]; then
        seacr_to_narrowpeak "$input_bed" "$output_peak"
    else
        echo "Warning: Expected output file not found: $input_bed"
    fi
done

echo "All SEACR analyses completed for ${SAMPLE}"

# Move AUC files to final location
mv *.auc.bed "${RESULTS_DIR}/peaks/seacr/auc/" 2>/dev/null || true