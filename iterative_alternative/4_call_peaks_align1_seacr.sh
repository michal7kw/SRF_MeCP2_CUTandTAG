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

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Create directories
mkdir -p ${RESULTS_DIR}/peaks/seacr/{auc,no_control} ${RESULTS_DIR}/bedgraph ${RESULTS_DIR}/peaks || exit 1

# Configuration
SEACR="$BASE_DIR/SEACR-master/SEACR_1.3.sh"
TIMEOUT=3600  # 1 hour timeout

# Convert SEACR output to narrowPeak format
seacr_to_narrowpeak() {
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "peak_"NR, $4, ".", $4, -1, -1, -1}' "$1" > "$2"
}

# Process IgM control
if [ ! -f "${INPUT_DIR}/bedgraph/IgM.bedgraph" ]; then
    echo "Processing IgM control..."
    bedtools genomecov -bg -ibam "${INPUT_DIR}/aligned/IgM.bam" > "${INPUT_DIR}/bedgraph/IgM.bedgraph"
fi

# Get sample names and current sample
SAMPLES=($(ls ../DATA/{EXOGENOUS,ENDOGENOUS}/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//' | grep -v "IgM"))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Create bedgraph if needed
if [ ! -f "${INPUT_DIR}/bedgraph/${SAMPLE}.bedgraph" ]; then
    echo "Creating bedgraph for ${SAMPLE}..."
    bedtools genomecov -bg -ibam "${INPUT_DIR}/aligned/${SAMPLE}.bam" > "${INPUT_DIR}/bedgraph/${SAMPLE}.bedgraph"
fi

# Create temporary directory
TEMP_DIR=$(mktemp -d -p ${RESULTS_DIR}/peaks/seacr/auc)
trap 'rm -rf "$TEMP_DIR"' EXIT

# Copy SEACR files to temp directory
cp -r SEACR-master/* "$TEMP_DIR/"
cd "$TEMP_DIR" || exit 1

# Run SEACR with IgM control
echo "Running SEACR with IgM control for ${SAMPLE}..."
timeout $TIMEOUT ./SEACR_1.3.sh \
    "${INPUT_DIR}/bedgraph/${SAMPLE}.bedgraph" \
    "${INPUT_DIR}/bedgraph/IgM.bedgraph" \
    non \
    stringent \
    "${RESULTS_DIR}/peaks/seacr/${SAMPLE}"

if [ $? -eq 0 ]; then
    # Convert to narrowPeak if successful
    seacr_to_narrowpeak \
        "${RESULTS_DIR}/peaks/seacr/${SAMPLE}.stringent.bed" \
        "${RESULTS_DIR}/peaks/seacr/${SAMPLE}.stringent.narrowPeak"
    
    # Copy to final location
    cp "${RESULTS_DIR}/peaks/seacr/${SAMPLE}.stringent.narrowPeak" \
       "${RESULTS_DIR}/peaks/${SAMPLE}_peaks.narrowPeak"
else
    echo "Error: SEACR failed for ${SAMPLE} with IgM control"
    exit 1
fi

# Run SEACR without control
echo "Running SEACR without control for ${SAMPLE}..."
timeout $TIMEOUT ./SEACR_1.3.sh \
    "${INPUT_DIR}/bedgraph/${SAMPLE}.bedgraph" \
    0.05 \
    non \
    stringent \
    "${RESULTS_DIR}/peaks/seacr/${SAMPLE}_no_control"

if [ $? -eq 0 ]; then
    # Convert to narrowPeak if successful
    seacr_to_narrowpeak \
        "${RESULTS_DIR}/peaks/seacr/${SAMPLE}_no_control.stringent.bed" \
        "${RESULTS_DIR}/peaks/seacr/${SAMPLE}_no_control.stringent.narrowPeak"
    
    # Move to no_control directory
    mv "${RESULTS_DIR}/peaks/seacr/${SAMPLE}_no_control.stringent."* \
       "${RESULTS_DIR}/peaks/seacr/no_control/"
else
    echo "Error: SEACR failed for ${SAMPLE} without control"
    exit 1
fi

# Move AUC files to final location
mv *.auc.bed "${RESULTS_DIR}/peaks/seacr/auc/" 2>/dev/null || true