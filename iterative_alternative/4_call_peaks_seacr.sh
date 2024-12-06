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
#SBATCH --error="logs/peaks_%a.err"
#SBATCH --output="logs/peaks_%a.out"
#SBATCH --array=0-9  # Excluding IgM controls

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Path to SEACR script
SEACR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/SEACR-master/SEACR_1.3.sh"

# Create output directories
mkdir -p results/peaks/seacr
mkdir -p results/bedgraph

# Add this near the top of the script, after the cd command
FORCE_REPROCESS=true  # Set to true to force reprocessing of all files

# First, process IgM control if it hasn't been done yet
if [ ! -f "results/bedgraph/IgM.bedgraph" ] || [ "$FORCE_REPROCESS" = true ]; then
    echo "Processing IgM control..."
    bedtools genomecov -bg -ibam results/aligned/IgM.bam > results/bedgraph/IgM.bedgraph
fi

# Get sample names (excluding IgM)
EXOGENOUS_SAMPLES=($(ls ../DATA/EXOGENOUS/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))
ENDOGENOUS_SAMPLES=($(ls ../DATA/ENDOGENOUS/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//' | grep -v "IgM"))
ALL_SAMPLES=("${EXOGENOUS_SAMPLES[@]}" "${ENDOGENOUS_SAMPLES[@]}")

# Get current sample
SAMPLE=${ALL_SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Convert BAM to bedgraph if not already done
if [ ! -f "results/bedgraph/${SAMPLE}.bedgraph" ] || [ "$FORCE_REPROCESS" = true ]; then
    echo "Processing ${SAMPLE}..."
    bedtools genomecov -bg -ibam results/aligned/${SAMPLE}.bam > results/bedgraph/${SAMPLE}.bedgraph
fi

# Function to convert SEACR output to narrowPeak format
seacr_to_narrowpeak() {
    input=$1
    output=$2
    # Convert SEACR output to narrowPeak format
    # narrowPeak format: chr start end name score strand signalValue pValue qValue peak
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "peak_"NR, $4, ".", $4, -1, -1, -1}' "$input" > "$output"
}

# Run SEACR peak calling with IgM control
# SEACR arguments:
# 1. Target data bedgraph file
# 2. IgM control bedgraph file
# 3. "non" - no normalization between target and control
# 4. "stringent" - more stringent threshold for peak calling
# 5. Output prefix for results
echo "Calling peaks for ${SAMPLE} with IgM control..."
bash $SEACR \
    results/bedgraph/${SAMPLE}.bedgraph \
    results/bedgraph/IgM.bedgraph \
    non \
    stringent \
    results/peaks/seacr/${SAMPLE}

# Convert to narrowPeak
seacr_to_narrowpeak "results/peaks/seacr/${SAMPLE}.stringent.bed" "results/peaks/seacr/${SAMPLE}.stringent.narrowPeak"

# Also run without control (top 1% of peaks)
# SEACR arguments:
# 1. Target data bedgraph file  
# 2. "0.01" - threshold representing top 1% of peaks
# 3. "non" - no normalization (not applicable without control)
# 4. "stringent" - more stringent threshold for peak calling
# 5. Output prefix for results
echo "Calling peaks for ${SAMPLE} without control..."
bash $SEACR \
    results/bedgraph/${SAMPLE}.bedgraph \
    0.01 \
    non \
    stringent \
    results/peaks/seacr/${SAMPLE}_no_control

# Convert to narrowPeak
seacr_to_narrowpeak "results/peaks/seacr/${SAMPLE}_no_control.stringent.bed" "results/peaks/seacr/${SAMPLE}_no_control.stringent.narrowPeak"

# Move .auc files to the desired directory
mv *.auc* results/peaks/seacr/auc/
mv *.txt* results/peaks/seacr/auc/