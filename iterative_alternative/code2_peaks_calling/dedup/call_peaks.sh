#!/bin/bash
#SBATCH --job-name=pc_dedup
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peaks_calling/dedup/peaks_calling.err"
#SBATCH --output="logs/peaks_calling/dedup/peaks_calling.out"
#SBATCH --array=0-10  # Excluding IgM controls

# Define paths
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"

# Define sample names
SAMPLES=(
    "NeuM2" "NeuM3" "NeuV1" "NeuV2" "NeuV3"
    "NSCM1" "NSCM2" "NSCM3" "NSCv1" "NSCv2" "NSCv3"
)
# Get current sample
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

INPUT_DIR="results_1b"
RESULTS_DIR="results/dedup"
QC_DIR="${RESULTS_DIR}/qc/${SAMPLE}"

INPUT_BAM="${INPUT_DIR}/aligned/${SAMPLE}.dedup.bam"
CONTROL_BAM="${INPUT_DIR}/aligned/IgM.dedup.bam"
BLACKLIST="${BASE_DIR}/DATA/mm10-blacklist.bed"


# Set working directory
WORK_DIR="${BASE_DIR}/iterative_alternative"
cd $WORK_DIR || exit 1

# Create required directories
mkdir -p ${RESULTS_DIR}/{peaks/{broad,narrow},qc,logs,reports}
mkdir -p "$QC_DIR"

# Load required modules and environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# # Source the module system
# if [ -f "/etc/profile.d/modules.sh" ]; then
#     source /etc/profile.d/modules.sh
# elif [ -f "/usr/share/Modules/init/bash" ]; then
#     source /usr/share/Modules/init/bash
# fi

# # Check if MACS2 is available
# if ! command -v macs2 &> /dev/null; then
#     echo "Error: macs2 command not found. Loading module..."
#     module load macs2/2.2.7.1
#     if ! command -v macs2 &> /dev/null; then
#         echo "Error: Failed to load macs2. Exiting."
#         exit 1
#     fi
# fi

# # Check if samtools is available
# if ! command -v samtools &> /dev/null; then
#     echo "Error: samtools not found. Loading module..."
#     module load samtools
#     if ! command -v samtools &> /dev/null; then
#         echo "Error: Failed to load samtools. Control quality checks will be skipped."
#         SKIP_CONTROL_QC=true
#     fi
# fi

# # Print environment info
# echo "Using MACS2 from: $(which macs2)"
# echo "MACS2 version: $(macs2 --version)"

# Validation functions
validate_inputs() {
    local files=("$INPUT_BAM" "$CONTROL_BAM" "$BLACKLIST")
    for file in "${files[@]}"; do
        if [ ! -f "$file" ]; then
            echo "Error: Required file $file not found"
            exit 1
        fi
    done
}

check_bam_integrity() {
    local bam=$1
    if ! samtools quickcheck "$bam"; then
        echo "Error: BAM file $bam is corrupted"
        exit 1
    fi
}

# QC functions
generate_fragment_dist() {
    local bam=$1
    local output="${QC_DIR}/fragment_dist.txt"
    
    # Use proper fragment length calculation for paired-end reads
    samtools view "$bam" | \
        awk '$9 > 0 {print $9}' > "$output"  # Only take positive TLEN values
    
    # Check if the output file has data
    if [ -s "$output" ]; then
        Rscript -e "
            data <- read.table('$output')\$V1
            png('${QC_DIR}/fragment_dist.png')
            hist(data, 
                breaks=100, 
                main='Fragment Size Distribution',
                xlab='Fragment Length (bp)',
                ylab='Frequency')
            dev.off()
        "
    else
        echo "Warning: No fragment length information found in BAM file"
    fi
}

calculate_frip() {
    local bam=$1
    local peaks=$2
    local total_reads=$(samtools view -c "$bam")
    local reads_in_peaks=$(bedtools intersect -a "$bam" -b "$peaks" -c | awk '{sum+=$NF} END {print sum}')
    echo "scale=4; $reads_in_peaks / $total_reads" | bc > "${QC_DIR}/frip_score.txt"
}

# Peak calling function with optimized parameters for Cut&Tag
call_peaks() {
    local peak_type=$1
    local outdir="${RESULTS_DIR}/peaks/$peak_type"
    
    # Common MACS2 parameters optimized for Cut&Tag
    local common_params=(
        -t "$INPUT_BAM"
        -c "$CONTROL_BAM"
        -n "${SAMPLE}_${peak_type}"
        --outdir "$outdir"
        -g mm
        -f BAMPE
        --nomodel
        --nolambda
        -q 0.05
        --keep-dup all
        --extsize 200
    )

    if [ "$peak_type" = "broad" ]; then
        macs2 callpeak "${common_params[@]}" \
            --broad \
            --broad-cutoff 0.05
    else
        macs2 callpeak "${common_params[@]}"
    fi

    # Filter blacklist regions
    bedtools intersect \
        -v -a "${outdir}/${SAMPLE}_${peak_type}_peaks.${peak_type}Peak" \
        -b "$BLACKLIST" \
        > "${outdir}/${SAMPLE}_${peak_type}_peaks.filtered.${peak_type}Peak"
}

# TSS enrichment analysis
calculate_tss_enrichment() {
    computeMatrix reference-point \
        -S "${QC_DIR}/${SAMPLE}_signal.bw" \
        -R ./DATA/mm10_TSS.bed \
        -o "${QC_DIR}/tss_matrix.gz" \
        --referencePoint TSS \
        -b 5000 -a 5000

    plotHeatmap \
        -m "${QC_DIR}/tss_matrix.gz" \
        -o "${QC_DIR}/tss_enrichment.png" \
        --colorMap RdYlBu \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 --zMax 10 \
        --plotWidth 12 \
        --plotHeight 8
}

# Main execution
main() {
    # Validate inputs
    validate_inputs
    check_bam_integrity "$INPUT_BAM"
    check_bam_integrity "$CONTROL_BAM"

    # Generate QC metrics
    generate_fragment_dist "$INPUT_BAM"

    # Call peaks
    call_peaks "broad"
    call_peaks "narrow"

    # Post-processing QC
    calculate_frip "$INPUT_BAM" "${RESULTS_DIR}/peaks/narrow/${SAMPLE}_narrow_peaks.filtered.narrowPeak"
    calculate_tss_enrichment

    # Generate signal tracks
    bamCoverage -b "$INPUT_BAM" \
        -o "${QC_DIR}/${SAMPLE}_signal.bw" \
        --binSize 10 \
        --normalizeUsing RPKM \
        --extendReads

    # Generate QC report
    Rscript ./generate_qc_report.R \
        --sample "$SAMPLE" \
        --qc-dir "$QC_DIR" \
        --output "${RESULTS_DIR}/reports/${SAMPLE}_qc_report.html"
}

# Execute main function with error handling
{
    main
} 2>&1 | tee "${RESULTS_DIR}/logs/${SAMPLE}_processing.log"

echo "Analysis completed for sample $SAMPLE"


#####################################################################
# This script performs peak calling and QC analysis on CUT&Tag data
#
# Input files:
# - BAM files in ${INPUT_DIR}/aligned/
#   - Sample BAM files (NeuM2,3, NeuV1-3, NSCM1-3, NSCv1-3)
#   - Control BAM file (IgM.bam)
# - Blacklist regions (mm10-blacklist.bed)
#
# Output files in ${RESULTS_DIR}:
# - peaks/narrow/: MACS2 narrow peak calls (filtered and unfiltered)
# - qc/${SAMPLE}/: QC metrics and plots
#   - Fragment size distribution
#   - TSS enrichment plots
#   - Signal tracks (.bw)
#   - FRiP scores
# - reports/: HTML QC reports per sample
# - logs/: Processing logs
#
# The script:
# 1. Validates input files
# 2. Calls peaks using MACS2 in narrow mode
# 3. Generates QC metrics (fragment dist, TSS enrichment, FRiP)
# 4. Creates normalized signal tracks
# 5. Compiles QC report
#####################################################################
