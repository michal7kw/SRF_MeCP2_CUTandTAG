#!/bin/bash
#SBATCH --job-name=peak_compare
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/peak_compare.err"
#SBATCH --output="logs/peak_compare.out"

# Set working directory
BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative"
cd $BASE_DIR || exit 1

source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

RESULTS_DIR="results_2_align2_005"

# Create analysis directories
mkdir -p ${RESULTS_DIR}/peak_analysis/{differential,overlaps,cpg_analysis,plots}

# Define sample groups
WILD_TYPE_NEU=("NeuV2" "NeuV3")
MUTANT_NEU=("NeuM2" "NeuM3")
WILD_TYPE_NSC=("NSCv1" "NSCv2" "NSCv3")
MUTANT_NSC=("NSCM1" "NSCM2" "NSCM3")

# Path to CpG islands
CpG_ISLANDS="${BASE_DIR}/DATA/cpg_islands.bed"

# Function to merge peaks from replicates
merge_peaks() {
    local group_name=$1
    shift
    local samples=("$@")
    local output="${RESULTS_DIR}/peak_analysis/${group_name}_merged.bed"
    
    # Concatenate peaks from all replicates
    for sample in "${samples[@]}"; do
        cat "${RESULTS_DIR}/peaks/narrow/${sample}_narrow_peaks.filtered.narrowPeak"
    done | \
    # Sort and merge overlapping peaks
    sort -k1,1 -k2,2n | \
    bedtools merge -c 7 -o max > "$output"
}

# Function to filter peaks overlapping CpG islands
filter_cpg_peaks() {
    local peaks=$1
    local output=$2
    
    # Find peaks overlapping CpG islands and keep original peak scores
    bedtools intersect -a "$peaks" -b "$CpG_ISLANDS" -wa \
        | sort -k1,1 -k2,2n \
        | uniq > "$output"
}

# Function to calculate peak statistics
calculate_peak_stats() {
    local group_name=$1
    local peaks=$2
    local output_prefix="${RESULTS_DIR}/peak_analysis/stats/${group_name}"
    
    # Calculate peak widths
    awk '{print $3-$2}' "$peaks" > "${output_prefix}_widths.txt"
    
    # Calculate basic statistics
    Rscript - <<EOF
    widths <- scan("${output_prefix}_widths.txt")
    stats <- c(mean=mean(widths), median=median(widths), 
               sd=sd(widths), total=sum(widths))
    write.table(t(stats), "${output_prefix}_stats.txt", 
                quote=FALSE, sep="\t")
EOF
}

# Function to normalize peak signals
normalize_peaks() {
    local input_bam=$1
    local peaks=$2
    local output=$3
    
    # Calculate RPKM normalized coverage
    bedtools coverage -a "$peaks" -b "$input_bam" \
        -g mm10.genome -mean -pc | \
    awk -v OFS="\t" '{
        # Normalize by peak width and library size
        width = $3-$2
        print $1, $2, $3, $4 * 1000 / width
    }' > "$output"
}

# Function to perform differential peak analysis
analyze_differential_peaks() {
    local wt_peaks=$1
    local mut_peaks=$2
    local output_prefix=$3
    local cpg_only=${4:-false}  # New parameter for CpG-only analysis
    
    local wt_input="$wt_peaks"
    local mut_input="$mut_peaks"
    
    if [ "$cpg_only" = true ]; then
        # Filter for CpG-overlapping peaks
        wt_input="${output_prefix}_wt_cpg.bed"
        mut_input="${output_prefix}_mut_cpg.bed"
        filter_cpg_peaks "$wt_peaks" "$wt_input"
        filter_cpg_peaks "$mut_peaks" "$mut_input"
    fi
    
    # Find common and unique peaks
    bedtools intersect -a "$wt_input" -b "$mut_input" -wo > "${output_prefix}_common.bed"
    bedtools intersect -a "$wt_input" -b "$mut_input" -v > "${output_prefix}_wt_specific.bed"
    bedtools intersect -a "$mut_input" -b "$wt_input" -v > "${output_prefix}_mut_specific.bed"
    
    # Analyze peak differences
    Rscript ../scripts/analyze_differential_peaks.R \
        --common "${RESULTS_DIR}/peak_analysis/${output_prefix}_common.bed" \
        --wt-specific "${RESULTS_DIR}/peak_analysis/${output_prefix}_wt_specific.bed" \
        --mut-specific "${RESULTS_DIR}/peak_analysis/${output_prefix}_mut_specific.bed" \
        --output "${RESULTS_DIR}/peak_analysis/${output_prefix}_analysis.pdf" \
        --cpg-only "$cpg_only"
}

# Function to analyze CpG island overlaps
analyze_cpg_overlaps() {
    local peaks=$1
    local group_name=$2
    local output_prefix="${RESULTS_DIR}/peak_analysis/cpg_analysis/${group_name}"
    
    # Find peaks overlapping CpG islands
    bedtools intersect -a "$peaks" -b "$CpG_ISLANDS" -wo > "${output_prefix}_cpg_overlaps.bed"
    
    # Calculate statistics
    Rscript ../scripts/analyze_cpg_overlaps.R \
        --overlaps "${output_prefix}_cpg_overlaps.bed" \
        --peaks "$peaks" \
        --output "${output_prefix}_cpg_analysis.pdf"
}

# Main execution
main() {
    # Merge peaks for each group
    merge_peaks "wt_neu" "${WILD_TYPE_NEU[@]}"
    merge_peaks "mut_neu" "${MUTANT_NEU[@]}"
    merge_peaks "wt_nsc" "${WILD_TYPE_NSC[@]}"
    merge_peaks "mut_nsc" "${MUTANT_NSC[@]}"
    
    # Calculate statistics for each group
    for group in wt_neu mut_neu wt_nsc mut_nsc; do
        calculate_peak_stats "$group" "${RESULTS_DIR}/peak_analysis/${group}_merged.bed"
    done
    
    # Perform differential analysis on all peaks
    analyze_differential_peaks \
        "${RESULTS_DIR}/peak_analysis/wt_neu_merged.bed" \
        "${RESULTS_DIR}/peak_analysis/mut_neu_merged.bed" \
        "${RESULTS_DIR}/peak_analysis/differential/neu" \
        false

    # Perform differential analysis on CpG-only peaks
    analyze_differential_peaks \
        "${RESULTS_DIR}/peak_analysis/wt_neu_merged.bed" \
        "${RESULTS_DIR}/peak_analysis/mut_neu_merged.bed" \
        "${RESULTS_DIR}/peak_analysis/differential/neu_cpg" \
        true
    
    # Perform differential analysis on NSC comparisons
    analyze_differential_peaks \
        "${RESULTS_DIR}/peak_analysis/wt_nsc_merged.bed" \
        "${RESULTS_DIR}/peak_analysis/mut_nsc_merged.bed" \
        "${RESULTS_DIR}/peak_analysis/differential/nsc" \
        false
    
    # Analyze CpG island overlaps
    for group in wt_neu mut_neu wt_nsc mut_nsc; do
        analyze_cpg_overlaps \
            "${RESULTS_DIR}/peak_analysis/${group}_merged.bed" \
            "$group"
    done
}

# Execute main function with error handling
{
    main
} 2>&1 | tee "${RESULTS_DIR}/logs/peak_analysis.log" 
