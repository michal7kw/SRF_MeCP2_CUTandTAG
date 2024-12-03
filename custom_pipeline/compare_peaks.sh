#!/bin/bash
#SBATCH --job-name=compare_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/compare_peaks_%j.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/compare_peaks_%j.out"

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Create directories for analysis
PEAKS_DIR="results/peaks"
ANALYSIS_DIR="results/peak_analysis"
mkdir -p "$ANALYSIS_DIR"

# Function to process bedGraph files
process_bedgraph() {
    local SAMPLE=$1
    echo "Processing $SAMPLE..."
    
    # Sort bedGraph
    sort -k1,1 -k2,2n "$PEAKS_DIR/${SAMPLE}_treat_pileup.bdg" > \
        "$ANALYSIS_DIR/${SAMPLE}_sorted_treat_pileup.bdg"
        
    # Calculate peak sizes and signal intensities
    awk -v sample="$SAMPLE" 'BEGIN{OFS="\t"} 
        {
            size=$3-$2
            signal=$4
            print sample, size, signal
        }' "$PEAKS_DIR/${SAMPLE}_peaks.broadPeak" > "$ANALYSIS_DIR/${SAMPLE}_peak_metrics.txt"
}

# Function to compare peaks between conditions
compare_conditions() {
    local COND1=$1
    local COND2=$2
    
    echo "Comparing $COND1 vs $COND2..."
    
    # Find overlapping peaks
    bedtools intersect -a "$PEAKS_DIR/${COND1}_peaks.broadPeak" \
        -b "$PEAKS_DIR/${COND2}_peaks.broadPeak" -wo > \
        "$ANALYSIS_DIR/${COND1}_vs_${COND2}_overlaps.bed"
    
    # Analyze overlapping peaks
    awk 'BEGIN{OFS="\t"} 
        {
            size1=$3-$2
            size2=$13-$12
            signal1=$7
            signal2=$17
            ratio=signal1/signal2
            print $1, $2, $3, size1, size2, signal1, signal2, ratio
        }' "$ANALYSIS_DIR/${COND1}_vs_${COND2}_overlaps.bed" > \
        "$ANALYSIS_DIR/${COND1}_vs_${COND2}_comparison.txt"
    
    # Generate summary statistics
    Rscript - <<EOF
    data <- read.table("$ANALYSIS_DIR/${COND1}_vs_${COND2}_comparison.txt", 
                      header=FALSE, 
                      col.names=c("chr", "start", "end", "size1", "size2", 
                                "signal1", "signal2", "ratio"))
    
    # Calculate statistics
    stats <- data.frame(
        mean_size_ratio = mean(data\$size1/data\$size2),
        median_size_ratio = median(data\$size1/data\$size2),
        mean_signal_ratio = mean(data\$ratio),
        median_signal_ratio = median(data\$ratio)
    )
    
    # Save statistics
    write.table(stats, 
                "$ANALYSIS_DIR/${COND1}_vs_${COND2}_stats.txt", 
                quote=FALSE, 
                sep="\t")
    
    # Generate plots
    pdf("$ANALYSIS_DIR/${COND1}_vs_${COND2}_plots.pdf")
    
    # Size comparison
    plot(data\$size1, data\$size2, 
         xlab="${COND1} peak size", 
         ylab="${COND2} peak size",
         main="Peak Size Comparison")
    abline(0,1, col="red")
    
    # Signal comparison
    plot(data\$signal1, data\$signal2, 
         xlab="${COND1} signal", 
         ylab="${COND2} signal",
         main="Peak Signal Comparison")
    abline(0,1, col="red")
    
    # Ratio distributions
    hist(log2(data\$ratio), 
         main="Log2 Signal Ratio Distribution",
         xlab="Log2(${COND1}/${COND2})")
    
    dev.off()
EOF
}

# Process all samples
for bdg in "$PEAKS_DIR"/*_treat_pileup.bdg; do
    SAMPLE=$(basename "$bdg" _treat_pileup.bdg)
    process_bedgraph "$SAMPLE"
done

# Compare conditions (adjust these based on your sample names)
# Example: compare_conditions "NSCv1" "NeuV1"
# Add your condition comparisons here

echo "Peak analysis completed. Results in $ANALYSIS_DIR" 