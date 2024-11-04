#!/bin/bash
#SBATCH --job-name=Heatmaps
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_ChipSeq/Visualization/logs/heatmaps.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_ChipSeq/Visualization/logs/heatmaps.out"

echo "Starting heatmap creation script..."

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

echo "Conda environment activated."

# Define input directory
INPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_ChipSeq/custom_pipeline/results"

# Define output directory
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_ChipSeq/Visualization/output"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR
echo "Output directory created: $OUTPUT_DIR"

# Define sample groups
endogenous_samples=(NeuM2 NeuM3 NSCM1 NSCM2 NSCM3 IgM)
exogenous_samples=(NeuV1 NeuV2 NeuV3 NSCv1 NSCv2 NSCv3)

echo "Sample groups defined."

# Combine all peak files
# This command concatenates all narrowPeak files, sorts them by chromosome and start position,
# then merges overlapping regions to create a single BED file of all unique peak regions
echo "Combining all peak files..."
cat ${INPUT_DIR}/peaks/*_peaks.narrowPeak | sort -k1,1 -k2,2n | bedtools merge > ${OUTPUT_DIR}/all_peaks_merged.bed
echo "All peaks merged into: ${OUTPUT_DIR}/all_peaks_merged.bed"

# Create a matrix of read coverages for endogenous samples
# This command uses multiBamSummary from deepTools to compute read coverages over the merged peak regions
# for all endogenous samples. It outputs a numpy array (.npz) file for further analysis.
echo "Creating matrix of read coverages for endogenous samples..."
multiBamSummary BED-file --BED ${OUTPUT_DIR}/all_peaks_merged.bed \
    --bamfiles ${INPUT_DIR}/aligned/NeuM2.bam ${INPUT_DIR}/aligned/NeuM3.bam \
    ${INPUT_DIR}/aligned/NSCM1.bam ${INPUT_DIR}/aligned/NSCM2.bam \
    ${INPUT_DIR}/aligned/NSCM3.bam ${INPUT_DIR}/aligned/IgM.bam \
    -o ${OUTPUT_DIR}/endogenous_results.npz \
    --labels ${endogenous_samples[@]}
echo "Endogenous samples matrix created: ${OUTPUT_DIR}/endogenous_results.npz"

# Create a matrix of read coverages for exogenous samples
echo "Creating matrix of read coverages for exogenous samples..."
multiBamSummary BED-file --BED ${OUTPUT_DIR}/all_peaks_merged.bed \
    --bamfiles ${INPUT_DIR}/aligned/NeuV1.bam ${INPUT_DIR}/aligned/NeuV2.bam \
    ${INPUT_DIR}/aligned/NeuV3.bam ${INPUT_DIR}/aligned/NSCv1.bam \
    ${INPUT_DIR}/aligned/NSCv2.bam ${INPUT_DIR}/aligned/NSCv3.bam \
    -o ${OUTPUT_DIR}/exogenous_results.npz \
    --labels ${exogenous_samples[@]}
echo "Exogenous samples matrix created: ${OUTPUT_DIR}/exogenous_results.npz"

# Plot heatmaps for endogenous samples
# This command uses plotHeatmap from deepTools to create a heatmap of peak intensities
# --colorMap RdYlBu sets the color scheme
# --whatToShow specifies to show both the heatmap and colorbar
# --zMin and --zMax set the range for color scaling
# --sortRegions descend sorts the regions by decreasing intensity
echo "Plotting heatmap for endogenous samples..."
plotHeatmap -m ${OUTPUT_DIR}/endogenous_results.npz -out ${OUTPUT_DIR}/endogenous_heatmap.png \
    --colorMap RdYlBu --whatToShow 'heatmap and colorbar' \
    --zMin -3 --zMax 3 --sortRegions descend \
    --plotTitle "Endogenous Samples Peak Intensity"
echo "Endogenous samples heatmap created: ${OUTPUT_DIR}/endogenous_heatmap.png"

# Plot heatmaps for exogenous samples
echo "Plotting heatmap for exogenous samples..."
plotHeatmap -m ${OUTPUT_DIR}/exogenous_results.npz -out ${OUTPUT_DIR}/exogenous_heatmap.png \
    --colorMap RdYlBu --whatToShow 'heatmap and colorbar' \
    --zMin -3 --zMax 3 --sortRegions descend \
    --plotTitle "Exogenous Samples Peak Intensity"
echo "Exogenous samples heatmap created: ${OUTPUT_DIR}/exogenous_heatmap.png"

# Create a matrix of read coverages for all samples
echo "Creating matrix of read coverages for all samples..."
multiBamSummary BED-file --BED ${OUTPUT_DIR}/all_peaks_merged.bed \
    --bamfiles ${INPUT_DIR}/aligned/*.bam \
    -o ${OUTPUT_DIR}/all_samples_results.npz \
    --labels ${endogenous_samples[@]} ${exogenous_samples[@]}
echo "All samples matrix created: ${OUTPUT_DIR}/all_samples_results.npz"

# Create a correlation heatmap for all samples
# This command uses plotCorrelation from deepTools to create a correlation heatmap
# --corMethod spearman specifies to use Spearman correlation
# --skipZeros ignores regions with no coverage
# --whatToPlot heatmap specifies to create a heatmap
# --colorMap RdYlBu sets the color scheme
# --plotNumbers adds correlation values to the heatmap
echo "Creating correlation heatmap for all samples..."
plotCorrelation -in ${OUTPUT_DIR}/all_samples_results.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o ${OUTPUT_DIR}/correlation_heatmap.png \
    --outFileCorMatrix ${OUTPUT_DIR}/correlation_matrix.tab
echo "Correlation heatmap created: ${OUTPUT_DIR}/correlation_heatmap.png"
echo "Correlation matrix saved: ${OUTPUT_DIR}/correlation_matrix.tab"

echo "Heatmap creation script completed successfully."