#!/bin/bash
#SBATCH --job-name=metaprofiles
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it

# Load required modules or activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Set working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative"
cd $WORKDIR

# Create output directory
mkdir -p results_1/metaprofiles

# Define file paths
GTF="/path/to/mm10.gtf"  # You'll need to specify the path to your mm10 GTF file
BLACKLIST="/path/to/mm10-blacklist.bed"  # Optional: specify path to blacklist regions

# Define sample groups
NEU_MECP2_SAMPLES=(
    "results_1/bigwig/NeuM2.bw"
    "results_1/bigwig/NeuM3.bw"
)

NEU_VECTOR_SAMPLES=(
    "results_1/bigwig/NeuV1.bw"
    "results_1/bigwig/NeuV2.bw"
    "results_1/bigwig/NeuV3.bw"
)

NSC_MECP2_SAMPLES=(
    "results_1/bigwig/NSCM1.bw"
    "results_1/bigwig/NSCM2.bw"
    "results_1/bigwig/NSCM3.bw"
)

NSC_VECTOR_SAMPLES=(
    "results_1/bigwig/NSCv1.bw"
    "results_1/bigwig/NSCv2.bw"
    "results_1/bigwig/NSCv3.bw"
)

# Step 1: Create regions file for analysis
echo "Creating regions file..."
computeMatrix scale-regions \
    -S "${NEU_MECP2_SAMPLES[@]}" "${NEU_VECTOR_SAMPLES[@]}" \
       "${NSC_MECP2_SAMPLES[@]}" "${NSC_VECTOR_SAMPLES[@]}" \
    -R $GTF \
    --beforeRegionStartLength 5000 \
    --regionBodyLength 5000 \
    --afterRegionStartLength 5000 \
    --skipZeros \
    --numberOfProcessors 16 \
    --transcriptID gene \
    --metagene \
    -o results_1/metaprofiles/matrix.gz \
    --outFileSortedRegions results_1/metaprofiles/regions.bed

# Step 2: Create profile plot
echo "Generating profile plot..."
plotProfile \
    -m results_1/metaprofiles/matrix.gz \
    -out results_1/metaprofiles/metaprofile.pdf \
    --perGroup \
    --colors "#FF0000" "#0000FF" "#00FF00" "#000000" \
    --plotTitle "Mecp2 binding profile around genes" \
    --samplesLabel "NEU-Mecp2" "NEU-Vector" "NSC-Mecp2" "NSC-Vector" \
    --regionsLabel "Genes" \
    --yAxisLabel "Average signal" \
    --xAxisLabel "Distance from TSS (bp)"

# Step 3: Create heatmap
echo "Generating heatmap..."
plotHeatmap \
    -m results_1/metaprofiles/matrix.gz \
    -out results_1/metaprofiles/heatmap.pdf \
    --colorMap RdYlBu \
    --whatToShow 'plot, heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 4 \
    --zMin 0 \
    --zMax 10

# Step 4: Generate TSS-centered profile
echo "Generating TSS-centered profile..."
computeMatrix reference-point \
    -S "${NEU_MECP2_SAMPLES[@]}" "${NEU_VECTOR_SAMPLES[@]}" \
       "${NSC_MECP2_SAMPLES[@]}" "${NSC_VECTOR_SAMPLES[@]}" \
    -R $GTF \
    --referencePoint TSS \
    --beforeRegionStartLength 5000 \
    --afterRegionStartLength 5000 \
    --skipZeros \
    --numberOfProcessors 16 \
    -o results_1/metaprofiles/tss_matrix.gz

plotProfile \
    -m results_1/metaprofiles/tss_matrix.gz \
    -out results_1/metaprofiles/tss_profile.pdf \
    --perGroup \
    --colors "#FF0000" "#0000FF" "#00FF00" "#000000" \
    --plotTitle "Mecp2 binding profile around TSS" \
    --samplesLabel "NEU-Mecp2" "NEU-Vector" "NSC-Mecp2" "NSC-Vector" \
    --regionsLabel "TSS" \
    --yAxisLabel "Average signal" \
    --xAxisLabel "Distance from TSS (bp)"

# Step 5: Generate statistics
echo "Generating statistics..."
plotProfileWithStderr \
    -m results_1/metaprofiles/matrix.gz \
    -out results_1/metaprofiles/profile_with_stderr.pdf \
    --perGroup \
    --plotTitle "Mecp2 binding profile with standard error" \
    --samplesLabel "NEU-Mecp2" "NEU-Vector" "NSC-Mecp2" "NSC-Vector" \
    --regionsLabel "Genes" \
    --yAxisLabel "Average signal" \
    --xAxisLabel "Distance from TSS (bp)"

echo "Analysis complete. Check results in results_1/metaprofiles/" 