#!/bin/bash
#SBATCH --job-name=tss_and_cpg_enrich
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/tss_and_cpg_enrich.err"
#SBATCH --output="logs/tss_and_cpg_enrich.out"

# This Bash script generates and plots the distribution of MeCP2 binding (both exogenous and endogenous) 
# around CpG islands and TSS (Transcription Start Sites) regions for two cell types (Neurons and NSCs) 
# using deepTools' computeMatrix and plotProfile commands.


# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Create output directories
RESULTS_DIR="results_1b"
OUTPUT_DIR="additional_scripts" 

mkdir -p ${OUTPUT_DIR}/distribution_plots

# Get sample lists for Neurons
NEURONS_EXO=($(ls ../DATA/EXOGENOUS/NeuV*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))
NEURONS_ENDO=($(ls ../DATA/ENDOGENOUS/NeuM*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))

# Get sample lists for NSCs
NSC_EXO=($(ls ../DATA/EXOGENOUS/NSCv*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))
NSC_ENDO=($(ls ../DATA/ENDOGENOUS/NSCM*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))

# Function to process samples and generate matrices
generate_matrices() {
    local cell_type=$1
    local region=$2
    local bed_file=$3
    local output_prefix=$4
    local exo_samples=("${!5}")
    local endo_samples=("${!6}")
    
    # Process exogenous samples
    exo_bigwigs=""
    for sample in "${exo_samples[@]}"; do
        exo_bigwigs="$exo_bigwigs ${RESULTS_DIR}/bigwig/${sample}.bw"
    done
    
    # Process endogenous samples
    endo_bigwigs=""
    for sample in "${endo_samples[@]}"; do
        endo_bigwigs="$endo_bigwigs ${RESULTS_DIR}/bigwig/${sample}.bw"
    done
    
    # Generate matrices
    computeMatrix reference-point \
        --referencePoint center \
        -b 3000 -a 3000 \
        -R ${bed_file} \
        -S $exo_bigwigs $endo_bigwigs \
        --skipZeros \
        --numberOfProcessors 8 \
        -o ${OUTPUT_DIR}/distribution_plots/${output_prefix}_matrix.gz
    
    # Plot profile
    plotProfile \
        -m ${OUTPUT_DIR}/distribution_plots/${output_prefix}_matrix.gz \
        -o ${OUTPUT_DIR}/distribution_plots/${output_prefix}_profile.png \
        --plotTitle "MeCP2 distribution at ${region} (${cell_type})" \
        --averageType mean \
        --perGroup \
        --colors blue lightblue \
        --samplesLabel "Exogenous" "Endogenous" \
        --legendLocation upper-right \
        --yMin 0 \
        --plotHeight 6 \
        --plotWidth 8
}

# Generate plots for Neurons
if [ -s "./DATA/cpg_islands.bed" ]; then
    generate_matrices "Neurons" "CpG Islands" "./DATA/cpg_islands.bed" "neurons_cpg_islands" NEURONS_EXO[@] NEURONS_ENDO[@]
else
    echo "Error: CpG islands bed file is missing or empty"
fi

if [ -s "./DATA/mm10_TSS.bed" ]; then
    generate_matrices "Neurons" "TSS" "./DATA/mm10_TSS.bed" "neurons_tss" NEURONS_EXO[@] NEURONS_ENDO[@]
else
    echo "Error: TSS bed file is missing or empty"
fi

# Generate plots for NSCs
if [ -s "./DATA/cpg_islands.bed" ]; then
    generate_matrices "NSCs" "CpG Islands" "./DATA/cpg_islands.bed" "nscs_cpg_islands" NSC_EXO[@] NSC_ENDO[@]
else
    echo "Error: CpG islands bed file is missing or empty"
fi

if [ -s "./DATA/mm10_TSS.bed" ]; then
    generate_matrices "NSCs" "TSS" "./DATA/mm10_TSS.bed" "nscs_tss" NSC_EXO[@] NSC_ENDO[@]
else
    echo "Error: TSS bed file is missing or empty"
fi