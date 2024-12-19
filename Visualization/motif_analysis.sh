#!/bin/bash
#SBATCH --job-name=Motif_Analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Visualization/logs/motif_analysis.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Visualization/logs/motif_analysis.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake


# Define input directory
INPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/results_5_new_005"

# Define output directory
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Visualization/results/"


DATA_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/DATA"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Define sample groups
endogenous_samples=(NeuM2 NeuM3 NSCM1 NSCM2 NSCM3)
exogenous_samples=(NeuV1 NeuV2 NeuV3 NSCv1 NSCv2 NSCv3)

# Function to run MEME on a given sample
run_meme() {
    sample=$1
    condition=$2
    
    # Extract sequences from peak regions
    bedtools getfasta -fi ${DATA_DIR}/mm10.fa -bed ${INPUT_DIR}/peaks/${sample}_peaks.narrowPeak \
        -fo ${OUTPUT_DIR}/${sample}_peak_sequences.fa

    # Run MEME
    meme ${OUTPUT_DIR}/${sample}_peak_sequences.fa -dna -oc ${OUTPUT_DIR}/meme_output_${condition}_${sample} \
        -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 6 -maxw 50 \
        -objfun classic -revcomp
}

# Run MEME for endogenous samples
for sample in "${endogenous_samples[@]}"; do
    run_meme $sample "endogenous"
done

# Run MEME for exogenous samples
for sample in "${exogenous_samples[@]}"; do
    run_meme $sample "exogenous"
done

# Combine MEME results
cat ${OUTPUT_DIR}/meme_output_endogenous_*/meme.txt > ${OUTPUT_DIR}/combined_endogenous_motifs.txt
cat ${OUTPUT_DIR}/meme_output_exogenous_*/meme.txt > ${OUTPUT_DIR}/combined_exogenous_motifs.txt

# Compare motifs between conditions using TOMTOM
tomtom -oc ${OUTPUT_DIR}/tomtom_output ${OUTPUT_DIR}/combined_endogenous_motifs.txt ${OUTPUT_DIR}/combined_exogenous_motifs.txt