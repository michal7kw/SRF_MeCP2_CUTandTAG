#!/bin/bash
#SBATCH --job-name=final_list
#SBATCH --account kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --mail-type=ALL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=kubacki.michal@hst.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/final_list.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/final_list.out"

echo "Starting peak analysis pipeline..."

# Load the appropriate conda environment (if needed)
echo "Loading conda environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

echo "Changing to results directory..."
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/results

# Create output directory for consensus peaks
echo "Creating consensus peaks directory..."
mkdir -p consensus_peaks

# Function to merge peaks from replicates
merge_peaks() {
    condition=$1
    tissue=$2
    files=$3
    
    echo "Processing ${tissue} ${condition} peaks..."
    echo "Input files: $files"
    
    # Merge all replicates
    echo "Sorting peaks..."
    cat $files | sort -k1,1 -k2,2n > tmp_sorted.bed
    
    # Merge overlapping peaks and require support from at least 2 replicates
    echo "Merging overlapping peaks with minimum 2 replicate support..."
    bedtools merge -i tmp_sorted.bed -c 1 -o count | awk '$4 >= 2' > consensus_peaks/${tissue}_${condition}_consensus.bed
    
    # Report number of peaks
    peak_count=$(wc -l < consensus_peaks/${tissue}_${condition}_consensus.bed)
    echo "Generated ${peak_count} consensus peaks for ${tissue} ${condition}"
    
    # Clean up
    echo "Cleaning up temporary files..."
    rm tmp_sorted.bed
}

echo "=== Processing Neuron Endogenous peaks ==="
merge_peaks "Endo" "Neuron" "./peaks/NeuM2_peaks.narrowPeak ./peaks/NeuM3_peaks.narrowPeak"

echo "=== Processing Neuron Exogenous peaks ==="
merge_peaks "Exo" "Neuron" "./peaks/NeuV1_peaks.narrowPeak ./peaks/NeuV2_peaks.narrowPeak ./peaks/NeuV3_peaks.narrowPeak"

echo "=== Processing NSC Endogenous peaks ==="
merge_peaks "Endo" "NSC" "./peaks/NSCM1_peaks.narrowPeak ./peaks/NSCM2_peaks.narrowPeak ./peaks/NSCM3_peaks.narrowPeak"

echo "=== Processing NSC Exogenous peaks ==="
merge_peaks "Exo" "NSC" "./peaks/NSCv1_peaks.narrowPeak ./peaks/NSCv2_peaks.narrowPeak ./peaks/NSCv3_peaks.narrowPeak"

# Generate statistics
echo "Generating peak statistics..."
echo "=== Peak Statistics ===" > peak_statistics.txt
echo "Generated on: $(date)" >> peak_statistics.txt
echo "----------------------------------------" >> peak_statistics.txt
for file in consensus_peaks/*.bed; do
    count=$(wc -l < "$file")
    echo "$file: $count peaks" >> peak_statistics.txt
    echo "$file: $count peaks"
done

echo "=== Generating intersection between conditions ==="
# Create intersection between conditions
for tissue in Neuron NSC; do
    echo "Processing ${tissue} intersections..."
    
    echo "Finding shared peaks between Endo and Exo for ${tissue}..."
    bedtools intersect -a ./consensus_peaks/${tissue}_Endo_consensus.bed \
                      -b ./consensus_peaks/${tissue}_Exo_consensus.bed \
                      > ./consensus_peaks/${tissue}_shared_peaks.bed
    shared_count=$(wc -l < ./consensus_peaks/${tissue}_shared_peaks.bed)
    echo "${tissue} shared peaks: ${shared_count}"
                      
    echo "Finding Endo-specific peaks for ${tissue}..."
    bedtools subtract -a ./consensus_peaks/${tissue}_Endo_consensus.bed \
                     -b ./consensus_peaks/${tissue}_Exo_consensus.bed \
                     > ./consensus_peaks/${tissue}_Endo_specific.bed
    endo_specific_count=$(wc -l < ./consensus_peaks/${tissue}_Endo_specific.bed)
    echo "${tissue} Endo-specific peaks: ${endo_specific_count}"
                     
    echo "Finding Exo-specific peaks for ${tissue}..."
    bedtools subtract -a ./consensus_peaks/${tissue}_Exo_consensus.bed \
                     -b ./consensus_peaks/${tissue}_Endo_consensus.bed \
                     > ./consensus_peaks/${tissue}_Exo_specific.bed
    exo_specific_count=$(wc -l < ./consensus_peaks/${tissue}_Exo_specific.bed)
    echo "${tissue} Exo-specific peaks: ${exo_specific_count}"
done

echo "Peak analysis pipeline completed successfully!"