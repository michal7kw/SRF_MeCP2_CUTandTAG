#!/bin/bash
#SBATCH --job-name=convert_narrowPeak_to_BED
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Visualization/logs/genome_browser.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Visualization/logs/genome_browser.out"

# Create output directory if it doesn't exist
mkdir -p /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Visualization/output

# Convert narrowPeak to BED format
for file in /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/results/peaks/*_peaks.narrowPeak; do
    filename=$(basename "$file")
    output_file="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Visualization/output/${filename%.narrowPeak}.bed"
    cut -f 1-6 "$file" > "$output_file"
done