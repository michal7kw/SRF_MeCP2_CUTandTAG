#!/bin/bash
#SBATCH --job-name=to_bw
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/to_bw.err"
#SBATCH --output="logs/to_bw.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative/results_1b/aligned

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Create output directories
mkdir -p bigwig_files
mkdir -p logs

# Function to convert BAM to BigWig
bam_to_bigwig() {
    local bam_file=$1
    local sample_name=$(basename "$bam_file" .bam)
    
    echo "Processing $sample_name..."
    
    # Generate BigWig file
    bamCoverage \
        --bam "$bam_file" \
        --outFileName "bigwig_files/${sample_name}.bw" \
        --binSize 10 \
        --normalizeUsing RPKM \
        --numberOfProcessors 32 \
        --ignoreForNormalization chrX chrY chrM \
        2> "logs/${sample_name}.log"
}

# Process each BAM file (excluding dedup files)
for bam_file in *.bam; do
    if [[ ! $bam_file == *"dedup"* ]]; then
        bam_to_bigwig "$bam_file"
    fi
done

# Create a README file
cat > bigwig_files/README.txt << EOL
BigWig files prepared for IGV visualization:

File generation parameters:
- Bin size: 10 bp
- Normalization: RPKM
- Excluded from normalization: chrX, chrY, chrM

To use in IGV:
1. Open IGV on your local computer
2. Go to File -> Load from File
3. Select the BigWig (.bw) files you want to visualize

Note: Make sure to load the appropriate reference genome in IGV before loading the BigWig files.
EOL

echo "Conversion complete. BigWig files are in the bigwig_files directory"