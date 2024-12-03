#!/bin/bash
#SBATCH --job-name=call_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/peaks.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/logs/peaks.out"

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Source the module system and load required modules
source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/Modules/init/bash
module load macs2/2.2.7.1
module load parallel
module load samtools

# Print environment info
echo "Using MACS2 from: $(which macs2)"
echo "MACS2 version: $(macs2 --version)"

# Set paths and parameters
BAM_DIR="results/aligned"
OUTPUT_DIR="results/peaks"
CONTROL_BAM="results/aligned/IgM.bam"
GENOME_SIZE="hs"
QVALUE="0.01"
FORMAT="BAMPE"
THREADS=16

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Process each BAM file
find "$BAM_DIR" -name "*.bam" ! -name "$(basename $CONTROL_BAM)" -print0 | while IFS= read -r -d '' BAM_FILE; do
    SAMPLE=$(basename "$BAM_FILE" .bam)
    
    echo "Processing sample: $SAMPLE"
    
    # Run MACS2
    macs2 callpeak \
        -t "$BAM_FILE" \
        -c "$CONTROL_BAM" \
        -f "$FORMAT" \
        -g "$GENOME_SIZE" \
        -n "$SAMPLE" \
        --outdir "$OUTPUT_DIR" \
        -q "$QVALUE" \
        --nomodel \
        --call-summits \
        --broad \
        2>&1 | tee "$OUTPUT_DIR/${SAMPLE}_macs2.log"
        
    if [ $? -eq 0 ] && [ -f "$OUTPUT_DIR/${SAMPLE}_peaks.narrowPeak" ]; then
        echo "Successfully processed $SAMPLE"
    else
        echo "Error processing $SAMPLE"
    fi
done

echo "Peak calling completed for all samples"