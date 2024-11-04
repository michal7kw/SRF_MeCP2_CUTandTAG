#!/bin/bash
#SBATCH --job-name=Merge_replicates
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_ChipSeq/custom_pipeline/logs/merge_replicates.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_ChipSeq/custom_pipeline/logs/merge_replicates.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_ChipSeq/custom_pipeline/results


# Create output directories
mkdir -p consensus_peaks bigwig

# Function to combine BAM files and create bigwig
combine_bams_to_bigwig() {
    local condition=$1
    local tissue=$2
    shift 2
    local bams=("$@")  # Get all remaining arguments as BAM files
    
    local merged_bam="consensus_peaks/${condition}_${tissue}_merged.bam"
    local bigwig_file="bigwig/${condition}_${tissue}.bw"
    
    # Check if final bigwig exists
    if [[ -f "$bigwig_file" ]]; then
        echo "Bigwig file already exists for ${condition} ${tissue}, skipping..."
        return
    fi
    
    echo "Processing ${condition} ${tissue} with BAMs: ${bams[@]}"
    
    # Merge BAM files if needed
    if [[ ! -f "$merged_bam" ]]; then
        samtools merge -f "$merged_bam" "${bams[@]}"
        samtools index "$merged_bam"
    fi
    
    # Create bigwig file (normalized to CPM)
    bamCoverage --bam "$merged_bam" \
                --outFileName "$bigwig_file" \
                --binSize 10 \
                --normalizeUsing CPM \
                --numberOfProcessors 32
}

# Function to call consensus peaks
call_consensus_peaks() {
    local condition=$1
    local tissue=$2
    local control=$3
    shift 3
    local peak_files=("$@")  # Get all remaining arguments as peak files
    
    local consensus_prefix="consensus_peaks/${condition}_${tissue}_consensus"
    local final_peaks="${consensus_prefix}_peaks.narrowPeak"
    
    # Check if final peaks file exists
    if [[ -f "$final_peaks" ]]; then
        echo "Consensus peaks already exist for ${condition} ${tissue}, skipping..."
        return
    fi
    
    echo "Calling consensus peaks for ${condition} ${tissue}"
    echo "Peak files: ${peak_files[@]}"
    
    local all_peaks="consensus_peaks/${condition}_${tissue}_all_peaks.bed"
    local merged_peaks="consensus_peaks/${condition}_${tissue}_merged_peaks.bed"
    
    # Create merged peak file
    cat "${peak_files[@]}" > "$all_peaks"
    
    # Sort and merge overlapping peaks
    sortBed -i "$all_peaks" | \
    mergeBed -i stdin > "$merged_peaks"
    
    # Call peaks on merged BAM using merged peaks as candidates
    macs2 callpeak \
        -t "consensus_peaks/${condition}_${tissue}_merged.bam" \
        -c "$control" \
        -f BAMPE \
        -g hs \
        --keep-dup auto \
        -n "${condition}_${tissue}_consensus" \
        --outdir consensus_peaks
}

# Process Exogenous samples
## Neurons
combine_bams_to_bigwig "Exogenous" "Neuron" \
    "aligned/NeuV1.bam" \
    "aligned/NeuV2.bam" \
    "aligned/NeuV3.bam"

call_consensus_peaks "Exogenous" "Neuron" "aligned/IgM.bam" \
    "peaks/NeuV1_peaks.narrowPeak" \
    "peaks/NeuV2_peaks.narrowPeak" \
    "peaks/NeuV3_peaks.narrowPeak"

## NSCs
combine_bams_to_bigwig "Exogenous" "NSC" \
    "aligned/NSCv1.bam" \
    "aligned/NSCv2.bam" \
    "aligned/NSCv3.bam"

call_consensus_peaks "Exogenous" "NSC" "aligned/IgM.bam" \
    "peaks/NSCv1_peaks.narrowPeak" \
    "peaks/NSCv2_peaks.narrowPeak" \
    "peaks/NSCv3_peaks.narrowPeak"

# Process Endogenous samples
## Neurons (only replicates 2 and 3)
combine_bams_to_bigwig "Endogenous" "Neuron" \
    "aligned/NeuM2.bam" \
    "aligned/NeuM3.bam"

call_consensus_peaks "Endogenous" "Neuron" "aligned/IgM.bam" \
    "peaks/NeuM2_peaks.narrowPeak" \
    "peaks/NeuM3_peaks.narrowPeak"

## NSCs
combine_bams_to_bigwig "Endogenous" "NSC" \
    "aligned/NSCM1.bam" \
    "aligned/NSCM2.bam" \
    "aligned/NSCM3.bam"

call_consensus_peaks "Endogenous" "NSC" "aligned/IgM.bam" \
    "peaks/NSCM1_peaks.narrowPeak" \
    "peaks/NSCM2_peaks.narrowPeak" \
    "peaks/NSCM3_peaks.narrowPeak"

# Process Control (IgG)
combine_bams_to_bigwig "Control" "IgG" "aligned/IgM.bam"