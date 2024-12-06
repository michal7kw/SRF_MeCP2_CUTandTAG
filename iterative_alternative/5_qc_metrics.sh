#!/bin/bash
#SBATCH --job-name=qc_metrics
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/qc_%a.err"
#SBATCH --output="logs/qc_%a.out"
#SBATCH --array=0-11

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline/iterative_processing
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Source modules
# source /etc/profile.d/modules.sh
# module load samtools/1.13
# module load bedtools/2.29.1
# module load R/4.1.0
# module load deeptools/3.5.1

# Get sample names
EXOGENOUS_SAMPLES=($(ls ../DATA/EXOGENOUS/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))
ENDOGENOUS_SAMPLES=($(ls ../DATA/ENDOGENOUS/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))
ALL_SAMPLES=("${EXOGENOUS_SAMPLES[@]}" "${ENDOGENOUS_SAMPLES[@]}")

# Get current sample
SAMPLE=${ALL_SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Create output directories
mkdir -p results/qc/fragment_sizes
mkdir -p results/qc/tss_enrichment
mkdir -p results/qc/frip
mkdir -p results/bigwig

# 1. Fragment size distribution
samtools view -f 0x2 results/aligned/${SAMPLE}.bam | awk -F'\t' '{print sqrt($9^2)}' | \
    sort | uniq -c | awk '{print $2,$1}' > results/qc/fragment_sizes/${SAMPLE}_sizes.txt

# 2. Generate bigWig for visualization and TSS enrichment
bamCoverage --bam results/aligned/${SAMPLE}.bam \
    --outFileName results/bigwig/${SAMPLE}.bw \
    --binSize 10 \
    --normalizeUsing RPKM \
    --numberOfProcessors 8

# 4. Calculate FRiP score
if [[ -f "results/peaks/seacr/${SAMPLE}.stringent.bed" ]]; then
    total_reads=$(samtools view -c -F 4 results/aligned/${SAMPLE}.bam)
    reads_in_peaks=$(bedtools intersect -a results/aligned/${SAMPLE}.bam \
        -b results/peaks/seacr/${SAMPLE}.stringent.bed -u -f 0.20 | samtools view -c)
    echo -e "${SAMPLE}\t${reads_in_peaks}\t${total_reads}" > results/qc/frip/${SAMPLE}_frip.txt
fi 