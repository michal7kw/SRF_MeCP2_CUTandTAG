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
#SBATCH --error="logs/qc_1.err"
#SBATCH --output="logs/qc_1.out"
#SBATCH --array=0-11

# Set working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Get sample names
EXOGENOUS_SAMPLES=($(ls ../DATA/EXOGENOUS/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))
ENDOGENOUS_SAMPLES=($(ls ../DATA/ENDOGENOUS/*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//'))
ALL_SAMPLES=("${EXOGENOUS_SAMPLES[@]}" "${ENDOGENOUS_SAMPLES[@]}")

# Get current sample
SAMPLE=${ALL_SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Create output directories
mkdir -p results_1/qc/fragment_sizes
mkdir -p results_1/qc/tss_enrichment
mkdir -p results_1/qc/frip
mkdir -p results_1/bigwig

# 1. Fragment size distribution
samtools view -f 0x2 results_1/aligned/${SAMPLE}.bam | awk -F'\t' '{print sqrt($9^2)}' | \
    sort | uniq -c | awk '{print $2,$1}' > results_1/qc/fragment_sizes/${SAMPLE}_sizes.txt

# 2. Generate bigWig for visualization and TSS enrichment
bamCoverage --bam results_1/aligned/${SAMPLE}.bam \
    --outFileName results_1/bigwig/${SAMPLE}.bw \
    --binSize 10 \
    --normalizeUsing RPKM \
    --numberOfProcessors 8

# 4. Calculate FRiP score
if [[ -f "results_1/peaks/seacr/${SAMPLE}.stringent.bed" ]]; then
    total_reads=$(samtools view -c -F 4 results_1/aligned/${SAMPLE}.bam)
    reads_in_peaks=$(bedtools intersect -a results_1/aligned/${SAMPLE}.bam \
        -b results_1/peaks/seacr/${SAMPLE}.stringent.bed -u -f 0.20 | samtools view -c)
    echo -e "${SAMPLE}\t${reads_in_peaks}\t${total_reads}" > results_1/qc/frip/${SAMPLE}_frip.txt
fi 