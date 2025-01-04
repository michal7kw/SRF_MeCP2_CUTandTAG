#!/bin/bash
#SBATCH --job-name=2_mecp2_meth
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/methylation_analysis_2.err"
#SBATCH --output="logs/methylation_analysis_2.out"

WORKING_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/Methylation"
cd $WORKING_DIR || exit 1

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "conda could not be found"
    exit 1
fi

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Check if required Python packages are installed
python3 -c "import pandas, numpy, seaborn, matplotlib, pyBigWig, pysam, pyranges" || {
    echo "Missing required Python packages"
    exit 1
}

# Create directories
mkdir -p logs plots/methylation plots/expression plots/integrated

# Set environment variables for parallel processing
export MKL_NUM_THREADS=16
export NUMEXPR_NUM_THREADS=16
export OMP_NUM_THREADS=16
export PYTHONUNBUFFERED=1
export PYTHONPATH="${PYTHONPATH}:/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG"

# Check if required input files exist
required_files=(
    "../DATA/mm10.fa"
    "../DATA/cpg_islands.bed"
    "../DATA/gencode.vM10.annotation.gtf"
    "../iterative_alternative/DATA/DEA_NEU.csv"
    "../iterative_alternative/DATA/DEA_NSC.csv"
)

for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Required file not found: $file"
        exit 1
    fi
done

# Run the analysis
echo "Starting methylation analysis at $(date)"
# LOG_FILE="logs/methylation_analysis_2_$(date +%Y%m%d_%H%M%S).log"

#  2>&1 | tee "$LOG_FILE"
python3 analyze_methylation_expression.py
# PYTHON_EXIT=${PIPESTATUS[0]}

# if [ $PYTHON_EXIT -eq 0 ]; then
#     echo "Analysis completed successfully at $(date)"
    
#     # Create results archive
#     ARCHIVE_NAME="plots_$(date +%Y%m%d_%H%M%S).tar.gz"
#     tar -czf "$ARCHIVE_NAME" plots/
    
#     # Generate summary
#     echo -e "\nSummary of results:"
#     echo "===================="
#     echo "Number of plots generated: $(find plots/ -name "*.pdf" | wc -l)"
    
#     if [ -f "plots/analysis_results.txt" ]; then
#         echo "Analysis results file size: $(du -h plots/analysis_results.txt | cut -f1)"
#         echo "Results have been archived to: $ARCHIVE_NAME"
#     else
#         echo "Warning: analysis_results.txt not found"
#     fi
    
#     # Check if any plots were generated
#     if [ $(find plots/ -name "*.pdf" | wc -l) -eq 0 ]; then
#         echo "Warning: No plots were generated"
#         exit 1
#     fi
# else
#     echo "Analysis failed at $(date)"
#     echo "Check log file for details: $LOG_FILE"
#     exit 1
# fi
