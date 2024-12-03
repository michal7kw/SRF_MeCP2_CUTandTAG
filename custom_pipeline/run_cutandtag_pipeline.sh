#!/bin/bash

# Master script to run the entire CUT&Tag analysis pipeline

# 1. Quality Control
sbatch fastqc_array.sh
wait_for_jobs "fastqc"

# 2. Trimming
sbatch trim_reads_array.sh
wait_for_jobs "trim_reads"

# 3. Alignment
sbatch align_array.sh
wait_for_jobs "align"

# 4. Peak Calling
sbatch call_peaks_seacr.sh
wait_for_jobs "peaks_seacr"

# 5. QC Metrics
sbatch qc_metrics.sh
wait_for_jobs "qc_metrics"

# 6. Generate QC Report
Rscript generate_qc_report.R

# 7. Differential Analysis
Rscript differential_analysis.R

# 8. Visualization
Rscript visualization_workflow.R

# Function to wait for jobs to complete
wait_for_jobs() {
    while squeue -u $USER | grep -q "$1"; do
        sleep 60
    done
} 