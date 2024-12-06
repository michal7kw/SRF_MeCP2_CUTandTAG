# For generating final HTML report
library(rmarkdown)
library(knitr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
qc_files <- args[1]
peak_metrics <- args[2]
output_file <- args[3]

# Read QC metrics
qc_data <- read_qc_files(qc_files)
peak_data <- read.table(peak_metrics)

# Generate report
rmarkdown::render(
    "report_template.Rmd",
    output_file = output_file,
    params = list(
        qc_data = qc_data,
        peak_data = peak_data
    )
)