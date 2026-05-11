# For analyzing read cross-correlation
library(phantompeakqualtools)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
bam_file <- args[1]
plot_file <- args[2]
metrics_file <- args[3]

# Run cross-correlation analysis
cc_results <- run_spp(bam_file, plot=plot_file)

# Save metrics
write.table(cc_results, metrics_file, sep="\t", quote=FALSE)