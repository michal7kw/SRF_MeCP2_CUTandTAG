library(rtracklayer)
library(GenomicRanges)
library(genomation)
library(ComplexHeatmap)
library(circlize)

args <- commandArgs(trailingOnly = TRUE)
bigwig_file <- args[1]
peaks_file <- args[2]
output_plot <- args[3]
window_size <- as.numeric(args[4])
bin_size <- as.numeric(args[5])

# Read data
peaks <- import(peaks_file)
signal <- import(bigwig_file)

# Create ScoreMatrix
sm <- ScoreMatrix(target=signal,
                  windows=resize(peaks, width=window_size, fix="center"),
                  bin.num=window_size/bin_size)

# Sort regions by mean signal
region_means <- rowMeans(sm, na.rm=TRUE)
sm <- sm[order(region_means, decreasing=TRUE),]

# Create heatmap
pdf(output_plot, width=8, height=10)
col_fun <- colorRamp2(c(min(sm), mean(sm), max(sm)), 
                      c("blue", "white", "red"))
Heatmap(sm,
        name="Signal",
        col=col_fun,
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        show_row_names=FALSE,
        column_title="Distance from peak center")
dev.off() 