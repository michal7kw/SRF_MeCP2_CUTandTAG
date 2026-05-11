library(rtracklayer)
library(GenomicRanges)
library(genomation)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
bigwig_file <- args[1]
peaks_file <- args[2]
output_plot <- args[3]
upstream <- as.numeric(args[4])
downstream <- as.numeric(args[5])

# Read data
peaks <- import(peaks_file)
signal <- import(bigwig_file)

# Create ScoreMatrix
sm <- ScoreMatrix(target=signal, 
                  windows=resize(peaks, width=upstream+downstream, fix="center"),
                  strand.aware=TRUE)

# Calculate average profile
profile <- colMeans(sm, na.rm=TRUE)
positions <- seq(-upstream, downstream, length.out=length(profile))

# Create plot
pdf(output_plot)
plot(positions, profile, type="l", 
     xlab="Distance from peak center (bp)",
     ylab="Average signal",
     main="Metagene plot")
dev.off() 