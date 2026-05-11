library(optparse)
library(ggplot2)
library(dplyr)
library(GenomicRanges)

# Parse command line arguments
option_list <- list(
    make_option("--overlaps", type="character", help="CpG overlap file"),
    make_option("--peaks", type="character", help="All peaks file"),
    make_option("--output", type="character", help="Output PDF file")
)
opts <- parse_args(OptionParser(option_list=option_list))

# Read data
overlaps <- read.table(opts$overlaps, header=FALSE)
all_peaks <- read.table(opts$peaks, header=FALSE)

# Calculate basic statistics
total_peaks <- nrow(all_peaks)
peaks_with_cpg <- length(unique(overlaps$V4))
percent_with_cpg <- peaks_with_cpg / total_peaks * 100

# Analyze overlap characteristics
overlap_stats <- overlaps %>%
    group_by(V4) %>%
    summarize(
        total_overlap = sum(V13),
        peak_width = first(V3) - first(V2),
        overlap_fraction = total_overlap / peak_width,
        cpg_count = n()
    )

# Create visualizations
pdf(opts$output, width=10, height=12)

# 1. Overall overlap statistics
par(mfrow=c(2,2))
pie(c(peaks_with_cpg, total_peaks - peaks_with_cpg),
    labels=c("With CpG", "Without CpG"),
    main="Peaks Overlapping CpG Islands")

# 2. Overlap fraction distribution
ggplot(overlap_stats, aes(x=overlap_fraction)) +
    geom_histogram(bins=50) +
    labs(title="Distribution of CpG Overlap Fractions",
         x="Fraction of Peak Overlapping CpG Islands",
         y="Count")

# 3. CpG count per peak
ggplot(overlap_stats, aes(x=cpg_count)) +
    geom_histogram(bins=30) +
    labs(title="Number of CpG Islands per Peak",
         x="Number of CpG Islands",
         y="Count")

# 4. Peak width vs overlap
ggplot(overlap_stats, aes(x=peak_width, y=total_overlap)) +
    geom_point(alpha=0.5) +
    geom_smooth(method="lm") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title="Peak Width vs CpG Overlap",
         x="Peak Width (log10)",
         y="Total CpG Overlap (log10)")

dev.off()

# Generate summary statistics
write.table(
    data.frame(
        TotalPeaks = total_peaks,
        PeaksWithCpG = peaks_with_cpg,
        PercentWithCpG = percent_with_cpg,
        MedianOverlapFraction = median(overlap_stats$overlap_fraction),
        MedianCpGCount = median(overlap_stats$cpg_count)
    ),
    file=sub(".pdf$", "_stats.txt", opts$output),
    quote=FALSE, row.names=FALSE, sep="\t"
) 