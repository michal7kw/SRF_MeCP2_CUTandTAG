library(optparse)
library(ggplot2)
library(dplyr)
library(GenomicRanges)

# Parse command line arguments
option_list <- list(
    make_option("--common", type="character", help="Common peaks file"),
    make_option("--wt-specific", type="character", help="WT-specific peaks file"),
    make_option("--mut-specific", type="character", help="Mutant-specific peaks file"),
    make_option("--output", type="character", help="Output PDF file"),
    make_option("--cpg-only", type="logical", default=FALSE, 
                help="Whether analysis is for CpG-only peaks")
)
opts <- parse_args(OptionParser(option_list=option_list))

# Read and process peak data
read_peaks <- function(file) {
    read.table(file, header=FALSE) %>%
        select(V1, V2, V3, V7) %>%
        setNames(c("chr", "start", "end", "signal"))
}

common_peaks <- read_peaks(opts$common)
wt_specific <- read_peaks(opts$wt_specific)
mut_specific <- read_peaks(opts$mut_specific)

# Calculate peak width-normalized signals
analyze_signals <- function(df) {
    df %>% mutate(
        width = end - start,
        total_signal = signal,  # Keep raw signal
        signal_density = signal / width,  # Signal density
        log2_total_signal = log2(total_signal + 1),
        log2_signal_density = log2(signal_density + 1)
    )
}

common_peaks <- analyze_signals(common_peaks)
wt_specific <- analyze_signals(wt_specific)
mut_specific <- analyze_signals(mut_specific)

# Create enhanced visualizations
create_signal_plots <- function(common, wt_specific, mut_specific) {
    # Raw signal distribution
    p1 <- ggplot() +
        geom_density(data=common, aes(log2_total_signal, color="Common")) +
        geom_density(data=wt_specific, aes(log2_total_signal, color="WT-specific")) +
        geom_density(data=mut_specific, aes(log2_total_signal, color="Mutant-specific")) +
        labs(title="Total Signal Distribution",
             x="log2(Total Signal)", y="Density")

    # Signal density distribution
    p2 <- ggplot() +
        geom_density(data=common, aes(log2_signal_density, color="Common")) +
        geom_density(data=wt_specific, aes(log2_signal_density, color="WT-specific")) +
        geom_density(data=mut_specific, aes(log2_signal_density, color="Mutant-specific")) +
        labs(title="Signal Density Distribution",
             x="log2(Signal/Width)", y="Density")

    # Scatter plot of width vs signal
    p3 <- ggplot(common, aes(x=width, y=total_signal)) +
        geom_point(alpha=0.3) +
        geom_smooth(method="lm") +
        scale_x_log10() +
        scale_y_log10() +
        labs(title="Peak Width vs Total Signal",
             x="Peak Width (bp)", y="Total Signal")

    # Return list of plots
    list(total_signal=p1, signal_density=p2, width_vs_signal=p3)
}

# In the main execution:
plots <- create_signal_plots(common_peaks, wt_specific, mut_specific)

pdf(opts$output, width=12, height=12)
grid.arrange(
    plots$total_signal, plots$signal_density,
    plots$width_vs_signal, peak_numbers_plot,
    ncol=2
)
dev.off()

# Enhanced statistics output
write.table(
    data.frame(
        Category = c("Common", "WT-specific", "Mutant-specific"),
        Count = c(nrow(common_peaks), nrow(wt_specific), nrow(mut_specific)),
        MedianWidth = c(
            median(common_peaks$width),
            median(wt_specific$width),
            median(mut_specific$width)
        ),
        MedianTotalSignal = c(
            median(common_peaks$total_signal),
            median(wt_specific$total_signal),
            median(mut_specific$total_signal)
        ),
        MedianSignalDensity = c(
            median(common_peaks$signal_density),
            median(wt_specific$signal_density),
            median(mut_specific$signal_density)
        )
    ),
    file=sub(".pdf$", "_stats.txt", opts$output),
    quote=FALSE, row.names=FALSE, sep="\t"
) 