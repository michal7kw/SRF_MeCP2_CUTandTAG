#!/usr/bin/env Rscript

library(ChIPseeker)
library(GenomicFeatures)
library(rtracklayer)
library(ggplot2)
library(dplyr)

# Set working directory and paths
setwd("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/paper_finals")
OUTPUT_DIR <- "outputs"
PEAKS_DIR <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/peaks/narrow"
GENCODE_GTF <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA/gencode.vM25.annotation.gtf"

# Create output directory if it doesn't exist
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Function to process peaks and create annotation plot
analyze_peaks <- function(peak_files, sample_types, conditions) {
    # Import GTF file
    txdb <- makeTxDbFromGFF(GENCODE_GTF, format="gtf")
    
    # Process each group of peaks
    peak_annos <- list()
    for(i in seq_along(peak_files)) {
        # Read peaks
        peaks <- readPeakFile(peak_files[i])
        
        # Annotate peaks with parameters optimized for CUT&Tag
        peakAnno <- annotatePeak(peaks, TxDb=txdb,
                                tssRegion=c(-3000, 3000), # Narrower TSS region due to CUT&Tag precision
                                level="transcript", # More detailed annotation
                                verbose=FALSE)
        
        # Store annotation
        peak_annos[[paste(sample_types[i], conditions[i], sep="_")]] <- peakAnno
    }
    
    # Combine annotation statistics
    anno_stats <- lapply(peak_annos, function(x) {
        df <- as.data.frame(x@annoStat)
        data.frame(Feature=df$Feature,
                  Frequency=df$Frequency)
    })
    
    # Calculate average frequencies across samples
    all_stats <- do.call(rbind, Map(function(df, name) {
        df$Group <- name
        df
    }, anno_stats, names(anno_stats)))
    
    # Aggregate and calculate means
    mean_stats <- all_stats %>%
        group_by(Feature) %>%
        summarise(Mean_Frequency = mean(Frequency)) %>%
        arrange(desc(Mean_Frequency))
    
    # Create color palette
    # Simplified categories for clearer visualization
    mean_stats$Feature <- gsub(" \\(.*\\)", "", mean_stats$Feature)
    mean_stats$Feature <- gsub("Distal Intergenic", "Intergenic", mean_stats$Feature)
    
    # Aggregate similar categories
    mean_stats <- mean_stats %>%
        group_by(Feature) %>%
        summarise(Mean_Frequency = sum(Mean_Frequency)) %>%
        arrange(desc(Mean_Frequency))
    
    # Define colors matching the example pie chart
    colors <- c(
        "Promoter" = "#FFF7BC",
        "1st Intron" = "#E6E6E6",
        "Intron" = "#2C7FB8",
        "Exon" = "#7FCDBB",
        "Intergenic" = "#225EA8"
    )[mean_stats$Feature]
    
    # Create pie chart
    p <- ggplot(mean_stats, aes(x="", y=Mean_Frequency, fill=Feature)) +
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) +
        scale_fill_manual(values=colors) +
        theme_minimal() +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              plot.title = element_text(hjust=0.5)) +
        labs(fill="Feature")
    
    # Save plot
    ggsave(file.path(OUTPUT_DIR, "peak_annotation_pie.pdf"), p, width=8, height=6)
    ggsave(file.path(OUTPUT_DIR, "peak_annotation_pie.png"), p, width=8, height=6, dpi=300)
    
    return(mean_stats)
}

# Main execution
main <- function() {
    # Define sample files and metadata
    peak_files <- c(
        file.path(PEAKS_DIR, "NeuM2_narrow_peaks.filtered.narrowPeak"),
        file.path(PEAKS_DIR, "NeuM3_narrow_peaks.filtered.narrowPeak"),
        file.path(PEAKS_DIR, "NSCM1_narrow_peaks.filtered.narrowPeak"),
        file.path(PEAKS_DIR, "NSCM2_narrow_peaks.filtered.narrowPeak"),
        file.path(PEAKS_DIR, "NSCM3_narrow_peaks.filtered.narrowPeak"),
        file.path(PEAKS_DIR, "NeuV1_narrow_peaks.filtered.narrowPeak"),
        file.path(PEAKS_DIR, "NeuV2_narrow_peaks.filtered.narrowPeak"),
        file.path(PEAKS_DIR, "NeuV3_narrow_peaks.filtered.narrowPeak"),
        file.path(PEAKS_DIR, "NSCv1_narrow_peaks.filtered.narrowPeak"),
        file.path(PEAKS_DIR, "NSCv2_narrow_peaks.filtered.narrowPeak"),
        file.path(PEAKS_DIR, "NSCv3_narrow_peaks.filtered.narrowPeak")
    )
    
    sample_types <- c(
        rep("Neuron", 2), rep("NSC", 3),
        rep("Neuron", 3), rep("NSC", 3)
    )
    
    conditions <- c(
        rep("Endogenous", 5),
        rep("Exogenous", 6)
    )
    
    # Run analysis
    results <- analyze_peaks(peak_files, sample_types, conditions)
    print(results)
}

# Run main function
main()
