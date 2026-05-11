#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
    library(optparse)
    library(ChIPseeker)
    library(GenomicFeatures)
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    library(rtracklayer)
    library(ggplot2)
})

# Parse command line arguments
option_list <- list(
    make_option("--peaks-dir", type="character", default=NULL, dest="peaks_dir",
                help="Directory containing narrowPeak files"),
    make_option("--output-dir", type="character", default=NULL, dest="output_dir",
                help="Output directory for plots and results"),
    make_option("--peak-type", type="character", default="narrow", dest="peak_type",
                help="Type of peaks to process: narrow or broad")
)
opts <- parse_args(OptionParser(option_list=option_list))

# Add diagnostic prints
print("Command line arguments received:")
print(opts)
print("Output directory path:")
print(opts$output_dir)
print("Peak type:")
print(opts$peak_type)

# Ensure the peak type is valid
if (!opts$peak_type %in% c("narrow", "broad")) {
    stop("Invalid peak type. Must be 'narrow' or 'broad'.")
}

# Ensure the output directory path is not NULL and is a character string
if (is.null(opts$output_dir) || !is.character(opts$output_dir)) {
    stop("Invalid or missing output directory path")
}

# Create output directory
dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

# Main function
main <- function() {
    # Load TxDb
    message("Loading mm10 annotations...")
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    
    # Process peak files
    peak_files <- list(
        endo = file.path(opts$peaks_dir, paste0("NPCs_endo_combined.", opts$peak_type, "Peak")),
        exo = file.path(opts$peaks_dir, paste0("NPCs_exo_combined.", opts$peak_type, "Peak"))
    )
    
    # Store annotations
    annotations <- list()
    
    # Process each peak file
    for (condition in names(peak_files)) {
        if (file.exists(peak_files[[condition]])) {
            message(sprintf("Processing %s peaks...", condition))
            
            # Read the peaks file manually
            peaks_df <- read.table(peak_files[[condition]], 
                                 col.names=c("chrom", "start", "end", "name", 
                                           "score", "strand", "signalValue",
                                           "pValue", "qValue", "peak"),
                                 colClasses=c("character", "numeric", "numeric", 
                                            "character", "numeric", "character",
                                            "numeric", "numeric", "numeric", "numeric"))
            
            # Convert "." in strand to "*"
            peaks_df$strand[peaks_df$strand == "."] <- "*"
            
            # Convert to GRanges
            peaks <- GRanges(seqnames=peaks_df$chrom,
                           ranges=IRanges(start=peaks_df$start, end=peaks_df$end),
                           name=peaks_df$name,
                           score=peaks_df$score,
                           strand=peaks_df$strand,
                           signalValue=peaks_df$signalValue,
                           pValue=peaks_df$pValue,
                           qValue=peaks_df$qValue,
                           peak=peaks_df$peak)
            
            # Annotate peaks
            annotations[[condition]] <- annotatePeak(peaks,
                                                   TxDb = txdb,
                                                   annoDb = "org.Mm.eg.db",
                                                   tssRegion = c(-3000, 3000),
                                                   level = "transcript")
            
            # Generate individual plots
            pdf(file.path(opts$output_dir, sprintf("%s_annotation_plots.pdf", condition)), 
                width = 10, height = 12)
            
            # Plot genomic annotation
            print(plotAnnoBar(annotations[[condition]]))
            
            # Plot distance to TSS
            print(plotDistToTSS(annotations[[condition]]))
            
            # Plot annotation pie chart
            print(plotAnnoPie(annotations[[condition]]))
            
            dev.off()
            
            # Save detailed annotation results
            anno_df <- as.data.frame(annotations[[condition]]@anno)
            write.table(anno_df,
                       file = file.path(opts$output_dir, 
                                      sprintf("%s_peak_annotations.tsv", condition)),
                       sep = "\t", quote = FALSE, row.names = FALSE)
        } else {
            warning(sprintf("File not found: %s", peak_files[[condition]]))
        }
    }
    
    # Create combined visualization if both conditions are present
    if (length(annotations) == 2) {
        pdf(file.path(opts$output_dir, "combined_annotation_plots.pdf"), 
            width = 12, height = 8)
        
        # Combined bar plot
        print(plotAnnoBar(annotations))
        
        # Combined distance to TSS
        print(plotDistToTSS(annotations))
        
        dev.off()
    }
}

# Run main function
main() 