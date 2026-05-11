#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(ggplot2)
  library(GenomicFeatures)
  library(rtracklayer)
})

# Set up logging
log_info <- function(msg) {
  cat(sprintf("[INFO] %s: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

# Main function to process peaks
process_peaks <- function(peaks_file, output_dir) {
  log_info("Starting peak annotation process")
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load peaks and convert is_higher_than_bg to logical
  log_info("Loading peaks data")
  peaks <- read.csv(peaks_file)
  peaks$is_higher_than_bg <- as.logical(peaks$is_higher_than_bg)
  
  # Convert to GRanges object
  peaks_gr <- GRanges(
    seqnames = peaks$chromosome,
    ranges = IRanges(start = peaks$tss - 500, end = peaks$tss + 500),
    strand = peaks$strand,
    gene = peaks$gene,
    fold_change = peaks$fold_change,
    bm3_signal = peaks$bm3_signal,
    avg_bg_signal = peaks$avg_bg_signal,
    log2FoldChange = peaks$log2FoldChange,
    padj = peaks$padj,
    is_higher_than_bg = peaks$is_higher_than_bg
  )
  
  # Load annotation database
  log_info("Loading annotation database")
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  
  # Annotate peaks
  log_info("Annotating peaks")
  peakAnno <- annotatePeak(peaks_gr, 
                          TxDb = txdb,
                          annoDb = "org.Mm.eg.db",
                          tssRegion = c(-3000, 3000))
  
  # Save detailed annotation results
  log_info("Saving annotation results")
  anno_df <- as.data.frame(peakAnno)
  anno_df$is_higher_than_bg <- as.logical(anno_df$is_higher_than_bg)
  write.csv(anno_df, 
            file = file.path(output_dir, "peak_annotations.csv"), 
            row.names = FALSE)
  
  # Separate and save gene lists based on background comparison
  log_info("Separating genes based on background signal comparison")
  
  # Create separate dataframes for higher and lower than background
  higher_than_bg <- anno_df[anno_df$is_higher_than_bg == TRUE, ]
  lower_than_bg <- anno_df[anno_df$is_higher_than_bg == FALSE, ]
  
  # Save separated gene lists
  write.csv(higher_than_bg,
            file = file.path(output_dir, "peaks_higher_than_bg.csv"),
            row.names = FALSE)
  write.csv(lower_than_bg,
            file = file.path(output_dir, "peaks_lower_than_bg.csv"),
            row.names = FALSE)
  
  # Update summary statistics
  summary_stats <- data.frame(
    total_peaks = length(peaks_gr),
    promoter_peaks = sum(abs(anno_df$distanceToTSS) <= 1000),
    distal_peaks = sum(abs(anno_df$distanceToTSS) > 1000),
    avg_distance_to_tss = mean(abs(anno_df$distanceToTSS)),
    correlation = cor(anno_df$log2FoldChange, anno_df$fold_change, 
                     use = "complete.obs"),
    higher_than_bg_count = sum(anno_df$is_higher_than_bg),
    lower_than_bg_count = sum(!anno_df$is_higher_than_bg)
  )
  
  # Generate peak annotation pie plot
  log_info("Generating peak annotation pie plot")
  pdf(file.path(output_dir, "peak_annotation_pie.pdf"))
  plotAnnoPie(peakAnno)
  dev.off()
  
  write.csv(summary_stats,
            file = file.path(output_dir, "summary_statistics.csv"),
            row.names = FALSE)
  
  log_info("Analysis completed successfully")
}

# Main execution
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  base_dir <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Cross_final"
  
  # Input/output paths
  peaks_file <- file.path(base_dir, "results/9analyze_smarcb1_dea_genes/smarcb1_tss_analysis.csv")
  output_dir <- file.path(base_dir, "results/10annotate_smarcb1_peaks")
  
  # Run analysis
  process_peaks(peaks_file, output_dir)
}

# Run main function
if (!interactive()) {
  main()
}
