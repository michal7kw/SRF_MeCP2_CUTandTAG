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
  library(GenomicRanges)
})

# Set up logging
log_info <- function(msg) {
  cat(sprintf("[INFO] %s: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

# Function to read narrowPeak files
read_narrow_peaks <- function(file_path) {
  peaks <- read.table(file_path, 
                     col.names = c("chr", "start", "end", "name", "score", 
                                 "strand", "signalValue", "pValue", 
                                 "qValue", "peak"))
  GRanges(
    seqnames = peaks$chr,
    ranges = IRanges(start = peaks$start, end = peaks$end),
    strand = "*",
    score = peaks$score
  )
}

# Function to check if gene has peaks nearby
find_genes_with_peaks <- function(tss_data, peak_files, window_size = 1000) {
  log_info("Finding genes with peaks in their vicinity")
  
  # Create GRanges for TSS regions
  tss_ranges <- GRanges(
    seqnames = tss_data$chromosome,
    ranges = IRanges(start = tss_data$tss - window_size,
                    end = tss_data$tss + window_size),
    strand = tss_data$strand,
    gene = tss_data$gene
  )
  
  # Read and process each peak file
  has_peaks <- rep(FALSE, length(tss_ranges))
  for (peak_file in peak_files) {
    log_info(sprintf("Processing peak file: %s", basename(peak_file)))
    peaks <- read_narrow_peaks(peak_file)
    
    # Find overlaps
    overlaps <- findOverlaps(tss_ranges, peaks)
    has_peaks[unique(queryHits(overlaps))] <- TRUE
  }
  
  # Return genes that have peaks
  return(tss_data$gene[has_peaks])
}

# Main function to process peaks
process_peaks <- function(peaks_file, peak_files, output_dir) {
  log_info("Starting peak annotation process")
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load peaks and convert is_higher_than_bg to logical
  log_info("Loading peaks data")
  peaks <- read.csv(peaks_file)
  
  # Find genes with peaks nearby
  genes_with_peaks <- find_genes_with_peaks(peaks, peak_files)
  log_info(sprintf("Found %d genes with peaks nearby", length(genes_with_peaks)))
  
  # Filter peaks data
  peaks_filtered <- peaks[peaks$gene %in% genes_with_peaks, ]
  peaks_filtered$is_higher_than_bg <- as.logical(peaks_filtered$is_higher_than_bg)
  
  # Save filtered gene list
  write.csv(peaks_filtered,
            file = file.path(output_dir, "filtered_genes.csv"),
            row.names = FALSE)
  
  # Create a unique identifier for each peak
  peak_ids <- paste0("peak_", seq_len(nrow(peaks_filtered)))
  
  # Convert to GRanges object with unique names
  peaks_gr <- GRanges(
    seqnames = peaks_filtered$chromosome,
    ranges = IRanges(start = peaks_filtered$tss - 500, 
                    end = peaks_filtered$tss + 500),
    strand = peaks_filtered$strand,
    gene = peaks_filtered$gene,
    fold_change = peaks_filtered$fold_change,
    bm3_signal = peaks_filtered$bm3_signal,
    avg_bg_signal = peaks_filtered$avg_bg_signal,
    log2FoldChange = peaks_filtered$log2FoldChange,
    padj = peaks_filtered$padj,
    is_higher_than_bg = peaks_filtered$is_higher_than_bg,
    peak_id = peak_ids  # Add unique identifier
  )
  names(peaks_gr) <- peak_ids
  
  # Load annotation database
  log_info("Loading annotation database")
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  
  # Annotate peaks
  log_info("Annotating peaks")
  peakAnno <- suppressWarnings(
    annotatePeak(peaks_gr, 
                 TxDb = txdb,
                 annoDb = "org.Mm.eg.db",
                 tssRegion = c(-3000, 3000))
  )
  
  # Convert to data frame and ensure proper column types
  log_info("Processing annotation results")
  anno_df <- as.data.frame(peakAnno)
  
  # Match back to original peaks using peak_id
  anno_df$is_higher_than_bg <- peaks_gr$is_higher_than_bg[match(anno_df$peak, names(peaks_gr))]
  
  # Save detailed annotation results
  write.csv(anno_df, 
            file = file.path(output_dir, "peak_annotations.csv"), 
            row.names = FALSE)
  
  # Separate and save gene lists based on background comparison
  log_info("Separating genes based on background signal comparison")
  higher_than_bg <- anno_df[anno_df$is_higher_than_bg, ]
  lower_than_bg <- anno_df[!anno_df$is_higher_than_bg, ]
  
  write.csv(higher_than_bg,
            file = file.path(output_dir, "peaks_higher_than_bg.csv"),
            row.names = FALSE)
  write.csv(lower_than_bg,
            file = file.path(output_dir, "peaks_lower_than_bg.csv"),
            row.names = FALSE)
  
  # Generate plots and statistics
  generate_plots(anno_df, peakAnno, txdb, output_dir, peaks_gr)
  generate_summary_stats(peaks_gr, anno_df, output_dir)
  
  log_info("Analysis completed successfully")
}

# Function to generate plots
generate_plots <- function(anno_df, peakAnno, txdb, output_dir, peaks_gr) {
  log_info("Generating visualization plots")
  
  # Standard ChIPseeker plots
  log_info("Generating ChIPseeker plots")
  
  # Peak annotation pie chart
  pdf(file.path(output_dir, "peak_annotation_pie.pdf"))
  tryCatch({
    plotAnnoPie(peakAnno)
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    log_info(sprintf("Warning: Could not generate pie chart: %s", e$message))
  })
  
  # Peak distribution
  tryCatch({
    # Create the plot object first
    p_dist <- plotDistToTSS(peakAnno,
                           title = "Distribution of peaks relative to TSS",
                           ylab = "Peak Count")
    
    # Save to PDF with specific dimensions
    pdf(file.path(output_dir, "peak_distribution.pdf"), 
        width = 8, 
        height = 6)
    print(p_dist)
    dev.off()
    
    log_info("Successfully generated peak distribution plot")
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    log_info(sprintf("Warning: Could not generate peak distribution plot: %s", e$message))
  })
  
  # Average profile plot
  tryCatch({
    # Check if peakAnno object contains valid data
    if (!is.null(peakAnno) && length(peakAnno) > 0) {
      output_file <- file.path(output_dir, "average_profile.pdf")
      
      # Create tagMatrix first
      tagMatrix <- getTagMatrix(peaks_gr, 
                              TxDb = txdb,
                              upstream = 3000, 
                              downstream = 3000,
                              type = "start")
      
      if (!is.null(tagMatrix) && ncol(tagMatrix) > 0 && nrow(tagMatrix) > 0) {
        # Open PDF device
        pdf(output_file, width = 10, height = 7)
        
        # Plot using tagMatrix
        plotAvgProf(tagMatrix,
                   xlim = c(-3000, 3000),
                   xlab = "Distance from TSS (bp)",
                   ylab = "Read Count Frequency",
                   main = "Average Peak Profile",
                   by = "gene",
                   conf = 0.95)
        
        # Explicitly close the device
        dev.off()
        
        # Verify the file was created and has content
        if (file.exists(output_file) && file.size(output_file) > 0) {
          log_info("Successfully generated average profile plot")
        } else {
          log_info("Warning: Average profile plot file was not created properly")
        }
      } else {
        log_info("Warning: Could not generate tag matrix for profile plot")
      }
    } else {
      log_info("Warning: No valid peak data available for profile plot")
    }
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    log_info(sprintf("Warning: Could not generate average profile plot: %s", e$message))
  })
  
  log_info("Plot generation completed")
}

# Function to generate summary statistics
generate_summary_stats <- function(peaks_gr, anno_df, output_dir) {
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
  
  write.csv(summary_stats,
            file = file.path(output_dir, "summary_statistics.csv"),
            row.names = FALSE)
}

# Main execution
main <- function() {
  base_dir <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/Cross_final"
  peaks_dir <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1/results/peaks"
  
  # Input/output paths
  peaks_file <- file.path(base_dir, "results/9analyze_smarcb1_dea_genes/smarcb1_tss_analysis.csv")
  output_dir <- file.path(base_dir, "results/11annotate_genes_filtered_from_9")
  
  # Get peak files
  peak_files <- list.files(peaks_dir, pattern = "*.narrowPeak$", full.names = TRUE)
  
  # Run analysis
  process_peaks(peaks_file, peak_files, output_dir)
}

# Run main function
if (!interactive()) {
  main()
} 