library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)  # For gene annotations
library(clusterProfiler)
library(ggplot2)
library(gridExtra)

# Initialize txdb object
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Parse command line arguments
library(optparse)

option_list <- list(
    make_option("--work-dir", 
                type="character", 
                default=NULL, 
                dest="work_dir",
                help="Working directory path"),
    make_option("--data-file",
                type="character",
                default=NULL,
                dest="data_file",
                help="CSV file to analyze")
)

######################################### Function definitions #########################################
# Function to read peak file
read_peak_file <- function(file_path) {
    # Read CSV file
    peaks_df <- read.csv(file_path, check.names = FALSE)
    
    # Print column names and their lengths for debugging
    cat("Column lengths:\n")
    for (col in names(peaks_df)) {
        cat(sprintf("%s: %d\n", col, length(peaks_df[[col]])))
    }
    
    # Check if all columns have the same length
    lengths <- sapply(peaks_df, length)
    if (length(unique(lengths)) > 1) {
        stop("Inconsistent column lengths in input file:\n",
             paste(names(peaks_df), lengths, sep=": ", collapse="\n"))
    }
    
    # Convert to GRanges object
    gr <- tryCatch({
        GRanges(
            seqnames = peaks_df[["chr"]],
            ranges = IRanges(
                start = peaks_df[["start"]],
                end = peaks_df[["end"]]
            ),
            # Add metadata columns
            exo_signal = peaks_df[["exo_signal"]],
            endo_signal = peaks_df[["endo_signal"]],
            enrichment = peaks_df[["enrichment"]],
            pvalue = peaks_df[["pvalue"]],
            binding_type = peaks_df[["binding_type"]],
            binding_type_by_peaks = peaks_df[["binding_type_by_peaks"]],
            significant = peaks_df[["significant"]],
            cpg_score = peaks_df[["cpg_score"]],
            cpg_name = peaks_df[["cpg_name"]],
            cpg_length = peaks_df[["cpg_length"]],
            region_length = peaks_df[["region_length"]],
            region_start = peaks_df[["region_start"]],
            region_end = peaks_df[["region_end"]]
        )
    }, error = function(e) {
        # Print the first few rows of the dataframe for debugging
        print("First few rows of the input data:")
        print(head(peaks_df))
        stop("Error creating GRanges object: ", e$message)
    })
    
    return(gr)
}

# Function to perform comprehensive peak annotation
annotate_peaks_comprehensive <- function(peaks) {
    # Create peaks_annotation directory
    dir.create("peaks_annotation_combined", showWarnings = FALSE)
    
    # Annotate peaks
    peakAnno <- annotatePeak(peaks, 
                            tssRegion = c(-3000, 3000),
                            TxDb = txdb,
                            annoDb = "org.Mm.eg.db",
                            level = "transcript",
                            verbose = FALSE)
    
    # Generate plots
    pdf("peaks_annotation_combined/peak_distribution.pdf")
    plotAnnoPie(peakAnno)
    plotDistToTSS(peakAnno)
    dev.off()
    
    # Create detailed annotation dataframe
    anno_df <- as.data.frame(peakAnno@anno)
    
    # Create a GRanges object from annotation data to ensure proper matching
    anno_gr <- GRanges(
        seqnames = anno_df$seqnames,
        ranges = IRanges(start = anno_df$start, end = anno_df$end)
    )
    
    # Find overlaps between annotation and original peaks
    overlaps <- findOverlaps(anno_gr, peaks)
    
    # Add metadata columns using the overlaps
    anno_df$binding_type <- NA
    anno_df$binding_type_by_peaks <- NA
    anno_df$exo_signal <- NA
    anno_df$endo_signal <- NA
    anno_df$enrichment <- NA
    anno_df$cpg_score <- NA
    anno_df$region_start <- NA
    anno_df$region_end <- NA
    
    # Transfer metadata using overlaps
    anno_df$binding_type[queryHits(overlaps)] <- peaks$binding_type[subjectHits(overlaps)]
    anno_df$binding_type_by_peaks[queryHits(overlaps)] <- peaks$binding_type_by_peaks[subjectHits(overlaps)]
    anno_df$exo_signal[queryHits(overlaps)] <- peaks$exo_signal[subjectHits(overlaps)]
    anno_df$endo_signal[queryHits(overlaps)] <- peaks$endo_signal[subjectHits(overlaps)]
    anno_df$enrichment[queryHits(overlaps)] <- peaks$enrichment[subjectHits(overlaps)]
    anno_df$cpg_score[queryHits(overlaps)] <- peaks$cpg_score[subjectHits(overlaps)]
    anno_df$region_start[queryHits(overlaps)] <- peaks$region_start[subjectHits(overlaps)]
    anno_df$region_end[queryHits(overlaps)] <- peaks$region_end[subjectHits(overlaps)]
    
    # Save complete annotation
    write.csv(anno_df, "peaks_annotation_combined/complete_peak_annotation.csv", row.names = FALSE)
    
    # Generate summary statistics
    summary_stats <- data.frame(
        Total_Peaks = length(peaks),
        TSS_Proximal = sum(abs(anno_df$distanceToTSS) <= 3000),
        Promoter = sum(grepl("Promoter", anno_df$annotation)),
        Intron = sum(grepl("Intron", anno_df$annotation)),
        Exon = sum(grepl("Exon", anno_df$annotation)),
        Intergenic = sum(grepl("Intergenic", anno_df$annotation))
    )
    
    write.csv(summary_stats, "peaks_annotation_combined/annotation_summary.csv", row.names = FALSE)
    
    # Perform GO enrichment analysis
    genes <- unique(anno_df$geneId[!is.na(anno_df$geneId)])
    ego <- enrichGO(gene = genes,
                   universe = keys(org.Mm.eg.db),
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)
    
    if (!is.null(ego) && nrow(ego) > 0) {
        pdf("peaks_annotation_combined/GO_enrichment.pdf")
        print(dotplot(ego, showCategory = 20))
        dev.off()
        
        write.csv(as.data.frame(ego), 
                 "peaks_annotation_combined/GO_enrichment_results.csv",
                 row.names = FALSE)
    }
    
    return(peakAnno)
}

#############################################################################################################
# Parse arguments
opts <- parse_args(OptionParser(option_list=option_list))

# Check if required arguments are provided
if (is.null(opts$work_dir) || !is.character(opts$work_dir)) {
    stop("Invalid or missing work directory path")
}
if (is.null(opts$data_file) || !is.character(opts$data_file)) {
    stop("Invalid or missing data file path")
}

# Set working directory
setwd(opts$work_dir)

# Read peak file
file_path <- file.path(opts$data_file)
message("Reading file: ", file_path)
peaks <- read_peak_file(file_path)

# Perform comprehensive peak annotation
peaks_anno <- annotate_peaks_comprehensive(peaks)

# Create summary of peaks by binding type
binding_type_summary <- data.frame(
    binding_type = mcols(peaks)$binding_type,
    enrichment = mcols(peaks)$enrichment,
    cpg_score = mcols(peaks)$cpg_score
) %>%
    group_by(binding_type) %>%
    summarise(
        count = n(),
        avg_enrichment = mean(enrichment),
        avg_cpg_score = mean(cpg_score),
        tss_proximal = sum(abs(peaks_anno@anno$distanceToTSS[binding_type == .data$binding_type]) <= 3000)
    )

write.csv(binding_type_summary, 
         "peaks_annotation_combined/binding_type_summary.csv",
         row.names = FALSE)

# Extract TSS-proximal peaks for each binding type
anno_df <- as.data.frame(peaks_anno@anno)

# Create a GRanges object from annotation data
anno_gr <- GRanges(
    seqnames = anno_df$seqnames,
    ranges = IRanges(start = anno_df$start, end = anno_df$end)
)

# Find overlaps between annotation and original peaks
overlaps <- findOverlaps(anno_gr, peaks)

# Add metadata columns
anno_df$binding_type <- NA
anno_df$binding_type_by_peaks <- NA
anno_df$exo_signal <- NA
anno_df$endo_signal <- NA
anno_df$enrichment <- NA
anno_df$cpg_score <- NA
anno_df$region_start <- NA
anno_df$region_end <- NA

# Transfer metadata using overlaps
anno_df$binding_type[queryHits(overlaps)] <- peaks$binding_type[subjectHits(overlaps)]
anno_df$binding_type_by_peaks[queryHits(overlaps)] <- peaks$binding_type_by_peaks[subjectHits(overlaps)]
anno_df$exo_signal[queryHits(overlaps)] <- peaks$exo_signal[subjectHits(overlaps)]
anno_df$endo_signal[queryHits(overlaps)] <- peaks$endo_signal[subjectHits(overlaps)]
anno_df$enrichment[queryHits(overlaps)] <- peaks$enrichment[subjectHits(overlaps)]
anno_df$cpg_score[queryHits(overlaps)] <- peaks$cpg_score[subjectHits(overlaps)]
anno_df$region_start[queryHits(overlaps)] <- peaks$region_start[subjectHits(overlaps)]
anno_df$region_end[queryHits(overlaps)] <- peaks$region_end[subjectHits(overlaps)]

# Function to get TSS peaks
get_tss_peaks <- function(anno_df, binding_type = NULL) {
    df <- anno_df
    if (!is.null(binding_type)) {
        df <- df[df$binding_type == binding_type, ]
    }
    
    tss_peaks <- df %>%
        filter(abs(distanceToTSS) <= 3000) %>%
        mutate(tss_position = ifelse(strand == "+", 
                                   start - distanceToTSS, 
                                   end + distanceToTSS)) %>%
        select(
            chr = seqnames,
            start,
            end,
            exo_signal,
            endo_signal,
            enrichment,
            pvalue,
            binding_type,
            binding_type_by_peaks,
            significant,
            cpg_score,
            cpg_name,
            cpg_length,
            region_start,
            region_end,
            region_length,
            distanceToTSS,
            tss_position,
            strand,
            geneId,
            SYMBOL
        )
    return(tss_peaks)
}

# Get TSS peaks for all peaks combined
all_tss_peaks <- get_tss_peaks(anno_df)
write.csv(all_tss_peaks, 
         "peaks_annotation_combined/all_tss_peaks.csv", 
         row.names = FALSE, 
         quote = FALSE)

# Get TSS peaks for each binding type
binding_types <- unique(anno_df$binding_type)
for (bt in binding_types) {
    tss_peaks <- get_tss_peaks(anno_df, bt)
    safe_name <- gsub(" ", "_", bt)  # Replace spaces with underscores
    write.csv(tss_peaks, 
             sprintf("peaks_annotation_combined/%s_tss_peaks.csv", safe_name),
             row.names = FALSE,
             quote = FALSE)
}

# Create summary of TSS peaks by binding type
tss_summary <- data.frame(
    binding_type = binding_types,
    total_peaks = sapply(binding_types, function(bt) sum(anno_df$binding_type == bt)),
    tss_peaks = sapply(binding_types, function(bt) sum(anno_df$binding_type == bt & abs(anno_df$distanceToTSS) <= 3000))
) %>%
    mutate(tss_percentage = (tss_peaks / total_peaks) * 100)

write.csv(tss_summary,
         "peaks_annotation_combined/tss_peaks_summary.csv",
         row.names = FALSE,
         quote = FALSE) 