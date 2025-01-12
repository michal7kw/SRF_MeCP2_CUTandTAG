#!/usr/bin/env Rscript

# Function to install Bioconductor packages
install_bioc_packages <- function(packages) {
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    # Install all required Bioconductor packages
    BiocManager::install(packages, update = FALSE)
}

# Function to install CRAN packages
install_cran_packages <- function(packages) {
    new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
    if (length(new_packages) > 0) {
        install.packages(new_packages)
    }
}

# First, install required CRAN packages
cran_packages <- c(
    "optparse",
    "ggplot2",
    "dplyr"
)
install_cran_packages(cran_packages)

# Then install Bioconductor packages in the correct order
bioc_packages <- c(
    "GenomicRanges",
    "rtracklayer",
    "GenomicFeatures",
    "TxDb.Mmusculus.UCSC.mm10.knownGene",
    "org.Mm.eg.db",
    "BSgenome.Mmusculus.UCSC.mm10",
    "TxDb.Hsapiens.UCSC.hg19.knownGene",
    "ChIPseeker"
)
install_bioc_packages(bioc_packages)

# Load packages with suppressed messages
suppressPackageStartupMessages({
    # Load base packages first
    library(GenomicRanges)
    library(rtracklayer)
    library(GenomicFeatures)
    
    # Load annotation packages
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    library(BSgenome.Mmusculus.UCSC.mm10)
    
    # Load ChIPseeker after its dependencies
    library(ChIPseeker)
    
    # Load other packages
    library(optparse)
    library(ggplot2)
    library(dplyr)
})

# Parse command line arguments (simplified)
option_list <- list(
    make_option("--peaks-dir", type="character", default=NULL, dest="peaks_dir",
                help="Directory containing narrowPeak files"),
    make_option("--output-dir", type="character", default=NULL, dest="output_dir",
                help="Output directory for plots")
)
opts <- parse_args(OptionParser(option_list=option_list))

# Add diagnostic prints
print("Command line arguments received:")
print(opts)
print("Output directory path:")
print(opts$output_dir)

# Ensure the output directory path is not NULL and is a character string
if (is.null(opts$output_dir) || !is.character(opts$output_dir)) {
    stop("Invalid or missing output directory path")
}

# Create output directory with error handling
tryCatch({
    dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
}, error = function(e) {
    print("Error creating directory:")
    print(e)
    stop("Failed to create output directory")
})

# Get CpG islands from UCSC
get_cpg_islands <- function() {
    session <- browserSession("UCSC")
    genome(session) <- "mm10"
    query <- ucscTableQuery(session, table="cpgIslandExt")
    cpg <- track(query)
    
    # Convert to GRanges and ensure proper chromosome naming
    cpg_gr <- GRanges(
        seqnames = seqnames(cpg),
        ranges = ranges(cpg),
        strand = strand(cpg)
    )
    
    # Filter out non-standard chromosomes
    standard_chroms <- paste0("chr", c(1:19, "X", "Y", "M"))
    cpg_gr <- cpg_gr[seqnames(cpg_gr) %in% standard_chroms]
    
    return(cpg_gr)
}

# Load and process peaks
load_peaks <- function(peak_file) {
    # Read the narrowPeak file with more flexible column specifications
    peaks_df <- read.table(peak_file, 
                          col.names = c("chrom", "start", "end", "name", "score",
                                      "strand", "signalValue", "pValue", "qValue", "peak"),
                          colClasses = c("character", "numeric", "numeric", "character",
                                       "numeric", "character", "numeric", "numeric",
                                       "numeric", "numeric"))
    
    # Fix strand values - replace empty or invalid strands with '*'
    peaks_df$strand[!peaks_df$strand %in% c("+", "-")] <- "*"
    
    # Convert to GRanges object
    peaks_gr <- GRanges(
        seqnames = peaks_df$chrom,
        ranges = IRanges(start = peaks_df$start, end = peaks_df$end),
        strand = peaks_df$strand,
        score = peaks_df$score,
        signalValue = peaks_df$signalValue,
        pValue = peaks_df$pValue,
        qValue = peaks_df$qValue,
        peak = peaks_df$peak
    )
    
    return(peaks_gr)
}

# Analyze peak distribution
analyze_peak_distribution <- function(peaks, cpg_islands, txdb) {
    # Create promoter regions (2kb upstream)
    promoter <- promoters(txdb, upstream=2000, downstream=0)
    
    # Filter peaks to include only standard chromosomes
    standard_chroms <- paste0("chr", c(1:19, "X", "Y", "M"))
    peaks <- peaks[seqnames(peaks) %in% standard_chroms]
    
    # Create genomic state priority list
    genomic_regions <- list(
        promoter_cpg = subsetByOverlaps(promoter, cpg_islands),
        other_cpg = cpg_islands,
        promoter = promoter,
        genes = genes(txdb)
    )
    
    # Initialize counters
    categories <- c(
        "Promoter CpG Island" = 0,
        "Other CpG Island" = 0,
        "Promoter non-CpG" = 0,
        "Gene body" = 0,
        "Intergenic" = 0
    )
    
    # Categorize each peak
    peak_categories <- rep("Intergenic", length(peaks))
    
    # Check overlaps in priority order
    overlaps_promoter_cpg <- overlapsAny(peaks, genomic_regions$promoter_cpg)
    peak_categories[overlaps_promoter_cpg] <- "Promoter CpG Island"
    
    remaining_peaks <- !overlaps_promoter_cpg
    overlaps_other_cpg <- overlapsAny(peaks[remaining_peaks], genomic_regions$other_cpg)
    peak_categories[remaining_peaks][overlaps_other_cpg] <- "Other CpG Island"
    
    remaining_peaks <- remaining_peaks & !overlaps_other_cpg
    overlaps_promoter <- overlapsAny(peaks[remaining_peaks], genomic_regions$promoter)
    peak_categories[remaining_peaks][overlaps_promoter] <- "Promoter non-CpG"
    
    remaining_peaks <- remaining_peaks & !overlaps_promoter
    overlaps_genes <- overlapsAny(peaks[remaining_peaks], genomic_regions$genes)
    peak_categories[remaining_peaks][overlaps_genes] <- "Gene body"
    
    # Count categories
    for (cat in names(categories)) {
        categories[cat] <- sum(peak_categories == cat)
    }
    
    return(categories)
}

# Plot distribution (unchanged)
plot_distributions <- function(endo_data, exo_data, output_dir) {
    # Prepare data for plotting
    plot_data <- data.frame(
        Category = rep(names(endo_data), 2),
        Count = c(endo_data, exo_data),
        Type = rep(c("Endogenous", "Exogenous"), each=length(endo_data))
    )
    
    # Calculate percentages
    plot_data <- plot_data %>%
        group_by(Type) %>%
        mutate(Percentage = Count / sum(Count) * 100)
    
    # Define colors
    colors <- c(
        "Promoter CpG Island" = "#e5e5b5",
        "Other CpG Island" = "#b5d7e5",
        "Promoter non-CpG" = "#3182bd",
        "Gene body" = "#fd8d3c",
        "Intergenic" = "#e7969c"
    )
    
    # Create plots
    p1 <- ggplot(plot_data, aes(x="", y=Percentage, fill=Category)) +
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) +
        facet_wrap(~Type) +
        scale_fill_manual(values=colors) +
        theme_minimal() +
        labs(title="Peak Distribution in NPCs",
             fill="Genomic Region") +
        theme(axis.text.x=element_blank(),
              axis.ticks=element_blank(),
              panel.grid=element_blank())
    
    # Save plots
    ggsave(file.path(output_dir, "peak_distributions_NPCs_combined_cpg_specific.pdf"),
           p1, width=12, height=6)
    
    # Also save as PNG for easy viewing
    ggsave(file.path(output_dir, "peak_distributions_NPCs_combined_cpg_specific.png"),
           p1, width=12, height=6, dpi=300)
    
    # Print statistics
    for (type in unique(plot_data$Type)) {
        cat(sprintf("\nDistribution for %s:\n", type))
        type_data <- plot_data %>% filter(Type == type)
        for (i in 1:nrow(type_data)) {
            cat(sprintf("%s: %d peaks (%.1f%%)\n",
                       type_data$Category[i],
                       type_data$Count[i],
                       type_data$Percentage[i]))
        }
    }
}

# Main execution
main <- function() {
    message("Loading mm10 annotations...")
    # Load annotations directly from Bioconductor packages
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    
    message("Fetching CpG islands from UCSC...")
    cpg_islands <- get_cpg_islands()
    
    message("Processing peak files...")
    # Process peak files
    peak_files <- list(
        endo = file.path(opts$peaks_dir, "NPCs_endo_combined.narrowPeak"),
        exo = file.path(opts$peaks_dir, "NPCs_exo_combined.narrowPeak")
    )
    
    results <- list()
    for (condition in names(peak_files)) {
        if (file.exists(peak_files[[condition]])) {
            message(sprintf("Processing %s peaks...", condition))
            peaks <- load_peaks(peak_files[[condition]])
            results[[condition]] <- analyze_peak_distribution(peaks, cpg_islands, txdb)
        } else {
            warning(sprintf("File not found: %s", peak_files[[condition]]))
        }
    }
    
    if (length(results) == 2) {
        message("Creating plots...")
        plot_distributions(results$endo, results$exo, opts$output_dir)
    }
}

# Run main function
main() 