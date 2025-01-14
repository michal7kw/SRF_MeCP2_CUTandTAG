library(tidyverse)
library(optparse)
library(GenomicRanges)
library(org.Mm.eg.db)
library(clusterProfiler)

# Add package installation checks and installation if needed
required_packages <- c("VennDiagram", "tidyverse", "optparse", "GenomicRanges", 
                      "org.Mm.eg.db", "clusterProfiler")

for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
        message(sprintf("Installing package: %s", package))
        if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install(package, update = FALSE)
    }
    # Load package with error handling
    tryCatch({
        library(package, character.only = TRUE)
    }, error = function(e) {
        message(sprintf("Warning: Failed to load package %s: %s", package, e$message))
    })
}

# Add error handling for biomaRt connections
options(timeout = 300)  # Increase timeout for downloads
options(connectionObserver = NULL)  # Disable connection observer which can cause issues

# Parse command line arguments
option_list <- list(
    make_option("--work-dir", 
                type="character", 
                default=NULL, 
                dest="work_dir",
                help="Working directory containing the analysis results")
)

# Parse arguments
opts <- parse_args(OptionParser(option_list=option_list))

# Add diagnostic prints
print("Command line arguments received:")
print(opts)
print("Working directory path:")
print(opts$work_dir)

# Check if work-dir argument is provided
if (is.null(opts$work_dir) || !is.character(opts$work_dir)) {
    stop("Invalid or missing work directory path")
}

# Set working directory
message(sprintf("Setting working directory to: %s", opts$work_dir))
setwd(opts$work_dir)

# Create output directory for gene analysis
dir.create("gene_analysis", showWarnings = FALSE)

# Function to read and process gene data
process_gene_data <- function(file_path, category) {
    data <- read_csv(file_path) %>%
        mutate(category = category) %>%
        # Rename seqnames to chr for consistency
        rename(chr = seqnames)
    
    # Convert ENTREZ IDs to gene symbols
    genes_with_symbols <- data %>%
        filter(!is.na(geneId)) %>%
        mutate(
            geneId = as.character(geneId),
            category = category
        )
    
    return(genes_with_symbols)
}

# Read the TSS-proximal peaks data
endo_genes <- process_gene_data("peaks_annotation/endo_only_tss_peaks.csv", "Endo Only")
exo_genes <- process_gene_data("peaks_annotation/exo_only_tss_peaks.csv", "Exo Only")
both_genes <- process_gene_data("peaks_annotation/both_tss_peaks.csv", "Both")

# Combine all gene data
all_genes <- bind_rows(endo_genes, exo_genes, both_genes)

# Calculate basic statistics
gene_stats <- all_genes %>%
    group_by(category) %>%
    summarise(
        unique_genes = n_distinct(geneId),
        mean_distance = mean(abs(distanceToTSS)),
        median_distance = median(abs(distanceToTSS))
    )

# Save gene statistics
write_csv(gene_stats, "gene_analysis/gene_statistics.csv")

# Analyze genes targeted by both Exo and Endo
both_detailed <- both_genes %>%
    arrange(abs(distanceToTSS))

# Save detailed analysis of dual-targeted genes
write_csv(both_detailed, "gene_analysis/both_targeted_genes_detailed.csv")

# Create Venn diagram data
venn_data <- list(
    Endo = unique(endo_genes$geneId),
    Exo = unique(exo_genes$geneId),
    Both = unique(both_genes$geneId)
)

# Modify the Venn diagram section to include error handling
tryCatch({
    venn.plot <- venn.diagram(
        venn_data,
        filename = "gene_analysis/gene_overlaps_venn.png",
        imagetype = "png",
        height = 3000,
        width = 3000,
        resolution = 300,
        compression = "lzw",
        lwd = 2,
        col = c("#440154FF", "#21908CFF", "#FDE725FF"),
        fill = c(alpha("#440154FF", 0.3), 
                alpha("#21908CFF", 0.3), 
                alpha("#FDE725FF", 0.3)),
        main = "Gene Overlaps Between Categories"
    )
}, error = function(e) {
    warning("Failed to create Venn diagram: ", e$message)
})

# Calculate enrichment scores for dual-targeted regions
calculate_enrichment_scores <- function(both_data) {
    # Read the original peak files to get signal values
    both_peaks <- read_csv("cpg_enrichment_both.csv")
    
    enrichment_scores <- both_peaks %>%
        select(chr, start, end, exo_signal, endo_signal, enrichment) %>%
        mutate(
            log2_enrichment = log2(enrichment),
            normalized_enrichment = (exo_signal - endo_signal) / (exo_signal + endo_signal)
        )
    
    # Add diagnostic print to check column names
    message("Enrichment scores columns: ", paste(colnames(enrichment_scores), collapse=", "))
    
    return(enrichment_scores)
}

# Before the join, add diagnostic prints
message("Genes columns: ", paste(colnames(both_genes), collapse=", "))

enrichment_scores <- calculate_enrichment_scores(both_genes)

# Combine enrichment scores with gene information
genes_with_enrichment <- both_genes %>%
    left_join(
        enrichment_scores,
        by = c("chr", "start", "end")
    ) %>%
    arrange(desc(normalized_enrichment))

# Save detailed enrichment analysis
write_csv(genes_with_enrichment, "gene_analysis/genes_with_enrichment_scores.csv")

# Create distance vs enrichment plot with improved smoothing
ggplot(genes_with_enrichment, 
       aes(x = abs(distanceToTSS), y = normalized_enrichment)) +
    geom_point(alpha = 0.5, color = "#21908CFF") +
    geom_smooth(method = "gam",  # Changed from loess to gam
                formula = y ~ s(x, bs = "cs"),  # Using cubic spline smoothing
                color = "#FDE725FF") +
    theme_minimal() +
    scale_x_continuous(trans = "log1p") +  # Log transform x-axis to handle wide range
    labs(
        title = "Enrichment Score vs Distance to TSS",
        x = "Distance to TSS (bp)",
        y = "Normalized Enrichment Score"
    )
ggsave("gene_analysis/distance_vs_enrichment.pdf", width = 10, height = 6)

# Create enrichment score visualization with improved binning
ggplot(genes_with_enrichment, aes(x = normalized_enrichment)) +
    geom_histogram(bins = 30,  # Reduced number of bins
                  fill = "#21908CFF", 
                  alpha = 0.7,
                  boundary = 0) +  # Align bins on zero
    theme_minimal() +
    labs(
        title = "Distribution of Normalized Enrichment Scores",
        x = "Normalized Enrichment Score",
        y = "Count"
    )
ggsave("gene_analysis/enrichment_score_distribution.pdf", width = 10, height = 6)

# Generate summary statistics for high-enrichment genes
high_enrichment_threshold <- quantile(genes_with_enrichment$normalized_enrichment, 0.75)
high_enrichment_genes <- genes_with_enrichment %>%
    filter(normalized_enrichment >= high_enrichment_threshold) %>%
    arrange(desc(normalized_enrichment))

# Save high-enrichment genes
write_csv(high_enrichment_genes, "gene_analysis/high_enrichment_genes.csv")

# Create summary report
report <- c(
    "Gene Analysis Summary",
    "==================",
    "",
    paste("Total unique genes analyzed:", nrow(all_genes)),
    paste("Genes with both Exo and Endo binding:", nrow(both_genes)),
    paste("High-enrichment genes:", nrow(high_enrichment_genes)),
    "",
    "Category-wise statistics:",
    capture.output(print(gene_stats)),
    "",
    "Enrichment Score Statistics:",
    paste("Median enrichment score:", median(genes_with_enrichment$normalized_enrichment)),
    paste("Mean enrichment score:", mean(genes_with_enrichment$normalized_enrichment)),
    paste("High-enrichment threshold:", high_enrichment_threshold)
)

# Save summary report
writeLines(report, "gene_analysis/analysis_summary.txt")

message("Gene analysis complete. Results saved in gene_analysis directory.") 