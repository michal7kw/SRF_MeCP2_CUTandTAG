# Check and install required packages if missing
required_packages <- c("TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db")

for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message(sprintf("Package %s not found. Installing...", pkg))
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        BiocManager::install(pkg)
    }
}

# Load required packages
library(tidyverse)
library(optparse)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(biomaRt)

# Parse command line arguments
option_list <- list(
    make_option("--work-dir", 
                type="character", 
                default=NULL, 
                dest="work_dir",
                help="Working directory for output files")
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

# Create output directory
dir.create("coding_sequence_analysis", showWarnings = FALSE)

# Read TSS peaks file first
tss_peaks <- read.csv("./peaks_annotation/exo_only_tss_peaks.csv")

# Function to get coding sequence lengths using local databases as fallback
get_coding_lengths <- function() {
    message("Attempting to get coding lengths...")
    
    # First try biomaRt
    tryCatch({
        message("Trying Ensembl biomaRt...")
        
        # Try different Ensembl mirrors with HTTPS
        mirrors <- c(
            "https://useast.ensembl.org",
            "https://uswest.ensembl.org",
            "https://asia.ensembl.org",
            "https://www.ensembl.org"
        )
        
        for (mirror in mirrors) {
            tryCatch({
                message(sprintf("Trying mirror: %s", mirror))
                mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = mirror)
                
                coding_lengths <- getBM(
                    attributes = c("external_gene_name", "cds_length"),
                    filters = "external_gene_name",
                    values = tss_peaks$SYMBOL,
                    mart = mouse
                )
                
                if (nrow(coding_lengths) > 0) {
                    message(sprintf("Successfully connected to %s", mirror))
                    return(coding_lengths)
                }
            }, error = function(e) {
                message(sprintf("Failed to connect to %s: %s", mirror, e$message))
            })
        }
        
        # If we get here, try local database
        stop("All Ensembl mirrors failed, trying local database...")
        
    }, error = function(e) {
        message("Falling back to local UCSC database...")
        
        # Get transcript lengths from local UCSC database
        txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
        
        # Get gene symbols - handle many:1 mappings
        gene_symbols <- suppressWarnings(
            select(org.Mm.eg.db, 
                   keys = tss_peaks$SYMBOL,
                   columns = c("ENTREZID", "SYMBOL"),
                   keytype = "SYMBOL") %>%
            # Remove duplicates by keeping the first mapping for each symbol
            distinct(SYMBOL, .keep_all = TRUE)
        )
        
        # Get CDS by transcript
        cds_by_tx <- cdsBy(txdb, by="tx", use.names=TRUE)
        
        # Calculate CDS lengths
        cds_lengths <- lapply(cds_by_tx, function(x) sum(width(x)))
        
        # Convert to data frame
        coding_lengths <- data.frame(
            entrez_id = names(cds_lengths),
            cds_length = unlist(cds_lengths)
        )
        
        # Join with gene symbols and handle multiple mappings
        coding_lengths <- coding_lengths %>%
            left_join(gene_symbols, by = c("entrez_id" = "ENTREZID")) %>%
            filter(!is.na(SYMBOL)) %>%
            group_by(SYMBOL) %>%
            # For each gene symbol, keep the longest CDS
            slice_max(cds_length, n = 1) %>%
            ungroup() %>%
            select(external_gene_name = SYMBOL, cds_length)
        
        message(sprintf("Found CDS lengths for %d genes", nrow(coding_lengths)))
        return(coding_lengths)
    })
}

# Get coding sequence lengths
coding_lengths <- get_coding_lengths()

# Clean and process the data
coding_analysis <- coding_lengths %>%
    filter(!is.na(cds_length)) %>%
    group_by(external_gene_name) %>%
    # Take the longest CDS for each gene
    slice_max(cds_length, n = 1) %>%
    ungroup()

# Calculate basic statistics
total_genes <- nrow(coding_analysis)
genes_over_4kb <- sum(coding_analysis$cds_length > 4000)
fraction_over_4kb <- genes_over_4kb / total_genes

# Create summary statistics
summary_stats <- list(
    total_genes = total_genes,
    genes_over_4kb = genes_over_4kb,
    fraction_over_4kb = fraction_over_4kb,
    mean_length = mean(coding_analysis$cds_length),
    median_length = median(coding_analysis$cds_length),
    min_length = min(coding_analysis$cds_length),
    max_length = max(coding_analysis$cds_length)
)

# Save summary statistics
writeLines(
    c(
        "Coding Sequence Length Analysis",
        "============================",
        "",
        paste("Total protein-coding genes analyzed:", summary_stats$total_genes),
        paste("Genes with CDS > 4kb:", summary_stats$genes_over_4kb),
        paste("Fraction of genes > 4kb:", round(summary_stats$fraction_over_4kb * 100, 2), "%"),
        "",
        "Length Statistics (bp):",
        paste("Mean length:", round(summary_stats$mean_length, 2)),
        paste("Median length:", round(summary_stats$median_length, 2)),
        paste("Minimum length:", summary_stats$min_length),
        paste("Maximum length:", summary_stats$max_length)
    ),
    "coding_sequence_analysis/summary_statistics.txt"
)

# Create detailed output of genes > 4kb
genes_over_4kb_details <- coding_analysis %>%
    filter(cds_length > 4000) %>%
    arrange(desc(cds_length))

# Save detailed gene list
write_csv(genes_over_4kb_details, 
          "coding_sequence_analysis/genes_over_4kb.csv")

# Create visualization of length distribution
p1 <- ggplot(coding_analysis, aes(x = cds_length)) +
    geom_histogram(bins = 100, fill = "#21908CFF", alpha = 0.7) +
    geom_vline(xintercept = 4000, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(
        title = "Distribution of Coding Sequence Lengths",
        subtitle = paste0(round(fraction_over_4kb * 100, 2), "% of genes > 4kb"),
        x = "Coding Sequence Length (bp)",
        y = "Count"
    )

ggsave("coding_sequence_analysis/cds_length_distribution.pdf", 
       p1, width = 10, height = 6)

# Create log-scaled version for better visualization
p2 <- ggplot(coding_analysis, aes(x = cds_length)) +
    geom_histogram(bins = 100, fill = "#21908CFF", alpha = 0.7) +
    geom_vline(xintercept = 4000, color = "red", linetype = "dashed") +
    scale_x_log10() +
    theme_minimal() +
    labs(
        title = "Distribution of Coding Sequence Lengths (Log Scale)",
        subtitle = paste0(round(fraction_over_4kb * 100, 2), "% of genes > 4kb"),
        x = "Coding Sequence Length (bp, log scale)",
        y = "Count"
    )

ggsave("coding_sequence_analysis/cds_length_distribution_log.pdf", 
       p2, width = 10, height = 6)

# Create cumulative distribution plot
p3 <- coding_analysis %>%
    arrange(cds_length) %>%
    mutate(cumulative_fraction = row_number() / n()) %>%
    ggplot(aes(x = cds_length, y = cumulative_fraction)) +
    geom_line(color = "#21908CFF") +
    geom_vline(xintercept = 4000, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(
        title = "Cumulative Distribution of Coding Sequence Lengths",
        x = "Coding Sequence Length (bp)",
        y = "Cumulative Fraction"
    )

ggsave("coding_sequence_analysis/cds_length_cumulative.pdf", 
       p3, width = 10, height = 6)

message("Analysis complete. Results saved in coding_sequence_analysis directory.") 