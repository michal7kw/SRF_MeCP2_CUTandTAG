# Check and install required packages if missing
required_packages <- c(
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Hs.eg.db",
    "tidyverse",
    "GenomicFeatures",
    "biomaRt",
    "ggplot2",
    "scales"  # For better number formatting
)

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
suppressPackageStartupMessages({
    library(tidyverse)
    library(GenomicFeatures)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(biomaRt)
    library(scales)
})

# Create output directory
output_dir <- "human_coding_sequence_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to get coding sequence lengths using multiple sources
get_human_coding_lengths <- function() {
    message("Retrieving human coding sequence lengths...")
    
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
                human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = mirror)
                
                coding_lengths <- getBM(
                    attributes = c(
                        "external_gene_name",
                        "cds_length",
                        "gene_biotype",
                        "chromosome_name"
                    ),
                    filters = "gene_biotype",
                    values = "protein_coding",
                    mart = human
                )
                
                if (nrow(coding_lengths) > 0) {
                    message(sprintf("Successfully retrieved %d genes from %s", 
                                  nrow(coding_lengths), mirror))
                    return(coding_lengths)
                }
            }, error = function(e) {
                message(sprintf("Failed to connect to %s: %s", mirror, e$message))
            })
        }
        
        stop("All Ensembl mirrors failed, trying local database...")
        
    }, error = function(e) {
        message("Falling back to local UCSC database...")
        
        # Get transcript lengths from local UCSC database
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
        
        # Get CDS ranges
        cds_by_gene <- cdsBy(txdb, by="gene")
        
        if (length(cds_by_gene) == 0) {
            stop("No CDS data found in local database")
        }
        
        # Calculate CDS lengths for each gene
        cds_lengths <- vapply(cds_by_gene, function(x) sum(width(x)), numeric(1))
        
        # Get gene symbols
        gene_ids <- names(cds_by_gene)
        
        # Map Entrez IDs to gene symbols
        gene_symbols <- suppressWarnings({
            AnnotationDbi::select(org.Hs.eg.db,
                                keys = gene_ids,
                                columns = c("SYMBOL"),
                                keytype = "ENTREZID")
        })
        
        # Create data frame
        coding_lengths <- data.frame(
            entrez_id = gene_ids,
            cds_length = cds_lengths
        ) %>%
            dplyr::left_join(gene_symbols, by = c("entrez_id" = "ENTREZID")) %>%
            dplyr::filter(!is.na(SYMBOL)) %>%
            dplyr::group_by(SYMBOL) %>%
            dplyr::slice_max(cds_length, n = 1) %>%
            dplyr::ungroup() %>%
            dplyr::rename(external_gene_name = SYMBOL)
        
        if (nrow(coding_lengths) == 0) {
            stop("Failed to retrieve any gene data from local database")
        }
        
        message(sprintf("Found CDS lengths for %d genes from local database", 
                       nrow(coding_lengths)))
        return(coding_lengths)
    })
}

# Get coding sequence lengths
coding_lengths <- get_human_coding_lengths()

# Clean and process the data
coding_analysis <- coding_lengths %>%
    filter(!is.na(cds_length)) %>%
    filter(cds_length > 0) %>%  # Remove entries with zero length
    group_by(external_gene_name) %>%
    slice_max(cds_length, n = 1) %>%  # Take the longest CDS for each gene
    ungroup()

# Calculate statistics
summary_stats <- list(
    total_genes = nrow(coding_analysis),
    genes_over_4kb = sum(coding_analysis$cds_length > 4000),
    fraction_over_4kb = sum(coding_analysis$cds_length > 4000) / nrow(coding_analysis),
    quantiles = quantile(coding_analysis$cds_length, probs = c(0, 0.25, 0.5, 0.75, 1)),
    mean_length = mean(coding_analysis$cds_length),
    sd_length = sd(coding_analysis$cds_length)
)

# Save summary statistics with improved formatting
writeLines(
    c(
        "Human Coding Sequence Length Analysis",
        "===================================",
        "",
        paste("Total protein-coding genes analyzed:", 
              format(summary_stats$total_genes, big.mark = ",")),
        paste("Genes with CDS > 4kb:", 
              format(summary_stats$genes_over_4kb, big.mark = ",")),
        paste("Fraction of genes > 4kb:", 
              percent(summary_stats$fraction_over_4kb, accuracy = 0.1)),
        "",
        "Length Statistics (bp):",
        paste("Mean length:", format(round(summary_stats$mean_length), big.mark = ",")),
        paste("Standard deviation:", format(round(summary_stats$sd_length), big.mark = ",")),
        "",
        "Quantiles:",
        paste("Minimum:", format(round(summary_stats$quantiles[1]), big.mark = ",")),
        paste("25th percentile:", format(round(summary_stats$quantiles[2]), big.mark = ",")),
        paste("Median:", format(round(summary_stats$quantiles[3]), big.mark = ",")),
        paste("75th percentile:", format(round(summary_stats$quantiles[4]), big.mark = ",")),
        paste("Maximum:", format(round(summary_stats$quantiles[5]), big.mark = ","))
    ),
    file.path(output_dir, "summary_statistics.txt")
)

# Save detailed gene list of sequences > 4kb
genes_over_4kb_details <- coding_analysis %>%
    filter(cds_length > 4000) %>%
    arrange(desc(cds_length))

write_csv(genes_over_4kb_details, 
          file.path(output_dir, "genes_over_4kb.csv"))

# Create visualizations with improved aesthetics
# 1. Regular histogram
p1 <- ggplot(coding_analysis, aes(x = cds_length)) +
    geom_histogram(bins = 100, fill = "#21908CFF", alpha = 0.7) +
    geom_vline(xintercept = 4000, color = "red", linetype = "dashed") +
    theme_minimal(base_size = 12) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = comma) +
    labs(
        title = "Distribution of Human Coding Sequence Lengths",
        subtitle = sprintf("%s of genes > 4kb", 
                         percent(summary_stats$fraction_over_4kb, accuracy = 0.1)),
        x = "Coding Sequence Length (bp)",
        y = "Count"
    )

# 2. Log-scaled histogram
p2 <- ggplot(coding_analysis, aes(x = cds_length)) +
    geom_histogram(bins = 100, fill = "#21908CFF", alpha = 0.7) +
    geom_vline(xintercept = 4000, color = "red", linetype = "dashed") +
    scale_x_log10(labels = comma) +
    scale_y_continuous(labels = comma) +
    theme_minimal(base_size = 12) +
    labs(
        title = "Distribution of Human Coding Sequence Lengths (Log Scale)",
        subtitle = sprintf("%s of genes > 4kb", 
                         percent(summary_stats$fraction_over_4kb, accuracy = 0.1)),
        x = "Coding Sequence Length (bp, log scale)",
        y = "Count"
    )

# 3. Cumulative distribution
p3 <- coding_analysis %>%
    arrange(cds_length) %>%
    mutate(cumulative_fraction = row_number() / n()) %>%
    ggplot(aes(x = cds_length, y = cumulative_fraction)) +
    geom_line(color = "#21908CFF", size = 1) +
    geom_vline(xintercept = 4000, color = "red", linetype = "dashed") +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = percent) +
    theme_minimal(base_size = 12) +
    labs(
        title = "Cumulative Distribution of Human Coding Sequence Lengths",
        x = "Coding Sequence Length (bp)",
        y = "Cumulative Fraction"
    )

# Save plots with higher resolution
ggsave(file.path(output_dir, "cds_length_distribution.pdf"), 
       p1, width = 10, height = 6, dpi = 300)
ggsave(file.path(output_dir, "cds_length_distribution_log.pdf"), 
       p2, width = 10, height = 6, dpi = 300)
ggsave(file.path(output_dir, "cds_length_cumulative.pdf"), 
       p3, width = 10, height = 6, dpi = 300)

message("Analysis complete. Results saved in human_coding_sequence_analysis directory.") 