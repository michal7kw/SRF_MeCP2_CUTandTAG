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

# Function to get coding sequence lengths using biomaRt
get_coding_lengths <- function() {
    message("Connecting to Ensembl...")
    
    # Connect to Ensembl
    mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    
    # Get coding sequence information
    message("Retrieving coding sequence information...")
    coding_seq_info <- getBM(
        attributes = c(
            "external_gene_name",
            "ensembl_gene_id",
            "cds_length",
            "gene_biotype"
        ),
        filters = "biotype",
        values = "protein_coding",
        mart = mart
    )
    
    return(coding_seq_info)
}

# Get coding sequence lengths
coding_lengths <- get_coding_lengths()

# Clean and process the data
coding_analysis <- coding_lengths %>%
    filter(!is.na(cds_length)) %>%
    group_by(ensembl_gene_id) %>%
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