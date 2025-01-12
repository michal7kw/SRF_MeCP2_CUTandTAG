#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(rmarkdown)

# Set working directory
setwd("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative")

# Function to read fragment sizes
read_fragment_sizes <- function(sample) {
    read.table(paste0("results_1b/qc/fragment_sizes/", sample, "_sizes.txt")) %>%
        mutate(sample = sample)
}

# Function to read FRiP scores
read_frip <- function(sample) {
    read.table(paste0("results_1b/qc/frip/", sample, "_frip.txt")) %>%
        mutate(frip = V2/V3 * 100,
               sample = sample) # Add sample column to FRiP data
}

# Get sample names
exo_samples <- list.files("../DATA/EXOGENOUS", "*_R1_001.fastq.gz") %>%
    gsub("_R1_001.fastq.gz", "", .)
endo_samples <- list.files("../DATA/ENDOGENOUS", "*_R1_001.fastq.gz") %>%
    gsub("_R1_001.fastq.gz", "", .)
all_samples <- c(exo_samples, endo_samples)
all_samples <- all_samples[!grepl("IgM", all_samples)]

# Read and combine data
fragment_sizes <- lapply(all_samples, read_fragment_sizes) %>%
    bind_rows()

frip_scores <- lapply(all_samples, read_frip) %>%
    bind_rows()

# Generate plots
p1 <- ggplot(fragment_sizes, aes(x=V1, y=V2, color=sample)) +
    geom_line() +
    theme_minimal() +
    labs(x="Fragment Size", y="Count", title="Fragment Size Distribution")

# Convert sample to factor to avoid the error with x aesthetic
frip_scores$sample <- as.factor(frip_scores$sample)

p2 <- ggplot(frip_scores, aes(x=sample, y=frip)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x="Sample", y="FRiP %", title="Fraction of Reads in Peaks")

# Save plots
pdf_dir <- "results_1b/qc"
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)  # Ensure directory exists

ggsave(file.path(pdf_dir, "fragment_size_dist.pdf"), p1, width=10, height=6)
ggsave(file.path(pdf_dir, "frip_scores.pdf"), p2, width=10, height=6)

# Generate HTML report
rmarkdown::render("qc_report.Rmd",
                 output_file = "results_1b/qc/qc_report.html",
                 knit_root_dir = getwd())  # Explicitly set working directory for knitting