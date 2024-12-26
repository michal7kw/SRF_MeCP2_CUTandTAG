#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(rtracklayer)
library(ComplexHeatmap)
library(circlize)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


RESULTS_DIR="results_1"
# Set working directory
setwd("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/iterative_alternative")

# Create visualization directories
dir.create("${RESULTS_DIR}/visualization", recursive = TRUE)
dir.create("${RESULTS_DIR}/visualization/heatmaps", recursive = TRUE)
dir.create("${RESULTS_DIR}/visualization/profiles", recursive = TRUE)
dir.create("${RESULTS_DIR}/visualization/peak_analysis", recursive = TRUE)

# 1. Fragment Size Distribution Analysis
plot_fragment_distribution <- function() {
    # Read fragment sizes
    frag_sizes <- read.table("${RESULTS_DIR}/qc/fragment_sizes_combined.txt", header=TRUE)
    
    # Create nucleosome periodicity plot
    p <- ggplot(frag_sizes, aes(x=size, y=count, color=sample)) +
        geom_line(size=1) +
        theme_minimal() +
        scale_x_continuous(breaks=seq(0, 1000, 100)) +
        labs(x="Fragment Size (bp)", 
             y="Count",
             title="Fragment Size Distribution",
             subtitle="Shows nucleosome periodicity") +
        theme(legend.position="bottom")
    
    ggsave("${RESULTS_DIR}/visualization/fragment_size_distribution.pdf", p, width=10, height=6)
}

# 2. Peak Profile Heatmaps
generate_peak_heatmaps <- function() {
    # Read bigWig files and peak files
    samples <- list.files("${RESULTS_DIR}/bigwig", pattern="*.bw$")
    peaks <- list.files("${RESULTS_DIR}/peaks/seacr", pattern="*.stringent.bed$")
    
    # Generate heatmap matrix
    computeMatrix <- function(bw, peaks, outfile) {
        system(paste("computeMatrix reference-point",
                    "-S", bw,
                    "-R", peaks,
                    "--referencePoint center",
                    "-b 2000 -a 2000",
                    "-o", outfile))
    }
    
    # Plot heatmap
    plotHeatmap <- function(matrix, outfile) {
        system(paste("plotHeatmap",
                    "-m", matrix,
                    "-o", outfile,
                    "--colorMap viridis",
                    "--whatToShow 'heatmap and colorbar'",
                    "--zMin 0 --zMax 100"))
    }
}

# 3. TSS Enrichment Analysis
plot_tss_enrichment <- function() {
    # Read TSS enrichment scores
    tss_scores <- read.table("${RESULTS_DIR}/qc/tss_enrichment_combined.txt", header=TRUE)
    
    p <- ggplot(tss_scores, aes(x=position, y=signal, color=sample)) +
        geom_line(size=1) +
        theme_minimal() +
        labs(x="Distance from TSS (bp)",
             y="Signal",
             title="TSS Enrichment Profile") +
        theme(legend.position="bottom")
    
    ggsave("${RESULTS_DIR}/visualization/tss_enrichment_profile.pdf", p, width=10, height=6)
}

# 4. Peak Annotation
annotate_peaks <- function() {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
    # Read peaks
    peaks <- list.files("${RESULTS_DIR}/peaks/seacr", pattern="*.stringent.bed$")
    
    for(peak in peaks) {
        # Convert to GRanges
        peak_gr <- import(peak)
        
        # Annotate peaks
        peak_annot <- annotatePeak(peak_gr, 
                                 TxDb=txdb,
                                 tssRegion=c(-3000, 3000))
        
        # Plot annotation
        pdf(paste0("${RESULTS_DIR}/visualization/peak_analysis/", 
                  gsub(".bed", "_annotation.pdf", basename(peak))))
        plotAnnoPie(peak_annot)
        plotDistToTSS(peak_annot)
        dev.off()
    }
}

# 5. Differential Binding Analysis Visualization
plot_differential_binding <- function() {
    # Read differential binding results
    diff_peaks <- read.table("${RESULTS_DIR}/differential/differential_peaks.txt", header=TRUE)
    
    # MA plot
    p1 <- ggplot(diff_peaks, aes(x=baseMean, y=log2FoldChange, color=padj < 0.05)) +
        geom_point(size=1) +
        scale_x_log10() +
        theme_minimal() +
        labs(x="Mean normalized counts",
             y="log2 fold change",
             title="MA Plot")
    
    # Volcano plot
    p2 <- ggplot(diff_peaks, aes(x=log2FoldChange, y=-log10(padj), color=padj < 0.05)) +
        geom_point(size=1) +
        theme_minimal() +
        labs(x="log2 fold change",
             y="-log10 adjusted p-value",
             title="Volcano Plot")
    
    ggsave("${RESULTS_DIR}/visualization/differential_MA_plot.pdf", p1, width=8, height=6)
    ggsave("${RESULTS_DIR}/visualization/differential_volcano_plot.pdf", p2, width=8, height=6)
}

# Main execution
main <- function() {
    plot_fragment_distribution()
    generate_peak_heatmaps()
    plot_tss_enrichment()
    annotate_peaks()
    plot_differential_binding()
}

main() 