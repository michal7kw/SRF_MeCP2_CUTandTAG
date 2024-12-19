#!/usr/bin/env Rscript

library(DESeq2)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(clusterProfiler)
library(org.Hs.eg.db)

# Set working directory
setwd("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_CUTandTAG/custom_pipeline")

# Create output directory
dir.create("results_1/differential", recursive = TRUE)

# Function to create count matrix
create_count_matrix <- function(peaks, bam_files, sample_names) {
    # Convert peaks to GRanges
    peaks_gr <- import(peaks)
    
    # Count reads in peaks
    counts <- matrix(0, nrow=length(peaks_gr), ncol=length(bam_files))
    colnames(counts) <- sample_names
    
    for(i in seq_along(bam_files)) {
        counts[,i] <- countOverlaps(peaks_gr, import(bam_files[i]))
    }
    
    return(counts)
}

# Function to perform differential analysis
run_deseq2 <- function(counts, sample_info) {
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = sample_info,
                                 design = ~ condition)
    
    # Run DESeq2
    dds <- DESeq(dds)
    
    # Get results
    res <- results(dds)
    
    return(list(dds=dds, res=res))
}

# Function to perform pathway analysis
run_pathway_analysis <- function(diff_peaks, output_prefix) {
    # Convert peak coordinates to gene names
    peak_annot <- annotatePeak(diff_peaks,
                              TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                              annoDb="org.Hs.eg.db")
    
    # Get gene list
    genes <- peak_annot@anno$ENTREZID
    
    # Run GO enrichment
    ego <- enrichGO(gene = genes,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05)
    
    # Save results
    write.csv(as.data.frame(ego), 
              file=paste0(output_prefix, "_GO_enrichment.csv"))
    
    return(ego)
}

# Main execution
main <- function() {
    # Get sample information
    samples <- data.frame(
        condition = factor(rep(c("EXOGENOUS", "ENDOGENOUS"), each=3)),
        sample = c(paste0("NeuV", 1:3), paste0("NSCV", 1:3))
    )
    
    # Get BAM files
    bam_files <- file.path("results_1/aligned", 
                          paste0(samples$sample, ".bam"))
    
    # Create consensus peak set
    consensus_peaks <- "results_1/peaks/consensus_peaks.bed"
    
    # Create count matrix
    counts <- create_count_matrix(consensus_peaks, bam_files, samples$sample)
    
    # Run differential analysis
    diff_results <- run_deseq2(counts, samples)
    
    # Save results
    write.csv(as.data.frame(diff_results$res),
              file="results_1/differential/differential_peaks.csv")
    
    # Run pathway analysis for significant peaks
    sig_peaks <- diff_results$res[which(diff_results$res$padj < 0.05),]
    pathway_results <- run_pathway_analysis(sig_peaks,
                                          "results_1/differential/pathway_analysis")
}

main() 