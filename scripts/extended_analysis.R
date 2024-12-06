# Extended analysis of Cut&Tag data
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(clusterProfiler)
library(GenomicFeatures)
library(ggplot2)
library(ComplexHeatmap)
library(rtracklayer)
library(chromstaR)
library(DiffBind)

# Set up parameters
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
output_dir <- "results/extended_analysis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Helper function to get sample groups
get_sample_groups <- function() {
    list(
        endogenous = list(
            neu = c("NeuM2", "NeuM3"),
            nsc = c("NSCM1", "NSCM2", "NSCM3"),
            control = "IgM"
        ),
        exogenous = list(
            neu = c("NeuV1", "NeuV2", "NeuV3"),
            nsc = c("NSCv1", "NSCv2", "NSCv3")
        )
    )
}

# 1. Peak Annotation and Genomic Distribution Analysis
analyze_peak_distribution <- function(peak_files, output_prefix) {
    # Read peaks
    peak_list <- lapply(peak_files, readPeakFile)
    names(peak_list) <- basename(peak_files)
    
    # Annotate peaks
    peakAnno_list <- lapply(peak_list, annotatePeak, 
                           TxDb=txdb,
                           annoDb="org.Mm.eg.db")
    
    # Plot genomic distribution
    pdf(file.path(output_dir, paste0(output_prefix, "_genomic_distribution.pdf")))
    plotAnnoBar(peakAnno_list)
    plotDistToTSS(peakAnno_list)
    dev.off()
    
    return(peakAnno_list)
}

# 2. Chromatin State Analysis using ChromstaR
analyze_chromatin_states <- function(bam_files, peak_files) {
    # Prepare ChromstaR input
    chromstar_data <- preparePeakInput(bam_files, peak_files)
    
    # Run chromatin state analysis
    chromstar_results <- runChromstaR(chromstar_data)
    
    # Plot chromatin state heatmap
    pdf(file.path(output_dir, "chromatin_states_heatmap.pdf"))
    plotChromStates(chromstar_results)
    dev.off()
    
    return(chromstar_results)
}

# 3. GO and Pathway Enrichment Analysis
perform_enrichment_analysis <- function(peak_annotations) {
    # Extract genes near peaks
    genes <- lapply(peak_annotations, function(x) {
        unique(x@anno$geneId)
    })
    
    # Perform GO enrichment
    go_results <- lapply(genes, function(g) {
        enrichGO(gene = g,
                OrgDb = org.Mm.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)
    })
    
    # Plot results
    pdf(file.path(output_dir, "go_enrichment.pdf"))
    for(i in seq_along(go_results)) {
        if(nrow(go_results[[i]]) > 0) {
            dotplot(go_results[[i]], showCategory=20, title=names(go_results)[i])
        }
    }
    dev.off()
    
    return(go_results)
}

# 4. Differential Binding Analysis Visualization
visualize_differential_binding <- function(dba_object) {
    # Create PCA plot
    pdf(file.path(output_dir, "differential_binding_pca.pdf"))
    dba.plotPCA(dba_object, DBA_CONDITION, label=DBA_ID)
    dev.off()
    
    # Create correlation heatmap
    pdf(file.path(output_dir, "sample_correlation_heatmap.pdf"))
    dba.plotHeatmap(dba_object)
    dev.off()
    
    # Create MA plots for contrasts
    pdf(file.path(output_dir, "ma_plots.pdf"))
    dba.plotMA(dba_object)
    dev.off()
}

# Main execution
main <- function() {
    # Get sample groups
    samples <- get_sample_groups()
    
    # Process endogenous samples
    endogenous_peaks <- list.files("results/peaks", 
                                  pattern="*_peaks.narrowPeak", 
                                  full.names=TRUE)
    
    # 1. Peak Distribution Analysis
    cat("Performing peak distribution analysis...\n")
    peak_annotations <- analyze_peak_distribution(endogenous_peaks, "endogenous")
    
    # 2. Chromatin State Analysis
    cat("Performing chromatin state analysis...\n")
    bam_files <- list.files("results/aligned", 
                           pattern="*.bam$", 
                           full.names=TRUE)
    chromatin_states <- analyze_chromatin_states(bam_files, endogenous_peaks)
    
    # 3. GO Enrichment Analysis
    cat("Performing GO enrichment analysis...\n")
    enrichment_results <- perform_enrichment_analysis(peak_annotations)
    
    # 4. Create summary report
    cat("Creating summary report...\n")
    rmarkdown::render("scripts/analysis_report.Rmd",
                     output_file=file.path(output_dir, "analysis_report.html"),
                     params=list(
                         peak_annotations=peak_annotations,
                         chromatin_states=chromatin_states,
                         enrichment_results=enrichment_results
                     ))
    
    cat("Analysis complete! Results are in:", output_dir, "\n")
}

# Run the analysis
main() 