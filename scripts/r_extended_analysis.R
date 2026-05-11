# Extended analysis of Cut&Tag data

# Load necessary libraries
library(ChIPseeker) # For peak annotation and genomic feature analysis
library(TxDb.Mmusculus.UCSC.mm10.knownGene) # Transcript database for Mus musculus (mouse)
library(org.Mm.eg.db) # Annotation database for mouse
library(clusterProfiler) # For GO and pathway enrichment analysis
library(GenomicFeatures) # For genomic feature operations
library(ggplot2) # For plotting
library(ComplexHeatmap) # For creating complex heatmaps
library(rtracklayer) # For reading and writing genomic files
library(chromstaR) # For chromatin state analysis
library(DiffBind) # For differential binding analysis

# Set up parameters
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene # Assign the mouse transcript database
output_dir <- "results/extended_analysis" # Define the output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE) # Create the output directory if it doesn't exist

# Helper function to get sample groups
get_sample_groups <- function() {
    # Define sample groups for endogenous and exogenous conditions
    list(
        endogenous = list(
            neu = c("NeuM2", "NeuM3"), # Neuronal samples (endogenous)
            nsc = c("NSCM1", "NSCM2", "NSCM3"), # Neural stem cell samples (endogenous)
            control = "IgM" # Control sample (endogenous)
        ),
        exogenous = list(
            neu = c("NeuV1", "NeuV2", "NeuV3"), # Neuronal samples (exogenous)
            nsc = c("NSCv1", "NSCv2", "NSCv3") # Neural stem cell samples (exogenous)
        )
    )
}

# 1. Peak Annotation and Genomic Distribution Analysis
analyze_peak_distribution <- function(peak_files, output_prefix) {
    # This function analyzes the genomic distribution of peaks

    # Read peaks using ChIPseeker's readPeakFile function
    peak_list <- lapply(peak_files, readPeakFile)

    # Assign names to the peak lists based on the basename of the peak files
    names(peak_list) <- basename(peak_files)
    
    # Annotate peaks using ChIPseeker's annotatePeak function
    peakAnno_list <- lapply(peak_list, annotatePeak, 
                           TxDb=txdb, # Use the mouse transcript database
                           annoDb="org.Mm.eg.db") # Use the mouse annotation database
    
    # Plot genomic distribution using ChIPseeker's plotting functions
    pdf(file.path(output_dir, paste0(output_prefix, "_genomic_distribution.pdf"))) # Open a PDF file for plotting
    plotAnnoBar(peakAnno_list) # Plot the annotation bar plot
    plotDistToTSS(peakAnno_list) # Plot the distance to transcription start site (TSS)
    dev.off() # Close the PDF file
    
    return(peakAnno_list) # Return the list of annotated peaks
}

# 2. Chromatin State Analysis using ChromstaR
analyze_chromatin_states <- function(bam_files, peak_files) {
    # This function analyzes chromatin states using the chromstaR package

    # Prepare input for chromstaR using preparePeakInput function
    chromstar_data <- preparePeakInput(bam_files, peak_files)
    
    # Run chromatin state analysis using runChromstaR function
    chromstar_results <- runChromstaR(chromstar_data)
    
    # Plot chromatin state heatmap using plotChromStates function
    pdf(file.path(output_dir, "chromatin_states_heatmap.pdf")) # Open a PDF file for plotting
    plotChromStates(chromstar_results) # Plot the chromatin states
    dev.off() # Close the PDF file
    
    return(chromstar_results) # Return the results of the chromatin state analysis
}

# 3. GO and Pathway Enrichment Analysis
perform_enrichment_analysis <- function(peak_annotations) {
    # This function performs GO enrichment analysis on the peaks

    # Extract genes near peaks from the peak annotations
    genes <- lapply(peak_annotations, function(x) {
        unique(x@anno$geneId) # Get unique gene IDs
    })
    
    # Perform GO enrichment analysis using clusterProfiler's enrichGO function
    go_results <- lapply(genes, function(g) {
        enrichGO(gene = g, # Use the extracted gene IDs
                OrgDb = org.Mm.eg.db, # Use the mouse annotation database
                ont = "BP", # Perform enrichment for biological processes
                pAdjustMethod = "BH", # Use Benjamini-Hochberg method for p-value adjustment
                pvalueCutoff = 0.05) # Set p-value cutoff to 0.05
    })
    
    # Plot the GO enrichment results using dotplot
    pdf(file.path(output_dir, "go_enrichment.pdf")) # Open a PDF file for plotting
    for(i in seq_along(go_results)) {
        if(nrow(go_results[[i]]) > 0) {
            dotplot(go_results[[i]], showCategory=20, title=names(go_results)[i]) # Create dotplot for each GO result
        }
    }
    dev.off() # Close the PDF file
    
    return(go_results) # Return the GO enrichment results
}

# 4. Differential Binding Analysis Visualization
visualize_differential_binding <- function(dba_object) {
    # This function visualizes the results of differential binding analysis

    # Create PCA plot using DiffBind's dba.plotPCA function
    pdf(file.path(output_dir, "differential_binding_pca.pdf")) # Open a PDF file for plotting
    dba.plotPCA(dba_object, DBA_CONDITION, label=DBA_ID) # Plot PCA
    dev.off() # Close the PDF file
    
    # Create correlation heatmap using DiffBind's dba.plotHeatmap function
    pdf(file.path(output_dir, "sample_correlation_heatmap.pdf")) # Open a PDF file for plotting
    dba.plotHeatmap(dba_object) # Plot heatmap
    dev.off() # Close the PDF file
    
    # Create MA plots for contrasts using DiffBind's dba.plotMA function
    pdf(file.path(output_dir, "ma_plots.pdf")) # Open a PDF file for plotting
    dba.plotMA(dba_object) # Plot MA plots
    dev.off() # Close the PDF file
}

# Main execution
main <- function() {
    # This is the main function that orchestrates the analysis

    # Get sample groups using the helper function
    samples <- get_sample_groups()
    
    # Process endogenous samples
    endogenous_peaks <- list.files("results/peaks", 
                                  pattern="*_peaks.narrowPeak", # Look for files ending with "_peaks.narrowPeak"
                                  full.names=TRUE) # Get the full file paths
    
    # 1. Peak Distribution Analysis
    cat("Performing peak distribution analysis...\n") # Print a message to the console
    peak_annotations <- analyze_peak_distribution(endogenous_peaks, "endogenous") # Analyze peak distribution
    
    # 2. Chromatin State Analysis
    cat("Performing chromatin state analysis...\n") # Print a message to the console
    bam_files <- list.files("results/aligned", 
                           pattern="*.bam$", # Look for files ending with ".bam"
                           full.names=TRUE) # Get the full file paths
    chromatin_states <- analyze_chromatin_states(bam_files, endogenous_peaks) # Analyze chromatin states
    
    # 3. GO Enrichment Analysis
    cat("Performing GO enrichment analysis...\n") # Print a message to the console
    enrichment_results <- perform_enrichment_analysis(peak_annotations) # Perform GO enrichment analysis
    
    # 4. Create summary report
    cat("Creating summary report...\n") # Print a message to the console
    rmarkdown::render("scripts/analysis_report.Rmd", # Path to the R Markdown report
                     output_file=file.path(output_dir, "analysis_report.html"), # Output file path
                     params=list( # Parameters to pass to the R Markdown report
                         peak_annotations=peak_annotations, # Peak annotations
                         chromatin_states=chromatin_states, # Chromatin states
                         enrichment_results=enrichment_results # Enrichment results
                     ))
    
    cat("Analysis complete! Results are in:", output_dir, "\n") # Print a completion message
}

# Run the analysis
main() # Call the main function to start the analysis