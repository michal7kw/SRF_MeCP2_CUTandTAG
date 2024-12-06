#!/usr/bin/env Rscript

# Load libraries in correct order
library(BiocGenerics)  # Load first to avoid conflicts
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
library(dplyr, warn.conflicts = FALSE)  # Load last to ensure its functions take precedence

# Parse command line arguments (optional)
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
    input_dir <- args[1]
    output_dir <- args[2]
} else {
    # Default directories
    input_dir <- "results"
    output_dir <- "results/extended_analysis"
}

# Create output directory immediately and check if successful
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) {
    stop("Failed to create output directory: ", output_dir)
}

# Set up parameters
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Helper function to get sample groups
get_sample_groups <- function() {
    # Add validation for sample names
    validate_sample <- function(sample_name) {
        if(is.na(sample_name) || !is.character(sample_name)) {
            warning("Invalid sample name detected")
            return(FALSE)
        }
        return(TRUE)
    }
    
    groups <- list(
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
    
    # Validate all sample names
    all_samples <- unlist(groups, recursive = TRUE)
    valid_samples <- sapply(all_samples, validate_sample)
    
    if(!all(valid_samples)) {
        warning("Some sample names are invalid")
    }
    
    return(groups)
}

# 1. Peak Annotation and Genomic Distribution Analysis
analyze_peak_distribution <- function(peak_files, output_prefix) {
    cat("Analyzing peak distribution...\n")
    
    if(length(peak_files) == 0) {
        stop("No peak files provided")
    }
    
    # Read peaks with error checking
    peak_list <- lapply(peak_files, function(file) {
        tryCatch({
            readPeakFile(file)
        }, error = function(e) {
            warning("Failed to read peak file: ", file)
            return(NULL)
        })
    })
    
    peak_list <- peak_list[!sapply(peak_list, is.null)]
    if(length(peak_list) == 0) {
        stop("No valid peak files could be read")
    }
    
    names(peak_list) <- basename(peak_files)[1:length(peak_list)]
    
    # Annotate peaks
    peakAnno_list <- lapply(peak_list, annotatePeak, 
                           TxDb=txdb,
                           annoDb="org.Mm.eg.db")
    
    # Plot with error checking
    pdf_file <- file.path(output_dir, paste0(output_prefix, "_genomic_distribution.pdf"))
    tryCatch({
        pdf(pdf_file)
        plotAnnoBar(peakAnno_list)
        plotDistToTSS(peakAnno_list)
        dev.off()
        
        # Verify PDF was created and has content
        if(!file.exists(pdf_file) || file.size(pdf_file) < 100) {
            warning("PDF file may be empty or invalid: ", pdf_file)
        }
    }, error = function(e) {
        if(dev.cur() > 1) dev.off()  # Ensure device is closed on error
        warning("Failed to create plot: ", e$message)
    })
    
    return(peakAnno_list)
}

# 2. Chromatin State Analysis using binned approach
analyze_chromatin_states <- function(bam_files, peak_files) {
    cat("Analyzing chromatin states...\n")
    
    # Filter out unsorted BAM files
    bam_files <- bam_files[!grepl("unsorted", bam_files)]
    
    # Validate BAM files and their indices
    valid_bams <- sapply(bam_files, function(bam) {
        bai_file <- paste0(bam, ".bai")
        if (!file.exists(bai_file)) {
            warning("BAM index missing for: ", bam)
            return(FALSE)
        }
        tryCatch({
            header <- Rsamtools::scanBamHeader(bam)
            return(TRUE)
        }, error = function(e) {
            warning("Invalid or corrupted BAM file: ", bam)
            return(FALSE)
        })
    })
    
    bam_files <- bam_files[valid_bams]
    if(length(bam_files) == 0) {
        stop("No valid BAM files found")
    }
    
    # Create experiment list with validation
    experiments <- data.frame(
        file = bam_files,
        condition = ifelse(grepl("Neu", basename(bam_files)), "Neuron", "NSC"),
        replicate = as.integer(gsub(".*([0-9]+).*", "\\1", basename(bam_files))),
        stringsAsFactors = FALSE
    )
    
    # Remove rows with NA values
    experiments <- experiments[complete.cases(experiments), ]
    
    # Increase bin size to reduce memory usage
    binsize <- 1000  # Changed from 200 to 1000
    
    # Get proper genome information
    library(GenomeInfoDb)
    # Get chromosome information from the first BAM file
    first_bam <- Rsamtools::BamFile(bam_files[1])
    si <- seqinfo(first_bam)
    
    # Create genomic bins
    bins <- tileGenome(si, tilewidth=binsize, cut.last.tile.in.chrom=TRUE)
    
    # Process BAM files and get counts
    bin_counts <- lapply(bam_files, function(bam) {
        tryCatch({
            reads <- GenomicAlignments::readGAlignments(bam)
            counts <- countOverlaps(bins, reads)
            return(counts)
        }, error = function(e) {
            message("Error processing BAM file: ", bam)
            message(e$message)
            return(NULL)
        })
    })
    
    # Remove any NULL results from failed processing
    bin_counts <- bin_counts[!sapply(bin_counts, is.null)]
    
    if (length(bin_counts) == 0) {
        stop("No valid BAM files could be processed")
    }
    
    # Convert to matrix
    count_matrix <- do.call(cbind, bin_counts)
    colnames(count_matrix) <- basename(bam_files)[1:ncol(count_matrix)]
    
    # Generate heatmap with memory management
    tryCatch({
        pdf(file.path(output_dir, "chromatin_states_heatmap.pdf"))
        # Subsample data if too large
        if(ncol(count_matrix) > 1000) {
            selected_bins <- sample(1:ncol(count_matrix), 1000)
            heatmap_data <- t(scale(count_matrix[, selected_bins]))
        } else {
            heatmap_data <- t(scale(count_matrix))
        }
        # Use ComplexHeatmap instead of base heatmap
        Heatmap(heatmap_data,
                name = "Z-score",
                show_row_names = TRUE,
                show_column_names = FALSE,
                column_title = "Genomic Bins",
                row_title = "Samples")
        dev.off()
    }, error = function(e) {
        message("Error generating heatmap: ", e$message)
    })
    
    # Create results object with similar structure to original
    chromstar_results <- list(
        counts = count_matrix,
        bins = bins,
        experiments = experiments,
        genome = genome(si)[1]  # Get genome info from seqinfo
    )
    
    return(chromstar_results)
}

# 3. GO and Pathway Enrichment Analysis
perform_enrichment_analysis <- function(peak_annotations) {
    cat("Performing GO enrichment analysis...\n")
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
visualize_differential_binding <- function(peak_files, bam_files, sample_groups) {
    cat("Performing differential binding analysis...\n")
    # Create sample sheet
    samples <- data.frame(
        SampleID = basename(bam_files),
        Condition = ifelse(grepl("Neu", basename(bam_files)), "Neuron", "NSC"),
        bamReads = bam_files,
        Peaks = peak_files,
        PeakCaller = "narrow"
    )
    
    # Create DBA object
    dba_obj <- dba(sampleSheet=samples)
    
    # Create PCA plot
    pdf(file.path(output_dir, "differential_binding_pca.pdf"))
    dba.plotPCA(dba_obj, attributes=DBA_CONDITION)
    dev.off()
    
    # Create correlation heatmap
    pdf(file.path(output_dir, "sample_correlation_heatmap.pdf"))
    dba.plotHeatmap(dba_obj)
    dev.off()
    
    return(dba_obj)
}

# Main execution
main <- function() {
    cat("Starting extended analysis...\n")
    
    # Verify input directories exist
    if (!dir.exists(file.path(input_dir, "peaks")) || !dir.exists(file.path(input_dir, "aligned"))) {
        stop("Input directories not found. Please check your input directory structure.")
    }
    
    # Get sample groups
    samples <- get_sample_groups()
    
    # Get input files
    peak_files <- list.files(file.path(input_dir, "peaks"), 
                            pattern="*_peaks.narrowPeak", 
                            full.names=TRUE)
    bam_files <- list.files(file.path(input_dir, "aligned"), 
                           pattern="*.bam$", 
                           full.names=TRUE)
    
    if(length(peak_files) == 0 || length(bam_files) == 0) {
        stop("No input files found. Please check your input directory.")
    }
    
    # Create results directory if it doesn't exist
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # 1. Peak Distribution Analysis
    peak_annotations <- analyze_peak_distribution(peak_files, "all_samples")
    
    # 2. Chromatin State Analysis
    chromatin_states <- analyze_chromatin_states(bam_files, peak_files)
    
    # 3. GO Enrichment Analysis
    enrichment_results <- perform_enrichment_analysis(peak_annotations)
    
    # 4. Differential Binding Analysis
    diff_binding <- visualize_differential_binding(peak_files, bam_files, samples)
    
    # 5. Create summary report
    cat("Creating summary report...\n")
    
    # Check if report template exists
    report_template <- "scripts/analysis_report.Rmd"
    if (!file.exists(report_template)) {
        warning("Report template not found. Skipping report generation.")
    } else {
        tryCatch({
            rmarkdown::render(report_template,
                         output_file=file.path(output_dir, "analysis_report.html"),
                         params=list(
                             peak_annotations=peak_annotations,
                             chromatin_states=chromatin_states,
                             enrichment_results=enrichment_results,
                             diff_binding=diff_binding
                         ))
        }, error = function(e) {
            warning("Failed to generate report: ", e$message)
        })
    }
    
    cat("Analysis complete! Results are in:", output_dir, "\n")
}

# Run the analysis if script is executed directly
if (!interactive()) {
    main()
} 