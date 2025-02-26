#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(rtracklayer)
library(ComplexHeatmap)
library(circlize)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(EnrichedHeatmap)
library(gridExtra)
library(RColorBrewer)

#### SAMPLES ####
# /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/peaks/narrow$ ls -1 *.filtered.narrowPeak
# NeuM2_narrow_peaks.filtered.narrowPeak
# NeuM3_narrow_peaks.filtered.narrowPeak
# NeuV1_narrow_peaks.filtered.narrowPeak
# NeuV2_narrow_peaks.filtered.narrowPeak
# NeuV3_narrow_peaks.filtered.narrowPeak
# NSCM1_narrow_peaks.filtered.narrowPeak
# NSCM2_narrow_peaks.filtered.narrowPeak
# NSCM3_narrow_peaks.filtered.narrowPeak
# NSCv1_narrow_peaks.filtered.narrowPeak
# NSCv2_narrow_peaks.filtered.narrowPeak
# NSCv3_narrow_peaks.filtered.narrowPeak

### WHERE ###
# NeuM2,Neuron,Endogenous
# NeuM3,Neuron,Endogenous
# NSCM1,NSC,Endogenous
# NSCM2,NSC,Endogenous
# NSCM3,NSC,Endogenous

# NeuV1,Neuron,Exogenous
# NeuV2,Neuron,Exogenous
# NeuV3,Neuron,Exogenous
# NSCv1,NSC,Exogenous
# NSCv2,NSC,Exogenous
# NSCv3,NSC,Exogenous

setwd("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/paper_finals")
OUTPUT_DIR <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/paper_finals/outputs"
BIGWIG_DIR <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results_1b/bigwig"
PEAKS_DIR <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/peaks/narrow"
TSS_BED <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA/mm10_TSS.bed"

# Create output directory if it doesn't exist
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Function to generate heatmaps and metaprofiles for Mecp2 Endo and Exo around TSS
generate_mecp2_tss_heatmaps <- function() {
  # Define sample groups
  endo_samples <- c("NeuM2", "NeuM3", "NSCM1", "NSCM2", "NSCM3")
  exo_samples <- c("NeuV1", "NeuV2", "NeuV3", "NSCv1", "NSCv2", "NSCv3")
  
  # Group by cell type and MeCP2 type
  neuron_endo <- c("NeuM2", "NeuM3")
  nsc_endo <- c("NSCM1", "NSCM2", "NSCM3")
  neuron_exo <- c("NeuV1", "NeuV2", "NeuV3")
  nsc_exo <- c("NSCv1", "NSCv2", "NSCv3")
  
  # Create temporary matrix files
  temp_dir <- file.path(OUTPUT_DIR, "temp")
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # First, let's prepare a clean TSS file to avoid the skipping issues
  prepare_tss_file <- function() {
    # Define the output file path
    clean_tss_file <- file.path(temp_dir, "clean_tss.bed")
    
    # Check if the file already exists
    if (file.exists(clean_tss_file)) {
      message(paste("Clean TSS file already exists, reusing:", clean_tss_file))
      return(clean_tss_file)
    }
    
    # Check if TSS file exists
    if (!file.exists(TSS_BED)) {
      stop("TSS BED file not found: ", TSS_BED)
    }
    
    # Read the TSS file
    tss <- tryCatch({
      import(TSS_BED)
    }, error = function(e) {
      stop("Error importing TSS BED file: ", e$message)
    })
    
    # Create a simplified version with unique names
    tss_clean <- tss
    names(tss_clean) <- paste0("gene_", seq_along(tss_clean))
    
    # Write to a temporary file
    export(tss_clean, clean_tss_file)
    
    return(clean_tss_file)
  }
  
  # Prepare clean TSS file
  clean_tss_file <- prepare_tss_file()
  
  # Function to compute matrix using deepTools
  compute_matrix <- function(bw_files, regions, output_file, upstream=5000, downstream=5000) {
    # Check if output file already exists
    if (file.exists(output_file)) {
      message(paste("Matrix file already exists, skipping:", output_file))
      return(output_file)
    }
    
    bw_files_str <- paste(bw_files, collapse = " ")
    cmd <- paste0("computeMatrix reference-point",
                 " -S ", bw_files_str,
                 " -R ", regions,
                 " --referencePoint TSS",
                 " -b ", upstream, " -a ", downstream,
                 " -o ", output_file,
                 " --skipZeros",
                 " --missingDataAsZero")
    print(cmd)
    system(cmd)
  }
  
  # Function to plot heatmap using deepTools
  plot_heatmap <- function(matrix_file, output_file, title="", zmin=0, zmax=0.3) {
    # Check if output file already exists
    if (file.exists(output_file)) {
      message(paste("Output file already exists, skipping:", output_file))
      return(output_file)
    }
    
    cmd <- paste0("plotHeatmap -m ", matrix_file,
                 " -o ", output_file,
                 " --colorMap 'Blues'",
                 " --whatToShow 'heatmap and colorbar'",
                 " --zMin ", zmin, " --zMax ", zmax,
                 " --heatmapHeight 15",
                 " --heatmapWidth 7.5",
                 " --xAxisLabel \"Distance from TSS (bp)\"",
                 " --refPointLabel \"TSS\"",
                 " --regionsLabel \"Genes\"",
                 " --plotTitle \"", title, "\"")
    print(cmd)
    system(cmd)
  }
  
  # Function to plot profile using deepTools - fixed parameters
  plot_profile <- function(matrix_file, output_file, title="") {
    # Check if output file already exists
    if (file.exists(output_file)) {
      message(paste("Output file already exists, skipping:", output_file))
      return(output_file)
    }
    
    # Get the number of samples in the matrix
    # Use computeMatrixOperations to get info about the matrix
    temp_info_file <- tempfile()
    cmd_info <- paste0("computeMatrixOperations info -m ", matrix_file, " > ", temp_info_file)
    system(cmd_info)
    
    # Read the first few lines to determine number of samples
    info_lines <- readLines(temp_info_file, n=20)
    file.remove(temp_info_file)
    
    # Find the line with sample info
    sample_line <- grep("sample_labels", info_lines, value=TRUE)
    if (length(sample_line) > 0) {
      # Extract sample count
      sample_count <- length(strsplit(gsub(".*\\[|\\].*", "", sample_line), ",")[[1]])
    } else {
      # Default to 2 if we can't determine
      sample_count <- 2
    }
    
    # Generate enough colors
    colors <- paste(rep(c("blue", "red", "green", "orange", "purple", "brown"), length.out=sample_count), collapse=" ")
    
    cmd <- paste0("plotProfile -m ", matrix_file,
                 " -o ", output_file,
                 " --colors ", colors,
                 " --plotTitle \"", title, "\"")
    print(cmd)
    system(cmd)
  }
  
  # Generate matrices and heatmaps for Neuron Endogenous vs Exogenous
  neuron_endo_bw <- file.path(BIGWIG_DIR, paste0(neuron_endo, ".bw"))
  neuron_exo_bw <- file.path(BIGWIG_DIR, paste0(neuron_exo, ".bw"))
  
  # Combine bigwig files for each group
  neuron_endo_matrix <- file.path(temp_dir, "neuron_endo_tss_matrix.gz")
  neuron_exo_matrix <- file.path(temp_dir, "neuron_exo_tss_matrix.gz")
  
  # Compute matrices
  compute_matrix(neuron_endo_bw, clean_tss_file, neuron_endo_matrix)
  compute_matrix(neuron_exo_bw, clean_tss_file, neuron_exo_matrix)
  
  # Plot heatmaps
  neuron_endo_heatmap <- file.path(OUTPUT_DIR, "neuron_endo_tss_heatmap.png")
  if (!file.exists(neuron_endo_heatmap)) {
    plot_heatmap(neuron_endo_matrix, neuron_endo_heatmap, "Neuron Endogenous MeCP2 around TSS")
  }
  
  neuron_exo_heatmap <- file.path(OUTPUT_DIR, "neuron_exo_tss_heatmap.png")
  if (!file.exists(neuron_exo_heatmap)) {
    plot_heatmap(neuron_exo_matrix, neuron_exo_heatmap, "Neuron Exogenous MeCP2 around TSS")
  }
  
  # Plot profiles
  neuron_endo_profile <- file.path(OUTPUT_DIR, "neuron_endo_tss_profile.png")
  if (!file.exists(neuron_endo_profile)) {
    plot_profile(neuron_endo_matrix, neuron_endo_profile, "Neuron Endogenous MeCP2 around TSS")
  }
  
  neuron_exo_profile <- file.path(OUTPUT_DIR, "neuron_exo_tss_profile.png")
  if (!file.exists(neuron_exo_profile)) {
    plot_profile(neuron_exo_matrix, neuron_exo_profile, "Neuron Exogenous MeCP2 around TSS")
  }
  
  # Generate matrices and heatmaps for NSC Endogenous vs Exogenous
  nsc_endo_bw <- file.path(BIGWIG_DIR, paste0(nsc_endo, ".bw"))
  nsc_exo_bw <- file.path(BIGWIG_DIR, paste0(nsc_exo, ".bw"))
  
  # Combine bigwig files for each group
  nsc_endo_matrix <- file.path(temp_dir, "nsc_endo_tss_matrix.gz")
  nsc_exo_matrix <- file.path(temp_dir, "nsc_exo_tss_matrix.gz")
  
  # Compute matrices
  compute_matrix(nsc_endo_bw, clean_tss_file, nsc_endo_matrix)
  compute_matrix(nsc_exo_bw, clean_tss_file, nsc_exo_matrix)
  
  # Plot heatmaps
  nsc_endo_heatmap <- file.path(OUTPUT_DIR, "nsc_endo_tss_heatmap.png")
  if (!file.exists(nsc_endo_heatmap)) {
    plot_heatmap(nsc_endo_matrix, nsc_endo_heatmap, "NSC Endogenous MeCP2 around TSS")
  }
  
  nsc_exo_heatmap <- file.path(OUTPUT_DIR, "nsc_exo_tss_heatmap.png")
  if (!file.exists(nsc_exo_heatmap)) {
    plot_heatmap(nsc_exo_matrix, nsc_exo_heatmap, "NSC Exogenous MeCP2 around TSS")
  }
  
  # Plot profiles
  nsc_endo_profile <- file.path(OUTPUT_DIR, "nsc_endo_tss_profile.png")
  if (!file.exists(nsc_endo_profile)) {
    plot_profile(nsc_endo_matrix, nsc_endo_profile, "NSC Endogenous MeCP2 around TSS")
  }
  
  nsc_exo_profile <- file.path(OUTPUT_DIR, "nsc_exo_tss_profile.png")
  if (!file.exists(nsc_exo_profile)) {
    plot_profile(nsc_exo_matrix, nsc_exo_profile, "NSC Exogenous MeCP2 around TSS")
  }
  
  # Create combined plots for comparison
  # Combine Endo vs Exo for Neurons
  neuron_combined_matrix <- file.path(temp_dir, "neuron_combined_tss_matrix.gz")
  compute_matrix(c(neuron_endo_bw[1], neuron_exo_bw[1]), clean_tss_file, neuron_combined_matrix)
  
  neuron_combined_heatmap <- file.path(OUTPUT_DIR, "neuron_endo_vs_exo_tss_heatmap.png")
  if (!file.exists(neuron_combined_heatmap)) {
    plot_heatmap(neuron_combined_matrix, neuron_combined_heatmap, "Neuron Endogenous vs Exogenous MeCP2 around TSS")
  }
  
  neuron_combined_profile <- file.path(OUTPUT_DIR, "neuron_endo_vs_exo_tss_profile.png")
  if (!file.exists(neuron_combined_profile)) {
    plot_profile(neuron_combined_matrix, neuron_combined_profile, "Neuron Endogenous vs Exogenous MeCP2 around TSS")
  }
  
  # Combine Endo vs Exo for NSCs
  nsc_combined_matrix <- file.path(temp_dir, "nsc_combined_tss_matrix.gz")
  compute_matrix(c(nsc_endo_bw[1], nsc_exo_bw[1]), clean_tss_file, nsc_combined_matrix)
  
  nsc_combined_heatmap <- file.path(OUTPUT_DIR, "nsc_endo_vs_exo_tss_heatmap.png")
  if (!file.exists(nsc_combined_heatmap)) {
    plot_heatmap(nsc_combined_matrix, nsc_combined_heatmap, "NSC Endogenous vs Exogenous MeCP2 around TSS")
  }
  
  nsc_combined_profile <- file.path(OUTPUT_DIR, "nsc_endo_vs_exo_tss_profile.png")
  if (!file.exists(nsc_combined_profile)) {
    plot_profile(nsc_combined_matrix, nsc_combined_profile, "NSC Endogenous vs Exogenous MeCP2 around TSS")
  }
  
  # Add CpG islands analysis
  generate_cpg_island_profiles <- function() {
    # Define CpG islands file
    cpg_islands_file <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA/cpg_islands.bed"
    
    # Check if CpG islands file exists
    if (!file.exists(cpg_islands_file)) {
      stop("CpG islands BED file not found: ", cpg_islands_file)
    }
    
    # Create a clean CpG islands file with unique names
    prepare_cpg_islands_file <- function() {
      # Define the output file path
      clean_cpg_file <- file.path(temp_dir, "clean_cpg_islands.bed")
      
      # Check if the file already exists
      if (file.exists(clean_cpg_file)) {
        message(paste("Clean CpG islands file already exists, reusing:", clean_cpg_file))
        return(clean_cpg_file)
      }
      
      # Read the CpG islands file
      cpg_islands <- tryCatch({
        import(cpg_islands_file)
      }, error = function(e) {
        stop("Error importing CpG islands BED file: ", e$message)
      })
      
      # Create a simplified version with unique names
      cpg_islands_clean <- cpg_islands
      names(cpg_islands_clean) <- paste0("cpg_", seq_along(cpg_islands_clean))
      
      # Write to a temporary file
      export(cpg_islands_clean, clean_cpg_file)
      
      return(clean_cpg_file)
    }
    
    # Prepare clean CpG islands file
    clean_cpg_file <- prepare_cpg_islands_file()
    
    # Generate matrices for Endo vs Exo around CpG islands
    # Combine all endogenous samples
    all_endo_bw <- file.path(BIGWIG_DIR, paste0(endo_samples, ".bw"))
    all_exo_bw <- file.path(BIGWIG_DIR, paste0(exo_samples, ".bw"))
    
    # Create matrices
    endo_cpg_matrix <- file.path(temp_dir, "endo_cpg_matrix.gz")
    exo_cpg_matrix <- file.path(temp_dir, "exo_cpg_matrix.gz")
    
    # Compute matrices with 3kb upstream and downstream
    if (!file.exists(endo_cpg_matrix)) {
      compute_matrix(all_endo_bw, clean_cpg_file, endo_cpg_matrix, upstream=3000, downstream=3000)
    } else {
      message(paste("Matrix file already exists, reusing:", endo_cpg_matrix))
    }
    
    if (!file.exists(exo_cpg_matrix)) {
      compute_matrix(all_exo_bw, clean_cpg_file, exo_cpg_matrix, upstream=3000, downstream=3000)
    } else {
      message(paste("Matrix file already exists, reusing:", exo_cpg_matrix))
    }
    
    # Combined matrix for Endo vs Exo
    combined_cpg_matrix <- file.path(temp_dir, "combined_cpg_matrix.gz")
    if (!file.exists(combined_cpg_matrix)) {
      # Use representative samples from each group
      compute_matrix(c(all_endo_bw[1], all_exo_bw[1]), clean_cpg_file, combined_cpg_matrix, upstream=3000, downstream=3000)
    } else {
      message(paste("Matrix file already exists, reusing:", combined_cpg_matrix))
    }
    
    # Plot heatmaps
    endo_cpg_heatmap <- file.path(OUTPUT_DIR, "endo_cpg_heatmap.png")
    if (!file.exists(endo_cpg_heatmap)) {
      plot_heatmap(endo_cpg_matrix, endo_cpg_heatmap, "Endogenous MeCP2 around CpG Islands")
    } else {
      message(paste("Heatmap file already exists, skipping:", endo_cpg_heatmap))
    }
    
    exo_cpg_heatmap <- file.path(OUTPUT_DIR, "exo_cpg_heatmap.png")
    if (!file.exists(exo_cpg_heatmap)) {
      plot_heatmap(exo_cpg_matrix, exo_cpg_heatmap, "Exogenous MeCP2 around CpG Islands")
    } else {
      message(paste("Heatmap file already exists, skipping:", exo_cpg_heatmap))
    }
    
    # Plot combined heatmap
    combined_cpg_heatmap <- file.path(OUTPUT_DIR, "endo_vs_exo_cpg_heatmap.png")
    if (!file.exists(combined_cpg_heatmap)) {
      plot_heatmap(combined_cpg_matrix, combined_cpg_heatmap, "Endogenous vs Exogenous MeCP2 around CpG Islands")
    } else {
      message(paste("Heatmap file already exists, skipping:", combined_cpg_heatmap))
    }
    
    # Plot profiles
    endo_cpg_profile <- file.path(OUTPUT_DIR, "endo_cpg_profile.png")
    if (!file.exists(endo_cpg_profile)) {
      plot_profile(endo_cpg_matrix, endo_cpg_profile, "Endogenous MeCP2 around CpG Islands")
    } else {
      message(paste("Profile file already exists, skipping:", endo_cpg_profile))
    }
    
    exo_cpg_profile <- file.path(OUTPUT_DIR, "exo_cpg_profile.png")
    if (!file.exists(exo_cpg_profile)) {
      plot_profile(exo_cpg_matrix, exo_cpg_profile, "Exogenous MeCP2 around CpG Islands")
    } else {
      message(paste("Profile file already exists, skipping:", exo_cpg_profile))
    }
    
    # Plot combined profile
    combined_cpg_profile <- file.path(OUTPUT_DIR, "endo_vs_exo_cpg_profile.png")
    if (!file.exists(combined_cpg_profile)) {
      plot_profile(combined_cpg_matrix, combined_cpg_profile, "Endogenous vs Exogenous MeCP2 around CpG Islands")
    } else {
      message(paste("Profile file already exists, skipping:", combined_cpg_profile))
    }
    
    # Create custom R plot with ggplot2 similar to the image provided
    create_custom_cpg_profile <- function() {
      custom_cpg_profile <- file.path(OUTPUT_DIR, "custom_endo_vs_exo_cpg_profile.png")
      
      if (file.exists(custom_cpg_profile)) {
        message(paste("Custom profile file already exists, skipping:", custom_cpg_profile))
        return(custom_cpg_profile)
      }
      
      # Extract data from matrices
      extract_matrix_data <- function(matrix_file) {
        # Use system command to extract data from the matrix
        temp_data_file <- tempfile(fileext = ".tab")
        cmd <- paste0("computeMatrixOperations info -m ", matrix_file, " > ", temp_data_file)
        system(cmd)
        
        # Read the data
        data <- tryCatch({
          read.table(temp_data_file, header = TRUE, sep = "\t", skip = 1)
        }, error = function(e) {
          warning("Error reading matrix data: ", e$message)
          return(NULL)
        }, finally = {
          # Clean up
          if (file.exists(temp_data_file)) file.remove(temp_data_file)
        })
        
        return(data)
      }
      
      # Extract data from both matrices
      endo_data <- extract_matrix_data(endo_cpg_matrix)
      exo_data <- extract_matrix_data(exo_cpg_matrix)
      
      if (is.null(endo_data) || is.null(exo_data)) {
        warning("Could not extract data from matrices for custom plot")
        return(NULL)
      }
      
      # Create a data frame for plotting
      endo_mean <- colMeans(as.matrix(endo_data), na.rm = TRUE)
      exo_mean <- colMeans(as.matrix(exo_data), na.rm = TRUE)
      
      # Create position vector (from -3kb to +3kb)
      num_points <- length(endo_mean)
      positions <- seq(-3000, 3000, length.out = num_points)
      
      plot_data <- data.frame(
        Position = rep(positions, 2),
        Enrichment = c(endo_mean, exo_mean),
        Type = factor(rep(c("Endogenous", "Exogenous"), each = num_points))
      )
      
      # Create the plot
      p <- ggplot(plot_data, aes(x = Position, y = Enrichment, color = Type)) +
        geom_line(size = 1) +
        scale_color_manual(values = c("Endogenous" = "lightblue", "Exogenous" = "blue")) +
        labs(
          x = "Distance from CpG island (bp)",
          y = "MeCP2 enrichment",
          title = "MeCP2 enrichment around CpG islands"
        ) +
        theme_minimal() +
        theme(
          legend.position = "top",
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
        ) +
        annotate("text", x = 0, y = max(c(endo_mean, exo_mean)) * 0.9, 
                 label = "CpG island", size = 4)
      
      # Save the plot
      ggsave(custom_cpg_profile, p, width = 10, height = 6, dpi = 300)
      
      return(custom_cpg_profile)
    }
    
    # Create custom profile
    custom_profile <- create_custom_cpg_profile()
    
    return("Completed generating CpG island profiles")
  }
  
  # Run CpG island analysis
  tryCatch({
    cpg_result <- generate_cpg_island_profiles()
    message(cpg_result)
  }, error = function(e) {
    warning("Error in CpG island analysis: ", e$message)
  })
  
  # Generate genome-wide distribution analysis
  analyze_genomic_distribution <- function() {
    # Load TxDb for mouse genome
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    
    # Function to analyze peaks
    analyze_peaks <- function(peak_file, sample_name) {
      # Import peaks
      peaks <- tryCatch({
        import(peak_file)
      }, error = function(e) {
        warning("Error importing peak file: ", e$message)
        return(NULL)
      })
      
      if (is.null(peaks)) {
        warning("Skipping peak analysis for ", sample_name, " due to import error")
        return(NULL)
      }
      
      # Annotate peaks
      peak_annot <- tryCatch({
        annotatePeak(peaks, 
                   tssRegion=c(-3000, 3000),
                   TxDb=txdb,
                   annoDb="org.Mm.eg.db")
      }, error = function(e) {
        warning("Error annotating peaks for ", sample_name, ": ", e$message)
        return(NULL)
      })
      
      if (is.null(peak_annot)) {
        warning("Skipping annotation plotting for ", sample_name, " due to annotation error")
        return(NULL)
      }
      
      # Create output directory for genomic distribution plots
      genomic_dir <- file.path(OUTPUT_DIR, "genomic_distribution")
      dir.create(genomic_dir, showWarnings = FALSE, recursive = TRUE)
      
      # Plot genomic annotation
      pdf(file.path(genomic_dir, paste0(sample_name, "_genomic_annotation.pdf")), width=10, height=8)
      tryCatch({
        print(plotAnnoPie(peak_annot))
        print(plotDistToTSS(peak_annot))
      }, error = function(e) {
        warning("Error plotting annotation for ", sample_name, ": ", e$message)
      }, finally = {
        dev.off()
      })
      
      return(peak_annot)
    }
    
    # Analyze peaks for each group
    neuron_endo_peaks <- file.path(PEAKS_DIR, "NeuM2_narrow_peaks.filtered.narrowPeak")
    neuron_exo_peaks <- file.path(PEAKS_DIR, "NeuV1_narrow_peaks.filtered.narrowPeak")
    nsc_endo_peaks <- file.path(PEAKS_DIR, "NSCM1_narrow_peaks.filtered.narrowPeak")
    nsc_exo_peaks <- file.path(PEAKS_DIR, "NSCv1_narrow_peaks.filtered.narrowPeak")
    
    neuron_endo_annot <- analyze_peaks(neuron_endo_peaks, "Neuron_Endo")
    neuron_exo_annot <- analyze_peaks(neuron_exo_peaks, "Neuron_Exo")
    nsc_endo_annot <- analyze_peaks(nsc_endo_peaks, "NSC_Endo")
    nsc_exo_annot <- analyze_peaks(nsc_exo_peaks, "NSC_Exo")
    
    # Return annotations for further analysis if needed
    return(list(neuron_endo=neuron_endo_annot, 
                neuron_exo=neuron_exo_annot,
                nsc_endo=nsc_endo_annot,
                nsc_exo=nsc_exo_annot))
  }
  
  # Run genomic distribution analysis
  tryCatch({
    genomic_annotations <- analyze_genomic_distribution()
  }, error = function(e) {
    warning("Error in genomic distribution analysis: ", e$message)
  })
  
  # Return a message indicating completion
  return("Completed generating MeCP2 Endo and Exo heatmaps, metaprofiles, and genomic distribution analysis")
}

# Main execution
main <- function() {
  tryCatch({
    result <- generate_mecp2_tss_heatmaps()
    print(result)
  }, error = function(e) {
    print(paste("Error in main execution:", e$message))
  })
}

main()