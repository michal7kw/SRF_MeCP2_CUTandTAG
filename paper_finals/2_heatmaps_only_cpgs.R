#!/usr/bin/env Rscript

# Load required libraries
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
library(R6)

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
CPG_ISLANDS_BED <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA/cpg_islands.bed"

# Create output directory if it doesn't exist
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Set this to TRUE to force recomputation of all matrices
FORCE_RECOMPUTE <- TRUE

# Enhanced logging system
log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0("[", timestamp, "] [", level, "] ", message, "\n"))
}

# Execution timer class
Timer <- R6::R6Class("Timer",
  public = list(
    start_time = NULL,
    initialize = function() {
      self$start_time <- Sys.time()
      log_message("Execution timer initialized")
    },
    elapsed = function() {
      difftime(Sys.time(), self$start_time, units = "secs")
    },
    format_duration = function() {
      dur <- as.numeric(self$elapsed())
      sprintf("%02d:%02d:%02d", dur %/% 3600, (dur %% 3600) %/% 60, dur %% 60)
    }
  )
)

# Memory monitoring
get_memory_usage <- function() {
  tryCatch({
    mem <- system("grep VmRSS /proc/$PPID/status | awk '{print $2}'", intern = TRUE)
    paste0(round(as.numeric(mem)/1024, 1), " GB")
  }, error = function(e) "N/A")
}

# Progress tracking
ProgressTracker <- R6::R6Class("ProgressTracker",
  public = list(
    total_steps = 0,
    current_step = 0,
    initialize = function(total) {
      self$total_steps <- total
      log_message(paste("Initialized progress tracker with", total, "steps"))
    },
    increment = function() {
      self$current_step <- self$current_step + 1
      mem_usage <- get_memory_usage()
      log_message(sprintf("Progress: %d/%d (%.1f%%) | Memory: %s",
                        self$current_step, self$total_steps,
                        (self$current_step/self$total_steps)*100,
                        mem_usage))
    }
  )
)

# Function to compute matrix using deepTools
compute_matrix <- function(bw_files, regions, output_file, upstream=3000, downstream=3000, reference_point="center", force_recompute=FORCE_RECOMPUTE) {
  # Check if output file already exists
  if (file.exists(output_file) && !force_recompute) {
    log_message(paste("Matrix file already exists, skipping:", output_file))
    return(output_file)
  } else if (file.exists(output_file) && force_recompute) {
    log_message(paste("Matrix file exists but force_recompute=TRUE, regenerating:", output_file))
    # Remove the existing file to ensure clean regeneration
    file.remove(output_file)
  }
  
  log_message("Computing matrix for CpG islands analysis...")
  log_message(paste0("Using bigwig files: ", paste(basename(bw_files), collapse=", ")))
  log_message(paste0("Using regions file: ", basename(regions)))
  log_message(paste0("Output will be saved to: ", basename(output_file)))
  
  bw_files_str <- paste(bw_files, collapse = " ")
  
  # Reference-point approach for CpG islands
  cmd <- paste0("computeMatrix reference-point",
               " -S ", bw_files_str,
               " -R ", regions,
               " --referencePoint ", reference_point,
               " -b ", upstream, " -a ", downstream,
               " -o ", output_file,
               " --skipZeros",
               " --missingDataAsZero")
  
  log_message("Executing command:")
  log_message(cmd)
  system_result <- system(cmd)
  
  if (system_result == 0) {
    log_message("Matrix computation completed successfully.")
  } else {
    log_message(paste("Matrix computation failed with exit code:", system_result), "ERROR")
  }
  
  return(output_file)
}

# Function to prepare CpG islands file
prepare_cpg_islands_file <- function(temp_dir, force_recompute=FORCE_RECOMPUTE) {
  log_message("=== Preparing CpG islands file ===")
  
  # Define output file
  cpg_islands_file <- file.path(temp_dir, "cpg_islands.bed")
  
  # Check if file already exists
  if (file.exists(cpg_islands_file) && !force_recompute) {
    log_message("CpG islands file already exists, reusing")
    return(cpg_islands_file)
  } else if (file.exists(cpg_islands_file) && force_recompute) {
    log_message("CpG islands file exists but force_recompute=TRUE, regenerating")
    file.remove(cpg_islands_file)
  }
  
  # Check if source file exists
  if (!file.exists(CPG_ISLANDS_BED)) {
    stop("CpG islands BED file not found: ", CPG_ISLANDS_BED)
  }
  
  log_message("Importing CpG islands from BED file...")
  
  # Read the custom format BED file manually
  tryCatch({
    # Read the file with read.table to handle the custom format
    cpg_data <- read.table(CPG_ISLANDS_BED, sep="\t", stringsAsFactors=FALSE, 
                          col.names=c("chrom", "start", "end", "score", "cpg_label", "cpg_count"))
    
    log_message(paste("Read", nrow(cpg_data), "CpG islands from file"))
    
    # Create a GRanges object from the data
    cpg_islands <- GRanges(
      seqnames = cpg_data$chrom,
      ranges = IRanges(start = cpg_data$start, end = cpg_data$end),
      score = cpg_data$score,
      cpg_count = cpg_data$cpg_count
    )
    
    log_message(paste("Created GRanges object with", length(cpg_islands), "CpG islands"))
    
  }, error = function(e) {
    stop("Error processing CpG islands BED file: ", e$message)
  })
  
  # Export CpG islands in standard BED format
  export(cpg_islands, cpg_islands_file)
  log_message(paste("CpG islands exported to:", cpg_islands_file))
  
  return(cpg_islands_file)
}

# Function to plot heatmap using deepTools
plot_heatmap <- function(matrix_file, output_file, title="", use_scale_regions=FALSE, zmin=0, zmax=20, force_recompute=FORCE_RECOMPUTE) {
  # Check if output file already exists
  if (file.exists(output_file) && !force_recompute) {
    log_message(paste("Output file already exists, skipping:", output_file))
    return(output_file)
  } else if (file.exists(output_file) && force_recompute) {
    log_message(paste("Output file exists but force_recompute=TRUE, regenerating:", output_file))
    # Remove the existing file to ensure clean regeneration
    file.remove(output_file)
  }
  
  log_message("Plotting heatmap...")
  log_message(paste0("Using matrix file: ", basename(matrix_file)))
  log_message(paste0("Output will be saved to: ", basename(output_file)))
  log_message(paste0("Title: ", as.character(title)))
  
  if (use_scale_regions) {
    # For TSS to TES heatmaps
    cmd <- paste0("plotHeatmap -m ", matrix_file,
                 " -o ", output_file,
                 " --colorMap 'Blues'",
                 " --whatToShow 'heatmap and colorbar'",
                 " --zMin ", zmin, " --zMax ", zmax,
                 " --heatmapHeight 15",
                 " --heatmapWidth 7.5",
                 " --xAxisLabel \"\"",
                 " --regionsLabel \"Genes\"",
                 " --plotTitle \"", title, "\"",
                 " --startLabel \"TSS\"",
                 " --endLabel \"TES\"")
  } else {
    # Original TSS-centered heatmaps
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
  }
  log_message("Executing command:")
  log_message(cmd)
  system_result <- system(cmd)
  
  if (system_result == 0) {
    log_message("Heatmap plotting completed successfully.")
  } else {
    log_message(paste("Heatmap plotting failed with exit code:", system_result), "ERROR")
  }
  
  return(output_file)
}

# Function to plot profile using deepTools
plot_profile <- function(matrix_file, output_file, title="", force_recompute=FORCE_RECOMPUTE) {
  # Check if output file already exists
  if (file.exists(output_file) && !force_recompute) {
    log_message(paste("Output file already exists, skipping:", output_file))
    return(output_file)
  } else if (file.exists(output_file) && force_recompute) {
    log_message(paste("Output file exists but force_recompute=TRUE, regenerating:", output_file))
    # Remove the existing file to ensure clean regeneration
    file.remove(output_file)
  }
  
  log_message("Plotting profile...")
  log_message(paste0("Using matrix file: ", basename(matrix_file)))
  log_message(paste0("Output will be saved to: ", basename(output_file)))
  log_message(paste0("Title: ", as.character(title)))
  
  # Get the number of samples in the matrix
  # Use computeMatrixOperations to get info about the matrix
  temp_info_file <- tempfile()
  cmd_info <- paste0("computeMatrixOperations info -m ", matrix_file, " > ", temp_info_file)
  log_message("Getting matrix info:")
  log_message(cmd_info)
  system(cmd_info)
  
  # Read the first few lines to determine number of samples
  info_lines <- readLines(temp_info_file, n=20)
  file.remove(temp_info_file)
  
  # Find the line with sample info
  sample_line <- grep("sample_labels", info_lines, value=TRUE)
  if (length(sample_line) > 0) {
    # Extract sample count
    sample_count <- length(strsplit(gsub(".*\\[|\\].*", "", sample_line), ",")[[1]])
    log_message(paste("Detected", sample_count, "samples in matrix"))
  } else {
    # Default to 2 if we can't determine
    sample_count <- 2
    log_message("Could not determine sample count, defaulting to 2")
  }
  
  # Generate enough colors
  colors <- paste(rep(c("blue", "red", "green", "orange", "purple", "brown"), length.out=sample_count), collapse=" ")
  
  cmd <- paste0("plotProfile -m ", matrix_file,
               " -o ", output_file,
               " --colors ", colors,
               " --plotTitle \"", as.character(title), "\"",
               " --refPointLabel \"CpG Island Center\"",
               " --regionsLabel \"CpG Islands\"",
               " --yAxisLabel \"Signal\"")
  
  log_message("Executing command:")
  log_message(cmd)
  system_result <- system(cmd)
  
  if (system_result == 0) {
    log_message("Profile plotting completed successfully.")
  } else {
    log_message(paste("Profile plotting failed with exit code:", system_result), "ERROR")
  }
  
  return(output_file)
}

# Function to generate CpG island profiles for Mecp2 Endo and Exo
generate_mecp2_cpg_profiles <- function(temp_dir, cpg_islands_file) {
  tracker <- ProgressTracker$new(8)  # Total number of main steps
  
  log_message("=== Starting MeCP2 CpG island profile generation ===")
  tracker$increment()
  
  # Define sample groups
  endo_samples <- c("NeuM2", "NeuM3", "NSCM1", "NSCM2", "NSCM3")
  exo_samples <- c("NeuV1", "NeuV2", "NeuV3", "NSCv1", "NSCv2", "NSCv3")
  
  log_message(paste("Endogenous samples:", paste(endo_samples, collapse=", ")))
  log_message(paste("Exogenous samples:", paste(exo_samples, collapse=", ")))
  
  # Group by cell type and MeCP2 type
  neuron_endo <- c("NeuM2", "NeuM3")
  nsc_endo <- c("NSCM1", "NSCM2", "NSCM3")
  neuron_exo <- c("NeuV1", "NeuV2", "NeuV3")
  nsc_exo <- c("NSCv1", "NSCv2", "NSCv3")
  
  # Create temporary matrix files
  log_message(paste("Temporary directory:", temp_dir))
  
  # Generate matrices for Neuron Endogenous vs Exogenous
  log_message("=== Processing Neuron Endogenous vs Exogenous ===")
  neuron_endo_bw <- file.path(BIGWIG_DIR, paste0(neuron_endo, ".bw"))
  neuron_exo_bw <- file.path(BIGWIG_DIR, paste0(neuron_exo, ".bw"))
  
  log_message(paste("Neuron Endogenous bigwig files:", paste(basename(neuron_endo_bw), collapse=", ")))
  log_message(paste("Neuron Exogenous bigwig files:", paste(basename(neuron_exo_bw), collapse=", ")))
  
  # Combine bigwig files for each group
  neuron_endo_matrix <- file.path(temp_dir, "neuron_endo_cpg_matrix.gz")
  neuron_exo_matrix <- file.path(temp_dir, "neuron_exo_cpg_matrix.gz")
  
  # Compute matrices for CpG islands
  log_message("=== Computing CpG island matrices for Neurons ===")
  compute_matrix(neuron_endo_bw, cpg_islands_file, neuron_endo_matrix)
  tracker$increment()
  compute_matrix(neuron_exo_bw, cpg_islands_file, neuron_exo_matrix)
  tracker$increment()
  
  # Plot CpG island profiles for Neurons
  log_message("=== Plotting CpG island profiles for Neurons ===")
  neuron_endo_profile <- file.path(OUTPUT_DIR, "neuron_endo_cpg_profile.png")
  neuron_exo_profile <- file.path(OUTPUT_DIR, "neuron_exo_cpg_profile.png")
  
  plot_profile(neuron_endo_matrix, neuron_endo_profile, "Neuron Endogenous MeCP2 at CpG Islands")
  tracker$increment()
  plot_profile(neuron_exo_matrix, neuron_exo_profile, "Neuron Exogenous MeCP2 at CpG Islands")
  tracker$increment()
  
  # Generate matrices for NSC Endogenous vs Exogenous
  log_message("=== Processing NSC Endogenous vs Exogenous ===")
  nsc_endo_bw <- file.path(BIGWIG_DIR, paste0(nsc_endo, ".bw"))
  nsc_exo_bw <- file.path(BIGWIG_DIR, paste0(nsc_exo, ".bw"))
  
  log_message(paste("NSC Endogenous bigwig files:", paste(basename(nsc_endo_bw), collapse=", ")))
  log_message(paste("NSC Exogenous bigwig files:", paste(basename(nsc_exo_bw), collapse=", ")))
  
  # Combine bigwig files for each group
  nsc_endo_matrix <- file.path(temp_dir, "nsc_endo_cpg_matrix.gz")
  nsc_exo_matrix <- file.path(temp_dir, "nsc_exo_cpg_matrix.gz")
  
  # Compute matrices for CpG islands
  log_message("=== Computing CpG island matrices for NSCs ===")
  compute_matrix(nsc_endo_bw, cpg_islands_file, nsc_endo_matrix)
  tracker$increment()
  compute_matrix(nsc_exo_bw, cpg_islands_file, nsc_exo_matrix)
  tracker$increment()
  
  # Plot CpG island profiles for NSCs
  log_message("=== Plotting CpG island profiles for NSCs ===")
  nsc_endo_profile <- file.path(OUTPUT_DIR, "nsc_endo_cpg_profile.png")
  nsc_exo_profile <- file.path(OUTPUT_DIR, "nsc_exo_cpg_profile.png")
  
  plot_profile(nsc_endo_matrix, nsc_endo_profile, "NSC Endogenous MeCP2 at CpG Islands")
  tracker$increment()
  plot_profile(nsc_exo_matrix, nsc_exo_profile, "NSC Exogenous MeCP2 at CpG Islands")
  tracker$increment()
  
  log_message("=== MeCP2 CpG island profile generation completed ===")
}

# Main execution
main <- function() {
  timer <- Timer$new()
  
  # Central directory management
  temp_dir <- file.path(OUTPUT_DIR, "temp")
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  tryCatch({
    log_message("=== Starting MeCP2 CpG islands analysis pipeline ===")
    log_message(paste("System time:", Sys.time()))
    log_message(paste("R version:", R.version.string))
    log_message(paste("Platform:", R.version$platform))
    
    log_message("Loading required packages")
    print(sessionInfo())
    
    # Prepare CpG islands file
    cpg_islands_file <- prepare_cpg_islands_file(temp_dir)
    
    log_message("Generating CpG island profiles")
    generate_mecp2_cpg_profiles(temp_dir, cpg_islands_file)
    
    log_message(paste("Completed successfully in", timer$format_duration()))
    log_message("Profile plots have been generated in the outputs directory")
  },
  error = function(e) {
    log_message(paste("Fatal error:", e$message), "ERROR")
    log_message(paste("Stack trace:", paste(deparse(e$call), collapse = "\n")), "ERROR")
    log_message(paste("Execution failed after", timer$format_duration()), "ERROR")
    stop(e)
  })
}

main()