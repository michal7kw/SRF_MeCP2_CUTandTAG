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
library(patchwork) # For combining plots

# Color palette for different sample types
COLORS <- list(
  "NSC_ENDO" = "#83CBEB",
  "NSC_EXO" = "#0070C0",
  "NEU_ENDO" = "#FF9999",
  "NEU_EXO" = "#FF3300"
)

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
FORCE_RECOMPUTE <- FALSE

#' Format duration in seconds to HH:MM:SS
format_duration <- function(seconds) {
  # Ensure seconds is numeric and convert to integer
  seconds <- as.integer(seconds)
  hours <- seconds %/% 3600
  minutes <- (seconds %% 3600) %/% 60
  secs <- seconds %% 60
  sprintf("%02d:%02d:%02d", hours, minutes, secs)
}

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
    initialize = function(total_steps) {
      self$total_steps <- total_steps
      self$current_step <- 0
      log_message(paste("Initialized progress tracker with", total_steps, "steps"))
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

# Constants and configuration
CACHE_DIR <- file.path(dirname(OUTPUT_DIR), "cache")
MATRIX_CACHE_PATTERN <- "matrix_%s_%s.gz"  # Format: matrix_[type]_[condition].gz

#' Enhanced cache management
CacheManager <- R6::R6Class("CacheManager",
  public = list(
    cache_dir = NULL,
    initialize = function(cache_dir = CACHE_DIR) {
      self$cache_dir <- cache_dir
      dir.create(self$cache_dir, recursive = TRUE, showWarnings = FALSE)
    },
    
    get_cache_path = function(key, type) {
      file.path(self$cache_dir, sprintf(MATRIX_CACHE_PATTERN, type, key))
    },
    
    exists = function(key, type) {
      file.exists(self$get_cache_path(key, type))
    },
    
    is_valid = function(key, type, source_files) {
      cache_file <- self$get_cache_path(key, type)
      if (!file.exists(cache_file)) return(FALSE)
      
      cache_mtime <- file.info(cache_file)$mtime
      source_mtimes <- sapply(source_files, function(f) file.info(f)$mtime)
      all(cache_mtime > source_mtimes)
    }
  )
)

#' Improved matrix computation with caching and dimension validation
compute_matrix <- function(bw_files, regions, output_file, cache_manager = NULL,
                         upstream = 3000, downstream = 3000, bin_size = 50,
                         reference_point = "center", force_recompute = FALSE) {
  
  # Ensure cache_manager exists
  if (is.null(cache_manager)) {
    cache_manager <- CacheManager$new()
  }
  
  # Generate cache key based on input parameters
  cache_key <- digest::digest(list(
    bw_files = bw_files,
    regions = regions,
    upstream = upstream,
    downstream = downstream,
    reference_point = reference_point,
    bin_size = bin_size
  ))
  
  # Check cache validity
  if (!force_recompute && 
      cache_manager$is_valid(cache_key, "matrix", c(bw_files, regions))) {
    cached_file <- cache_manager$get_cache_path(cache_key, "matrix")
    log_message(sprintf("Using cached matrix: %s", basename(cached_file)))
    file.copy(cached_file, output_file, overwrite = TRUE)
    return(output_file)
  }
  
  # Calculate expected number of bins
  total_length <- upstream + downstream
  expected_bins <- total_length / bin_size
  log_message(sprintf("Expected bins: %d (upstream: %d, downstream: %d, bin_size: %d)",
                     expected_bins, upstream, downstream, bin_size))
  
  # Compute matrix if not cached or force recompute
  log_message(sprintf("Computing matrix for %d bigwig files...", length(bw_files)))
  
  # Split computation for large datasets
  chunk_size <- 5  # Process 5 bigwig files at a time
  chunks <- split(bw_files, ceiling(seq_along(bw_files)/chunk_size))
  
  temp_matrices <- vector("character", length(chunks))
  
  for (i in seq_along(chunks)) {
    chunk_files <- chunks[[i]]
    temp_output <- tempfile(pattern = "matrix_chunk", fileext = ".gz")
    
    cmd <- sprintf(
      "computeMatrix reference-point -S %s -R %s --referencePoint %s -b %d -a %d --binSize %d -o %s --skipZeros --missingDataAsZero --numberOfProcessors %d",
      paste(chunk_files, collapse = " "),
      regions,
      reference_point,
      upstream,
      downstream,
      bin_size,
      temp_output,
      parallel::detectCores() - 1
    )
    
    system_result <- system(cmd)
    if (system_result != 0) {
      stop(sprintf("Matrix computation failed for chunk %d", i))
    }
    
    # Validate matrix dimensions
    cmd_info <- sprintf("computeMatrixOperations info -m %s", temp_output)
    matrix_info <- system(cmd_info, intern = TRUE)
    matrix_shape <- grep("shape", matrix_info, value = TRUE)
    current_bins <- as.numeric(strsplit(gsub(".*shape: \\((\\d+), (\\d+)\\).*", "\\2", matrix_shape), " ")[[1]])
    
    if (current_bins != expected_bins) {
      stop(sprintf("Matrix dimension mismatch in chunk %d. Expected %d bins, got %d bins",
                  i, expected_bins, current_bins))
    }
    
    temp_matrices[i] <- temp_output
  }
  
  # Merge chunks if necessary
  if (length(temp_matrices) > 1) {
    log_message("Merging matrix chunks...")
    cmd_merge <- sprintf(
      "computeMatrixOperations rbind -m %s -o %s",
      paste(temp_matrices, collapse = " "),
      output_file
    )
    merge_result <- system(cmd_merge)
    if (merge_result != 0) {
      stop("Failed to merge matrix chunks")
    }
  } else {
    file.copy(temp_matrices[1], output_file, overwrite = TRUE)
  }
  
  # Validate final matrix dimensions
  cmd_info <- sprintf("computeMatrixOperations info -m %s", output_file)
  matrix_info <- system(cmd_info, intern = TRUE)
  matrix_shape <- grep("shape", matrix_info, value = TRUE)
  final_bins <- as.numeric(strsplit(gsub(".*shape: \\((\\d+), (\\d+)\\).*", "\\2", matrix_shape), " ")[[1]])
  
  if (final_bins != expected_bins) {
    stop(sprintf("Final matrix dimension mismatch. Expected %d bins, got %d bins",
                expected_bins, final_bins))
  }
  
  # Cache the result
  file.copy(output_file, cache_manager$get_cache_path(cache_key, "matrix"), 
            overwrite = TRUE)
  
  # Cleanup
  unlink(temp_matrices)
  
  output_file
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

# Function to merge bigwig files
merge_bigwigs <- function(bw_files, output_file, force_recompute=FORCE_RECOMPUTE) {
  # Check if output file already exists
  if (file.exists(output_file) && !force_recompute) {
    log_message(paste("Merged bigwig file already exists, skipping:", output_file))
    return(output_file)
  } else if (file.exists(output_file) && force_recompute) {
    log_message(paste("Merged bigwig file exists but force_recompute=TRUE, regenerating:", output_file))
    # Remove the existing file to ensure clean regeneration
    file.remove(output_file)
  }
  
  log_message("Merging bigwig files...")
  log_message(paste0("Input files: ", paste(basename(bw_files), collapse=", ")))
  log_message(paste0("Output will be saved to: ", basename(output_file)))
  
  bw_files_str <- paste(bw_files, collapse = " ")
  
  # Use wiggletools to merge bigwig files (mean)
  cmd <- paste0("wiggletools mean ", bw_files_str, " | wigToBigWig stdin mm10.chrom.sizes ", output_file)
  
  log_message("Executing command:")
  log_message(cmd)
  system_result <- system(cmd)
  
  if (system_result == 0) {
    log_message("Bigwig merging completed successfully.")
  } else {
    log_message(paste("Bigwig merging failed with exit code:", system_result), "ERROR")
    # Try alternative method using deepTools
    log_message("Trying alternative method using deepTools...")
    cmd_alt <- paste0("bigwigCompare -b1 ", bw_files[1], " -b2 ", bw_files[2], 
                     " --operation mean -o ", output_file, 
                     " --binSize 50 --verbose")
    log_message("Executing command:")
    log_message(cmd_alt)
    system_result_alt <- system(cmd_alt)
    
    if (system_result_alt == 0) {
      log_message("Bigwig merging completed successfully using deepTools.")
    } else {
      log_message(paste("Bigwig merging failed with both methods. Exit code:", system_result_alt), "ERROR")
      stop("Failed to merge bigwig files")
    }
  }
  
  return(output_file)
}

# Function to create combined profile plot
create_combined_profile_plot <- function(endo_matrix, exo_matrix, output_file, title="MeCP2 at CpG Islands", cell_type="", force_recompute=FORCE_RECOMPUTE) {
  # Check if output file already exists
  if (file.exists(output_file) && !force_recompute) {
    log_message(paste("Combined profile plot already exists, skipping:", output_file))
    return(output_file)
  } else if (file.exists(output_file) && force_recompute) {
    log_message(paste("Combined profile plot exists but force_recompute=TRUE, regenerating:", output_file))
    file.remove(output_file)
  }
  
  log_message("Creating combined profile plot...")
  log_message(paste0("Using matrices: ", basename(endo_matrix), " and ", basename(exo_matrix)))
  log_message(paste0("Output will be saved to: ", basename(output_file)))
  
  # Determine colors based on cell type
  endo_color <- if(cell_type == "neuron") COLORS$NEU_ENDO else COLORS$NSC_ENDO
  exo_color <- if(cell_type == "neuron") COLORS$NEU_EXO else COLORS$NSC_EXO
  
  # Create temporary profile files
  temp_dir <- dirname(output_file)
  endo_profile_temp <- file.path(temp_dir, paste0("temp_endo_profile_", cell_type, ".png"))
  exo_profile_temp <- file.path(temp_dir, paste0("temp_exo_profile_", cell_type, ".png"))
  
  # Generate individual profiles
  cmd_endo <- paste0("plotProfile -m ", endo_matrix,
                    " -o ", endo_profile_temp,
                    " --colors ", endo_color,
                    " --plotTitle \"Endogenous MeCP2\"",
                    " --refPointLabel \"CpG Island Center\"",
                    " --regionsLabel \"CpG Islands\"",
                    " --yAxisLabel \"Signal\"",
                    " --startLabel \"-3kb\"",
                    " --endLabel \"+3kb\"")
  
  cmd_exo <- paste0("plotProfile -m ", exo_matrix,
                   " -o ", exo_profile_temp,
                   " --colors ", exo_color,
                   " --plotTitle \"Exogenous MeCP2\"",
                   " --refPointLabel \"CpG Island Center\"",
                   " --regionsLabel \"CpG Islands\"",
                   " --yAxisLabel \"Signal\"",
                   " --startLabel \"-3kb\"",
                   " --endLabel \"+3kb\"")
  
  log_message("Generating endogenous profile:")
  log_message(cmd_endo)
  system(cmd_endo)
  
  log_message("Generating exogenous profile:")
  log_message(cmd_exo)
  system(cmd_exo)
  
  # Combine the profiles
  cmd_combine <- paste0("montage ", endo_profile_temp, " ", exo_profile_temp, 
                       " -geometry +0+0 -tile 2x1 ", output_file)
  
  log_message("Combining profiles:")
  log_message(cmd_combine)
  system_result <- system(cmd_combine)
  
  if (system_result == 0) {
    log_message("Combined profile plot created successfully.")
    # Clean up temporary files
    file.remove(endo_profile_temp, exo_profile_temp)
  } else {
    log_message(paste("Failed to create combined profile plot. Exit code:", system_result), "ERROR")
  }
  
  return(output_file)
}

# Function to generate merged CpG island profiles for Mecp2 Endo and Exo
generate_merged_mecp2_cpg_profiles <- function(temp_dir, cpg_islands_file) {
  tracker <- ProgressTracker$new(8)
  
  log_message("=== Starting merged MeCP2 CpG island profile generation ===")
  tracker$increment()
  
  # Define sample groups
  sample_groups <- list(
    neuron_endo = c("NeuM2", "NeuM3"),
    neuron_exo = c("NeuV1", "NeuV2", "NeuV3"),
    nsc_endo = c("NSCM1", "NSCM2", "NSCM3"),
    nsc_exo = c("NSCv1", "NSCv2", "NSCv3")
  )
  
  # Create a function to get bigwig files for a group
  get_bigwig_files <- function(group) {
    file.path(BIGWIG_DIR, paste0(sample_groups[[group]], ".bw"))
  }
  
  # Define consistent parameters for all matrices
  upstream_val <- 3000
  downstream_val <- 3000
  bin_size_val <- 50
  cache_manager <- CacheManager$new()
  
  # Process each cell type
  cell_types <- c("neuron", "nsc")
  for (cell_type in cell_types) {
    endo_group <- paste0(cell_type, "_endo")
    exo_group <- paste0(cell_type, "_exo")
    
    # Get bigwig files
    endo_bw <- get_bigwig_files(endo_group)
    exo_bw <- get_bigwig_files(exo_group)
    
    # Compute matrices
    endo_matrix <- file.path(temp_dir, paste0(cell_type, "_endo_cpg_matrix.gz"))
    exo_matrix <- file.path(temp_dir, paste0(cell_type, "_exo_cpg_matrix.gz"))
    
    compute_matrix(endo_bw, cpg_islands_file, endo_matrix, cache_manager, 
                  upstream = upstream_val, downstream = downstream_val, bin_size = bin_size_val)
    tracker$increment()
    compute_matrix(exo_bw, cpg_islands_file, exo_matrix, cache_manager,
                  upstream = upstream_val, downstream = downstream_val, bin_size = bin_size_val)
    tracker$increment()
    
    # Create combined profile
    combined_profile <- file.path(OUTPUT_DIR, paste0(cell_type, "_endo_exo_cpg_profile.png"))
    create_combined_profile_plot(
      endo_matrix, 
      exo_matrix, 
      combined_profile,
      paste(toupper(substr(cell_type, 1, 1)), substr(cell_type, 2, nchar(cell_type)), "MeCP2 at CpG Islands"),
      cell_type
    )
    tracker$increment()
  }
  
  # Create all-samples profile
  all_samples_matrix <- file.path(temp_dir, "all_samples_cpg_matrix.gz")
  all_bw_files <- c(
    get_bigwig_files("neuron_endo")[1],
    get_bigwig_files("neuron_exo")[1],
    get_bigwig_files("nsc_endo")[1],
    get_bigwig_files("nsc_exo")[1]
  )
  
  compute_matrix(all_bw_files, cpg_islands_file, all_samples_matrix, cache_manager,
                upstream = upstream_val, downstream = downstream_val, bin_size = bin_size_val)
  tracker$increment()
  
  # Create final combined profile with all samples
  all_samples_profile <- file.path(OUTPUT_DIR, "all_samples_cpg_profile.png")
  plot_profile(
    all_samples_matrix,
    all_samples_profile,
    "MeCP2 at CpG Islands",
    sample_types = c("NEU_ENDO", "NEU_EXO", "NSC_ENDO", "NSC_EXO")
  )
  tracker$increment()
  
  log_message("=== Merged MeCP2 CpG island profile generation completed ===")
}

# Add memory management utilities
gc_and_report <- function() {
  gc_result <- gc(full = TRUE)
  mem_used <- gc_result[2, 2] # Used memory in MB
  log_message(sprintf("Memory after GC: %.2f GB", mem_used/1024))
}

#' Memory-efficient plot generation
generate_plots <- function(matrices, output_dir, tracker) {
  # Clean up memory before starting
  gc_and_report()
  
  # First merge the matrices
  merged_matrix <- file.path(dirname(matrices[[1]]), "merged_matrix.gz")
  
  # Check if all matrices exist and have the same dimensions
  log_message("Checking matrix dimensions before merging...")
  
  for (matrix_file in matrices) {
    if (!file.exists(matrix_file)) {
      log_message(sprintf("Matrix file does not exist: %s", matrix_file), "ERROR")
      stop("Missing matrix file")
    }
    
    # Get matrix info using computeMatrixOperations info
    cmd_info <- sprintf("computeMatrixOperations info -m %s", matrix_file)
    matrix_info <- system(cmd_info, intern = TRUE)
    log_message(sprintf("Matrix %s info: %s", basename(matrix_file), 
                       paste(matrix_info[grep("shape", matrix_info)], collapse=", ")))
  }
  
  # Use computeMatrixOperations to merge matrices
  cmd_merge <- sprintf(
    "computeMatrixOperations rbind -m %s -o %s",
    paste(matrices, collapse=" "),
    merged_matrix
  )
  
  log_message("Merging matrices...")
  log_message(sprintf("Command: %s", cmd_merge))
  
  result_merge <- system(cmd_merge)
  
  if(result_merge != 0) {
    log_message("Failed to merge matrices", "ERROR")
    stop("Matrix merging failed")
  }
  
  # Now create plot using the merged matrix
  cmd_plot <- sprintf(
    'plotProfile --verbose -m %s -o %s --colors "%s" "%s" "%s" "%s" --plotTitle "MeCP2 at CpG Islands" --refPointLabel "CpG Island Center" --regionsLabel "CpG Islands" --yAxisLabel "Signal" --samplesLabel "Neuron Endo" "Neuron Exo" "NSC Endo" "NSC Exo" --legendLocation "upper-right"',
    merged_matrix,
    file.path(output_dir, "combined_profile.png"),
    COLORS$NEU_ENDO,
    COLORS$NEU_EXO,
    COLORS$NSC_ENDO,
    COLORS$NSC_EXO
  )
  
  log_message("Generating plot...")
  log_message(sprintf("Command: %s", cmd_plot))
  
  result_plot <- system(cmd_plot)
  
  if(result_plot == 0) {
    log_message("Successfully generated combined profile plot")
    tracker$increment()
  } else {
    # Try with simpler options if the first attempt fails
    log_message("First attempt failed, trying with simplified options...", "WARNING")
    
    cmd_simple <- sprintf(
      'plotProfile --verbose -m %s -o %s --plotTitle "MeCP2 at CpG Islands" --refPointLabel "CpG Island Center" --regionsLabel "CpG Islands" --yAxisLabel "Signal"',
      merged_matrix,
      file.path(output_dir, "combined_profile.png")
    )
    
    log_message(sprintf("Trying alternative command: %s", cmd_simple))
    
    result_simple <- system(cmd_simple)
    
    if(result_simple == 0) {
      log_message("Successfully generated combined profile plot with simplified options")
      tracker$increment()
    } else {
      log_message("Both plotting attempts failed", "ERROR")
      log_message("Checking environment and matrix file...", "WARNING")
      
      # Check if matrix file exists and is readable
      if(file.exists(merged_matrix)) {
        log_message(sprintf("Merged matrix file exists: %s", merged_matrix))
        if(file.access(merged_matrix, mode=4) == 0) {
          log_message("Merged matrix file is readable")
        } else {
          log_message("Merged matrix file is not readable", "WARNING")
        }
      } else {
        log_message("Merged matrix file does not exist", "ERROR")
      }
      
      # Check Python environment
      python_version <- system("python --version", intern=TRUE)
      log_message(sprintf("Python version: %s", python_version))
      
      # Check deepTools version
      deeptools_version <- system("plotProfile --version", intern=TRUE)
      log_message(sprintf("deepTools version: %s", deeptools_version))
      
      stop("Plot generation failed. Check logs for details.")
    }
  }
  
  # Clean up merged matrix file
  if(file.exists(merged_matrix)) {
    unlink(merged_matrix)
  }
  
  # Clean up memory after plotting
  gc_and_report()
}

#' Main execution function with progress tracking
main <- function() {
  start_time <- Sys.time()
  cache_manager <- CacheManager$new()
  tracker <- ProgressTracker$new(total_steps = 8)
  
  tryCatch({
    # Create temporary directory
    temp_dir <- tempfile("heatmap_analysis_")
    dir.create(temp_dir)
    on.exit(unlink(temp_dir, recursive = TRUE))
    
    # Define sample groups
    sample_groups <- list(
      neuron_endo = c("NeuM2", "NeuM3"),
      neuron_exo = c("NeuV1", "NeuV2", "NeuV3"),
      nsc_endo = c("NSCM1", "NSCM2", "NSCM3"),
      nsc_exo = c("NSCv1", "NSCv2", "NSCv3")
    )
    
    # Function to get bigwig files for a group
    get_bigwig_files <- function(group) {
      file.path(BIGWIG_DIR, paste0(sample_groups[[group]], ".bw"))
    }
    
    # Process each condition with proper caching
    conditions <- list(
      neuron_endo = get_bigwig_files("neuron_endo"),
      neuron_exo = get_bigwig_files("neuron_exo"),
      nsc_endo = get_bigwig_files("nsc_endo"),
      nsc_exo = get_bigwig_files("nsc_exo")
    )
    
    # Define consistent parameters for all matrices
    upstream_val <- 3000
    downstream_val <- 3000
    bin_size_val <- 50
    
    matrices <- lapply(names(conditions), function(condition) {
      output_file <- file.path(temp_dir, sprintf("%s_matrix.gz", condition))
      compute_matrix(
        conditions[[condition]], 
        CPG_ISLANDS_BED,
        output_file,
        cache_manager,
        upstream = upstream_val,
        downstream = downstream_val,
        bin_size = bin_size_val,
        force_recompute = FORCE_RECOMPUTE
      )
      tracker$increment()
      output_file
    })
    
    # Generate final plots
    generate_plots(matrices, OUTPUT_DIR, tracker)
    
    # Report completion
    end_time <- Sys.time()
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
    log_message(sprintf("Analysis completed in %s", format_duration(duration)))
    
  }, error = function(e) {
    log_message(sprintf("Fatal error: %s", e$message), "ERROR")
    log_message(sprintf("Stack trace: %s", paste(deparse(e$call), collapse = "\n")), "ERROR")
    quit(status = 1)
  })
}

main()
