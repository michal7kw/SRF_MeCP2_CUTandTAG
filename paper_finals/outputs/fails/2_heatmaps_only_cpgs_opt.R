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

# Define sample groups and their corresponding bigwig files
SAMPLE_GROUPS <- list(
  neuron_endo = c("NeuM2", "NeuM3"),
  neuron_exo = c("NeuV1", "NeuV2", "NeuV3"),
  nsc_endo = c("NSCM1", "NSCM2", "NSCM3"),
  nsc_exo = c("NSCv1", "NSCv2", "NSCv3")
)

# Function to get bigwig files for a group
get_bigwig_files <- function(group) {
  file.path(BIGWIG_DIR, paste0(SAMPLE_GROUPS[[group]], ".bw"))
}

# Define bigwig file paths
neuron_endo_bw <- get_bigwig_files("neuron_endo")
neuron_exo_bw <- get_bigwig_files("neuron_exo")
nsc_endo_bw <- get_bigwig_files("nsc_endo")
nsc_exo_bw <- get_bigwig_files("nsc_exo")

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

#' Improved matrix computation with caching
compute_matrix <- function(bw_files, regions, output_file, cache_manager,
                         upstream = 3000, downstream = 3000,
                         reference_point = "center", force_recompute = FALSE) {
  
  # Generate cache key based on input parameters
  cache_key <- digest::digest(list(
    bw_files = bw_files,
    regions = regions,
    upstream = upstream,
    downstream = downstream,
    reference_point = reference_point
  ))
  
  # Check cache validity
  if (!force_recompute && 
      cache_manager$is_valid(cache_key, "matrix", c(bw_files, regions))) {
    cached_file <- cache_manager$get_cache_path(cache_key, "matrix")
    log_message(sprintf("Using cached matrix: %s", basename(cached_file)))
    file.copy(cached_file, output_file, overwrite = TRUE)
    return(output_file)
  }
  
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
      "computeMatrix reference-point -S %s -R %s --referencePoint %s -b %d -a %d -o %s --skipZeros --missingDataAsZero --numberOfProcessors %d",
      paste(chunk_files, collapse = " "),
      regions,
      reference_point,
      upstream,
      downstream,
      temp_output,
      parallel::detectCores() - 1
    )
    
    system_result <- system(cmd)
    if (system_result != 0) {
      stop(sprintf("Matrix computation failed for chunk %d", i))
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
    system(cmd_merge)
  } else {
    file.copy(temp_matrices[1], output_file, overwrite = TRUE)
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
    # Remove the existing file to ensure clean regeneration
    file.remove(output_file)
  }
  
  log_message("Creating combined profile plot...")
  log_message(paste0("Using matrices: ", basename(endo_matrix), " and ", basename(exo_matrix)))
  log_message(paste0("Output will be saved to: ", basename(output_file)))
  
  # Create temporary profile files
  temp_dir <- dirname(output_file)
  endo_profile_temp <- file.path(temp_dir, paste0("temp_endo_profile_", cell_type, ".png"))
  exo_profile_temp <- file.path(temp_dir, paste0("temp_exo_profile_", cell_type, ".png"))
  
  # Generate individual profiles
  cmd_endo <- paste0("plotProfile -m ", endo_matrix,
                    " -o ", endo_profile_temp,
                    " --colors blue",
                    " --plotTitle 'Endogenous MeCP2'",
                    " --refPointLabel 'CpG Island Center'",
                    " --regionsLabel 'CpG Islands'",
                    " --yAxisLabel 'Signal'",
                    " --startLabel -3kb",
                    " --endLabel +3kb")
  
  cmd_exo <- paste0("plotProfile -m ", exo_matrix,
                   " -o ", exo_profile_temp,
                   " --colors red",
                   " --plotTitle 'Exogenous MeCP2'",
                   " --refPointLabel 'CpG Island Center'",
                   " --regionsLabel 'CpG Islands'",
                   " --yAxisLabel 'Signal'",
                   " --startLabel -3kb",
                   " --endLabel +3kb")
  
  log_message("Generating endogenous profile:")
  log_message(cmd_endo)
  system(cmd_endo)
  
  log_message("Generating exogenous profile:")
  log_message(cmd_exo)
  system(cmd_exo)
  
  # Now combine the two profiles into one figure using ImageMagick
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
    # Try alternative method using R's grid package
    log_message("Trying alternative method using R...")
    
    # Create a combined plot using R's grid package
    tryCatch({
      # Create a new PNG device
      png(output_file, width=1200, height=600)
      
      # Create a layout with 1 row and 2 columns
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(1, 2)))
      
      # Add a title
      grid.text(title, y=0.95, gp=gpar(fontsize=16, fontface="bold"))
      
      # Add the images
      if (file.exists(endo_profile_temp) && file.exists(exo_profile_temp)) {
        endo_img <- png::readPNG(endo_profile_temp)
        exo_img <- png::readPNG(exo_profile_temp)
        
        # Display the images
        grid.raster(endo_img, vp=viewport(layout.pos.row=1, layout.pos.col=1))
        grid.raster(exo_img, vp=viewport(layout.pos.row=1, layout.pos.col=2))
      } else {
        grid.text("Failed to generate profile images", vp=viewport(layout.pos.row=1, layout.pos.col=1:2))
      }
      
      # Close the device
      dev.off()
      
      # Clean up temporary files
      if (file.exists(endo_profile_temp)) file.remove(endo_profile_temp)
      if (file.exists(exo_profile_temp)) file.remove(exo_profile_temp)
      
      log_message("Combined profile plot created successfully using R.")
    }, error = function(e) {
      log_message(paste("Failed to create combined profile plot using R:", e$message), "ERROR")
    })
  }
  
  return(output_file)
}

# Function to generate merged CpG island profiles for Mecp2 Endo and Exo
generate_merged_mecp2_cpg_profiles <- function(temp_dir, cpg_islands_file) {
  tracker <- ProgressTracker$new(total_steps = 8)
  
  # Check if all output files already exist
  output_files <- c(
    file.path(OUTPUT_DIR, "neuron_endo_exo_cpg_profile.png"),
    file.path(OUTPUT_DIR, "nsc_endo_exo_cpg_profile.png"),
    file.path(OUTPUT_DIR, "all_samples_cpg_profile.png")
  )
  
  if (all(file.exists(output_files)) && !FORCE_RECOMPUTE) {
    log_message("All profile plots already exist and FORCE_RECOMPUTE is FALSE. Skipping generation.")
    return(invisible(NULL))
  }
  
  log_message("=== Starting merged MeCP2 CpG island profile generation ===")
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
  
  # Create merged bigwig files
  log_message("=== Creating merged bigwig files ===")
  
  # Create a directory for merged bigwigs
  merged_bw_dir <- file.path(temp_dir, "merged_bigwigs")
  dir.create(merged_bw_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get bigwig file paths
  neuron_endo_bw <- file.path(BIGWIG_DIR, paste0(neuron_endo, ".bw"))
  neuron_exo_bw <- file.path(BIGWIG_DIR, paste0(neuron_exo, ".bw"))
  nsc_endo_bw <- file.path(BIGWIG_DIR, paste0(nsc_endo, ".bw"))
  nsc_exo_bw <- file.path(BIGWIG_DIR, paste0(nsc_exo, ".bw"))
  
  # Create merged bigwig files
  log_message("=== Merging bigwig files ===")
  
  # First, create a chromosome sizes file for mm10
  chrom_sizes_file <- file.path(temp_dir, "mm10.chrom.sizes")
  if (!file.exists(chrom_sizes_file)) {
    log_message("Creating mm10 chromosome sizes file...")
    # Use fetchChromSizes from UCSC tools or create manually
    chrom_sizes <- data.frame(
      chrom = paste0("chr", c(1:19, "X", "Y", "M")),
      size = c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 
              129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 
              104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698, 16299)
    )
    write.table(chrom_sizes, chrom_sizes_file, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    log_message("Created mm10 chromosome sizes file")
  }
  
  # Use computeMatrix directly on multiple bigwig files instead of merging
  # This approach is more reliable than trying to merge bigwigs
  
  # Compute matrices for CpG islands
  log_message("=== Computing CpG island matrices ===")
  
  # Matrices for Neurons
  neuron_endo_matrix <- file.path(temp_dir, "neuron_endo_cpg_matrix.gz")
  neuron_exo_matrix <- file.path(temp_dir, "neuron_exo_cpg_matrix.gz")
  
  compute_matrix(neuron_endo_bw, cpg_islands_file, neuron_endo_matrix)
  tracker$increment()
  compute_matrix(neuron_exo_bw, cpg_islands_file, neuron_exo_matrix)
  tracker$increment()
  
  # Matrices for NSCs
  nsc_endo_matrix <- file.path(temp_dir, "nsc_endo_cpg_matrix.gz")
  nsc_exo_matrix <- file.path(temp_dir, "nsc_exo_cpg_matrix.gz")
  
  compute_matrix(nsc_endo_bw, cpg_islands_file, nsc_endo_matrix)
  tracker$increment()
  compute_matrix(nsc_exo_bw, cpg_islands_file, nsc_exo_matrix)
  tracker$increment()
  
  # Create combined profile plots
  log_message("=== Creating combined profile plots ===")
  
  # Combined plot for Neurons
  neuron_combined_profile <- file.path(OUTPUT_DIR, "neuron_endo_exo_cpg_profile.png")
  create_combined_profile_plot(neuron_endo_matrix, neuron_exo_matrix, neuron_combined_profile, 
                              "Neuron MeCP2 at CpG Islands", "neuron")
  tracker$increment()
  
  # Combined plot for NSCs
  nsc_combined_profile <- file.path(OUTPUT_DIR, "nsc_endo_exo_cpg_profile.png")
  create_combined_profile_plot(nsc_endo_matrix, nsc_exo_matrix, nsc_combined_profile, 
                              "NSC MeCP2 at CpG Islands", "nsc")
  tracker$increment()
  
  # Create a combined plot for all samples
  log_message("=== Creating all-in-one profile plot ===")
  
  # Create a matrix with all samples
  all_samples_matrix <- file.path(temp_dir, "all_samples_cpg_matrix.gz")
  all_bw_files <- c(neuron_endo_bw, neuron_exo_bw, nsc_endo_bw, nsc_exo_bw)
  compute_matrix(all_bw_files, cpg_islands_file, all_samples_matrix)
  
  # Create a profile plot with all samples
  all_samples_profile <- file.path(OUTPUT_DIR, "all_samples_cpg_profile.png")
  cmd_all <- sprintf(
    "plotProfile -m %s -o %s --colors blue,red,green,orange --plotTitle 'MeCP2 at CpG Islands' --refPointLabel 'CpG Island Center' --regionsLabel 'CpG Islands' --yAxisLabel Signal --startLabel -3kb --endLabel +3kb --samplesLabel Neuron_Endo,Neuron_Exo,NSC_Endo,NSC_Exo --legendLocation upper-right",
    all_samples_matrix,
    all_samples_profile
  )
  
  log_message("Generating all-samples profile:")
  log_message(cmd_all)
  system(cmd_all)
  tracker$increment()
  
  log_message("=== Merged MeCP2 CpG island profile generation completed ===")
}

#' Format duration in seconds to HH:MM:SS
format_duration <- function(seconds) {
  # Ensure seconds is numeric and convert to integer
  seconds <- as.integer(seconds)
  hours <- seconds %/% 3600
  minutes <- (seconds %% 3600) %/% 60
  secs <- seconds %% 60
  sprintf("%02d:%02d:%02d", hours, minutes, secs)
}

#' Main execution function with progress tracking
main <- function() {
  start_time <- Sys.time()
  cache_manager <- CacheManager$new()
  tracker <- ProgressTracker$new(total_steps = 8)
  
  tryCatch({
    # Check if final outputs already exist
    final_outputs <- c(
      file.path(OUTPUT_DIR, "neuron_endo_exo_cpg_profile.png"),
      file.path(OUTPUT_DIR, "nsc_endo_exo_cpg_profile.png"),
      file.path(OUTPUT_DIR, "all_samples_cpg_profile.png"),
      file.path(OUTPUT_DIR, "combined_profile.png")
    )
    
    # Check if all outputs exist and FORCE_RECOMPUTE is FALSE
    if (all(file.exists(final_outputs)) && !FORCE_RECOMPUTE) {
      log_message("All output files already exist and FORCE_RECOMPUTE is FALSE. Skipping processing.")
      return(invisible(NULL))
    }
    
    # Create temporary directory
    temp_dir <- tempfile("heatmap_analysis_")
    dir.create(temp_dir)
    on.exit(unlink(temp_dir, recursive = TRUE))
    
    # Prepare CpG islands file first
    log_message("Preparing CpG islands file...")
    cpg_islands_file <- prepare_cpg_islands_file(temp_dir, force_recompute=FORCE_RECOMPUTE)
    
    # Process each condition with proper caching
    conditions <- list(
      neuron_endo = neuron_endo_bw,
      neuron_exo = neuron_exo_bw,
      nsc_endo = nsc_endo_bw,
      nsc_exo = nsc_exo_bw
    )
    
    matrices <- lapply(names(conditions), function(condition) {
      output_file <- file.path(temp_dir, sprintf("%s_matrix.gz", condition))
      if (file.exists(output_file) && !FORCE_RECOMPUTE) {
        log_message(sprintf("Matrix for %s already exists, skipping computation", condition))
        tracker$increment()
        return(output_file)
      }
      compute_matrix(
        conditions[[condition]], 
        cpg_islands_file,
        output_file,
        cache_manager,
        force_recompute = FORCE_RECOMPUTE
      )
      tracker$increment()
      output_file
    })
    
    # Generate final plots only if they don't exist or FORCE_RECOMPUTE is TRUE
    if (!all(file.exists(final_outputs)) || FORCE_RECOMPUTE) {
      generate_plots(matrices, OUTPUT_DIR, tracker)
    } else {
      log_message("All plots already exist, skipping plot generation")
      tracker$increment()
    }
    
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

#' Memory-efficient plot generation
generate_plots <- function(matrices, output_dir, tracker) {
  output_file <- file.path(output_dir, "combined_profile.png")
  
  if (file.exists(output_file) && !FORCE_RECOMPUTE) {
    log_message("Combined profile plot already exists and FORCE_RECOMPUTE is FALSE. Skipping generation.")
    tracker$increment()
    return(invisible(NULL))
  }
  
  # Generate combined plot
  cmd <- sprintf(
    "plotProfile -m %s -o %s --colors blue,red,green,orange --plotTitle 'MeCP2 at CpG Islands' --refPointLabel 'CpG Island Center' --regionsLabel 'CpG Islands' --yAxisLabel Signal --startLabel -3kb --endLabel +3kb --samplesLabel Neuron_Endo,Neuron_Exo,NSC_Endo,NSC_Exo --legendLocation upper-right",
    paste(matrices, collapse = " "),
    output_file
  )
  
  log_message("Executing plotProfile command:")
  log_message(cmd)
  system(cmd)
  tracker$increment()
}

main()
