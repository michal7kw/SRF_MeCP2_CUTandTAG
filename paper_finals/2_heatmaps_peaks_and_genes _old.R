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
library(R6)
library(reshape2) # Added for matrix manipulation

# Color palette for different sample types
COLORS <- list(
  "NSC_ENDO" = "#83CBEB",
  "NSC_EXO" = "#0070C0",
  "NEU_ENDO" = "#FF9999",
  "NEU_EXO" = "#FF3300"
)

# Enhanced logging system
log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0("[", timestamp, "] [", level, "] ", message, "\n"))
}

# Cache management system
CacheManager <- R6::R6Class("CacheManager",
  public = list(
    cache_dir = NULL,
    force_recompute = FALSE,
    
    initialize = function(cache_dir, force_recompute = FALSE) {
      self$cache_dir <- cache_dir
      self$force_recompute <- force_recompute
      dir.create(self$cache_dir, showWarnings = FALSE, recursive = TRUE)
      log_message(paste("Initialized cache manager at:", cache_dir))
    },
    
    get_cache_path = function(key, extension = "") {
      file.path(self$cache_dir, paste0(digest::digest(key, algo = "sha1"), extension))
    },
    
    exists = function(key, extension = "") {
      cache_path <- self$get_cache_path(key, extension)
      file.exists(cache_path)
    },
    
    save = function(key, extension = "", data = NULL, file_path = NULL) {
      cache_path <- self$get_cache_path(key, extension)
      if (!is.null(file_path)) {
        file.copy(file_path, cache_path, overwrite = TRUE)
      } else if (!is.null(data)) {
        saveRDS(data, cache_path)
      }
      cache_path
    },
    
    load = function(key, extension = "") {
      cache_path <- self$get_cache_path(key, extension)
      if (!file.exists(cache_path)) return(NULL)
      if (grepl("\\.rds$", extension)) {
        readRDS(cache_path)
      } else {
        cache_path
      }
    },
    
    clear = function() {
      unlink(file.path(self$cache_dir, "*"), recursive = TRUE)
      log_message("Cache cleared")
    }
  )
)

# Sample definitions
setwd("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/paper_finals")
OUTPUT_DIR <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/paper_finals/outputs_TSS_TES_R"
BIGWIG_DIR <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results_1b/bigwig"
GENCODE_GTF <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA/gencode.vM25.annotation.gtf"

# Set to TRUE to use GENCODE GTF instead of TxDb
USE_GENCODE <- FALSE

# Create output directory if it doesn't exist
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Set this to TRUE to force recomputation of all matrices and heatmaps
FORCE_RECOMPUTE <- TRUE

# Initialize cache manager
CACHE_DIR <- file.path(OUTPUT_DIR, "cache")
cache_manager <- CacheManager$new(CACHE_DIR, FORCE_RECOMPUTE)

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
      hours <- floor(dur / 3600)
      minutes <- floor((dur %% 3600) / 60)
      seconds <- floor(dur %% 60)
      sprintf("%02d:%02d:%02d", as.integer(hours), as.integer(minutes), as.integer(seconds))
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

# Modified parse_deeptools_matrix function
parse_deeptools_matrix <- function(matrix_file) {
  log_message(paste("Parsing matrix file:", matrix_file))
  
  # Check if this is a deepTools matrix or our R-native matrix
  meta_file <- paste0(matrix_file, ".meta")
  if (file.exists(meta_file)) {
    # This is an R-native matrix
    log_message("Detected R-native matrix format")
    meta_data <- readRDS(meta_file)
    
    # Load the corresponding RDS file
    rds_file <- sub("\\.gz$", ".rds", matrix_file)
    if (file.exists(rds_file)) {
      matrix_data <- readRDS(rds_file)
      
      # Extract matrix dimensions
      num_regions <- nrow(matrix_data)
      num_samples <- length(meta_data$samples)
      
      if (meta_data$use_scale_regions) {
        # For scale-regions matrices
        upstream_bins <- meta_data$upstream / meta_data$bin_size
        body_bins <- meta_data$region_body_length / meta_data$bin_size
        downstream_bins <- meta_data$downstream / meta_data$bin_size
        
        return(list(
          matrix_data = matrix_data,
          is_scale_regions = TRUE,
          upstream_bins = upstream_bins,
          body_bins = body_bins,
          downstream_bins = downstream_bins,
          samples = meta_data$samples,
          num_regions = num_regions
        ))
      } else {
        # For reference-point matrices
        total_bins <- (meta_data$upstream + meta_data$downstream) / meta_data$bin_size
        
        return(list(
          matrix_data = matrix_data,
          is_scale_regions = FALSE,
          total_bins = total_bins,
          samples = meta_data$samples,
          num_regions = num_regions
        ))
      }
    } else {
      log_message(paste("RDS file not found:", rds_file), "ERROR")
      return(NULL)
    }
  }
  
  # Original deepTools parsing code (with error handling)
  tryCatch({
    # Try to use deepTools to extract info
    cmd <- paste0("computeMatrixOperations info -m ", matrix_file)
    info <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    
    if (length(info) == 0 || grepl("error|exception", tolower(paste(info, collapse = " ")))) {
      log_message("Failed to extract info using computeMatrixOperations", "ERROR")
      
      # Fall back to tab file if it exists
      tab_file <- paste0(matrix_file, ".tab")
      if (file.exists(tab_file)) {
        log_message("Trying to parse the tab file directly")
        matrix_data <- read.table(tab_file, header = TRUE, sep = "\t", check.names = FALSE)
        
        # Make a best guess about the structure
        # Typically first 6 columns are metadata, rest are signal values
        meta_cols <- 6
        signal_cols <- ncol(matrix_data) - meta_cols
        
        # Guess if it's scale-regions based on column names
        is_scale_regions <- any(grepl("body|upstream|downstream", colnames(matrix_data)))
        
        if (is_scale_regions) {
          # Guess the bin distribution (this is approximate)
          total_bins <- signal_cols
          upstream_bins <- total_bins / 3
          body_bins <- total_bins / 3
          downstream_bins <- total_bins - upstream_bins - body_bins
          
          return(list(
            matrix_data = matrix_data,
            is_scale_regions = TRUE,
            upstream_bins = upstream_bins,
            body_bins = body_bins,
            downstream_bins = downstream_bins,
            total_bins = total_bins
          ))
        } else {
          return(list(
            matrix_data = matrix_data,
            is_scale_regions = FALSE,
            total_bins = signal_cols
          ))
        }
      } else {
        log_message("Tab file not found, cannot parse matrix", "ERROR")
        return(NULL)
      }
    }
    
    # Continue with the original parsing logic if deepTools worked
    log_message("Matrix information extracted, now parsing")
    
    # Get number of samples, regions, and bins
    samples_line <- grep("samples:", info, value = TRUE)
    regions_line <- grep("regions:", info, value = TRUE)
    bins_line <- grep("bins:", info, value = TRUE)
    
    # Extract numbers
    num_samples <- as.numeric(gsub(".*samples: ([0-9]+).*", "\\1", samples_line))
    num_regions <- as.numeric(gsub(".*regions: ([0-9]+).*", "\\1", regions_line))
    num_bins <- as.numeric(gsub(".*bins: ([0-9]+).*", "\\1", bins_line))
    
    # Determine if it's scale-regions
    is_scale_regions <- any(grepl("body", info))
    
    # Extract the matrix data (if possible)
    tab_file <- paste0(matrix_file, ".tab")
    if (file.exists(tab_file)) {
      matrix_data <- read.table(tab_file, header = TRUE, sep = "\t", check.names = FALSE)
    } else {
      # Try to extract using computeMatrixOperations
      temp_file <- tempfile(fileext = ".tab")
      cmd <- paste0("computeMatrixOperations dump -m ", matrix_file, " -o ", temp_file)
      system(cmd)
      
      if (file.exists(temp_file) && file.size(temp_file) > 0) {
        matrix_data <- read.table(temp_file, header = TRUE, sep = "\t", check.names = FALSE)
        unlink(temp_file)
      } else {
        log_message("Failed to extract matrix data", "ERROR")
        return(NULL)
      }
    }
    
    if (is_scale_regions) {
      # For scale-regions matrices
      # Try to determine the bin distribution
      scale_regions_info <- grep("upstream|body|downstream", info, value = TRUE)
      upstream_bins <- as.numeric(gsub(".*upstream: ([0-9]+).*", "\\1", 
                                      grep("upstream", scale_regions_info, value = TRUE)))
      body_bins <- as.numeric(gsub(".*body: ([0-9]+).*", "\\1", 
                                  grep("body", scale_regions_info, value = TRUE)))
      downstream_bins <- as.numeric(gsub(".*downstream: ([0-9]+).*", "\\1", 
                                        grep("downstream", scale_regions_info, value = TRUE)))
      
      if (is.na(upstream_bins) || is.na(body_bins) || is.na(downstream_bins)) {
        # If we couldn't parse the bin distribution, make a best guess
        upstream_bins <- num_bins / 3
        body_bins <- num_bins / 3
        downstream_bins <- num_bins - upstream_bins - body_bins
      }
      
      return(list(
        matrix_data = matrix_data,
        is_scale_regions = TRUE,
        upstream_bins = upstream_bins,
        body_bins = body_bins,
        downstream_bins = downstream_bins,
        num_samples = num_samples,
        num_regions = num_regions
      ))
    } else {
      # For reference-point matrices
      return(list(
        matrix_data = matrix_data,
        is_scale_regions = FALSE,
        total_bins = num_bins,
        num_samples = num_samples,
        num_regions = num_regions
      ))
    }
  }, error = function(e) {
    log_message(paste("Error parsing matrix file:", e$message), "ERROR")
    return(NULL)
  })
}

# IMPROVED: Function to compute the average profile for each sample
compute_average_profile <- function(matrix_data, sample_idx, 
                                   is_scale_regions = FALSE,
                                   upstream_bins = NULL,
                                   body_bins = NULL,
                                   downstream_bins = NULL,
                                   total_bins = NULL) {
  
  log_message(paste("Computing average profile for sample index", sample_idx))
  
  # Select the appropriate columns for this sample
  if (is_scale_regions) {
    # For scale-regions matrices
    start_col <- 6 + (sample_idx - 1) * (upstream_bins + body_bins + downstream_bins)
    end_col <- start_col + upstream_bins + body_bins + downstream_bins - 1
  } else {
    # For reference-point matrices
    start_col <- 6 + (sample_idx - 1) * total_bins
    end_col <- start_col + total_bins - 1
  }
  
  signal_cols <- start_col:end_col
  
  # Compute the average signal across all regions for each bin
  avg_profile <- colMeans(matrix_data[, signal_cols], na.rm = TRUE)
  
  # Mark TSS and TES positions for scale-regions
  if (is_scale_regions) {
    tss_position <- upstream_bins
    tes_position <- upstream_bins + body_bins
    attr(avg_profile, "tss_position") <- tss_position
    attr(avg_profile, "tes_position") <- tes_position
  }
  
  return(avg_profile)
}

# NEW: Function to verify matrix structure and extract important parameters
verify_matrix_structure <- function(matrix_file) {
  log_message(paste("Verifying matrix structure for", matrix_file))
  
  # Use computeMatrixOperations to get detailed info
  cmd <- paste0("computeMatrixOperations info -m ", matrix_file)
  info <- system(cmd, intern = TRUE)
  
  # Extract key parameters
  group_boundary_info <- grep("group_boundaries", info, value = TRUE)
  sample_info <- grep("samples:", info, value = TRUE)
  region_info <- grep("regions:", info, value = TRUE)
  bin_info <- grep("bins:", info, value = TRUE)
  
  # Parse and return the structure information
  log_message("Matrix structure information:")
  for (line in c(group_boundary_info, sample_info, region_info, bin_info)) {
    log_message(paste("  ", line))
  }
  
  # Determine if it's a scale-regions or reference-point matrix
  is_scale_regions <- any(grepl("body", info))
  
  log_message(paste("Matrix type:", ifelse(is_scale_regions, "scale-regions", "reference-point")))
  
  return(list(
    is_scale_regions = is_scale_regions,
    info = info
  ))
}

# Modified compute_matrix function to fall back to EnrichedHeatmap when deepTools fails
compute_matrix <- function(bw_files, regions, output_file, upstream=5000, downstream=5000, 
                          use_scale_regions=FALSE, force_recompute=FORCE_RECOMPUTE,
                          bin_size=50, region_body_length=5000, sort_regions=TRUE) {
  # Try deepTools approach first
  log_message(paste0("Attempting matrix computation using deepTools..."))
  
  # Create cache key based on input parameters
  cache_key <- list(
    bw_files = bw_files,
    regions = regions,
    upstream = upstream,
    downstream = downstream,
    use_scale_regions = use_scale_regions,
    bin_size = bin_size,
    region_body_length = region_body_length,
    sort_regions = sort_regions,
    modification_times = file.info(c(bw_files, regions))$mtime
  )
  
  # Check cache first
  if (!force_recompute) {
    cached_file <- cache_manager$load(cache_key, ".matrix.gz")
    if (!is.null(cached_file)) {
      log_message("Using cached matrix file")
      file.copy(cached_file, output_file, overwrite = TRUE)
      return(output_file)
    }
  }
  
  # Try deepTools first
  deeptools_success <- FALSE
  tryCatch({
    # Existing deepTools code...
    bw_files_str <- paste(bw_files, collapse = " ")
    
    if (use_scale_regions) {
      cmd <- paste0("computeMatrix scale-regions",
                   " -S ", bw_files_str,
                   " -R ", regions,
                   " -b ", upstream, " -a ", downstream,
                   " --regionBodyLength ", region_body_length,
                   " --binSize ", bin_size,
                   " --sortRegions ", ifelse(sort_regions, "descending", "no"),
                   " --sortUsing mean",
                   " --startLabel 'TSS'",
                   " --endLabel 'TES'",
                   " --skipZeros",
                   " --missingDataAsZero",
                   " -o ", output_file,
                   " --outFileNameMatrix ", paste0(output_file, ".tab"))
    } else {
      cmd <- paste0("computeMatrix reference-point",
                   " -S ", bw_files_str,
                   " -R ", regions,
                   " --referencePoint TSS",
                   " -b ", upstream, " -a ", downstream,
                   " --binSize ", bin_size,
                   " --sortRegions ", ifelse(sort_regions, "descending", "no"),
                   " --sortUsing mean",
                   " --skipZeros",
                   " --missingDataAsZero",
                   " -o ", output_file,
                   " --outFileNameMatrix ", paste0(output_file, ".tab"))
    }
    
    log_message("Executing command:")
    log_message(cmd)
    system_result <- system(cmd)
    
    if (system_result == 0) {
      log_message("Matrix computation with deepTools completed successfully.")
      deeptools_success <- TRUE
    } else {
      log_message(paste("Matrix computation with deepTools failed with exit code:", system_result), "ERROR")
    }
  }, error = function(e) {
    log_message(paste("Error in deepTools execution:", e$message), "ERROR")
  })
  
  # If deepTools failed, use R-native approach with EnrichedHeatmap
  if (!deeptools_success) {
    log_message("Falling back to R-native approach using rtracklayer and EnrichedHeatmap")
    
    # Prepare output directory for R-computed matrix
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Read regions
    regions_gr <- rtracklayer::import(regions)
    
    # Create list to store matrices
    combined_matrix <- NULL
    
    # Process each bigWig file
    for (i in seq_along(bw_files)) {
      bw_file <- bw_files[i]
      log_message(paste("Processing bigWig file:", basename(bw_file)))
      
      # Import bigWig using rtracklayer
      bw_data <- NULL
      tryCatch({
        bw_data <- rtracklayer::import(bw_file, as = "RleList")
      }, error = function(e) {
        log_message(paste("Error importing bigWig file:", e$message), "ERROR")
      })
      
      if (!is.null(bw_data)) {
        # Process depending on mode
        if (use_scale_regions) {
          # For TSS to TES (scale-regions)
          mat <- EnrichedHeatmap::normalizeToMatrix(
            bw_data, 
            regions_gr,
            value_column = "score",
            extend = c(upstream, downstream),
            w = bin_size,
            target_ratio = region_body_length / mean(width(regions_gr)),
            mean_mode = "absolute"
          )
        } else {
          # For TSS-centered (reference-point)
          mat <- EnrichedHeatmap::normalizeToMatrix(
            bw_data,
            regions_gr,
            value_column = "score",
            extend = c(upstream, downstream),
            w = bin_size,
            mean_mode = "absolute"
          )
        }
        
        # Add to combined matrix
        if (is.null(combined_matrix)) {
          combined_matrix <- mat
        } else {
          combined_matrix <- cbind(combined_matrix, mat)
        }
      }
    }
    
    if (!is.null(combined_matrix)) {
      # Convert to data frame and save
      mat_df <- as.data.frame(as.matrix(combined_matrix))
      
      # Add region information
      if (length(regions_gr) == nrow(mat_df)) {
        mat_df$seqnames <- as.character(seqnames(regions_gr))
        mat_df$start <- start(regions_gr)
        mat_df$end <- end(regions_gr)
        mat_df$strand <- as.character(strand(regions_gr))
        
        # If regions have names, add them
        if (!is.null(names(regions_gr))) {
          mat_df$name <- names(regions_gr)
        } else {
          mat_df$name <- paste0("region_", 1:length(regions_gr))
        }
      }
      
      # Save as RDS (native R format)
      rds_file <- sub("\\.gz$", ".rds", output_file)
      saveRDS(mat_df, rds_file)
      log_message(paste("Saved matrix data to:", rds_file))
      
      # Create a tab-delimited version for compatibility
      tab_file <- paste0(output_file, ".tab")
      write.table(mat_df, tab_file, sep = "\t", quote = FALSE, row.names = FALSE)
      log_message(paste("Saved matrix data to tab format:", tab_file))
      
      # Create a dummy deepTools-compatible matrix file
      file.create(output_file)
      log_message(paste("Created placeholder matrix file:", output_file))
      
      # Save custom metadata for later use
      meta_file <- paste0(output_file, ".meta")
      meta_data <- list(
        is_deeptools = FALSE,
        use_scale_regions = use_scale_regions,
        upstream = upstream,
        downstream = downstream,
        bin_size = bin_size,
        region_body_length = region_body_length,
        samples = basename(bw_files)
      )
      saveRDS(meta_data, meta_file)
      
      # Cache the result if needed
      cache_manager$save(cache_key, ".matrix.gz", file_path = output_file)
      
      return(output_file)
    } else {
      log_message("Failed to generate matrix using R-native approach", "ERROR")
      return(NULL)
    }
  }
  
  # Cache the result if deepTools was successful
  if (deeptools_success) {
    cache_manager$save(cache_key, ".matrix.gz", file_path = output_file)
  }
  
  return(output_file)
}

# IMPROVED: Function to plot composite profile directly from matrix data
plot_composite_profile <- function(matrix_file, output_file, title="", 
                                 sample_labels=NULL, colors=NULL, 
                                 ylim=NULL, force_recompute=FORCE_RECOMPUTE) {
  
  # Check if output file already exists
  if (file.exists(output_file) && !force_recompute) {
    log_message(paste("Output file already exists, skipping:", output_file))
    return(output_file)
  } else if (file.exists(output_file) && force_recompute) {
    log_message(paste("Output file exists but force_recompute=TRUE, regenerating:", output_file))
    file.remove(output_file)
  }
  
  log_message(paste("Creating composite profile plot:", output_file))
  
  # Parse the matrix file
  matrix_info <- parse_deeptools_matrix(matrix_file)
  if (is.null(matrix_info)) {
    log_message("Failed to parse matrix file", "ERROR")
    return(NULL)
  }
  
  # Extract necessary information
  matrix_data <- matrix_info$matrix_data
  is_scale_regions <- matrix_info$is_scale_regions
  
  # Get sample count from the matrix file
  cmd <- paste0("computeMatrixOperations info -m ", matrix_file, " | grep 'samples:' | cut -d ':' -f 2")
  sample_count <- as.numeric(trimws(system(cmd, intern = TRUE)))
  log_message(paste("Matrix contains", sample_count, "samples"))
  
  # If sample labels not provided, generate default ones
  if (is.null(sample_labels)) {
    sample_labels <- paste0("Sample ", 1:sample_count)
  }
  
  # If colors not provided, use defaults
  if (is.null(colors)) {
    # Use our predefined colors if we have the right number of samples
    if (sample_count <= length(unlist(COLORS))) {
      colors <- unlist(COLORS)[1:sample_count]
    } else {
      # Generate colors from RColorBrewer
      colors <- brewer.pal(min(9, sample_count), "Set1")
      if (sample_count > 9) {
        # Interpolate if we need more than 9 colors
        colors <- colorRampPalette(colors)(sample_count)
      }
    }
  }
  
  # Compute average profile for each sample
  profiles <- list()
  for (i in 1:sample_count) {
    if (is_scale_regions) {
      profiles[[i]] <- compute_average_profile(matrix_data, i, 
                                             is_scale_regions = TRUE,
                                             upstream_bins = matrix_info$upstream_bins,
                                             body_bins = matrix_info$body_bins,
                                             downstream_bins = matrix_info$downstream_bins)
    } else {
      profiles[[i]] <- compute_average_profile(matrix_data, i, 
                                             is_scale_regions = FALSE,
                                             total_bins = matrix_info$total_bins)
    }
  }
  
  # Create profile data frame for ggplot
  profile_df <- data.frame(
    Position = 1:length(profiles[[1]]),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:length(profiles)) {
    profile_df[[sample_labels[i]]] <- profiles[[i]]
  }
  
  # Melt the data frame for ggplot
  melted_df <- reshape2::melt(profile_df, id.vars = "Position", 
                            variable.name = "Sample", value.name = "Signal")
  
  # Add x-axis labels for scale-regions
  if (is_scale_regions) {
    # Mark TSS and TES positions
    tss_pos <- matrix_info$upstream_bins
    tes_pos <- tss_pos + matrix_info$body_bins
    
    # Create breaks and labels for the x-axis
    x_breaks <- c(1, tss_pos, tes_pos, length(profiles[[1]]))
    x_labels <- c(paste0("-", matrix_info$upstream_bins * 50, "bp"), "TSS", "TES", 
                 paste0("+", matrix_info$downstream_bins * 50, "bp"))
    
    # Create the plot
    p <- ggplot(melted_df, aes(x = Position, y = Signal, color = Sample)) +
      geom_line(size = 1) +
      scale_color_manual(values = colors) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels) +
      labs(title = title, x = "", y = "Signal") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()
      ) +
      geom_vline(xintercept = tss_pos, linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = tes_pos, linetype = "dashed", alpha = 0.5)
    
  } else {
    # For reference-point (TSS-centered) plots
    # Calculate distance from reference point
    bin_size <- 50  # Default bin size
    total_bins <- length(profiles[[1]])
    reference_point_index <- total_bins / 2  # Center point
    
    # Create breaks and labels
    x_breaks <- c(1, total_bins / 4, reference_point_index, 3 * total_bins / 4, total_bins)
    x_labels <- c(
      paste0("-", matrix_info$total_bins * bin_size / 2, "bp"),
      paste0("-", matrix_info$total_bins * bin_size / 4, "bp"),
      "TSS",
      paste0("+", matrix_info$total_bins * bin_size / 4, "bp"),
      paste0("+", matrix_info$total_bins * bin_size / 2, "bp")
    )
    
    # Create the plot
    p <- ggplot(melted_df, aes(x = Position, y = Signal, color = Sample)) +
      geom_line(size = 1) +
      scale_color_manual(values = colors) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels) +
      labs(title = title, x = "Distance from TSS", y = "Signal") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "right",
        panel.grid.minor = element_blank()
      ) +
      geom_vline(xintercept = reference_point_index, linetype = "dashed", alpha = 0.5)
  }
  
  # Set y-axis limits if provided
  if (!is.null(ylim)) {
    p <- p + ylim(ylim)
  }
  
  # Save the plot
  ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  log_message(paste("Profile plot saved to:", output_file))
  
  # Save the underlying data for debugging
  write.csv(melted_df, paste0(output_file, ".csv"), row.names = FALSE)
  
  return(output_file)
}

# Modified plot_heatmap function to handle deepTools failures
plot_heatmap <- function(matrix_file, output_file, title="", zmin=0, zmax=0.3, 
                        use_scale_regions=FALSE, force_recompute=FORCE_RECOMPUTE, 
                        sample_type=NULL) {
  # Check if output file already exists
  if (file.exists(output_file) && !force_recompute) {
    log_message(paste("Output file already exists, skipping:", output_file))
    return(output_file)
  } else if (file.exists(output_file) && force_recompute) {
    log_message(paste("Output file exists but force_recompute=TRUE, regenerating:", output_file))
    file.remove(output_file)
  }
  
  log_message(paste0("Plotting heatmap for ", ifelse(use_scale_regions, "TSS-to-TES", "TSS-centered"), " analysis..."))
  log_message(paste0("Using matrix file: ", basename(matrix_file)))
  log_message(paste0("Output will be saved to: ", basename(output_file)))
  log_message(paste0("Title: ", title))
  
  # Try deepTools plotHeatmap first
  deeptools_success <- FALSE
  tryCatch({
    if (use_scale_regions) {
      # For TSS to TES heatmaps
      cmd <- paste0("plotHeatmap -m ", matrix_file,
                  " -o ", output_file,
                  " --colorMap 'viridis'",
                  " --whatToShow 'heatmap and colorbar'",
                  " --zMin ", zmin, " --zMax ", zmax,
                  " --heatmapHeight 15",
                  " --heatmapWidth 7.5",
                  " --xAxisLabel \"\"",
                  " --regionsLabel \"Genes\"",
                  " --plotTitle \"", title, "\"",
                  " --startLabel \"TSS\"",
                  " --endLabel \"TES\"",
                  " --plotType 'lines'")
    } else {
      # Original TSS-centered heatmaps
      cmd <- paste0("plotHeatmap -m ", matrix_file,
                  " -o ", output_file,
                  " --colorMap 'viridis'",
                  " --whatToShow 'heatmap and colorbar'",
                  " --zMin ", zmin, " --zMax ", zmax,
                  " --heatmapHeight 15",
                  " --heatmapWidth 7.5",
                  " --xAxisLabel \"Distance from TSS (bp)\"",
                  " --refPointLabel \"TSS\"",
                  " --regionsLabel \"Genes\"",
                  " --plotTitle \"", title, "\"",
                  " --plotType 'lines'")
    }
    log_message("Executing command:")
    log_message(cmd)
    system_result <- system(cmd)
    
    if (system_result == 0) {
      log_message("Heatmap plotting completed successfully using deepTools.")
      deeptools_success <- TRUE
    } else {
      log_message(paste("deepTools heatmap plotting failed with exit code:", system_result), "ERROR")
    }
  }, error = function(e) {
    log_message(paste("Error in deepTools execution:", e$message), "ERROR")
  })
  
  # If deepTools failed, generate the heatmap using R's ComplexHeatmap
  if (!deeptools_success) {
    log_message("Falling back to R-native approach using ComplexHeatmap")
    
    # Parse the matrix file
    matrix_info <- parse_deeptools_matrix(matrix_file)
    if (is.null(matrix_info)) {
      log_message("Failed to parse matrix file, cannot create heatmap", "ERROR")
      return(NULL)
    }
    
    # Extract matrix data
    matrix_data <- matrix_info$matrix_data
    
    # Extract the signal columns (typically columns 7+)
    signal_cols <- 7:ncol(matrix_data)
    signal_matrix <- as.matrix(matrix_data[, signal_cols])
    
    # Set up colors
    col_fun <- circlize::colorRamp2(c(zmin, (zmin + zmax)/2, zmax), 
                                   c("blue", "white", "red"))
    
    # Create the heatmap using ComplexHeatmap
    ht <- ComplexHeatmap::Heatmap(
      signal_matrix,
      name = "Signal",
      col = col_fun,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      column_title = title,
      use_raster = TRUE,
      raster_quality = 4
    )
    
    # Save the heatmap to file
    png(output_file, width = 7.5, height = 15, units = "in", res = 300)
    draw(ht)
    dev.off()
    
    log_message(paste("R-native heatmap saved to:", output_file))
  }
  
  return(output_file)
}

# Function to prepare gene regions and TSS files
prepare_annotation_files <- function(temp_dir, force_recompute=FORCE_RECOMPUTE) {
  log_message("=== Preparing annotation files ===")
  
  # Define output files
  gene_regions_file <- file.path(temp_dir, "gene_regions.bed")
  tss_file <- file.path(temp_dir, "tss_regions.bed")
  
  if (USE_GENCODE) {
    log_message("Using GENCODE GTF annotation")
    
    # Check if GENCODE GTF file exists
    if (!file.exists(GENCODE_GTF)) {
      stop("GENCODE GTF file not found: ", GENCODE_GTF)
    }
    
    # Check if files already exist
    if (file.exists(gene_regions_file) && file.exists(tss_file) && !force_recompute) {
      log_message("Annotation files already exist, reusing")
      return(list(gene_regions=gene_regions_file, tss=tss_file))
    } else if ((file.exists(gene_regions_file) || file.exists(tss_file)) && force_recompute) {
      log_message("Annotation files exist but force_recompute=TRUE, regenerating")
      if (file.exists(gene_regions_file)) file.remove(gene_regions_file)
      if (file.exists(tss_file)) file.remove(tss_file)
    }
    
    log_message("Importing GENCODE GTF annotation...")
    # Import GTF file
    gencode <- tryCatch({
      import(GENCODE_GTF)
    }, error = function(e) {
      stop("Error importing GENCODE GTF file: ", e$message)
    })
    
    log_message(paste("Imported", length(gencode), "features from GENCODE"))
    
    # Extract genes
    genes <- gencode[gencode$type == "gene"]
    log_message(paste("Extracted", length(genes), "genes from GENCODE"))
    
    # Export gene regions
    export(genes, gene_regions_file)
    log_message(paste("Gene regions exported to:", gene_regions_file))
    
    # Create TSS regions
    tss_regions <- promoters(genes, upstream=1, downstream=1)
    log_message(paste("Created", length(tss_regions), "TSS regions from GENCODE"))
    
    # Export TSS regions
    export(tss_regions, tss_file)
    log_message(paste("TSS regions exported to:", tss_file))
    
  } else {
    log_message("Using TxDb.Mmusculus.UCSC.mm10.knownGene annotation")
    
    # Check if files already exist
    if (file.exists(gene_regions_file) && file.exists(tss_file) && !force_recompute) {
      log_message("Annotation files already exist, reusing")
      return(list(gene_regions=gene_regions_file, tss=tss_file))
    } else if ((file.exists(gene_regions_file) || file.exists(tss_file)) && force_recompute) {
      log_message("Annotation files exist but force_recompute=TRUE, regenerating")
      if (file.exists(gene_regions_file)) file.remove(gene_regions_file)
      if (file.exists(tss_file)) file.remove(tss_file)
    }
    
    log_message("Loading TxDb for mouse genome...")
    # Load TxDb for mouse genome
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

    log_message("Extracting gene regions...")
    # Extract gene regions
    genes <- genes(txdb)
    log_message(paste("Extracted", length(genes), "gene regions"))
    
    # Export gene regions
    export(genes, gene_regions_file)
    log_message(paste("Gene regions exported to:", gene_regions_file))
    
    log_message("Extracting TSS regions...")
    # Extract TSS regions
    tss_regions <- promoters(genes, upstream=1, downstream=1)
    log_message(paste("Created", length(tss_regions), "TSS regions"))
    
    # Export TSS regions
    export(tss_regions, tss_file)
    log_message(paste("TSS regions exported to:", tss_file))
  }
  
  return(list(gene_regions=gene_regions_file, tss=tss_file))
}

# IMPROVED: Function to generate heatmaps and metaprofiles for Mecp2 Endo and Exo around TSS
generate_mecp2_tss_heatmaps <- function(temp_dir, annotation_files) {
  tracker <- ProgressTracker$new(18)
  
  # Define sample groups
  sample_groups <- list(
    neuron_endo = c("NeuM2", "NeuM3"),
    neuron_exo = c("NeuV1", "NeuV2", "NeuV3"),
    nsc_endo = c("NSCM1", "NSCM2", "NSCM3"),
    nsc_exo = c("NSCv1", "NSCv2", "NSCv3")
  )
  
  # NEW: Sample labels for plots
  sample_labels <- list(
    neuron_endo = "Neuron Endogenous",
    neuron_exo = "Neuron Exogenous",
    nsc_endo = "NSC Endogenous",
    nsc_exo = "NSC Exogenous"
  )
  
  # Pre-compute all matrices
  matrices <- list()
  for (group_name in names(sample_groups)) {
    bw_files <- file.path(BIGWIG_DIR, paste0(sample_groups[[group_name]], ".bw"))
    
    # TSS matrices
    matrices[[paste0(group_name, "_tss")]] <- compute_matrix(
      bw_files,
      annotation_files$tss,
      file.path(temp_dir, paste0(group_name, "_tss_matrix.gz")),
      bin_size = 50,  # Explicitly set bin size
      sort_regions = TRUE  # Sort regions by signal strength
    )
    tracker$increment()
    
    # TSS to TES matrices with improved parameters
    matrices[[paste0(group_name, "_tss_tes")]] <- compute_matrix(
      bw_files,
      annotation_files$gene_regions,
      file.path(temp_dir, paste0(group_name, "_tss_tes_matrix.gz")),
      use_scale_regions = TRUE,
      bin_size = 50,  # Explicitly set bin size
      region_body_length = 5000,  # Fixed body length for better scaling
      sort_regions = TRUE  # Sort regions by signal strength
    )
    tracker$increment()
  }
  
  # Generate individual plots for each group
  for (group_name in names(sample_groups)) {
    # Get sample type for colors
    if (group_name == "neuron_endo") {
      sample_type <- "NEU_ENDO"
    } else if (group_name == "neuron_exo") {
      sample_type <- "NEU_EXO"
    } else if (group_name == "nsc_endo") {
      sample_type <- "NSC_ENDO"
    } else {
      sample_type <- "NSC_EXO"
    }
    
    # TSS to TES heatmap
    plot_heatmap(
      matrices[[paste0(group_name, "_tss_tes")]],
      file.path(OUTPUT_DIR, paste0(group_name, "_tss_tes_heatmap.png")),
      title = paste(sample_labels[[group_name]], "MeCP2 gene body coverage"),
      use_scale_regions = TRUE,
      sample_type = sample_type
    )
    tracker$increment()
    
    # TSS heatmap
    plot_heatmap(
      matrices[[paste0(group_name, "_tss")]],
      file.path(OUTPUT_DIR, paste0(group_name, "_tss_heatmap.png")),
      title = paste(sample_labels[[group_name]], "MeCP2 around TSS"),
      sample_type = sample_type
    )
    tracker$increment()
    
    # IMPROVED: Custom TSS profile with accurate bin positioning
    plot_composite_profile(
      matrices[[paste0(group_name, "_tss")]],
      file.path(OUTPUT_DIR, paste0(group_name, "_tss_profile.png")),
      title = paste(sample_labels[[group_name]], "MeCP2 around TSS"),
      sample_labels = c(sample_labels[[group_name]]),
      colors = c(COLORS[[sample_type]])
    )
    tracker$increment()
    
    # IMPROVED: Custom gene body profile with accurate bin positioning
    plot_composite_profile(
      matrices[[paste0(group_name, "_tss_tes")]],
      file.path(OUTPUT_DIR, paste0(group_name, "_tss_tes_profile.png")),
      title = paste(sample_labels[[group_name]], "MeCP2 gene body coverage"),
      sample_labels = c(sample_labels[[group_name]]),
      colors = c(COLORS[[sample_type]])
    )
    tracker$increment()
  }
  
  # Generate combined plots for each cell type
  generate_combined_plots <- function(cell_type) {
    log_message(paste("Generating combined plots for", cell_type))
    
    endo_group <- paste0(cell_type, "_endo")
    exo_group <- paste0(cell_type, "_exo")
    endo_type <- paste0(toupper(substring(cell_type, 1, 3)), "_ENDO")
    exo_type <- paste0(toupper(substring(cell_type, 1, 3)), "_EXO")
    
    # IMPROVED: Combine existing matrices rather than creating new ones
    # Create a combined TSS matrix (reference-point) from the individual sample bigWig files
    combined_files_tss <- c(
      file.path(BIGWIG_DIR, paste0(sample_groups[[endo_group]], ".bw")),
      file.path(BIGWIG_DIR, paste0(sample_groups[[exo_group]], ".bw"))
    )
    
    combined_tss_matrix <- compute_matrix(
      combined_files_tss,
      annotation_files$tss,
      file.path(temp_dir, paste0(cell_type, "_combined_tss_matrix.gz")),
      bin_size = 50,
      sort_regions = TRUE
    )
    tracker$increment()
    
    # Create a combined gene body (scale-regions) matrix
    combined_tss_tes_matrix <- compute_matrix(
      combined_files_tss,
      annotation_files$gene_regions,
      file.path(temp_dir, paste0(cell_type, "_combined_tss_tes_matrix.gz")),
      use_scale_regions = TRUE,
      bin_size = 50,
      region_body_length = 5000,
      sort_regions = TRUE
    )
    tracker$increment()
    
    # Get sample labels for the combined plots
    all_samples <- c(
      sample_groups[[endo_group]],
      sample_groups[[exo_group]]
    )
    
    sample_labels <- c(
      rep(paste(cell_type, "Endogenous"), length(sample_groups[[endo_group]])),
      rep(paste(cell_type, "Exogenous"), length(sample_groups[[exo_group]]))
    )
    
    # Get colors for the combined plots
    colors <- c(
      rep(COLORS[[endo_type]], length(sample_groups[[endo_group]])),
      rep(COLORS[[exo_type]], length(sample_groups[[exo_group]]))
    )
    
    # Create simplified labels for the composite plots
    simple_labels <- c(paste(cell_type, "Endogenous"), paste(cell_type, "Exogenous"))
    simple_colors <- c(COLORS[[endo_type]], COLORS[[exo_type]])
    
    # Plot combined profiles
    # TSS profile
    plot_composite_profile(
      combined_tss_matrix,
      file.path(OUTPUT_DIR, paste0(cell_type, "_endo_vs_exo_tss_profile.png")),
      title = paste(cell_type, "Endogenous vs Exogenous MeCP2 around TSS"),
      sample_labels = simple_labels,
      colors = simple_colors
    )
    tracker$increment()
    
    # Gene body profile
    plot_composite_profile(
      combined_tss_tes_matrix,
      file.path(OUTPUT_DIR, paste0(cell_type, "_endo_vs_exo_tss_tes_profile.png")),
      title = paste(cell_type, "Endogenous vs Exogenous MeCP2 gene body coverage"),
      sample_labels = simple_labels,
      colors = simple_colors
    )
    tracker$increment()
    
    # Combined heatmaps
    plot_heatmap(
      combined_tss_matrix,
      file.path(OUTPUT_DIR, paste0(cell_type, "_endo_vs_exo_tss_heatmap.png")),
      title = paste(cell_type, "Endogenous vs Exogenous MeCP2 around TSS")
    )
    tracker$increment()
    
    plot_heatmap(
      combined_tss_tes_matrix,
      file.path(OUTPUT_DIR, paste0(cell_type, "_endo_vs_exo_tss_tes_heatmap.png")),
      title = paste(cell_type, "Endogenous vs Exogenous MeCP2 gene body coverage"),
      use_scale_regions = TRUE
    )
    tracker$increment()
  }
  
  # Generate combined plots for both cell types
  generate_combined_plots("neuron")
  generate_combined_plots("nsc")
  
  # NEW: Generate diagnostic plots to compare TSS and TES signal heights
  generate_diagnostic_plots <- function() {
    log_message("Generating diagnostic plots to verify TSS and TES signal heights")
    
    # Create a data frame to store TSS and TES signal values
    diagnostics <- data.frame(
      Sample = character(),
      TSS_Signal = numeric(),
      TES_Signal = numeric(),
      TSS_TES_Ratio = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Extract signal values for each sample group
    for (group_name in names(sample_groups)) {
      # Parse the matrix to get the profile
      matrix_file <- matrices[[paste0(group_name, "_tss_tes")]]
      matrix_info <- parse_deeptools_matrix(matrix_file)
      
      if (!is.null(matrix_info) && matrix_info$is_scale_regions) {
        # Get the profile for this group
        profile <- compute_average_profile(
          matrix_info$matrix_data, 
          1,  # First sample in the matrix
          is_scale_regions = TRUE,
          upstream_bins = matrix_info$upstream_bins,
          body_bins = matrix_info$body_bins,
          downstream_bins = matrix_info$downstream_bins
        )
        
        # Get TSS and TES signal values
        tss_idx <- attr(profile, "tss_position")
        tes_idx <- attr(profile, "tes_position")
        
        tss_signal <- profile[tss_idx]
        tes_signal <- profile[tes_idx]
        
        # Add to diagnostics data frame
        diagnostics <- rbind(diagnostics, data.frame(
          Sample = group_name,
          TSS_Signal = tss_signal,
          TES_Signal = tes_signal,
          TSS_TES_Ratio = tss_signal / tes_signal,
          stringsAsFactors = FALSE
        ))
        
        # Create a diagnostic plot showing the raw signal
        raw_signal_df <- data.frame(
          Position = 1:length(profile),
          Signal = profile
        )
        
        p <- ggplot(raw_signal_df, aes(x = Position, y = Signal)) +
          geom_line(color = COLORS[[ifelse(group_name == "neuron_endo", "NEU_ENDO",
                                        ifelse(group_name == "neuron_exo", "NEU_EXO",
                                              ifelse(group_name == "nsc_endo", "NSC_ENDO", "NSC_EXO")))]]) +
          geom_vline(xintercept = tss_idx, color = "red", linetype = "dashed") +
          geom_vline(xintercept = tes_idx, color = "green", linetype = "dashed") +
          annotate("text", x = tss_idx, y = max(profile) * 0.9, label = "TSS", color = "red") +
          annotate("text", x = tes_idx, y = max(profile) * 0.8, label = "TES", color = "green") +
          labs(title = paste(group_name, "Raw Signal"),
               subtitle = paste("TSS/TES Ratio:", round(tss_signal / tes_signal, 2)),
               x = "Position (bins)", y = "Signal") +
          theme_minimal()
        
        ggsave(file.path(OUTPUT_DIR, paste0(group_name, "_raw_signal_diagnostic.png")), p, width = 10, height = 6, dpi = 300)
      }
    }
    
    # Save diagnostics to CSV
    write.csv(diagnostics, file.path(OUTPUT_DIR, "tss_tes_diagnostics.csv"), row.names = FALSE)
    
    # Create a comparison plot of TSS/TES ratios
    p <- ggplot(diagnostics, aes(x = Sample, y = TSS_TES_Ratio, fill = Sample)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = round(TSS_TES_Ratio, 2)), vjust = -0.5) +
      labs(title = "TSS/TES Signal Ratio Comparison",
           x = "", y = "TSS/TES Ratio") +
      scale_fill_manual(values = c(
        "neuron_endo" = COLORS$NEU_ENDO,
        "neuron_exo" = COLORS$NEU_EXO,
        "nsc_endo" = COLORS$NSC_ENDO,
        "nsc_exo" = COLORS$NSC_EXO
      )) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(OUTPUT_DIR, "tss_tes_ratio_comparison.png"), p, width = 10, height = 6, dpi = 300)
  }
  
  # Run the diagnostic function
  generate_diagnostic_plots()
  
  log_message("=== Analysis completed ===")
}

# Main execution
main <- function() {
  timer <- Timer$new()
  
  # Central directory management
  temp_dir <- file.path(OUTPUT_DIR, "temp")
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  tryCatch({
    log_message("=== Starting MeCP2 analysis pipeline ===")
    log_message(paste("System time:", Sys.time()))
    log_message(paste("R version:", R.version.string))
    log_message(paste("Platform:", R.version$platform))
    
    log_message("Loading required packages")
    print(sessionInfo())
    
    # Prepare files once and pass to downstream functions
    annotation_files <- prepare_annotation_files(temp_dir)
    
    # IMPROVED: Run our enhanced analysis function
    log_message("Generating heatmaps and profiles with improved methodology")
    result <- generate_mecp2_tss_heatmaps(temp_dir, annotation_files)
    
    log_message(paste("Completed successfully in", timer$format_duration()))
    return(result)
  },
  error = function(e) {
    log_message(paste("Fatal error:", e$message), "ERROR")
    log_message(paste("Stack trace:", paste(deparse(e$call), collapse = "\n")), "ERROR")
    log_message(paste("Execution failed after", timer$format_duration()), "ERROR")
    stop(e)
  })
}

main()
main()