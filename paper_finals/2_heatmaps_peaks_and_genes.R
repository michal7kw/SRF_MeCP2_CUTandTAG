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

# Modified compute_matrix function to handle strand issues properly
compute_matrix <- function(bw_files, regions, output_file, upstream=5000, downstream=5000, 
                          use_scale_regions=FALSE, force_recompute=FORCE_RECOMPUTE,
                          bin_size=50, region_body_length=5000, sort_regions=TRUE) {
  
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
    cached_file <- cache_manager$load(cache_key, ".rds")
    if (!is.null(cached_file)) {
      log_message("Using cached matrix data")
      return(cached_file)
    }
  }
  
  log_message(paste0("Computing matrix for ", ifelse(use_scale_regions, "TSS-to-TES", "TSS-centered"), " analysis..."))
  log_message(paste0("Using bigwig files: ", paste(basename(bw_files), collapse=", ")))
  log_message(paste0("Using regions file: ", basename(regions)))
  log_message(paste0("Output will be saved to: ", basename(output_file)))
  
  # Import regions as GRanges object
  regions_gr <- rtracklayer::import(regions)
  log_message(paste("Imported", length(regions_gr), "regions from", basename(regions)))
  
  # Ensure strand is properly set (if missing, set to '*')
  if (!("strand" %in% names(mcols(regions_gr)))) {
    strand(regions_gr) <- "*"
    log_message("Added default strand information to regions")
  }
  
  # Store results for each bigWig file
  matrix_list <- list()
  sample_names <- basename(bw_files)
  
  has_valid_matrices <- FALSE
  
  for (i in seq_along(bw_files)) {
    bw_file <- bw_files[i]
    sample_name <- sample_names[i]
    log_message(paste("Processing", sample_name))
    
    # Import bigWig file with proper error handling
    tryCatch({
      # Import as RleList and set options to avoid strand problems
      # The key fix here: import.bw without strand information
      bw_data <- rtracklayer::import.bw(bw_file, as = "RleList")
      
      # Generate matrix based on mode
      if (use_scale_regions) {
        # For TSS to TES (scale-regions mode)
        # First make sure we're using a non-stranded view of the RleList
        matrix <- tryCatch({
          # Make a copy of regions_gr without strand information for normalization
          noStrand_regions <- regions_gr
          strand(noStrand_regions) <- "*"
          
          # Use the non-stranded regions
          EnrichedHeatmap::normalizeToMatrix(
            bw_data,
            noStrand_regions,
            extend = c(upstream, downstream),
            target_ratio = region_body_length / mean(width(regions_gr)),
            w = bin_size,
            mean_mode = "absolute"
          )
        }, error = function(e) {
          log_message(paste("Error in scale-regions matrix generation:", e$message), "ERROR")
          # Try an alternative method without target_ratio that might be causing problems
          # Again, use non-stranded version
          noStrand_regions <- regions_gr
          strand(noStrand_regions) <- "*"
          
          EnrichedHeatmap::normalizeToMatrix(
            bw_data,
            noStrand_regions,
            extend = c(upstream, downstream),
            w = bin_size,
            mean_mode = "absolute"
          )
        })
      } else {
        # For TSS-centered (reference-point mode)
        # Create a non-stranded version of the regions
        noStrand_regions <- GRanges(
          seqnames = seqnames(regions_gr),
          ranges = IRanges(start = start(regions_gr), width = 1),
          strand = "*"  # Set all strands to unspecified
        )
        
        matrix <- EnrichedHeatmap::normalizeToMatrix(
          bw_data,
          noStrand_regions,
          extend = c(upstream, downstream),
          w = bin_size,
          mean_mode = "absolute"
        )
      }
      
      matrix_list[[sample_name]] <- matrix
      has_valid_matrices <- TRUE
      
    }, error = function(e) {
      log_message(paste("Error processing", sample_name, ":", e$message), "ERROR")
    })
  }
  
  # Check if we have any valid matrices
  if (!has_valid_matrices) {
    log_message("No valid matrices were generated", "ERROR")
    return(NULL)
  }
  
  # Combine matrices and metadata into a single object
  result <- list(
    matrices = matrix_list,
    regions = regions_gr,
    sample_names = names(matrix_list),
    parameters = list(
      use_scale_regions = use_scale_regions,
      upstream = upstream,
      downstream = downstream,
      bin_size = bin_size,
      region_body_length = region_body_length
    )
  )
  
  # Sort regions if requested
  if (sort_regions) {
    log_message("Sorting regions by signal strength")
    
    # Calculate average signal across all samples for each region
    avg_signal <- matrix(0, nrow = nrow(matrix_list[[1]]), ncol = 1)
    
    for (m in matrix_list) {
      avg_signal <- avg_signal + rowMeans(m, na.rm = TRUE)
    }
    
    avg_signal <- avg_signal / length(matrix_list)
    
    # Sort indices
    sort_idx <- order(avg_signal, decreasing = TRUE)
    
    # Reorder all matrices
    for (i in seq_along(matrix_list)) {
      result$matrices[[i]] <- matrix_list[[i]][sort_idx, ]
    }
    
    # Reorder regions
    result$regions <- regions_gr[sort_idx]
  }
  
  # Save to RDS file
  saveRDS(result, output_file)
  log_message(paste("Matrix data saved to", basename(output_file)))
  
  # Cache the result
  cache_manager$save(cache_key, ".rds", file_path = output_file)
  
  return(output_file)
}

# IMPROVED: Function to plot composite profile directly from matrix data
plot_composite_profile <- function(matrix_file, output_file, title="", 
                                 sample_labels=NULL, colors=NULL, 
                                 ylim=NULL, force_recompute=FORCE_RECOMPUTE) {
  
  # Check if matrix_file is NULL (indicating error in matrix generation)
  if (is.null(matrix_file)) {
    log_message(paste("Cannot create profile plot, matrix file is NULL"), "ERROR")
    return(NULL)
  }
  
  # Check if output file already exists
  if (file.exists(output_file) && !force_recompute) {
    log_message(paste("Output file already exists, skipping:", output_file))
    return(output_file)
  } else if (file.exists(output_file) && force_recompute) {
    log_message(paste("Output file exists but force_recompute=TRUE, regenerating:", output_file))
    file.remove(output_file)
  }
  
  log_message(paste("Creating composite profile plot:", output_file))
  
  # Load matrix data - detect format by file extension
  is_rds <- FALSE
  if (!is.null(matrix_file)) {
    is_rds <- grepl("\\.rds$", matrix_file)
  }
  
  if (is_rds) {
    # For R-native matrices
    matrix_data <- readRDS(matrix_file)
    
    # Get the matrices
    matrices <- matrix_data$matrices
    sample_names <- matrix_data$sample_names
    
    # Use provided sample labels or default to sample names
    if (is.null(sample_labels)) {
      sample_labels <- sample_names
    }
    
    # Use provided colors or default colors
    if (is.null(colors)) {
      if (length(sample_names) <= length(unlist(COLORS))) {
        colors <- unlist(COLORS)[1:length(sample_names)]
      } else {
        # Generate colors from RColorBrewer
        colors <- brewer.pal(min(9, length(sample_names)), "Set1")
        if (length(sample_names) > 9) {
          # Interpolate if we need more than 9 colors
          colors <- colorRampPalette(colors)(length(sample_names))
        }
      }
    }
    
    # Compute average profile for each sample
    profiles <- lapply(matrices, function(m) colMeans(m, na.rm = TRUE))
    
    # Create profile data frame for plotting
    profile_df <- data.frame(
      Position = 1:ncol(matrices[[1]])
    )
    
    for (i in seq_along(profiles)) {
      profile_df[[sample_labels[i]]] <- profiles[[i]]
    }
    
    # Melt the data frame for ggplot
    melted_df <- reshape2::melt(profile_df, id.vars = "Position", 
                              variable.name = "Sample", value.name = "Signal")
    
    # Determine plot type based on matrix mode
    if (matrix_data$parameters$use_scale_regions) {
      # For scale-regions (TSS to TES)
      upstream_bins <- matrix_data$parameters$upstream / matrix_data$parameters$bin_size
      body_bins <- matrix_data$parameters$region_body_length / matrix_data$parameters$bin_size
      downstream_bins <- matrix_data$parameters$downstream / matrix_data$parameters$bin_size
      
      tss_pos <- upstream_bins
      tes_pos <- tss_pos + body_bins
      
      # Create x-axis breaks and labels
      x_breaks <- c(1, tss_pos, tes_pos, ncol(matrices[[1]]))
      x_labels <- c(paste0("-", matrix_data$parameters$upstream, "bp"), 
                   "TSS", "TES",
                   paste0("+", matrix_data$parameters$downstream, "bp"))
      
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
      # For reference-point (TSS-centered)
      total_bins <- (matrix_data$parameters$upstream + matrix_data$parameters$downstream) / 
                    matrix_data$parameters$bin_size
      middle_bin <- total_bins / 2
      
      # Create x-axis breaks and labels
      x_breaks <- c(1, middle_bin/2, middle_bin, middle_bin*1.5, total_bins)
      x_labels <- c(
        paste0("-", matrix_data$parameters$upstream, "bp"),
        paste0("-", matrix_data$parameters$upstream/2, "bp"),
        "TSS",
        paste0("+", matrix_data$parameters$downstream/2, "bp"),
        paste0("+", matrix_data$parameters$downstream, "bp")
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
        geom_vline(xintercept = middle_bin, linetype = "dashed", alpha = 0.5)
    }
    
  } else {
    # Legacy format handling (not needed for new version)
    log_message("Legacy format detected - this should not happen with pure R-native implementation", "ERROR")
    return(NULL)
  }
  
  # Set y-axis limits if provided
  if (!is.null(ylim)) {
    p <- p + ylim(ylim)
  }
  
  # Save the plot
  ggsave(output_file, p, width = 10, height = 6, dpi = 300)
  log_message(paste("Profile plot saved to:", output_file))
  
  # Save the underlying data for reference
  write.csv(melted_df, paste0(output_file, ".csv"), row.names = FALSE)
  
  return(output_file)
}

# Modified plot_heatmap function to handle NULL matrix files
plot_heatmap <- function(matrix_file, output_file, title="", zmin=0, zmax=0.3, 
                        use_scale_regions=FALSE, force_recompute=FORCE_RECOMPUTE, 
                        sample_type=NULL) {
  # Check if matrix_file is NULL (indicating error in matrix generation)
  if (is.null(matrix_file)) {
    log_message(paste("Cannot create heatmap, matrix file is NULL"), "ERROR")
    return(NULL)
  }
  
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
  
  # Load matrix data
  matrix_data <- readRDS(matrix_file)
  
  # Get the matrices
  matrices <- matrix_data$matrices
  sample_names <- matrix_data$sample_names
  
  # Set up color function
  col_fun <- circlize::colorRamp2(c(zmin, (zmin + zmax)/2, zmax), 
                                c("blue", "yellow", "red"))
  
  # Determine heatmap type based on mode
  if (use_scale_regions) {
    # For TSS to TES (scale-regions)
    upstream_bins <- matrix_data$parameters$upstream / matrix_data$parameters$bin_size
    body_bins <- matrix_data$parameters$region_body_length / matrix_data$parameters$bin_size
    downstream_bins <- matrix_data$parameters$downstream / matrix_data$parameters$bin_size
    
    # Create column annotations to mark TSS and TES positions
    column_split <- rep(c("Upstream", "Gene Body", "Downstream"), 
                       c(upstream_bins, body_bins, downstream_bins))
  } else {
    # For TSS-centered (reference-point)
    total_bins <- (matrix_data$parameters$upstream + matrix_data$parameters$downstream) / 
                  matrix_data$parameters$bin_size
    middle_bin <- total_bins / 2
    
    # Create column split for visual reference
    column_split <- ifelse(1:ncol(matrices[[1]]) <= middle_bin, "Upstream", "Downstream")
  }
  
  # Create a list to hold the heatmaps
  ht_list <- lapply(seq_along(matrices), function(i) {
    mat <- matrices[[i]]
    sample <- sample_names[i]
    
    Heatmap(
      mat,
      name = sample,
      col = col_fun,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      column_title = ifelse(i == 1, title, sample),
      column_split = column_split,
      use_raster = TRUE,
      raster_quality = 4,
      heatmap_legend_param = list(
        title = "Signal",
        at = c(zmin, (zmin + zmax)/2, zmax),
        labels = c(zmin, (zmin + zmax)/2, zmax)
      )
    )
  })
  
  # Combine heatmaps
  if (length(ht_list) > 1) {
    combined_ht <- ht_list[[1]] + ht_list[[2]]
    if (length(ht_list) > 2) {
      for (i in 3:length(ht_list)) {
        combined_ht <- combined_ht + ht_list[[i]]
      }
    }
  } else {
    combined_ht <- ht_list[[1]]
  }
  
  # Create a profile above the heatmap
  profile_data <- lapply(matrices, function(m) colMeans(m, na.rm = TRUE))
  profile_df <- data.frame(
    Position = 1:ncol(matrices[[1]])
  )
  for (i in seq_along(profile_data)) {
    profile_df[[sample_names[i]]] <- profile_data[[i]]
  }
  
  # Convert to long format for plotting
  melted_df <- reshape2::melt(profile_df, id.vars = "Position", 
                            variable.name = "Sample", value.name = "Signal")
  
  # Create the profile plot
  if (use_scale_regions) {
    # For scale-regions, mark TSS and TES
    tss_pos <- upstream_bins
    tes_pos <- tss_pos + body_bins
    
    p <- ggplot(melted_df, aes(x = Position, y = Signal, color = Sample)) +
      geom_line(size = 1) +
      scale_color_manual(values = COLORS[1:length(unique(melted_df$Sample))]) +
      labs(title = paste("Profile:", title), x = "", y = "Signal") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      ) +
      geom_vline(xintercept = tss_pos, linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = tes_pos, linetype = "dashed", alpha = 0.5) +
      annotate("text", x = tss_pos - 5, y = max(melted_df$Signal) * 0.9, label = "TSS") +
      annotate("text", x = tes_pos + 5, y = max(melted_df$Signal) * 0.9, label = "TES")
  } else {
    # For reference-point (TSS-centered)
    middle_bin <- (matrix_data$parameters$upstream + matrix_data$parameters$downstream) / 
                  (2 * matrix_data$parameters$bin_size)
    
    p <- ggplot(melted_df, aes(x = Position, y = Signal, color = Sample)) +
      geom_line(size = 1) +
      scale_color_manual(values = COLORS[1:length(unique(melted_df$Sample))]) +
      labs(title = paste("Profile:", title), x = "Distance from TSS", y = "Signal") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      ) +
      geom_vline(xintercept = middle_bin, linetype = "dashed", alpha = 0.5) +
      annotate("text", x = middle_bin - 5, y = max(melted_df$Signal) * 0.9, label = "TSS")
  }
  
  # Save plot and heatmap together
  png(output_file, width = 10, height = 8, units = "in", res = 300)
  # Draw the profile plot at the top
  profile_grob <- ggplotGrob(p)
  # Draw the heatmap below
  pushViewport(viewport(x = 0.5, y = 0.3, width = 0.8, height = 0.6))
  draw(combined_ht)
  upViewport()
  # Draw the profile at the top
  pushViewport(viewport(x = 0.5, y = 0.8, width = 0.8, height = 0.3))
  grid.draw(profile_grob)
  upViewport()
  dev.off()
  
  log_message(paste("Heatmap saved to:", output_file))
  
  # Save profile data for reference
  profile_output <- paste0(sub("\\.png$", "", output_file), "_profile.csv")
  write.csv(melted_df, profile_output, row.names = FALSE)
  
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
    # Extract gene regions - passing FALSE to include all genes
    genes <- genes(txdb, single.strand.genes.only=FALSE)
    log_message(paste("Extracted", length(genes), "gene regions"))
    
    # For GRangesList from single.strand.genes.only=FALSE, we need to unlist to get a GRanges
    if (is(genes, "GRangesList")) {
      log_message("Flattening GRangesList to GRanges for compatibility")
      genes <- unlist(genes)
    }
    
    # Export gene regions
    export(genes, gene_regions_file)
    log_message(paste("Gene regions exported to:", gene_regions_file))
    
    log_message("Extracting TSS regions...")
    # Extract TSS regions - ensure we handle both GRanges and GRangesList
    tss_regions <- promoters(genes, upstream=1, downstream=1)
    log_message(paste("Created", length(tss_regions), "TSS regions"))
    
    # Export TSS regions
    export(tss_regions, tss_file)
    log_message(paste("TSS regions exported to:", tss_file))
  }
  
  return(list(gene_regions=gene_regions_file, tss=tss_file))
}

# Function to generate heatmaps and metaprofiles for Mecp2 Endo and Exo around TSS
generate_mecp2_tss_heatmaps <- function(temp_dir, annotation_files) {
  log_message("Starting R-native analysis of MeCP2 binding patterns")
  tracker <- ProgressTracker$new(18)
  
  # Define sample groups
  sample_groups <- list(
    neuron_endo = c("NeuM2", "NeuM3"),
    neuron_exo = c("NeuV1", "NeuV2", "NeuV3"),
    nsc_endo = c("NSCM1", "NSCM2", "NSCM3"),
    nsc_exo = c("NSCv1", "NSCv2", "NSCv3")
  )
  
  # Sample labels for plots
  sample_labels <- list(
    neuron_endo = "Neuron Endogenous",
    neuron_exo = "Neuron Exogenous",
    nsc_endo = "NSC Endogenous",
    nsc_exo = "NSC Exogenous"
  )
  
  # Pre-compute all matrices using our R-native approach
  matrices <- list()
  for (group_name in names(sample_groups)) {
    bw_files <- file.path(BIGWIG_DIR, paste0(sample_groups[[group_name]], ".bw"))
    
    # TSS matrices
    matrix_file <- file.path(temp_dir, paste0(group_name, "_tss_matrix.rds"))
    matrices[[paste0(group_name, "_tss")]] <- compute_matrix(
      bw_files,
      annotation_files$tss,
      matrix_file,
      bin_size = 50,
      sort_regions = TRUE
    )
    tracker$increment()
    
    # TSS to TES matrices
    matrix_file <- file.path(temp_dir, paste0(group_name, "_tss_tes_matrix.rds"))
    matrices[[paste0(group_name, "_tss_tes")]] <- compute_matrix(
      bw_files,
      annotation_files$gene_regions,
      matrix_file,
      use_scale_regions = TRUE,
      bin_size = 50,
      region_body_length = 5000,
      sort_regions = TRUE
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
    
    # TSS profile
    plot_composite_profile(
      matrices[[paste0(group_name, "_tss")]],
      file.path(OUTPUT_DIR, paste0(group_name, "_tss_profile.png")),
      title = paste(sample_labels[[group_name]], "MeCP2 around TSS"),
      sample_labels = c(sample_labels[[group_name]]),
      colors = c(COLORS[[sample_type]])
    )
    tracker$increment()
    
    # Gene body profile
    plot_composite_profile(
      matrices[[paste0(group_name, "_tss_tes")]],
      file.path(OUTPUT_DIR, paste0(group_name, "_tss_tes_profile.png")),
      title = paste(sample_labels[[group_name]], "MeCP2 gene body coverage"),
      sample_labels = c(sample_labels[[group_name]]),
      colors = c(COLORS[[sample_type]])
    )
    tracker$increment()
  }
  
  # Function to generate combined plots for Endo vs Exo in a cell type
  generate_combined_plots <- function(cell_type) {
    log_message(paste("Generating combined plots for", cell_type))
    
    endo_group <- paste0(cell_type, "_endo")
    exo_group <- paste0(cell_type, "_exo")
    endo_type <- paste0(toupper(substring(cell_type, 1, 3)), "_ENDO")
    exo_type <- paste0(toupper(substring(cell_type, 1, 3)), "_EXO")
    
    # Create combined matrices from individual files
    # For TSS (reference-point)
    combined_bw_files <- c(
      file.path(BIGWIG_DIR, paste0(sample_groups[[endo_group]], ".bw")),
      file.path(BIGWIG_DIR, paste0(sample_groups[[exo_group]], ".bw"))
    )
    
    combined_tss_matrix_file <- file.path(temp_dir, paste0(cell_type, "_combined_tss_matrix.rds"))
    combined_tss_matrix <- compute_matrix(
      combined_bw_files,
      annotation_files$tss,
      combined_tss_matrix_file,
      bin_size = 50,
      sort_regions = TRUE
    )
    tracker$increment()
    
    # For gene body (scale-regions)
    combined_tss_tes_matrix_file <- file.path(temp_dir, paste0(cell_type, "_combined_tss_tes_matrix.rds"))
    combined_tss_tes_matrix <- compute_matrix(
      combined_bw_files,
      annotation_files$gene_regions,
      combined_tss_tes_matrix_file,
      use_scale_regions = TRUE,
      bin_size = 50,
      region_body_length = 5000,
      sort_regions = TRUE
    )
    tracker$increment()
    
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
  
  # Generate diagnostic plots to compare TSS and TES signal heights
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
      # Get the matrix data
      matrix_file <- matrices[[paste0(group_name, "_tss_tes")]]
      
      # Properly handle both formats (R-native RDS or deepTools matrix)
      is_rds <- grepl("\\.rds$", matrix_file)
      
      if (is_rds) {
        # For R-native matrices
        matrix_data <- readRDS(matrix_file)
        
        if (!is.null(matrix_data) && matrix_data$parameters$use_scale_regions) {
          # Get the first matrix for this group (assuming similar patterns across replicates)
          matrix <- matrix_data$matrices[[1]]
          profile <- colMeans(matrix, na.rm = TRUE)
          
          # Calculate positions for TSS and TES
          upstream_bins <- matrix_data$parameters$upstream / matrix_data$parameters$bin_size
          body_bins <- matrix_data$parameters$region_body_length / matrix_data$parameters$bin_size
          
          tss_idx <- upstream_bins
          tes_idx <- tss_idx + body_bins
          
          # Extract signal values at TSS and TES
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
      } else {
        # For legacy deepTools matrices (which we aren't using anymore, but keeping for compatibility)
        # Use the old parsing method...
        log_message(paste("Legacy format detected for", basename(matrix_file)), "WARNING")
        next
      }
    }
    
    # Proceed only if we have diagnostic data
    if (nrow(diagnostics) > 0) {
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
    } else {
      log_message("No diagnostic data was generated", "WARNING")
    }
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
    log_message("=== Starting MeCP2 analysis pipeline (R-native approach) ===")
    log_message(paste("System time:", Sys.time()))
    log_message(paste("R version:", R.version.string))
    log_message(paste("Platform:", R.version$platform))
    
    log_message("Loading required packages")
    print(sessionInfo())
    
    # Prepare files once and pass to downstream functions
    annotation_files <- prepare_annotation_files(temp_dir)
    
    # Run R-native analysis
    log_message("Generating heatmaps and profiles using R-native methodology")
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