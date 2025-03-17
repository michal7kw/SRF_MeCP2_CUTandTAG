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
# NSCv3,Exogenous

setwd("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/paper_finals")
OUTPUT_DIR <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/paper_finals/outputs"
BIGWIG_DIR <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results_1b/bigwig"
PEAKS_DIR <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/peaks/broad"
TSS_BED <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA/mm10_TSS.bed"
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
compute_matrix <- function(bw_files, regions, output_file, upstream=5000, downstream=5000, use_scale_regions=FALSE, force_recompute=FORCE_RECOMPUTE) {
  # Create cache key based on input parameters
  cache_key <- list(
    bw_files = bw_files,
    regions = regions,
    upstream = upstream,
    downstream = downstream,
    use_scale_regions = use_scale_regions,
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
  
  log_message(paste0("Computing matrix for ", ifelse(use_scale_regions, "TSS-to-TES", "TSS-centered"), " analysis..."))
  log_message(paste0("Using bigwig files: ", paste(basename(bw_files), collapse=", ")))
  log_message(paste0("Using regions file: ", basename(regions)))
  log_message(paste0("Output will be saved to: ", basename(output_file)))
  
  bw_files_str <- paste(bw_files, collapse = " ")
  
  if (use_scale_regions) {
    # For TSS to TES scaled regions
    cmd <- paste0("computeMatrix scale-regions",
                 " -S ", bw_files_str,
                 " -R ", regions,
                 " -b ", upstream, " -a ", downstream,
                 " --regionBodyLength 5000",
                 " -o ", output_file,
                 " --skipZeros",
                 " --missingDataAsZero")
  } else {
    # Original reference-point approach
    cmd <- paste0("computeMatrix reference-point",
                 " -S ", bw_files_str,
                 " -R ", regions,
                 " --referencePoint TSS",
                 " -b ", upstream, " -a ", downstream,
                 " -o ", output_file,
                 " --skipZeros",
                 " --missingDataAsZero")
  }
  log_message("Executing command:")
  log_message(cmd)
  system_result <- system(cmd)
  
  if (system_result == 0) {
    log_message("Matrix computation completed successfully.")
  } else {
    log_message(paste("Matrix computation failed with exit code:", system_result), "ERROR")
  }
  
  # Cache the result
  if (system_result == 0) {
    cache_manager$save(cache_key, ".matrix.gz", file_path = output_file)
  }
  
  return(output_file)
}

# Function to plot heatmap using deepTools
plot_heatmap <- function(matrix_file, output_file, title="", zmin=0, zmax=0.3, use_scale_regions=FALSE, force_recompute=FORCE_RECOMPUTE, sample_type=NULL) {
  # Check if output file already exists
  if (file.exists(output_file) && !force_recompute) {
    log_message(paste("Output file already exists, skipping:", output_file))
    return(output_file)
  } else if (file.exists(output_file) && force_recompute) {
    log_message(paste("Output file exists but force_recompute=TRUE, regenerating:", output_file))
    # Remove the existing file to ensure clean regeneration
    file.remove(output_file)
  }
  
  log_message(paste0("Plotting heatmap for ", ifelse(use_scale_regions, "TSS-to-TES", "TSS-centered"), " analysis..."))
  log_message(paste0("Using matrix file: ", basename(matrix_file)))
  log_message(paste0("Output will be saved to: ", basename(output_file)))
  log_message(paste0("Title: ", title))
  
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
                 " --endLabel \"TES\"")
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

# Function to plot profile using deepTools - fixed parameters
plot_profile <- function(matrix_file, output_file, title="", force_recompute=FORCE_RECOMPUTE, sample_types=NULL) {
  # Check if output file already exists
  if (file.exists(output_file) && !force_recompute) {
    log_message(paste("Output file already exists, skipping:", output_file))
    return(output_file)
  } else if (file.exists(output_file) && force_recompute) {
    log_message(paste("Output file exists but force_recompute=TRUE, regenerating:", output_file))
    # Remove the existing file to ensure clean regeneration
    file.remove(output_file)
  }
  
  log_message(paste0("Plotting profile..."))
  log_message(paste0("Using matrix file: ", basename(matrix_file)))
  log_message(paste0("Output will be saved to: ", basename(output_file)))
  log_message(paste0("Title: ", title))
  
  # Get colors based on sample types
  colors <- NULL
  if (!is.null(sample_types)) {
    colors <- sapply(sample_types, function(type) COLORS[[type]])
    colors <- paste(colors, collapse=" ")
  } else {
    # Default colors if no sample types specified
    colors <- paste(COLORS$NEU_ENDO, COLORS$NEU_EXO, COLORS$NSC_ENDO, COLORS$NSC_EXO, collapse=" ")
  }
  
  cmd <- paste0("plotProfile -m ", matrix_file,
               " -o ", output_file,
               " --colors ", colors,
               " --plotTitle \"", title, "\"")
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

# Function to generate heatmaps and metaprofiles for Mecp2 Endo and Exo around TSS
generate_mecp2_tss_heatmaps <- function(temp_dir, annotation_files) {
  tracker <- ProgressTracker$new(18)
  
  # Create a function to get or compute matrix
  get_or_compute_matrix <- function(bw_files, regions, output_file, use_scale_regions = FALSE, ...) {
    if (file.exists(output_file) && !FORCE_RECOMPUTE) {
      log_message(paste("Reusing existing matrix:", basename(output_file)))
      return(output_file)
    }
    compute_matrix(bw_files, regions, output_file, use_scale_regions = use_scale_regions, ...)
  }
  
  # Define sample groups
  sample_groups <- list(
    neuron_endo = c("NeuM2", "NeuM3"),
    neuron_exo = c("NeuV1", "NeuV2", "NeuV3"),
    nsc_endo = c("NSCM1", "NSCM2", "NSCM3"),
    nsc_exo = c("NSCv1", "NSCv2", "NSCv3")
  )
  
  # Pre-compute all matrices in parallel if possible
  matrices <- list()
  for (group_name in names(sample_groups)) {
    bw_files <- file.path(BIGWIG_DIR, paste0(sample_groups[[group_name]], ".bw"))
    
    # TSS matrices
    matrices[[paste0(group_name, "_tss")]] <- get_or_compute_matrix(
      bw_files,
      annotation_files$tss,
      file.path(temp_dir, paste0(group_name, "_tss_matrix.gz"))
    )
    
    # TSS to TES matrices
    matrices[[paste0(group_name, "_tss_tes")]] <- get_or_compute_matrix(
      bw_files,
      annotation_files$gene_regions,
      file.path(temp_dir, paste0(group_name, "_tss_tes_matrix.gz")),
      use_scale_regions = TRUE
    )
  }
  
  # Function to generate plots for a group
  generate_group_plots <- function(group_name, sample_type) {
    # TSS to TES heatmap
    plot_heatmap(
      matrices[[paste0(group_name, "_tss_tes")]],
      file.path(OUTPUT_DIR, paste0(group_name, "_tss_tes_heatmap.png")),
      title = paste(group_name, "MeCP2 gene body coverage"),
      use_scale_regions = TRUE,
      sample_type = sample_type
    )
    
    # TSS heatmap
    plot_heatmap(
      matrices[[paste0(group_name, "_tss")]],
      file.path(OUTPUT_DIR, paste0(group_name, "_tss_heatmap.png")),
      title = paste(group_name, "MeCP2 around TSS"),
      sample_type = sample_type
    )
    
    # TSS profile
    plot_profile(
      matrices[[paste0(group_name, "_tss")]],
      file.path(OUTPUT_DIR, paste0(group_name, "_tss_profile.png")),
      title = paste(group_name, "MeCP2 around TSS"),
      sample_types = c(sample_type)
    )

    # TSS to TES profile
    plot_profile(
      matrices[[paste0(group_name, "_tss_tes")]],
      file.path(OUTPUT_DIR, paste0(group_name, "_tss_tes_profile.png")),
      title = paste(group_name, "MeCP2 gene body coverage"),
      sample_types = c(sample_type)
    )
  }
  
  # Generate all individual plots
  generate_group_plots("neuron_endo", "NEU_ENDO")
  generate_group_plots("neuron_exo", "NEU_EXO")
  generate_group_plots("nsc_endo", "NSC_ENDO")
  generate_group_plots("nsc_exo", "NSC_EXO")
  
  # Generate combined plots
  generate_combined_plots <- function(cell_type) {
    endo_group <- paste0(cell_type, "_endo")
    exo_group <- paste0(cell_type, "_exo")
    endo_type <- paste0(toupper(cell_type), "_ENDO")
    exo_type <- paste0(toupper(cell_type), "_EXO")
    
    # Combined TSS matrix
    combined_matrix <- get_or_compute_matrix(
      c(file.path(BIGWIG_DIR, paste0(sample_groups[[endo_group]][1], ".bw")),
        file.path(BIGWIG_DIR, paste0(sample_groups[[exo_group]][1], ".bw"))),
      annotation_files$tss,
      file.path(temp_dir, paste0(cell_type, "_combined_tss_matrix.gz"))
    )

    # Combined TSS to TES matrix
     combined_tss_tes_matrix <- get_or_compute_matrix(
      c(file.path(BIGWIG_DIR, paste0(sample_groups[[endo_group]][1], ".bw")),
        file.path(BIGWIG_DIR, paste0(sample_groups[[exo_group]][1], ".bw"))),
      annotation_files$gene_regions,
      file.path(temp_dir, paste0(cell_type, "_combined_tss_tes_matrix.gz")),
      use_scale_regions = TRUE
    )
    
    # Combined plots
    plot_heatmap(
      combined_matrix,
      file.path(OUTPUT_DIR, paste0(cell_type, "_endo_vs_exo_tss_heatmap.png")),
      title = paste(cell_type, "Endogenous vs Exogenous MeCP2 around TSS")
    )
    
    plot_profile(
      combined_matrix,
      file.path(OUTPUT_DIR, paste0(cell_type, "_endo_vs_exo_tss_profile.png")),
      title = paste(cell_type, "Endogenous vs Exogenous MeCP2 around TSS"),
      sample_types = c(endo_type, exo_type)
    )

    # Combined TSS to TES plots
    plot_profile(
      combined_tss_tes_matrix,
      file.path(OUTPUT_DIR, paste0(cell_type, "_endo_vs_exo_tss_tes_profile.png")),
      title = paste(cell_type, "Endogenous vs Exogenous MeCP2 gene body coverage"),
      sample_types = c(endo_type, exo_type)
    )
  }
  
  # Generate combined plots for both cell types
  generate_combined_plots("neuron")
  generate_combined_plots("nsc")
  
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
    
    log_message("Generating heatmaps")
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