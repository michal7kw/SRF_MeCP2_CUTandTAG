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
GENCODE_GTF <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/DATA/gencode.vM25.annotation.gtf"

# Set to TRUE to use GENCODE GTF instead of TxDb
USE_GENCODE <- FALSE

# Create output directory if it doesn't exist
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Set this to TRUE to force recomputation of all matrices and heatmaps
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
compute_matrix <- function(bw_files, regions, output_file, upstream=5000, downstream=5000, use_scale_regions=FALSE, force_recompute=FORCE_RECOMPUTE) {
  # Check if output file already exists
  if (file.exists(output_file) && !force_recompute) {
    log_message(paste("Matrix file already exists, skipping:", output_file))
    return(output_file)
  } else if (file.exists(output_file) && force_recompute) {
    log_message(paste("Matrix file exists but force_recompute=TRUE, regenerating:", output_file))
    # Remove the existing file to ensure clean regeneration
    file.remove(output_file)
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
  
  return(output_file)
}

# Function to plot heatmap using deepTools
plot_heatmap <- function(matrix_file, output_file, title="", zmin=0, zmax=0.3, use_scale_regions=FALSE, force_recompute=FORCE_RECOMPUTE) {
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

# Function to plot profile using deepTools - fixed parameters
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
  
  log_message(paste0("Plotting profile..."))
  log_message(paste0("Using matrix file: ", basename(matrix_file)))
  log_message(paste0("Output will be saved to: ", basename(output_file)))
  log_message(paste0("Title: ", title))
  
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
  tracker <- ProgressTracker$new(18)  # Total number of main steps
  
  log_message("=== Starting MeCP2 heatmap generation ===")
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
  
  # Generate matrices and heatmaps for Neuron Endogenous vs Exogenous
  log_message("=== Processing Neuron Endogenous vs Exogenous ===")
  neuron_endo_bw <- file.path(BIGWIG_DIR, paste0(neuron_endo, ".bw"))
  neuron_exo_bw <- file.path(BIGWIG_DIR, paste0(neuron_exo, ".bw"))
  
  log_message(paste("Neuron Endogenous bigwig files:", paste(basename(neuron_endo_bw), collapse=", ")))
  log_message(paste("Neuron Exogenous bigwig files:", paste(basename(neuron_exo_bw), collapse=", ")))
  
  # Combine bigwig files for each group
  neuron_endo_matrix <- file.path(temp_dir, "neuron_endo_tss_matrix.gz")
  neuron_exo_matrix <- file.path(temp_dir, "neuron_exo_tss_matrix.gz")
  
  # Compute matrices for TSS-centered analysis
  log_message("=== Computing TSS-centered matrices for Neurons ===")
  compute_matrix(neuron_endo_bw, annotation_files$tss, neuron_endo_matrix)
  tracker$increment()
  compute_matrix(neuron_exo_bw, annotation_files$tss, neuron_exo_matrix)
  tracker$increment()
  
  # Create TSS to TES matrices
  log_message("=== Computing TSS-to-TES matrices for Neurons ===")
  neuron_endo_tss_tes_matrix <- file.path(temp_dir, "neuron_endo_tss_tes_matrix.gz")
  neuron_exo_tss_tes_matrix <- file.path(temp_dir, "neuron_exo_tss_tes_matrix.gz")
  
  # Compute matrices for TSS to TES analysis
  compute_matrix(neuron_endo_bw, annotation_files$gene_regions, neuron_endo_tss_tes_matrix, upstream=5000, downstream=5000, use_scale_regions=TRUE)
  tracker$increment()
  compute_matrix(neuron_exo_bw, annotation_files$gene_regions, neuron_exo_tss_tes_matrix, upstream=5000, downstream=5000, use_scale_regions=TRUE)
  tracker$increment()
  
  # Plot TSS to TES heatmaps
  log_message("=== Plotting TSS-to-TES heatmaps for Neurons ===")
  neuron_endo_tss_tes_heatmap <- file.path(OUTPUT_DIR, "neuron_endo_tss_tes_heatmap.png")
  plot_heatmap(neuron_endo_tss_tes_matrix, neuron_endo_tss_tes_heatmap, "Endogenous", use_scale_regions=TRUE)
  tracker$increment()
  
  neuron_exo_tss_tes_heatmap <- file.path(OUTPUT_DIR, "neuron_exo_tss_tes_heatmap.png")
  plot_heatmap(neuron_exo_tss_tes_matrix, neuron_exo_tss_tes_heatmap, "Exogenous", use_scale_regions=TRUE)
  tracker$increment()
  
  # Plot heatmaps for TSS-centered analysis
  log_message("=== Plotting TSS-centered heatmaps for Neurons ===")
  neuron_endo_heatmap <- file.path(OUTPUT_DIR, "neuron_endo_tss_heatmap.png")
  plot_heatmap(neuron_endo_matrix, neuron_endo_heatmap, "Neuron Endogenous MeCP2 around TSS")
  tracker$increment()
  
  neuron_exo_heatmap <- file.path(OUTPUT_DIR, "neuron_exo_tss_heatmap.png")
  plot_heatmap(neuron_exo_matrix, neuron_exo_heatmap, "Neuron Exogenous MeCP2 around TSS")
  tracker$increment()
  
  # Plot profiles
  log_message("=== Plotting TSS-centered profiles for Neurons ===")
  neuron_endo_profile <- file.path(OUTPUT_DIR, "neuron_endo_tss_profile.png")
  plot_profile(neuron_endo_matrix, neuron_endo_profile, "Neuron Endogenous MeCP2 around TSS")
  tracker$increment()
  
  neuron_exo_profile <- file.path(OUTPUT_DIR, "neuron_exo_tss_profile.png")
  plot_profile(neuron_exo_matrix, neuron_exo_profile, "Neuron Exogenous MeCP2 around TSS")
  tracker$increment()
  
  # Generate matrices and heatmaps for NSC Endogenous vs Exogenous
  log_message("=== Processing NSC Endogenous vs Exogenous ===")
  nsc_endo_bw <- file.path(BIGWIG_DIR, paste0(nsc_endo, ".bw"))
  nsc_exo_bw <- file.path(BIGWIG_DIR, paste0(nsc_exo, ".bw"))
  
  log_message(paste("NSC Endogenous bigwig files:", paste(basename(nsc_endo_bw), collapse=", ")))
  log_message(paste("NSC Exogenous bigwig files:", paste(basename(nsc_exo_bw), collapse=", ")))
  
  # Combine bigwig files for each group
  nsc_endo_matrix <- file.path(temp_dir, "nsc_endo_tss_matrix.gz")
  nsc_exo_matrix <- file.path(temp_dir, "nsc_exo_tss_matrix.gz")
  
  # Compute matrices for TSS-centered analysis
  log_message("=== Computing TSS-centered matrices for NSCs ===")
  compute_matrix(nsc_endo_bw, annotation_files$tss, nsc_endo_matrix)
  tracker$increment()
  compute_matrix(nsc_exo_bw, annotation_files$tss, nsc_exo_matrix)
  tracker$increment()
  
  # Create TSS to TES matrices for NSC
  log_message("=== Computing TSS-to-TES matrices for NSCs ===")
  nsc_endo_tss_tes_matrix <- file.path(temp_dir, "nsc_endo_tss_tes_matrix.gz")
  nsc_exo_tss_tes_matrix <- file.path(temp_dir, "nsc_exo_tss_tes_matrix.gz")
  
  # Compute matrices for TSS to TES analysis
  compute_matrix(nsc_endo_bw, annotation_files$gene_regions, nsc_endo_tss_tes_matrix, upstream=5000, downstream=5000, use_scale_regions=TRUE)
  tracker$increment()
  compute_matrix(nsc_exo_bw, annotation_files$gene_regions, nsc_exo_tss_tes_matrix, upstream=5000, downstream=5000, use_scale_regions=TRUE)
  tracker$increment()
  
  # Plot TSS to TES heatmaps for NSC
  log_message("=== Plotting TSS-to-TES heatmaps for NSCs ===")
  nsc_endo_tss_tes_heatmap <- file.path(OUTPUT_DIR, "nsc_endo_tss_tes_heatmap.png")
  plot_heatmap(nsc_endo_tss_tes_matrix, nsc_endo_tss_tes_heatmap, "Endogenous", use_scale_regions=TRUE)
  tracker$increment()
  
  nsc_exo_tss_tes_heatmap <- file.path(OUTPUT_DIR, "nsc_exo_tss_tes_heatmap.png")
  plot_heatmap(nsc_exo_tss_tes_matrix, nsc_exo_tss_tes_heatmap, "Exogenous", use_scale_regions=TRUE)
  tracker$increment()
  
  # Plot heatmaps
  log_message("=== Plotting TSS-centered heatmaps for NSCs ===")
  nsc_endo_heatmap <- file.path(OUTPUT_DIR, "nsc_endo_tss_heatmap.png")
  plot_heatmap(nsc_endo_matrix, nsc_endo_heatmap, "NSC Endogenous MeCP2 around TSS")
  tracker$increment()
  
  nsc_exo_heatmap <- file.path(OUTPUT_DIR, "nsc_exo_tss_heatmap.png")
  plot_heatmap(nsc_exo_matrix, nsc_exo_heatmap, "NSC Exogenous MeCP2 around TSS")
  tracker$increment()
  
  # Plot profiles
  log_message("=== Plotting TSS-centered profiles for NSCs ===")
  nsc_endo_profile <- file.path(OUTPUT_DIR, "nsc_endo_tss_profile.png")
  plot_profile(nsc_endo_matrix, nsc_endo_profile, "NSC Endogenous MeCP2 around TSS")
  tracker$increment()
  
  nsc_exo_profile <- file.path(OUTPUT_DIR, "nsc_exo_tss_profile.png")
  plot_profile(nsc_exo_matrix, nsc_exo_profile, "NSC Exogenous MeCP2 around TSS")
  tracker$increment()
  
  # Create combined plots for comparison
  # Combine Endo vs Exo for Neurons
  log_message("=== Creating combined comparison plots ===")
  neuron_combined_matrix <- file.path(temp_dir, "neuron_combined_tss_matrix.gz")
  compute_matrix(c(neuron_endo_bw[1], neuron_exo_bw[1]), annotation_files$tss, neuron_combined_matrix)
  tracker$increment()
  
  neuron_combined_heatmap <- file.path(OUTPUT_DIR, "neuron_endo_vs_exo_tss_heatmap.png")
  plot_heatmap(neuron_combined_matrix, neuron_combined_heatmap, "Neuron Endogenous vs Exogenous MeCP2 around TSS")
  tracker$increment()
  
  neuron_combined_profile <- file.path(OUTPUT_DIR, "neuron_endo_vs_exo_tss_profile.png")
  plot_profile(neuron_combined_matrix, neuron_combined_profile, "Neuron Endogenous vs Exogenous MeCP2 around TSS")
  tracker$increment()
  
  # Combine Endo vs Exo for NSCs
  nsc_combined_matrix <- file.path(temp_dir, "nsc_combined_tss_matrix.gz")
  compute_matrix(c(nsc_endo_bw[1], nsc_exo_bw[1]), annotation_files$tss, nsc_combined_matrix)
  tracker$increment()
  
  nsc_combined_heatmap <- file.path(OUTPUT_DIR, "nsc_endo_vs_exo_tss_heatmap.png")
  plot_heatmap(nsc_combined_matrix, nsc_combined_heatmap, "NSC Endogenous vs Exogenous MeCP2 around TSS")
  tracker$increment()
  
  nsc_combined_profile <- file.path(OUTPUT_DIR, "nsc_endo_vs_exo_tss_profile.png")
  plot_profile(nsc_combined_matrix, nsc_combined_profile, "NSC Endogenous vs Exogenous MeCP2 around TSS")
  tracker$increment()
  
  # Create combined TSS to TES plots
  neuron_combined_tss_tes_matrix <- file.path(temp_dir, "neuron_combined_tss_tes_matrix.gz")
  compute_matrix(c(neuron_endo_bw[1], neuron_exo_bw[1]), annotation_files$gene_regions, neuron_combined_tss_tes_matrix, upstream=5000, downstream=5000, use_scale_regions=TRUE)
  tracker$increment()
  
  neuron_combined_tss_tes_heatmap <- file.path(OUTPUT_DIR, "neuron_endo_vs_exo_tss_tes_heatmap.png")
  plot_heatmap(neuron_combined_tss_tes_matrix, neuron_combined_tss_tes_heatmap, "Neuron Endogenous vs Exogenous", use_scale_regions=TRUE)
  tracker$increment()
  
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