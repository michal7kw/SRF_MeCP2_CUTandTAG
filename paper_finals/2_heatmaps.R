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
    clean_tss_file <- file.path(temp_dir, "clean_tss.bed")
    export(tss_clean, clean_tss_file)
    
    return(clean_tss_file)
  }
  
  # Prepare clean TSS file
  clean_tss_file <- prepare_tss_file()
  
  # Function to compute matrix using deepTools
  compute_matrix <- function(bw_files, regions, output_file, upstream=5000, downstream=5000) {
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
    cmd <- paste0("plotProfile -m ", matrix_file,
                 " -o ", output_file,
                 " --colors blue red",
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
  plot_heatmap(neuron_endo_matrix, 
              file.path(OUTPUT_DIR, "neuron_endo_tss_heatmap.png"), 
              "Neuron Endogenous MeCP2 around TSS")
  
  plot_heatmap(neuron_exo_matrix, 
              file.path(OUTPUT_DIR, "neuron_exo_tss_heatmap.png"), 
              "Neuron Exogenous MeCP2 around TSS")
  
  # Plot profiles
  plot_profile(neuron_endo_matrix, 
              file.path(OUTPUT_DIR, "neuron_endo_tss_profile.png"), 
              "Neuron Endogenous MeCP2 around TSS")
  
  plot_profile(neuron_exo_matrix, 
              file.path(OUTPUT_DIR, "neuron_exo_tss_profile.png"), 
              "Neuron Exogenous MeCP2 around TSS")
  
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
  plot_heatmap(nsc_endo_matrix, 
              file.path(OUTPUT_DIR, "nsc_endo_tss_heatmap.png"), 
              "NSC Endogenous MeCP2 around TSS")
  
  plot_heatmap(nsc_exo_matrix, 
              file.path(OUTPUT_DIR, "nsc_exo_tss_heatmap.png"), 
              "NSC Exogenous MeCP2 around TSS")
  
  # Plot profiles
  plot_profile(nsc_endo_matrix, 
              file.path(OUTPUT_DIR, "nsc_endo_tss_profile.png"), 
              "NSC Endogenous MeCP2 around TSS")
  
  plot_profile(nsc_exo_matrix, 
              file.path(OUTPUT_DIR, "nsc_exo_tss_profile.png"), 
              "NSC Exogenous MeCP2 around TSS")
  
  # Create combined plots for comparison
  # Combine Endo vs Exo for Neurons
  neuron_combined_matrix <- file.path(temp_dir, "neuron_combined_tss_matrix.gz")
  compute_matrix(c(neuron_endo_bw[1], neuron_exo_bw[1]), clean_tss_file, neuron_combined_matrix)
  
  plot_heatmap(neuron_combined_matrix, 
              file.path(OUTPUT_DIR, "neuron_endo_vs_exo_tss_heatmap.png"), 
              "Neuron Endogenous vs Exogenous MeCP2 around TSS")
  
  plot_profile(neuron_combined_matrix, 
              file.path(OUTPUT_DIR, "neuron_endo_vs_exo_tss_profile.png"), 
              "Neuron Endogenous vs Exogenous MeCP2 around TSS")
  
  # Combine Endo vs Exo for NSCs
  nsc_combined_matrix <- file.path(temp_dir, "nsc_combined_tss_matrix.gz")
  compute_matrix(c(nsc_endo_bw[1], nsc_exo_bw[1]), clean_tss_file, nsc_combined_matrix)
  
  plot_heatmap(nsc_combined_matrix, 
              file.path(OUTPUT_DIR, "nsc_endo_vs_exo_tss_heatmap.png"), 
              "NSC Endogenous vs Exogenous MeCP2 around TSS")
  
  plot_profile(nsc_combined_matrix, 
              file.path(OUTPUT_DIR, "nsc_endo_vs_exo_tss_profile.png"), 
              "NSC Endogenous vs Exogenous MeCP2 around TSS")
  
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