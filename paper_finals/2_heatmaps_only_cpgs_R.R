# Pure R implementation for MeCP2 CpG island profiling

# Add command line argument parsing
args <- commandArgs(trailingOnly = TRUE)

# Parse command line arguments
arg_list <- list()
for (i in seq_along(args)) {
  if (startsWith(args[i], "--")) {
    arg_name <- sub("^--", "", strsplit(args[i], "=")[[1]][1])
    arg_value <- strsplit(args[i], "=")[[1]][2]
    arg_list[[arg_name]] <- arg_value
  }
}

# Check required arguments
required_args <- c("bw-dir", "cpg-file", "output-dir")
missing_args <- required_args[!required_args %in% names(arg_list)]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

# Set number of cores if provided
if ("cores" %in% names(arg_list)) {
  cores <- as.numeric(arg_list[["cores"]])
} else {
  cores <- 8
}

library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(data.table)
library(parallel)
library(patchwork)
library(BiocParallel)

# Function to import CpG islands BED file
import_cpg_islands <- function(cpg_file) {
  # Read the BED file with custom column names to handle the format
  cpg_data <- read.table(cpg_file, col.names = c("chrom", "start", "end", "id", "type", "count"), 
                         colClasses = c("character", "integer", "integer", "integer", "character", "integer"))
  
  # Create a GRanges object
  cpg_islands <- GRanges(
    seqnames = cpg_data$chrom,
    ranges = IRanges(start = cpg_data$start, end = cpg_data$end),
    id = cpg_data$id,
    count = cpg_data$count
  )
  
  return(cpg_islands)
}

# Function to import bigWig files and extract signal
import_bigwig <- function(bw_file, regions, upstream = 3000, downstream = 3000, bin_size = 50) {
  # Import the bigWig file
  bw <- import(bw_file, as = "RleList")
  
  # Process in parallel for larger datasets
  profile_matrix <- bplapply(seq_along(regions), function(i) {
    region <- regions[i]
    chrom <- as.character(seqnames(region))
    center <- start(region) + (end(region) - start(region)) / 2
    
    # Define the region to extract
    start_pos <- center - upstream
    end_pos <- center + downstream
    
    # Skip if chromosome not in bigWig file
    if(!chrom %in% names(bw)) {
      return(rep(0, (upstream + downstream) / bin_size))
    }
    
    # Extract signal
    signal <- tryCatch({
      # Get signal values
      signal_range <- bw[[chrom]][start_pos:end_pos]
      
      # Bin the signal
      binned_signal <- sapply(seq(1, length(signal_range), bin_size), function(i) {
        end_idx <- min(i + bin_size - 1, length(signal_range))
        mean(signal_range[i:end_idx], na.rm = TRUE)
      })
      
      # Replace NAs with zeros
      binned_signal[is.na(binned_signal)] <- 0
      binned_signal
    }, error = function(e) {
      # Return zeros if an error occurs
      rep(0, (upstream + downstream) / bin_size)
    })
    
    return(signal)
  }, BPPARAM = MulticoreParam(min(8, detectCores())))
  
  # Convert to matrix
  profile_matrix <- do.call(rbind, profile_matrix)
  return(profile_matrix)
}

# Function to merge multiple bigWig profiles
merge_profiles <- function(bw_files, regions, upstream = 3000, downstream = 3000, bin_size = 50) {
  # Process each bigWig file
  all_profiles <- lapply(bw_files, function(bw_file) {
    import_bigwig(bw_file, regions, upstream, downstream, bin_size)
  })
  
  # Calculate mean profile
  n_profiles <- length(all_profiles)
  if(n_profiles == 0) return(NULL)
  
  # Initialize with first profile
  merged_profile <- all_profiles[[1]]
  
  # Add remaining profiles
  if(n_profiles > 1) {
    for(i in 2:n_profiles) {
      merged_profile <- merged_profile + all_profiles[[i]]
    }
  }
  
  # Calculate mean
  merged_profile <- merged_profile / n_profiles
  return(merged_profile)
}

# Function to generate profile plot
plot_profile <- function(profile_matrix, x_values = NULL, title = "MeCP2 Profile", sample_label = "Sample") {
  # Create x-axis values if not provided
  if(is.null(x_values)) {
    n_bins <- ncol(profile_matrix)
    x_values <- seq(-3000, 3000, length.out = n_bins)
  }
  
  # Calculate mean profile
  mean_profile <- colMeans(profile_matrix, na.rm = TRUE)
  
  # Calculate standard error
  se_profile <- apply(profile_matrix, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  
  # Create data frame for plotting
  plot_data <- data.frame(
    Position = x_values,
    Signal = mean_profile,
    SE_lower = mean_profile - se_profile,
    SE_upper = mean_profile + se_profile
  )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = Position, y = Signal)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = SE_lower, ymax = SE_upper), alpha = 0.2) +
    labs(
      title = title,
      x = "Distance from CpG Island Center (bp)",
      y = "Signal"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  return(p)
}

# Main function to run the analysis
run_mecp2_analysis <- function(bw_dir, cpg_file, output_dir, upstream = 3000, downstream = 3000, bin_size = 50) {
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Import CpG islands using custom function
  cpg_islands <- import_cpg_islands(cpg_file)
  
  # Define sample groups
  sample_groups <- list(
    neuron_endo = c("NeuM2.bw", "NeuM3.bw"),
    neuron_exo = c("NeuV1.bw", "NeuV2.bw", "NeuV3.bw"),
    nsc_endo = c("NSCM1.bw", "NSCM2.bw", "NSCM3.bw"),
    nsc_exo = c("NSCv1.bw", "NSCv2.bw", "NSCv3.bw")
  )
  
  # Process each group
  profiles <- list()
  
  for(group_name in names(sample_groups)) {
    # Get file paths
    bw_files <- file.path(bw_dir, sample_groups[[group_name]])
    
    # Generate profile
    profile_matrix <- merge_profiles(bw_files, cpg_islands, upstream, downstream, bin_size)
    profiles[[group_name]] <- profile_matrix
    
    # Create individual plot
    p <- plot_profile(profile_matrix, title = paste("MeCP2", group_name), sample_label = group_name)
    ggsave(file.path(output_dir, paste0(group_name, "_profile.png")), p, width = 8, height = 6)
  }
  
  # Create combined plots
  
  # Neuron comparison
  x_values <- seq(-upstream, downstream, length.out = ncol(profiles$neuron_endo))
  neuron_data <- data.frame(
    Position = rep(x_values, 2),
    Signal = c(colMeans(profiles$neuron_endo), colMeans(profiles$neuron_exo)),
    Group = rep(c("Endogenous", "Exogenous"), each = length(x_values))
  )
  
  p_neuron <- ggplot(neuron_data, aes(x = Position, y = Signal, color = Group)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("Endogenous" = "#FF9999", "Exogenous" = "#FF3300")) +
    labs(
      title = "Neuron MeCP2 at CpG Islands",
      x = "Distance from CpG Island Center (bp)",
      y = "Signal"
    ) +
    theme_minimal()
  
  # NSC comparison
  nsc_data <- data.frame(
    Position = rep(x_values, 2),
    Signal = c(colMeans(profiles$nsc_endo), colMeans(profiles$nsc_exo)),
    Group = rep(c("Endogenous", "Exogenous"), each = length(x_values))
  )
  
  p_nsc <- ggplot(nsc_data, aes(x = Position, y = Signal, color = Group)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("Endogenous" = "#83CBEB", "Exogenous" = "#0070C0")) +
    labs(
      title = "NSC MeCP2 at CpG Islands",
      x = "Distance from CpG Island Center (bp)",
      y = "Signal"
    ) +
    theme_minimal()
  
  # Combined cell type comparison
  combined_data <- data.frame(
    Position = rep(x_values, 4),
    Signal = c(colMeans(profiles$neuron_endo), colMeans(profiles$neuron_exo), 
               colMeans(profiles$nsc_endo), colMeans(profiles$nsc_exo)),
    Group = rep(c("Neuron Endo", "Neuron Exo", "NSC Endo", "NSC Exo"), each = length(x_values))
  )
  
  p_combined <- ggplot(combined_data, aes(x = Position, y = Signal, color = Group)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("Neuron Endo" = "#FF9999", "Neuron Exo" = "#FF3300",
                                  "NSC Endo" = "#83CBEB", "NSC Exo" = "#0070C0")) +
    labs(
      title = "MeCP2 at CpG Islands",
      x = "Distance from CpG Island Center (bp)",
      y = "Signal"
    ) +
    theme_minimal()
  
  # Save combined plots
  combined_plot <- p_neuron + p_nsc
  ggsave(file.path(output_dir, "neuron_nsc_comparison.png"), combined_plot, width = 12, height = 6)
  ggsave(file.path(output_dir, "all_samples_comparison.png"), p_combined, width = 8, height = 6)
  
  return(list(profiles = profiles, plots = list(neuron = p_neuron, nsc = p_nsc, combined = p_combined)))
}

# Execute the analysis with provided arguments
result <- run_mecp2_analysis(
  bw_dir = arg_list[["bw-dir"]],
  cpg_file = arg_list[["cpg-file"]],
  output_dir = arg_list[["output-dir"]],
  upstream = 3000,
  downstream = 3000,
  bin_size = 50
)