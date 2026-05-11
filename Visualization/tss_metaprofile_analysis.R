# Load required libraries
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(ggplot2)
library(patchwork)

# Read sample sheets
endo_samples <- read.csv("DATA/ENDOGENOUS_sample_sheet_1.csv")
exo_samples <- read.csv("DATA/EXOGENOUS_sample_sheet_1.csv")
all_samples <- rbind(endo_samples, exo_samples)

# Get TxDb object for mouse mm10
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Function to read peaks and create GRanges object
read_peaks <- function(peak_file, sample_id) {
    message(paste("Reading peaks for:", sample_id, "from:", peak_file))
    peaks <- readPeakFile(peak_file)
    mcols(peaks)$name <- paste0(sample_id, "_peak_", seq_len(length(peaks)))
    return(peaks)
}

# Create list of peak files with proper sample IDs
peak_list <- list()
for(i in 1:nrow(all_samples)) {
    tryCatch({
        peak_list[[all_samples$SampleID[i]]] <- read_peaks(
            all_samples$Peaks[i],
            all_samples$SampleID[i]
        )
    }, error = function(e) {
        message(paste("Error reading peaks for:", all_samples$SampleID[i]))
        message(e)
    })
}

# Generate promoter regions for TSS analysis
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

# Function to create TSS profile plots
create_tss_profiles <- function(endo_peaks, exo_peaks, tissue_name, color_scheme) {
    message(paste("Calculating TSS profiles for", tissue_name))
    
    # Calculate matrices for endogenous peaks
    endo_matrices <- lapply(endo_peaks, function(x) {
        getTagMatrix(x, windows=promoter)
    })
    endo_mean <- colMeans(do.call(rbind, lapply(endo_matrices, colMeans)))
    
    # Calculate matrices for exogenous peaks
    exo_matrices <- lapply(exo_peaks, function(x) {
        getTagMatrix(x, windows=promoter)
    })
    exo_mean <- colMeans(do.call(rbind, lapply(exo_matrices, colMeans)))
    
    # Create data frame for plotting
    pos <- seq(-3000, 3000, length.out = length(endo_mean))
    df <- data.frame(
        Position = rep(pos, 2),
        Count = c(endo_mean, exo_mean),
        Type = factor(rep(c("Endogenous", "Exogenous"), each = length(pos)))
    )
    
    # Set colors
    if (color_scheme == "blue") {
        colors <- c("Endogenous" = "lightblue", "Exogenous" = "darkblue")
    } else {
        colors <- c("Endogenous" = "pink", "Exogenous" = "red")
    }
    
    # Create plot
    p <- ggplot(df, aes(x = Position, y = Count, color = Type)) +
        geom_line(linewidth = 1) +
        scale_color_manual(values = colors) +
        theme_minimal() +
        labs(
            title = paste("MeCP2 distribution at TSS\n", tissue_name),
            x = "Distance from TSS (bp)",
            y = "Average Signal"
        ) +
        theme(
            plot.title = element_text(hjust = 0.5, size = 14),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12),
            legend.position = "top",
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()
        )
    
    return(p)
}

# Separate samples by tissue type and condition
neuron_endo <- peak_list[all_samples$Tissue == "Neuron" & all_samples$Condition == "Endogenous"]
neuron_exo <- peak_list[all_samples$Tissue == "Neuron" & all_samples$Condition == "Exogenous"]
nsc_endo <- peak_list[all_samples$Tissue == "NSC" & all_samples$Condition == "Endogenous"]
nsc_exo <- peak_list[all_samples$Tissue == "NSC" & all_samples$Condition == "Exogenous"]

# Create output directory if it doesn't exist
dir.create("results", showWarnings = FALSE)

# Create plots for each tissue type
message("Creating plots for Neurons")
neuron_tss_plot <- create_tss_profiles(neuron_endo, neuron_exo, "Neurons", "blue")
message("Creating plots for NSCs")
nsc_tss_plot <- create_tss_profiles(nsc_endo, nsc_exo, "NSCs", "red")

# Save individual plots
ggsave("results/neuron_tss_profile.pdf", neuron_tss_plot, width = 6, height = 4)
ggsave("results/nsc_tss_profile.pdf", nsc_tss_plot, width = 6, height = 4)

# Create combined plot
combined_plot <- neuron_tss_plot / nsc_tss_plot
ggsave("results/combined_tss_profiles.pdf", combined_plot, width = 8, height = 8) 