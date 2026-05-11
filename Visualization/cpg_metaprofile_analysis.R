# Load required libraries
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(ggplot2)
library(patchwork)
library(GenomicRanges)

# Read sample sheets
endo_samples <- read.csv("DATA/ENDOGENOUS_sample_sheet_1.csv")
exo_samples <- read.csv("DATA/EXOGENOUS_sample_sheet_1.csv")
all_samples <- rbind(endo_samples, exo_samples)

# Get TxDb object for mouse mm10
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Read CpG Islands data
cpg_islands <- read.table("../DATA/cpg_islands.bed", 
                         col.names=c("chr", "start", "end", "id", "cpg_label", "cpg_count"),
                         header=FALSE)

# Filter for standard chromosomes only
standard_chroms <- paste0("chr", c(1:19, "X", "Y", "M"))
cpg_islands <- cpg_islands[cpg_islands$chr %in% standard_chroms,]

# Create GRanges for CpG islands
cpg_ranges <- GRanges(
    seqnames = cpg_islands$chr,
    ranges = IRanges(start = cpg_islands$start, 
                    end = cpg_islands$end),
    strand = "*"
)

# Set and keep only standard chromosomes
seqlevels(cpg_ranges) <- standard_chroms
seqlevels(cpg_ranges, pruning.mode="coarse") <- standard_chroms

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

# Function to create CpG profile plots
create_cpg_profiles <- function(endo_peaks, exo_peaks, tissue_name, color_scheme) {
    message(paste("Calculating CpG Islands profiles for", tissue_name))
    
    # Create bioRegion from CpG islands
    cpg_bioregion <- makeBioRegionFromGranges(
        gr = cpg_ranges,
        upstream = 3000,
        downstream = 3000,
        type = "body",
        by = "center"
    )
    
    # Define standard chromosomes
    standard_chroms <- paste0("chr", c(1:19, "X", "Y", "M"))
    
    # Calculate matrices for endogenous peaks
    endo_matrices <- lapply(endo_peaks, function(x) {
        # Keep only standard chromosomes
        x <- keepSeqlevels(x, standard_chroms, pruning.mode="coarse")
        seqlevels(x) <- seqlevels(cpg_bioregion)
        getTagMatrix(x, windows = cpg_bioregion, nbin = 100)
    })
    endo_mean <- colMeans(do.call(rbind, lapply(endo_matrices, colMeans)))
    
    # Calculate matrices for exogenous peaks
    exo_matrices <- lapply(exo_peaks, function(x) {
        # Keep only standard chromosomes
        x <- keepSeqlevels(x, standard_chroms, pruning.mode="coarse")
        seqlevels(x) <- seqlevels(cpg_bioregion)
        getTagMatrix(x, windows = cpg_bioregion, nbin = 100)
    })
    exo_mean <- colMeans(do.call(rbind, lapply(exo_matrices, colMeans)))
    
    # Create data frame for plotting
    pos <- seq(-3000, 3000, length.out = 100)
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
            title = paste("MeCP2 distribution at CpG Islands\n", tissue_name),
            x = "Distance from CpG Island Center (bp)",
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
neuron_cpg_plot <- create_cpg_profiles(neuron_endo, neuron_exo, "Neurons", "blue")
message("Creating plots for NSCs")
nsc_cpg_plot <- create_cpg_profiles(nsc_endo, nsc_exo, "NSCs", "red")

# Save individual plots
ggsave("results/neuron_cpg_profile.pdf", neuron_cpg_plot, width = 6, height = 4)
ggsave("results/nsc_cpg_profile.pdf", nsc_cpg_plot, width = 6, height = 4)

# Create combined plot
combined_plot <- neuron_cpg_plot / nsc_cpg_plot
ggsave("results/combined_cpg_profiles.pdf", combined_plot, width = 8, height = 8) 