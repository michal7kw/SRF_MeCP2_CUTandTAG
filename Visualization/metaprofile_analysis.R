# Load required libraries
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(ggplot2)
library(patchwork)
library(GenomicFeatures)

# Read sample sheets
endo_samples <- read.csv("DATA/ENDOGENOUS_sample_sheet_1.csv")
exo_samples <- read.csv("DATA/EXOGENOUS_sample_sheet_1.csv")
all_samples <- rbind(endo_samples, exo_samples)

# Get TxDb object for mouse mm10
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# After getting TxDb object, add CpG Islands data
cpg_islands <- read.table("../DATA/cpg_islands.bed", 
                         col.names=c("chr", "start", "end", "id", "cpg_label", "cpg_count"),
                         header=FALSE)

# Filter for standard chromosomes only and ensure proper chromosome naming
standard_chroms <- paste0("chr", c(1:19, "X", "Y", "M"))
cpg_islands <- cpg_islands[cpg_islands$chr %in% standard_chroms,]

# Function to create windows around regions
createWindowsAroundRegions <- function(gr, upstream, downstream) {
    # Calculate centers
    centers <- (start(gr) + end(gr)) %/% 2
    
    # Create new ranges centered on the middle of each region
    new_ranges <- GRanges(
        seqnames = seqnames(gr),
        ranges = IRanges(
            start = centers - upstream,
            end = centers + downstream
        ),
        strand = strand(gr)
    )
    return(new_ranges)
}

# Create standardized windows around CpG islands
# Create GRanges for CpG islands with standardized windows
cpg_ranges <- GRanges(
    seqnames = cpg_islands$chr,
    ranges = IRanges(start = cpg_islands$start, 
                    end = cpg_islands$end),
    strand = "*",
    type = "CpG_island",  # Add type metadata
    gene_id = paste0("CpG_", seq_len(nrow(cpg_islands))),  # Required for TxDb
    tx_id = paste0("CpG_", seq_len(nrow(cpg_islands))),    # Required for TxDb
    tx_name = paste0("CpG_", seq_len(nrow(cpg_islands)))   # Add tx_name to avoid warning
)

# Set and keep only standard chromosomes
seqlevels(cpg_ranges) <- standard_chroms
seqlevels(cpg_ranges, pruning.mode="coarse") <- standard_chroms

# Create a mock TxDb object for CpG islands
# First create transcript-like features
cpg_transcripts <- GRanges(
    seqnames = seqnames(cpg_ranges),
    ranges = ranges(cpg_ranges),
    strand = strand(cpg_ranges),
    type = "transcript",
    gene_id = mcols(cpg_ranges)$gene_id,
    tx_id = mcols(cpg_ranges)$tx_id,
    tx_name = mcols(cpg_ranges)$tx_name
)

# Combine features
all_features <- c(cpg_ranges, cpg_transcripts)

# Create a mock TxDb object for CpG islands
cpg_txdb <- makeTxDbFromGRanges(all_features)

# Calculate centers of CpG islands and create fixed-width windows
centers <- (start(cpg_ranges) + end(cpg_ranges)) %/% 2
window_size <- 6000  # Total window size (±3000 from center)
cpg_windows <- GRanges(
    seqnames = seqnames(cpg_ranges),
    ranges = IRanges(
        start = centers - (window_size/2),
        width = rep(window_size, length(centers))
    ),
    strand = strand(cpg_ranges)
)

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

# Function to create combined profile plots
create_tissue_profiles <- function(endo_peaks, exo_peaks, tissue_name, color_scheme, feature_type = "TSS") {
    message(paste("Calculating", feature_type, "profiles for", tissue_name))
    
    # Select the appropriate windows based on feature type
    if (feature_type == "TSS") {
        windows <- promoter
        x_label <- "Distance from TSS (bp)"
    } else {
        # Create fixed-width windows centered on CpG islands
        cpg_centers <- (start(cpg_ranges) + end(cpg_ranges)) %/% 2
        temp_windows <- GRanges(
            seqnames = seqnames(cpg_ranges),
            ranges = IRanges(
                start = cpg_centers - 3000,
                end = cpg_centers + 3000
            ),
            strand = strand(cpg_ranges)
        )
        # Use makeBioRegionFromGranges to create compatible windows
        windows <- makeBioRegionFromGranges(temp_windows, upstream = 0, downstream = 0)
        x_label <- "Distance from CpG Island Center (bp)"
    }
    
    # Ensure chromosome naming is consistent
    seqlevels(windows) <- seqlevels(txdb)
    
    # Calculate matrices for endogenous peaks
    endo_matrices <- lapply(endo_peaks, function(x) {
        # Ensure peak chromosomes match window chromosomes
        seqlevels(x) <- seqlevels(windows)
        getTagMatrix(x, windows=windows)
    })
    endo_mean <- colMeans(do.call(rbind, lapply(endo_matrices, colMeans)))
    
    # Calculate matrices for exogenous peaks
    exo_matrices <- lapply(exo_peaks, function(x) {
        # Ensure peak chromosomes match window chromosomes
        seqlevels(x) <- seqlevels(windows)
        getTagMatrix(x, windows=windows)
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
            title = paste("MeCP2 distribution at", feature_type, "\n", tissue_name),
            x = x_label,
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

# Create plots for each tissue type and feature
message("Creating plots for Neurons")
neuron_tss_plot <- create_tissue_profiles(neuron_endo, neuron_exo, "Neurons", "blue", "TSS")
neuron_cpg_plot <- create_tissue_profiles(neuron_endo, neuron_exo, "Neurons", "blue", "CpG Islands")

message("Creating plots for NSCs")
nsc_tss_plot <- create_tissue_profiles(nsc_endo, nsc_exo, "NSCs", "red", "TSS")
nsc_cpg_plot <- create_tissue_profiles(nsc_endo, nsc_exo, "NSCs", "red", "CpG Islands")

# Save individual plots
ggsave("results/neuron_tss_profile.pdf", neuron_tss_plot, width = 6, height = 4)
ggsave("results/neuron_cpg_profile.pdf", neuron_cpg_plot, width = 6, height = 4)
ggsave("results/nsc_tss_profile.pdf", nsc_tss_plot, width = 6, height = 4)
ggsave("results/nsc_cpg_profile.pdf", nsc_cpg_plot, width = 6, height = 4)

# Create combined plot with both TSS and CpG Islands
combined_plot <- (neuron_tss_plot + neuron_cpg_plot) / 
                 (nsc_tss_plot + nsc_cpg_plot)
ggsave("results/combined_profiles.pdf", combined_plot, width = 12, height = 8)

# Perform peak annotation
annotate_peaks <- function(peaks, name) {
    message(paste("Annotating peaks for:", name))
    # Convert peaks to GRanges if needed and annotate
    if (is.list(peaks)) {
        # If it's a list of peaks, merge them first
        merged_peaks <- do.call(c, peaks)
        anno <- annotatePeak(merged_peaks, 
                           TxDb = txdb,
                           annoDb = "org.Mm.eg.db",
                           tssRegion = c(-3000, 3000))
    } else {
        # If it's already a GRanges object
        anno <- annotatePeak(peaks, 
                           TxDb = txdb,
                           annoDb = "org.Mm.eg.db",
                           tssRegion = c(-3000, 3000))
    }
    return(anno)
}

# Create annotation plots
message("Creating annotation plots...")

# Create lists for annotation
neuron_peaks <- list(
    Endogenous = do.call(c, neuron_endo),
    Exogenous = do.call(c, neuron_exo)
)

nsc_peaks <- list(
    Endogenous = do.call(c, nsc_endo),
    Exogenous = do.call(c, nsc_exo)
)

# Annotate peaks
neuron_anno <- lapply(neuron_peaks, function(x) {
    annotatePeak(x, TxDb = txdb, annoDb = "org.Mm.eg.db", tssRegion = c(-3000, 3000))
})

nsc_anno <- lapply(nsc_peaks, function(x) {
    annotatePeak(x, TxDb = txdb, annoDb = "org.Mm.eg.db", tssRegion = c(-3000, 3000))
})

# Create and save annotation plots
neuron_anno_plot <- plotAnnoBar(neuron_anno)
nsc_anno_plot <- plotAnnoBar(nsc_anno)

# Save annotation plots
ggsave("results/neuron_peak_annotation.pdf", neuron_anno_plot, width = 8, height = 6)
ggsave("results/nsc_peak_annotation.pdf", nsc_anno_plot, width = 8, height = 6)

# Save processed data
saveRDS(list(
    neuron_annotation = neuron_anno,
    nsc_annotation = nsc_anno
), "results/chipseeker_analysis.rds")