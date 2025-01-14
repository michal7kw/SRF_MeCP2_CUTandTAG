library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)  # For gene annotations
library(clusterProfiler)
library(ggplot2)
library(gridExtra)

# Parse command line arguments
library(optparse)

# Create option parser with explicit dest parameter
option_list <- list(
    make_option("--work-dir", 
                type="character", 
                default=NULL, 
                dest="work_dir",
                help="Working directory path")
)

# Parse arguments
opts <- parse_args(OptionParser(option_list=option_list))

# Add diagnostic prints
print("Command line arguments received:")
print(opts)
print("Working directory path:")
print(opts$work_dir)

# Check if work-dir argument is provided
if (is.null(opts$work_dir) || !is.character(opts$work_dir)) {
    stop("Invalid or missing work directory path")
}

# Set working directory
message(sprintf("Setting working directory to: %s", opts$work_dir))
setwd(opts$work_dir)


# Read peak files
read_peak_file <- function(file_path) {
  peaks <- read_csv(file_path) %>%
    # Convert to GRanges object
    makeGRangesFromDataFrame(
      keep.extra.columns = TRUE,
      seqnames.field = "chr",
      start.field = "start",
      end.field = "end"
    )
  return(peaks)
}

# Read the three peak files
endo_only_peaks <- read_peak_file("cpg_enrichment_endo_only.csv")
exo_only_peaks <- read_peak_file("cpg_enrichment_exo_only.csv")
both_peaks <- read_peak_file("cpg_enrichment_both.csv")

# Get TSS annotations
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Create peaks_annotation directory if it doesn't exist
dir.create("peaks_annotation", showWarnings = FALSE)

# Analyze peak distribution relative to TSS
analyze_peaks <- function(peaks, name) {
  # Annotate peaks
  peakAnno <- annotatePeak(peaks, 
                          tssRegion = c(-3000, 3000),  # Define TSS region as ±3kb
                          TxDb = txdb,
                          annoDb = "org.Mm.eg.db")
  
  # Plot peak distribution
  pdf(paste0("peaks_annotation/", name, "_distribution.pdf"))
  plotAnnoPie(peakAnno)
  plotDistToTSS(peakAnno)
  dev.off()
  
  # Return detailed annotation
  return(peakAnno)
}

# Analyze each peak set
endo_anno <- analyze_peaks(endo_only_peaks, "endo_only")
exo_anno <- analyze_peaks(exo_only_peaks, "exo_only")
both_anno <- analyze_peaks(both_peaks, "both")

# Extract peaks near TSS (within ±3kb)
get_tss_peaks <- function(peaks, anno) {
  anno_df <- as.data.frame(anno@anno)
  tss_peaks <- anno_df %>%
    filter(abs(distanceToTSS) <= 3000) %>%
    select(seqnames, start, end, distanceToTSS, geneId, SYMBOL)
  return(tss_peaks)
}

# Get TSS-proximal peaks for each category
endo_tss <- get_tss_peaks(endo_only_peaks, endo_anno)
exo_tss <- get_tss_peaks(exo_only_peaks, exo_anno)
both_tss <- get_tss_peaks(both_peaks, both_anno)

# Save TSS-proximal peaks
write_csv(endo_tss, "peaks_annotation/endo_only_tss_peaks.csv")
write_csv(exo_tss, "peaks_annotation/exo_only_tss_peaks.csv")
write_csv(both_tss, "peaks_annotation/both_tss_peaks.csv")

# Generate summary statistics
summary_stats <- data.frame(
  Category = c("Endo Only", "Exo Only", "Both"),
  Total_Peaks = c(length(endo_only_peaks), length(exo_only_peaks), length(both_peaks)),
  TSS_Proximal = c(nrow(endo_tss), nrow(exo_tss), nrow(both_tss))
)

summary_stats <- summary_stats %>%
  mutate(Percent_TSS = (TSS_Proximal / Total_Peaks) * 100)

write_csv(summary_stats, "peaks_annotation/peak_location_summary.csv")

# Plot distribution of distances to TSS
plot_tss_dist <- function(anno_list, names) {
  distances <- map2_dfr(anno_list, names, ~{
    as.data.frame(.x@anno) %>%
      select(distanceToTSS) %>%
      mutate(Category = .y)
  })
  
  ggplot(distances, aes(x = distanceToTSS, fill = Category)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    xlim(-3000, 3000) +
    labs(x = "Distance to TSS (bp)",
         y = "Density",
         title = "Distribution of Peak Distances to TSS")
  
  ggsave("peaks_annotation/tss_distance_distribution.pdf", width = 10, height = 6)
}

plot_tss_dist(
  list(endo_anno, exo_anno, both_anno),
  c("Endo Only", "Exo Only", "Both")
)

#####################################################################################################################################

# Enhanced genomic annotation analysis
analyze_genomic_regions <- function(peaks, name) {
    # Create output directory for plots
    dir.create("peaks_annotation_adv", showWarnings = FALSE)
    
    # 1. Detailed genomic annotation
    detailed_anno <- annotatePeak(peaks,
                                TxDb = txdb,
                                annoDb = "org.Mm.eg.db",
                                tssRegion = c(-3000, 3000),
                                level = "transcript",
                                verbose = FALSE)
    
    # 2. Generate comprehensive plots
    
    # Pie chart of genomic annotations
    pdf(paste0("peaks_annotation_adv/", name, "_pie_chart.pdf"))
    plotAnnoPie(detailed_anno)
    dev.off()
    
    # 3. Generate detailed statistics
    anno_stats <- as.data.frame(detailed_anno@anno) %>%
        group_by(annotation) %>%
        summarise(
            count = n(),
            percentage = n() / nrow(.) * 100
        )
    
    # 4. Get promoter vs non-promoter stats
    promoter_stats <- data.frame(
        region = c("Promoter", "Non-promoter"),
        count = c(
            sum(anno_stats$count[grep("Promoter", anno_stats$annotation)]),
            sum(anno_stats$count[!grepl("Promoter", anno_stats$annotation)])
        )
    )
    promoter_stats$percentage <- promoter_stats$count / sum(promoter_stats$count) * 100
    
    # Return statistics
    return(list(
        detailed_stats = anno_stats,
        promoter_stats = promoter_stats,
        full_annotation = detailed_anno
    ))
}

# Analyze each peak set
endo_analysis <- analyze_genomic_regions(endo_only_peaks, "endo_only")
exo_analysis <- analyze_genomic_regions(exo_only_peaks, "exo_only")
both_analysis <- analyze_genomic_regions(both_peaks, "both")

# Create comparative analysis
create_comparative <- function(endo_analysis, exo_analysis, both_analysis) {
    # Combine statistics for comparison
    combined_stats <- bind_rows(
        endo_analysis$detailed_stats %>% mutate(category = "Endo Only"),
        exo_analysis$detailed_stats %>% mutate(category = "Exo Only"),
        both_analysis$detailed_stats %>% mutate(category = "Both")
    )
    
    # Create summary table
    summary_table <- combined_stats %>%
        spread(category, percentage) %>%
        arrange(desc(`Both`))
    
    write_csv(summary_table, "peaks_annotation_adv/region_distribution_summary.csv")
}

create_comparative(endo_analysis, exo_analysis, both_analysis)

# Modified feature enrichment analysis
calculate_feature_enrichment <- function(peaks, name) {
    # Get peak annotations
    peak_anno <- annotatePeak(peaks,
                            TxDb = txdb,
                            annoDb = "org.Mm.eg.db",
                            tssRegion = c(-3000, 3000))
    
    # Extract gene IDs from peaks
    peak_genes <- unique(peak_anno@anno$geneId)
    
    # Get all genes as background
    all_genes <- keys(org.Mm.eg.db)
    
    # Perform GO enrichment analysis
    ego <- enrichGO(gene = peak_genes,
                   universe = all_genes,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)
    
    # Plot results
    if (!is.null(ego) && nrow(ego) > 0) {
        pdf(paste0("peaks_annotation_adv/", name, "_GO_enrichment.pdf"))
        print(dotplot(ego, showCategory = 15, title = paste(name, "GO Enrichment")))
        dev.off()
        
        # Save results to CSV
        write.csv(as.data.frame(ego), 
                 file = paste0("peaks_annotation_adv/", name, "_GO_enrichment.csv"))
    }
    
    return(ego)
}

# Calculate enrichment for each peak set
endo_enrichment <- calculate_feature_enrichment(endo_only_peaks, "endo_only")
exo_enrichment <- calculate_feature_enrichment(exo_only_peaks, "exo_only")
both_enrichment <- calculate_feature_enrichment(both_peaks, "both")

# Save all statistics
write_csv(endo_analysis$detailed_stats, "peaks_annotation_adv/endo_only_detailed_stats.csv")
write_csv(exo_analysis$detailed_stats, "peaks_annotation_adv/exo_only_detailed_stats.csv")
write_csv(both_analysis$detailed_stats, "peaks_annotation_adv/both_detailed_stats.csv")