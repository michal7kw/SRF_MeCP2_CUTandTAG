library(optparse)
library(GenomicRanges)
library(rtracklayer)

# Parse command line arguments
option_list <- list(
    make_option("--peaks", type="character", help="Comma-separated list of peak files"),
    make_option("--output", type="character", help="Output file path")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Read peak files
peak_files <- unlist(strsplit(opt$peaks, " "))
peak_list <- lapply(peak_files, import)
names(peak_list) <- basename(peak_files)

# Calculate overlap between peaks
calculate_overlap <- function(gr1, gr2) {
    overlap <- findOverlaps(gr1, gr2)
    n_overlap <- length(unique(queryHits(overlap)))
    return(n_overlap / length(gr1))
}

# Calculate pairwise reproducibility
results <- data.frame()
for(i in 1:(length(peak_list)-1)) {
    for(j in (i+1):length(peak_list)) {
        overlap <- calculate_overlap(peak_list[[i]], peak_list[[j]])
        results <- rbind(results, data.frame(
            Sample1 = names(peak_list)[i],
            Sample2 = names(peak_list)[j],
            Reproducibility = overlap
        ))
    }
}

# Write results
write.table(results, opt$output, sep="\t", quote=FALSE, row.names=FALSE) 