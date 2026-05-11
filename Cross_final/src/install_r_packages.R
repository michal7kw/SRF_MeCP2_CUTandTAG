#!/usr/bin/env Rscript

# Function to install packages if not already installed
install_if_missing <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    if (package_name %in% c("ChIPseeker", "TxDb.Mmusculus.UCSC.mm10.knownGene", 
                           "org.Mm.eg.db", "clusterProfiler", "GenomicFeatures")) {
      BiocManager::install(package_name)
    } else {
      install.packages(package_name)
    }
  }
}

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of required packages
packages <- c(
  "ChIPseeker",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "org.Mm.eg.db",
  "clusterProfiler",
  "GenomicFeatures",
  "rtracklayer",
  "ggplot2"
)

# Install missing packages
for (pkg in packages) {
  install_if_missing(pkg)
}
