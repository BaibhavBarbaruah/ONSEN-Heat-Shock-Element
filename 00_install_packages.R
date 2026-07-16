# Install/check packages required by the complete flat repository.

cran_packages <- c(
  "data.table", "dplyr", "tidyr", "readr", "stringr", "purrr", "tibble",
  "ggplot2", "ggseqlogo", "forcats", "scales", "patchwork", "matrixStats",
  "openxlsx", "readxl"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_packages <- c(
  "Biostrings", "GenomicRanges", "IRanges", "S4Vectors",
  "rtracklayer", "DESeq2", "Rsubread"
)
missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("Package installation/check completed.")
