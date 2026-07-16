# Reproduce Fig. 7C sequence identity and Fig. 7D representative motif logos.
# This script recreates the final Fig. 7C and Fig. 7D panels; journal-ready assembled artwork remains a submission asset.

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "ggplot2", "ggseqlogo", "patchwork"))
message_config()

variant_file <- repo_file("source_data/Figure7_variant_metrics.tsv")
hsf_file <- repo_file("motifs/MA1667.2_HSFC1.jaspar")
dof_file <- repo_file("motifs/MA0981.2_DOF1.8.jaspar")
metadata_file <- repo_file("source_data/Figure7_logo_model_metadata.csv")

variants <- data.table::fread(variant_file, data.table = FALSE, check.names = FALSE)
required_variant_columns <- c("variant_id", "mismatch_count", "sequence_identity_percent")
if (!all(required_variant_columns %in% names(variants))) {
  stop("Figure 7 variant source is missing: ", paste(setdiff(required_variant_columns, names(variants)), collapse = ", "), call. = FALSE)
}
if (nrow(variants) != 53L || sum(variants$mismatch_count == 0) != 1L) {
  stop("Unexpected Figure 7 variant source dimensions.", call. = FALSE)
}

variants$variant_id <- factor(variants$variant_id, levels = rev(variants$variant_id))
identity_min <- min(variants$sequence_identity_percent, na.rm = TRUE)

p7c <- ggplot2::ggplot(variants, ggplot2::aes(sequence_identity_percent, variant_id)) +
  ggplot2::geom_segment(ggplot2::aes(xend = 100, yend = variant_id), colour = "grey78", linewidth = 0.45) +
  ggplot2::geom_point(shape = 21, size = 3.0, stroke = 0.55, fill = "#C77DD4", colour = "black") +
  ggplot2::geom_text(ggplot2::aes(label = paste0(mismatch_count, " nt")), hjust = -0.15, size = 3.1) +
  ggplot2::scale_x_continuous(
    breaks = seq(floor(identity_min / 2) * 2, 100, 2), limits = c(floor(identity_min / 2) * 2, 101.8),
    labels = function(x) paste0(x, "%"), expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(x = "Sequence identity relative to Col-0", y = NULL, tag = "C") +
  theme_onsen(12) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 8),
    plot.tag = ggplot2::element_text(face = "bold", size = 18, margin = ggplot2::margin(r = 14, b = 10)),
    plot.tag.position = "topleft", plot.tag.location = "margin",
    plot.margin = ggplot2::margin(t = 26, r = 18, b = 18, l = 34)
  )

parse_one_jaspar <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header <- sub("^>", "", lines[grepl("^>", lines)][1])
  parts <- strsplit(header, "[\\t ]+")[[1]]
  rows <- lapply(c("A", "C", "G", "T"), function(base) {
    line <- lines[grepl(paste0("^", base, "[[:space:]]+\\["), lines)][1]
    as.numeric(strsplit(trimws(gsub(".*\\[|\\].*", "", line)), "[[:space:]]+")[[1]])
  })
  names(rows) <- c("A", "C", "G", "T")
  pfm <- do.call(rbind, rows)
  rownames(pfm) <- names(rows)
  list(id = parts[[1]], name = parts[[2]], pfm = pfm)
}

hsf <- parse_one_jaspar(hsf_file)
dof <- parse_one_jaspar(dof_file)
if (!identical(c(hsf$id, dof$id), c("MA1667.2", "MA0981.2")) || !identical(c(hsf$name, dof$name), c("HSFC1", "DOF1.8"))) {
  stop("Representative logo metadata does not match the final Fig. 7D models.", call. = FALSE)
}
logo_metadata <- read.csv(metadata_file, check.names = FALSE)
if (!all(c("MA1667.2", "MA0981.2") %in% logo_metadata$JASPAR_ID)) stop("Figure 7 logo metadata is incomplete.", call. = FALSE)

p_hsf <- ggseqlogo::ggseqlogo(hsf$pfm, method = "bits", seq_type = "dna") +
  ggplot2::labs(title = paste0(hsf$name, " (", hsf$id, ")"), x = NULL, y = "Information content (bits)", tag = "D") +
  theme_onsen(11) + ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 18), plot.tag.position = "topleft")
p_dof <- ggseqlogo::ggseqlogo(dof$pfm, method = "bits", seq_type = "dna") +
  ggplot2::labs(title = paste0(dof$name, " (", dof$id, ")"), x = "Motif position", y = "Information content (bits)") + theme_onsen(11)
p7d <- (p_hsf / p_dof) & ggplot2::theme(plot.margin = ggplot2::margin(10, 12, 10, 16))

if (ONSEN_MAKE_FIGURES) {
  save_plot_pair(p7c, "Fig7C_sequence_identity", 8.0, 10.0)
  save_plot_pair(p7d, "Fig7D_representative_JASPAR_logos", 7.0, 7.6)
}

data.table::fwrite(
  data.frame(
    family = c("HSF", "DOF"), model_name = c(hsf$name, dof$name), JASPAR_ID = c(hsf$id, dof$id),
    motif_width = c(ncol(hsf$pfm), ncol(dof$pfm)), stringsAsFactors = FALSE
  ),
  out_file("Figure7_logo_model_metadata_reproduced.csv")
)

message("Figure 7C identity and Figure 7D logo workflow completed.")
