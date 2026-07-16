# 13. SAVE DATA
###############################################################################

data.table::fwrite(region_metrics, file.path(OUT_DIR, "Revision_R1_3_region_level_HSF_metrics_all_thresholds.csv"))
data.table::fwrite(region_metrics, file.path(OUT_DIR, "S8C_region_metrics.tsv.gz"), sep = "\t")
data.table::fwrite(continuous_metrics, file.path(OUT_DIR, "Revision_R1_3_region_level_continuous_HSF_scores.csv"))
data.table::fwrite(threshold_stats, file.path(OUT_DIR, "Revision_R1_3_threshold_sensitivity_statistics.csv"))
data.table::fwrite(continuous_stats, file.path(OUT_DIR, "Revision_R1_3_continuous_score_statistics.csv"))

motif_inventory <- data.frame(HSF_model = names(pwm_models), JASPAR_ID = vapply(pwm_models, `[[`, character(1), "id"),
                              motif_width = vapply(pwm_models, `[[`, integer(1), "width"), pseudocount = PSEUDOCOUNT)

qc <- data.frame(item = c("Genome file", "JASPAR file", "Background coordinate file", "Chromosome source columns available",
                          "ONSEN regions", "Background regions before exclusion", "Direct overlaps removed", "Final background regions",
                          "HSF models", "Thresholds", "Pseudocount", "Both strands scanned"),
                 value = c(basename(GENOME_FILE), basename(JASPAR_FILE), basename(BG_FILE), paste(intersect(c("chr_clean", "seqid"), names(bg_raw)), collapse = ", "),
                           nrow(onsen), n_before, n_removed, nrow(bg), length(pwm_models), paste(THRESHOLDS, collapse = ", "),
                           PSEUDOCOUNT, "Yes"))

data.table::fwrite(qc, file.path(OUT_DIR, "Revision_R1_3_analysis_QC.csv"))
data.table::fwrite(motif_inventory, file.path(OUT_DIR, "Revision_R1_3_HSF_model_inventory.csv"))

###############################################################################
# 14. RECONSTRUCT FINAL TABLE S8
###############################################################################

TABLE_FILE <- file.path(OUT_DIR, "Table_S8.xlsx")
wb <- openxlsx::createWorkbook(creator = "Baibhav R. Barbaruah")

write_sheet <- function(name, title, data) {
  openxlsx::addWorksheet(wb, name); openxlsx::writeData(wb, name, title, startRow = 1)
  openxlsx::writeDataTable(wb, name, data, startRow = 3, tableStyle = "TableStyleMedium2")
  openxlsx::freezePane(wb, name, firstActiveRow = 4); openxlsx::setColWidths(wb, name, 1:ncol(data), widths = "auto")
}

write_sheet("S8A_threshold_stats", "Table S8A. HSF-family motif-density threshold sensitivity.", threshold_stats)
write_sheet("S8B_continuous_stats", "Table S8B. Continuous HSF relative-score comparisons.", continuous_stats)
write_sheet("S8C_region_metrics", "Table S8C. Region-level threshold and continuous-score metrics.", region_metrics)
write_sheet("S8D_QC", "Table S8D. Analysis quality-control information.", dplyr::bind_rows(
  qc |> dplyr::mutate(section = "Analysis QC"),
  motif_inventory |> dplyr::transmute(section = "HSF motif inventory", item = HSF_model,
                                      value = paste0(JASPAR_ID, "; width=", motif_width, " bp; pseudocount=", pseudocount))
) |> dplyr::select(section, item, value))

openxlsx::saveWorkbook(wb, TABLE_FILE, overwrite = TRUE)

###############################################################################
# 15. RECONSTRUCT FINAL FIGURE S3
###############################################################################

plot_summary <- region_metrics |>
  dplyr::group_by(class, threshold) |>
  dplyr::summarise(median = median(HSF_hits_per_kb), Q1 = quantile(HSF_hits_per_kb, 0.25),
                   Q3 = quantile(HSF_hits_per_kb, 0.75), .groups = "drop")

theme_revision <- function() {
  ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = "black", size = 12),
                   axis.title = ggplot2::element_text(colour = "black", size = 14, face = "bold"),
                   legend.text = ggplot2::element_text(size = 11), legend.title = ggplot2::element_blank(),
                   plot.tag = ggplot2::element_text(size = 16, face = "bold"))
}

panel_a <- ggplot2::ggplot(plot_summary, ggplot2::aes(threshold, median, colour = class, fill = class, group = class)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = Q1, ymax = Q3), alpha = 0.18, colour = NA) +
  ggplot2::geom_line(linewidth = 1) + ggplot2::geom_point(shape = 21, size = 3.3, colour = "black") +
  ggplot2::scale_colour_manual(values = CLASS_COLOURS) + ggplot2::scale_fill_manual(values = CLASS_COLOURS) +
  ggplot2::scale_x_continuous(breaks = THRESHOLDS) +
  ggplot2::labs(x = "Relative PWM-score threshold", y = "Median HSF motif-position hits per kb") +
  theme_revision() + ggplot2::theme(legend.position = "top")

set.seed(20260715)

panel_b <- ggplot2::ggplot(continuous_metrics, ggplot2::aes(class, maximum_HSF_relative_score, fill = class)) +
  ggplot2::geom_violin(trim = TRUE, alpha = 0.60, colour = "black") +
  ggplot2::geom_boxplot(width = 0.22, outlier.shape = NA) +
  ggplot2::geom_jitter(data = continuous_metrics |> dplyr::filter(class == ONSEN_CLASS), width = 0.08, shape = 21, size = 2.4) +
  ggplot2::scale_fill_manual(values = CLASS_COLOURS) +
  ggplot2::scale_x_discrete(labels = setNames(c("Strict TE\nbackground", "ONSEN\nwindows"), c(BG_CLASS, ONSEN_CLASS))) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(x = NULL, y = "Maximum HSF relative PWM score") +
  theme_revision() + ggplot2::theme(legend.position = "none")

panel_c <- ggplot2::ggplot(continuous_metrics, ggplot2::aes(class, mean_top5_HSF_relative_score, fill = class)) +
  ggplot2::geom_violin(trim = TRUE, alpha = 0.60, colour = "black") +
  ggplot2::geom_boxplot(width = 0.22, outlier.shape = NA) +
  ggplot2::geom_jitter(data = continuous_metrics |> dplyr::filter(class == ONSEN_CLASS), width = 0.08, shape = 21, size = 2.4) +
  ggplot2::scale_fill_manual(values = CLASS_COLOURS) +
  ggplot2::scale_x_discrete(labels = setNames(c("Strict TE\nbackground", "ONSEN\nwindows"), c(BG_CLASS, ONSEN_CLASS))) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(x = NULL, y = "Mean of top five HSF relative scores") +
  theme_revision() + ggplot2::theme(legend.position = "none")

figure_s3 <- (panel_a / (panel_b | panel_c)) + patchwork::plot_annotation(tag_levels = "A")

# Mandatory for the user's RStudio workflow: display before saving.
print(figure_s3)

FIG_PDF <- file.path(FIG_DIR, "FigS3_HSF_threshold_and_continuous_scores.pdf")
FIG_PNG <- file.path(FIG_DIR, "FigS3_HSF_threshold_and_continuous_scores.png")

ggplot2::ggsave(FIG_PDF, figure_s3, width = 10.5, height = 10, units = "in")
ggplot2::ggsave(FIG_PNG, figure_s3, width = 10.5, height = 10, units = "in", dpi = 600, bg = "white")

###############################################################################
# 16. VALIDATION AND FINAL OUTPUT
###############################################################################

capture.output(sessionInfo(), file = file.path(OUT_DIR, "Revision_R1_3_sessionInfo.txt"))

if (!file.exists(TABLE_FILE) || !file.exists(FIG_PDF) || !file.exists(FIG_PNG)) stop("One or more final output files were not created.")
if (nrow(region_metrics[region_metrics$class == ONSEN_CLASS & region_metrics$threshold == 0.85, ]) != 16L) stop("ONSEN row validation failed.")
if (!all(sort(unique(region_metrics$threshold)) == THRESHOLDS)) stop("Threshold validation failed.")

cat("\n============================================================\n")
cat("REVIEWER 1.3 ANALYSIS COMPLETED\n")
cat("============================================================\n\n")

cat("Threshold statistics:\n"); print(threshold_stats, row.names = FALSE)
cat("\nContinuous-score statistics:\n"); print(continuous_stats, row.names = FALSE)

cat("\nOutputs:\n")
cat("Table S8 reconstructed:\n  ", TABLE_FILE, "\n", sep = "")
cat("Figure S3 PDF:\n  ", FIG_PDF, "\n", sep = "")
cat("Figure S3 PNG:\n  ", FIG_PNG, "\n", sep = "")
cat("Source-data directory:\n  ", OUT_DIR, "\n", sep = "")
