# 13. SAVE DATA
###############################################################################

data.table::fwrite(region_metrics, file.path(OUT_DIR, "S7_region_level_HSF_metrics_all_thresholds.csv"))
data.table::fwrite(region_metrics, file.path(OUT_DIR, "S7D_region_threshold.tsv.gz"), sep = "\t")
data.table::fwrite(continuous_metrics, file.path(OUT_DIR, "S7_region_level_continuous_HSF_scores.csv"))
data.table::fwrite(threshold_stats, file.path(OUT_DIR, "S7_threshold_sensitivity_statistics.csv"))
data.table::fwrite(continuous_stats, file.path(OUT_DIR, "S7_continuous_score_statistics.csv"))

motif_inventory <- data.frame(HSF_model = names(pwm_models), JASPAR_ID = vapply(pwm_models, `[[`, character(1), "id"),
                              motif_width = vapply(pwm_models, `[[`, integer(1), "width"), pseudocount = PSEUDOCOUNT)

qc <- data.frame(
  item = c(
    "Primary HSF count definition", "PWM-score thresholds",
    "Final ONSEN terminal windows", "Final strict non-ONSEN TE background regions",
    "Continuous-score ONSEN handling", "Continuous-score background handling",
    "Archived background regions removed", "Final continuous-score universe"
  ),
  value = c(
    paste(
      "Exact coordinate-unique placements: predictions with identical 1-based",
      "start-end coordinates are collapsed across HSF PWM models and strands;",
      "different coordinates remain separate."
    ),
    paste(THRESHOLDS, collapse = ", "),
    nrow(onsen), nrow(bg),
    "The same curated set of 16 ONSEN terminal windows is used for threshold and continuous-score analyses.",
    "Background continuous-score rows are restricted by region ID to the final 1,930-region universe.",
    n_harmonisation_removed,
    "16 ONSEN terminal windows + 1,930 strict non-ONSEN TE background regions."
  ),
  stringsAsFactors = FALSE
)

original_qc <- data.frame(
  item = c(
    "Genome file", "JASPAR file", "Background coordinate file",
    "Chromosome source columns available", "ONSEN regions",
    "Background regions in original scan", "Direct coordinate overlaps removed in original scan",
    "Regions excluded by final annotation and ONSEN/ATCOPIA78 harmonisation",
    "Final background regions used in reported analyses", "HSF models",
    "Thresholds", "Pseudocount", "Both strands scanned", "Reconciliation"
  ),
  value = c(
    paste0(basename(GENOME_FILE), " (external input; see repository INPUT_PROVENANCE.tsv)"),
    paste0(basename(JASPAR_FILE), " (external input; see repository INPUT_PROVENANCE.tsv)"),
    paste0(basename(BG_FILE), " (external input; see repository INPUT_PROVENANCE.tsv)"),
    paste(intersect(c("chr_clean", "seqid"), names(bg_raw)), collapse = ", "),
    nrow(onsen), n_before, n_removed, n_harmonisation_removed, nrow(bg),
    length(pwm_models), paste(THRESHOLDS, collapse = ", "), PSEUDOCOUNT, "Yes",
    paste(
      "Reported threshold and continuous-score statistics use the final 1,930-region universe;",
      "the 1,942-region value is retained only as archived scan provenance."
    )
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(qc, file.path(OUT_DIR, "S7_analysis_QC.csv"))
data.table::fwrite(motif_inventory, file.path(OUT_DIR, "S7_HSF_model_inventory.csv"))

###############################################################################
# 14. RECONSTRUCT FINAL TABLE S7
###############################################################################

TABLE_FILE <- file.path(OUT_DIR, "Table_S7.xlsx")
wb <- openxlsx::createWorkbook(creator = "Baibhav R. Barbaruah")

write_sheet <- function(name, title, data) {
  openxlsx::addWorksheet(wb, name); openxlsx::writeData(wb, name, title, startRow = 1)
  openxlsx::writeDataTable(wb, name, data, startRow = 3, tableStyle = "TableStyleMedium2")
  openxlsx::freezePane(wb, name, firstActiveRow = 4); openxlsx::setColWidths(wb, name, 1:ncol(data), widths = "auto")
}

threshold_table <- threshold_stats |>
  dplyr::transmute(
    `Relative PWM-score threshold` = relative_score_threshold,
    `n ONSEN terminal windows` = n_ONSEN_regions,
    `n strict non-ONSEN TE regions` = n_background_regions,
    `ONSEN median density (placements/kb)` = ONSEN_median_hits_per_kb,
    `Background median density (placements/kb)` = background_median_hits_per_kb,
    `Wilcoxon W` = Wilcoxon_W,
    `Two-sided Wilcoxon P` = Wilcoxon_P,
    `BH-adjusted P` = Wilcoxon_P_BH,
    `Cliff's delta` = Cliffs_delta,
    Significance = dplyr::if_else(Wilcoxon_P_BH < 0.0001, "****", "ns")
  )

threshold_summary <- region_metrics |>
  dplyr::group_by(threshold, class) |>
  dplyr::summarise(
    n = dplyr::n(), median = median(HSF_hits_per_kb),
    Q1 = unname(stats::quantile(HSF_hits_per_kb, 0.25)),
    Q3 = unname(stats::quantile(HSF_hits_per_kb, 0.75)),
    mean = mean(HSF_hits_per_kb), .groups = "drop"
  ) |>
  dplyr::transmute(
    `Relative PWM-score threshold` = threshold,
    `Region class` = as.character(class), n,
    `Median density (placements/kb)` = median,
    `Q1 density (placements/kb)` = Q1,
    `Q3 density (placements/kb)` = Q3,
    `Mean density (placements/kb)` = mean
  )

region_threshold <- region_metrics |>
  dplyr::transmute(
    `Region ID` = region_id, `Region class` = as.character(class),
    `Region length (bp)` = width_bp,
    `Relative PWM-score threshold` = threshold,
    `Non-redundant HSF motif-coordinate placements` = nonredundant_HSF_motif_coordinate_placements,
    `Density (placements/kb)` = HSF_hits_per_kb
  )
region_continuous <- continuous_metrics |>
  dplyr::transmute(
    `Region ID` = region_id, `Region class` = as.character(class),
    Chromosome = chromosome, Start = start, End = end, `Width (bp)` = width_bp,
    `Maximum HSF relative PWM score` = maximum_HSF_relative_score,
    `Mean top-five HSF relative PWM score` = mean_top5_HSF_relative_score,
    `Median per-model maximum HSF score` = median_per_model_maximum_score
  )
continuous_table <- continuous_stats |>
  dplyr::mutate(metric_column = metric)

write_sheet("S7A_threshold_stats", "Table S7A. Threshold sensitivity of non-redundant HSF motif-coordinate placement density.", threshold_table)
write_sheet("S7B_threshold_summary", "Table S7B. Distribution summaries for non-redundant HSF motif-coordinate placement density.", threshold_summary)
write_sheet("S7C_continuous_stats", "Table S7C. Threshold-independent continuous HSF relative-score comparisons using the revised region universe.", continuous_table)
write_sheet("S7D_region_threshold", "Table S7D. Region-level non-redundant HSF motif-coordinate placement metrics across PWM-score thresholds.", region_threshold)
write_sheet("S7E_region_continuous", "Table S7E. Region-level continuous HSF relative-score metrics for the revised comparison universe.", region_continuous)
write_sheet("S7F_QC", "Table S7F. Revised analysis definitions and quality-control information.", qc)
write_sheet("S7G_HSF_models", "Table S7G. Arabidopsis HSF PWM-model inventory.", motif_inventory)
write_sheet("S7H_original_QC", "Table S7H. Original continuous-score scan QC and reconciliation to the final background universe.", original_qc)

openxlsx::saveWorkbook(wb, TABLE_FILE, overwrite = TRUE)

###############################################################################
# 15. RECONSTRUCT FINAL FIGURE S2
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
  ggplot2::labs(x = "Relative PWM-score threshold", y = "Median non-redundant HSF motif-coordinate placements per kb") +
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

figure_s2 <- (panel_a / (panel_b | panel_c)) + patchwork::plot_annotation(tag_levels = "A")

# Mandatory for the user's RStudio workflow: display before saving.
print(figure_s2)

FIG_PDF <- file.path(FIG_DIR, "FigS2_HSF_threshold_and_continuous_scores.pdf")
FIG_PNG <- file.path(FIG_DIR, "FigS2_HSF_threshold_and_continuous_scores.png")

ggplot2::ggsave(FIG_PDF, figure_s2, width = 10.5, height = 10, units = "in")
ggplot2::ggsave(FIG_PNG, figure_s2, width = 10.5, height = 10, units = "in", dpi = 600, bg = "white")

###############################################################################
# 16. VALIDATION AND FINAL OUTPUT
###############################################################################

capture.output(sessionInfo(), file = file.path(OUT_DIR, "S7_sessionInfo.txt"))

if (!file.exists(TABLE_FILE) || !file.exists(FIG_PDF) || !file.exists(FIG_PNG)) stop("One or more final output files were not created.")
if (nrow(region_metrics[region_metrics$class == ONSEN_CLASS & region_metrics$threshold == 0.85, ]) != 16L) stop("ONSEN row validation failed.")
if (!all(sort(unique(region_metrics$threshold)) == THRESHOLDS)) stop("Threshold validation failed.")

cat("\n============================================================\n")
cat("FINAL HSF THRESHOLD AND CONTINUOUS-SCORE ANALYSIS COMPLETED\n")
cat("============================================================\n\n")

cat("Threshold statistics:\n"); print(threshold_stats, row.names = FALSE)
cat("\nContinuous-score statistics:\n"); print(continuous_stats, row.names = FALSE)

cat("\nOutputs:\n")
cat("Table S7 reconstructed:\n  ", TABLE_FILE, "\n", sep = "")
cat("Figure S2 PDF:\n  ", FIG_PDF, "\n", sep = "")
cat("Figure S2 PNG:\n  ", FIG_PNG, "\n", sep = "")
cat("Source-data directory:\n  ", OUT_DIR, "\n", sep = "")
