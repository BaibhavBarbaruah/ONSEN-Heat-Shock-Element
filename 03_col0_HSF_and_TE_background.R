# Col-0 ONSEN HSF terminal-window and strict non-ONSEN TE comparison.
# Covers Fig. 3 and final Tables S4 and S6.
#
# The principal metric is the number of unique physical 1-based start/end HSF
# motif-coordinate placements after collapsing duplicate PWM-model/strand hits.
# The full genome-wide exact-coordinate rescan is implemented in
# 03B_threshold_and_continuous_sensitivity.R. This script assembles/validates
# the Fig. 3 quantities from the final deposited source tables and never falls
# back to raw PWM-model-position hit counts.

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "tidyr", "ggplot2"))
message_config()

read_source <- function(filename, skip) {
  data.table::fread(
    repo_file(file.path("supplementary_table_source", filename)),
    skip = skip, data.table = FALSE, check.names = FALSE
  )
}

s4 <- read_source("Table_S4__S4_Col0_HSF.tsv", skip = 1L)
s6 <- read_source("Table_S6__S6_TE_background.tsv", skip = 3L)
s7_regions <- read_source("Table_S7__S7D_region_threshold.tsv", skip = 2L)

expected_s4 <- c(
  "ONSEN copy", "Terminal candidate window",
  "Non-redundant HSF motif-coordinate placements", "Distinct HSF PWM models",
  "Maximum relative score", "Terminal-window width (bp)",
  "Non-redundant HSF motif-coordinate density (placements/kb)"
)
assert_columns(s4, expected_s4, "Table S4 source")
if (nrow(s4) != 16L) stop("Table S4 must contain sixteen terminal windows.")
if (!all(s4[["Terminal-window width (bp)"]] == 800L)) stop("All final ONSEN terminal windows must be 800 bp.")
if (min(s4[["Non-redundant HSF motif-coordinate placements"]]) != 33L ||
    max(s4[["Non-redundant HSF motif-coordinate placements"]]) != 49L) {
  stop("Table S4 exact-coordinate placement range must be 33-49.")
}

assert_columns(
  s6,
  c(
    "Relative PWM-score threshold", "n ONSEN terminal windows",
    "n strict non-ONSEN TE background regions",
    "ONSEN median non-redundant HSF placement density (placements/kb)",
    "Background median non-redundant HSF placement density (placements/kb)",
    "Wilcoxon P", "BH-adjusted P", "Cliff's delta"
  ),
  "Table S6 source"
)
key85 <- s6[s6[["Relative PWM-score threshold"]] == 0.85, , drop = FALSE]
if (nrow(key85) != 1L ||
    key85[["n ONSEN terminal windows"]] != 16L ||
    key85[["n strict non-ONSEN TE background regions"]] != 1930L ||
    abs(key85[["ONSEN median non-redundant HSF placement density (placements/kb)"]] - 49.375) > 1e-12 ||
    abs(key85[["Background median non-redundant HSF placement density (placements/kb)"]] - 6.81238953796413) > 1e-12) {
  stop("Table S6 does not match the final 0.85 Fig. 3 comparison.")
}

assert_columns(
  s7_regions,
  c(
    "Region ID", "Region class", "Region length (bp)",
    "Relative PWM-score threshold",
    "Non-redundant HSF motif-coordinate placements", "Density (placements/kb)"
  ),
  "Table S7 region source"
)
fig3c_data <- s7_regions[s7_regions[["Relative PWM-score threshold"]] == 0.85, , drop = FALSE]
if (sum(fig3c_data[["Region class"]] == "ONSEN terminal candidate windows") != 16L ||
    sum(fig3c_data[["Region class"]] == "Strict non-ONSEN TE background") != 1930L) {
  stop("Fig. 3C source must contain 16 ONSEN and 1,930 background regions.")
}

safe_write_csv(s4, "Table_S4_Col0_HSF_repository.csv")
safe_write_csv(s6, "Table_S6_TE_background_repository.csv")
safe_write_csv(fig3c_data, "Fig3C_exact_coordinate_density_repository.csv")

# When large external inputs are available, invoke the exact-coordinate rescan
# rather than any raw-hit approximation.
if (ONSEN_RUN_LARGE_STEPS || ONSEN_FORCE_RESCAN) {
  source(repo_file("03B_threshold_and_continuous_sensitivity.R"), echo = FALSE)
}

if (ONSEN_MAKE_FIGURES) {
  s4$terminal_side <- ifelse(grepl("^Left", s4[["Terminal candidate window"]]), "Left terminal window", "Right terminal window")
  s4$copy_order <- factor(s4[["ONSEN copy"]], levels = paste0("ONSEN", 1:8))

  p3a <- ggplot2::ggplot(
    s4,
    ggplot2::aes(copy_order, .data[["Non-redundant HSF motif-coordinate placements"]], fill = terminal_side)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.68, colour = "black") +
    ggplot2::labs(x = "ONSEN copy", y = "Non-redundant HSF\nmotif-coordinate placements", fill = NULL) +
    theme_onsen(12) +
    ggplot2::theme(legend.position = "top", axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
  save_plot_pair(p3a, "Fig3A_terminal_window_HSF_placements", 7.2, 5.0)

  paired <- s4 |>
    dplyr::select(copy_order, terminal_side, count = .data[["Non-redundant HSF motif-coordinate placements"]]) |>
    tidyr::pivot_wider(names_from = terminal_side, values_from = count)
  paired_long <- paired |>
    tidyr::pivot_longer(cols = c("Left terminal window", "Right terminal window"), names_to = "terminal_side", values_to = "count")
  p3b <- ggplot2::ggplot(paired_long, ggplot2::aes(terminal_side, count, group = copy_order)) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::labs(x = "Terminal candidate window", y = "Non-redundant HSF\nmotif-coordinate placements") +
    theme_onsen(12)
  save_plot_pair(p3b, "Fig3B_paired_terminal_windows", 6.2, 5.2)

  p3c <- ggplot2::ggplot(
    fig3c_data,
    ggplot2::aes(.data[["Region class"]], .data[["Density (placements/kb)"]], fill = .data[["Region class"]])
  ) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.65) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA) +
    ggplot2::labs(x = NULL, y = "Non-redundant HSF motif-coordinate\nplacements per kb", fill = NULL) +
    theme_onsen(12) +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 15, hjust = 1))
  save_plot_pair(p3c, "Fig3C_ONSEN_vs_strict_TE_background", 6.4, 5.2)

  # Fig. 3D requires the final annotated outlier table listed in
  # INPUT_PROVENANCE.tsv. It is intentionally not guessed from family names.
  outlier_file <- find_input(
    "strict_TE_background_annotated_with_HSF_outlier_classes.csv",
    required = FALSE
  )
  if (!is.na(outlier_file)) {
    outliers <- read_table_auto(outlier_file)
    family_col <- c("TE_family", "te_family", "family", "Alias", "alias")[
      c("TE_family", "te_family", "family", "Alias", "alias") %in% names(outliers)
    ][1]
    class_col <- names(outliers)[grepl("outlier", names(outliers), ignore.case = TRUE)][1]
    if (!is.na(family_col) && !is.na(class_col)) {
      z <- outliers[as.logical(outliers[[class_col]]) %in% TRUE, , drop = FALSE] |>
        dplyr::count(.data[[family_col]], sort = TRUE, name = "n_outlier_regions")
      p3d <- ggplot2::ggplot(z, ggplot2::aes(n_outlier_regions, reorder(.data[[family_col]], n_outlier_regions))) +
        ggplot2::geom_col() +
        ggplot2::labs(x = "Number of HSF-rich TE outlier regions", y = "TE family / alias") +
        theme_onsen(11)
      save_plot_pair(p3d, "Fig3D_HSF_rich_TE_outlier_families", 6.4, 5.2)
    }
  }
}

message("Fig. 3 / Tables S4 and S6 exact-coordinate source validation completed.")
