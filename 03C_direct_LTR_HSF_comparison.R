# Direct element-level HSF-density comparison for Fig. 4C.
# Principal units: eight ONSEN elements and 794 annotation-defined non-ONSEN
# LTR-retrotransposon elements. This is distinct from the broader 1,930-region
# strict TE-background sensitivity analysis in Fig. 3 and Tables S6-S7.

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "ggplot2"))

message_config()

input_file <- find_any_input(c(
  "Fig4C_direct_ONSEN_nonONSEN_LTR_element_HSF_density.tsv",
  "Fig4C_direct_ONSEN_nonONSEN_LTR_element_HSF_density.csv",
  "Revision_Fig4C_direct_LTR_element_HSF_density.csv",
  "direct_ONSEN_vs_nonONSEN_LTR_HSF_element_density.csv"
), required = FALSE)

rename_first <- function(x, target, alternatives) {
  hit <- alternatives[alternatives %in% names(x)][1]
  if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <- target
  x
}

standardize_class <- function(x) {
  result <- as.character(x)
  result[grepl("non[- _]?ONSEN|comparator", result, ignore.case = TRUE)] <-
    "non-ONSEN LTR-retrotransposon"
  result[grepl("ONSEN", result, ignore.case = TRUE) &
           !grepl("non[- _]?ONSEN", result, ignore.case = TRUE)] <- "ONSEN"
  result
}

if (is.na(input_file)) {
  summary_file <- repo_file(file.path(
    "source_data", "Figure4C_direct_LTR_HSF_summary.tsv"
  ))
  summary_only <- read.delim(summary_file, check.names = FALSE, stringsAsFactors = FALSE)
  safe_write_csv(summary_only, "Fig4C_direct_LTR_HSF_summary_repository.csv")
  warning(
    "Element-level Fig. 4C scan input was not found. The deposited summary was ",
    "validated and copied, but element-level statistics/plotting were skipped. ",
    "See INPUT_PROVENANCE.tsv for the required processed input."
  )
} else {
  element_data <- read_table_auto(input_file)
  element_data <- rename_first(
    element_data, "element_id",
    c("element_id", "TE_id", "locus_id", "feature_id", "candidate_id")
  )
  element_data <- rename_first(
    element_data, "element_class",
    c("element_class", "locus_class", "region_class", "class", "group")
  )
  element_data <- rename_first(
    element_data, "HSF_density_placements_per_kb",
    c(
      "HSF_density_placements_per_kb", "HSF_placements_per_kb",
      "HSF_hits_per_kb", "hsf_density", "density_placements_per_kb"
    )
  )
  element_data <- rename_first(
    element_data, "HSF_placements",
    c("HSF_placements", "nonredundant_HSF_placements", "HSF_hits")
  )
  element_data <- rename_first(
    element_data, "sequence_length_bp",
    c("sequence_length_bp", "length_bp", "width_bp", "region_length_bp")
  )
  assert_columns(element_data, c("element_id", "element_class"), "Fig. 4C input")
  element_data$element_class <- standardize_class(element_data$element_class)
  element_data <- element_data[element_data$element_class %in% c(
    "ONSEN", "non-ONSEN LTR-retrotransposon"
  ), , drop = FALSE]

  duplicate_units <- anyDuplicated(element_data$element_id) > 0L
  if (duplicate_units) {
    assert_columns(
      element_data,
      c("HSF_placements", "sequence_length_bp"),
      "multi-row Fig. 4C element input"
    )
    element_data <- element_data |>
      dplyr::group_by(element_id, element_class) |>
      dplyr::summarise(
        HSF_placements = sum(as.numeric(HSF_placements), na.rm = TRUE),
        sequence_length_bp = sum(as.numeric(sequence_length_bp), na.rm = TRUE),
        HSF_density_placements_per_kb =
          1000 * HSF_placements / sequence_length_bp,
        .groups = "drop"
      )
  } else if (!"HSF_density_placements_per_kb" %in% names(element_data)) {
    assert_columns(
      element_data,
      c("HSF_placements", "sequence_length_bp"),
      "Fig. 4C element input"
    )
    element_data$HSF_density_placements_per_kb <-
      1000 * as.numeric(element_data$HSF_placements) /
      as.numeric(element_data$sequence_length_bp)
  }

  element_data$HSF_density_placements_per_kb <-
    as.numeric(element_data$HSF_density_placements_per_kb)
  element_data <- element_data[is.finite(
    element_data$HSF_density_placements_per_kb
  ), , drop = FALSE]

  n_by_class <- table(element_data$element_class)
  if (!identical(unname(n_by_class["ONSEN"]), 8L) ||
      !identical(unname(n_by_class["non-ONSEN LTR-retrotransposon"]), 794L)) {
    stop(
      "Fig. 4C input must contain exactly 8 ONSEN and 794 non-ONSEN ",
      "LTR-retrotransposon elements after filtering. Observed: ",
      paste(names(n_by_class), n_by_class, collapse = "; "),
      call. = FALSE
    )
  }

  onsen <- element_data$HSF_density_placements_per_kb[
    element_data$element_class == "ONSEN"
  ]
  comparator <- element_data$HSF_density_placements_per_kb[
    element_data$element_class == "non-ONSEN LTR-retrotransposon"
  ]
  test <- stats::wilcox.test(onsen, comparator, exact = FALSE)
  stats_out <- data.frame(
    comparison = "ONSEN versus annotation-defined non-ONSEN LTR-retrotransposons",
    metric = "non-redundant HSF motif-coordinate placement density (placements/kb)",
    ONSEN_n = length(onsen),
    ONSEN_median = stats::median(onsen),
    comparator_n = length(comparator),
    comparator_median = stats::median(comparator),
    Wilcoxon_P = test$p.value,
    Cliffs_delta = cliffs_delta(onsen, comparator),
    stringsAsFactors = FALSE
  )

  safe_write_csv(element_data, "Fig4C_direct_LTR_HSF_element_data_repository.csv")
  safe_write_csv(stats_out, "Fig4C_direct_LTR_HSF_statistics_repository.csv")

  if (ONSEN_MAKE_FIGURES) {
    element_data$element_class <- factor(
      element_data$element_class,
      levels = c("non-ONSEN LTR-retrotransposon", "ONSEN")
    )
    p4c <- ggplot2::ggplot(
      element_data,
      ggplot2::aes(element_class, HSF_density_placements_per_kb, fill = element_class)
    ) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.55, colour = "grey35") +
      ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85) +
      ggplot2::scale_fill_manual(values = c(
        "non-ONSEN LTR-retrotransposon" = "#83AEE8",
        "ONSEN" = "#E989AE"
      )) +
      ggplot2::labs(
        x = NULL,
        y = "Non-redundant HSF placements per kb",
        fill = NULL
      ) +
      theme_onsen(13) +
      ggplot2::theme(legend.position = "none")
    save_plot_pair(p4c, "Fig4C_direct_LTR_HSF_density", 6.4, 5.2)
  }
}

message("Direct ONSEN/non-ONSEN LTR HSF-density analysis completed.")
