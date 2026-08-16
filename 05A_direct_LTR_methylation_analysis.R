# Principal direct element/copy-level methylation analysis for Fig. 4A and
# Table S8A. The public methylome is basal/unstressed Col-0 leaf material.

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "ggplot2"))

message_config()

input_file <- find_any_input(c(
  "Fig4A_direct_ONSEN_nonONSEN_LTR_element_methylation.tsv",
  "Fig4A_direct_ONSEN_nonONSEN_LTR_element_methylation.csv",
  "Revision_Fig4_direct_LTR_element_methylation.csv",
  "direct_ONSEN_vs_nonONSEN_LTR_element_methylation.csv"
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

effect_magnitude <- function(delta) {
  absolute <- abs(delta)
  ifelse(absolute < 0.147, "negligible",
    ifelse(absolute < 0.33, "small",
      ifelse(absolute < 0.474, "medium", "large")))
}

significance_label <- function(adjusted_p) {
  ifelse(adjusted_p < 1e-4, "****",
    ifelse(adjusted_p < 1e-3, "***",
      ifelse(adjusted_p < 1e-2, "**",
        ifelse(adjusted_p < 0.05, "*", "ns"))))
}

if (is.na(input_file)) {
  summary_file <- repo_file(file.path(
    "source_data", "Figure4A_direct_LTR_methylation_summary.tsv"
  ))
  summary_only <- read.delim(summary_file, check.names = FALSE, stringsAsFactors = FALSE)
  safe_write_csv(summary_only, "Fig4A_direct_LTR_methylation_summary_repository.csv")
  warning(
    "Element-level Fig. 4A methylation input was not found. The deposited ",
    "Table S8A summary was validated and copied, but element-level statistics/",
    "plotting were skipped. See INPUT_PROVENANCE.tsv."
  )
} else {
  methylation <- read_table_auto(input_file)
  methylation <- rename_first(
    methylation, "element_id",
    c("element_id", "ONSEN_copy", "TE_id", "locus_id", "feature_id")
  )
  methylation <- rename_first(
    methylation, "element_class",
    c("element_class", "locus_class", "region_class", "class", "group")
  )
  methylation <- rename_first(
    methylation, "context",
    c("context", "methylation_context", "mc_context")
  )
  methylation <- rename_first(
    methylation, "weighted_methylation_percent",
    c(
      "weighted_methylation_percent", "weighted_methylation",
      "weighted_methylation_pct", "methylation_percent"
    )
  )
  methylation <- rename_first(
    methylation, "methylated_bases_sum",
    c("methylated_bases_sum", "methylated_bases", "methylated_count")
  )
  methylation <- rename_first(
    methylation, "total_bases_sum",
    c("total_bases_sum", "total_bases", "coverage_sum", "total_count")
  )
  assert_columns(
    methylation,
    c("element_id", "element_class", "context"),
    "Fig. 4A element methylation input"
  )
  methylation$element_class <- standardize_class(methylation$element_class)
  methylation$context <- toupper(as.character(methylation$context))
  methylation <- methylation[
    methylation$element_class %in% c("ONSEN", "non-ONSEN LTR-retrotransposon") &
      methylation$context %in% c("CG", "CHG", "CHH"),
    , drop = FALSE
  ]

  if (all(c("methylated_bases_sum", "total_bases_sum") %in% names(methylation))) {
    methylation <- methylation |>
      dplyr::group_by(element_id, element_class, context) |>
      dplyr::summarise(
        methylated_bases_sum = sum(as.numeric(methylated_bases_sum), na.rm = TRUE),
        total_bases_sum = sum(as.numeric(total_bases_sum), na.rm = TRUE),
        weighted_methylation_percent = weighted_methylation(
          methylated_bases_sum, total_bases_sum
        ),
        .groups = "drop"
      )
  } else {
    assert_columns(
      methylation,
      "weighted_methylation_percent",
      "Fig. 4A element methylation input"
    )
    if (anyDuplicated(methylation[c("element_id", "context")])) {
      stop(
        "Multiple terminal rows per element/context require methylated and total ",
        "base counts for weighted aggregation.",
        call. = FALSE
      )
    }
    methylation$weighted_methylation_percent <-
      as.numeric(methylation$weighted_methylation_percent)
  }
  methylation <- methylation[is.finite(
    methylation$weighted_methylation_percent
  ), , drop = FALSE]

  expected_onsen <- c(CG = 6L, CHG = 6L, CHH = 8L)
  contexts <- c("CG", "CHG", "CHH")
  stats_out <- do.call(rbind, lapply(contexts, function(context_name) {
    subset <- methylation[methylation$context == context_name, , drop = FALSE]
    onsen <- subset$weighted_methylation_percent[
      subset$element_class == "ONSEN"
    ]
    comparator <- subset$weighted_methylation_percent[
      subset$element_class == "non-ONSEN LTR-retrotransposon"
    ]
    if (length(onsen) != expected_onsen[[context_name]] || length(comparator) != 779L) {
      stop(
        "Unexpected Fig. 4A sample size for ", context_name,
        ": ONSEN=", length(onsen), ", comparator=", length(comparator),
        call. = FALSE
      )
    }
    test <- stats::wilcox.test(onsen, comparator, exact = FALSE)
    data.frame(
      context = context_name,
      ONSEN_unit = "ONSEN copy",
      ONSEN_n = length(onsen),
      ONSEN_median_weighted_methylation_percent = stats::median(onsen),
      comparator_unit = "non-ONSEN LTR-retrotransposon element",
      comparator_n = length(comparator),
      comparator_median_weighted_methylation_percent = stats::median(comparator),
      Wilcoxon_P = test$p.value,
      Cliffs_delta = cliffs_delta(onsen, comparator),
      stringsAsFactors = FALSE
    )
  }))
  stats_out$BH_adjusted_P <- p.adjust(stats_out$Wilcoxon_P, method = "BH")
  stats_out$effect_magnitude <- effect_magnitude(stats_out$Cliffs_delta)
  stats_out$significance <- significance_label(stats_out$BH_adjusted_P)

  safe_write_csv(methylation, "Fig4A_direct_LTR_methylation_element_data_repository.csv")
  safe_write_csv(stats_out, "Fig4A_direct_LTR_methylation_statistics_repository.csv")

  if (ONSEN_MAKE_FIGURES) {
    methylation$context <- factor(methylation$context, levels = contexts)
    methylation$element_class <- factor(
      methylation$element_class,
      levels = c("non-ONSEN LTR-retrotransposon", "ONSEN")
    )
    p4a <- ggplot2::ggplot(
      methylation,
      ggplot2::aes(context, weighted_methylation_percent, fill = element_class)
    ) +
      ggplot2::geom_boxplot(
        position = ggplot2::position_dodge(width = 0.75),
        width = 0.62, outlier.shape = NA, linewidth = 0.6
      ) +
      ggplot2::geom_point(
        ggplot2::aes(group = element_class),
        position = ggplot2::position_jitterdodge(
          jitter.width = 0.15, dodge.width = 0.75
        ),
        size = 1.15, alpha = 0.35
      ) +
      ggplot2::scale_fill_manual(values = c(
        "non-ONSEN LTR-retrotransposon" = "#83AEE8",
        "ONSEN" = "#E989AE"
      )) +
      ggplot2::labs(
        x = "Methylation context",
        y = "Weighted methylation (%)",
        fill = NULL
      ) +
      theme_onsen(13) +
      ggplot2::theme(legend.position = "top")
    save_plot_pair(p4a, "Fig4A_direct_LTR_methylation", 6.8, 5.4)
  }
}

message("Direct ONSEN/non-ONSEN LTR methylation analysis completed.")
