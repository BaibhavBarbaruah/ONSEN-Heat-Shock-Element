# Non-redundant HSF-compatible sequence locations.
# Covers Fig. S1 and source data for final Table S5.

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "tidyr", "readr", "ggplot2"))

message_config()
threshold <- 0.85

normalize_49bp_hits <- function(x) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("sequence_id", c("sequence_id", "sequence_type", "sequence"))
  rename_first("motif_id", c("motif_id", "tf_id", "jaspar_id"))
  rename_first("motif_name", c("motif_name", "tf_name", "model_name"))
  rename_first("tf_family", c("tf_family", "family"))
  rename_first("strand", c("strand"))
  rename_first("forward_start", c("forward_start", "motif_start", "start", "position"))
  rename_first("forward_end", c("forward_end", "motif_end", "end"))
  rename_first("motif_width", c("motif_width", "motif_length", "width"))
  rename_first("matched_sequence", c("matched_sequence", "match_sequence"))
  rename_first("relative_score", c("relative_score", "score"))
  assert_columns(x, c("sequence_id", "motif_id", "motif_name", "relative_score"), "49-bp hit table")
  if (!"tf_family" %in% names(x)) x$tf_family <- classify_tf_family(x$motif_name)
  if (!"motif_width" %in% names(x) && "matched_sequence" %in% names(x)) {
    x$motif_width <- nchar(x$matched_sequence)
  }
  if (!"forward_end" %in% names(x) && all(c("forward_start", "motif_width") %in% names(x))) {
    x$forward_end <- as.integer(x$forward_start) + as.integer(x$motif_width) - 1L
  }
  x
}

normalize_copy_hits <- function(x) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("sequence_id", c("sequence_id", "candidate_id"))
  rename_first("motif_id", c("motif_id", "tf_id", "jaspar_id"))
  rename_first("motif_name", c("motif_name", "tf_name", "model_name"))
  rename_first("strand", c("strand"))
  rename_first("scan_start", c("scan_start", "position", "start", "motif_start"))
  rename_first("motif_width", c("motif_width", "motif_length", "width"))
  rename_first("matched_sequence", c("matched_sequence", "match_sequence"))
  rename_first("relative_score", c("relative_score", "score"))
  assert_columns(
    x,
    c("sequence_id", "motif_id", "motif_name", "strand", "scan_start", "relative_score"),
    "copy-wide HSF hit table"
  )
  if (!"motif_width" %in% names(x) && "matched_sequence" %in% names(x)) {
    x$motif_width <- nchar(x$matched_sequence)
  }
  x
}

hit49_file <- find_any_input(c(
  "native_vs_mutated_49bp_TF_motif_hits_annotated.csv",
  "native_vs_mutated_49bp_TF_motif_hits_repository.csv",
  "native_vs_mutated_49bp_TF_motif_hits.csv"
))
copy_hit_file <- find_any_input(c(
  "Col0_ONSEN_LTRcandidate_JASPAR2026_Arabidopsis_HSF_hits.csv",
  "Col0_ONSEN_HSF_hits_repository.csv",
  "Col0_ONSEN_LTRcandidate_high_confidence_motif_hits_COMBINED.csv"
))

h49 <- normalize_49bp_hits(read_table_auto(hit49_file))
h49 <- h49[
  (h49$tf_family == "HSF" | grepl("^HSF", h49$motif_name, ignore.case = TRUE)) &
    as.numeric(h49$relative_score) >= threshold,
  ,
  drop = FALSE
]
h49$forward_start <- as.integer(h49$forward_start)
h49$forward_end <- as.integer(h49$forward_end)

exact49 <- h49 |>
  dplyr::distinct(sequence_id, forward_start, forward_end, .keep_all = TRUE)
clustered49 <- merge_overlapping_intervals(h49)
locations49 <- summarize_interval_clusters(clustered49)

raw49 <- h49 |>
  dplyr::count(sequence_id, name = "HSF_PWM_model_position_hits")
exact49_summary <- exact49 |>
  dplyr::count(sequence_id, name = "exact_coordinate_HSF_placements")
overlap49_summary <- locations49 |>
  dplyr::count(sequence_id, name = "strict_overlap_merged_HSF_regions")
summary49 <- raw49 |>
  dplyr::full_join(exact49_summary, by = "sequence_id") |>
  dplyr::full_join(overlap49_summary, by = "sequence_id") |>
  dplyr::mutate(
    HSF_PWM_model_position_hits = tidyr::replace_na(HSF_PWM_model_position_hits, 0L),
    exact_coordinate_HSF_placements = tidyr::replace_na(exact_coordinate_HSF_placements, 0L),
    strict_overlap_merged_HSF_regions = tidyr::replace_na(strict_overlap_merged_HSF_regions, 0L),
    analysis_space = "49-bp sequence"
  )

copy_hits <- normalize_copy_hits(read_table_auto(copy_hit_file))
copy_hits <- copy_hits[as.numeric(copy_hits$relative_score) >= threshold, , drop = FALSE]

windows <- read.csv(repo_file("ONSEN_Col0_terminal_candidate_windows.csv"))
window_widths <- setNames(windows$width_bp, windows$window_id)
copy_hits$sequence_width <- window_widths[copy_hits$sequence_id]
copy_hits$sequence_width[is.na(copy_hits$sequence_width)] <- 800L

# Convert reverse-strand scan positions to forward-sequence coordinates:
# forward_start = L - scan_start - motif_width + 2
copy_hits$forward_start <- ifelse(
  copy_hits$strand == "+",
  as.integer(copy_hits$scan_start),
  as.integer(copy_hits$sequence_width) -
    as.integer(copy_hits$scan_start) -
    as.integer(copy_hits$motif_width) + 2L
)
copy_hits$forward_end <- copy_hits$forward_start + as.integer(copy_hits$motif_width) - 1L

exact_copy <- copy_hits |>
  dplyr::distinct(sequence_id, forward_start, forward_end, .keep_all = TRUE)
clustered_copy <- merge_overlapping_intervals(copy_hits)
locations_copy <- summarize_interval_clusters(clustered_copy)

raw_copy <- copy_hits |>
  dplyr::count(sequence_id, name = "HSF_PWM_model_position_hits")
exact_copy_summary <- exact_copy |>
  dplyr::count(sequence_id, name = "exact_coordinate_HSF_placements")
overlap_copy_summary <- locations_copy |>
  dplyr::count(sequence_id, name = "strict_overlap_merged_HSF_regions")
summary_copy <- raw_copy |>
  dplyr::full_join(exact_copy_summary, by = "sequence_id") |>
  dplyr::full_join(overlap_copy_summary, by = "sequence_id") |>
  dplyr::mutate(
    HSF_PWM_model_position_hits = tidyr::replace_na(HSF_PWM_model_position_hits, 0L),
    exact_coordinate_HSF_placements = tidyr::replace_na(exact_coordinate_HSF_placements, 0L),
    strict_overlap_merged_HSF_regions = tidyr::replace_na(strict_overlap_merged_HSF_regions, 0L),
    strict_overlap_fraction_of_exact = ifelse(
      exact_coordinate_HSF_placements > 0,
      strict_overlap_merged_HSF_regions / exact_coordinate_HSF_placements,
      NA_real_
    ),
    analysis_space = "Col-0 ONSEN 800-bp terminal candidate window"
  )

combined_summary <- dplyr::bind_rows(summary49, summary_copy)

safe_write_csv(h49, "Revision_R1_2_49bp_HSF_raw_hits_threshold_0p85_repository.csv")
safe_write_csv(exact49, "Revision_R1_2_49bp_HSF_exact_coordinate_placements_threshold_0p85_repository.csv")
safe_write_csv(clustered49, "Revision_R1_2_49bp_HSF_strict_overlap_clusters_threshold_0p85_repository.csv")
safe_write_csv(locations49, "Revision_R1_2_49bp_HSF_strict_overlap_regions_threshold_0p85_repository.csv")
safe_write_csv(summary49, "Revision_R1_2_49bp_HSF_exact_vs_overlap_summary_threshold_0p85_repository.csv")

safe_write_csv(copy_hits, "Revision_R1_2_Col0_ONSEN_HSF_raw_hits_forward_coordinates_threshold_0p85_repository.csv")
safe_write_csv(exact_copy, "Revision_R1_2_Col0_ONSEN_HSF_exact_coordinate_placements_threshold_0p85_repository.csv")
safe_write_csv(clustered_copy, "Revision_R1_2_Col0_ONSEN_HSF_strict_overlap_clusters_threshold_0p85_repository.csv")
safe_write_csv(locations_copy, "Revision_R1_2_Col0_ONSEN_HSF_strict_overlap_regions_threshold_0p85_repository.csv")
safe_write_csv(summary_copy, "Revision_R1_2_Col0_ONSEN_HSF_exact_vs_overlap_summary_threshold_0p85_repository.csv")
safe_write_csv(combined_summary, "Revision_R1_2_combined_HSF_exact_vs_overlap_summary_threshold_0p85_repository.csv")

if (ONSEN_MAKE_FIGURES) {
  figure_data <- summary_copy |>
    tidyr::pivot_longer(
      cols = c("exact_coordinate_HSF_placements", "strict_overlap_merged_HSF_regions"),
      names_to = "metric", values_to = "count"
    ) |>
    dplyr::mutate(
      metric = dplyr::recode(
        metric,
        exact_coordinate_HSF_placements = "Exact-coordinate HSF placements",
        strict_overlap_merged_HSF_regions = "Strict overlap-merged HSF regions"
      ),
      sequence_id = factor(sequence_id, levels = rev(windows$window_id))
    )

  p_s1 <- ggplot2::ggplot(
    figure_data,
    ggplot2::aes(count, sequence_id, group = sequence_id)
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(shape = metric), size = 3.3) +
    ggplot2::scale_shape_manual(values = c(
      "Exact-coordinate HSF placements" = 16,
      "Strict overlap-merged HSF regions" = 17
    )) +
    ggplot2::labs(
      x = "HSF motif count",
      y = "Col-0 ONSEN terminal candidate window",
      shape = NULL
    ) +
    theme_onsen(13) +
    ggplot2::theme(legend.position = "top")
  save_plot_pair(p_s1, "FigS1_exact_coordinate_vs_strict_overlap_HSF", 7.2, 7.4)
}

message("Non-redundant HSF-location analysis completed.")
