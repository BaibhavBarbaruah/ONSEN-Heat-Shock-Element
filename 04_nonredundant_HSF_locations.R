# Non-redundant HSF-compatible sequence locations.
# Covers Fig. S4 and source data for Table S13.

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

clustered49 <- merge_overlapping_intervals(h49)
locations49 <- summarize_interval_clusters(clustered49)

raw49 <- h49 |>
  dplyr::count(sequence_id, name = "HSF_motif_model_position_hits")
nr49 <- locations49 |>
  dplyr::count(sequence_id, name = "nonredundant_HSF_locations")
summary49 <- dplyr::full_join(raw49, nr49, by = "sequence_id") |>
  dplyr::mutate(
    HSF_motif_model_position_hits = tidyr::replace_na(HSF_motif_model_position_hits, 0L),
    nonredundant_HSF_locations = tidyr::replace_na(nonredundant_HSF_locations, 0L),
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

clustered_copy <- merge_overlapping_intervals(copy_hits)
locations_copy <- summarize_interval_clusters(clustered_copy)

raw_copy <- copy_hits |>
  dplyr::count(sequence_id, name = "HSF_motif_model_position_hits")
nr_copy <- locations_copy |>
  dplyr::count(sequence_id, name = "nonredundant_HSF_locations")
summary_copy <- dplyr::full_join(raw_copy, nr_copy, by = "sequence_id") |>
  dplyr::mutate(
    HSF_motif_model_position_hits = tidyr::replace_na(HSF_motif_model_position_hits, 0L),
    nonredundant_HSF_locations = tidyr::replace_na(nonredundant_HSF_locations, 0L),
    proportion_collapsed = 1 -
      nonredundant_HSF_locations / HSF_motif_model_position_hits,
    analysis_space = "Col-0 ONSEN 800-bp terminal candidate window"
  )

combined_summary <- dplyr::bind_rows(summary49, summary_copy)

safe_write_csv(h49, "Revision_R1_2_49bp_HSF_raw_hits_threshold_0p85_repository.csv")
safe_write_csv(clustered49, "Revision_R1_2_49bp_HSF_clustered_hits_threshold_0p85_repository.csv")
safe_write_csv(locations49, "Revision_R1_2_49bp_HSF_nonredundant_locations_threshold_0p85_repository.csv")
safe_write_csv(summary49, "Revision_R1_2_49bp_HSF_raw_vs_nonredundant_summary_threshold_0p85_repository.csv")

safe_write_csv(copy_hits, "Revision_R1_2_Col0_ONSEN_HSF_raw_hits_forward_coordinates_threshold_0p85_repository.csv")
safe_write_csv(clustered_copy, "Revision_R1_2_Col0_ONSEN_HSF_clustered_hits_threshold_0p85_repository.csv")
safe_write_csv(locations_copy, "Revision_R1_2_Col0_ONSEN_HSF_nonredundant_locations_threshold_0p85_repository.csv")
safe_write_csv(summary_copy, "Revision_R1_2_Col0_ONSEN_HSF_raw_vs_nonredundant_summary_threshold_0p85_repository.csv")
safe_write_csv(combined_summary, "Revision_R1_2_combined_HSF_raw_vs_nonredundant_summary_threshold_0p85_repository.csv")

if (ONSEN_MAKE_FIGURES) {
  figure_data <- summary_copy |>
    tidyr::pivot_longer(
      cols = c("HSF_motif_model_position_hits", "nonredundant_HSF_locations"),
      names_to = "metric", values_to = "count"
    ) |>
    dplyr::mutate(
      metric = dplyr::recode(
        metric,
        HSF_motif_model_position_hits = "PWM model-position hits",
        nonredundant_HSF_locations = "Non-redundant HSF locations"
      ),
      sequence_id = factor(sequence_id, levels = rev(windows$window_id))
    )

  p_s4 <- ggplot2::ggplot(
    figure_data,
    ggplot2::aes(count, sequence_id, group = sequence_id)
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(shape = metric), size = 3.3) +
    ggplot2::scale_shape_manual(values = c(
      "PWM model-position hits" = 16,
      "Non-redundant HSF locations" = 17
    )) +
    ggplot2::labs(
      x = "HSF motif count",
      y = "Col-0 ONSEN terminal candidate window",
      shape = NULL
    ) +
    theme_onsen(13) +
    ggplot2::theme(legend.position = "top")
  save_plot_pair(p_s4, "FigS4_HSF_model_hits_vs_nonredundant_locations", 7.2, 7.4)
}

message("Non-redundant HSF-location analysis completed.")
