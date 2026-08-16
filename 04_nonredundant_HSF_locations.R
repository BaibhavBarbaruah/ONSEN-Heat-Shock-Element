# Non-redundant HSF-compatible sequence locations.
# Covers final Fig. S1 and Table S5.

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

h49_all <- normalize_49bp_hits(read_table_auto(hit49_file))
sequence_ids49 <- unique(h49_all$sequence_id)
h49 <- h49_all[
  (h49_all$tf_family == "HSF" | grepl("^HSF", h49_all$motif_name, ignore.case = TRUE)) &
    as.numeric(h49_all$relative_score) >= threshold,
  ,
  drop = FALSE
]
h49$forward_start <- as.integer(h49$forward_start)
h49$forward_end <- as.integer(h49$forward_end)

clustered49 <- merge_overlapping_intervals(h49)
locations49 <- summarize_interval_clusters(clustered49)

raw49 <- h49 |>
  dplyr::count(sequence_id, name = "HSF_PWM_model_position_predictions")
principal49 <- h49 |>
  dplyr::distinct(sequence_id, forward_start, forward_end)
principal49_counts <- principal49 |>
  dplyr::count(
    sequence_id,
    name = "principal_nonredundant_HSF_coordinate_placements"
  )
strict49_counts <- locations49 |>
  dplyr::count(sequence_id, name = "strict_overlap_merged_HSF_regions")
model49 <- h49 |>
  dplyr::group_by(sequence_id) |>
  dplyr::summarise(
    distinct_HSF_PWM_models = dplyr::n_distinct(motif_id),
    maximum_relative_PWM_score = max(as.numeric(relative_score), na.rm = TRUE),
    .groups = "drop"
  )
summary49 <- data.frame(sequence_id = sequence_ids49) |>
  dplyr::left_join(raw49, by = "sequence_id") |>
  dplyr::left_join(principal49_counts, by = "sequence_id") |>
  dplyr::left_join(strict49_counts, by = "sequence_id") |>
  dplyr::left_join(model49, by = "sequence_id") |>
  dplyr::mutate(
    HSF_PWM_model_position_predictions = tidyr::replace_na(
      HSF_PWM_model_position_predictions, 0L
    ),
    principal_nonredundant_HSF_coordinate_placements = tidyr::replace_na(
      principal_nonredundant_HSF_coordinate_placements, 0L
    ),
    strict_overlap_merged_HSF_regions = tidyr::replace_na(
      strict_overlap_merged_HSF_regions, 0L
    ),
    distinct_HSF_PWM_models = tidyr::replace_na(distinct_HSF_PWM_models, 0L),
    sequence_length_bp = 49L,
    principal_placement_density_per_kb =
      principal_nonredundant_HSF_coordinate_placements / sequence_length_bp * 1000,
    strict_overlap_density_per_kb =
      strict_overlap_merged_HSF_regions / sequence_length_bp * 1000,
    strict_overlap_reduction_percent = dplyr::if_else(
      principal_nonredundant_HSF_coordinate_placements > 0,
      100 * (1 - strict_overlap_merged_HSF_regions /
        principal_nonredundant_HSF_coordinate_placements),
      NA_real_
    ),
    analysis_space = "49-bp sequence"
  )

copy_hits <- normalize_copy_hits(read_table_auto(copy_hit_file))
copy_hits <- copy_hits[as.numeric(copy_hits$relative_score) >= threshold, , drop = FALSE]

windows <- read.csv(repo_file("ONSEN_Col0_terminal_candidate_windows.csv"))
source_sequence_id <- as.character(copy_hits$sequence_id)
needs_window_id <- !source_sequence_id %in% windows$window_id
if (any(needs_window_id)) {
  parsed_copy <- sub(
    ".*(ONSEN[1-8]).*", "\\1", source_sequence_id[needs_window_id],
    perl = TRUE
  )
  parsed_side <- ifelse(
    grepl("left|-L$", source_sequence_id[needs_window_id], ignore.case = TRUE),
    "-L", "-R"
  )
  copy_hits$sequence_id[needs_window_id] <- paste0(parsed_copy, parsed_side)
}
if (any(!copy_hits$sequence_id %in% windows$window_id)) {
  stop("Could not map all Col-0 hit identifiers to the sixteen terminal windows.", call. = FALSE)
}
copy_hits$source_sequence_id <- source_sequence_id
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
  dplyr::count(sequence_id, name = "HSF_PWM_model_position_predictions")
principal_copy <- copy_hits |>
  dplyr::distinct(sequence_id, forward_start, forward_end)
principal_copy_counts <- principal_copy |>
  dplyr::count(
    sequence_id,
    name = "principal_nonredundant_HSF_coordinate_placements"
  )
strict_copy_counts <- locations_copy |>
  dplyr::count(sequence_id, name = "strict_overlap_merged_HSF_regions")
model_copy <- copy_hits |>
  dplyr::group_by(sequence_id) |>
  dplyr::summarise(
    distinct_HSF_PWM_models = dplyr::n_distinct(motif_id),
    maximum_relative_PWM_score = max(as.numeric(relative_score), na.rm = TRUE),
    .groups = "drop"
  )
summary_copy <- windows |>
  dplyr::transmute(
    sequence_id = window_id,
    ONSEN_copy = copy_id,
    terminal_side,
    sequence_length_bp = width_bp
  ) |>
  dplyr::left_join(raw_copy, by = "sequence_id") |>
  dplyr::left_join(principal_copy_counts, by = "sequence_id") |>
  dplyr::left_join(strict_copy_counts, by = "sequence_id") |>
  dplyr::left_join(model_copy, by = "sequence_id") |>
  dplyr::mutate(
    HSF_PWM_model_position_predictions = tidyr::replace_na(
      HSF_PWM_model_position_predictions, 0L
    ),
    principal_nonredundant_HSF_coordinate_placements = tidyr::replace_na(
      principal_nonredundant_HSF_coordinate_placements, 0L
    ),
    strict_overlap_merged_HSF_regions = tidyr::replace_na(
      strict_overlap_merged_HSF_regions, 0L
    ),
    distinct_HSF_PWM_models = tidyr::replace_na(distinct_HSF_PWM_models, 0L),
    principal_placement_density_per_kb =
      principal_nonredundant_HSF_coordinate_placements / sequence_length_bp * 1000,
    strict_overlap_density_per_kb =
      strict_overlap_merged_HSF_regions / sequence_length_bp * 1000,
    strict_overlap_reduction_percent = dplyr::if_else(
      principal_nonredundant_HSF_coordinate_placements > 0,
      100 * (1 - strict_overlap_merged_HSF_regions /
        principal_nonredundant_HSF_coordinate_placements),
      NA_real_
    ),
    analysis_space = "Col-0 ONSEN 800-bp terminal candidate window"
  )

combined_summary <- dplyr::bind_rows(summary49, summary_copy)

safe_write_csv(h49, "S5_49bp_HSF_PWM_model_predictions_threshold_0p85.csv")
safe_write_csv(principal49, "S5_49bp_HSF_principal_exact_coordinates_threshold_0p85.csv")
safe_write_csv(clustered49, "S5_49bp_HSF_strict_overlap_clusters_threshold_0p85.csv")
safe_write_csv(locations49, "S5_49bp_HSF_strict_overlap_regions_threshold_0p85.csv")
safe_write_csv(summary49, "S5_49bp_HSF_principal_vs_strict_summary_threshold_0p85.csv")

safe_write_csv(copy_hits, "S5_Col0_ONSEN_HSF_PWM_model_predictions_threshold_0p85.csv")
safe_write_csv(principal_copy, "S5_Col0_ONSEN_HSF_principal_exact_coordinates_threshold_0p85.csv")
safe_write_csv(clustered_copy, "S5_Col0_ONSEN_HSF_strict_overlap_clusters_threshold_0p85.csv")
safe_write_csv(locations_copy, "S5_Col0_ONSEN_HSF_strict_overlap_regions_threshold_0p85.csv")
safe_write_csv(summary_copy, "S5_Col0_ONSEN_HSF_principal_vs_strict_summary_threshold_0p85.csv")
safe_write_csv(combined_summary, "S5_combined_HSF_principal_vs_strict_summary_threshold_0p85.csv")

if (ONSEN_MAKE_FIGURES) {
  figure_data <- summary_copy |>
    tidyr::pivot_longer(
      cols = c(
        "principal_nonredundant_HSF_coordinate_placements",
        "strict_overlap_merged_HSF_regions"
      ),
      names_to = "metric", values_to = "count"
    ) |>
    dplyr::mutate(
      metric = dplyr::recode(
        metric,
        principal_nonredundant_HSF_coordinate_placements =
          "Principal exact-coordinate placements",
        strict_overlap_merged_HSF_regions =
          "Strict overlap-merged regions"
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
      "Principal exact-coordinate placements" = 16,
      "Strict overlap-merged regions" = 17
    )) +
    ggplot2::labs(
      x = "HSF motif-region count",
      y = "Col-0 ONSEN terminal candidate window",
      shape = NULL
    ) +
    theme_onsen(13) +
    ggplot2::theme(legend.position = "top")
  save_plot_pair(p_s4, "FigS1_principal_vs_strict_overlap_HSF_regions", 7.2, 7.4)
}

message("Non-redundant HSF-location analysis completed.")
