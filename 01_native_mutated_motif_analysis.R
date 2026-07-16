# Native versus in-silico-designed 49-bp ONSEN motif analysis.
# Covers Fig. 1, Fig. 2 and source data for Tables S1-S3.

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "tidyr", "readr", "stringr",
                   "ggplot2", "forcats", "scales", "patchwork", "Biostrings"))

message_config()

sequences <- read_fasta_simple(repo_file("ONSEN_49bp_sequences.fasta"))
sequence_table <- data.frame(
  sequence_id = c("Native HSE window", "HSE-mutated window"),
  sequence = unname(sequences[c("native_49bp_ONSEN_HSE", "in_silico_designed_HSE_disrupted_49bp")]),
  stringsAsFactors = FALSE
)
sequence_table$GC_percent <- vapply(sequence_table$sequence, gc_percent, numeric(1))
safe_write_csv(sequence_table, "native_vs_mutated_49bp_sequences_repository.csv")

native_units <- find_canonical_hse_units(sequence_table$sequence[[1]])
native_units$unit_id <- paste0("HSE", seq_len(nrow(native_units)))
safe_write_csv(native_units, "native_49bp_HSE_units_repository.csv")

# Prefer exact processed scan outputs used in the manuscript.
processed_hits_file <- find_any_input(c(
  "native_vs_mutated_49bp_TF_motif_hits_annotated.csv",
  "native_vs_mutated_49bp_TF_motif_hits.csv"
), required = FALSE)

processed_scores_file <- find_any_input(c(
  "native_vs_mutated_49bp_ALL_window_scores_unfiltered.csv",
  "native_vs_mutated_49bp_best_score_per_motif.csv"
), required = FALSE)

normalize_hits <- function(x) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("sequence_id", c("sequence_id", "sequence_type", "sequence"))
  rename_first("motif_id", c("motif_id", "tf_id", "jaspar_id"))
  rename_first("motif_name", c("motif_name", "tf_name", "model_name"))
  rename_first("tf_family", c("tf_family", "family"))
  rename_first("relative_score", c("relative_score", "rel_score", "score"))
  rename_first("strand", c("strand"))
  rename_first("forward_start", c("forward_start", "motif_start", "start", "position"))
  rename_first("forward_end", c("forward_end", "motif_end", "end"))
  rename_first("motif_width", c("motif_width", "width", "motif_length"))
  rename_first("matched_sequence", c("matched_sequence", "match_sequence", "best_match"))

  assert_columns(x, c("sequence_id", "motif_id", "motif_name", "relative_score"), "49-bp motif hits")
  if (!"tf_family" %in% names(x)) x$tf_family <- classify_tf_family(x$motif_name)
  if (!"motif_width" %in% names(x) && "matched_sequence" %in% names(x)) {
    x$motif_width <- nchar(x$matched_sequence)
  }
  if (!"forward_end" %in% names(x) && all(c("forward_start", "motif_width") %in% names(x))) {
    x$forward_end <- as.integer(x$forward_start) + as.integer(x$motif_width) - 1L
  }
  x
}

normalize_scores <- function(x) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("sequence_id", c("sequence_id", "sequence_type", "sequence"))
  rename_first("motif_id", c("motif_id", "tf_id", "jaspar_id"))
  rename_first("motif_name", c("motif_name", "tf_name", "model_name"))
  rename_first("tf_family", c("tf_family", "family"))
  rename_first("relative_score", c("relative_score", "best_relative_score", "score"))
  rename_first("forward_start", c("forward_start", "motif_start", "best_start", "start", "position"))
  rename_first("strand", c("strand", "best_strand"))
  rename_first("matched_sequence", c("matched_sequence", "best_sequence", "best_match"))
  assert_columns(x, c("sequence_id", "motif_id", "motif_name", "relative_score"), "49-bp motif scores")
  if (!"tf_family" %in% names(x)) x$tf_family <- classify_tf_family(x$motif_name)
  x
}

if (!is.na(processed_hits_file) && !ONSEN_FORCE_RESCAN) {
  message("Using exact processed motif-hit file: ", processed_hits_file)
  hits <- normalize_hits(read_table_auto(processed_hits_file))
} else {
  jaspar_file <- find_any_input(c(
    "JASPAR2026_CORE_plants_nonredundant_pfms_jaspar.txt",
    "JASPAR2024_CORE_plants_nonredundant_pfms_jaspar.txt"
  ))
  motifs <- parse_jaspar_pfms(jaspar_file)
  hits <- as.data.frame(scan_sequences_against_motifs(
    sequence_table, motifs, threshold = 0.85, pseudocount = 0.8, retain_all = FALSE
  ))
}

if (!is.na(processed_scores_file) && !ONSEN_FORCE_RESCAN) {
  message("Using exact processed motif-score file: ", processed_scores_file)
  scores_raw <- normalize_scores(read_table_auto(processed_scores_file))
  # If the input contains every scanned position, retain the best occurrence per sequence/model.
  scores <- scores_raw |>
    dplyr::group_by(sequence_id, motif_id, motif_name, tf_family) |>
    dplyr::slice_max(relative_score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
} else {
  if (!exists("motifs")) {
    jaspar_file <- find_any_input(c(
      "JASPAR2026_CORE_plants_nonredundant_pfms_jaspar.txt",
      "JASPAR2024_CORE_plants_nonredundant_pfms_jaspar.txt"
    ))
    motifs <- parse_jaspar_pfms(jaspar_file)
  }
  score_rows <- list()
  k <- 0L
  for (i in seq_len(nrow(sequence_table))) {
    for (motif in motifs) {
      k <- k + 1L
      score_rows[[k]] <- best_score_one_motif(
        sequence_table$sequence_id[[i]],
        sequence_table$sequence[[i]],
        motif,
        pseudocount = 0.8
      )
    }
  }
  scores <- data.table::rbindlist(score_rows, fill = TRUE) |> as.data.frame()
  scores$tf_family <- classify_tf_family(scores$motif_name)
}

hits$relative_score <- as.numeric(hits$relative_score)
scores$relative_score <- as.numeric(scores$relative_score)

safe_write_csv(hits, "native_vs_mutated_49bp_TF_motif_hits_repository.csv")
safe_write_csv(scores, "native_vs_mutated_49bp_best_score_per_motif_repository.csv")

family_summary <- hits |>
  dplyr::filter(relative_score >= 0.85) |>
  dplyr::count(sequence_id, tf_family, name = "high_confidence_motif_position_hits") |>
  tidyr::complete(
    sequence_id = sequence_table$sequence_id,
    tf_family,
    fill = list(high_confidence_motif_position_hits = 0L)
  ) |>
  dplyr::arrange(tf_family, sequence_id)

safe_write_csv(family_summary, "native_vs_mutated_49bp_family_summary_repository.csv")

best_wide <- scores |>
  dplyr::select(sequence_id, motif_id, motif_name, tf_family, relative_score,
                dplyr::any_of(c("forward_start", "strand", "matched_sequence"))) |>
  tidyr::pivot_wider(
    names_from = sequence_id,
    values_from = c(relative_score, dplyr::any_of(c("forward_start", "strand", "matched_sequence"))),
    names_sep = "__"
  )

native_score_col <- grep("^relative_score__Native HSE window$", names(best_wide), value = TRUE)
designed_score_col <- grep("^relative_score__HSE-mutated window$", names(best_wide), value = TRUE)
native_pos_col <- grep("^forward_start__Native HSE window$", names(best_wide), value = TRUE)
designed_pos_col <- grep("^forward_start__HSE-mutated window$", names(best_wide), value = TRUE)

if (length(native_score_col) != 1L || length(designed_score_col) != 1L) {
  stop("Could not construct native/designed best-score columns.")
}

best_wide$native_best_relative_score <- best_wide[[native_score_col]]
best_wide$designed_best_relative_score <- best_wide[[designed_score_col]]
best_wide$delta_designed_minus_native <- best_wide$designed_best_relative_score -
  best_wide$native_best_relative_score

best_wide$effect_class <- vapply(seq_len(nrow(best_wide)), function(i) {
  classify_motif_effect(
    best_wide$native_best_relative_score[[i]],
    best_wide$designed_best_relative_score[[i]],
    threshold = 0.85,
    effect_margin = 0.05,
    native_position = if (length(native_pos_col)) best_wide[[native_pos_col]][[i]] else NA_integer_,
    designed_position = if (length(designed_pos_col)) best_wide[[designed_pos_col]][[i]] else NA_integer_
  )
}, character(1))

effect_summary <- best_wide |>
  dplyr::count(effect_class, name = "motif_models")

safe_write_csv(best_wide, "native_vs_mutated_49bp_motif_effects_repository.csv")
safe_write_csv(effect_summary, "native_vs_mutated_49bp_effect_summary_repository.csv")

# ------------------------------- Fig. 1A -------------------------------------
if (ONSEN_MAKE_FIGURES) {
  native_bases <- split_bases(sequence_table$sequence[[1]])
  designed_bases <- split_bases(sequence_table$sequence[[2]])
  substitutions <- which(native_bases != designed_bases)
  hse_positions <- unlist(Map(seq, native_units$start, native_units$end))

  seq_plot <- do.call(rbind, lapply(seq_len(nrow(sequence_table)), function(i) {
    data.frame(
      sequence_id = sequence_table$sequence_id[[i]],
      position = seq_along(native_bases),
      base = split_bases(sequence_table$sequence[[i]]),
      stringsAsFactors = FALSE
    )
  }))
  seq_plot$tile_class <- "Other sequence"
  seq_plot$tile_class[
    seq_plot$sequence_id == "Native HSE window" & seq_plot$position %in% hse_positions
  ] <- "Native HSE-like unit"
  seq_plot$tile_class[
    seq_plot$sequence_id == "HSE-mutated window" & seq_plot$position %in% substitutions
  ] <- "Substituted base"
  seq_plot$sequence_id <- factor(
    seq_plot$sequence_id,
    levels = c("Native HSE window", "HSE-mutated window")
  )

  p1a <- ggplot2::ggplot(seq_plot, ggplot2::aes(position, sequence_id)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = tile_class), colour = "black",
      linewidth = 0.25, width = 0.95, height = 0.70
    ) +
    ggplot2::geom_text(ggplot2::aes(label = base), size = 3.1, fontface = "bold") +
    ggplot2::scale_fill_manual(values = c(
      "Native HSE-like unit" = "#F6D36F",
      "Substituted base" = "#E99AB8",
      "Other sequence" = "#F2F2F2"
    )) +
    ggplot2::scale_x_continuous(
      breaks = c(1, 10, 20, 30, 40, 49),
      limits = c(0.5, 49.5), expand = c(0, 0)
    ) +
    ggplot2::labs(x = "Position in 49-bp window", y = NULL, fill = NULL) +
    theme_onsen(13) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.y = ggplot2::element_text(face = "bold"),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  save_plot_pair(p1a, "Fig1A_native_designed_sequence", 10, 3.5)

  # ----------------------------- Fig. 1B -------------------------------------
  main_families <- c("HSF", "AP2/ERF", "DOF", "GATA", "MYB", "NAC",
                     "bZIP", "C2H2/ZAT", "WRKY", "ARF", "HD-ZIP", "LBD")
  p1b_data <- family_summary |>
    dplyr::filter(tf_family %in% main_families) |>
    dplyr::mutate(
      tf_family = factor(tf_family, levels = main_families),
      sequence_id = factor(sequence_id, levels = sequence_table$sequence_id)
    )

  p1b <- ggplot2::ggplot(
    p1b_data,
    ggplot2::aes(tf_family, high_confidence_motif_position_hits, fill = sequence_id)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.76),
      width = 0.68, colour = "black", linewidth = 0.25
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = high_confidence_motif_position_hits),
      position = ggplot2::position_dodge(width = 0.76),
      vjust = -0.25, size = 3.3
    ) +
    ggplot2::scale_fill_manual(values = c(
      "Native HSE window" = "#3A2D8F",
      "HSE-mutated window" = "#0B9A6B"
    )) +
    ggplot2::labs(
      x = "TF family", y = "High-confidence motif-position hits", fill = "Sequence"
    ) +
    theme_onsen(13) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  save_plot_pair(p1b, "Fig1B_family_hit_counts", 8.2, 5.2)

  # Shared heatmap function for Fig. 2 panels.
  best_long <- scores |>
    dplyr::select(sequence_id, motif_id, motif_name, tf_family, relative_score)

  make_score_heatmap <- function(data, stem, x_label, width, height) {
    data$sequence_id <- factor(
      data$sequence_id,
      levels = c("HSE-mutated window", "Native HSE window")
    )
    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(motif_name, sequence_id, fill = relative_score)
    ) +
      ggplot2::geom_tile(
        ggplot2::aes(colour = relative_score >= 0.85),
        linewidth = 0.8
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.3f", relative_score)),
        size = 2.8, fontface = "bold"
      ) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA), guide = "none") +
      ggplot2::scale_fill_gradientn(
        colours = c("#F3A6A6", "#F6D07A", "#D9D8E5", "#7662F1"),
        limits = c(0.60, 1.00), oob = scales::squish,
        name = "Best relative\nPWM score"
      ) +
      ggplot2::labs(x = x_label, y = NULL) +
      theme_onsen(12) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )
    save_plot_pair(p, stem, width, height)
  }

  hsf_models <- read.csv(repo_file("Arabidopsis_HSF_models_JASPAR2026.csv"))
  fig2a <- best_long |>
    dplyr::filter(toupper(motif_name) %in% toupper(hsf_models$HSF_model))
  make_score_heatmap(fig2a, "Fig2A_HSF_best_scores", "Arabidopsis HSF-family motif model", 8.2, 3.8)

  changed_ap2 <- best_wide |>
    dplyr::filter(tf_family == "AP2/ERF", effect_class %in% c("gained", "strengthened")) |>
    dplyr::arrange(dplyr::desc(delta_designed_minus_native)) |>
    dplyr::slice_head(n = 15) |>
    dplyr::pull(motif_id)
  fig2b <- best_long |> dplyr::filter(motif_id %in% changed_ap2)
  make_score_heatmap(fig2b, "Fig2B_AP2ERF_best_scores", "AP2/ERF-family motif model", 10, 3.8)

  dof_models <- best_long |>
    dplyr::filter(tf_family == "DOF") |>
    dplyr::group_by(motif_id, motif_name) |>
    dplyr::filter(max(relative_score, na.rm = TRUE) >= 0.85) |>
    dplyr::ungroup()
  make_score_heatmap(dof_models, "Fig2C_DOF_best_scores", "DOF-family motif model", 7.5, 3.8)

  other_ids <- best_wide |>
    dplyr::filter(
      !tf_family %in% c("HSF", "AP2/ERF", "DOF", "Other"),
      pmax(native_best_relative_score, designed_best_relative_score, na.rm = TRUE) >= 0.85
    ) |>
    dplyr::pull(motif_id)
  fig2d <- best_long |>
    dplyr::filter(motif_id %in% other_ids) |>
    dplyr::mutate(motif_name = paste(tf_family, motif_name, sep = " | "))
  make_score_heatmap(fig2d, "Fig2D_additional_TF_best_scores", "TF motif model", 13, 4.2)
}

message("Native/designed 49-bp motif analysis completed.")
