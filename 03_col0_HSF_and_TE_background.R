# Col-0 ONSEN HSF scan and strict TE-only background analysis.
# Covers Fig. 3 and source data for Tables S4-S5.

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "tidyr", "readr", "stringr",
                   "purrr", "ggplot2", "forcats", "Biostrings", "readxl"))

message_config()

thresholds <- c(0.85, 0.90)
window_coords <- read.csv(repo_file("ONSEN_Col0_terminal_candidate_windows.csv"))
window_coords$chromosome <- standard_chr(window_coords$chromosome)

normalize_hsf_summary <- function(x) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("sequence_id", c("sequence_id", "candidate_id", "LTR_candidate", "ltr_candidate"))
  rename_first("HSF_hits", c("HSF_hits", "hsf_hits", "HSF_motif_position_hits",
                            "HSF motif-position hits", "motif_position_hits"))
  rename_first("unique_HSF_models", c("unique_HSF_models", "unique_hsf_models",
                                     "Unique HSF models"))
  rename_first("maximum_relative_score", c("maximum_relative_score", "max_relative_score",
                                          "Maximum relative score"))
  rename_first("width_bp", c("width_bp", "LTR_window_width_bp", "window_width_bp",
                            "LTR window width (bp)"))
  rename_first("HSF_hits_per_kb", c("HSF_hits_per_kb", "hsf_hits_per_kb",
                                   "HSF_motif_density_hits_per_kb",
                                   "HSF motif density (hits/kb)"))
  rename_first("threshold", c("threshold", "relative_score_threshold"))
  assert_columns(x, c("sequence_id", "HSF_hits"), "HSF summary")
  if (!"width_bp" %in% names(x)) x$width_bp <- 800L
  if (!"HSF_hits_per_kb" %in% names(x)) {
    x$HSF_hits_per_kb <- as.numeric(x$HSF_hits) / as.numeric(x$width_bp) * 1000
  }
  x
}

normalize_hits <- function(x, sequence_width_lookup = NULL) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("sequence_id", c("sequence_id", "candidate_id"))
  rename_first("motif_id", c("motif_id", "tf_id", "jaspar_id"))
  rename_first("motif_name", c("motif_name", "tf_name"))
  rename_first("strand", c("strand"))
  rename_first("scan_start", c("scan_start", "position", "start", "motif_start"))
  rename_first("motif_width", c("motif_width", "motif_length", "width"))
  rename_first("matched_sequence", c("matched_sequence", "match_sequence"))
  rename_first("relative_score", c("relative_score", "score"))
  assert_columns(x, c("sequence_id", "motif_id", "motif_name", "strand",
                      "scan_start", "relative_score"), "HSF hit table")
  if (!"motif_width" %in% names(x) && "matched_sequence" %in% names(x)) {
    x$motif_width <- nchar(x$matched_sequence)
  }
  x
}

copy_summary_file <- find_any_input(c(
  "Col0_ONSEN_LTRcandidate_JASPAR2026_Arabidopsis_HSF_summary.csv",
  "Col0_ONSEN_LTRcandidate_HSF_summary_complete.csv"
), required = FALSE)

copy_hits_file <- find_any_input(c(
  "Col0_ONSEN_LTRcandidate_JASPAR2026_Arabidopsis_HSF_hits.csv",
  "Col0_ONSEN_LTRcandidate_high_confidence_motif_hits_COMBINED.csv"
), required = FALSE)

background_coords_file <- find_input("strict_TE_only_background_coordinates.csv")
background_summary_file <- find_any_input(c(
  "ONSEN_vs_strict_TE_only_background_threshold_0p85_0p90_summary.csv",
  "ONSEN_vs_strict_TE_only_background_threshold_0p85_0p90_STATS.csv"
), required = FALSE)

# A scanner is used only when exact processed hit summaries are unavailable or
# ONSEN_FORCE_RESCAN=true.
scan_from_coordinates <- function(coordinates, sequence_id_col, thresholds) {
  if (!ONSEN_RUN_LARGE_STEPS && !ONSEN_FORCE_RESCAN) {
    stop(
      "Large HSF rescanning is disabled and the exact processed summary was not found.\n",
      "Set ONSEN_RUN_LARGE_STEPS=true or supply the processed HSF summary."
    )
  }

  genome_file <- find_input("Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz")
  jaspar_file <- find_input("JASPAR2026_CORE_plants_nonredundant_pfms_jaspar.txt")
  genome <- Biostrings::readDNAStringSet(genome_file)
  names(genome) <- standard_chr(names(genome))

  hsf_names <- read.csv(repo_file("Arabidopsis_HSF_models_JASPAR2026.csv"))$HSF_model
  motifs_all <- parse_jaspar_pfms(jaspar_file)
  motif_names <- toupper(vapply(motifs_all, `[[`, character(1), "name"))
  hsf_motifs <- motifs_all[motif_names %in% toupper(hsf_names)]
  if (!length(hsf_motifs)) {
    hsf_motifs <- motifs_all[grepl("^HSF", motif_names)]
  }
  if (!length(hsf_motifs)) stop("No Arabidopsis HSF models found in JASPAR file.")

  extract_sequence <- function(chr, start, end) {
    chr <- standard_chr(chr)
    if (!chr %in% names(genome)) stop("Chromosome absent from FASTA: ", chr)
    as.character(Biostrings::subseq(genome[[chr]], start = start, end = end))
  }

  sequences <- data.frame(
    sequence_id = coordinates[[sequence_id_col]],
    sequence = purrr::pmap_chr(
      list(coordinates$chromosome, coordinates$start, coordinates$end),
      extract_sequence
    ),
    stringsAsFactors = FALSE
  )

  all_summaries <- list()
  all_hits <- list()
  for (cutoff in thresholds) {
    hits <- as.data.frame(scan_sequences_against_motifs(
      sequences, hsf_motifs, threshold = cutoff, pseudocount = 0.8
    ))
    summary <- hits |>
      dplyr::group_by(sequence_id) |>
      dplyr::summarise(
        HSF_hits = dplyr::n(),
        unique_HSF_models = dplyr::n_distinct(motif_id),
        maximum_relative_score = max(relative_score, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::right_join(
        data.frame(
          sequence_id = coordinates[[sequence_id_col]],
          width_bp = coordinates$end - coordinates$start + 1L
        ),
        by = "sequence_id"
      ) |>
      dplyr::mutate(
        HSF_hits = tidyr::replace_na(HSF_hits, 0L),
        unique_HSF_models = tidyr::replace_na(unique_HSF_models, 0L),
        HSF_hits_per_kb = HSF_hits / width_bp * 1000,
        threshold = cutoff
      )
    hits$threshold <- cutoff
    all_summaries[[as.character(cutoff)]] <- summary
    all_hits[[as.character(cutoff)]] <- hits
  }

  list(
    summary = dplyr::bind_rows(all_summaries),
    hits = dplyr::bind_rows(all_hits)
  )
}

# ---------------------- Col-0 sixteen terminal windows -----------------------
if (!is.na(copy_summary_file) && !ONSEN_FORCE_RESCAN) {
  copy_summary_raw <- normalize_hsf_summary(read_table_auto(copy_summary_file))
  if (!"threshold" %in% names(copy_summary_raw)) copy_summary_raw$threshold <- 0.85
  copy_summary <- copy_summary_raw |>
    dplyr::left_join(
      window_coords |>
        dplyr::select(window_id, copy_id, terminal_side, chromosome,
                      start_1based, end_1based, width_bp),
      by = c("sequence_id" = "window_id")
    )
  if (any(copy_summary$threshold == 0.90)) {
    copy_summary_all <- copy_summary
  } else {
    # Use the exact project threshold-0.90 summary when available.
    strict_copy_file <- find_input(
      "Col0_ONSEN_LTRcandidate_HSF_summary_threshold_0p90.csv",
      required = FALSE
    )
    if (!is.na(strict_copy_file)) {
      strict_copy <- normalize_hsf_summary(read_table_auto(strict_copy_file))
      strict_copy$threshold <- 0.90
      strict_copy <- strict_copy |>
        dplyr::left_join(
          window_coords |>
            dplyr::select(window_id, copy_id, terminal_side, chromosome,
                          start_1based, end_1based, width_bp),
          by = c("sequence_id" = "window_id")
        )
      copy_summary_all <- dplyr::bind_rows(copy_summary, strict_copy)
    } else {
      copy_summary_all <- copy_summary
    }
  }
} else {
  scan_coords <- window_coords |>
    dplyr::transmute(
      sequence_id = window_id,
      chromosome,
      start = start_1based,
      end = end_1based
    )
  scanned <- scan_from_coordinates(scan_coords, "sequence_id", thresholds)
  copy_summary_all <- scanned$summary |>
    dplyr::left_join(
      window_coords |>
        dplyr::select(window_id, copy_id, terminal_side, chromosome,
                      start_1based, end_1based, width_bp),
      by = c("sequence_id" = "window_id")
    )
  copy_hits <- scanned$hits
}

if (!is.na(copy_hits_file) && !ONSEN_FORCE_RESCAN) {
  copy_hits <- normalize_hits(read_table_auto(copy_hits_file))
  if (!"threshold" %in% names(copy_hits)) copy_hits$threshold <- 0.85
}
safe_write_csv(copy_summary_all, "Col0_ONSEN_HSF_summary_repository.csv")
if (exists("copy_hits")) safe_write_csv(copy_hits, "Col0_ONSEN_HSF_hits_repository.csv")

# ------------------------- Strict TE-only background -------------------------
background_coords <- read_table_auto(background_coords_file)

rename_coordinate_columns <- function(x) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("sequence_id", c("sequence_id", "candidate_id", "region_id", "feature_id", "GeneID"))
  rename_first("chromosome", c("chromosome", "chr", "Chr", "seqnames"))
  rename_first("start", c("start", "Start", "region_start"))
  rename_first("end", c("end", "End", "region_end"))
  rename_first("TE_family", c("TE_family", "te_family", "family", "Alias", "alias"))
  assert_columns(x, c("chromosome", "start", "end"), "strict TE coordinates")
  if (!"sequence_id" %in% names(x)) x$sequence_id <- paste0("strict_TE_", seq_len(nrow(x)))
  x$chromosome <- standard_chr(x$chromosome)
  x$start <- as.integer(x$start)
  x$end <- as.integer(x$end)
  x
}
background_coords <- rename_coordinate_columns(background_coords)

if (nrow(background_coords) != 1942L) {
  warning("The final manuscript background contained 1,942 regions; current coordinate file contains ",
          nrow(background_coords), ".")
}

processed_bg_summary_file <- find_any_input(c(
  "ONSEN_vs_strict_TE_only_background_threshold_0p85_0p90_summary.csv",
  "strict_TE_background_HSF_summary_thresholds_0p85_0p90.csv"
), required = FALSE)

if (!is.na(processed_bg_summary_file) && !ONSEN_FORCE_RESCAN) {
  combined_processed <- read_table_auto(processed_bg_summary_file)
  # This file may contain class-level statistics rather than region-level data.
  region_bg_file <- find_any_input(c(
    "strict_TE_background_annotated_with_HSF_outlier_classes.csv",
    "strict_TE_background_HSF_counts_threshold_0p85.csv",
    "ONSEN_vs_strict_TE_background_HSF_summary_thresholds_0p85_0p90.csv"
  ), required = FALSE)

  if (!is.na(region_bg_file)) {
    region_bg <- read_table_auto(region_bg_file)
    region_bg <- normalize_hsf_summary(region_bg)
    if (!"threshold" %in% names(region_bg)) region_bg$threshold <- 0.85
    background_summary_all <- region_bg
  } else {
    background_summary_all <- data.frame()
  }
} else {
  scanned_bg <- scan_from_coordinates(background_coords, "sequence_id", thresholds)
  background_summary_all <- scanned_bg$summary
  background_hits <- scanned_bg$hits
  safe_write_csv(background_hits, "strict_TE_background_HSF_hits_repository.csv")
}

# If exact region-level background counts are not found, the published class
# statistics are retained from Table S5 and the processed project summaries.
if (nrow(background_summary_all)) {
  safe_write_csv(background_summary_all, "strict_TE_background_HSF_summary_repository.csv")
}

# ------------------------------- Statistics ----------------------------------
# Prefer exact class-level project statistics when available.
stats_file <- find_any_input(c(
  "ONSEN_vs_strict_TE_only_background_threshold_0p85_0p90_STATS.csv",
  "ONSEN_vs_strict_TE_background_HSF_stats.csv"
), required = FALSE)

if (!is.na(stats_file) && !ONSEN_FORCE_RESCAN) {
  stats_table <- read_table_auto(stats_file)
} else if (nrow(background_summary_all)) {
  stats_rows <- list()
  for (cutoff in intersect(thresholds, unique(background_summary_all$threshold))) {
    onsen <- copy_summary_all$HSF_hits_per_kb[copy_summary_all$threshold == cutoff]
    background <- background_summary_all$HSF_hits_per_kb[
      background_summary_all$threshold == cutoff
    ]
    test <- stats::wilcox.test(onsen, background, exact = FALSE)
    stats_rows[[as.character(cutoff)]] <- data.frame(
      threshold = cutoff,
      n_ONSEN = length(onsen),
      n_background = length(background),
      ONSEN_median_hits_per_kb = median(onsen, na.rm = TRUE),
      background_median_hits_per_kb = median(background, na.rm = TRUE),
      Wilcoxon_P = test$p.value,
      Cliffs_delta = cliffs_delta(onsen, background),
      background_regions_at_or_above_ONSEN_median =
        sum(background >= median(onsen, na.rm = TRUE), na.rm = TRUE),
      background_regions_at_or_above_ONSEN_maximum =
        sum(background >= max(onsen, na.rm = TRUE), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  stats_table <- dplyr::bind_rows(stats_rows)
} else {
  # Table S5 is included in the repository and remains the authoritative
  # processed summary when region-level background values are unavailable.
  stats_table <- readxl::read_excel(repo_file("Table_S5.xlsx"), skip = 1)
}
safe_write_csv(stats_table, "ONSEN_vs_strict_TE_HSF_statistics_repository.csv")

# --------------------------- Outlier family counts ---------------------------
outlier_file <- find_any_input(c(
  "strict_TE_background_annotated_with_HSF_outlier_classes.csv",
  "top_strict_TE_only_background_HSF_rich_regions_threshold_0p85.csv"
), required = FALSE)

if (!is.na(outlier_file)) {
  outliers <- read_table_auto(outlier_file)
  family_col <- c("TE_family", "te_family", "family", "Alias", "alias")[
    c("TE_family", "te_family", "family", "Alias", "alias") %in% names(outliers)
  ][1]
  class_col <- c("HSF_outlier_class", "outlier_class", "is_HSF_rich_outlier")[
    c("HSF_outlier_class", "outlier_class", "is_HSF_rich_outlier") %in% names(outliers)
  ][1]
  if (!is.na(class_col)) {
    keep <- !is.na(outliers[[class_col]]) &
      !tolower(as.character(outliers[[class_col]])) %in% c("no", "false", "not_outlier", "background")
    outliers <- outliers[keep, , drop = FALSE]
  }
  if (!is.na(family_col)) {
    outlier_family <- outliers |>
      dplyr::count(TE_family = .data[[family_col]], name = "n_HSF_rich_outlier_regions") |>
      dplyr::arrange(dplyr::desc(n_HSF_rich_outlier_regions))
    safe_write_csv(outlier_family, "HSF_rich_TE_outlier_family_counts_repository.csv")
  }
}

# -------------------------------- Figures ------------------------------------
if (ONSEN_MAKE_FIGURES) {
  copy_085 <- copy_summary_all |>
    dplyr::filter(threshold == 0.85) |>
    dplyr::mutate(
      copy_id = ifelse(is.na(copy_id), sub("-[LR]$", "", sequence_id), copy_id),
      terminal_side = ifelse(
        is.na(terminal_side),
        ifelse(grepl("-L$", sequence_id), "Left LTR", "Right LTR"),
        ifelse(tolower(terminal_side) == "left", "Left LTR", "Right LTR")
      ),
      copy_id = factor(copy_id, levels = paste0("ONSEN", 1:8)),
      terminal_side = factor(terminal_side, levels = c("Left LTR", "Right LTR"))
    )

  p3a <- ggplot2::ggplot(
    copy_085,
    ggplot2::aes(copy_id, HSF_hits, fill = terminal_side)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.76),
      width = 0.68, colour = "black", linewidth = 0.25
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = HSF_hits),
      position = ggplot2::position_dodge(width = 0.76),
      vjust = -0.25, size = 3.2
    ) +
    ggplot2::scale_fill_manual(values = c("Left LTR" = "#A67BE8", "Right LTR" = "#ED82BD")) +
    ggplot2::labs(
      x = "ONSEN copy", y = "HSF-family motif-position hits", fill = "LTR candidate"
    ) +
    theme_onsen(13) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  save_plot_pair(p3a, "Fig3A_Col0_ONSEN_HSF_hits", 7.5, 5.0)

  paired <- copy_085 |>
    dplyr::select(copy_id, terminal_side, HSF_hits) |>
    tidyr::pivot_wider(names_from = terminal_side, values_from = HSF_hits) |>
    tidyr::pivot_longer(
      cols = c("Left LTR", "Right LTR"),
      names_to = "terminal_side", values_to = "HSF_hits"
    )
  p3b <- ggplot2::ggplot(
    paired,
    ggplot2::aes(terminal_side, HSF_hits, group = copy_id)
  ) +
    ggplot2::geom_line(colour = "grey40", linewidth = 0.6) +
    ggplot2::geom_point(
      ggplot2::aes(fill = terminal_side),
      shape = 21, size = 3.4, colour = "black"
    ) +
    ggplot2::scale_fill_manual(values = c("Left LTR" = "#A67BE8", "Right LTR" = "#ED82BD")) +
    ggplot2::labs(x = "LTR candidate", y = "HSF-family motif-position hits") +
    theme_onsen(13) +
    ggplot2::theme(legend.position = "none")
  save_plot_pair(p3b, "Fig3B_left_right_HSF_hits", 5.5, 5.0)

  if (nrow(background_summary_all)) {
    bg_085 <- background_summary_all |>
      dplyr::filter(threshold == 0.85) |>
      dplyr::transmute(class = "Strict TE-only background", HSF_hits_per_kb)
    onsen_085 <- copy_085 |>
      dplyr::transmute(class = "ONSEN LTR candidates", HSF_hits_per_kb)
    density_plot <- dplyr::bind_rows(bg_085, onsen_085)
    density_plot$class <- factor(
      density_plot$class,
      levels = c("Strict TE-only background", "ONSEN LTR candidates")
    )

    p3c <- ggplot2::ggplot(
      density_plot,
      ggplot2::aes(class, HSF_hits_per_kb, fill = class)
    ) +
      ggplot2::geom_violin(trim = TRUE, alpha = 0.65, colour = "black") +
      ggplot2::geom_boxplot(width = 0.24, outlier.shape = NA, colour = "black") +
      ggplot2::geom_jitter(width = 0.12, size = 1.1, alpha = 0.45) +
      ggplot2::scale_fill_manual(values = c(
        "Strict TE-only background" = "#83AEE8",
        "ONSEN LTR candidates" = "#E989AE"
      )) +
      ggplot2::labs(x = NULL, y = "HSF motif-position hits per kb") +
      theme_onsen(13) +
      ggplot2::theme(legend.position = "none")
    save_plot_pair(p3c, "Fig3C_ONSEN_vs_strict_TE_background", 6.2, 5.2)
  }

  if (exists("outlier_family")) {
    p3d <- ggplot2::ggplot(
      outlier_family,
      ggplot2::aes(
        n_HSF_rich_outlier_regions,
        forcats::fct_reorder(TE_family, n_HSF_rich_outlier_regions)
      )
    ) +
      ggplot2::geom_col(colour = "black", fill = "#9C2DB4", linewidth = 0.25) +
      ggplot2::labs(
        x = "Number of HSF-rich TE outlier regions",
        y = "TE family / alias"
      ) +
      theme_onsen(13)
    save_plot_pair(p3d, "Fig3D_HSF_rich_TE_outlier_families", 6.2, 4.8)
  }
}

message("Col-0 ONSEN and strict TE-only HSF analysis completed.")
