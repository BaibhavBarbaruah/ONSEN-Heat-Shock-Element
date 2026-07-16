# Public Col-0 leaf methylation analysis.
# Covers Fig. 4 and source data for Table S6.
# The methylome is basal/unstressed and is not interpreted as heat-induced change.

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "tidyr", "readr", "stringr",
                   "purrr", "ggplot2", "forcats", "scales",
                   "GenomicRanges", "IRanges", "S4Vectors"))

message_config()

normalize_coordinates <- function(x, default_class = NA_character_) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("locus_id", c("locus_id", "sequence_id", "candidate_id", "region_id",
                            "window_id", "feature_id", "GeneID", "plot_label"))
  rename_first("chromosome", c("chromosome", "chr", "Chr", "seqnames", "chrom"))
  rename_first("start", c("start", "Start", "start_1based", "region_start"))
  rename_first("end", c("end", "End", "end_1based", "region_end"))
  rename_first("locus_class", c("locus_class", "candidate_class", "class", "region_class"))
  rename_first("HSF_hits_per_kb", c("HSF_hits_per_kb", "hsf_hits_per_kb",
                                   "HSF_motif_density_hits_per_kb",
                                   "HSF motif density (hits/kb)"))
  assert_columns(x, c("chromosome", "start", "end"), "coordinate table")
  if (!"locus_id" %in% names(x)) x$locus_id <- paste0("region_", seq_len(nrow(x)))
  if (!"locus_class" %in% names(x)) x$locus_class <- default_class
  x$chromosome <- standard_chr(x$chromosome)
  x$start <- as.integer(x$start)
  x$end <- as.integer(x$end)
  x
}

normalize_locus_methylation <- function(x) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("locus_id", c("locus_id", "candidate_id", "region_id", "plot_label"))
  rename_first("locus_class", c("locus_class", "candidate_class", "class"))
  rename_first("methylation_context", c("methylation_context", "context", "mc_context"))
  rename_first("weighted_methylation_percent", c(
    "weighted_methylation_percent", "weighted_methylation",
    "weighted_methylation_pct", "methylation_percent"
  ))
  rename_first("HSF_hits_per_kb", c("HSF_hits_per_kb", "hsf_hits_per_kb",
                                   "HSF motif density (hits/kb)"))
  assert_columns(
    x,
    c("locus_id", "locus_class", "methylation_context", "weighted_methylation_percent"),
    "locus methylation summary"
  )
  x$methylation_context <- toupper(x$methylation_context)
  x$weighted_methylation_percent <- as.numeric(x$weighted_methylation_percent)
  x
}

read_methylome <- function(path) {
  mc <- data.table::fread(path, data.table = TRUE)
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(mc)][1]
    if (!is.na(hit) && target != hit) data.table::setnames(mc, hit, target)
  }
  rename_first("chromosome", c("chromosome", "chrom", "chr", "Chr"))
  rename_first("position", c("position", "pos", "start"))
  rename_first("mc_class", c("mc_class", "class", "context"))
  rename_first("methylated_bases", c("methylated_bases", "methylated", "mc_count", "num_c"))
  rename_first("total_bases", c("total_bases", "total", "coverage", "cov", "num_total"))
  assert_columns(
    mc,
    c("chromosome", "position", "mc_class", "methylated_bases", "total_bases"),
    "methylome"
  )
  mc[, chromosome := standard_chr(chromosome)]
  mc[, position := as.integer(position)]
  mc[, methylated_bases := as.numeric(methylated_bases)]
  mc[, total_bases := as.numeric(total_bases)]
  mc <- mc[is.finite(methylated_bases) & is.finite(total_bases) & total_bases > 0]
  mc[, methylation_context := data.table::fifelse(
    substr(mc_class, 2, 2) == "G", "CG",
    data.table::fifelse(substr(mc_class, 3, 3) == "G", "CHG", "CHH")
  )]
  mc
}

summarize_calls_over_regions <- function(mc, regions) {
  region_dt <- data.table::as.data.table(regions)
  region_dt[, start := as.integer(start)]
  region_dt[, end := as.integer(end)]
  data.table::setkey(region_dt, chromosome, start, end)

  calls <- mc[, .(
    chromosome, start = position, end = position,
    methylation_context, methylated_bases, total_bases
  )]
  data.table::setkey(calls, chromosome, start, end)

  matched <- data.table::foverlaps(
    calls, region_dt,
    by.x = c("chromosome", "start", "end"),
    by.y = c("chromosome", "start", "end"),
    type = "within", nomatch = 0L
  )

  summary <- matched[, .(
    methylated_bases_sum = sum(methylated_bases, na.rm = TRUE),
    total_bases_sum = sum(total_bases, na.rm = TRUE),
    weighted_methylation_percent =
      weighted_methylation(methylated_bases, total_bases),
    n_methylation_call_rows = .N
  ), by = .(locus_id, locus_class, methylation_context)]

  list(matched = matched, summary = summary)
}

methylome_file <- find_input("GSM1085222_mC_calls_Col_0.tsv.gz", required = FALSE)

# --------------------------- Fig. 4A broad control ----------------------------
broad_locus_file <- find_input(
  "Revision_Fig4_ONSEN_vs_ordinary_TE_Col0_leaf_methylation_locus_summary.csv",
  required = FALSE
)
broad_stats_file <- find_input(
  "Revision_Fig4_ONSEN_vs_ordinary_TE_Wilcoxon_BH_statistics.csv",
  required = FALSE
)

if (!is.na(broad_locus_file) && !ONSEN_FORCE_RESCAN) {
  broad_locus <- normalize_locus_methylation(read_table_auto(broad_locus_file))
  message("Using completed ordinary-TE methylation control: ", broad_locus_file)
} else {
  if (is.na(methylome_file)) stop("Public methylome file is required for raw recalculation.")
  onsen_regions <- read.csv(repo_file("ONSEN_Col0_terminal_candidate_windows.csv")) |>
    dplyr::transmute(
      locus_id = window_id,
      locus_class = "ONSEN LTR candidate",
      chromosome = standard_chr(chromosome),
      start = start_1based,
      end = end_1based
    )

  ordinary_regions <- normalize_coordinates(
    read_table_auto(find_input("strict_TE_only_background_coordinates.csv")),
    default_class = "Ordinary background TE"
  )
  ordinary_regions$locus_class <- "Ordinary background TE"

  # Exclude any direct overlap with ONSEN coordinates.
  gr_onsen <- GenomicRanges::GRanges(
    seqnames = onsen_regions$chromosome,
    ranges = IRanges::IRanges(onsen_regions$start, onsen_regions$end)
  )
  gr_bg <- GenomicRanges::GRanges(
    seqnames = ordinary_regions$chromosome,
    ranges = IRanges::IRanges(ordinary_regions$start, ordinary_regions$end)
  )
  overlap_bg <- unique(S4Vectors::subjectHits(
    GenomicRanges::findOverlaps(gr_onsen, gr_bg)
  ))
  if (length(overlap_bg)) ordinary_regions <- ordinary_regions[-overlap_bg, , drop = FALSE]

  if (nrow(ordinary_regions) != 1942L) {
    warning("Expected 1,942 ordinary TE regions after exclusion; found ", nrow(ordinary_regions))
  }

  mc <- read_methylome(methylome_file)
  broad_result <- summarize_calls_over_regions(mc, dplyr::bind_rows(
    onsen_regions, ordinary_regions |>
      dplyr::select(locus_id, locus_class, chromosome, start, end)
  ))
  broad_locus <- as.data.frame(broad_result$summary)
}

contexts <- c("CG", "CHG", "CHH")
if (!is.na(broad_stats_file) && !ONSEN_FORCE_RESCAN) {
  broad_stats <- read_table_auto(broad_stats_file)
} else {
  broad_stats <- do.call(rbind, lapply(contexts, function(context) {
    z <- broad_locus[broad_locus$methylation_context == context, , drop = FALSE]
    ordinary <- z$weighted_methylation_percent[
      z$locus_class == "Ordinary background TE"
    ]
    onsen <- z$weighted_methylation_percent[
      z$locus_class == "ONSEN LTR candidate"
    ]
    test <- stats::wilcox.test(ordinary, onsen, exact = FALSE)
    data.frame(
      methylation_context = context,
      ordinary_TE_n = length(ordinary),
      ONSEN_n = length(onsen),
      ordinary_TE_median = median(ordinary, na.rm = TRUE),
      ONSEN_median = median(onsen, na.rm = TRUE),
      Wilcoxon_P = test$p.value,
      stringsAsFactors = FALSE
    )
  }))
  broad_stats$BH_adjusted_P <- p.adjust(broad_stats$Wilcoxon_P, method = "BH")
  broad_stats$significance_label <- cut(
    broad_stats$BH_adjusted_P,
    breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
    labels = c("****", "***", "**", "*", "ns")
  )
}

broad_class <- broad_locus |>
  dplyr::group_by(locus_class, methylation_context) |>
  dplyr::summarise(
    n_loci = dplyr::n(),
    median_weighted_methylation_percent =
      median(weighted_methylation_percent, na.rm = TRUE),
    mean_weighted_methylation_percent =
      mean(weighted_methylation_percent, na.rm = TRUE),
    minimum_weighted_methylation_percent =
      min(weighted_methylation_percent, na.rm = TRUE),
    maximum_weighted_methylation_percent =
      max(weighted_methylation_percent, na.rm = TRUE),
    .groups = "drop"
  )

safe_write_csv(broad_locus, "Revision_Fig4_ordinary_TE_methylation_locus_summary_repository.csv")
safe_write_csv(broad_class, "Revision_Fig4_ordinary_TE_methylation_class_summary_repository.csv")
safe_write_csv(broad_stats, "Revision_Fig4_ordinary_TE_methylation_statistics_repository.csv")

# ---------------------- Fig. 4B-D selected HSF-rich loci ---------------------
selected_summary_file <- find_any_input(c(
  "Figure5_candidate_loci_Col0_leaf_methylation_summary_by_context.csv",
  "Figure5_candidate_loci_Col0_leaf_methylation_summary_wide.csv"
), required = FALSE)

selected_coord_file <- find_any_input(c(
  "Figure5_candidate_ONSEN_and_HSF_rich_TE_loci.csv",
  "Figure5_HSF_rich_non_ONSEN_TE_outlier_coordinates_ANNOTATED.csv"
), required = FALSE)

if (!is.na(selected_summary_file) && !ONSEN_FORCE_RESCAN) {
  selected_locus <- normalize_locus_methylation(read_table_auto(selected_summary_file))
} else if (!is.na(selected_coord_file) && !is.na(methylome_file)) {
  if (!exists("mc")) mc <- read_methylome(methylome_file)
  selected_regions <- normalize_coordinates(read_table_auto(selected_coord_file))
  selected_result <- summarize_calls_over_regions(mc, selected_regions)
  selected_locus <- as.data.frame(selected_result$summary)
} else {
  selected_locus <- data.frame()
  warning("Selected HSF-rich locus methylation inputs were not found; Fig. 4B-D may be skipped.")
}

# Join HSF density where possible.
if (nrow(selected_locus) && !"HSF_hits_per_kb" %in% names(selected_locus)) {
  hsf_sources <- c(
    "Figure5_candidate_ONSEN_and_HSF_rich_TE_loci.csv",
    "strict_TE_background_annotated_with_HSF_outlier_classes.csv",
    "Col0_ONSEN_HSF_summary_repository.csv"
  )
  density_tables <- list()
  for (source_name in hsf_sources) {
    path <- find_input(source_name, required = FALSE)
    if (!is.na(path)) {
      z <- tryCatch(normalize_coordinates(read_table_auto(path)), error = function(e) NULL)
      if (!is.null(z) && "HSF_hits_per_kb" %in% names(z)) {
        density_tables[[length(density_tables) + 1L]] <- z |>
          dplyr::select(locus_id, HSF_hits_per_kb)
      }
    }
  }
  if (length(density_tables)) {
    density_table <- dplyr::bind_rows(density_tables) |>
      dplyr::distinct(locus_id, .keep_all = TRUE)
    selected_locus <- selected_locus |>
      dplyr::left_join(density_table, by = "locus_id")
  }
}

if (nrow(selected_locus)) {
  safe_write_csv(selected_locus, "Fig4_selected_loci_methylation_repository.csv")
}

# Aggregate CHH profile. Prefer an existing exact profile table.
profile_candidates <- if (dir.exists(ONSEN_DATA_ROOT)) {
  list.files(
    ONSEN_DATA_ROOT,
    pattern = "CHH.*profile|profile.*CHH",
    recursive = TRUE, full.names = TRUE, ignore.case = TRUE
  )
} else character()
profile_file <- if (length(profile_candidates)) profile_candidates[[1]] else NA_character_

normalize_profile <- function(x) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("locus_class", c("locus_class", "candidate_class", "class"))
  rename_first("bin_index", c("bin_index", "bin", "normalized_bin"))
  rename_first("normalized_position", c("normalized_position", "position_label", "position"))
  rename_first("mean_CHH_methylation", c("mean_CHH_methylation", "mean_methylation",
                                       "CHH_methylation", "mean"))
  rename_first("sem_CHH_methylation", c("sem_CHH_methylation", "sem", "se", "SEM"))
  assert_columns(x, c("locus_class", "bin_index", "mean_CHH_methylation"), "CHH profile")
  x
}

if (!is.na(profile_file) && !ONSEN_FORCE_RESCAN) {
  chh_profile <- normalize_profile(read_table_auto(profile_file))
} else if (!is.na(selected_coord_file) && !is.na(methylome_file) && ONSEN_RUN_LARGE_STEPS) {
  if (!exists("mc")) mc <- read_methylome(methylome_file)
  selected_regions <- normalize_coordinates(read_table_auto(selected_coord_file))

  flank_bp <- 250L
  flank_bins <- 5L
  body_bins <- 20L

  create_bins_one <- function(locus_id, locus_class, chromosome, start, end) {
    width <- end - start + 1L
    make_segments <- function(segment_start, segment_end, n_bins, offset, region_part) {
      cuts <- floor(seq(segment_start, segment_end + 1, length.out = n_bins + 1L))
      data.frame(
        locus_id = locus_id,
        locus_class = locus_class,
        chromosome = chromosome,
        start = cuts[-length(cuts)],
        end = pmax(cuts[-1] - 1L, cuts[-length(cuts)]),
        bin_index = offset + seq_len(n_bins),
        region_part = region_part,
        stringsAsFactors = FALSE
      )
    }
    dplyr::bind_rows(
      make_segments(max(1L, start - flank_bp), start - 1L, flank_bins, 0L, "upstream"),
      make_segments(start, end, body_bins, flank_bins, "body"),
      make_segments(end + 1L, end + flank_bp, flank_bins, flank_bins + body_bins, "downstream")
    )
  }

  bins <- purrr::pmap_dfr(
    selected_regions[, c("locus_id", "locus_class", "chromosome", "start", "end")],
    create_bins_one
  )
  bin_dt <- data.table::as.data.table(bins)
  data.table::setkey(bin_dt, chromosome, start, end)
  chh_calls <- mc[methylation_context == "CHH", .(
    chromosome, start = position, end = position,
    methylated_bases, total_bases
  )]
  data.table::setkey(chh_calls, chromosome, start, end)
  matched_bins <- data.table::foverlaps(
    chh_calls, bin_dt,
    by.x = c("chromosome", "start", "end"),
    by.y = c("chromosome", "start", "end"),
    type = "within", nomatch = 0L
  )
  locus_profile <- matched_bins[, .(
    CHH_methylation = weighted_methylation(methylated_bases, total_bases)
  ), by = .(locus_id, locus_class, bin_index, region_part)]

  chh_profile <- locus_profile |>
    dplyr::group_by(locus_class, bin_index, region_part) |>
    dplyr::summarise(
      mean_CHH_methylation = mean(CHH_methylation, na.rm = TRUE),
      sem_CHH_methylation = sem(CHH_methylation),
      n_loci = sum(is.finite(CHH_methylation)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      normalized_position = dplyr::case_when(
        bin_index == 1L ~ "-250 bp",
        bin_index == flank_bins + 1L ~ "5' edge",
        bin_index == flank_bins + body_bins ~ "3' edge",
        bin_index == flank_bins + body_bins + flank_bins ~ "+250 bp",
        TRUE ~ ""
      )
    )
} else {
  chh_profile <- data.frame()
}

if (nrow(chh_profile)) {
  safe_write_csv(chh_profile, "Fig4D_aggregate_CHH_profile_repository.csv")
}

# -------------------------------- Figures ------------------------------------
if (ONSEN_MAKE_FIGURES) {
  plot_broad <- broad_locus |>
    dplyr::mutate(
      methylation_context = factor(methylation_context, levels = contexts),
      locus_class = factor(
        locus_class,
        levels = c("Ordinary background TE", "ONSEN LTR candidate")
      )
    )
  stat_label_col <- c("significance_label", "label")[
    c("significance_label", "label") %in% names(broad_stats)
  ][1]
  if (is.na(stat_label_col)) {
    p_col <- c("BH_adjusted_P", "BH_adjusted_p", "adjusted_P")[
      c("BH_adjusted_P", "BH_adjusted_p", "adjusted_P") %in% names(broad_stats)
    ][1]
    labels <- ifelse(broad_stats[[p_col]] <= 1e-4, "****",
      ifelse(broad_stats[[p_col]] <= 0.05, "*", "ns"))
  } else {
    labels <- as.character(broad_stats[[stat_label_col]])
  }
  annotation <- data.frame(
    methylation_context = factor(contexts, levels = contexts),
    y = 106, label = labels[match(contexts, broad_stats$methylation_context)]
  )

  p4a <- ggplot2::ggplot(
    plot_broad,
    ggplot2::aes(methylation_context, weighted_methylation_percent, fill = locus_class)
  ) +
    ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.62, outlier.shape = NA, linewidth = 0.6
    ) +
    ggplot2::geom_point(
      ggplot2::aes(group = locus_class),
      position = ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
      size = 1.15, alpha = 0.35
    ) +
    ggplot2::geom_text(
      data = annotation,
      ggplot2::aes(methylation_context, y, label = label),
      inherit.aes = FALSE, size = 4.4, fontface = "bold"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 110), clip = "on") +
    ggplot2::scale_fill_manual(values = c(
      "Ordinary background TE" = "#83AEE8",
      "ONSEN LTR candidate" = "#E989AE"
    )) +
    ggplot2::labs(
      x = "Methylation context",
      y = "Weighted methylation (%)",
      fill = NULL
    ) +
    theme_onsen(13) +
    ggplot2::theme(legend.position = "top")
  save_plot_pair(p4a, "Fig4A_ONSEN_vs_ordinary_TE_methylation", 6.8, 5.4)

  if (nrow(selected_locus)) {
    heat_data <- selected_locus |>
      dplyr::mutate(
        methylation_context = factor(methylation_context, levels = contexts),
        locus_id = factor(locus_id, levels = rev(unique(locus_id)))
      )
    p4b <- ggplot2::ggplot(
      heat_data,
      ggplot2::aes(methylation_context, locus_id, fill = weighted_methylation_percent)
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
      ggplot2::scale_fill_gradientn(
        colours = c("#F0F5F9", "#F7D18A", "#D8738A", "#7B2B83"),
        limits = c(0, 100), oob = scales::squish,
        name = "Weighted\nmethylation (%)"
      ) +
      ggplot2::labs(x = "Methylation context", y = NULL) +
      theme_onsen(11)
    save_plot_pair(p4b, "Fig4B_locus_methylation_heatmap", 6.0, 7.2)

    chh_scatter <- selected_locus |>
      dplyr::filter(methylation_context == "CHH", is.finite(HSF_hits_per_kb))
    if (nrow(chh_scatter)) {
      p4c <- ggplot2::ggplot(
        chh_scatter,
        ggplot2::aes(HSF_hits_per_kb, weighted_methylation_percent, colour = locus_class)
      ) +
        ggplot2::geom_point(size = 3, alpha = 0.85) +
        ggplot2::geom_text(
          ggplot2::aes(label = locus_id),
          size = 2.5, check_overlap = TRUE, vjust = -0.6
        ) +
        ggplot2::labs(
          x = "HSF motif density (hits per kb)",
          y = "CHH methylation (%)",
          colour = NULL
        ) +
        theme_onsen(12) +
        ggplot2::theme(legend.position = "top")
      save_plot_pair(p4c, "Fig4C_HSF_density_vs_CHH", 6.5, 5.2)
    }
  }

  if (nrow(chh_profile)) {
    p4d <- ggplot2::ggplot(
      chh_profile,
      ggplot2::aes(bin_index, mean_CHH_methylation, colour = locus_class, fill = locus_class)
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = pmax(0, mean_CHH_methylation - sem_CHH_methylation),
          ymax = pmin(100, mean_CHH_methylation + sem_CHH_methylation)
        ),
        alpha = 0.18, colour = NA
      ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point(size = 1.8) +
      ggplot2::labs(
        x = "Normalized position",
        y = "CHH methylation (%)",
        colour = NULL, fill = NULL
      ) +
      theme_onsen(12) +
      ggplot2::theme(legend.position = "top")
    save_plot_pair(p4d, "Fig4D_aggregate_CHH_profile", 7.2, 4.8)
  }
}

message("Basal Col-0 leaf methylation analysis completed.")
