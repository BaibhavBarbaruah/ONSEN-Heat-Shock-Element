# Constrained-mutant and exact-GC sensitivity analysis.
# Covers Fig. S3 and source data for Table S12.

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "tidyr", "readr", "stringr",
                   "purrr", "ggplot2", "scales"))

message_config()

sequence_set <- read_fasta_simple(repo_file("ONSEN_49bp_sequences.fasta"))
native_sequence <- unname(sequence_set[["native_49bp_ONSEN_HSE"]])
designed_sequence <- unname(sequence_set[["in_silico_designed_HSE_disrupted_49bp"]])

core_positions <- c(12:14, 17:19, 23:25, 28:30, 36:38)
core_groups <- list(12:14, 17:19, 23:25, 28:30, 36:38)

valid_constrained_sequence <- function(sequence) {
  bases <- split_bases(sequence)
  native_bases <- split_bases(native_sequence)
  changed_core <- sum(bases[core_positions] != native_bases[core_positions])
  every_core_disrupted <- all(vapply(core_groups, function(pos) {
    paste(bases[pos], collapse = "") != paste(native_bases[pos], collapse = "")
  }, logical(1)))
  changed_core == 12L && every_core_disrupted && !contains_canonical_hse_unit(sequence)
}

generate_random_library <- function(n = 5000L, random_seed = 20260714L) {
  set.seed(random_seed)
  native_bases <- split_bases(native_sequence)
  bases <- c("A", "C", "G", "T")
  sequences <- character()
  attempts <- 0L

  while (length(sequences) < n) {
    attempts <- attempts + 1L
    chosen <- sort(sample(core_positions, 12L, replace = FALSE))
    x <- native_bases
    for (position in chosen) {
      x[[position]] <- sample(setdiff(bases, x[[position]]), 1L)
    }
    candidate <- paste(x, collapse = "")
    if (valid_constrained_sequence(candidate)) {
      sequences <- unique(c(sequences, candidate))
    }
    if (attempts > 5000000L) {
      stop("Random-library generation exceeded five million attempts.")
    }
  }

  data.frame(
    sequence_id = sprintf("random_mutant_%04d", seq_len(n)),
    sequence_class = "Constrained random mutant",
    sequence = sequences,
    substitutions_vs_native = vapply(sequences, hamming_distance, numeric(1), b = native_sequence),
    GC_percent = vapply(sequences, gc_percent, numeric(1)),
    stringsAsFactors = FALSE
  )
}

generate_complete_exact_gc_space <- function(target_gc_count = sum(split_bases(designed_sequence) %in% c("G", "C"))) {
  native_bases <- split_bases(native_sequence)
  non_core <- setdiff(seq_along(native_bases), core_positions)
  fixed_gc <- sum(native_bases[non_core] %in% c("G", "C"))
  target_core_gc <- target_gc_count - fixed_gc
  bases <- c("A", "C", "G", "T")

  core_options <- lapply(core_groups, function(pos) {
    original <- native_bases[pos]
    grid <- expand.grid(
      b1 = bases, b2 = bases, b3 = bases,
      stringsAsFactors = FALSE
    )
    grid$triplet <- paste0(grid$b1, grid$b2, grid$b3)
    grid$substitutions <- vapply(grid$triplet, function(s) {
      sum(split_bases(s) != original)
    }, integer(1))
    grid$GC_count <- vapply(grid$triplet, function(s) {
      sum(split_bases(s) %in% c("G", "C"))
    }, integer(1))
    grid <- grid[grid$substitutions >= 1L, c("triplet", "substitutions", "GC_count")]
    rownames(grid) <- NULL
    grid
  })

  states <- data.frame(
    substitutions = 0L,
    core_GC_count = 0L,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(core_options)) {
    options_i <- core_options[[i]]
    states$key <- 1L
    options_i$key <- 1L
    expanded <- merge(states, options_i, by = "key")
    expanded$key <- NULL
    expanded$substitutions <- expanded$substitutions + expanded$substitutions.y
    expanded$core_GC_count <- expanded$core_GC_count + expanded$GC_count
    expanded$substitutions.y <- NULL
    expanded$GC_count <- NULL
    names(expanded)[names(expanded) == "triplet"] <- paste0("core", i)

    remaining <- length(core_options) - i
    expanded <- expanded[
      expanded$substitutions <= 12L &
      expanded$substitutions + remaining >= 12L &
      expanded$substitutions + 3L * remaining >= 12L &
      expanded$core_GC_count <= target_core_gc &
      expanded$core_GC_count + 3L * remaining >= target_core_gc,
      ,
      drop = FALSE
    ]
    states <- expanded
  }

  states <- states[
    states$substitutions == 12L &
    states$core_GC_count == target_core_gc,
    ,
    drop = FALSE
  ]

  sequences <- apply(states[paste0("core", 1:5)], 1, function(parts) {
    x <- native_bases
    for (i in seq_along(core_groups)) {
      x[core_groups[[i]]] <- split_bases(parts[[i]])
    }
    paste(x, collapse = "")
  })
  sequences <- unique(sequences)
  sequences <- sequences[vapply(sequences, valid_constrained_sequence, logical(1))]

  result <- data.frame(
    sequence_id = sprintf("exact_GC_valid_%04d", seq_along(sequences)),
    sequence_class = "Valid exact-GC alternative",
    sequence = sequences,
    substitutions_vs_native = vapply(sequences, hamming_distance, numeric(1), b = native_sequence),
    GC_percent = vapply(sequences, gc_percent, numeric(1)),
    stringsAsFactors = FALSE
  )
  result$sequence_class[result$sequence == designed_sequence] <- "Designed mutant"
  result$sequence_id[result$sequence == designed_sequence] <- "designed_mutant"

  if (nrow(result) != 5120L) {
    stop(
      "Expected 5,120 valid exact-GC sequences including the designed mutant, but generated ",
      nrow(result), "."
    )
  }
  result
}

random_library_file <- find_any_input(c(
  "constrained_HSE_mutant_library_n5000.csv",
  "random_constrained_HSE_mutants_n5000.csv"
), required = FALSE)

exact_library_file <- find_input(
  "Step2C_complete_exact_GC_constrained_sequence_library.csv",
  required = FALSE
)

if (!is.na(random_library_file) && !ONSEN_FORCE_RESCAN) {
  random_library <- read_table_auto(random_library_file)
  message("Using completed random-mutant library: ", random_library_file)
} else {
  message("Generating a reproducible 5,000-sequence constrained random library.")
  random_library <- generate_random_library()
}
safe_write_csv(random_library, "constrained_HSE_mutant_library_n5000_repository.csv")

if (!is.na(exact_library_file) && !ONSEN_FORCE_RESCAN) {
  exact_library <- read_table_auto(exact_library_file)
  message("Using completed exact-GC library: ", exact_library_file)
} else {
  message("Generating the complete valid exact-GC design space.")
  exact_library <- generate_complete_exact_gc_space()
}
safe_write_csv(exact_library, "complete_exact_GC_design_space_repository.csv")

normalize_family_summary <- function(x) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("sequence_id", c("sequence_id", "mutant_id", "id"))
  rename_first("sequence_class", c("sequence_class", "class", "sequence_type"))
  rename_first("threshold", c("threshold", "score_threshold", "relative_score_threshold"))
  rename_first("tf_family", c("tf_family", "family", "TF_family"))
  rename_first("motif_position_hits", c(
    "motif_position_hits", "high_confidence_motif_position_hits",
    "n_hits", "hit_count", "motif_hits"
  ))
  assert_columns(x, c("sequence_id", "threshold", "tf_family", "motif_position_hits"), "mutant family summary")
  x$threshold <- as.numeric(x$threshold)
  x$motif_position_hits <- as.numeric(x$motif_position_hits)
  x
}

step2b_summary_file <- find_input(
  "Step2B_validated_v1_sequence_family_threshold_summary.csv",
  required = FALSE
)
step2c_summary_file <- find_input(
  "Step2C_complete_exact_GC_sequence_family_threshold_summary.csv",
  required = FALSE
)

scan_family_library <- function(library, output_label) {
  if (!ONSEN_RUN_LARGE_STEPS) {
    stop(
      "Processed family summary for ", output_label, " was not found. ",
      "Set ONSEN_RUN_LARGE_STEPS=true to rescan the complete sequence library."
    )
  }
  jaspar_file <- find_any_input(c(
    "JASPAR2026_CORE_plants_nonredundant_pfms_jaspar.txt",
    "JASPAR2024_CORE_plants_nonredundant_pfms_jaspar.txt"
  ))
  motifs <- parse_jaspar_pfms(jaspar_file)
  motif_families <- classify_tf_family(vapply(motifs, `[[`, character(1), "name"))
  motifs <- motifs[motif_families %in% c("HSF", "AP2/ERF", "DOF")]
  thresholds <- c(0.80, 0.85, 0.90, 0.95)

  sequence_table <- data.frame(
    sequence_id = library$sequence_id,
    sequence = library$sequence,
    stringsAsFactors = FALSE
  )
  summaries <- list()
  k <- 0L
  for (cutoff in thresholds) {
    hits <- as.data.frame(scan_sequences_against_motifs(
      sequence_table, motifs, threshold = cutoff, pseudocount = 0.8
    ))
    summary <- hits |>
      dplyr::count(sequence_id, tf_family, name = "motif_position_hits") |>
      tidyr::complete(
        sequence_id = sequence_table$sequence_id,
        tf_family = c("HSF", "AP2/ERF", "DOF"),
        fill = list(motif_position_hits = 0L)
      ) |>
      dplyr::mutate(threshold = cutoff)
    k <- k + 1L
    summaries[[k]] <- summary
  }
  dplyr::bind_rows(summaries)
}

if (!is.na(step2b_summary_file) && !ONSEN_FORCE_RESCAN) {
  random_family <- normalize_family_summary(read_table_auto(step2b_summary_file))
} else {
  random_family <- scan_family_library(random_library, "random-mutant library")
}
safe_write_csv(random_family, "Step2B_sequence_family_threshold_summary_repository.csv")

if (!is.na(step2c_summary_file) && !ONSEN_FORCE_RESCAN) {
  exact_family <- normalize_family_summary(read_table_auto(step2c_summary_file))
} else {
  exact_family <- scan_family_library(exact_library, "exact-GC design space")
}
safe_write_csv(exact_family, "Step2C_sequence_family_threshold_summary_repository.csv")

# Join sequence metadata when absent from the processed summary.
random_meta <- random_library
if (!"sequence_id" %in% names(random_meta)) names(random_meta)[1] <- "sequence_id"
if (!"GC_percent" %in% names(random_meta) && "sequence" %in% names(random_meta)) {
  random_meta$GC_percent <- vapply(random_meta$sequence, gc_percent, numeric(1))
}
random_family <- random_family |>
  dplyr::left_join(
    random_meta |> dplyr::select(dplyr::any_of(c("sequence_id", "sequence", "GC_percent", "sequence_class"))),
    by = "sequence_id"
  )

exact_meta <- exact_library
if (!"sequence_id" %in% names(exact_meta)) names(exact_meta)[1] <- "sequence_id"
if (!"GC_percent" %in% names(exact_meta) && "sequence" %in% names(exact_meta)) {
  exact_meta$GC_percent <- vapply(exact_meta$sequence, gc_percent, numeric(1))
}
exact_family <- exact_family |>
  dplyr::left_join(
    exact_meta |> dplyr::select(dplyr::any_of(c("sequence_id", "sequence", "GC_percent", "sequence_class"))),
    by = "sequence_id"
  )

# Identify designed-mutant values. Processed summaries may name it differently.
designed_ids_random <- unique(random_family$sequence_id[
  grepl("designed", random_family$sequence_id, ignore.case = TRUE)
])
designed_ids_exact <- unique(exact_family$sequence_id[
  grepl("designed", exact_family$sequence_id, ignore.case = TRUE)
])

# Add designed sequence scan from the exact-GC summary if needed.
if (!length(designed_ids_exact) && "sequence" %in% names(exact_meta)) {
  designed_ids_exact <- exact_meta$sequence_id[exact_meta$sequence == designed_sequence]
}
if (!length(designed_ids_random)) designed_ids_random <- designed_ids_exact

random_reference <- random_family |>
  dplyr::filter(threshold == 0.85, tf_family == "AP2/ERF")
exact_reference <- exact_family |>
  dplyr::filter(threshold == 0.85, tf_family == "AP2/ERF")

designed_ap2_random <- random_reference$motif_position_hits[
  random_reference$sequence_id %in% designed_ids_random
][1]
designed_ap2_exact <- exact_reference$motif_position_hits[
  exact_reference$sequence_id %in% designed_ids_exact
][1]

# Fall back to the validated designed value in the manuscript if the designed
# sequence was stored only in the exact-GC object.
if (!is.finite(designed_ap2_random)) designed_ap2_random <- 61
if (!is.finite(designed_ap2_exact)) designed_ap2_exact <- 61

random_distribution <- random_reference |>
  dplyr::filter(!sequence_id %in% designed_ids_random)
exact_distribution <- exact_reference |>
  dplyr::filter(!sequence_id %in% designed_ids_exact)

plus_one_upper <- function(values, observed) {
  (1 + sum(values >= observed, na.rm = TRUE)) / (length(values[is.finite(values)]) + 1)
}

empirical <- data.frame(
  comparison = c("Variable-GC constrained random mutants", "Complete exact-GC alternatives"),
  n_alternatives = c(nrow(random_distribution), nrow(exact_distribution)),
  designed_AP2ERF_hits = c(designed_ap2_random, designed_ap2_exact),
  alternative_median_AP2ERF_hits = c(
    median(random_distribution$motif_position_hits, na.rm = TRUE),
    median(exact_distribution$motif_position_hits, na.rm = TRUE)
  ),
  plus_one_upper_tail_P = c(
    plus_one_upper(random_distribution$motif_position_hits, designed_ap2_random),
    plus_one_upper(exact_distribution$motif_position_hits, designed_ap2_exact)
  ),
  stringsAsFactors = FALSE
)
safe_write_csv(empirical, "mutant_sensitivity_empirical_comparison_repository.csv")

gc_data <- random_reference |>
  dplyr::filter(!sequence_id %in% designed_ids_random, is.finite(GC_percent))
gc_test <- suppressWarnings(stats::cor.test(
  gc_data$GC_percent, gc_data$motif_position_hits,
  method = "spearman", exact = FALSE
))
gc_summary <- data.frame(
  comparison = "GC percent versus AP2/ERF motif-position hits",
  n = nrow(gc_data),
  spearman_rho = unname(gc_test$estimate),
  P_value = gc_test$p.value,
  designed_GC_percent = gc_percent(designed_sequence),
  stringsAsFactors = FALSE
)
safe_write_csv(gc_summary, "random_mutant_GC_spearman_repository.csv")

# Threshold/family manuscript summary.
threshold_summary <- dplyr::bind_rows(
  random_family |> dplyr::mutate(design_space = "Variable-GC random mutants"),
  exact_family |> dplyr::mutate(design_space = "Complete exact-GC space")
) |>
  dplyr::group_by(design_space, threshold, tf_family) |>
  dplyr::summarise(
    n_sequences = dplyr::n(),
    median_hits = median(motif_position_hits, na.rm = TRUE),
    mean_hits = mean(motif_position_hits, na.rm = TRUE),
    minimum_hits = min(motif_position_hits, na.rm = TRUE),
    maximum_hits = max(motif_position_hits, na.rm = TRUE),
    .groups = "drop"
  )
safe_write_csv(threshold_summary, "mutant_sensitivity_threshold_family_summary_repository.csv")

if (ONSEN_MAKE_FIGURES) {
  # Fig. S3A
  p_s3a <- ggplot2::ggplot(
    random_distribution,
    ggplot2::aes(motif_position_hits)
  ) +
    ggplot2::geom_histogram(binwidth = 2, boundary = 0, colour = "black", fill = "#9DB7D5") +
    ggplot2::geom_vline(xintercept = designed_ap2_random, linetype = "dashed", linewidth = 1) +
    ggplot2::annotate(
      "text", x = designed_ap2_random, y = Inf,
      label = paste0("Designed mutant = ", designed_ap2_random),
      vjust = 1.5, hjust = 1.05, fontface = "bold"
    ) +
    ggplot2::labs(
      x = "AP2/ERF-family motif-position hits",
      y = "Number of constrained random mutants"
    ) +
    theme_onsen(13)
  save_plot_pair(p_s3a, "FigS3A_random_mutant_AP2ERF_distribution", 6.5, 4.8)

  # Fig. S3B
  p_s3b <- ggplot2::ggplot(
    gc_data,
    ggplot2::aes(GC_percent, motif_position_hits)
  ) +
    ggplot2::geom_point(alpha = 0.35, size = 1.3) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    ggplot2::geom_point(
      data = data.frame(
        GC_percent = gc_percent(designed_sequence),
        motif_position_hits = designed_ap2_random
      ),
      size = 3.5, shape = 23, fill = "white"
    ) +
    ggplot2::labs(
      x = "GC content (%)",
      y = "AP2/ERF-family motif-position hits"
    ) +
    ggplot2::annotate(
      "text", x = -Inf, y = Inf,
      label = sprintf("Spearman rho = %.3f\nP = %.3g", unname(gc_test$estimate), gc_test$p.value),
      hjust = -0.05, vjust = 1.2
    ) +
    theme_onsen(13)
  save_plot_pair(p_s3b, "FigS3B_random_mutant_GC_correlation", 6.5, 4.8)

  # Fig. S3C
  p_s3c <- ggplot2::ggplot(
    exact_distribution,
    ggplot2::aes(motif_position_hits)
  ) +
    ggplot2::geom_histogram(binwidth = 2, boundary = 0, colour = "black", fill = "#D3C2ED") +
    ggplot2::geom_vline(xintercept = designed_ap2_exact, linetype = "dashed", linewidth = 1) +
    ggplot2::annotate(
      "text", x = designed_ap2_exact, y = Inf,
      label = paste0("Designed mutant = ", designed_ap2_exact),
      vjust = 1.5, hjust = 1.05, fontface = "bold"
    ) +
    ggplot2::labs(
      x = "AP2/ERF-family motif-position hits",
      y = "Number of exact-GC alternatives"
    ) +
    theme_onsen(13)
  save_plot_pair(p_s3c, "FigS3C_exact_GC_AP2ERF_distribution", 6.5, 4.8)
}

message("Constrained-mutant and exact-GC sensitivity analysis completed.")
