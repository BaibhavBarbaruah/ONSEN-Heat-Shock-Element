# Natural-accession ONSEN-like HSE analysis.
# Covers Fig. 6, Fig. 7, Fig. S1, Fig. S2 and source data for Tables S8-S10.

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "tidyr", "readr", "stringr",
                   "purrr", "tibble", "ggplot2", "forcats", "scales",
                   "Biostrings"))

message_config()

sequence_set <- read_fasta_simple(repo_file("ONSEN_49bp_sequences.fasta"))
seed <- unname(sequence_set[["native_49bp_ONSEN_HSE"]])
seed_rc <- rev_comp(seed)
seed_width <- nchar(seed)
max_mismatches <- 4L
flank_bp <- 400L

main_chromosome <- function(name) {
  name <- as.character(name)
  grepl(
    "(^Chr[1-5]$)|(^[1-5]$)|chromosome[_ -]?[1-5]$|(^|[._-])Chr?[1-5]($|[._-])",
    name,
    ignore.case = TRUE
  )
}

# -------------------------- Optional raw genome scan --------------------------
scan_accession_genomes <- function(genome_files) {
  all_hits <- list()
  hit_index <- 0L

  for (file in genome_files) {
    accession <- basename(file)
    accession <- sub("\\.gz$", "", accession, ignore.case = TRUE)
    accession <- sub("\\.(fa|fasta|fna)$", "", accession, ignore.case = TRUE)
    accession <- sub("\\.chr\\.all.*$", "", accession, ignore.case = TRUE)

    genome <- Biostrings::readDNAStringSet(file)
    for (i in seq_along(genome)) {
      chromosome <- names(genome)[[i]]
      chromosome_sequence <- genome[[i]]
      chromosome_width <- Biostrings::width(chromosome_sequence)

      plus <- Biostrings::matchPattern(
        seed, chromosome_sequence,
        max.mismatch = max_mismatches,
        with.indels = FALSE, fixed = TRUE
      )
      minus <- Biostrings::matchPattern(
        seed_rc, chromosome_sequence,
        max.mismatch = max_mismatches,
        with.indels = FALSE, fixed = TRUE
      )

      append_hits <- function(matches, strand) {
        if (!length(matches)) return()
        for (j in seq_along(matches)) {
          start <- BiocGenerics::start(matches)[[j]]
          end <- BiocGenerics::end(matches)[[j]]
          raw_sequence <- as.character(Biostrings::subseq(
            chromosome_sequence, start = start, end = end
          ))
          oriented <- if (strand == "+") raw_sequence else rev_comp(raw_sequence)
          window_start <- max(1L, start - flank_bp)
          window_end <- min(chromosome_width, end + flank_bp)
          raw_window <- as.character(Biostrings::subseq(
            chromosome_sequence, start = window_start, end = window_end
          ))
          oriented_window <- if (strand == "+") raw_window else rev_comp(raw_window)

          hit_index <<- hit_index + 1L
          all_hits[[hit_index]] <<- data.frame(
            accession = accession,
            chromosome = chromosome,
            main_chromosome = main_chromosome(chromosome),
            seed_start = start,
            seed_end = end,
            strand = strand,
            mismatches_to_Col0 = hamming_distance(oriented, seed, ignore_n = TRUE),
            seed_sequence = oriented,
            window_start = window_start,
            window_end = window_end,
            candidate_window_sequence = oriented_window,
            stringsAsFactors = FALSE
          )
        }
      }

      append_hits(plus, "+")
      append_hits(minus, "-")
    }
  }

  hits <- dplyr::bind_rows(all_hits) |>
    dplyr::arrange(accession, chromosome, seed_start, strand) |>
    dplyr::group_by(accession) |>
    dplyr::mutate(
      candidate_id = paste0(
        accession, "_HSEwin_",
        stringr::str_pad(dplyr::row_number(), 3, pad = "0")
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(candidate_id, dplyr::everything())
  hits
}

candidate_metadata_file <- find_any_input(c(
  "accession_ONSEN_like_HSE_candidate_windows_mismatch_leq4_metadata.csv",
  "accession_ONSEN_like_49bp_HSE_seed_hits_mismatch_leq4.csv",
  "accession_ONSEN_like_49bp_HSE_seed_hits.csv"
), required = FALSE)

if (!is.na(candidate_metadata_file) && !ONSEN_FORCE_RESCAN) {
  candidates <- read_table_auto(candidate_metadata_file)
  message("Using completed accession candidate metadata: ", candidate_metadata_file)
} else {
  if (!ONSEN_RUN_LARGE_STEPS) {
    stop(
      "Accession candidate metadata was not found and raw genome scanning is disabled.\n",
      "Set ONSEN_RUN_LARGE_STEPS=true to scan accession assemblies."
    )
  }
  genome_files <- list.files(
    ONSEN_DATA_ROOT,
    pattern = "\\.(fa|fasta|fna)(\\.gz)?$",
    recursive = TRUE, full.names = TRUE, ignore.case = TRUE
  )
  genome_files <- genome_files[
    grepl("An-1|C24|Cvi|Eri|Kyo|Ler|Sha", basename(genome_files), ignore.case = TRUE)
  ]
  if (!length(genome_files)) stop("No accession genome FASTA files were found.")
  candidates <- scan_accession_genomes(genome_files)
}

normalize_candidates <- function(x) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("candidate_id", c("candidate_id", "sequence_id", "hit_id"))
  rename_first("accession", c("accession", "Accession"))
  rename_first("chromosome", c("chromosome", "chr", "seqnames", "Chromosome"))
  rename_first("seed_start", c("seed_start", "start", "hit_start"))
  rename_first("seed_end", c("seed_end", "end", "hit_end"))
  rename_first("strand", c("strand"))
  rename_first("mismatches_to_Col0", c(
    "mismatches_to_Col0", "mismatches_to_col0_seed", "mismatch_count", "mismatches"
  ))
  rename_first("seed_sequence", c(
    "seed_sequence", "candidate_sequence", "49bp_sequence", "sequence_49bp"
  ))
  rename_first("candidate_window_sequence", c(
    "candidate_window_sequence", "window_sequence", "sequence", "plusminus400_sequence"
  ))
  rename_first("main_chromosome", c(
    "main_chromosome", "is_main_chromosome", "mainchr"
  ))
  assert_columns(x, c("accession", "chromosome"), "accession candidate table")
  if (!"candidate_id" %in% names(x)) {
    x$candidate_id <- paste0(x$accession, "_candidate_", seq_len(nrow(x)))
  }
  if (!"main_chromosome" %in% names(x)) {
    x$main_chromosome <- main_chromosome(x$chromosome)
  }
  x$main_chromosome <- as.logical(x$main_chromosome)
  if (!"mismatches_to_Col0" %in% names(x) && "seed_sequence" %in% names(x)) {
    x$mismatches_to_Col0 <- vapply(
      x$seed_sequence, hamming_distance, numeric(1), b = seed, ignore_n = TRUE
    )
  }
  if (!"seed_sequence" %in% names(x) && "candidate_window_sequence" %in% names(x)) {
    # The 49-bp seed is centered in an 849-bp +/-400 window.
    x$seed_sequence <- substring(
      x$candidate_window_sequence,
      flank_bp + 1L,
      flank_bp + seed_width
    )
  }
  x
}
candidates <- normalize_candidates(candidates)
main_candidates <- candidates[candidates$main_chromosome %in% TRUE, , drop = FALSE]

safe_write_csv(candidates, "accession_ONSEN_like_candidates_all_repository.csv")
safe_write_csv(main_candidates, "accession_ONSEN_like_candidates_main_chromosomes_repository.csv")

# HSE architecture within the 49-bp seed.
if ("seed_sequence" %in% names(main_candidates)) {
  architecture <- main_candidates |>
    dplyr::rowwise() |>
    dplyr::mutate(
      canonical_HSE_units = nrow(find_canonical_hse_units(seed_sequence)),
      HSE_architecture_class = dplyr::case_when(
        canonical_HSE_units >= 5L ~ "Canonical five-unit",
        canonical_HSE_units == 4L ~ "Partial four-unit",
        canonical_HSE_units == 3L ~ "Weak three-unit",
        TRUE ~ "Low HSE architecture"
      ),
      seed_identity_percent = 100 * (seed_width - mismatches_to_Col0) / seed_width
    ) |>
    dplyr::ungroup()
} else {
  architecture_file <- find_input(
    "accession_ONSEN_like_HSE_candidate_windows_mismatch_leq4_HSE_architecture_scores.csv"
  )
  architecture <- read_table_auto(architecture_file)
}
safe_write_csv(architecture, "accession_HSE_architecture_repository.csv")

# -------------------------- HSF scan and summaries ----------------------------
hsf_candidate_file <- find_input(
  "FIXED_accession_ONSEN_like_mainchr_candidate_windows_JASPAR2026_Arabidopsis_HSF_candidate_summary_threshold_0.85.csv",
  required = FALSE
)
hsf_accession_file <- find_input(
  "FIXED_accession_ONSEN_like_mainchr_candidate_windows_JASPAR2026_Arabidopsis_HSF_accession_summary_threshold_0.85.csv",
  required = FALSE
)
hsf_model_summary_file <- find_input(
  "FIXED_accession_by_HSF_motif_model_summary_threshold_0.85.csv",
  required = FALSE
)

normalize_hsf_candidate <- function(x) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("candidate_id", c("candidate_id", "sequence_id"))
  rename_first("accession", c("accession", "Accession"))
  rename_first("HSF_hits", c("HSF_hits", "hsf_hits", "motif_position_hits"))
  rename_first("HSF_hits_per_kb", c("HSF_hits_per_kb", "hsf_density", "motif_density"))
  rename_first("unique_HSF_models", c("unique_HSF_models", "unique_models"))
  assert_columns(x, c("candidate_id", "accession"), "accession HSF summary")
  x
}

if (!is.na(hsf_candidate_file) && !ONSEN_FORCE_RESCAN) {
  hsf_candidates <- normalize_hsf_candidate(read_table_auto(hsf_candidate_file))
} else {
  if (!ONSEN_RUN_LARGE_STEPS) {
    stop(
      "Fixed accession HSF summary was not found and large motif scanning is disabled."
    )
  }
  jaspar_file <- find_input("JASPAR2026_CORE_plants_nonredundant_pfms_jaspar.txt")
  motifs <- parse_jaspar_pfms(jaspar_file)
  hsf_names <- read.csv(repo_file("Arabidopsis_HSF_models_JASPAR2026.csv"))$HSF_model
  motif_names <- toupper(vapply(motifs, `[[`, character(1), "name"))
  hsf_motifs <- motifs[motif_names %in% toupper(hsf_names)]

  seq_col <- if ("candidate_window_sequence" %in% names(main_candidates)) {
    "candidate_window_sequence"
  } else {
    "seed_sequence"
  }
  sequence_table <- main_candidates |>
    dplyr::transmute(
      sequence_id = candidate_id,
      sequence = .data[[seq_col]]
    )
  hits <- as.data.frame(scan_sequences_against_motifs(
    sequence_table, hsf_motifs, threshold = 0.85, pseudocount = 0.8
  ))
  hsf_candidates <- hits |>
    dplyr::group_by(sequence_id) |>
    dplyr::summarise(
      HSF_hits = dplyr::n(),
      unique_HSF_models = dplyr::n_distinct(motif_id),
      maximum_relative_score = max(relative_score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::right_join(
      main_candidates |>
        dplyr::transmute(
          sequence_id = candidate_id,
          accession,
          width_bp = nchar(.data[[seq_col]])
        ),
      by = "sequence_id"
    ) |>
    dplyr::mutate(
      HSF_hits = tidyr::replace_na(HSF_hits, 0L),
      unique_HSF_models = tidyr::replace_na(unique_HSF_models, 0L),
      HSF_hits_per_kb = HSF_hits / width_bp * 1000
    ) |>
    dplyr::rename(candidate_id = sequence_id)
}
safe_write_csv(hsf_candidates, "accession_HSF_candidate_summary_repository.csv")

if (!is.na(hsf_accession_file) && !ONSEN_FORCE_RESCAN) {
  hsf_accession <- read_table_auto(hsf_accession_file)
} else {
  hsf_accession <- hsf_candidates |>
    dplyr::group_by(accession) |>
    dplyr::summarise(
      candidate_windows = dplyr::n(),
      median_HSF_hits_per_kb = median(HSF_hits_per_kb, na.rm = TRUE),
      mean_HSF_hits_per_kb = mean(HSF_hits_per_kb, na.rm = TRUE),
      .groups = "drop"
    )
}
safe_write_csv(hsf_accession, "accession_HSF_accession_summary_repository.csv")

# ----------------------- Paired HSE/LTR structural proxies -------------------
proxy_file <- find_input(
  "putative_ONSEN_like_copy_proxy_summary_mismatch_leq4.csv",
  required = FALSE
)

if (!is.na(proxy_file) && !ONSEN_FORCE_RESCAN) {
  proxy_summary <- read_table_auto(proxy_file)
} else {
  pair_min <- as.integer(Sys.getenv("ONSEN_PAIR_MIN_BP", "4000"))
  pair_max <- as.integer(Sys.getenv("ONSEN_PAIR_MAX_BP", "6000"))

  paired_ids <- character()
  pair_rows <- list()
  pair_index <- 0L

  for (accession in unique(main_candidates$accession)) {
    accession_data <- main_candidates[main_candidates$accession == accession, , drop = FALSE]
    for (chromosome in unique(accession_data$chromosome)) {
      z <- accession_data[accession_data$chromosome == chromosome, , drop = FALSE]
      z <- z[order(z$seed_start), , drop = FALSE]
      used <- rep(FALSE, nrow(z))
      if (nrow(z) >= 2L) {
        for (i in seq_len(nrow(z) - 1L)) {
          if (used[[i]]) next
          distances <- z$seed_start[(i + 1L):nrow(z)] - z$seed_start[[i]]
          possible <- which(distances >= pair_min & distances <= pair_max)
          if (length(possible)) {
            j <- i + possible[[1]]
            if (!used[[j]]) {
              pair_index <- pair_index + 1L
              pair_rows[[pair_index]] <- data.frame(
                accession = accession,
                chromosome = chromosome,
                left_candidate = z$candidate_id[[i]],
                right_candidate = z$candidate_id[[j]],
                HSE_to_HSE_distance_bp = z$seed_start[[j]] - z$seed_start[[i]],
                stringsAsFactors = FALSE
              )
              used[c(i, j)] <- TRUE
              paired_ids <- c(paired_ids, z$candidate_id[c(i, j)])
            }
          }
        }
      }
    }
  }

  paired_table <- dplyr::bind_rows(pair_rows)
  safe_write_csv(paired_table, "putative_ONSEN_like_paired_HSE_candidates_repository.csv")

  proxy_summary <- main_candidates |>
    dplyr::group_by(accession) |>
    dplyr::summarise(
      candidate_windows = dplyr::n(),
      paired_HSE_LTR_copy_proxy = sum(candidate_id %in% paired_ids) / 2,
      unpaired_HSE_candidate = sum(!candidate_id %in% paired_ids),
      .groups = "drop"
    )
}
safe_write_csv(proxy_summary, "putative_ONSEN_like_copy_proxy_summary_repository.csv")

# --------------------------- Accession summary -------------------------------
candidate_summary <- architecture |>
  dplyr::group_by(accession) |>
  dplyr::summarise(
    candidate_windows = dplyr::n(),
    exact_seed_matches = sum(mismatches_to_Col0 == 0, na.rm = TRUE),
    one_to_two_mismatch_seeds = sum(mismatches_to_Col0 %in% 1:2, na.rm = TRUE),
    three_to_four_mismatch_seeds = sum(mismatches_to_Col0 %in% 3:4, na.rm = TRUE),
    median_seed_identity = median(seed_identity_percent, na.rm = TRUE),
    median_seed_HSE_score = median(canonical_HSE_units, na.rm = TRUE),
    maximum_seed_HSE_score = max(canonical_HSE_units, na.rm = TRUE),
    .groups = "drop"
  )
safe_write_csv(candidate_summary, "accession_candidate_architecture_summary_repository.csv")

# ---------------------------- Natural variants -------------------------------
variant_accession_file <- find_input(
  "natural_accession_49bp_HSE_seed_variant_accession_summary.csv",
  required = FALSE
)
variant_summary_file <- find_input(
  "natural_accession_49bp_HSE_seed_variant_summary.csv",
  required = FALSE
)
family_long_file <- find_input(
  "natural_accession_variant_seed_TF_family_summary_long.csv",
  required = FALSE
)
family_delta_file <- find_input(
  "natural_accession_variant_seed_TF_family_delta_vs_Col0_long.csv",
  required = FALSE
)

if (!is.na(variant_summary_file) && !ONSEN_FORCE_RESCAN) {
  variant_summary <- read_table_auto(variant_summary_file)
} else {
  assert_columns(architecture, c("accession", "seed_sequence"), "accession architecture table")
  unique_variants <- architecture |>
    dplyr::filter(accession != "Col-0", seed_sequence != seed) |>
    dplyr::distinct(accession, seed_sequence, .keep_all = TRUE) |>
    dplyr::arrange(accession, mismatches_to_Col0, seed_sequence) |>
    dplyr::group_by(accession) |>
    dplyr::mutate(
      variant_id = paste0(
        accession, "-v",
        stringr::str_pad(dplyr::row_number(), 2, pad = "0")
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::rowwise() |>
    dplyr::mutate(
      substitution_positions = paste(
        which(split_bases(seed_sequence) != split_bases(seed)),
        collapse = ","
      ),
      substitution_count = hamming_distance(seed_sequence, seed),
      sequence_identity_percent = 100 * (seed_width - substitution_count) / seed_width
    ) |>
    dplyr::ungroup()
  variant_summary <- unique_variants
}
safe_write_csv(variant_summary, "natural_accession_49bp_variant_summary_repository.csv")

if (!is.na(variant_accession_file) && !ONSEN_FORCE_RESCAN) {
  variant_accession <- read_table_auto(variant_accession_file)
} else {
  variant_accession <- variant_summary |>
    dplyr::group_by(accession) |>
    dplyr::summarise(
      unique_variant_sequences = dplyr::n(),
      median_substitution_count = median(substitution_count, na.rm = TRUE),
      minimum_sequence_identity_percent = min(sequence_identity_percent, na.rm = TRUE),
      .groups = "drop"
    )
}
safe_write_csv(variant_accession, "natural_accession_49bp_variant_accession_summary_repository.csv")

normalize_family_long <- function(x) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("variant_id", c("variant_id", "sequence_id", "variant"))
  rename_first("accession", c("accession", "Accession"))
  rename_first("tf_family", c("tf_family", "family", "TF_family"))
  rename_first("motif_ID_count", c("motif_ID_count", "motif_ids", "unique_motif_ids",
                                  "motif_count", "high_confidence_motif_ids"))
  rename_first("delta_vs_Col0", c("delta_vs_Col0", "delta_vs_col0", "delta"))
  x
}

if (!is.na(family_long_file) && !ONSEN_FORCE_RESCAN) {
  family_long <- normalize_family_long(read_table_auto(family_long_file))
} else {
  if (!ONSEN_RUN_LARGE_STEPS) {
    stop(
      "Natural variant TF-family summary was not found and large motif scanning is disabled."
    )
  }
  jaspar_file <- find_input("JASPAR2026_CORE_plants_nonredundant_pfms_jaspar.txt")
  motifs <- parse_jaspar_pfms(jaspar_file)

  sequence_table <- dplyr::bind_rows(
    data.frame(sequence_id = "Col-0 reference", sequence = seed),
    variant_summary |>
      dplyr::transmute(sequence_id = variant_id, sequence = seed_sequence)
  )
  hits <- as.data.frame(scan_sequences_against_motifs(
    sequence_table, motifs, threshold = 0.85, pseudocount = 0.8
  ))
  family_long <- hits |>
    dplyr::distinct(sequence_id, tf_family, motif_id) |>
    dplyr::count(sequence_id, tf_family, name = "motif_ID_count") |>
    tidyr::complete(
      sequence_id = sequence_table$sequence_id,
      tf_family,
      fill = list(motif_ID_count = 0L)
    ) |>
    dplyr::rename(variant_id = sequence_id)
}
safe_write_csv(family_long, "natural_variant_TF_family_absolute_repository.csv")

if (!is.na(family_delta_file) && !ONSEN_FORCE_RESCAN) {
  family_delta <- normalize_family_long(read_table_auto(family_delta_file))
} else {
  reference <- family_long |>
    dplyr::filter(grepl("Col-0", variant_id, ignore.case = TRUE)) |>
    dplyr::select(tf_family, Col0_motif_ID_count = motif_ID_count)
  family_delta <- family_long |>
    dplyr::left_join(reference, by = "tf_family") |>
    dplyr::mutate(delta_vs_Col0 = motif_ID_count - Col0_motif_ID_count)
}
safe_write_csv(family_delta, "natural_variant_TF_family_delta_repository.csv")

# -------------------------------- Figures ------------------------------------
if (ONSEN_MAKE_FIGURES) {
  accession_order <- c("Col-0", "An-1", "C24", "Cvi", "Eri", "Kyo", "Ler", "Sha")

  # Fig. S5A
  s1a <- candidate_summary |>
    dplyr::select(
      accession, exact_seed_matches,
      one_to_two_mismatch_seeds, three_to_four_mismatch_seeds
    ) |>
    tidyr::pivot_longer(
      -accession, names_to = "seed_class", values_to = "count"
    ) |>
    dplyr::mutate(
      accession = factor(accession, levels = accession_order),
      seed_class = factor(
        seed_class,
        levels = c(
          "exact_seed_matches",
          "one_to_two_mismatch_seeds",
          "three_to_four_mismatch_seeds"
        ),
        labels = c("Exact seed", "1-2 mismatches", "3-4 mismatches")
      )
    )
  p_s1a <- ggplot2::ggplot(
    s1a, ggplot2::aes(accession, count, fill = seed_class)
  ) +
    ggplot2::geom_col(colour = "black", linewidth = 0.2) +
    ggplot2::labs(
      x = "Accession", y = "Main-chromosome candidate windows",
      fill = "Seed class"
    ) +
    theme_onsen(12) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  save_plot_pair(p_s1a, "FigS5A_accession_candidate_abundance", 7.0, 5.2)

  # Fig. S5B
  hsf_candidates$accession <- factor(hsf_candidates$accession, levels = accession_order)
  p_s1b <- ggplot2::ggplot(
    hsf_candidates,
    ggplot2::aes(accession, HSF_hits_per_kb, fill = accession)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, colour = "black") +
    ggplot2::geom_jitter(width = 0.15, size = 1.3, alpha = 0.55) +
    ggplot2::labs(
      x = "Accession", y = "HSF-family motif-position hits per kb"
    ) +
    theme_onsen(12) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  save_plot_pair(p_s1b, "FigS5B_accession_HSF_density", 7.0, 5.2)

  # Fig. S5C
  proxy_plot <- proxy_summary
  names(proxy_plot) <- gsub(" ", "_", names(proxy_plot))
  paired_col <- c(
    "paired_HSE_LTR_copy_proxy", "putative_paired_copy_proxy",
    "paired_candidates", "paired_proxy"
  )[c(
    "paired_HSE_LTR_copy_proxy", "putative_paired_copy_proxy",
    "paired_candidates", "paired_proxy"
  ) %in% names(proxy_plot)][1]
  unpaired_col <- c(
    "unpaired_HSE_candidate", "unpaired_candidates", "unpaired"
  )[c(
    "unpaired_HSE_candidate", "unpaired_candidates", "unpaired"
  ) %in% names(proxy_plot)][1]
  if (!is.na(paired_col)) {
    if (is.na(unpaired_col)) {
      candidate_col <- c("candidate_windows", "candidate_count")[
        c("candidate_windows", "candidate_count") %in% names(proxy_plot)
      ][1]
      proxy_plot$unpaired_HSE_candidate <- proxy_plot[[candidate_col]] -
        2 * proxy_plot[[paired_col]]
      unpaired_col <- "unpaired_HSE_candidate"
    }
    s1c <- proxy_plot |>
      dplyr::transmute(
        accession,
        `Paired HSE/LTR copy proxy` = .data[[paired_col]],
        `Unpaired HSE candidate` = .data[[unpaired_col]]
      ) |>
      tidyr::pivot_longer(
        -accession, names_to = "proxy_class", values_to = "count"
      ) |>
      dplyr::mutate(accession = factor(accession, levels = accession_order))
    p_s1c <- ggplot2::ggplot(
      s1c, ggplot2::aes(count, accession, fill = proxy_class)
    ) +
      ggplot2::geom_col(position = "stack", colour = "black", linewidth = 0.2) +
      ggplot2::labs(x = "Count", y = "Accession", fill = NULL) +
      theme_onsen(12) +
      ggplot2::theme(legend.position = "top")
    save_plot_pair(p_s1c, "FigS5C_paired_unpaired_proxies", 7.0, 5.2)
  }

  # Fig. S4A
  architecture$accession <- factor(architecture$accession, levels = accession_order)
  p_s2a <- ggplot2::ggplot(
    architecture,
    ggplot2::aes(accession, canonical_HSE_units, fill = accession)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, colour = "black") +
    ggplot2::geom_jitter(width = 0.15, size = 1.3, alpha = 0.55) +
    ggplot2::geom_hline(yintercept = 5, linetype = "dashed") +
    ggplot2::labs(x = "Accession", y = "Canonical HSE-like units in 49-bp seed") +
    theme_onsen(12) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  save_plot_pair(p_s2a, "FigS4A_accession_HSE_architecture", 7.0, 5.2)

  # Fig. S4B scaled accession summary
  s2b_metrics <- candidate_summary |>
    dplyr::mutate(
      accession = factor(accession, levels = accession_order)
    ) |>
    tidyr::pivot_longer(
      -accession, names_to = "metric", values_to = "value"
    ) |>
    dplyr::group_by(metric) |>
    dplyr::mutate(
      scaled_value = ifelse(
        max(value, na.rm = TRUE) == min(value, na.rm = TRUE),
        0,
        (value - min(value, na.rm = TRUE)) /
          (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))
      )
    ) |>
    dplyr::ungroup()
  p_s2b <- ggplot2::ggplot(
    s2b_metrics,
    ggplot2::aes(metric, accession, fill = scaled_value)
  ) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_gradient(
      low = "#F7F7F7", high = "#5A4AA8",
      limits = c(0, 1), name = "Scaled\nwithin metric"
    ) +
    ggplot2::labs(x = NULL, y = "Accession") +
    theme_onsen(10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_plot_pair(p_s2b, "FigS4B_scaled_accession_summary", 9.0, 5.2)

  # Fig. S4C HSF-model heatmap from exact processed model summary.
  if (!is.na(hsf_model_summary_file)) {
    model_summary <- read_table_auto(hsf_model_summary_file)
    accession_col <- c("accession", "Accession")[
      c("accession", "Accession") %in% names(model_summary)
    ][1]
    model_col <- c("HSF_model", "motif_name", "tf_name", "model")[
      c("HSF_model", "motif_name", "tf_name", "model") %in% names(model_summary)
    ][1]
    density_col <- names(model_summary)[
      grepl("median.*density|motif.*density|hits_per_kb", names(model_summary), ignore.case = TRUE)
    ][1]
    if (!is.na(accession_col) && !is.na(model_col) && !is.na(density_col)) {
      s2c <- model_summary |>
        dplyr::transmute(
          accession = factor(.data[[accession_col]], levels = accession_order),
          HSF_model = .data[[model_col]],
          density = as.numeric(.data[[density_col]])
        ) |>
        dplyr::group_by(HSF_model) |>
        dplyr::mutate(
          scaled_density = ifelse(
            max(density, na.rm = TRUE) == min(density, na.rm = TRUE),
            0,
            (density - min(density, na.rm = TRUE)) /
              (max(density, na.rm = TRUE) - min(density, na.rm = TRUE))
          )
        ) |>
        dplyr::ungroup()
      p_s2c <- ggplot2::ggplot(
        s2c,
        ggplot2::aes(HSF_model, accession, fill = scaled_density)
      ) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::scale_fill_gradient(
          low = "#F7F7F7", high = "#5A4AA8",
          limits = c(0, 1), name = "Scaled median\nmotif density"
        ) +
        ggplot2::labs(x = "Arabidopsis HSF motif model", y = "Accession") +
        theme_onsen(10) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      save_plot_pair(p_s2c, "FigS4C_HSF_model_compatibility", 8.5, 5.2)
    }
  }

  # Fig. 6 natural variant nucleotide matrix.
  variant_id_col <- c("variant_id", "sequence_id", "variant")[
    c("variant_id", "sequence_id", "variant") %in% names(variant_summary)
  ][1]
  sequence_col <- c("seed_sequence", "sequence", "variant_sequence")[
    c("seed_sequence", "sequence", "variant_sequence") %in% names(variant_summary)
  ][1]
  substitution_col <- c("substitution_count", "n_variants", "mismatch_count")[
    c("substitution_count", "n_variants", "mismatch_count") %in% names(variant_summary)
  ][1]
  if (!is.na(variant_id_col) && !is.na(sequence_col)) {
    variant_plot <- variant_summary |>
      dplyr::transmute(
        variant_id = .data[[variant_id_col]],
        accession,
        sequence = .data[[sequence_col]],
        substitution_count = if (!is.na(substitution_col)) {
          as.integer(.data[[substitution_col]])
        } else {
          vapply(.data[[sequence_col]], hamming_distance, numeric(1), b = seed)
        }
      )
    plot_sequences <- dplyr::bind_rows(
      data.frame(
        variant_id = "Col-0 reference",
        accession = "Col-0",
        sequence = seed,
        substitution_count = 0L
      ),
      variant_plot
    )
    hse_positions <- unlist(Map(
      seq,
      find_canonical_hse_units(seed)$start,
      find_canonical_hse_units(seed)$end
    ))
    seed_bases <- split_bases(seed)
    fig6 <- plot_sequences |>
      tidyr::expand_grid(position = seq_along(seed_bases)) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        base = substring(sequence, position, position),
        reference_base = seed_bases[[position]],
        display_base = ifelse(
          variant_id == "Col-0 reference" || base != reference_base,
          base, "."
        ),
        tile_class = dplyr::case_when(
          variant_id == "Col-0 reference" & position %in% hse_positions ~ "Col-0 HSE unit",
          variant_id == "Col-0 reference" ~ "Col-0 reference",
          base == reference_base & position %in% hse_positions ~ "Same inside HSE",
          base == reference_base ~ "Same as Col-0",
          base != reference_base & position %in% hse_positions ~ "Variant inside HSE",
          TRUE ~ "Variant outside HSE"
        ),
        row_label = paste0(variant_id, " | ", substitution_count, " variants")
      ) |>
      dplyr::ungroup()
    p6 <- ggplot2::ggplot(
      fig6,
      ggplot2::aes(position, forcats::fct_rev(row_label), fill = tile_class)
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.1) +
      ggplot2::geom_text(
        ggplot2::aes(label = display_base),
        size = 2.0, fontface = "bold"
      ) +
      ggplot2::scale_fill_manual(values = c(
        "Col-0 reference" = "#D8C7FF",
        "Col-0 HSE unit" = "#F6D36F",
        "Same as Col-0" = "#F6F6F6",
        "Same inside HSE" = "#FFF2B8",
        "Variant inside HSE" = "#E76F9E",
        "Variant outside HSE" = "#6AA5D9"
      )) +
      ggplot2::scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 49)) +
      ggplot2::labs(
        x = "Position in 49-bp ONSEN HSE seed",
        y = NULL, fill = NULL
      ) +
      theme_onsen(9) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 6.5),
        legend.position = "bottom"
      )
    save_plot_pair(p6, "Fig6_natural_accession_HSE_seed_variants", 11, 10)
  }

  # Fig. 7A absolute HSF/AP2/ERF/DOF family counts.
  main_families <- c("HSF", "AP2/ERF", "DOF")
  fig7a <- family_long |>
    dplyr::filter(tf_family %in% main_families) |>
    dplyr::mutate(tf_family = factor(tf_family, levels = main_families))
  p7a <- ggplot2::ggplot(
    fig7a,
    ggplot2::aes(tf_family, motif_ID_count, group = variant_id, colour = tf_family)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.45, size = 1.3) +
    ggplot2::labs(
      x = "TF family",
      y = "High-confidence motif IDs per variant"
    ) +
    theme_onsen(12) +
    ggplot2::theme(legend.position = "none")
  save_plot_pair(p7a, "Fig7A_natural_variant_family_counts", 6.5, 5.0)

  # Fig. 7B delta heatmap. Show only changed variants.
  delta_col <- c("delta_vs_Col0", "delta_vs_col0", "delta")[
    c("delta_vs_Col0", "delta_vs_col0", "delta") %in% names(family_delta)
  ][1]
  family_delta$delta_vs_Col0 <- family_delta[[delta_col]]
  changed_ids <- family_delta |>
    dplyr::group_by(variant_id) |>
    dplyr::summarise(any_change = any(delta_vs_Col0 != 0, na.rm = TRUE), .groups = "drop") |>
    dplyr::filter(any_change) |>
    dplyr::pull(variant_id)
  plot_families <- c("HSF", "AP2/ERF", "DOF", "MYB", "NAC", "bZIP", "WRKY")
  fig7b <- family_delta |>
    dplyr::filter(
      variant_id %in% changed_ids,
      tf_family %in% plot_families
    ) |>
    dplyr::mutate(
      tf_family = factor(tf_family, levels = plot_families),
      variant_id = factor(variant_id, levels = rev(changed_ids))
    )
  p7b <- ggplot2::ggplot(
    fig7b,
    ggplot2::aes(tf_family, variant_id, fill = delta_vs_Col0)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(delta_vs_Col0 == 0, "", sprintf("%+d", delta_vs_Col0))),
      size = 2.0
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#4C78A8", mid = "white", high = "#E45786",
      midpoint = 0, name = "Delta motif IDs\nvs Col-0"
    ) +
    ggplot2::labs(x = "TF family", y = NULL) +
    theme_onsen(9) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6.5))
  save_plot_pair(p7b, "Fig7B_natural_variant_TF_family_delta", 8.5, 10)
}

message("Natural-accession ONSEN-like HSE analysis completed.")
