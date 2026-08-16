# Natural-accession ONSEN-like candidate analysis.
# Covers Fig. 6-7 accession summaries, Fig. S3-S4 and Tables S11-S13.
#
# The deposited final source tables are authoritative for accession-level
# abundance, exact-coordinate HSF density and seed-variant PWM-model summaries.
# This script validates those outputs and, when the processed chromosome-level
# candidate files are supplied, checks their consistency without substituting a
# raw PWM-model-position hit count for the manuscript's exact-coordinate metric.

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "tidyr", "ggplot2"))
message_config()

ACCESSIONS <- c("Col-0", "An-1", "C24", "Cvi", "Eri", "Kyo", "Ler", "Sha")
EXPECTED_CANDIDATES <- c(`Col-0` = 19L, `An-1` = 7L, C24 = 11L, Cvi = 9L,
                         Eri = 16L, Kyo = 17L, Ler = 15L, Sha = 41L)
EXPECTED_PROXIES <- c(`Col-0` = 8L, `An-1` = 2L, C24 = 5L, Cvi = 4L,
                      Eri = 4L, Kyo = 7L, Ler = 6L, Sha = 16L)

read_final <- function(filename, skip = 1L) {
  data.table::fread(
    repo_file(file.path("supplementary_table_source", filename)),
    skip = skip, data.table = FALSE, check.names = FALSE
  )
}

s11 <- read_final("Table_S11__S11_accessions.tsv")
s12 <- read_final("Table_S12__S12_seed_variants.tsv")
s13 <- read_final("Table_S13__S13_TF_family.tsv")

assert_columns(
  s11,
  c(
    "Accession", "Candidate windows",
    "Median non-redundant HSF motif-coordinate placements per kb",
    "Median distinct HSF PWM models", "Putative paired HSE/LTR proxies"
  ),
  "Table S11 source"
)
assert_columns(
  s12,
  c(
    "Accession", "Candidate windows", "Exact seeds", "1-2 mismatch seeds",
    "3-4 mismatch seeds", "Median mismatch", "Maximum mismatch",
    "Median HSE-like units retained"
  ),
  "Table S12 source"
)
assert_columns(
  s13,
  c(
    "Accession", "Variant seed", "Mismatch count", "Mismatch positions",
    "Distinct HSF PWM models", "Distinct DOF PWM models",
    "Distinct AP2/ERF PWM models", "Δ HSF PWM models vs Col-0",
    "Δ DOF PWM models vs Col-0"
  ),
  "Table S13 source"
)

if (!identical(as.character(s11$Accession), ACCESSIONS)) {
  stop("Table S11 accession order/set does not match the final manuscript.")
}
observed_candidates <- setNames(as.integer(s11[["Candidate windows"]]), s11$Accession)
observed_proxies <- setNames(as.integer(s11[["Putative paired HSE/LTR proxies"]]), s11$Accession)
if (!identical(observed_candidates[ACCESSIONS], EXPECTED_CANDIDATES[ACCESSIONS])) {
  stop("Table S11 candidate-window abundance does not match the final values.")
}
if (!identical(observed_proxies[ACCESSIONS], EXPECTED_PROXIES[ACCESSIONS])) {
  stop("Table S11 paired HSE/LTR proxy counts do not match the final values.")
}
if (!identical(as.integer(s12[["Candidate windows"]]), as.integer(s11[["Candidate windows"]]))) {
  stop("Tables S11 and S12 disagree on accession candidate abundance.")
}
if (any(s12[["Maximum mismatch"]] > 4L, na.rm = TRUE)) {
  stop("Detected candidate set exceeds the <=4-mismatch ascertainment rule.")
}
if (!any(s13$Accession == "Col-0" & s13[["Mismatch count"]] == 0L &
         s13[["Distinct HSF PWM models"]] == 15L &
         s13[["Distinct DOF PWM models"]] == 8L)) {
  stop("Table S13 lacks the Col-0 reference PWM-model row.")
}

# Manuscript-level recurrent pattern checks. These are distinct PWM-model
# identities, not independent binding sites.
if (min(s13[["Δ HSF PWM models vs Col-0"]], na.rm = TRUE) > -5L) {
  stop("Expected Cvi-derived reduction of up to five distinct HSF PWM models was not found.")
}
if (max(s13[["Δ DOF PWM models vs Col-0"]], na.rm = TRUE) != 7L) {
  stop("Expected maximum natural-variant DOF gain of seven distinct PWM models was not found.")
}

safe_write_csv(s11, "Table_S11_accession_summary_repository.csv")
safe_write_csv(s12, "Table_S12_seed_variant_summary_repository.csv")
safe_write_csv(s13, "Table_S13_TF_family_summary_repository.csv")

# Optional consistency checks against the processed chromosome-level inputs
# listed in INPUT_PROVENANCE.tsv. These inputs encode the exact final candidate
# discovery and exact-coordinate HSF scan; no alternative raw-hit metric is
# substituted when they are absent.
if (ONSEN_RUN_LARGE_STEPS || ONSEN_FORCE_RESCAN) {
  candidate_file <- find_input(
    "accession_ONSEN_like_HSE_candidate_windows_mismatch_leq4_metadata.csv"
  )
  hsf_file <- find_input(
    "FIXED_accession_ONSEN_like_mainchr_candidate_windows_JASPAR2026_Arabidopsis_HSF_candidate_summary_threshold_0.85.csv"
  )
  proxy_file <- find_input("putative_ONSEN_like_copy_proxy_summary_mismatch_leq4.csv")

  candidates <- read_table_auto(candidate_file)
  hsf <- read_table_auto(hsf_file)
  proxies <- read_table_auto(proxy_file)

  acc_col <- c("accession", "Accession")[c("accession", "Accession") %in% names(candidates)][1]
  main_col <- c("main_chromosome", "main_chr", "is_main_chromosome")[
    c("main_chromosome", "main_chr", "is_main_chromosome") %in% names(candidates)
  ][1]
  mismatch_col <- c("mismatches_to_Col0", "mismatch_count", "mismatches")[
    c("mismatches_to_Col0", "mismatch_count", "mismatches") %in% names(candidates)
  ][1]
  if (is.na(acc_col) || is.na(mismatch_col)) stop("Candidate metadata lacks accession/mismatch columns.")
  main_candidates <- if (!is.na(main_col)) candidates[candidates[[main_col]] %in% TRUE, , drop = FALSE] else candidates
  counts <- table(factor(main_candidates[[acc_col]], levels = ACCESSIONS))
  if (!identical(as.integer(counts), as.integer(EXPECTED_CANDIDATES[ACCESSIONS]))) {
    stop("Processed accession candidate metadata does not reproduce Table S11 candidate counts.")
  }
  if (any(as.integer(main_candidates[[mismatch_col]]) > 4L, na.rm = TRUE)) {
    stop("Processed candidate metadata violates the <=4-mismatch rule.")
  }

  hsf_acc_col <- c("accession", "Accession")[c("accession", "Accession") %in% names(hsf)][1]
  density_col <- names(hsf)[grepl("HSF.*(per_kb|density)|placements.*per.*kb", names(hsf), ignore.case = TRUE)][1]
  if (is.na(hsf_acc_col) || is.na(density_col)) {
    stop("Processed accession HSF summary lacks accession or exact-coordinate density columns.")
  }
  hsf_medians <- hsf |>
    dplyr::group_by(.data[[hsf_acc_col]]) |>
    dplyr::summarise(median_density = median(as.numeric(.data[[density_col]]), na.rm = TRUE), .groups = "drop")
  names(hsf_medians)[1] <- "Accession"
  check <- s11 |>
    dplyr::select(Accession, expected = .data[["Median non-redundant HSF motif-coordinate placements per kb"]]) |>
    dplyr::left_join(hsf_medians, by = "Accession")
  if (any(abs(check$expected - check$median_density) > 1e-8, na.rm = TRUE)) {
    stop("Processed exact-coordinate HSF summary does not reproduce Table S11 medians.")
  }

  proxy_acc_col <- c("accession", "Accession")[c("accession", "Accession") %in% names(proxies)][1]
  proxy_count_col <- names(proxies)[grepl("paired.*proxy|paired_HSE_LTR", names(proxies), ignore.case = TRUE)][1]
  if (!is.na(proxy_acc_col) && !is.na(proxy_count_col)) {
    proxy_map <- setNames(as.integer(proxies[[proxy_count_col]]), proxies[[proxy_acc_col]])
    if (!identical(proxy_map[ACCESSIONS], EXPECTED_PROXIES[ACCESSIONS])) {
      stop("Processed structural-proxy summary does not reproduce Table S11.")
    }
  }
}

if (ONSEN_MAKE_FIGURES) {
  s11$Accession <- factor(s11$Accession, levels = ACCESSIONS)

  p_s4a <- ggplot2::ggplot(s11, ggplot2::aes(Accession, .data[["Candidate windows"]])) +
    ggplot2::geom_col(colour = "black") +
    ggplot2::labs(x = NULL, y = "Candidate windows") + theme_onsen(12)
  save_plot_pair(p_s4a, "FigS4A_accession_candidate_abundance", 6.4, 4.8)

  p_s4b <- ggplot2::ggplot(
    s11,
    ggplot2::aes(Accession, .data[["Median non-redundant HSF motif-coordinate placements per kb"]])
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = NULL, y = "Median HSF motif-coordinate\nplacements per kb") + theme_onsen(12)
  save_plot_pair(p_s4b, "FigS4B_accession_HSF_density_summary", 6.4, 4.8)

  proxy_long <- s11 |>
    dplyr::transmute(
      Accession,
      `Paired HSE/LTR copy proxy` = .data[["Putative paired HSE/LTR proxies"]],
      `Unpaired HSE candidate` = .data[["Candidate windows"]] - 2 * .data[["Putative paired HSE/LTR proxies"]]
    ) |>
    tidyr::pivot_longer(-Accession, names_to = "candidate_type", values_to = "count")
  p_s4c <- ggplot2::ggplot(proxy_long, ggplot2::aes(count, Accession, shape = candidate_type)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = "Count", y = NULL, shape = NULL) + theme_onsen(12) +
    ggplot2::theme(legend.position = "top")
  save_plot_pair(p_s4c, "FigS4C_accession_structural_proxies", 6.4, 4.8)
}

message("Natural-accession final-source validation completed.")
