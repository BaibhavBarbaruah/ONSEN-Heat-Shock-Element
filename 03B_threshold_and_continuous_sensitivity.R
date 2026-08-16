# HSF threshold-sensitivity and continuous-score analysis for Fig. S2/Table S7.
#
# Principal threshold metric: unique 1-based forward-sequence (start, end)
# coordinates after collapsing duplicate HSF PWM-model/strand predictions.
# Different coordinates remain separate. Continuous-score summaries are
# threshold-independent and use all HSF PWM scores.

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "tidyr", "ggplot2", "patchwork", "Biostrings"))
message_config()

THRESHOLDS <- c(0.80, 0.85, 0.90, 0.95)
PSEUDOCOUNT <- 0.8
HSF_NAMES <- c("HSFC1", "HSFA6B", "HSFB3", "HSFA6A", "HSFA4A",
               "HSFA1E", "HSFB4", "HSFB2B", "HSFB2A", "HSFA1B")
BG_CLASS <- "Strict non-ONSEN TE background"
ONSEN_CLASS <- "ONSEN terminal candidate windows"

read_final_s7 <- function(sheet_stem) {
  path <- repo_file(file.path("supplementary_table_source", sheet_stem))
  x <- data.table::fread(path, skip = 2L, data.table = FALSE, check.names = FALSE)
  x
}

# The repository ships the final region-level source tables, so the lightweight
# default run validates/copies those outputs. Set ONSEN_RUN_LARGE_STEPS=true to
# rescan TAIR10 and reconstruct the metrics from sequence and JASPAR PFMs.
if (!ONSEN_RUN_LARGE_STEPS && !ONSEN_FORCE_RESCAN) {
  region_threshold <- read_final_s7("Table_S7__S7D_region_threshold.tsv")
  region_continuous <- read_final_s7("Table_S7__S7E_region_continuous.tsv")
  threshold_stats <- read_final_s7("Table_S7__S7A_threshold_stats.tsv")
  continuous_stats <- read_final_s7("Table_S7__S7C_continuous_stats.tsv")

  stopifnot(
    nrow(region_threshold) == (16L + 1930L) * length(THRESHOLDS),
    sum(region_threshold[["Region class"]] == ONSEN_CLASS &
          region_threshold[["Relative PWM-score threshold"]] == 0.85) == 16L,
    sum(region_threshold[["Region class"]] == BG_CLASS &
          region_threshold[["Relative PWM-score threshold"]] == 0.85) == 1930L
  )

  key85 <- threshold_stats[threshold_stats[["Relative PWM-score threshold"]] == 0.85, ]
  stopifnot(
    nrow(key85) == 1L,
    abs(key85[["ONSEN median density (placements/kb)"]] - 49.375) < 1e-10,
    abs(key85[["Background median density (placements/kb)"]] - 6.81238953796413) < 1e-10
  )

  data.table::fwrite(region_threshold, out_file("S7D_region_threshold_metrics_repository.tsv"), sep = "\t")
  data.table::fwrite(region_continuous, out_file("S7E_region_continuous_metrics_repository.tsv"), sep = "\t")
  data.table::fwrite(threshold_stats, out_file("S7A_threshold_statistics_repository.tsv"), sep = "\t")
  data.table::fwrite(continuous_stats, out_file("S7C_continuous_statistics_repository.tsv"), sep = "\t")
  message("Validated deposited final Table S7 source data (16 ONSEN + 1,930 background regions).")
} else {
  BG_FILE <- find_input("strict_TE_only_background_coordinates.csv")
  GENOME_FILE <- find_input("Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz")
  JASPAR_FILE <- find_input("JASPAR2026_CORE_plants_nonredundant_pfms_jaspar.txt")

  onsen_source <- read.csv(repo_file("ONSEN_Col0_terminal_candidate_windows.csv"), check.names = FALSE)
  onsen <- data.frame(
    region_id = onsen_source$window_id,
    chromosome = standard_chr(onsen_source$chromosome),
    start = as.integer(onsen_source$start_1based),
    end = as.integer(onsen_source$end_1based),
    width_bp = as.integer(onsen_source$width_bp),
    class = ONSEN_CLASS,
    stringsAsFactors = FALSE
  )
  stopifnot(nrow(onsen) == 16L, all(onsen$width_bp == 800L))

  bg_raw <- data.table::fread(BG_FILE, data.table = FALSE, check.names = FALSE)
  pick <- function(candidates, required = TRUE) {
    hit <- candidates[candidates %in% names(bg_raw)][1]
    if (is.na(hit) && required) stop("Missing required strict-TE coordinate column: ", paste(candidates, collapse = ", "))
    hit
  }
  chr_col <- pick(c("chr_clean", "seqid", "chromosome", "chr", "Chr", "seqnames", "chrom"))
  start_col <- pick(c("start", "Start", "start_1based", "region_start", "feature_start"))
  end_col <- pick(c("end", "End", "end_1based", "region_end", "feature_end"))
  id_col <- pick(c("background_id", "fasta_id", "region_id", "sequence_id", "candidate_id", "feature_id", "GeneID", "ID"), FALSE)
  family_col <- pick(c("TE_family", "te_family", "family", "Alias", "alias", "transposon_family"), FALSE)

  bg <- data.frame(
    region_id = if (is.na(id_col)) paste0("strict_TE_", seq_len(nrow(bg_raw))) else as.character(bg_raw[[id_col]]),
    chromosome = standard_chr(bg_raw[[chr_col]]),
    start = as.integer(bg_raw[[start_col]]),
    end = as.integer(bg_raw[[end_col]]),
    TE_family = if (is.na(family_col)) "" else as.character(bg_raw[[family_col]]),
    stringsAsFactors = FALSE
  )
  bg$region_id[is.na(bg$region_id) | bg$region_id == ""] <- paste0("strict_TE_", which(is.na(bg$region_id) | bg$region_id == ""))
  bg$region_id <- make.unique(bg$region_id)
  bg$width_bp <- bg$end - bg$start + 1L

  direct_overlap <- vapply(seq_len(nrow(bg)), function(i) {
    same <- onsen$chromosome == bg$chromosome[i]
    any(same & onsen$start <= bg$end[i] & onsen$end >= bg$start[i])
  }, logical(1))
  bg <- bg[!direct_overlap, , drop = FALSE]
  residual_onsen <- grepl("ONSEN|ATCOPIA[ _-]?78|COPIA[ _-]?78",
                          paste(bg$region_id, bg$TE_family), ignore.case = TRUE)
  bg <- bg[!residual_onsen, , drop = FALSE]
  if (nrow(bg) != 1930L) stop("Expected 1,930 strict non-ONSEN TE regions after harmonisation; observed ", nrow(bg), ".")
  bg$class <- BG_CLASS

  genome <- Biostrings::readDNAStringSet(GENOME_FILE)
  names(genome) <- standard_chr(names(genome))
  coordinates <- dplyr::bind_rows(
    onsen |> dplyr::select(region_id, chromosome, start, end, width_bp, class),
    bg |> dplyr::select(region_id, chromosome, start, end, width_bp, class)
  )
  if (any(!coordinates$chromosome %in% names(genome))) stop("Coordinate chromosome missing from TAIR10 FASTA.")
  coordinates$sequence <- mapply(
    function(chr, start, end) as.character(Biostrings::subseq(genome[[chr]], start = start, end = end)),
    coordinates$chromosome, coordinates$start, coordinates$end, USE.NAMES = FALSE
  )

  motifs_all <- parse_jaspar_pfms(JASPAR_FILE)
  motif_names <- toupper(sub("\\s+.*$", "", trimws(vapply(motifs_all, `[[`, character(1), "name"))))
  hsf_motifs <- lapply(HSF_NAMES, function(name) {
    idx <- which(motif_names == name)
    if (!length(idx)) stop("Missing HSF PWM model in JASPAR collection: ", name)
    motifs_all[[idx[1]]]
  })
  names(hsf_motifs) <- HSF_NAMES

  scan_region <- function(region_id, sequence, region_class, chromosome, start, end, width_bp) {
    coordinate_sets <- lapply(THRESHOLDS, function(x) data.frame(forward_start = integer(), forward_end = integer()))
    model_maxima <- numeric(length(hsf_motifs))
    all_scores <- numeric()

    for (m in seq_along(hsf_motifs)) {
      hits <- scan_one_motif(region_id, sequence, hsf_motifs[[m]],
                             threshold = -Inf, pseudocount = PSEUDOCOUNT, retain_all = TRUE)
      scores <- hits$relative_score[is.finite(hits$relative_score)]
      model_maxima[m] <- max(scores)
      all_scores <- c(all_scores, scores)
      for (t in seq_along(THRESHOLDS)) {
        z <- hits[is.finite(hits$relative_score) & hits$relative_score >= THRESHOLDS[t],
                  c("forward_start", "forward_end"), drop = FALSE]
        coordinate_sets[[t]] <- rbind(coordinate_sets[[t]], z)
      }
    }

    threshold_rows <- lapply(seq_along(THRESHOLDS), function(t) {
      unique_coords <- unique(coordinate_sets[[t]])
      n <- nrow(unique_coords)
      data.frame(
        region_id = region_id, class = region_class, chromosome = chromosome,
        start = start, end = end, width_bp = width_bp,
        threshold = THRESHOLDS[t],
        exact_coordinate_HSF_placements = n,
        HSF_hits_per_kb = n / width_bp * 1000,
        stringsAsFactors = FALSE
      )
    })
    top5 <- head(sort(all_scores, decreasing = TRUE), 5L)
    continuous <- data.frame(
      region_id = region_id, class = region_class, chromosome = chromosome,
      start = start, end = end, width_bp = width_bp,
      maximum_HSF_relative_score = max(all_scores),
      mean_top5_HSF_relative_score = mean(top5),
      median_per_model_maximum_score = median(model_maxima),
      stringsAsFactors = FALSE
    )
    list(threshold = dplyr::bind_rows(threshold_rows), continuous = continuous)
  }

  scanned <- lapply(seq_len(nrow(coordinates)), function(i) {
    scan_region(coordinates$region_id[i], coordinates$sequence[i], coordinates$class[i],
                coordinates$chromosome[i], coordinates$start[i], coordinates$end[i], coordinates$width_bp[i])
  })
  region_metrics <- dplyr::bind_rows(lapply(scanned, `[[`, "threshold"))
  continuous_metrics <- dplyr::bind_rows(lapply(scanned, `[[`, "continuous"))

  threshold_stats <- dplyr::bind_rows(lapply(THRESHOLDS, function(cutoff) {
    z <- region_metrics[region_metrics$threshold == cutoff, ]
    a <- z$HSF_hits_per_kb[z$class == ONSEN_CLASS]
    b <- z$HSF_hits_per_kb[z$class == BG_CLASS]
    test <- stats::wilcox.test(a, b, exact = FALSE)
    data.frame(
      threshold = cutoff, n_ONSEN = length(a), n_background = length(b),
      ONSEN_median = median(a), background_median = median(b),
      Wilcoxon_W = unname(test$statistic), Wilcoxon_P = test$p.value,
      Cliffs_delta = cliffs_delta(a, b), stringsAsFactors = FALSE
    )
  }))
  threshold_stats$BH_adjusted_P <- p.adjust(threshold_stats$Wilcoxon_P, method = "BH")

  metric_columns <- c("maximum_HSF_relative_score", "mean_top5_HSF_relative_score",
                      "median_per_model_maximum_score")
  continuous_stats <- dplyr::bind_rows(lapply(metric_columns, function(metric) {
    a <- continuous_metrics[[metric]][continuous_metrics$class == ONSEN_CLASS]
    b <- continuous_metrics[[metric]][continuous_metrics$class == BG_CLASS]
    test <- stats::wilcox.test(a, b, exact = FALSE)
    data.frame(metric = metric, n_ONSEN = length(a), n_background = length(b),
               ONSEN_median = median(a), background_median = median(b),
               Wilcoxon_P = test$p.value, Cliffs_delta = cliffs_delta(a, b),
               stringsAsFactors = FALSE)
  }))
  continuous_stats$BH_adjusted_P <- p.adjust(continuous_stats$Wilcoxon_P, method = "BH")

  data.table::fwrite(region_metrics, out_file("Revision_R1_3_region_level_HSF_metrics_all_thresholds.csv"))
  data.table::fwrite(continuous_metrics, out_file("Revision_R1_3_region_level_continuous_HSF_scores.csv"))
  data.table::fwrite(threshold_stats, out_file("Revision_R1_3_threshold_sensitivity_statistics.csv"))
  data.table::fwrite(continuous_stats, out_file("Revision_R1_3_continuous_score_statistics.csv"))

  key <- threshold_stats[threshold_stats$threshold == 0.85, ]
  if (abs(key$ONSEN_median - 49.375) > 1e-9 || abs(key$background_median - 6.81238953796413) > 1e-9) {
    stop("Full rescan did not reproduce the final 0.85 Table S7 medians; check input versions and motif collection.")
  }

  message("Full exact-coordinate HSF threshold/continuous rescan completed successfully.")
}
