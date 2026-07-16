# HSF threshold-sensitivity and continuous-score analysis.
# Reconstructs the analysis underlying final Fig. S3 and Table S8.

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "tidyr", "ggplot2", "patchwork", "openxlsx", "Biostrings"))
message_config()

THRESHOLDS <- c(0.80, 0.85, 0.90, 0.95)
PSEUDOCOUNT <- 0.8
HSF_NAMES <- c("HSFC1", "HSFA6B", "HSFB3", "HSFA6A", "HSFA4A", "HSFA1E", "HSFB4", "HSFB2B", "HSFB2A", "HSFA1B")
BG_CLASS <- "Strict non-ONSEN TE background"
ONSEN_CLASS <- "ONSEN terminal candidate windows"
CLASS_COLOURS <- setNames(c("#83AEE8", "#E989AE"), c(BG_CLASS, ONSEN_CLASS))

BG_FILE <- find_input("strict_TE_only_background_coordinates.csv")
GENOME_FILE <- find_input("Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz")
JASPAR_FILE <- find_input("JASPAR2026_CORE_plants_nonredundant_pfms_jaspar.txt")
OUT_DIR <- ONSEN_OUTPUT_ROOT
FIG_DIR <- ONSEN_OUTPUT_ROOT
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

standard_chr <- function(x) {
  x <- trimws(as.character(x)); x <- sub("\\s.*$", "", x); x <- sub("^chromosome", "", x, ignore.case = TRUE)
  x <- sub("^chr", "", x, ignore.case = TRUE); x[x %in% as.character(1:5)] <- paste0("Chr", x[x %in% as.character(1:5)])
  x
}

coalesce_columns <- function(data, candidates, label, required = TRUE) {
  columns <- candidates[candidates %in% names(data)]
  if (!length(columns)) {
    if (required) stop("Could not identify ", label, " column.\nExpected one of:\n  ", paste(candidates, collapse = ", "),
                       "\nAvailable columns:\n  ", paste(names(data), collapse = ", "), call. = FALSE)
    return(rep(NA_character_, nrow(data)))
  }
  output <- as.character(data[[columns[1]]])
  if (length(columns) > 1L) for (column in columns[-1]) {
    replacement <- as.character(data[[column]])
    use <- is.na(output) | trimws(output) == ""
    output[use] <- replacement[use]
  }
  output
}

###############################################################################
# 4. SIXTEEN FIXED COL-0 TERMINAL CANDIDATE WINDOWS
###############################################################################

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
if (nrow(onsen) != 16L || any(onsen$width_bp != 800L)) stop("ONSEN coordinate validation failed.", call. = FALSE)

###############################################################################
# 5. READ STRICT-TE BACKGROUND
###############################################################################

bg_raw <- data.table::fread(BG_FILE, data.table = FALSE, check.names = FALSE)

cat("\nBackground columns detected:\n", paste(names(bg_raw), collapse = ", "), "\n", sep = "")

bg_chr <- coalesce_columns(bg_raw, c("chr_clean", "seqid", "chromosome", "chr", "Chr", "seqnames", "sequence_name", "chrom"), "chromosome")
bg_start <- coalesce_columns(bg_raw, c("start", "Start", "start_1based", "region_start", "feature_start"), "start coordinate")
bg_end <- coalesce_columns(bg_raw, c("end", "End", "end_1based", "region_end", "feature_end"), "end coordinate")
bg_id <- coalesce_columns(bg_raw, c("background_id", "fasta_id", "region_id", "sequence_id", "candidate_id", "feature_id", "GeneID", "ID"), "region identifier", required = FALSE)
bg_family <- coalesce_columns(bg_raw, c("TE_family", "te_family", "family", "Alias", "alias", "transposon_family"), "TE family", required = FALSE)

bg <- data.frame(region_id = bg_id, chromosome = standard_chr(bg_chr), start = suppressWarnings(as.integer(bg_start)),
                 end = suppressWarnings(as.integer(bg_end)), TE_family = bg_family, stringsAsFactors = FALSE)

missing_id <- is.na(bg$region_id) | trimws(bg$region_id) == ""
bg$region_id[missing_id] <- paste0("strict_TE_", which(missing_id))
bg$region_id <- make.unique(bg$region_id)

invalid <- is.na(bg$chromosome) | is.na(bg$start) | is.na(bg$end) | bg$start < 1L | bg$end < bg$start
if (any(invalid)) {
  print(bg[invalid, , drop = FALSE])
  stop("Invalid strict-TE coordinates detected in ", sum(invalid), " rows.", call. = FALSE)
}

bg$width_bp <- bg$end - bg$start + 1L
bg$class <- BG_CLASS

cat("Strict-TE regions read: ", nrow(bg), "\n", sep = "")
cat("Chromosomes: ", paste(sort(unique(bg$chromosome)), collapse = ", "), "\n", sep = "")

###############################################################################
# 6. REMOVE DIRECT ONSEN OVERLAPS
###############################################################################

overlaps_onsen <- vapply(seq_len(nrow(bg)), function(i) {
  same_chr <- onsen$chromosome == bg$chromosome[i]
  any(same_chr & onsen$start <= bg$end[i] & onsen$end >= bg$start[i])
}, logical(1))

n_before <- nrow(bg); n_removed <- sum(overlaps_onsen)
bg <- bg[!overlaps_onsen, , drop = FALSE]

cat("Direct ONSEN-overlapping background regions removed: ", n_removed, "\n", sep = "")
cat("Final strict non-ONSEN TE regions: ", nrow(bg), "\n", sep = "")
if (nrow(bg) != 1942L) warning("Expected 1,942 background regions, but found ", nrow(bg), ".")

###############################################################################
# 7. READ GENOME AND EXTRACT SEQUENCES
###############################################################################

cat("\nReading genome...\n")
genome <- Biostrings::readDNAStringSet(GENOME_FILE)
names(genome) <- standard_chr(names(genome))

coordinates <- dplyr::bind_rows(
  onsen |> dplyr::select(region_id, chromosome, start, end, width_bp, class),
  bg |> dplyr::select(region_id, chromosome, start, end, width_bp, class)
)

missing_chr <- setdiff(unique(coordinates$chromosome), names(genome))
if (length(missing_chr)) stop("Chromosomes missing from genome FASTA: ", paste(missing_chr, collapse = ", "), call. = FALSE)

genome_lengths <- Biostrings::width(genome); names(genome_lengths) <- names(genome)
outside <- coordinates$end > genome_lengths[coordinates$chromosome]
if (any(outside)) stop("Coordinates extend beyond chromosome boundaries: ", paste(coordinates$region_id[outside], collapse = ", "))

extract_sequence <- function(chr, start, end) as.character(Biostrings::subseq(genome[[chr]], start = start, end = end))
coordinates$sequence <- mapply(extract_sequence, coordinates$chromosome, coordinates$start, coordinates$end, USE.NAMES = FALSE)

###############################################################################
# 8. READ JASPAR MOTIFS
###############################################################################

parse_jaspar <- function(path) {
  lines <- readLines(path, warn = FALSE); headers <- grep("^>", lines)
  if (!length(headers)) stop("No JASPAR headers found.")
  ends <- c(headers[-1] - 1L, length(lines)); motifs <- list()

  for (i in seq_along(headers)) {
    block <- lines[headers[i]:ends[i]]; header <- sub("^>\\s*", "", block[1]); parts <- strsplit(header, "\\s+")[[1]]
    motif_id <- parts[1]; motif_name <- if (length(parts) > 1L) paste(parts[-1], collapse = " ") else motif_id
    rows <- list()

    for (base in c("A", "C", "G", "T")) {
      line <- grep(paste0("^\\s*", base, "\\s*\\["), block, value = TRUE)
      if (!length(line)) next
      text <- sub(".*\\[", "", line[1]); text <- sub("\\].*", "", text); text <- gsub(",", " ", text)
      rows[[base]] <- as.numeric(strsplit(trimws(text), "\\s+")[[1]])
    }

    if (length(rows) == 4L) {
      pfm <- do.call(rbind, rows[c("A", "C", "G", "T")]); rownames(pfm) <- c("A", "C", "G", "T")
      motifs[[length(motifs) + 1L]] <- list(id = motif_id, name = motif_name, pfm = pfm)
    }
  }
  motifs
}

cat("Reading JASPAR motifs...\n")
motifs_all <- parse_jaspar(JASPAR_FILE)
motif_names <- toupper(vapply(motifs_all, function(x) sub("\\s+.*$", "", trimws(x$name)), character(1)))

missing_models <- setdiff(HSF_NAMES, motif_names)
if (length(missing_models)) stop("Missing HSF models: ", paste(missing_models, collapse = ", "), call. = FALSE)

hsf_motifs <- lapply(HSF_NAMES, function(x) motifs_all[[which(motif_names == x)[1]]])
names(hsf_motifs) <- HSF_NAMES

###############################################################################
# 9. PWM FUNCTIONS
###############################################################################

pfm_to_pwm <- function(pfm, pseudocount = 0.8) {
  adjusted <- pfm + pseudocount
  probabilities <- sweep(adjusted, 2, colSums(adjusted), "/")
  pwm <- log2(sweep(probabilities, 1, rep(0.25, 4), "/")); rownames(pwm) <- c("A", "C", "G", "T")
  list(pwm = pwm, width = ncol(pwm), minimum = sum(apply(pwm, 2, min)), maximum = sum(apply(pwm, 2, max)))
}

pwm_models <- lapply(hsf_motifs, function(x) {
  model <- pfm_to_pwm(x$pfm, PSEUDOCOUNT); model$id <- x$id; model$name <- x$name; model
})

score_one_strand <- function(sequence, model) {
  sequence <- toupper(sequence); bases <- strsplit(sequence, "", fixed = TRUE)[[1]]
  n_starts <- length(bases) - model$width + 1L
  if (n_starts < 1L) return(numeric())

  invalid_positions <- which(!bases %in% c("A", "C", "G", "T"))
  clean_bases <- bases; clean_bases[invalid_positions] <- "A"
  clean_sequence <- Biostrings::DNAString(paste0(clean_bases, collapse = ""))

  raw <- Biostrings::PWMscoreStartingAt(model$pwm, clean_sequence, starting.at = seq_len(n_starts))
  relative <- (raw - model$minimum) / (model$maximum - model$minimum)
  relative <- pmax(0, pmin(1, relative))

  if (length(invalid_positions)) {
    invalid_starts <- rep(FALSE, n_starts)
    for (position in invalid_positions) {
      first <- max(1L, position - model$width + 1L); last <- min(position, n_starts)
      if (first <= last) invalid_starts[first:last] <- TRUE
    }
    relative[invalid_starts] <- NA_real_
  }
  as.numeric(relative)
}

score_both_strands <- function(sequence, model) {
  reverse_sequence <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(toupper(sequence))))
  scores <- c(score_one_strand(sequence, model), score_one_strand(reverse_sequence, model))
  scores[is.finite(scores)]
}

top_n <- function(x, n = 5L) {
  x <- x[is.finite(x)]
  if (!length(x)) return(numeric())
  head(sort(x, decreasing = TRUE), n)
}

###############################################################################
# 10. SCAN ONE REGION
###############################################################################

scan_region <- function(region_id, sequence, region_class, chromosome, start, end, width_bp) {
  hit_counts <- integer(length(THRESHOLDS)); model_counts <- integer(length(THRESHOLDS))
  overall_top <- numeric(); model_maxima <- numeric(length(pwm_models))

  for (m in seq_along(pwm_models)) {
    scores <- score_both_strands(sequence, pwm_models[[m]])
    if (!length(scores)) { model_maxima[m] <- NA_real_; next }

    model_maxima[m] <- max(scores)
    overall_top <- top_n(c(overall_top, top_n(scores, 5L)), 5L)

    for (t in seq_along(THRESHOLDS)) {
      count <- sum(scores >= THRESHOLDS[t])
      hit_counts[t] <- hit_counts[t] + count
      if (count > 0L) model_counts[t] <- model_counts[t] + 1L
    }
  }

  data.frame(region_id = region_id, class = region_class, chromosome = chromosome, start = start, end = end, width_bp = width_bp,
             threshold = THRESHOLDS, HSF_motif_model_position_hits = hit_counts, unique_HSF_models = model_counts,
             HSF_hits_per_kb = hit_counts / width_bp * 1000,
             maximum_HSF_relative_score = if (length(overall_top)) max(overall_top) else NA_real_,
             mean_top5_HSF_relative_score = if (length(overall_top)) mean(overall_top) else NA_real_,
             median_per_model_maximum_score = if (any(is.finite(model_maxima))) median(model_maxima, na.rm = TRUE) else NA_real_)
}

###############################################################################
# 11. SCAN ALL REGIONS
###############################################################################

cat("\nScanning ", nrow(coordinates), " regions against ", length(pwm_models), " HSF models...\n", sep = "")

results <- vector("list", nrow(coordinates))
pb <- utils::txtProgressBar(min = 0, max = nrow(coordinates), style = 3)

for (i in seq_len(nrow(coordinates))) {
  results[[i]] <- scan_region(coordinates$region_id[i], coordinates$sequence[i], coordinates$class[i], coordinates$chromosome[i],
                              coordinates$start[i], coordinates$end[i], coordinates$width_bp[i])
  utils::setTxtProgressBar(pb, i)
}
close(pb)

region_metrics <- dplyr::bind_rows(results)
region_metrics$class <- factor(region_metrics$class, levels = c(BG_CLASS, ONSEN_CLASS))

continuous_metrics <- region_metrics |>
  dplyr::select(region_id, class, chromosome, start, end, width_bp, maximum_HSF_relative_score,
                mean_top5_HSF_relative_score, median_per_model_maximum_score) |>
  dplyr::distinct()

###############################################################################
# 12. STATISTICS
###############################################################################

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  sum(sign(outer(x, y, "-"))) / (length(x) * length(y))
}

wilcox_p <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value)
}

threshold_stats <- lapply(THRESHOLDS, function(cutoff) {
  current <- region_metrics[region_metrics$threshold == cutoff, ]
  x <- current$HSF_hits_per_kb[current$class == ONSEN_CLASS]
  y <- current$HSF_hits_per_kb[current$class == BG_CLASS]

  data.frame(relative_score_threshold = cutoff, n_ONSEN_regions = length(x), n_background_regions = length(y),
             ONSEN_median_hits_per_kb = median(x), ONSEN_Q1_hits_per_kb = unname(quantile(x, 0.25)),
             ONSEN_Q3_hits_per_kb = unname(quantile(x, 0.75)), background_median_hits_per_kb = median(y),
             background_Q1_hits_per_kb = unname(quantile(y, 0.25)), background_Q3_hits_per_kb = unname(quantile(y, 0.75)),
             Wilcoxon_P = wilcox_p(x, y), Cliffs_delta = cliffs_delta(x, y),
             background_regions_at_or_above_ONSEN_median = sum(y >= median(x)),
             background_regions_at_or_above_ONSEN_maximum = sum(y >= max(x)))
}) |> dplyr::bind_rows()

threshold_stats$Wilcoxon_P_BH <- p.adjust(threshold_stats$Wilcoxon_P, method = "BH")

metric_labels <- c(maximum_HSF_relative_score = "Maximum HSF relative PWM score",
                   mean_top5_HSF_relative_score = "Mean of top five HSF relative PWM scores",
                   median_per_model_maximum_score = "Median of ten per-model maximum scores")

continuous_stats <- lapply(names(metric_labels), function(metric) {
  x <- continuous_metrics[[metric]][continuous_metrics$class == ONSEN_CLASS]
  y <- continuous_metrics[[metric]][continuous_metrics$class == BG_CLASS]

  data.frame(metric = metric_labels[[metric]], metric_column = metric, n_ONSEN_regions = sum(is.finite(x)),
             n_background_regions = sum(is.finite(y)), ONSEN_median = median(x, na.rm = TRUE),
             ONSEN_Q1 = unname(quantile(x, 0.25, na.rm = TRUE)), ONSEN_Q3 = unname(quantile(x, 0.75, na.rm = TRUE)),
             background_median = median(y, na.rm = TRUE), background_Q1 = unname(quantile(y, 0.25, na.rm = TRUE)),
             background_Q3 = unname(quantile(y, 0.75, na.rm = TRUE)), Wilcoxon_P = wilcox_p(x, y),
             Cliffs_delta = cliffs_delta(x, y))
}) |> dplyr::bind_rows()

continuous_stats$Wilcoxon_P_BH <- p.adjust(continuous_stats$Wilcoxon_P, method = "BH")

###############################################################################
