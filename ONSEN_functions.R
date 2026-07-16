# Shared utilities for the flat ONSEN-HSE repository.

source("ONSEN_config.R")

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) {
    stop(
      "Missing required R package(s): ", paste(missing, collapse = ", "),
      "\nRun 00_install_packages.R first.",
      call. = FALSE
    )
  }
}

standard_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "Chr", x, ignore.case = TRUE)
  x <- sub("^([1-5])$", "Chr\\1", x)
  x
}

clean_dna <- function(x) {
  x <- toupper(gsub("\\s+", "", as.character(x)))
  if (any(!grepl("^[ACGTN]+$", x))) {
    warning("One or more sequences contain characters outside A/C/G/T/N.")
  }
  x
}

rev_comp <- function(sequence) {
  sequence <- clean_dna(sequence)
  chars <- strsplit(sequence, "", fixed = TRUE)[[1]]
  complement <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
  paste(rev(unname(complement[chars])), collapse = "")
}

split_bases <- function(sequence) {
  strsplit(clean_dna(sequence), "", fixed = TRUE)[[1]]
}

gc_percent <- function(sequence) {
  x <- split_bases(sequence)
  100 * sum(x %in% c("G", "C")) / length(x)
}

hamming_distance <- function(a, b, ignore_n = FALSE) {
  aa <- split_bases(a)
  bb <- split_bases(b)
  if (length(aa) != length(bb)) stop("Hamming distance requires equal-length sequences.")
  keep <- rep(TRUE, length(aa))
  if (ignore_n) keep <- aa != "N" & bb != "N"
  sum(aa[keep] != bb[keep])
}

contains_canonical_hse_unit <- function(sequence) {
  sequence <- clean_dna(sequence)
  if (nchar(sequence) < 5) return(FALSE)
  windows <- substring(sequence, seq_len(nchar(sequence) - 4L), seq_len(nchar(sequence) - 4L) + 4L)
  cores <- substring(windows, 2, 4)
  any(cores %in% c("TTC", "GAA"))
}

find_canonical_hse_units <- function(sequence) {
  sequence <- clean_dna(sequence)
  if (nchar(sequence) < 5) return(data.frame())
  starts <- seq_len(nchar(sequence) - 4L)
  windows <- substring(sequence, starts, starts + 4L)
  cores <- substring(windows, 2, 4)
  keep <- cores %in% c("TTC", "GAA")
  data.frame(
    start = starts[keep],
    end = starts[keep] + 4L,
    sequence = windows[keep],
    consensus_type = ifelse(cores[keep] == "TTC", "nTTCn", "nGAAn"),
    stringsAsFactors = FALSE
  )
}

read_fasta_simple <- function(path) {
  require_packages("Biostrings")
  dna <- Biostrings::readDNAStringSet(path)
  setNames(as.character(dna), names(dna))
}

parse_jaspar_pfms <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header_idx <- grep("^>", lines)
  if (!length(header_idx)) stop("No JASPAR headers found in: ", path)
  ends <- c(header_idx[-1] - 1L, length(lines))
  motifs <- list()

  for (i in seq_along(header_idx)) {
    block <- lines[header_idx[i]:ends[i]]
    header <- sub("^>\\s*", "", block[[1]])
    parts <- strsplit(header, "\\s+")[[1]]
    motif_id <- parts[[1]]
    motif_name <- if (length(parts) > 1L) paste(parts[-1], collapse = " ") else motif_id

    rows <- list()
    for (base in c("A", "C", "G", "T")) {
      line <- grep(paste0("^\\s*", base, "\\s*\\["), block, value = TRUE)
      if (!length(line)) next
      numbers <- sub(".*\\[", "", line[[1]])
      numbers <- sub("\\].*", "", numbers)
      values <- as.numeric(strsplit(trimws(gsub(",", " ", numbers)), "\\s+")[[1]])
      rows[[base]] <- values
    }
    if (length(rows) != 4L) next
    pfm <- do.call(rbind, rows[c("A", "C", "G", "T")])
    rownames(pfm) <- c("A", "C", "G", "T")
    motifs[[length(motifs) + 1L]] <- list(id = motif_id, name = motif_name, pfm = pfm)
  }
  motifs
}

pfm_to_pwm <- function(pfm, pseudocount = 0.8, background = rep(0.25, 4)) {
  background <- background / sum(background)
  adjusted <- pfm + pseudocount
  probabilities <- sweep(adjusted, 2, colSums(adjusted), "/")
  pwm <- log2(sweep(probabilities, 1, background, "/"))
  rownames(pwm) <- c("A", "C", "G", "T")
  list(
    pwm = pwm,
    min_score = sum(apply(pwm, 2, min)),
    max_score = sum(apply(pwm, 2, max))
  )
}

score_pwm_window <- function(window, pwm) {
  bases <- split_bases(window)
  idx <- match(bases, rownames(pwm))
  if (anyNA(idx)) return(NA_real_)
  sum(pwm[cbind(idx, seq_along(idx))])
}

scan_one_motif <- function(sequence_id, sequence, motif, threshold = 0.85,
                           pseudocount = 0.8, retain_all = FALSE) {
  sequence <- clean_dna(sequence)
  model <- pfm_to_pwm(motif$pfm, pseudocount = pseudocount)
  width <- ncol(model$pwm)
  sequence_width <- nchar(sequence)
  if (width > sequence_width) return(data.frame())

  rows <- list()
  counter <- 0L
  strands <- c("+" = sequence, "-" = rev_comp(sequence))

  for (strand in names(strands)) {
    scan_sequence <- strands[[strand]]
    for (scan_start in seq_len(sequence_width - width + 1L)) {
      matched <- substring(scan_sequence, scan_start, scan_start + width - 1L)
      raw <- score_pwm_window(matched, model$pwm)
      relative <- (raw - model$min_score) / (model$max_score - model$min_score)

      if (strand == "+") {
        forward_start <- scan_start
      } else {
        forward_start <- sequence_width - scan_start - width + 2L
      }
      forward_end <- forward_start + width - 1L

      if (retain_all || (!is.na(relative) && relative >= threshold)) {
        counter <- counter + 1L
        rows[[counter]] <- data.frame(
          sequence_id = sequence_id,
          motif_id = motif$id,
          motif_name = motif$name,
          strand = strand,
          scan_start = scan_start,
          forward_start = forward_start,
          forward_end = forward_end,
          motif_width = width,
          matched_sequence = matched,
          raw_score = raw,
          relative_score = relative,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

best_score_one_motif <- function(sequence_id, sequence, motif, pseudocount = 0.8) {
  x <- scan_one_motif(
    sequence_id, sequence, motif,
    threshold = -Inf, pseudocount = pseudocount, retain_all = TRUE
  )
  if (!nrow(x)) return(data.frame())
  x[which.max(x$relative_score), , drop = FALSE]
}

classify_tf_family <- function(motif_name) {
  x <- toupper(as.character(motif_name))
  ifelse(grepl("HSF", x), "HSF",
  ifelse(grepl("AP2|ERF|DREB|RAP2", x), "AP2/ERF",
  ifelse(grepl("DOF", x), "DOF",
  ifelse(grepl("GATA", x), "GATA",
  ifelse(grepl("MYB", x), "MYB",
  ifelse(grepl("NAC", x), "NAC",
  ifelse(grepl("BZIP|GBF|ABF|TGA", x), "bZIP",
  ifelse(grepl("WRKY", x), "WRKY",
  ifelse(grepl("C2H2|ZAT|AZF", x), "C2H2/ZAT",
  ifelse(grepl("HD-ZIP|ATHB", x), "HD-ZIP",
  ifelse(grepl("ARF", x), "ARF",
  ifelse(grepl("LBD", x), "LBD", "Other"))))))))))))
}

scan_sequences_against_motifs <- function(sequence_table, motifs,
                                          threshold = 0.85,
                                          pseudocount = 0.8,
                                          retain_all = FALSE) {
  require_packages(c("data.table"))
  all_rows <- list()
  index <- 0L
  for (i in seq_len(nrow(sequence_table))) {
    for (motif in motifs) {
      index <- index + 1L
      all_rows[[index]] <- scan_one_motif(
        sequence_table$sequence_id[[i]],
        sequence_table$sequence[[i]],
        motif,
        threshold = threshold,
        pseudocount = pseudocount,
        retain_all = retain_all
      )
    }
  }
  x <- data.table::rbindlist(all_rows, fill = TRUE)
  if (nrow(x)) x[, tf_family := classify_tf_family(motif_name)]
  x
}

classify_motif_effect <- function(native_score, designed_score,
                                  threshold = 0.85, effect_margin = 0.05,
                                  native_position = NA_integer_,
                                  designed_position = NA_integer_) {
  native_high <- is.finite(native_score) && native_score >= threshold
  designed_high <- is.finite(designed_score) && designed_score >= threshold
  if (native_high && !designed_high) return("lost")
  if (!native_high && designed_high) return("gained")
  if (native_high && designed_high) {
    delta <- designed_score - native_score
    if (is.finite(native_position) && is.finite(designed_position) &&
        native_position != designed_position && abs(delta) < effect_margin) return("shifted")
    if (delta >= effect_margin) return("strengthened")
    if (delta <= -effect_margin) return("weakened")
    return("retained")
  }
  "below_threshold_in_both"
}

assert_columns <- function(x, required, object_name = deparse(substitute(x))) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(
      "Missing columns in ", object_name, ": ", paste(missing, collapse = ", "),
      "\nAvailable columns: ", paste(names(x), collapse = ", "),
      call. = FALSE
    )
  }
}

read_table_auto <- function(path) {
  require_packages("data.table")
  data.table::fread(path, data.table = FALSE, check.names = FALSE)
}

safe_write_csv <- function(x, filename, subdir = NULL) {
  path <- out_file(filename, subdir)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  comparisons <- outer(x, y, FUN = "-")
  (sum(comparisons > 0) - sum(comparisons < 0)) / (length(x) * length(y))
}

merge_overlapping_intervals <- function(x, group_col = "sequence_id",
                                        start_col = "forward_start",
                                        end_col = "forward_end") {
  assert_columns(x, c(group_col, start_col, end_col), "interval table")
  if (!nrow(x)) return(x)
  pieces <- split(x, x[[group_col]])
  result <- lapply(pieces, function(z) {
    z <- z[order(z[[start_col]], z[[end_col]]), , drop = FALSE]
    cluster <- integer(nrow(z))
    current <- 1L
    current_end <- z[[end_col]][[1]]
    cluster[[1]] <- current
    if (nrow(z) > 1L) {
      for (i in 2:nrow(z)) {
        if (z[[start_col]][[i]] <= current_end) {
          cluster[[i]] <- current
          current_end <- max(current_end, z[[end_col]][[i]], na.rm = TRUE)
        } else {
          current <- current + 1L
          cluster[[i]] <- current
          current_end <- z[[end_col]][[i]]
        }
      }
    }
    z$cluster_id <- cluster
    z
  })
  do.call(rbind, result)
}

summarize_interval_clusters <- function(clustered,
                                        group_col = "sequence_id",
                                        start_col = "forward_start",
                                        end_col = "forward_end") {
  if (!nrow(clustered)) return(data.frame())
  key <- interaction(clustered[[group_col]], clustered$cluster_id, drop = TRUE)
  pieces <- split(clustered, key)
  result <- lapply(pieces, function(z) {
    data.frame(
      sequence_id = z[[group_col]][[1]],
      cluster_id = z$cluster_id[[1]],
      merged_start = min(z[[start_col]], na.rm = TRUE),
      merged_end = max(z[[end_col]], na.rm = TRUE),
      merged_width = max(z[[end_col]], na.rm = TRUE) - min(z[[start_col]], na.rm = TRUE) + 1L,
      n_model_position_hits = nrow(z),
      n_unique_models = length(unique(z$motif_id)),
      motif_models = paste(sort(unique(z$motif_name)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}

weighted_methylation <- function(methylated, total) {
  if (!sum(total, na.rm = TRUE)) return(NA_real_)
  100 * sum(methylated, na.rm = TRUE) / sum(total, na.rm = TRUE)
}

sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

theme_onsen <- function(base_size = 13) {
  require_packages("ggplot2")
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(face = "bold", colour = "black"),
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold", colour = "black")
    )
}

save_plot_pair <- function(plot, stem, width, height, dpi = 600) {
  require_packages("ggplot2")
  ggplot2::ggsave(out_file(paste0(stem, ".pdf")), plot, width = width, height = height, units = "in")
  ggplot2::ggsave(out_file(paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = dpi)
}

write_workbook_sheet <- function(workbook, sheet_name, title, data,
                                 p_value_columns = character()) {
  require_packages("openxlsx")
  openxlsx::addWorksheet(workbook, sheet_name)
  openxlsx::writeData(workbook, sheet_name, title, startRow = 1, startCol = 1)
  openxlsx::writeData(workbook, sheet_name, data, startRow = 2, startCol = 1, withFilter = TRUE)
  title_style <- openxlsx::createStyle(fontSize = 12, textDecoration = "bold", wrapText = TRUE)
  header_style <- openxlsx::createStyle(fontSize = 11, textDecoration = "bold",
                                        fgFill = "#D9EAF7", border = "Bottom",
                                        halign = "center", valign = "center", wrapText = TRUE)
  openxlsx::addStyle(workbook, sheet_name, title_style, rows = 1, cols = seq_len(max(1, ncol(data))), gridExpand = TRUE)
  openxlsx::addStyle(workbook, sheet_name, header_style, rows = 2, cols = seq_len(max(1, ncol(data))), gridExpand = TRUE)
  openxlsx::freezePane(workbook, sheet_name, firstActiveRow = 3)
  openxlsx::setColWidths(workbook, sheet_name, cols = seq_len(max(1, ncol(data))), widths = "auto")
  if (length(p_value_columns)) {
    for (column in intersect(p_value_columns, names(data))) {
      j <- match(column, names(data))
      openxlsx::addStyle(
        workbook, sheet_name,
        openxlsx::createStyle(numFmt = "0.000E+00"),
        rows = 3:(nrow(data) + 2L), cols = j, gridExpand = TRUE
      )
    }
  }
  invisible(workbook)
}

write_session_info <- function(filename = "sessionInfo.txt") {
  path <- out_file(filename)
  sink(path)
  on.exit(sink(), add = TRUE)
  print(sessionInfo())
  invisible(path)
}
