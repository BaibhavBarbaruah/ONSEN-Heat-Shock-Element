# Col-0 non-stressed versus 24-h 37C heat-stress RNA-seq analysis.
# Covers Fig. 5 and source data for Table S10.
# Candidate-window signal is fractional and multimapping-aware.

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "tidyr", "readr", "stringr",
                   "purrr", "tibble", "ggplot2", "forcats", "DESeq2", "Rsubread"))

message_config()

run_alignment <- identical(tolower(Sys.getenv("ONSEN_RUN_ALIGNMENT", "false")), "true")
run_counting <- identical(tolower(Sys.getenv("ONSEN_RUN_FEATURECOUNTS", "false")), "true")
threads <- as.integer(Sys.getenv("ONSEN_THREADS", "8"))
RNA_LOG2FC_PSEUDOCOUNT_CPM <- 0.05

metadata_file <- find_any_input(c(
  "Col0_NS_vs_24h_37C_HS_sample_metadata.csv",
  "RNAseq_sample_metadata_template.csv"
))
sample_meta <- read.csv(metadata_file, check.names = FALSE)

rename_meta <- function(x) {
  rename_first <- function(target, alternatives) {
    hit <- alternatives[alternatives %in% names(x)][1]
    if (!is.na(hit) && target != hit) names(x)[names(x) == hit] <<- target
  }
  rename_first("sample_id", c("sample_id", "sample", "Sample"))
  rename_first("treatment", c("treatment", "condition", "group"))
  rename_first("bam_filename", c("bam_filename", "bam", "BAM"))
  rename_first("fastq_R1", c("fastq_R1", "R1", "read1"))
  rename_first("fastq_R2", c("fastq_R2", "R2", "read2"))
  assert_columns(x, c("sample_id", "treatment"), "RNA-seq sample metadata")
  x$treatment <- ifelse(grepl("HS|heat|37", x$treatment, ignore.case = TRUE), "HS", "NS")
  x$treatment <- factor(x$treatment, levels = c("NS", "HS"))
  x
}
sample_meta <- rename_meta(sample_meta)

# All six libraries were generated in our laboratory. Validate the exact
# three-replicate-per-condition BioSample/BioProject map before any analysis.
assert_columns(
  sample_meta,
  c("sample_id", "treatment", "replicate", "ecotype",
    "biosample_accession", "bioproject", "laboratory_provenance"),
  "RNA-seq sample metadata with DDBJ provenance"
)
expected_accessions <- data.frame(
  sample_id = c(
    "Col0_NS_rep1", "Col0_NS_rep2", "Col0_NS_rep3",
    "Col0_HS_24h_37C_rep1", "Col0_HS_24h_37C_rep2", "Col0_HS_24h_37C_rep3"
  ),
  treatment = c("NS", "NS", "NS", "HS", "HS", "HS"),
  replicate = c(1L, 2L, 3L, 1L, 2L, 3L),
  biosample_accession = c(
    "SAMD01789795", "SAMD01789796", "SAMD01789797",
    "SAMD01943917", "SAMD01943918", "SAMD01943919"
  ),
  bioproject = c(
    "PRJDB39904", "PRJDB39904", "PRJDB39904",
    "PRJDB42759", "PRJDB42759", "PRJDB42759"
  ),
  stringsAsFactors = FALSE
)
idx <- match(expected_accessions$sample_id, sample_meta$sample_id)
if (nrow(sample_meta) != 6L || anyNA(idx)) {
  stop("RNA-seq metadata must contain exactly the six validated Col-0 libraries.", call. = FALSE)
}
for (field in c("treatment", "replicate", "biosample_accession", "bioproject")) {
  observed <- as.character(sample_meta[[field]][idx])
  expected <- as.character(expected_accessions[[field]])
  if (!identical(observed, expected)) {
    stop("Incorrect DDBJ Col-0 assignment in field: ", field, call. = FALSE)
  }
}
if (!all(tolower(sample_meta$ecotype[idx]) %in% c("col-0", "col0")) ||
    !all(sample_meta$laboratory_provenance[idx] == "Generated in our laboratory")) {
  stop("RNA-seq metadata contains a non-Col-0 library or incorrect laboratory provenance.", call. = FALSE)
}

# --------------------------- Optional alignment -------------------------------
# Laboratory-generated paired-end read files or derived BAMs are resolved from the validated metadata.
# Alignment is disabled unless explicitly
# requested through ONSEN_RUN_ALIGNMENT=true.
if (run_alignment) {
  genome_fasta <- find_input("Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz")
  index_base <- out_file("TAIR10_subread_index")
  if (!file.exists(paste0(index_base, ".00.b.array"))) {
    Rsubread::buildindex(
      basename = index_base,
      reference = genome_fasta,
      memory = 8000
    )
  }

  assert_columns(sample_meta, c("fastq_R1", "fastq_R2"), "RNA-seq metadata for alignment")
  for (i in seq_len(nrow(sample_meta))) {
    read1 <- find_input(sample_meta$fastq_R1[[i]])
    read2 <- find_input(sample_meta$fastq_R2[[i]])
    bam_out <- out_file(paste0(sample_meta$sample_id[[i]], ".bam"))
    Rsubread::subjunc(
      index = index_base,
      readfile1 = read1,
      readfile2 = read2,
      output_file = bam_out,
      nthreads = threads,
      unique = FALSE,
      nBestLocations = 20
    )
  }
}

# Resolve BAM files from metadata or exact project filenames.
resolve_bam <- function(sample_id, bam_filename = NA_character_) {
  alternatives <- unique(na.omit(c(
    bam_filename,
    paste0(sample_id, ".bam"),
    sub("_24h_", "_", paste0(sample_id, ".bam"))
  )))
  find_any_input(alternatives, required = FALSE)
}
sample_meta$bam_path <- mapply(
  resolve_bam,
  sample_meta$sample_id,
  if ("bam_filename" %in% names(sample_meta)) sample_meta$bam_filename else NA_character_,
  USE.NAMES = FALSE
)

# ------------------------- Optional featureCounts -----------------------------
gene_count_file <- find_any_input(c(
  "Col0_NS_vs_24h_37C_HS_gene_counts_raw.tsv",
  "gene_counts.tsv",
  "Col0_NS_vs_24h_37C_HS_gene_counts_annotation.tsv"
), required = FALSE)

candidate_count_file <- find_any_input(c(
  "Col0_NS_vs_24h_37C_HS_candidate_ONSEN_HSF_rich_TE_loci_fractional_multimap_counts_annotated.tsv",
  "Col0_NS_vs_24h_37C_HS_candidate_ONSEN_HSF_rich_TE_loci_fractional_multimap_counts.tsv",
  "candidate_window_fractional_counts_annotated.tsv"
), required = FALSE)

if (run_counting || is.na(gene_count_file) || is.na(candidate_count_file)) {
  if (!all(file.exists(sample_meta$bam_path))) {
    stop(
      "BAM-to-count processing was requested or count files are absent, but one or more BAMs were not found.\n",
      paste(sample_meta$sample_id[!file.exists(sample_meta$bam_path)], collapse = ", ")
    )
  }

  gtf_file <- find_any_input(c(
    "Arabidopsis_thaliana.TAIR10.54.gtf",
    "Arabidopsis_thaliana.TAIR10.gtf"
  ))
  candidate_saf <- find_any_input(c(
    "candidate_ONSEN_HSF_rich_TE_loci_for_featureCounts.saf",
    "ONSEN_Col0_terminal_candidate_windows.saf"
  ))

  gene_fc <- Rsubread::featureCounts(
    files = sample_meta$bam_path,
    annot.ext = gtf_file,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    useMetaFeatures = TRUE,
    isPairedEnd = TRUE,
    requireBothEndsMapped = TRUE,
    countMultiMappingReads = FALSE,
    allowMultiOverlap = FALSE,
    nthreads = threads
  )
  gene_counts <- as.data.frame(gene_fc$counts)
  names(gene_counts) <- sample_meta$sample_id
  gene_counts <- tibble::rownames_to_column(gene_counts, "gene_id")
  gene_count_file <- out_file("Col0_NS_vs_24h_37C_HS_gene_counts_raw_repository.tsv")
  readr::write_tsv(gene_counts, gene_count_file)

  candidate_fc <- Rsubread::featureCounts(
    files = sample_meta$bam_path,
    annot.ext = candidate_saf,
    isGTFAnnotationFile = FALSE,
    useMetaFeatures = FALSE,
    isPairedEnd = TRUE,
    requireBothEndsMapped = TRUE,
    countMultiMappingReads = TRUE,
    fraction = TRUE,
    allowMultiOverlap = TRUE,
    nthreads = threads
  )
  candidate_counts <- as.data.frame(candidate_fc$counts)
  names(candidate_counts) <- sample_meta$sample_id
  candidate_counts <- tibble::rownames_to_column(candidate_counts, "candidate_id")
  candidate_count_file <- out_file(
    "Col0_NS_vs_24h_37C_HS_candidate_fractional_counts_repository.tsv"
  )
  readr::write_tsv(candidate_counts, candidate_count_file)

  library_summary <- data.frame(
    sample_id = sample_meta$sample_id,
    assigned_gene_counts = colSums(gene_fc$counts),
    stringsAsFactors = FALSE
  )
  safe_write_csv(library_summary, "Col0_NS_vs_24h_37C_HS_gene_library_summary_repository.csv")
}

# ----------------------------- Gene-level DESeq2 ------------------------------
normalize_count_matrix <- function(path, id_name) {
  x <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  # featureCounts output may include annotation columns.
  id_col <- c(id_name, "Geneid", "gene_id", "candidate_id", names(x)[1])[
    c(id_name, "Geneid", "gene_id", "candidate_id", names(x)[1]) %in% names(x)
  ][1]
  sample_columns <- intersect(sample_meta$sample_id, names(x))
  if (!length(sample_columns)) {
    # Match BAM basenames when featureCounts retained full BAM paths.
    for (j in seq_along(names(x))) {
      base <- tools::file_path_sans_ext(basename(names(x)[[j]]))
      hit <- match(base, sample_meta$sample_id)
      if (!is.na(hit)) names(x)[[j]] <- sample_meta$sample_id[[hit]]
    }
    sample_columns <- intersect(sample_meta$sample_id, names(x))
  }
  if (!length(sample_columns)) {
    stop("No count columns matched RNA-seq sample IDs in: ", path)
  }
  result <- x[, c(id_col, sample_columns), drop = FALSE]
  names(result)[1] <- id_name
  result
}

gene_counts <- normalize_count_matrix(gene_count_file, "gene_id")
count_matrix <- as.matrix(gene_counts[, sample_meta$sample_id, drop = FALSE])
rownames(count_matrix) <- gene_counts$gene_id
mode(count_matrix) <- "integer"

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = data.frame(
    treatment = sample_meta$treatment,
    row.names = sample_meta$sample_id
  ),
  design = ~ treatment
)
dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, contrast = c("treatment", "HS", "NS"))
res_table <- as.data.frame(res) |>
  tibble::rownames_to_column("gene_id") |>
  dplyr::arrange(padj)
safe_write_csv(res_table, "DESeq2_Col0_HS_24h_37C_vs_NS_all_genes_repository.csv")

vst_object <- DESeq2::vst(dds, blind = FALSE)
vst_matrix <- DESeq2::assay(vst_object)
vst_table <- as.data.frame(vst_matrix) |>
  tibble::rownames_to_column("gene_id")
safe_write_csv(vst_table, "DESeq2_Col0_NS_vs_24h_37C_HS_VST_matrix_repository.csv")

# Prefer the exact five-marker plot data already generated for the manuscript.
exact_marker_plot_file <- find_input(
  "Fig5A_5gene_heat_marker_2deltaVST_plot_data.tsv",
  required = FALSE
)
selected_marker_file <- find_input(
  "Fig5A_selected_5_heat_markers.tsv",
  required = FALSE
)

if (!is.na(exact_marker_plot_file) && !ONSEN_FORCE_RESCAN) {
  marker_points <- read_table_auto(exact_marker_plot_file)
} else {
  if (is.na(selected_marker_file)) {
    stop(
      "Exact heat-marker selection file was not found. Expected: ",
      "Fig5A_selected_5_heat_markers.tsv"
    )
  }
  selected <- read_table_auto(selected_marker_file)
  gene_id_col <- c("gene_id", "GeneID", "TAIR_ID")[c("gene_id", "GeneID", "TAIR_ID") %in% names(selected)][1]
  symbol_col <- c("gene_symbol", "symbol", "gene_name")[c("gene_symbol", "symbol", "gene_name") %in% names(selected)][1]
  if (is.na(gene_id_col)) stop("Selected heat-marker table lacks a gene ID column.")
  if (is.na(symbol_col)) selected$gene_symbol <- selected[[gene_id_col]] else selected$gene_symbol <- selected[[symbol_col]]

  marker_points <- vst_table |>
    dplyr::filter(gene_id %in% selected[[gene_id_col]]) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(sample_meta$sample_id),
      names_to = "sample_id", values_to = "vst"
    ) |>
    dplyr::left_join(
      sample_meta |> dplyr::select(sample_id, treatment),
      by = "sample_id"
    ) |>
    dplyr::left_join(
      selected |>
        dplyr::transmute(gene_id = .data[[gene_id_col]], gene_symbol),
      by = "gene_id"
    ) |>
    dplyr::group_by(gene_id) |>
    dplyr::mutate(
      NS_mean_vst = mean(vst[treatment == "NS"], na.rm = TRUE),
      relative_expression_2deltaVST = 2^(vst - NS_mean_vst)
    ) |>
    dplyr::ungroup()
}
safe_write_csv(marker_points, "Fig5A_heat_marker_plot_data_repository.csv")

# Normalize marker columns.
if (!"gene_symbol" %in% names(marker_points)) {
  symbol_col <- c("symbol", "gene_name")[c("symbol", "gene_name") %in% names(marker_points)][1]
  marker_points$gene_symbol <- if (!is.na(symbol_col)) marker_points[[symbol_col]] else marker_points$gene_id
}
if (!"treatment" %in% names(marker_points)) {
  marker_points <- marker_points |>
    dplyr::left_join(sample_meta |> dplyr::select(sample_id, treatment), by = "sample_id")
}
value_col <- c(
  "relative_expression_2deltaVST", "relative_expression",
  "two_delta_VST", "expression_2deltaVST"
)[c(
  "relative_expression_2deltaVST", "relative_expression",
  "two_delta_VST", "expression_2deltaVST"
) %in% names(marker_points)][1]
if (is.na(value_col)) stop("Heat-marker plot data lacks a relative-expression column.")

marker_hs <- marker_points |>
  dplyr::filter(treatment == "HS") |>
  dplyr::group_by(gene_id, gene_symbol) |>
  dplyr::summarise(
    mean_relative_expression = mean(.data[[value_col]], na.rm = TRUE),
    sem_relative_expression = sem(.data[[value_col]]),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    res_table |> dplyr::select(gene_id, padj),
    by = "gene_id"
  ) |>
  dplyr::mutate(
    significance = dplyr::case_when(
      padj <= 1e-4 ~ "****",
      padj <= 1e-3 ~ "***",
      padj <= 1e-2 ~ "**",
      padj <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
safe_write_csv(marker_hs, "Fig5A_heat_marker_summary_repository.csv")

# ----------------------- Candidate-window signal analysis --------------------
candidate_counts <- normalize_count_matrix(candidate_count_file, "candidate_id")

candidate_key_file <- find_any_input(c(
  "candidate_ONSEN_HSF_rich_TE_loci_key.tsv",
  "Col0_NS_vs_24h_37C_HS_candidate_TE_loci_featureCounts_annotation.tsv"
), required = FALSE)
if (!is.na(candidate_key_file)) {
  candidate_key <- read_table_auto(candidate_key_file)
} else {
  candidate_key <- data.frame(candidate_id = candidate_counts$candidate_id)
}
if (!"candidate_id" %in% names(candidate_key)) {
  names(candidate_key)[1] <- "candidate_id"
}

class_col <- c("candidate_class", "locus_class", "class")[
  c("candidate_class", "locus_class", "class") %in% names(candidate_key)
][1]
if (is.na(class_col)) {
  candidate_key$candidate_class <- ifelse(
    grepl("ONSEN|ATCOPIA78", candidate_key$candidate_id, ignore.case = TRUE),
    "ONSEN LTR", "HSF-rich non-ONSEN TE"
  )
} else {
  candidate_key$candidate_class <- candidate_key[[class_col]]
}

library_file <- find_any_input(c(
  "Col0_NS_vs_24h_37C_HS_gene_assigned_count_library_summary.tsv",
  "Col0_NS_vs_24h_37C_HS_gene_library_summary_repository.csv",
  "gene_assigned_count_library_summary.tsv"
), required = FALSE)

if (!is.na(library_file)) {
  library_summary <- read_table_auto(library_file)
  sample_col <- c("sample_id", "sample", "Sample")[c("sample_id", "sample", "Sample") %in% names(library_summary)][1]
  assigned_col <- names(library_summary)[
    grepl("assigned.*gene|gene.*assigned|assigned_counts|assigned", names(library_summary), ignore.case = TRUE)
  ][1]
  if (is.na(sample_col) || is.na(assigned_col)) {
    stop("Could not identify sample/assigned-count columns in library summary.")
  }
  library_summary <- library_summary |>
    dplyr::transmute(
      sample_id = .data[[sample_col]],
      assigned_gene_counts = as.numeric(.data[[assigned_col]])
    )
} else {
  library_summary <- data.frame(
    sample_id = sample_meta$sample_id,
    assigned_gene_counts = colSums(count_matrix),
    stringsAsFactors = FALSE
  )
}

candidate_long <- candidate_counts |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_meta$sample_id),
    names_to = "sample_id", values_to = "fractional_count"
  ) |>
  dplyr::left_join(sample_meta |> dplyr::select(sample_id, treatment), by = "sample_id") |>
  dplyr::left_join(library_summary, by = "sample_id") |>
  dplyr::left_join(
    candidate_key |> dplyr::select(candidate_id, candidate_class),
    by = "candidate_id"
  ) |>
  dplyr::mutate(
    candidate_CPM = as.numeric(fractional_count) / assigned_gene_counts * 1e6,
    candidate_class = dplyr::coalesce(candidate_class, "Candidate")
  )
safe_write_csv(candidate_long, "Fig5_candidate_window_signal_long_repository.csv")

class_replicate <- candidate_long |>
  dplyr::group_by(candidate_class, treatment, sample_id) |>
  dplyr::summarise(
    mean_CPM_per_locus = mean(candidate_CPM, na.rm = TRUE),
    .groups = "drop"
  )
safe_write_csv(class_replicate, "Fig5B_candidate_window_class_replicate_summary_repository.csv")

class_stats <- class_replicate |>
  dplyr::group_by(candidate_class) |>
  dplyr::summarise(
    Welch_P = tryCatch(
      stats::t.test(mean_CPM_per_locus ~ treatment, var.equal = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    Wilcoxon_P = tryCatch(
      stats::wilcox.test(mean_CPM_per_locus ~ treatment, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    NS_mean = mean(mean_CPM_per_locus[treatment == "NS"], na.rm = TRUE),
    HS_mean = mean(mean_CPM_per_locus[treatment == "HS"], na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    log2FC_HS_vs_NS = log2(
      (HS_mean + RNA_LOG2FC_PSEUDOCOUNT_CPM) /
        (NS_mean + RNA_LOG2FC_PSEUDOCOUNT_CPM)
    ),
    fold_from_log2FC = 2^log2FC_HS_vs_NS
  )
safe_write_csv(class_stats, "Fig5B_candidate_window_class_statistics_repository.csv")

window_summary <- candidate_long |>
  dplyr::group_by(candidate_id, candidate_class, treatment) |>
  dplyr::summarise(mean_CPM = mean(candidate_CPM, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = treatment, values_from = mean_CPM) |>
  dplyr::mutate(
    HS_signal_log2 = log2(HS + 1),
    log2FC_HS_vs_NS = log2((HS + RNA_LOG2FC_PSEUDOCOUNT_CPM) /
                              (NS + RNA_LOG2FC_PSEUDOCOUNT_CPM))
  ) |>
  dplyr::arrange(dplyr::desc(HS_signal_log2))
safe_write_csv(window_summary, "Fig5C_individual_candidate_window_signal_repository.csv")

# -------------------------------- Figures ------------------------------------
if (ONSEN_MAKE_FIGURES) {
  marker_hs$gene_symbol <- factor(
    marker_hs$gene_symbol,
    levels = marker_hs$gene_symbol
  )
  p5a <- ggplot2::ggplot(
    marker_hs,
    ggplot2::aes(gene_symbol, mean_relative_expression, fill = gene_symbol)
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_col(width = 0.67, colour = "black", show.legend = FALSE) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_relative_expression - sem_relative_expression,
        ymax = mean_relative_expression + sem_relative_expression
      ),
      width = 0.18
    ) +
    ggplot2::geom_point(
      data = marker_points |> dplyr::filter(treatment == "HS"),
      ggplot2::aes(
        x = gene_symbol, y = .data[[value_col]]
      ),
      inherit.aes = FALSE,
      position = ggplot2::position_jitter(width = 0.08),
      shape = 21, fill = "white", size = 2.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        y = mean_relative_expression + sem_relative_expression +
          max(mean_relative_expression, na.rm = TRUE) * 0.06,
        label = significance
      ),
      fontface = "bold", size = 4
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Relative expression vs NS baseline (2^delta-VST)"
    ) +
    theme_onsen(13) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_plot_pair(p5a, "Fig5A_heat_marker_expression", 6.2, 4.8)

  p5b <- ggplot2::ggplot(
    class_replicate,
    ggplot2::aes(candidate_class, mean_CPM_per_locus, fill = treatment)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.72),
      width = 0.66, colour = "black"
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitterdodge(
        jitter.width = 0.08, dodge.width = 0.72
      ),
      shape = 21, fill = "white", size = 2.3
    ) +
    ggplot2::scale_fill_manual(values = c(NS = "#B9CEDF", HS = "#EAA0BF")) +
    ggplot2::labs(
      x = NULL, y = "Candidate-window signal", fill = "Treatment"
    ) +
    theme_onsen(13) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
    )
  save_plot_pair(p5b, "Fig5B_candidate_window_class_signal", 6.5, 4.8)

  p5c <- ggplot2::ggplot(
    window_summary,
    ggplot2::aes(
      HS_signal_log2,
      forcats::fct_reorder(candidate_id, HS_signal_log2),
      colour = candidate_class,
      size = pmax(0, log2FC_HS_vs_NS)
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::labs(
      x = "HS candidate-window signal [log2(CPM + 1)]",
      y = NULL, colour = "Candidate class",
      size = "HS/NS\nlog2FC"
    ) +
    theme_onsen(11)
  save_plot_pair(p5c, "Fig5C_individual_candidate_window_signal", 8.2, 6.0)
}

message("RNA-seq analysis completed.")
