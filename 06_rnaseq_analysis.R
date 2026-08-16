# Col-0 non-stressed versus 24-h 37C heat-stress RNA-seq analysis.
# Covers Fig. 5 and final Tables S9-S10.
#
# Gene-level differential expression: DESeq2.
# Global TE differential expression: fractional multimapping-aware counts,
# edgeR::filterByExpr expression filtering, then limma-voom.
# Candidate-window signal: fractional counts normalized to assigned gene reads.
# Fig. 5B compares the SUM of CPM across the sixteen curated ONSEN terminal
# windows within each biological replicate (three replicates per condition).

source("ONSEN_functions.R")
require_packages(c("data.table", "dplyr", "tidyr", "readr", "tibble", "ggplot2"))
message_config()

GENE_CALL_PADJ <- 0.05
LOG2FC_CUTOFF <- 1
TE_DESCRIPTIVE_PSEUDOCOUNT <- 0.5
CANDIDATE_PSEUDOCOUNT_CPM <- 0.05

read_final_sheet <- function(filename, skip = 2L) {
  data.table::fread(
    repo_file(file.path("supplementary_table_source", filename)),
    skip = skip, data.table = FALSE, check.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# Lightweight repository validation using deposited final-numbered source data.
# Full re-analysis from processed count matrices is enabled with
# ONSEN_RUN_LARGE_STEPS=true or ONSEN_FORCE_RESCAN=true.
# -----------------------------------------------------------------------------
if (!ONSEN_RUN_LARGE_STEPS && !ONSEN_FORCE_RESCAN) {
  s9_summary <- read_final_sheet("Table_S9__Summary.tsv")
  s10_stats <- read_final_sheet("Table_S10__S10D_Statistics.tsv", skip = 11L)

  gene_row <- s9_summary[s9_summary[[1]] == "Genes", , drop = FALSE]
  te_row <- s9_summary[s9_summary[[1]] == "Transposable elements", , drop = FALSE]
  stopifnot(
    nrow(gene_row) == 1L,
    nrow(te_row) == 1L,
    as.integer(gene_row[[2]]) == 32833L,
    as.integer(gene_row[[3]]) == 25912L,
    as.integer(gene_row[[4]]) == 23922L,
    as.integer(gene_row[[5]]) == 3912L,
    as.integer(gene_row[[6]]) == 4185L,
    as.integer(te_row[[2]]) == 31189L,
    as.integer(te_row[[3]]) == 7632L,
    as.integer(te_row[[4]]) == 2101L,
    as.integer(te_row[[5]]) == 692L,
    as.integer(te_row[[6]]) == 397L
  )

  # Parse the deposited statistics sheet explicitly because it also contains
  # replicate-sum rows above the statistics block.
  s10_text <- readLines(
    repo_file(file.path("supplementary_table_source", "Table_S10__S10D_Statistics.tsv")),
    warn = FALSE
  )
  stopifnot(
    any(grepl("0.000757270061088614", s10_text, fixed = TRUE)),
    any(grepl("0.08085559837005224", s10_text, fixed = TRUE)),
    any(grepl("Replicate-summed fractional CPM", s10_text, fixed = TRUE))
  )

  message(
    "Validated deposited final RNA-seq source data: 32,833 genes; 31,189 TEs; ",
    "2,101 tested TEs; Fig. 5B Welch P=7.57270061088614e-4."
  )
} else {
  require_packages(c("DESeq2", "edgeR", "limma", "rtracklayer"))

  # ---------------------------------------------------------------------------
  # Metadata and exact DDBJ provenance
  # ---------------------------------------------------------------------------
  metadata_file <- find_any_input(c(
    "Col0_NS_vs_24h_37C_HS_sample_metadata.csv",
    "RNAseq_sample_metadata_template.csv"
  ))
  sample_meta <- read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)

  rename_first <- function(x, target, candidates) {
    hit <- candidates[candidates %in% names(x)][1]
    if (!is.na(hit) && hit != target) names(x)[names(x) == hit] <- target
    x
  }
  sample_meta <- rename_first(sample_meta, "sample_id", c("sample_id", "sample", "Sample"))
  sample_meta <- rename_first(sample_meta, "treatment", c("treatment", "condition", "group"))
  assert_columns(
    sample_meta,
    c("sample_id", "treatment", "replicate", "ecotype", "biosample_accession",
      "bioproject", "laboratory_provenance"),
    "RNA-seq metadata"
  )
  sample_meta$treatment <- ifelse(
    grepl("HS|heat|37", sample_meta$treatment, ignore.case = TRUE), "HS", "NS"
  )

  expected <- data.frame(
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
    bioproject = c(rep("PRJDB39904", 3L), rep("PRJDB42759", 3L)),
    stringsAsFactors = FALSE
  )
  idx <- match(expected$sample_id, sample_meta$sample_id)
  if (nrow(sample_meta) != 6L || anyNA(idx)) stop("Expected exactly the six validated Col-0 libraries.")
  for (field in c("treatment", "replicate", "biosample_accession", "bioproject")) {
    if (!identical(as.character(sample_meta[[field]][idx]), as.character(expected[[field]]))) {
      stop("RNA-seq metadata mismatch in field: ", field)
    }
  }
  if (!all(tolower(sample_meta$ecotype[idx]) %in% c("col-0", "col0")) ||
      !all(sample_meta$laboratory_provenance[idx] == "Generated in our laboratory")) {
    stop("RNA-seq metadata contains incorrect ecotype or laboratory provenance.")
  }
  sample_meta <- sample_meta[idx, , drop = FALSE]
  sample_meta$treatment <- factor(sample_meta$treatment, levels = c("NS", "HS"))

  sample_ids <- sample_meta$sample_id
  ns_ids <- sample_meta$sample_id[sample_meta$treatment == "NS"]
  hs_ids <- sample_meta$sample_id[sample_meta$treatment == "HS"]

  # ---------------------------------------------------------------------------
  # Generic count-table reader. Annotation columns are retained separately.
  # ---------------------------------------------------------------------------
  read_count_table <- function(path, id_candidates) {
    x <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
    # featureCounts may retain BAM paths as column names; map their basenames.
    for (j in seq_along(names(x))) {
      base <- tools::file_path_sans_ext(basename(names(x)[j]))
      if (base %in% sample_ids) names(x)[j] <- base
    }
    id_col <- id_candidates[id_candidates %in% names(x)][1]
    if (is.na(id_col)) id_col <- names(x)[1]
    missing_samples <- setdiff(sample_ids, names(x))
    if (length(missing_samples)) stop("Missing count columns: ", paste(missing_samples, collapse = ", "))
    counts <- as.matrix(x[, sample_ids, drop = FALSE])
    storage.mode(counts) <- "numeric"
    rownames(counts) <- as.character(x[[id_col]])
    list(raw = x, id_col = id_col, counts = counts)
  }

  # ---------------------------------------------------------------------------
  # Gene-level DESeq2: model the 25,912 genes with non-zero signal, then rejoin
  # results to the complete 32,833-gene annotation universe.
  # ---------------------------------------------------------------------------
  gene_file <- find_any_input(c(
    "Col0_NS_vs_24h_37C_HS_gene_counts_raw.tsv",
    "gene_counts.tsv",
    "Col0_NS_vs_24h_37C_HS_gene_counts_annotation.tsv"
  ))
  gene <- read_count_table(gene_file, c("gene_id", "Geneid", "GeneID"))
  gene_ids <- rownames(gene$counts)
  if (length(gene_ids) != 32833L) {
    stop("Expected 32,833 annotated genes in the complete-TAIR10 count universe; observed ", length(gene_ids), ".")
  }
  gene_nonzero <- rowSums(gene$counts) > 0
  if (sum(gene_nonzero) != 25912L) {
    stop("Expected 25,912 genes with non-zero signal; observed ", sum(gene_nonzero), ".")
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(gene$counts[gene_nonzero, , drop = FALSE]),
    colData = data.frame(treatment = sample_meta$treatment, row.names = sample_ids),
    design = ~ treatment
  )
  dds <- DESeq2::DESeq(dds)
  gene_fit <- as.data.frame(DESeq2::results(dds, contrast = c("treatment", "HS", "NS"))) |>
    tibble::rownames_to_column("gene_id")

  gene_all <- data.frame(gene_id = gene_ids, stringsAsFactors = FALSE) |>
    dplyr::left_join(gene_fit, by = "gene_id") |>
    dplyr::mutate(
      padj_available = as.integer(!is.na(padj)),
      DE_call = dplyr::case_when(
        !is.na(padj) & padj < GENE_CALL_PADJ & log2FoldChange >= LOG2FC_CUTOFF ~ "Up",
        !is.na(padj) & padj < GENE_CALL_PADJ & log2FoldChange <= -LOG2FC_CUTOFF ~ "Down",
        TRUE ~ "NS"
      )
    )
  if (sum(!is.na(gene_all$pvalue)) != 25911L || sum(!is.na(gene_all$padj)) != 23922L ||
      sum(gene_all$DE_call == "Up") != 3912L || sum(gene_all$DE_call == "Down") != 4185L) {
    stop("Gene-level DE results do not reproduce final Table S9 counts; check count matrix/DESeq2 inputs.")
  }
  safe_write_csv(gene_all, "Table_S9A_genes_repository.csv")

  # ---------------------------------------------------------------------------
  # Global TE analysis: all 31,189 Araport11 transposable_element loci retained;
  # expression filtering selects the 2,101 statistically tested loci.
  # ---------------------------------------------------------------------------
  te_file <- find_any_input(c(
    "Col0_NS_vs_24h_37C_HS_TE_fractional_multimap_counts.tsv",
    "Col0_NS_vs_24h_37C_HS_TE_fractional_multimap_counts_annotated.tsv",
    "TE_fractional_multimap_counts.tsv"
  ))
  te <- read_count_table(te_file, c("TE_ID", "te_id", "Geneid", "feature_id", "transposable_element_id"))
  te_ids <- rownames(te$counts)
  if (length(te_ids) != 31189L) stop("Expected 31,189 Araport11 TE loci; observed ", length(te_ids), ".")

  te_nonzero <- rowSums(te$counts) > 0
  if (sum(te_nonzero) != 7632L) stop("Expected 7,632 TEs with non-zero signal; observed ", sum(te_nonzero), ".")

  dge <- edgeR::DGEList(counts = te$counts, group = sample_meta$treatment)
  keep_te <- edgeR::filterByExpr(dge, group = sample_meta$treatment)
  if (sum(keep_te) != 2101L) {
    stop("edgeR::filterByExpr did not reproduce the final 2,101 tested TEs; observed ", sum(keep_te), ".")
  }
  dge_f <- dge[keep_te, , keep.lib.sizes = FALSE]
  dge_f <- edgeR::calcNormFactors(dge_f)
  design <- stats::model.matrix(~0 + sample_meta$treatment)
  colnames(design) <- c("NS", "HS")
  v <- limma::voom(dge_f, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::contrasts.fit(fit, limma::makeContrasts(HS_minus_NS = HS - NS, levels = design))
  fit <- limma::eBayes(fit)
  te_fit <- limma::topTable(fit, coef = "HS_minus_NS", number = Inf, sort.by = "none") |>
    tibble::rownames_to_column("TE_ID")

  te_means <- data.frame(
    TE_ID = te_ids,
    mean_NS_signal = rowMeans(te$counts[, ns_ids, drop = FALSE]),
    mean_HS_signal = rowMeans(te$counts[, hs_ids, drop = FALSE]),
    mean_signal_all_samples = rowMeans(te$counts),
    statistically_tested = as.integer(te_ids %in% rownames(dge_f$counts)),
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      descriptive_log2FC = log2(
        (mean_HS_signal + TE_DESCRIPTIVE_PSEUDOCOUNT) /
          (mean_NS_signal + TE_DESCRIPTIVE_PSEUDOCOUNT)
      )
    )

  te_stats <- te_means |>
    dplyr::left_join(
      te_fit |>
        dplyr::transmute(
          TE_ID,
          model_log2FC = logFC,
          log2FC_SE = ifelse(is.finite(t) & t != 0, logFC / t, NA_real_),
          test_statistic = t,
          P_value = P.Value,
          BH_adjusted_P = adj.P.Val
        ),
      by = "TE_ID"
    ) |>
    dplyr::mutate(
      DE_call = dplyr::case_when(
        !is.na(BH_adjusted_P) & BH_adjusted_P < GENE_CALL_PADJ & model_log2FC >= LOG2FC_CUTOFF ~ "Up",
        !is.na(BH_adjusted_P) & BH_adjusted_P < GENE_CALL_PADJ & model_log2FC <= -LOG2FC_CUTOFF ~ "Down",
        TRUE ~ "NS"
      )
    )
  if (sum(te_stats$DE_call == "Up") != 692L || sum(te_stats$DE_call == "Down") != 397L ||
      sum(!is.na(te_stats$BH_adjusted_P)) != 2101L) {
    stop("TE limma-voom results do not reproduce final Table S9 counts; check count/filtering inputs.")
  }

  # Add Araport11 metadata. Prefer annotation columns carried with the processed
  # count table; otherwise import the public GFF3 listed in INPUT_PROVENANCE.tsv.
  te_meta_candidates <- setdiff(names(te$raw), sample_ids)
  te_meta <- te$raw[, te_meta_candidates, drop = FALSE]
  names(te_meta)[names(te_meta) == te$id_col] <- "TE_ID"
  if (!all(c("TE_alias", "chromosome", "start", "end", "strand") %in% names(te_meta))) {
    gff_file <- find_input("Araport11_GFF3_genes_transposons.20241001.gff.gz")
    gff <- rtracklayer::import(gff_file)
    gdf <- as.data.frame(gff)
    type_col <- c("type", "feature", "feature_type")[c("type", "feature", "feature_type") %in% names(gdf)][1]
    if (!is.na(type_col)) gdf <- gdf[gdf[[type_col]] == "transposable_element", , drop = FALSE]
    id_col <- c("ID", "gene_id", "Name")[c("ID", "gene_id", "Name") %in% names(gdf)][1]
    if (is.na(id_col)) stop("Could not identify TE ID in Araport11 GFF3.")
    alias_col <- c("Alias", "alias", "family")[c("Alias", "alias", "family") %in% names(gdf)][1]
    name_col <- c("Name", "ID")[c("Name", "ID") %in% names(gdf)][1]
    te_meta <- data.frame(
      TE_ID = as.character(gdf[[id_col]]),
      TE_alias = if (is.na(alias_col)) NA_character_ else as.character(gdf[[alias_col]]),
      TE_name = if (is.na(name_col)) as.character(gdf[[id_col]]) else as.character(gdf[[name_col]]),
      chromosome = as.character(gdf$seqnames),
      start = as.integer(gdf$start), end = as.integer(gdf$end),
      strand = as.character(gdf$strand),
      stringsAsFactors = FALSE
    )
  }

  canonical_onsen <- data.frame(
    canonical_ONSEN_copy = paste0("ONSEN", 1:8),
    TE_ID = c("AT1TE12295", "AT3TE92525", "AT5TE15240", "AT1TE71045",
              "AT1TE59755", "AT3TE89830", "AT1TE24850", "AT3TE54550"),
    stringsAsFactors = FALSE
  )

  te_final <- te_stats |>
    dplyr::left_join(te_meta, by = "TE_ID") |>
    dplyr::left_join(canonical_onsen, by = "TE_ID")
  safe_write_csv(te_final, "Table_S9B_TEs_repository.csv")
  onsen_te <- canonical_onsen |>
    dplyr::left_join(te_final, by = c("TE_ID", "canonical_ONSEN_copy"))
  if (nrow(onsen_te) != 8L || any(onsen_te$DE_call != "Up")) {
    stop("Canonical ONSEN extraction did not reproduce eight upregulated loci.")
  }
  safe_write_csv(onsen_te, "Table_S9C_ONSEN_repository.csv")

  # ---------------------------------------------------------------------------
  # Candidate-window analysis for Fig. 5B-C and Table S10.
  # ---------------------------------------------------------------------------
  candidate_file <- find_any_input(c(
    "Col0_NS_vs_24h_37C_HS_candidate_ONSEN_HSF_rich_TE_loci_fractional_multimap_counts_annotated.tsv",
    "Col0_NS_vs_24h_37C_HS_candidate_ONSEN_HSF_rich_TE_loci_fractional_multimap_counts.tsv",
    "candidate_window_fractional_counts_annotated.tsv"
  ))
  candidate <- read_count_table(candidate_file, c("candidate_id", "Candidate window", "window_id"))

  library_file <- find_any_input(c(
    "Col0_NS_vs_24h_37C_HS_gene_assigned_count_library_summary.tsv",
    "gene_assigned_count_library_summary.tsv"
  ), required = FALSE)
  if (!is.na(library_file)) {
    lib <- data.table::fread(library_file, data.table = FALSE, check.names = FALSE)
    sample_col <- c("sample_id", "sample", "Sample")[c("sample_id", "sample", "Sample") %in% names(lib)][1]
    assigned_col <- names(lib)[grepl("assigned.*gene|gene.*assigned|assigned_counts|assigned", names(lib), ignore.case = TRUE)][1]
    if (is.na(sample_col) || is.na(assigned_col)) stop("Could not identify assigned-gene library-size columns.")
    library_sizes <- setNames(as.numeric(lib[[assigned_col]]), as.character(lib[[sample_col]]))[sample_ids]
  } else {
    library_sizes <- colSums(gene$counts)[sample_ids]
  }
  if (any(!is.finite(library_sizes) | library_sizes <= 0)) stop("Invalid gene-assigned library sizes.")

  candidate_ids <- rownames(candidate$counts)
  candidate_cpm <- sweep(candidate$counts, 2, library_sizes, "/") * 1e6
  windows <- read.csv(repo_file("ONSEN_Col0_terminal_candidate_windows.csv"), stringsAsFactors = FALSE)
  curated_ids <- windows$window_id
  if (!all(curated_ids %in% candidate_ids)) {
    stop("Candidate count matrix does not contain all sixteen curated ONSEN terminal windows.")
  }

  candidate_meta <- candidate$raw[, setdiff(names(candidate$raw), sample_ids), drop = FALSE]
  names(candidate_meta)[names(candidate_meta) == candidate$id_col] <- "candidate_id"
  candidate_meta$candidate_id <- as.character(candidate_meta$candidate_id)
  if (!"candidate_class" %in% names(candidate_meta)) {
    candidate_meta$candidate_class <- ifelse(
      candidate_meta$candidate_id %in% curated_ids,
      "ONSEN terminal candidate windows",
      ifelse(grepl("ONSEN|ATCOPIA78", candidate_meta$candidate_id, ignore.case = TRUE),
             "Additional ATCOPIA78 candidate regions", "HSF-rich non-ONSEN TE")
    )
  }

  cpm_table <- data.frame(candidate_id = candidate_ids, candidate_cpm, check.names = FALSE) |>
    dplyr::left_join(candidate_meta, by = "candidate_id")
  safe_write_csv(cpm_table, "Table_S10C_candidate_replicate_CPM_repository.csv")

  candidate_long <- data.frame(candidate_id = candidate_ids, candidate_cpm, check.names = FALSE) |>
    tidyr::pivot_longer(cols = dplyr::all_of(sample_ids), names_to = "sample_id", values_to = "fractional_CPM") |>
    dplyr::left_join(sample_meta |> dplyr::select(sample_id, treatment, replicate), by = "sample_id") |>
    dplyr::left_join(candidate_meta |> dplyr::select(dplyr::any_of(c("candidate_id", "candidate_class"))), by = "candidate_id")

  window_summary <- candidate_long |>
    dplyr::group_by(candidate_id, candidate_class, treatment) |>
    dplyr::summarise(mean_CPM = mean(fractional_CPM), .groups = "drop") |>
    tidyr::pivot_wider(names_from = treatment, values_from = mean_CPM) |>
    dplyr::mutate(
      descriptive_log2FC = log2((HS + CANDIDATE_PSEUDOCOUNT_CPM) /
                                  (NS + CANDIDATE_PSEUDOCOUNT_CPM)),
      Fig5C_x = log2(HS + 1)
    )
  safe_write_csv(window_summary, "Table_S10B_individual_windows_repository.csv")

  onsen_replicate_sum <- candidate_long |>
    dplyr::filter(candidate_id %in% curated_ids) |>
    dplyr::group_by(treatment, sample_id, replicate) |>
    dplyr::summarise(
      summed_fractional_CPM = sum(fractional_CPM),
      mean_fractional_CPM_per_window = mean(fractional_CPM),
      .groups = "drop"
    )
  if (nrow(onsen_replicate_sum) != 6L) stop("Expected six replicate-summed ONSEN observations.")
  welch <- stats::t.test(summed_fractional_CPM ~ treatment, data = onsen_replicate_sum, var.equal = FALSE)
  wilcox <- stats::wilcox.test(summed_fractional_CPM ~ treatment, data = onsen_replicate_sum, exact = FALSE)
  if (abs(welch$p.value - 0.000757270061088614) > 1e-10 ||
      abs(wilcox$p.value - 0.08085559837005224) > 1e-10) {
    stop("Candidate-window tests do not reproduce final Table S10D; check count/library-size inputs.")
  }
  candidate_stats <- data.frame(
    test = c("Welch's t-test", "Wilcoxon rank-sum test"),
    comparison_unit = "Replicate-summed fractional CPM across 16 curated ONSEN terminal windows",
    n_NS = 3L, n_HS = 3L,
    statistic = c(unname(welch$statistic), unname(wilcox$statistic)),
    degrees_of_freedom = c(unname(welch$parameter), NA_real_),
    two_sided_P = c(welch$p.value, wilcox$p.value),
    stringsAsFactors = FALSE
  )
  safe_write_csv(onsen_replicate_sum, "Table_S10D_replicate_sums_repository.csv")
  safe_write_csv(candidate_stats, "Table_S10D_statistics_repository.csv")

  # ---------------------------------------------------------------------------
  # Reproducible figure data/plots. These do not replace the journal PDFs; they
  # reproduce the analysis quantities used in Fig. 5.
  # ---------------------------------------------------------------------------
  if (ONSEN_MAKE_FIGURES) {
    p5b <- ggplot2::ggplot(
      onsen_replicate_sum,
      ggplot2::aes(x = treatment, y = summed_fractional_CPM, fill = treatment)
    ) +
      ggplot2::geom_col(stat = "summary", fun = mean, width = 0.65, colour = "black") +
      ggplot2::geom_point(size = 2.5, position = ggplot2::position_jitter(width = 0.05)) +
      ggplot2::labs(x = NULL, y = "Candidate-window signal") +
      theme_onsen(13) + ggplot2::theme(legend.position = "none")
    save_plot_pair(p5b, "Fig5B_ONSEN_16window_replicate_sum", 5.2, 4.8)

    p5c <- ggplot2::ggplot(
      window_summary,
      ggplot2::aes(Fig5C_x, candidate_class, size = pmax(0, descriptive_log2FC), colour = candidate_class)
    ) +
      ggplot2::geom_point(alpha = 0.8) +
      ggplot2::labs(x = "HS candidate-window signal [log2(CPM + 1)]", y = NULL,
                    size = "HS/NS\nlog2FC", colour = "Candidate class") +
      theme_onsen(11)
    save_plot_pair(p5c, "Fig5C_candidate_window_signal", 8.0, 5.8)

    gene_plot <- gene_all |>
      dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) |>
      dplyr::mutate(neglog10_padj = -log10(pmax(padj, .Machine$double.xmin)))
    p5d <- ggplot2::ggplot(gene_plot, ggplot2::aes(log2FoldChange, neglog10_padj, colour = DE_call)) +
      ggplot2::geom_point(alpha = 0.55, size = 0.7) +
      ggplot2::labs(x = "log2 fold change (HS / NS)", y = "-log10 adjusted P") + theme_onsen(11)
    save_plot_pair(p5d, "Fig5D_gene_volcano", 6.4, 5.2)

    te_plot <- te_final |>
      dplyr::filter(!is.na(BH_adjusted_P), !is.na(model_log2FC)) |>
      dplyr::mutate(neglog10_padj = -log10(pmax(BH_adjusted_P, .Machine$double.xmin)))
    p5e <- ggplot2::ggplot(te_plot, ggplot2::aes(model_log2FC, neglog10_padj, colour = DE_call)) +
      ggplot2::geom_point(alpha = 0.6, size = 0.9) +
      ggplot2::labs(x = "log2 fold change (HS / NS)", y = "-log10 adjusted P") + theme_onsen(11)
    save_plot_pair(p5e, "Fig5E_TE_volcano", 6.4, 5.2)
  }

  message("Full gene, TE and candidate-window RNA-seq analysis completed successfully.")
}
