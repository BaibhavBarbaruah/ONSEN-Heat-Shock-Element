# Validate the public repository against the final Biology Open manuscript data
# and supplementary package (Fig. S1-Fig. S4; Table S1-Table S13).

source("ONSEN_config.R")

for (package in c("data.table", "dplyr", "openxlsx")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package '", package, "' is required for validation.", call. = FALSE)
  }
}
message_config()

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

read_text <- function(path) {
  paste(readLines(file.path(REPO_ROOT, path), warn = FALSE), collapse = "\n")
}

read_final_tsv <- function(path, skip = 2L) {
  data.table::fread(
    file.path(REPO_ROOT, "supplementary_table_source", path),
    skip = skip, data.table = FALSE, check.names = FALSE
  )
}

required_scripts <- c(
  "ONSEN_config.R", "ONSEN_functions.R", "00_install_packages.R",
  "00_run_pipeline.R", "01_native_mutated_motif_analysis.R",
  "02_constrained_mutant_sensitivity.R", "03_col0_HSF_and_TE_background.R",
  "03B_threshold_and_continuous_sensitivity.R",
  "03C_direct_LTR_HSF_comparison.R", "04_nonredundant_HSF_locations.R",
  "05A_direct_LTR_methylation_analysis.R", "05_methylation_analysis.R",
  "06_rnaseq_analysis.R", "07_accession_analysis.R",
  "07B_figure7_identity_and_logo_workflow.R",
  "08_write_supplementary_tables.R", "09_write_session_info.R",
  "10_validate_manuscript_outputs.R"
)
required_metadata <- c(
  "README.md", "LICENSE", "CITATION.cff", "DATA_AVAILABILITY.md",
  "REPRODUCIBILITY_NOTES.md", "INPUT_PROVENANCE.tsv",
  "REPRODUCIBILITY_MATRIX.tsv", "R_PACKAGE_REQUIREMENTS.txt",
  "FULL_PACKAGE_MANIFEST.txt", "RNAseq_sample_metadata_template.csv",
  "ONSEN_49bp_sequences.fasta", "ONSEN_HSE_units_and_substitutions.csv",
  "ONSEN_Col0_terminal_candidate_windows.csv",
  "Arabidopsis_HSF_models_JASPAR2026.csv",
  "source_data/Figure4A_direct_LTR_methylation_summary.tsv",
  "source_data/Figure4C_direct_LTR_HSF_summary.tsv",
  "source_data/Figure7_variant_metrics.tsv",
  "source_data/Figure7_logo_model_metadata.csv",
  "motifs/MA1667.2_HSFC1.jaspar", "motifs/MA0981.2_DOF1.8.jaspar",
  "supplementary_table_source/README.md",
  "supplementary_table_source/SOURCE_SHEET_MANIFEST.tsv"
)
required_tables <- sprintf("Table_S%d.xlsx", 1:13)
required_files <- c(required_scripts, required_metadata, required_tables)
missing_files <- required_files[!file.exists(file.path(REPO_ROOT, required_files))]
add_check(
  "Required public files are present",
  !length(missing_files),
  if (length(missing_files)) paste(missing_files, collapse = "; ") else "complete"
)

obsolete_public_files <- c(
  "03B_threshold_and_continuous_sensitivity_part1.R",
  "03B_threshold_and_continuous_sensitivity_part2.R",
  "DDBJ_Update_Change_Log.md", "FINAL_NUMBERING_MAP.tsv",
  "CHANGELOG.md", "FILE_CHECKSUMS_SHA256.tsv"
)
present_obsolete <- obsolete_public_files[
  file.exists(file.path(REPO_ROOT, obsolete_public_files))
]
add_check(
  "Obsolete/internal public files are absent",
  !length(present_obsolete),
  if (length(present_obsolete)) paste(present_obsolete, collapse = "; ") else "clean"
)

# Public-facing R files must not contain hard-coded Windows drive roots.
r_files <- list.files(REPO_ROOT, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
drive_specific <- character()
for (path in r_files) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (grepl("(^|[^A-Za-z0-9])([A-Za-z]):[/\\\\]", text, perl = TRUE)) {
    drive_specific <- c(drive_specific, basename(path))
  }
}
add_check(
  "No hard-coded Windows drive roots",
  !length(drive_specific),
  if (length(drive_specific)) paste(drive_specific, collapse = "; ") else "portable"
)

# Core documentation must state the same final analytical definitions.
readme <- read_text("README.md")
data_availability <- read_text("DATA_AVAILABILITY.md")
notes <- read_text("REPRODUCIBILITY_NOTES.md")
provenance <- read_text("INPUT_PROVENANCE.tsv")
packages <- read_text("R_PACKAGE_REQUIREMENTS.txt")

add_check(
  "README has final supplementary scope",
  grepl("Fig. S1-Fig. S4", readme, fixed = TRUE) &&
    grepl("Table S1-Table S13", readme, fixed = TRUE),
  "Figs S1-S4; Tables S1-S13"
)
add_check(
  "README defines exact-coordinate HSF metric",
  grepl("unique `(start, end)` coordinate placements", readme, fixed = TRUE) &&
    grepl("identical start and end coordinates are counted once", readme, fixed = TRUE),
  "principal metric documented"
)
add_check(
  "README has final HSF background counts",
  grepl("1,930", readme, fixed = TRUE) &&
    grepl("49.375", readme, fixed = TRUE) &&
    grepl("6.81238953796413", readme, fixed = TRUE),
  "16 ONSEN windows versus 1,930 strict non-ONSEN TE regions"
)
add_check(
  "README has final RNA-seq methods",
  grepl("32,833", readme, fixed = TRUE) &&
    grepl("31,189", readme, fixed = TRUE) &&
    grepl("2,101", readme, fixed = TRUE) &&
    grepl("limma-voom", readme, fixed = TRUE) &&
    grepl("summed across the sixteen curated ONSEN terminal windows", readme, fixed = TRUE),
  "gene DESeq2; TE limma-voom; replicate-summed candidate signal"
)
add_check(
  "DDBJ accessions are synchronized",
  all(vapply(c(
    "PRJDB39904", "PRJDB42759",
    "SAMD01789795", "SAMD01789796", "SAMD01789797",
    "SAMD01943917", "SAMD01943918", "SAMD01943919"
  ), function(x) grepl(x, paste(readme, data_availability, provenance), fixed = TRUE), logical(1))),
  "both BioProjects and all six BioSamples present"
)
add_check(
  "Reference provenance states Ensembl Plants release 63",
  grepl("Ensembl Plants release 63", provenance, fixed = TRUE) &&
    grepl("Arabidopsis_thaliana.TAIR10.63.gtf", provenance, fixed = TRUE),
  "release 63 gene annotation"
)
add_check(
  "TE-analysis dependencies are documented",
  grepl("edgeR", packages, fixed = TRUE) && grepl("limma", packages, fixed = TRUE),
  "edgeR and limma present"
)

# Source-sheet manifest and formatted workbooks.
source_manifest <- read.delim(
  file.path(REPO_ROOT, "supplementary_table_source", "SOURCE_SHEET_MANIFEST.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
add_check(
  "Final TSV source manifest has 36 worksheets",
  nrow(source_manifest) == 36L,
  paste("observed", nrow(source_manifest))
)
add_check(
  "Source manifest covers Tables S1-S13",
  all(sprintf("Table_S%d.xlsx", 1:13) %in% source_manifest$workbook),
  "all thirteen workbooks represented"
)

unreadable_tables <- character()
for (table_name in required_tables) {
  sheets <- tryCatch(
    openxlsx::getSheetNames(file.path(REPO_ROOT, table_name)),
    error = function(e) character()
  )
  if (!length(sheets)) unreadable_tables <- c(unreadable_tables, table_name)
}
add_check(
  "All thirteen final workbooks open",
  !length(unreadable_tables),
  if (length(unreadable_tables)) paste(unreadable_tables, collapse = "; ") else "13/13"
)

s10_path <- file.path(REPO_ROOT, "Table_S10.xlsx")
s10_sheets <- openxlsx::getSheetNames(s10_path)
add_check(
  "Table S10 contains the statistics sheet",
  identical(
    s10_sheets,
    c("S10A_Class_summary", "S10B_Individual_windows",
      "S10C_Replicate_CPM", "S10D_Statistics")
  ),
  paste(s10_sheets, collapse = "; ")
)
s10_text <- readLines(
  file.path(REPO_ROOT, "supplementary_table_source", "Table_S10__S10D_Statistics.tsv"),
  warn = FALSE
)
add_check(
  "Table S10 statistics are exact",
  any(grepl("0.000757270061088614", s10_text, fixed = TRUE)) &&
    any(grepl("0.08085559837005224", s10_text, fixed = TRUE)) &&
    any(grepl("Replicate-summed fractional CPM", s10_text, fixed = TRUE)),
  "Welch and Wilcoxon results on replicate sums"
)

# Table S7: exact-coordinate threshold analysis.
s7_stats <- read_final_tsv("Table_S7__S7A_threshold_stats.tsv")
s7_regions <- read_final_tsv("Table_S7__S7D_region_threshold.tsv")
s7_key <- s7_stats[s7_stats[["Relative PWM-score threshold"]] == 0.85, , drop = FALSE]
add_check(
  "Table S7 has final 0.85 comparison",
  nrow(s7_key) == 1L &&
    as.integer(s7_key[["ONSEN n"]]) == 16L &&
    as.integer(s7_key[["Background n"]]) == 1930L &&
    abs(s7_key[["ONSEN median density (placements/kb)"]] - 49.375) < 1e-12 &&
    abs(s7_key[["Background median density (placements/kb)"]] - 6.81238953796413) < 1e-12,
  "n=16/1930; medians=49.375/6.81238953796413"
)
add_check(
  "Table S7 region-level metric is exact-coordinate",
  all(c(
    "Non-redundant HSF motif-coordinate placements",
    "Density (placements/kb)"
  ) %in% names(s7_regions)) &&
    sum(s7_regions[["Region class"]] == "ONSEN terminal candidate windows" &
          s7_regions[["Relative PWM-score threshold"]] == 0.85) == 16L &&
    sum(s7_regions[["Region class"]] == "Strict non-ONSEN TE background" &
          s7_regions[["Relative PWM-score threshold"]] == 0.85) == 1930L,
  "exact-coordinate regional source data complete"
)

# Table S9: genome-wide differential expression.
s9_summary <- read_final_tsv("Table_S9__Summary.tsv")
gene_row <- s9_summary[s9_summary[[1]] == "Genes", , drop = FALSE]
te_row <- s9_summary[s9_summary[[1]] == "Transposable elements", , drop = FALSE]
add_check(
  "Table S9 gene totals are exact",
  nrow(gene_row) == 1L &&
    identical(as.integer(unlist(gene_row[1, 2:6], use.names = FALSE)),
              c(32833L, 25912L, 23922L, 3912L, 4185L)),
  "32833 total; 25912 modeled; 23922 padj; 3912 up; 4185 down"
)
add_check(
  "Table S9 TE totals are exact",
  nrow(te_row) == 1L &&
    identical(as.integer(unlist(te_row[1, 2:6], use.names = FALSE)),
              c(31189L, 7632L, 2101L, 692L, 397L)),
  "31189 total; 7632 nonzero; 2101 tested; 692 up; 397 down"
)

s9_onsen <- read_final_tsv("Table_S9__S9C_ONSEN.tsv")
call_col <- names(s9_onsen)[grepl("call$|DE_call|classification", names(s9_onsen), ignore.case = TRUE)][1]
add_check(
  "Table S9 contains eight upregulated canonical ONSEN loci",
  nrow(s9_onsen) == 8L && !is.na(call_col) && all(s9_onsen[[call_col]] == "Up"),
  "8/8 canonical ONSEN loci classified Up"
)

# RNA-seq metadata: exact library/accession map.
rna_meta <- read.csv(
  file.path(REPO_ROOT, "RNAseq_sample_metadata_template.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
expected_biosamples <- c(
  "SAMD01789795", "SAMD01789796", "SAMD01789797",
  "SAMD01943917", "SAMD01943918", "SAMD01943919"
)
expected_bioprojects <- c(rep("PRJDB39904", 3L), rep("PRJDB42759", 3L))
add_check(
  "RNA-seq metadata map is exact",
  nrow(rna_meta) == 6L &&
    identical(rna_meta$biosample_accession, expected_biosamples) &&
    identical(rna_meta$bioproject, expected_bioprojects) &&
    identical(as.integer(rna_meta$replicate), c(1L, 2L, 3L, 1L, 2L, 3L)) &&
    all(rna_meta$laboratory_provenance == "Generated in our laboratory"),
  "six study-generated libraries; three biological replicates/condition"
)

# Code-level guards against the two previously identified analytical mismatches.
s3b <- read_text("03B_threshold_and_continuous_sensitivity.R")
s6 <- read_text("06_rnaseq_analysis.R")
add_check(
  "Threshold script uses coordinate deduplication",
  grepl("unique(coordinate_sets[[t]])", s3b, fixed = TRUE) &&
    grepl("forward_start", s3b, fixed = TRUE) &&
    grepl("forward_end", s3b, fixed = TRUE),
  "deduplicates physical start/end coordinates across model/strand predictions"
)
add_check(
  "RNA-seq script implements final TE method",
  grepl("edgeR::filterByExpr", s6, fixed = TRUE) &&
    grepl("limma::voom", s6, fixed = TRUE) &&
    grepl("2101L", s6, fixed = TRUE),
  "fractional TE counts; expression filtering; limma-voom"
)
add_check(
  "RNA-seq script implements replicate-summed Fig. 5B",
  grepl("summed_fractional_CPM = sum(fractional_CPM)", s6, fixed = TRUE) &&
    grepl("stats::t.test(summed_fractional_CPM ~ treatment", s6, fixed = TRUE),
  "sum across 16 ONSEN windows within each biological replicate"
)

# Figure 7 logo identifiers and author spelling.
logo <- read.csv(
  file.path(REPO_ROOT, "source_data", "Figure7_logo_model_metadata.csv"),
  stringsAsFactors = FALSE
)
add_check(
  "Figure 7 logo models are exact",
  any(logo$model_name == "HSFC1" & logo$JASPAR_ID == "MA1667.2") &&
    any(logo$model_name == "DOF1.8" & logo$JASPAR_ID == "MA0981.2"),
  "HSFC1 MA1667.2; DOF1.8 MA0981.2"
)

metadata_text <- paste(
  readme, data_availability, read_text("CITATION.cff"),
  read_text("REPRODUCIBILITY_MATRIX.tsv"), notes
)
add_check(
  "Author surname is consistently Airlangga",
  grepl("Rahmadani P. Airlangga", metadata_text, fixed = TRUE) &&
    !grepl("Airalangga", metadata_text, fixed = TRUE),
  "Rahmadani P. Airlangga"
)

report <- dplyr::bind_rows(checks)
dir.create(ONSEN_OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
write.csv(
  report,
  file.path(ONSEN_OUTPUT_ROOT, "repository_validation_report.csv"),
  row.names = FALSE
)

if (any(!report$passed)) {
  print(report[!report$passed, ], row.names = FALSE)
  stop("Repository validation failed; see repository_validation_report.csv.", call. = FALSE)
}

message("\n============================================================")
message("FINAL PUBLIC REPOSITORY VALIDATION PASSED")
message("============================================================")
