# Validate the final flat repository and journal-numbered source package.

source("ONSEN_config.R")

if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required for validation.", call. = FALSE)
message_config()

required_scripts <- c(
  "ONSEN_config.R", "ONSEN_functions.R", "00_install_packages.R", "00_run_pipeline.R",
  "01_native_mutated_motif_analysis.R", "02_constrained_mutant_sensitivity.R",
  "03_col0_HSF_and_TE_background.R", "03B_threshold_and_continuous_sensitivity.R",
  "03B_threshold_and_continuous_sensitivity_part1.R", "03B_threshold_and_continuous_sensitivity_part2.R",
  "04_nonredundant_HSF_locations.R", "05_methylation_analysis.R",
  "06_rnaseq_analysis.R", "07_accession_analysis.R",
  "07B_figure7_identity_and_logo_workflow.R", "08_write_supplementary_tables.R",
  "09_write_session_info.R", "10_validate_manuscript_outputs.R"
)
required_metadata <- c(
  "README.md", "LICENSE", "CITATION.cff", "CHANGELOG.md", "DATA_AVAILABILITY.md",
  "INPUT_PROVENANCE.tsv", "REPRODUCIBILITY_MATRIX.tsv", "REPRODUCIBILITY_NOTES.md",
  "FINAL_NUMBERING_MAP.tsv", "FILE_CHECKSUMS_SHA256.tsv", "FULL_PACKAGE_MANIFEST.txt",
  "R_PACKAGE_REQUIREMENTS.txt", "ONSEN_49bp_sequences.fasta", "ONSEN_HSE_units_and_substitutions.csv",
  "ONSEN_Col0_terminal_candidate_windows.csv", "ONSEN_Col0_terminal_candidate_windows.saf",
  "Arabidopsis_HSF_models_JASPAR2026.csv", "source_data/Figure7_variant_metrics.tsv",
  "source_data/Figure7_logo_model_metadata.csv", "source_data/S8A_threshold_stats.tsv",
  "source_data/S8B_continuous_stats.tsv", "source_data/S8D_QC.tsv",
  "source_data/Table_S13_natural_variant_TF_family.tsv",
  "source_data/Table_S14_published_accession_evidence.tsv", "source_data/README.md",
  "motifs/MA1667.2_HSFC1.jaspar", "motifs/MA0981.2_DOF1.8.jaspar", "motifs/README.md",
  "supplementary_table_source/README.md"
)
required_table_sources <- list.files(
  file.path(REPO_ROOT, "supplementary_table_source"),
  pattern = "^Table_S[0-9]+__.*\\.tsv$", recursive = FALSE, full.names = FALSE
)
for (i in 1:14) {
  if (!any(grepl(sprintf("^Table_S%d__", i), required_table_sources))) {
    stop("No deposited source sheet for Table S", i, ".", call. = FALSE)
  }
}
required_tables <- file.path("supplementary_table_source", required_table_sources)
required_files <- c(required_scripts, required_metadata, required_tables)
missing_files <- required_files[!file.exists(file.path(REPO_ROOT, required_files))]
if (length(missing_files)) stop("Repository package is incomplete. Missing:\n", paste(missing_files, collapse = "\n"), call. = FALSE)

# Public-facing R files must not contain hard-coded Windows drive roots.
r_files <- list.files(REPO_ROOT, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
for (path in r_files) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (grepl("(^|[^A-Za-z0-9])([A-Za-z]):[/\\\\]", txt, perl = TRUE)) stop("Hard-coded drive-specific path in: ", basename(path), call. = FALSE)
}

# Numbering declarations.
readme <- paste(readLines(file.path(REPO_ROOT, "README.md"), warn = FALSE), collapse = "\n")
for (token in c("Fig. S1-Fig. S5", "Table S1-Table S14", "MA1667.2", "MA0981.2")) {
  if (!grepl(token, readme, fixed = TRUE)) stop("README missing required token: ", token, call. = FALSE)
}
matrix <- read.delim(file.path(REPO_ROOT, "REPRODUCIBILITY_MATRIX.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
if (!all(sprintf("Table S%d", 1:14) %in% matrix$display_item)) stop("Reproducibility matrix lacks final Tables S1-S14.", call. = FALSE)
for (item in sprintf("Fig. S%d", 1:5)) if (!any(startsWith(matrix$display_item, item))) stop("Reproducibility matrix lacks ", item, call. = FALSE)

# SHA-256 checksum verification using the system sha256sum command.
checksum_file <- file.path(REPO_ROOT, "FILE_CHECKSUMS_SHA256.tsv")
checksums <- read.delim(checksum_file, check.names = FALSE, stringsAsFactors = FALSE)
sha256_one <- function(relative_path) {
  command <- Sys.which("sha256sum")
  if (!nzchar(command)) return(NA_character_)
  output <- system2(command, shQuote(file.path(REPO_ROOT, relative_path)), stdout = TRUE, stderr = TRUE)
  sub("[[:space:]].*$", "", output[[1]])
}
if (nzchar(Sys.which("sha256sum"))) {
  observed <- vapply(checksums$path, sha256_one, character(1))
  if (!all(tolower(observed) == tolower(checksums$sha256))) {
    bad <- checksums$path[tolower(observed) != tolower(checksums$sha256)]
    stop("Checksum mismatch for:\n", paste(bad, collapse = "\n"), call. = FALSE)
  }
} else {
  warning("sha256sum is unavailable; checksum verification was skipped.")
}

read_source_values <- function(table_number) {
  paths <- list.files(file.path(REPO_ROOT, "supplementary_table_source"), pattern = sprintf("^Table_S%d__.*\\.tsv$", table_number), full.names = TRUE)
  unlist(lapply(paths, function(path) unlist(read.delim(path, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE), use.names = FALSE)), use.names = FALSE)
}

numericize_values <- function(values) {
  text_values <- trimws(as.character(values)); text_values[text_values %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  text_values <- gsub("\u00A0", "", text_values, fixed = TRUE); text_values <- gsub(",", "", text_values, fixed = TRUE); text_values <- sub("%$", "", text_values)
  suppressWarnings(as.numeric(text_values))
}
contains_number <- function(values, target, tolerance = max(1e-10, abs(target) * 1e-5)) {
  numeric_values <- numericize_values(values); any(is.finite(numeric_values) & abs(numeric_values - target) <= tolerance)
}
contains_text <- function(values, pattern) any(grepl(pattern, as.character(values), ignore.case = TRUE))
checks <- list()
add_check <- function(name, passed, detail) checks[[length(checks) + 1L]] <<- data.frame(check = name, passed = isTRUE(passed), detail = detail, stringsAsFactors = FALSE)

# Final-numbered source-sheet reference values.
v1 <- read_source_values(1)
add_check("Table S1 lost motif-model count", contains_number(v1, 60), "Expected 60")
add_check("Table S1 gained motif-model count", contains_number(v1, 50), "Expected 50")
v4 <- read_source_values(4)
add_check("Table S4 random-mutant library size", contains_number(v4, 5000), "Expected 5,000")
add_check("Table S4 exact-GC design space", contains_number(v4, 5120), "Expected 5,120")
add_check("Table S4 designed AP2/ERF count", contains_number(v4, 61), "Expected 61")
v5 <- read_source_values(5)
s5_lines <- unlist(lapply(
  list.files(file.path(REPO_ROOT, "supplementary_table_source"), pattern = "^Table_S5__.*\\.tsv$", full.names = TRUE),
  readLines, warn = FALSE
))
add_check("Table S5 has sixteen windows", sum(grepl("^ONSEN[1-8]\\t", s5_lines)) == 16L, "Expected sixteen terminal windows")
add_check("Table S5 includes raw HSF minimum 48", contains_number(v5, 48), "Expected 48")
add_check("Table S5 includes raw HSF maximum 75", contains_number(v5, 75), "Expected 75")
v6 <- read_source_values(6)
add_check("Table S6 native raw HSF count", contains_number(v6, 31), "Expected 31")
add_check("Table S6 native non-redundant count", contains_number(v6, 1), "Expected 1")
v7 <- read_source_values(7)
add_check("Table S7 threshold 0.85 ONSEN median", contains_number(v7, 75), "Expected 75")
add_check("Table S7 Wilcoxon P", contains_number(v7, 7.33e-12, 2e-13), "Expected approximately 7.33e-12")
v8 <- read_source_values(8)
for (x in c(154.4, 75, 37.5, 12.5)) add_check(paste("Table S8 contains ONSEN median", x), contains_number(v8, x, 0.1), paste("Expected", x))
add_check("Table S8 includes threshold 0.95", contains_number(v8, 0.95), "Expected cutoff 0.95")
v9 <- read_source_values(9)
add_check("Table S9 ordinary-TE CHH median", contains_number(v9, 15.83017, 0.02), "Expected approximately 15.83")
add_check("Table S9 ONSEN CHH median", contains_number(v9, 53.61375, 0.02), "Expected approximately 53.61")
v11 <- read_source_values(11)
for (x in c(19, 7, 11, 9, 16, 17, 15, 41)) add_check(paste("Table S11 contains accession count", x), contains_number(v11, x), paste("Expected", x))
v14 <- read_source_values(14)
add_check("Table S14 contains Xu et al. 2024 evidence", contains_text(v14, "Xu et al[.] [(]2024[)]") && contains_text(v14, "10[.]1093/gbe/evae242"), "Expected cited Xu et al. 2024 evidence")
logo <- read.csv(repo_file("source_data/Figure7_logo_model_metadata.csv"), stringsAsFactors = FALSE)
add_check("Figure 7 upper logo is HSFC1 MA1667.2", any(logo$model_name == "HSFC1" & logo$JASPAR_ID == "MA1667.2"), "Expected HSFC1/MA1667.2")
add_check("Figure 7 lower logo is DOF1.8 MA0981.2", any(logo$model_name == "DOF1.8" & logo$JASPAR_ID == "MA0981.2"), "Expected DOF1.8/MA0981.2")

variants <- read.delim(repo_file("source_data/Figure7_variant_metrics.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
add_check("Figure 7 variant source has 53 rows", nrow(variants) == 53L, paste("Rows found:", nrow(variants)))
add_check("Figure 7 variant source has one Col-0 reference", sum(variants$mismatch_count == 0) == 1L, "Expected one zero-mismatch reference")

# Final-audit provenance and numbering checks.
rna_meta <- read.csv(file.path(REPO_ROOT, "RNAseq_sample_metadata_template.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
expected_biosamples <- c("SAMD01789795", "SAMD01789796", "SAMD01789797",
                         "SAMD01943917", "SAMD01943918", "SAMD01943919")
expected_bioprojects <- c(rep("PRJDB39904", 3L), rep("PRJDB42759", 3L))
add_check("RNA-seq template has six Col-0 libraries",
          nrow(rna_meta) == 6L && all(rna_meta$ecotype == "Col-0") &&
            identical(as.integer(rna_meta$replicate), c(1L, 2L, 3L, 1L, 2L, 3L)),
          "Expected exactly six Col-0 rows and three replicates per condition")
add_check("RNA-seq DDBJ BioSample map is exact",
          identical(rna_meta$biosample_accession, expected_biosamples),
          "Unexpected DDBJ BioSample assignment")
add_check("RNA-seq DDBJ BioProject map is exact",
          identical(rna_meta$bioproject, expected_bioprojects),
          "Unexpected DDBJ BioProject assignment")
add_check("RNA-seq laboratory provenance is explicit",
          all(rna_meta$laboratory_provenance == "Generated in our laboratory"),
          "Expected laboratory-generated provenance for all six libraries")

accession_code <- paste(readLines(file.path(REPO_ROOT, "07_accession_analysis.R"),
                                  warn = FALSE), collapse = "\n")
add_check("Final Fig. S4 output naming",
          grepl("FigS4A_accession_HSE_architecture", accession_code, fixed = TRUE),
          "Fig. S4 code/output mapping is not synchronized")
add_check("Final Fig. S5 output naming",
          grepl("FigS5A_accession_candidate_abundance", accession_code, fixed = TRUE),
          "Fig. S5 code/output mapping is not synchronized")

s9d <- file.path(REPO_ROOT,
                 "supplementary_table_source/Table_S9__S9D_selected_outliers.tsv")
add_check("Table S9 exact downstream outlier set is deposited",
          file.exists(s9d) && sum(grepl("^ATHILA|^ATCOPIA", readLines(s9d))) == 8L,
          "Expected eight selected non-ONSEN TE rows")

report <- dplyr::bind_rows(checks)
write.csv(report, file.path(ONSEN_OUTPUT_ROOT, "repository_validation_report.csv"), row.names = FALSE)
if (any(!report$passed)) {
  print(report[!report$passed, ], row.names = FALSE)
  stop("Repository validation failed; see repository_validation_report.csv.", call. = FALSE)
}

message("\n============================================================")
message("REPOSITORY VALIDATION PASSED")
message("============================================================")
message("All final-numbered files, checksums and selected manuscript values passed.")
