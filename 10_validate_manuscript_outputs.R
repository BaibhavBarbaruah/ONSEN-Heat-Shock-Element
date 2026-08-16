# Validate repository metadata and deposited tables against the final August
# 2026 Biology Open package: Figs S1-S4 and Tables S1-S13.

source("ONSEN_config.R")

for (package in c("dplyr", "openxlsx")) {
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
    detail = detail,
    stringsAsFactors = FALSE
  )
}

required_scripts <- c(
  "ONSEN_config.R", "ONSEN_functions.R", "00_install_packages.R",
  "00_run_pipeline.R", "01_native_mutated_motif_analysis.R",
  "02_constrained_mutant_sensitivity.R", "03_col0_HSF_and_TE_background.R",
  "03B_threshold_and_continuous_sensitivity.R",
  "03B_threshold_and_continuous_sensitivity_part1.R",
  "03B_threshold_and_continuous_sensitivity_part2.R",
  "03C_direct_LTR_HSF_comparison.R", "04_nonredundant_HSF_locations.R",
  "05A_direct_LTR_methylation_analysis.R", "05_methylation_analysis.R",
  "06_rnaseq_analysis.R", "07_accession_analysis.R",
  "07B_figure7_identity_and_logo_workflow.R",
  "08_write_supplementary_tables.R", "09_write_session_info.R",
  "10_validate_manuscript_outputs.R"
)

required_metadata <- c(
  "README.md", "LICENSE", "CITATION.cff", "CHANGELOG.md",
  "DATA_AVAILABILITY.md", "REPRODUCIBILITY_NOTES.md",
  "INPUT_PROVENANCE.tsv", "REPRODUCIBILITY_MATRIX.tsv",
  "FINAL_NUMBERING_MAP.tsv", "R_PACKAGE_REQUIREMENTS.txt",
  "RNAseq_sample_metadata_template.csv", "ONSEN_49bp_sequences.fasta",
  "ONSEN_HSE_units_and_substitutions.csv",
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
  "Required final files are present",
  !length(missing_files),
  if (length(missing_files)) paste(missing_files, collapse = "; ") else "complete"
)

add_check(
  "No obsolete root Table S14",
  !file.exists(file.path(REPO_ROOT, "Table_S14.xlsx")),
  "Final package stops at Table S13"
)
add_check(
  "No obsolete AP2/ERF source sheet",
  !file.exists(file.path(
    REPO_ROOT, "supplementary_table_source", "Table_S3__S3_AP2ERF_gained.tsv"
  )),
  "AP2/ERF-only workbook was deleted"
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

readme <- paste(readLines(file.path(REPO_ROOT, "README.md"), warn = FALSE), collapse = "\n")
add_check("README has Figs S1-S4", grepl("Fig. S1-Fig. S4", readme, fixed = TRUE), "expected final figures")
add_check("README has Tables S1-S13", grepl("Table S1-Table S13", readme, fixed = TRUE), "expected final tables")
add_check("README states no final S14", grepl("no final Table S14", readme, ignore.case = TRUE), "expected explicit deletion")
add_check("README assigns global DE to S9", grepl("final \\*\\*Table S9\\*\\*", readme), "expected global DE Table S9")

numbering <- read.delim(
  file.path(REPO_ROOT, "FINAL_NUMBERING_MAP.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
final_tables <- numbering$final_item[
  numbering$item_type == "Supplementary table" & numbering$final_item != "DELETED"
]
final_figs <- numbering$final_item[
  numbering$item_type == "Supplementary figure" & numbering$final_item != "DELETED"
]
add_check(
  "Final map has exactly Tables S1-S13",
  setequal(final_tables, sprintf("Table S%d", 1:13)) && length(final_tables) == 13L,
  paste(final_tables, collapse = "; ")
)
add_check(
  "Final map has exactly Figs S1-S4",
  setequal(final_figs, sprintf("Fig. S%d", 1:4)) && length(final_figs) == 4L,
  paste(final_figs, collapse = "; ")
)
add_check(
  "Final map assigns global DE to S9",
  any(numbering$pre_final_item == "NEW global gene/TE DE" &
        numbering$final_item == "Table S9"),
  "expected new global DE Table S9"
)

matrix <- read.delim(
  file.path(REPO_ROOT, "REPRODUCIBILITY_MATRIX.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "Matrix contains Tables S1-S13",
  all(sprintf("Table S%d", 1:13) %in% matrix$display_item),
  "all final tables mapped"
)
add_check(
  "Matrix excludes Table S14",
  !any(matrix$display_item == "Table S14"),
  "no obsolete final table"
)
add_check(
  "Matrix contains Figs S1-S4 only",
  all(sprintf("Fig. S%d", 1:4) %in% matrix$display_item) &&
    !any(matrix$display_item == "Fig. S5"),
  "final supplementary figures mapped"
)

source_manifest <- read.delim(
  file.path(REPO_ROOT, "supplementary_table_source", "SOURCE_SHEET_MANIFEST.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "Final TSV source manifest has 36 worksheets",
  nrow(source_manifest) == 36L,
  paste("observed", nrow(source_manifest))
)
add_check(
  "Final TSV source manifest covers S1-S13",
  all(sprintf("Table_S%d.xlsx", 1:13) %in% source_manifest$workbook),
  "all final workbooks mirrored"
)

# Verify every final workbook opens and has at least one sheet.
unreadable_tables <- character()
for (table_name in required_tables) {
  sheets <- tryCatch(
    openxlsx::getSheetNames(file.path(REPO_ROOT, table_name)),
    error = function(e) character()
  )
  if (!length(sheets)) unreadable_tables <- c(unreadable_tables, table_name)
}
add_check(
  "All final workbooks open",
  !length(unreadable_tables),
  if (length(unreadable_tables)) paste(unreadable_tables, collapse = "; ") else "13/13"
)

s10_path <- file.path(REPO_ROOT, "Table_S10.xlsx")
s10_sheets <- openxlsx::getSheetNames(s10_path)
add_check(
  "Table S10 contains statistics sheet",
  "S10D_Statistics" %in% s10_sheets,
  paste(s10_sheets, collapse = "; ")
)
s10a <- openxlsx::read.xlsx(
  s10_path,
  sheet = "S10A_Class_summary",
  colNames = FALSE,
  skipEmptyRows = FALSE
)
expected_s10_headers <- c(
  "Candidate class", "n windows", "Mean NS CPM", "Mean HS CPM",
  "Median NS CPM", "Median HS CPM",
  "Descriptive log2FC [(mean HS CPM + 0.05)/(mean NS CPM + 0.05)]"
)
observed_s10_headers <- as.character(unlist(s10a[2, 1:7], use.names = FALSE))
add_check(
  "Table S10 class-summary headers restored",
  identical(observed_s10_headers, expected_s10_headers),
  paste(observed_s10_headers, collapse = "; ")
)
s10d <- openxlsx::read.xlsx(
  s10_path,
  sheet = "S10D_Statistics",
  colNames = FALSE,
  skipEmptyRows = FALSE
)
s10d_text <- paste(unlist(s10d, use.names = FALSE), collapse = " ")
add_check(
  "Table S10 reports Welch and Wilcoxon results",
  grepl("0.000757270061", s10d_text, fixed = TRUE) &&
    grepl("0.08085559837", s10d_text, fixed = TRUE),
  "expected two-sided P values"
)

rna_meta <- read.csv(
  file.path(REPO_ROOT, "RNAseq_sample_metadata_template.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
expected_biosamples <- c(
  "SAMD01789795", "SAMD01789796", "SAMD01789797",
  "SAMD01943917", "SAMD01943918", "SAMD01943919"
)
expected_bioprojects <- c(rep("PRJDB39904", 3L), rep("PRJDB42759", 3L))
add_check(
  "RNA-seq metadata has six Col-0 libraries",
  nrow(rna_meta) == 6L && all(rna_meta$ecotype == "Col-0"),
  "expected six validated libraries"
)
add_check(
  "RNA-seq BioSample map is exact",
  identical(rna_meta$biosample_accession, expected_biosamples),
  "validated DDBJ assignments"
)
add_check(
  "RNA-seq BioProject map is exact",
  identical(rna_meta$bioproject, expected_bioprojects),
  "validated DDBJ projects"
)

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

text_files <- c("README.md", "CITATION.cff", "DATA_AVAILABILITY.md")
wrong_spelling <- vapply(text_files, function(path) {
  any(grepl(
    "Airalangga",
    readLines(file.path(REPO_ROOT, path), warn = FALSE),
    fixed = TRUE
  ))
}, logical(1))
add_check(
  "Author surname is consistently Airlangga",
  !any(wrong_spelling),
  if (any(wrong_spelling)) paste(text_files[wrong_spelling], collapse = "; ") else "correct"
)

report <- dplyr::bind_rows(checks)
write.csv(
  report,
  file.path(ONSEN_OUTPUT_ROOT, "repository_validation_report.csv"),
  row.names = FALSE
)

if (any(!report$passed)) {
  print(report[!report$passed, ], row.names = FALSE)
  stop(
    "Repository validation failed; see repository_validation_report.csv.",
    call. = FALSE
  )
}

message("\n============================================================")
message("FINAL S1-S13 REPOSITORY VALIDATION PASSED")
message("============================================================")
