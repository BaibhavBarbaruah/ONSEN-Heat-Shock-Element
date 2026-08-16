# Validate the public repository metadata/provenance against the final August 2026
# Biology Open revision. Historical files under supplementary_table_source/ retain
# their pre-final-numbering names by design; FINAL_NUMBERING_MAP.tsv is authoritative.

source("ONSEN_config.R")

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required for validation.", call. = FALSE)
}
message_config()

required_scripts <- c(
  "ONSEN_config.R", "ONSEN_functions.R", "00_install_packages.R", "00_run_pipeline.R",
  "01_native_mutated_motif_analysis.R", "02_constrained_mutant_sensitivity.R",
  "03_col0_HSF_and_TE_background.R", "03B_threshold_and_continuous_sensitivity.R",
  "03B_threshold_and_continuous_sensitivity_part1.R", "03B_threshold_and_continuous_sensitivity_part2.R",
  "04_nonredundant_HSF_locations.R", "05_methylation_analysis.R", "06_rnaseq_analysis.R",
  "07_accession_analysis.R", "07B_figure7_identity_and_logo_workflow.R",
  "08_write_supplementary_tables.R", "09_write_session_info.R", "10_validate_manuscript_outputs.R"
)

required_metadata <- c(
  "README.md", "LICENSE", "CITATION.cff", "CHANGELOG.md", "DATA_AVAILABILITY.md",
  "INPUT_PROVENANCE.tsv", "REPRODUCIBILITY_MATRIX.tsv", "FINAL_NUMBERING_MAP.tsv",
  "R_PACKAGE_REQUIREMENTS.txt", "RNAseq_sample_metadata_template.csv",
  "ONSEN_49bp_sequences.fasta", "ONSEN_HSE_units_and_substitutions.csv",
  "ONSEN_Col0_terminal_candidate_windows.csv", "Arabidopsis_HSF_models_JASPAR2026.csv",
  "source_data/Figure7_variant_metrics.tsv", "source_data/Figure7_logo_model_metadata.csv",
  "motifs/MA1667.2_HSFC1.jaspar", "motifs/MA0981.2_DOF1.8.jaspar",
  "supplementary_table_source/README.md"
)

required_files <- c(required_scripts, required_metadata)
missing_files <- required_files[!file.exists(file.path(REPO_ROOT, required_files))]
if (length(missing_files)) {
  stop("Repository package is incomplete. Missing:\n", paste(missing_files, collapse = "\n"), call. = FALSE)
}

# Public-facing R files must not contain hard-coded Windows drive roots.
r_files <- list.files(REPO_ROOT, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
for (path in r_files) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (grepl("(^|[^A-Za-z0-9])([A-Za-z]):[/\\\\]", txt, perl = TRUE)) {
    stop("Hard-coded drive-specific path in: ", basename(path), call. = FALSE)
  }
}

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, passed = isTRUE(passed), detail = detail, stringsAsFactors = FALSE
  )
}

# README declarations.
readme <- paste(readLines(file.path(REPO_ROOT, "README.md"), warn = FALSE), collapse = "\n")
add_check("README final supplementary figures", grepl("Fig. S1-Fig. S4", readme, fixed = TRUE), "Expected Fig. S1-Fig. S4")
add_check("README final supplementary tables", grepl("Table S1-Table S14", readme, fixed = TRUE), "Expected Table S1-Table S14")
add_check("README marks legacy source numbering", grepl("Legacy source-sheet numbering", readme, fixed = TRUE), "Expected explicit legacy-source note")
add_check("README records principal exact-coordinate metric", grepl("physical-forward exact-coordinate", readme, fixed = TRUE), "Expected principal metric definition")
add_check("README records final global DE table", grepl("Table S14", readme, fixed = TRUE) && grepl("genome-wide gene", readme, ignore.case = TRUE), "Expected global gene/TE Table S14")

# Final numbering map.
numbering <- read.delim(file.path(REPO_ROOT, "FINAL_NUMBERING_MAP.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
final_tables <- numbering$final_item[numbering$item_type == "Supplementary table"]
final_figs <- numbering$final_item[numbering$item_type == "Supplementary figure"]
add_check("Final map has Tables S1-S14", all(sprintf("Table S%d", 1:14) %in% final_tables), "Expected every final Table S1-S14")
add_check("Final map deletes AP2-only Table S3", any(numbering$pre_final_item == "Table S3" & numbering$final_item == "DELETED"), "Expected old AP2-only S3 deletion")
add_check("Final map has Figures S1-S4", all(sprintf("Fig. S%d", 1:4) %in% final_figs), "Expected every final Fig. S1-S4")
add_check("Final map has no retained Fig. S5", !any(final_figs == "Fig. S5"), "Final supplement should stop at Fig. S4")
add_check("Final map assigns new global DE to S14", any(grepl("NEW global gene/TE DE", numbering$pre_final_item, fixed = TRUE) & numbering$final_item == "Table S14"), "Expected new global DE as final Table S14")

# Reproducibility matrix declarations.
matrix <- read.delim(file.path(REPO_ROOT, "REPRODUCIBILITY_MATRIX.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
add_check("Matrix has Tables S1-S14", all(sprintf("Table S%d", 1:14) %in% matrix$display_item), "Expected every final Table S1-S14")
for (item in sprintf("Fig. S%d", 1:4)) {
  add_check(paste("Matrix has", item), any(matrix$display_item == item), paste("Expected", item))
}
add_check("Matrix excludes final Fig. S5", !any(matrix$display_item == "Fig. S5"), "Expected no final Fig. S5")
add_check("Matrix includes Fig. 5D gene volcano", any(matrix$display_item == "Fig. 5D"), "Expected global gene volcano")
add_check("Matrix includes Fig. 5E TE volcano", any(matrix$display_item == "Fig. 5E"), "Expected global TE volcano")

# RNA-seq accession/provenance map.
rna_meta <- read.csv(file.path(REPO_ROOT, "RNAseq_sample_metadata_template.csv"), check.names = FALSE, stringsAsFactors = FALSE)
expected_biosamples <- c(
  "SAMD01789795", "SAMD01789796", "SAMD01789797",
  "SAMD01943917", "SAMD01943918", "SAMD01943919"
)
expected_bioprojects <- c(rep("PRJDB39904", 3L), rep("PRJDB42759", 3L))
add_check("RNA-seq template has six Col-0 libraries", nrow(rna_meta) == 6L && all(rna_meta$ecotype == "Col-0"), "Expected six Col-0 libraries")
add_check("RNA-seq DDBJ BioSample map is exact", identical(rna_meta$biosample_accession, expected_biosamples), "Unexpected BioSample assignment")
add_check("RNA-seq DDBJ BioProject map is exact", identical(rna_meta$bioproject, expected_bioprojects), "Unexpected BioProject assignment")

# Figure 7 logo metadata remains intentionally illustrative.
logo <- read.csv(file.path(REPO_ROOT, "source_data/Figure7_logo_model_metadata.csv"), stringsAsFactors = FALSE)
add_check("Figure 7 upper logo is HSFC1 MA1667.2", any(logo$model_name == "HSFC1" & logo$JASPAR_ID == "MA1667.2"), "Expected HSFC1/MA1667.2")
add_check("Figure 7 lower logo is DOF1.8 MA0981.2", any(logo$model_name == "DOF1.8" & logo$JASPAR_ID == "MA0981.2"), "Expected DOF1.8/MA0981.2")

report <- dplyr::bind_rows(checks)
write.csv(report, file.path(ONSEN_OUTPUT_ROOT, "repository_validation_report.csv"), row.names = FALSE)

if (any(!report$passed)) {
  print(report[!report$passed, ], row.names = FALSE)
  stop("Repository validation failed; see repository_validation_report.csv.", call. = FALSE)
}

message("\n============================================================")
message("REPOSITORY METADATA VALIDATION PASSED")
message("============================================================")
message("Final August 2026 numbering, provenance and public metadata passed.")
