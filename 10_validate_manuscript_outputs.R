# Validate the public repository against the final Biology Open manuscript data
# and supplementary package (Fig. S1-Fig. S4; Table S1-Table S13).

source("ONSEN_config.R")
needed <- c("data.table", "dplyr", "openxlsx")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing validation package(s): ", paste(missing, collapse = ", "))
message_config()

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, passed = isTRUE(passed), detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}
text_of <- function(path) {
  full_path <- if (grepl("^(/|[A-Za-z]:[/\\\\])", path, perl = TRUE)) path else file.path(REPO_ROOT, path)
  paste(readLines(full_path, warn = FALSE), collapse = "\n")
}
read_final <- function(filename, skip = 2L) {
  data.table::fread(
    file.path(REPO_ROOT, "supplementary_table_source", filename),
    skip = skip, data.table = FALSE, check.names = FALSE
  )
}

required_scripts <- c(
  "ONSEN_config.R", "ONSEN_functions.R", "00_install_packages.R", "00_run_pipeline.R",
  "01_native_mutated_motif_analysis.R", "02_constrained_mutant_sensitivity.R",
  "03_col0_HSF_and_TE_background.R", "03B_threshold_and_continuous_sensitivity.R",
  "03C_direct_LTR_HSF_comparison.R", "04_nonredundant_HSF_locations.R",
  "05A_direct_LTR_methylation_analysis.R", "05_methylation_analysis.R",
  "06_rnaseq_analysis.R", "07_accession_analysis.R",
  "07B_figure7_identity_and_logo_workflow.R", "08_write_supplementary_tables.R",
  "09_write_session_info.R", "10_validate_manuscript_outputs.R"
)
required_docs <- c(
  "README.md", "DATA_AVAILABILITY.md", "REPRODUCIBILITY_NOTES.md",
  "REPRODUCIBILITY_MATRIX.tsv", "INPUT_PROVENANCE.tsv", "R_PACKAGE_REQUIREMENTS.txt",
  "FULL_PACKAGE_MANIFEST.txt", "CITATION.cff", "LICENSE",
  "RNAseq_sample_metadata_template.csv", "ONSEN_49bp_sequences.fasta",
  "ONSEN_HSE_units_and_substitutions.csv", "ONSEN_Col0_terminal_candidate_windows.csv",
  "Arabidopsis_HSF_models_JASPAR2026.csv",
  "source_data/Figure4A_direct_LTR_methylation_summary.tsv",
  "source_data/Figure4C_direct_LTR_HSF_summary.tsv",
  "source_data/Figure7_variant_metrics.tsv", "source_data/Figure7_logo_model_metadata.csv",
  "motifs/MA1667.2_HSFC1.jaspar", "motifs/MA0981.2_DOF1.8.jaspar",
  "supplementary_table_source/SOURCE_SHEET_MANIFEST.tsv"
)
required_tables <- sprintf("Table_S%d.xlsx", 1:13)
required <- c(required_scripts, required_docs, required_tables)
missing_files <- required[!file.exists(file.path(REPO_ROOT, required))]
add_check("Required public files", !length(missing_files),
          if (length(missing_files)) paste(missing_files, collapse = "; ") else "complete")

obsolete <- c(
  "03B_threshold_and_continuous_sensitivity_part1.R",
  "03B_threshold_and_continuous_sensitivity_part2.R",
  "DDBJ_Update_Change_Log.md", "FINAL_NUMBERING_MAP.tsv", "CHANGELOG.md",
  "FILE_CHECKSUMS_SHA256.tsv"
)
present_obsolete <- obsolete[file.exists(file.path(REPO_ROOT, obsolete))]
add_check("No obsolete/internal tracked files", !length(present_obsolete),
          if (length(present_obsolete)) paste(present_obsolete, collapse = "; ") else "clean")

# Portability.
r_files <- list.files(REPO_ROOT, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
drive_specific <- vapply(r_files, function(path) {
  grepl("(^|[^A-Za-z0-9])([A-Za-z]):[/\\\\]", text_of(path), perl = TRUE)
}, logical(1))
add_check("No hard-coded Windows drive roots", !any(drive_specific),
          if (any(drive_specific)) paste(basename(r_files[drive_specific]), collapse = "; ") else "portable")

# Documentation/provenance synchronization.
readme <- text_of("README.md")
data_avail <- text_of("DATA_AVAILABILITY.md")
provenance <- text_of("INPUT_PROVENANCE.tsv")
packages <- text_of("R_PACKAGE_REQUIREMENTS.txt")
metadata <- paste(readme, data_avail, provenance, text_of("REPRODUCIBILITY_NOTES.md"),
                  text_of("REPRODUCIBILITY_MATRIX.tsv"), text_of("CITATION.cff"))
add_check("Final supplementary scope",
          grepl("Fig. S1-Fig. S4", readme, fixed = TRUE) &&
            grepl("Table S1-Table S13", readme, fixed = TRUE),
          "Figs S1-S4; Tables S1-S13")
add_check("Exact-coordinate HSF metric documented",
          grepl("unique `(start, end)` coordinate placements", readme, fixed = TRUE) &&
            grepl("identical start and end coordinates are counted once", readme, fixed = TRUE),
          "physical start/end coordinate deduplication")
add_check("Final strict-TE comparison documented",
          grepl("1,930", readme, fixed = TRUE) && grepl("49.375", readme, fixed = TRUE) &&
            grepl("6.81238953796413", readme, fixed = TRUE),
          "16 ONSEN windows versus 1,930 strict non-ONSEN TE regions")
add_check("Final RNA-seq methods documented",
          all(vapply(c("32,833", "31,189", "2,101", "limma-voom",
                       "summed across the sixteen curated ONSEN terminal windows"),
                     function(x) grepl(x, readme, fixed = TRUE), logical(1))),
          "DESeq2 genes; filtered limma-voom TEs; replicate-summed Fig. 5B")
add_check("Ensembl Plants release 63 provenance",
          grepl("Ensembl Plants release 63", provenance, fixed = TRUE) &&
            grepl("Arabidopsis_thaliana.TAIR10.63.gtf", provenance, fixed = TRUE),
          "TAIR10 release-63 gene annotation")
add_check("TE-analysis dependencies documented",
          grepl("edgeR", packages, fixed = TRUE) && grepl("limma", packages, fixed = TRUE),
          "edgeR and limma")
add_check("Author surname",
          grepl("Rahmadani P. Airlangga", metadata, fixed = TRUE) &&
            !grepl("Airalangga", metadata, fixed = TRUE),
          "Rahmadani P. Airlangga")

# DDBJ map.
rna_meta <- read.csv(file.path(REPO_ROOT, "RNAseq_sample_metadata_template.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
expected_samples <- c("SAMD01789795", "SAMD01789796", "SAMD01789797",
                      "SAMD01943917", "SAMD01943918", "SAMD01943919")
expected_projects <- c(rep("PRJDB39904", 3L), rep("PRJDB42759", 3L))
add_check("Six-library DDBJ map",
          nrow(rna_meta) == 6L && identical(rna_meta$biosample_accession, expected_samples) &&
            identical(rna_meta$bioproject, expected_projects) &&
            identical(as.integer(rna_meta$replicate), c(1L,2L,3L,1L,2L,3L)) &&
            all(rna_meta$laboratory_provenance == "Generated in our laboratory"),
          "three study-generated NS and three study-generated HS libraries")

# Workbooks/source-sheet manifest.
manifest <- read.delim(file.path(REPO_ROOT, "supplementary_table_source", "SOURCE_SHEET_MANIFEST.tsv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
add_check("Final source-sheet manifest", nrow(manifest) == 36L &&
            all(required_tables %in% manifest$workbook), "36 worksheets across S1-S13")
unreadable <- required_tables[!vapply(required_tables, function(x) {
  length(tryCatch(openxlsx::getSheetNames(file.path(REPO_ROOT, x)), error = function(e) character())) > 0L
}, logical(1))]
add_check("All final workbooks open", !length(unreadable),
          if (length(unreadable)) paste(unreadable, collapse = "; ") else "13/13")
s10_sheets <- openxlsx::getSheetNames(file.path(REPO_ROOT, "Table_S10.xlsx"))
add_check("Table S10 has all four final sheets",
          identical(s10_sheets, c("S10A_Class_summary", "S10B_Individual_windows",
                                  "S10C_Replicate_CPM", "S10D_Statistics")),
          paste(s10_sheets, collapse = "; "))
s10_text <- readLines(file.path(REPO_ROOT, "supplementary_table_source",
                                "Table_S10__S10D_Statistics.tsv"), warn = FALSE)
add_check("Table S10 statistics",
          any(grepl("0.000757270061088614", s10_text, fixed = TRUE)) &&
            any(grepl("0.08085559837005224", s10_text, fixed = TRUE)) &&
            any(grepl("Replicate-summed fractional CPM", s10_text, fixed = TRUE)),
          "Welch and Wilcoxon on three replicate sums per condition")

# Table S7 exact-coordinate statistics.
s7 <- read_final("Table_S7__S7A_threshold_stats.tsv")
s7_regions <- read_final("Table_S7__S7D_region_threshold.tsv")
k85 <- s7[s7[["Relative PWM-score threshold"]] == 0.85, , drop = FALSE]
add_check("Table S7 final 0.85 comparison",
          nrow(k85) == 1L &&
            as.integer(k85[["n ONSEN terminal windows"]]) == 16L &&
            as.integer(k85[["n strict non-ONSEN TE regions"]]) == 1930L &&
            abs(k85[["ONSEN median density (placements/kb)"]] - 49.375) < 1e-12 &&
            abs(k85[["Background median density (placements/kb)"]] - 6.81238953796413) < 1e-12,
          "n=16/1930; medians=49.375/6.81238953796413")
add_check("Table S7 regional source completeness",
          sum(s7_regions[["Region class"]] == "ONSEN terminal candidate windows" &
                s7_regions[["Relative PWM-score threshold"]] == 0.85) == 16L &&
            sum(s7_regions[["Region class"]] == "Strict non-ONSEN TE background" &
                  s7_regions[["Relative PWM-score threshold"]] == 0.85) == 1930L &&
            "Non-redundant HSF motif-coordinate placements" %in% names(s7_regions),
          "exact-coordinate rows complete")

# Table S9 exact totals and canonical ONSEN extraction.
s9 <- read_final("Table_S9__Summary.tsv")
g <- s9[s9[[1]] == "Genes", , drop = FALSE]
t <- s9[s9[[1]] == "Transposable elements", , drop = FALSE]
add_check("Table S9 gene totals",
          nrow(g) == 1L && identical(as.integer(unlist(g[1, 2:6], use.names = FALSE)),
                                    c(32833L,25912L,23922L,3912L,4185L)),
          "32833/25912/23922/3912/4185")
add_check("Table S9 TE totals",
          nrow(t) == 1L && identical(as.integer(unlist(t[1, 2:6], use.names = FALSE)),
                                    c(31189L,7632L,2101L,692L,397L)),
          "31189/7632/2101/692/397")
s9_onsen <- read_final("Table_S9__S9C_ONSEN.tsv")
add_check("Eight canonical ONSEN loci are Up",
          nrow(s9_onsen) == 8L && "DE call" %in% names(s9_onsen) && all(s9_onsen[["DE call"]] == "Up"),
          "8/8 canonical loci")

# Table S11/S12 accession abundance consistency.
s11 <- read_final("Table_S11__S11_accessions.tsv", skip = 1L)
s12 <- read_final("Table_S12__S12_seed_variants.tsv", skip = 1L)
expected_accession_counts <- c(19L,7L,11L,9L,16L,17L,15L,41L)
add_check("Accession candidate counts",
          identical(as.integer(s11[["Candidate windows"]]), expected_accession_counts) &&
            identical(as.integer(s12[["Candidate windows"]]), expected_accession_counts),
          "Col-0/An-1/C24/Cvi/Eri/Kyo/Ler/Sha = 19/7/11/9/16/17/15/41")

# Code-level guards for substantive analysis definitions.
s3b <- text_of("03B_threshold_and_continuous_sensitivity.R")
s6 <- text_of("06_rnaseq_analysis.R")
s7code <- text_of("07_accession_analysis.R")
add_check("03B uses exact-coordinate deduplication",
          grepl("unique(coordinate_sets[[t]])", s3b, fixed = TRUE) &&
            grepl("forward_start", s3b, fixed = TRUE) && grepl("forward_end", s3b, fixed = TRUE),
          "no raw model-position count substituted")
add_check("06 implements final global TE method",
          grepl("edgeR::filterByExpr", s6, fixed = TRUE) &&
            grepl("limma::voom", s6, fixed = TRUE) && grepl("2101L", s6, fixed = TRUE),
          "expression filtering + limma-voom")
add_check("06 implements replicate-summed Fig. 5B",
          grepl("summed_fractional_CPM = sum(fractional_CPM)", s6, fixed = TRUE) &&
            grepl("stats::t.test(summed_fractional_CPM ~ treatment", s6, fixed = TRUE),
          "sum across 16 windows within each biological replicate")
add_check("07 does not substitute raw PWM-model-position hits",
          !grepl("HSF_hits = dplyr::n()", s7code, fixed = TRUE) &&
            grepl("exact-coordinate", s7code, fixed = TRUE),
          "final accession exact-coordinate source is used")

logo <- read.csv(file.path(REPO_ROOT, "source_data", "Figure7_logo_model_metadata.csv"),
                 stringsAsFactors = FALSE)
add_check("Figure 7 logo models",
          any(logo$model_name == "HSFC1" & logo$JASPAR_ID == "MA1667.2") &&
            any(logo$model_name == "DOF1.8" & logo$JASPAR_ID == "MA0981.2"),
          "HSFC1 MA1667.2; DOF1.8 MA0981.2")

report <- dplyr::bind_rows(checks)
dir.create(ONSEN_OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
write.csv(report, file.path(ONSEN_OUTPUT_ROOT, "repository_validation_report.csv"), row.names = FALSE)
if (any(!report$passed)) {
  print(report[!report$passed, ], row.names = FALSE)
  stop("Repository validation failed; see repository_validation_report.csv.", call. = FALSE)
}
message("\n============================================================")
message("FINAL PUBLIC REPOSITORY VALIDATION PASSED")
message("============================================================")
