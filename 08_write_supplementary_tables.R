# Rebuild the deposited PRE-FINAL-NUMBERING supplementary source sheets for
# provenance checks only.
#
# IMPORTANT (August 2026): the final Biology Open submission deleted the old
# AP2/ERF-only Table S3, shifted subsequent table numbers, and added a new global
# gene/TE differential-expression Table S14. The historical TSV files under
# supplementary_table_source/ intentionally retain their earlier numbering.
# FINAL_NUMBERING_MAP.tsv is authoritative for the final journal package.
#
# This script therefore writes to a dedicated legacy_rebuilt_tables directory and
# must not be used to generate the final submission workbooks.

source("ONSEN_functions.R")
require_packages(c("openxlsx", "readr"))
message_config()

SOURCE_DIR <- file.path(REPO_ROOT, "supplementary_table_source")
if (!dir.exists(SOURCE_DIR)) stop("Missing supplementary_table_source directory.", call. = FALSE)

LEGACY_OUT <- file.path(ONSEN_OUTPUT_ROOT, "legacy_rebuilt_tables")
dir.create(LEGACY_OUT, recursive = TRUE, showWarnings = FALSE)

read_tsv_raw <- function(path) {
  readr::read_tsv(path, col_names = FALSE, show_col_types = FALSE,
                  name_repair = "minimal", progress = FALSE)
}

sheet_name_from_path <- function(path) {
  x <- sub("^Table_S[0-9]+__", "", tools::file_path_sans_ext(basename(path)))
  substr(x, 1, 31)
}

write_raw_sheet <- function(wb, sheet, x) {
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, x, colNames = FALSE, rowNames = FALSE, keepNA = TRUE)
  openxlsx::freezePane(wb, sheet, firstActiveRow = 3)
  openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(x)), widths = "auto")
}

message("Rebuilding legacy deposited source sheets only. Final journal numbering is defined in FINAL_NUMBERING_MAP.tsv.")

for (table_number in 1:14) {
  pattern <- sprintf("^Table_S%d__.*\\.tsv$", table_number)
  paths <- sort(list.files(SOURCE_DIR, pattern = pattern, full.names = TRUE))

  if (table_number == 8) {
    generated_candidates <- c(
      out_file("S8C_region_metrics.tsv.gz"),
      out_file("Revision_R1_3_region_metrics.tsv.gz"),
      file.path(ONSEN_DATA_ROOT, "S8C_region_metrics.tsv.gz")
    )
    generated <- generated_candidates[file.exists(generated_candidates)][1]
    if (!is.na(generated)) paths <- c(paths, generated)
  }

  if (!length(paths)) stop("No legacy source sheets found for pre-final Table S", table_number, ".", call. = FALSE)

  wb <- openxlsx::createWorkbook(creator = "Baibhav R. Barbaruah")
  split_location_paths <- paths[grepl("Table_S6__Col0_ONSEN_locations_part[0-9]+\\.tsv$", basename(paths))]
  if (length(split_location_paths)) {
    split_location_paths <- sort(split_location_paths)
    pieces <- lapply(split_location_paths, read_tsv_raw)
    combined <- do.call(rbind, pieces)
    write_raw_sheet(wb, "Col0_ONSEN_locations", combined)
    paths <- setdiff(paths, split_location_paths)
  }

  for (path in paths) {
    sheet <- if (grepl("region_metrics", basename(path), ignore.case = TRUE)) {
      "S8C_region_metrics"
    } else {
      sheet_name_from_path(path)
    }
    x <- read_tsv_raw(path)
    write_raw_sheet(wb, sheet, x)
  }

  destination <- file.path(LEGACY_OUT, sprintf("PRE_FINAL_Table_S%d.xlsx", table_number))
  openxlsx::saveWorkbook(wb, destination, overwrite = TRUE)
  message("Legacy rebuild: ", destination)
}

message("Legacy source-sheet rebuild complete. Do not submit these workbooks as the final journal tables.")
