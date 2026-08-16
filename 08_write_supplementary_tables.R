# Rebuild data-equivalent final-numbered supplementary workbooks from the
# deposited TSV mirrors. The formatted root Table_S1.xlsx-Table_S13.xlsx files
# remain the authoritative journal submission artifacts.

source("ONSEN_functions.R")
require_packages(c("openxlsx", "readr"))
message_config()

SOURCE_DIR <- file.path(REPO_ROOT, "supplementary_table_source")
if (!dir.exists(SOURCE_DIR)) {
  stop("Missing supplementary_table_source directory.", call. = FALSE)
}

REBUILD_OUT <- file.path(ONSEN_OUTPUT_ROOT, "rebuilt_final_data_workbooks")
dir.create(REBUILD_OUT, recursive = TRUE, showWarnings = FALSE)

read_tsv_raw <- function(path) {
  readr::read_tsv(
    path,
    col_names = FALSE,
    show_col_types = FALSE,
    name_repair = "minimal",
    progress = FALSE
  )
}

sheet_name_from_path <- function(path) {
  x <- sub("^Table_S[0-9]+__", "", tools::file_path_sans_ext(basename(path)))
  substr(x, 1, 31)
}

write_raw_sheet <- function(workbook, sheet_name, values) {
  openxlsx::addWorksheet(workbook, sheet_name)
  openxlsx::writeData(
    workbook,
    sheet_name,
    values,
    colNames = FALSE,
    rowNames = FALSE,
    keepNA = TRUE
  )
  openxlsx::setColWidths(
    workbook,
    sheet_name,
    cols = seq_len(ncol(values)),
    widths = "auto"
  )
}

for (table_number in 1:13) {
  pattern <- sprintf("^Table_S%d__.*\\.tsv$", table_number)
  paths <- sort(list.files(SOURCE_DIR, pattern = pattern, full.names = TRUE))
  if (!length(paths)) {
    stop("No final source sheets found for Table S", table_number, ".", call. = FALSE)
  }

  workbook <- openxlsx::createWorkbook(creator = "Baibhav R. Barbaruah")
  for (path in paths) {
    write_raw_sheet(workbook, sheet_name_from_path(path), read_tsv_raw(path))
  }

  destination <- file.path(
    REBUILD_OUT,
    sprintf("Table_S%d_data_rebuild.xlsx", table_number)
  )
  openxlsx::saveWorkbook(workbook, destination, overwrite = TRUE)
  message("Data-equivalent rebuild: ", destination)
}

message(
  "Final S1-S13 data rebuild complete. Use the formatted root workbooks for ",
  "journal submission."
)
