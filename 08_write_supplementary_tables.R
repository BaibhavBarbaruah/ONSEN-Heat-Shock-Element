# Rebuild content-equivalent final supplementary workbooks from the deposited
# worksheet-level TSV exports. The tracked root workbooks remain authoritative
# for final formatting, merged cells, widths and numeric cell types.

source("ONSEN_functions.R")
require_packages(c("openxlsx", "readr"))
message_config()

SOURCE_DIR <- file.path(REPO_ROOT, "supplementary_table_source")
OUTPUT_DIR <- file.path(ONSEN_OUTPUT_ROOT, "rebuilt_final_tables")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(SOURCE_DIR)) {
  stop("Missing supplementary_table_source directory.", call. = FALSE)
}

read_tsv_raw <- function(path) {
  readr::read_tsv(
    path,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character()),
    na = character(),
    name_repair = "minimal",
    progress = FALSE,
    skip_empty_rows = FALSE,
    show_col_types = FALSE
  )
}

sheet_name_from_path <- function(path) {
  value <- sub("^Table_S[0-9]+__", "", tools::file_path_sans_ext(basename(path)))
  substr(value, 1, 31)
}

write_source_sheet <- function(workbook, sheet_name, values) {
  openxlsx::addWorksheet(workbook, sheet_name)
  openxlsx::writeData(
    workbook, sheet_name, values,
    colNames = FALSE, rowNames = FALSE, keepNA = FALSE
  )

  n_columns <- max(1L, ncol(values))
  if (n_columns > 1L) openxlsx::mergeCells(workbook, sheet_name, cols = seq_len(n_columns), rows = 1L)

  title_style <- openxlsx::createStyle(
    fontSize = 12, textDecoration = "bold",
    fgFill = "#D9EAF7", wrapText = TRUE, valign = "center"
  )
  openxlsx::addStyle(
    workbook, sheet_name, title_style,
    rows = 1L, cols = seq_len(n_columns), gridExpand = TRUE, stack = TRUE
  )
  openxlsx::setRowHeights(workbook, sheet_name, rows = 1L, heights = 36)
  openxlsx::setColWidths(workbook, sheet_name, cols = seq_len(n_columns), widths = "auto")
  openxlsx::freezePane(workbook, sheet_name, firstActiveRow = min(4L, nrow(values) + 1L))
}

expected_tables <- sprintf("Table_S%d.xlsx", 1:13)
for (table_number in 1:13) {
  pattern <- sprintf("^Table_S%d__.*\\.tsv$", table_number)
  paths <- sort(list.files(SOURCE_DIR, pattern = pattern, full.names = TRUE))
  if (!length(paths)) {
    stop("No final source sheets found for Table S", table_number, ".", call. = FALSE)
  }

  workbook <- openxlsx::createWorkbook(creator = "Baibhav R. Barbaruah")
  for (path in paths) {
    write_source_sheet(workbook, sheet_name_from_path(path), read_tsv_raw(path))
  }

  destination <- file.path(OUTPUT_DIR, expected_tables[table_number])
  openxlsx::saveWorkbook(workbook, destination, overwrite = TRUE)
  message("Rebuilt final Table S", table_number, ": ", destination)
}

message("Final supplementary-table content rebuild complete: ", OUTPUT_DIR)
