# Rebuild final journal-numbered supplementary workbooks from deposited TSV sources.
#
# Run 03B_threshold_and_continuous_sensitivity.R before this script so the
# complete Table S8 region-level sheet is available in ONSEN_OUTPUT_ROOT.

source("ONSEN_functions.R")
require_packages(c("openxlsx", "readr"))
message_config()

SOURCE_DIR <- file.path(REPO_ROOT, "supplementary_table_source")
if (!dir.exists(SOURCE_DIR)) stop("Missing supplementary_table_source directory.", call. = FALSE)

read_tsv_raw <- function(path) {
  readr::read_tsv(path, col_names = FALSE, show_col_types = FALSE, name_repair = "minimal", progress = FALSE)
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

  if (!length(paths)) stop("No source sheets found for Table S", table_number, ".", call. = FALSE)
  if (table_number == 8 && !any(grepl("region_metrics", basename(paths), ignore.case = TRUE))) {
    stop("Table S8 requires its generated region-level sheet. Run 03B_threshold_and_continuous_sensitivity.R first.", call. = FALSE)
  }

  wb <- openxlsx::createWorkbook(creator = "Baibhav R. Barbaruah")
  split_location_paths <- paths[grepl("Table_S6__Col0_ONSEN_locations_part[0-9]+\\.tsv$", basename(paths))]
  if (length(split_location_paths)) {
    split_location_paths <- sort(split_location_paths)
    pieces <- lapply(split_location_paths, read_tsv_raw)
    # Part 1 contains title, blank row and header; later parts contain data only.
    combined <- do.call(rbind, pieces)
    write_raw_sheet(wb, "Col0_ONSEN_locations", combined)
    paths <- setdiff(paths, split_location_paths)
  }
  for (path in paths) {
    sheet <- if (grepl("region_metrics", basename(path), ignore.case = TRUE)) "S8C_region_metrics" else sheet_name_from_path(path)
    x <- read_tsv_raw(path)
    write_raw_sheet(wb, sheet, x)
  }

  destination <- out_file(sprintf("Table_S%d.xlsx", table_number))
  openxlsx::saveWorkbook(wb, destination, overwrite = TRUE)
  message("Rebuilt: ", destination)
}

message("Final supplementary workbooks S1-S14 rebuilt successfully.")
