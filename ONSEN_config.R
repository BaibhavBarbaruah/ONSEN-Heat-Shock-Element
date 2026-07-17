# Portable configuration for the ONSEN HSE repository.

options(stringsAsFactors = FALSE)

REPO_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

ONSEN_DATA_ROOT <- Sys.getenv("ONSEN_DATA_ROOT", unset = REPO_ROOT)
ONSEN_DATA_ROOT <- normalizePath(ONSEN_DATA_ROOT, winslash = "/", mustWork = FALSE)

ONSEN_OUTPUT_ROOT <- Sys.getenv(
  "ONSEN_OUTPUT_ROOT",
  unset = file.path(REPO_ROOT, "reproduced_outputs")
)
ONSEN_OUTPUT_ROOT <- normalizePath(ONSEN_OUTPUT_ROOT, winslash = "/", mustWork = FALSE)
dir.create(ONSEN_OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

ONSEN_FORCE_RESCAN <- identical(tolower(Sys.getenv("ONSEN_FORCE_RESCAN", "false")), "true")
ONSEN_RUN_LARGE_STEPS <- identical(tolower(Sys.getenv("ONSEN_RUN_LARGE_STEPS", "false")), "true")
ONSEN_MAKE_FIGURES <- !identical(tolower(Sys.getenv("ONSEN_MAKE_FIGURES", "true")), "false")

repo_file <- function(filename) {
  path <- file.path(REPO_ROOT, filename)
  if (!file.exists(path)) stop("Repository file not found: ", path, call. = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

find_input <- function(filename, required = TRUE, extra = character()) {
  candidates <- unique(c(
    file.path(REPO_ROOT, filename),
    file.path(ONSEN_DATA_ROOT, filename),
    extra
  ))
  direct <- candidates[file.exists(candidates)]
  if (length(direct)) {
    return(normalizePath(direct[[1]], winslash = "/", mustWork = TRUE))
  }

  if (dir.exists(ONSEN_DATA_ROOT)) {
    hits <- list.files(
      ONSEN_DATA_ROOT,
      pattern = paste0("^", escape_regex(filename), "$"),
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = FALSE
    )
    if (length(hits)) {
      hits <- hits[order(nchar(hits), hits)]
      return(normalizePath(hits[[1]], winslash = "/", mustWork = TRUE))
    }
  }

  if (required) {
    stop(
      "Required input not found: ", filename,
      "\nSet ONSEN_DATA_ROOT to the project root or place the file in the repository root.",
      call. = FALSE
    )
  }
  NA_character_
}

find_any_input <- function(filenames, required = TRUE) {
  for (f in filenames) {
    hit <- find_input(f, required = FALSE)
    if (!is.na(hit)) return(hit)
  }
  if (required) {
    stop("None of the required input alternatives were found:\n", paste(filenames, collapse = "\n"), call. = FALSE)
  }
  NA_character_
}

out_file <- function(filename, subdir = NULL) {
  d <- if (is.null(subdir)) ONSEN_OUTPUT_ROOT else file.path(ONSEN_OUTPUT_ROOT, subdir)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, filename)
}

message_config <- function() {
  message("Repository root: ", REPO_ROOT)
  message("Data root:       ", ONSEN_DATA_ROOT)
  message("Output root:     ", ONSEN_OUTPUT_ROOT)
  message("Force rescan:    ", ONSEN_FORCE_RESCAN)
  message("Large steps:     ", ONSEN_RUN_LARGE_STEPS)
}
