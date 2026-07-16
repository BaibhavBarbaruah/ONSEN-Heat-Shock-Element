# Master runner for the final Biology Open revision workflow.

source("ONSEN_config.R")
message_config()

steps <- c(
  "01_native_mutated_motif_analysis.R",
  "02_constrained_mutant_sensitivity.R",
  "03_col0_HSF_and_TE_background.R",
  "03B_threshold_and_continuous_sensitivity.R",
  "04_nonredundant_HSF_locations.R",
  "05_methylation_analysis.R",
  "06_rnaseq_analysis.R",
  "07_accession_analysis.R",
  "07B_figure7_identity_and_logo_workflow.R",
  "08_write_supplementary_tables.R",
  "09_write_session_info.R",
  "10_validate_manuscript_outputs.R"
)

message("Run order:")
for (i in seq_along(steps)) message(i, ". ", steps[[i]])

if (!ONSEN_RUN_LARGE_STEPS) {
  message(
    "\nONSEN_RUN_LARGE_STEPS is FALSE.\n",
    "Scripts will prefer deposited/processed outputs where supported.\n",
    "Set ONSEN_RUN_LARGE_STEPS='true' to permit genome-wide scans and other large raw-data steps."
  )
}

for (script in steps) {
  message("\n============================================================")
  message("Running: ", script)
  message("============================================================")
  source(repo_file(script), echo = FALSE)
}
