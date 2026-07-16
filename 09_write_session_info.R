# Capture the exact R and package environment used for reproduction.

source("ONSEN_functions.R")
write_session_info("ONSEN_HSE_repository_sessionInfo.txt")

installed <- as.data.frame(installed.packages()[, c("Package", "Version", "LibPath")])
safe_write_csv(installed, "ONSEN_HSE_installed_packages.csv")

message("Session information written to: ", ONSEN_OUTPUT_ROOT)
