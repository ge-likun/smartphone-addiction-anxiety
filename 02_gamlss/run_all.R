if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.UTF-8"))
}
script_directory <- function() {
  sourced_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (length(sourced_file) && file.exists(sourced_file)) {
    return(dirname(sourced_file))
  }
  getwd()
}

analysis_directory <- script_directory()
setwd(analysis_directory)

source("01_setup_and_data.R", encoding = "UTF-8")
source("00_zanbi_order_safe.R", encoding = "UTF-8")
source("02_fit_models.R", encoding = "UTF-8")
source("03_residual_diagnostics.R", encoding = "UTF-8")
source("04_figures.R", encoding = "UTF-8")
source("05_tables.R", encoding = "UTF-8")
source("06_covariate_sensitivity.R", encoding = "UTF-8")

writeLines(capture.output(sessionInfo()), file.path(output_dirs$logs, "sessionInfo.txt"))
message("Analysis complete. Results are in: ", normalizePath(output_dirs$root))
