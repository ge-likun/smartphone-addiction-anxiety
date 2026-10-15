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

source("00_functions.R", encoding = "UTF-8")
source("01_descriptive_tables.R", encoding = "UTF-8")
source("02_moderated_mediation.R", encoding = "UTF-8")
source("03_retest_stability.R", encoding = "UTF-8")
