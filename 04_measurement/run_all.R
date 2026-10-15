if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.UTF-8"))
}
script_directory <- function() {
  sourced_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (length(sourced_file) && file.exists(sourced_file)) {
    return(dirname(normalizePath(sourced_file)))
  }
  normalizePath(getwd())
}

module_root <- script_directory()
setwd(module_root)

source(file.path("R", "01_social_support_cfa.R"), encoding = "UTF-8")
source(file.path("R", "03_measurement_tables.R"), encoding = "UTF-8")
