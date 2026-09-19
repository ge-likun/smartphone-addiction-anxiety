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
source("05_combine_binary_runs.R", encoding = "UTF-8")
