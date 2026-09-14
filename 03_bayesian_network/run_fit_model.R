arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments, value = TRUE)
frame_files <- vapply(
  sys.frames(),
  function(frame) if (is.null(frame$ofile)) NA_character_ else frame$ofile,
  character(1)
)
script_file <- if (length(file_argument)) {
  sub("^--file=", "", file_argument[1])
} else if (any(!is.na(frame_files))) {
  tail(frame_files[!is.na(frame_files)], 1)
} else {
  NA_character_
}
project_root <- if (is.na(script_file)) "." else dirname(script_file)

source(file.path(project_root, "00_packages.R"), encoding = "UTF-8")
source(file.path(project_root, "01_config.R"), encoding = "UTF-8")
source(file.path(project_root, "02_data.R"), encoding = "UTF-8")
source(file.path(project_root, "03_model.R"), encoding = "UTF-8")

if (file.exists(model_file)) {
  stop("The model output already exists: ", model_file, call. = FALSE)
}

network_data <- read_network_data(data_file)
fit_network_model(network_data, output_file = model_file)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_model_fit.txt"))
