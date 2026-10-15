if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.UTF-8"))
}
args <- commandArgs(FALSE)
file <- grep("^--file=", args, value = TRUE)
root <- if (length(file)) dirname(normalizePath(sub("^--file=", "", file[1]), winslash = "/")) else getwd()
files <- list.files(root, pattern = "\\.R$|\\.R\\.example$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("[/\\\\](results|data_raw|data)[/\\\\]", files)]
for (file in files) parse(file, encoding = "UTF-8")
cat("Parsed", length(files), "R files.\n")
artifacts <- list.files(root, pattern = "\\.(rds|RData|rda|xlsx|csv|docx|png|pdf|tiff)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
artifacts <- artifacts[!grepl("[/\\\\](results|data_raw|data)[/\\\\]", artifacts)]
if (length(artifacts)) stop("Unexpected generated artifacts in source directories: ", paste(artifacts, collapse = ", "))
cat("No participant files or generated artifacts in source directories.\n")
