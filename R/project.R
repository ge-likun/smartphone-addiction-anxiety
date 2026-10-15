project_root_from_script <- function() {
  frames <- lapply(sys.frames(), function(x) x$ofile)
  files <- Filter(function(x) length(x) == 1L && file.exists(x), frames)
  if (length(files)) return(dirname(normalizePath(files[[1]], winslash = "/")))
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) return(dirname(normalizePath(sub("^--file=", "", arg[1]), winslash = "/")))
  normalizePath(getwd(), winslash = "/")
}
absolute_path <- function(path, root) {
  if (grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", path)) path else file.path(root, path)
}
read_project_config <- function(root, file = "local_config.R") {
  path <- absolute_path(file, root)
  if (!file.exists(path)) {
    if (basename(path) != "local_config.R") stop("Configuration file not found: ", path)
    path <- file.path(root, "config.R.example")
  }
  e <- new.env(parent = baseenv())
  sys.source(path, e, keep.source = FALSE)
  if (!is.list(e$config)) stop("The configuration must define a list named config.")
  e$config
}
object_hash <- function(x) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(x, path, version = 3, compress = FALSE)
  unname(tools::md5sum(path))
}
model_data_signature <- function(data, columns) {
  setNames(lapply(columns, function(column) {
    x <- data[[column]]
    if (is.factor(x)) list(values = as.character(x), levels = levels(x), ordered = is.ordered(x))
    else if (is.numeric(x)) as.numeric(x) else x
  }), columns)
}
software_versions <- function(packages = c("brms", "cmdstanr", "posterior")) {
  x <- setNames(lapply(packages, function(p) as.character(utils::packageVersion(p))), packages)
  x$R <- paste(R.version$major, R.version$minor, sep = ".")
  if ("cmdstanr" %in% packages) x$CmdStan <- as.character(cmdstanr::cmdstan_version())
  x
}
write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}
run_module <- function(root, directory, runner, arguments = character()) {
  old <- getwd()
  on.exit(setwd(old))
  setwd(file.path(root, directory))
  executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  status <- system2(executable, c(shQuote(runner), vapply(arguments, shQuote, character(1))))
  if (status != 0L) stop("Stage failed: ", directory, "/", runner, " (exit ", status, ").")
  invisible(status)
}
