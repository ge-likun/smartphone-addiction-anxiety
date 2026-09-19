required_packages(c("flextable", "officer"))

input_dir <- Sys.getenv("BINARY_PARALLEL_BASE_DIR", unset = "")
output_dir <- Sys.getenv("BINARY_PATH_OUTPUT_DIR", unset = "")
if (!nzchar(input_dir) || !nzchar(output_dir)) {
  stop("Set BINARY_PARALLEL_BASE_DIR and BINARY_PATH_OUTPUT_DIR.")
}

sources <- c(
  kinship = "kinship_support",
  offline = "offline_friend_support",
  online = "online_friend_support"
)
read_result <- function(folder, filename) {
  path <- file.path(input_dir, folder, filename)
  if (!file.exists(path)) stop("Missing result: ", path)
  read.csv(path, check.names = FALSE)
}

status <- do.call(rbind, lapply(names(sources), function(folder) {
  read_result(folder, "run_status.csv")
}))
if (nrow(status) != 3L ||
    any(!status$converged) ||
    any(status$failed_bootstraps != 0L) ||
    length(unique(status$n)) != 1L ||
    length(unique(status$bootstrap_replicates)) != 1L) {
  stop("The three model runs do not pass the completion checks.")
}

path_summary <- do.call(rbind, lapply(names(sources), function(folder) {
  read_result(folder, "binary_path_coefficient_summary.csv")
}))
defined_summary <- do.call(rbind, lapply(names(sources), function(folder) {
  read_result(folder, "binary_defined_parameter_summary.csv")
}))
publication <- do.call(rbind, lapply(names(sources), function(folder) {
  read_result(folder, "binary_sensitivity_publication_source_data.csv")
}))
for (object_name in c("path_summary", "defined_summary")) {
  object <- get(object_name)
  names(object)[names(object) == "two_sided_p"] <- "normal_theory_p"
  assign(object_name, object)
}
if (nrow(path_summary) != 24L || nrow(defined_summary) != 48L) {
  stop("The combined coefficient files do not contain all three support sources.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(status, file.path(output_dir, "run_status.csv"), row.names = FALSE)
write.csv(
  path_summary,
  file.path(output_dir, "binary_path_coefficient_summary.csv"),
  row.names = FALSE
)
write.csv(
  defined_summary,
  file.path(output_dir, "binary_defined_parameter_summary.csv"),
  row.names = FALSE
)
write.csv(
  publication,
  file.path(output_dir, "binary_sensitivity_publication_source_data.csv"),
  row.names = FALSE
)
for (folder in names(sources)) {
  filename <- paste0(sources[[folder]], "_parameter_estimates.csv")
  from <- file.path(input_dir, folder, filename)
  if (!file.copy(from, file.path(output_dir, filename))) {
    stop("Could not copy model estimates: ", from)
  }
}

publication$row_type <- ifelse(
  publication[["B or product"]] == "", "header", "body"
)
publication <- publication[, c("row_type", setdiff(names(publication), "row_type"))]
save_table_docx(
  publication,
  "Binary addiction-status sensitivity analysis of sex-grouped path models",
  paste(
    "Addiction status was coded 0 for MPIQ scores of 32 or lower and 1 for scores greater than 32.",
    "Support scores were centred at their full-sample means; survey year was included in both equations.",
    "Path b1 is the support coefficient in the non-addicted group; b1 + b2 is the corresponding coefficient in the addicted group.",
    "Path b2 is the addicted-minus-non-addicted difference in the support coefficient.",
    "P values use two-sided normal-approximation tests with bootstrap standard errors; 95% intervals are bias-corrected bootstrap intervals.",
    "Product terms summarize associations measured at the same assessment."
  ),
  file.path(output_dir, "binary_addiction_status_sensitivity_publication.docx")
)

file.copy(
  file.path(input_dir, "kinship", "sessionInfo.txt"),
  file.path(output_dir, "sessionInfo.txt")
)
message("Combined binary sensitivity results written to ", output_dir)
