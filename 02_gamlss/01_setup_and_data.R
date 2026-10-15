required_packages <- c(
  "gamlss", "gamlss.dist", "dplyr", "tibble", "ggplot2", "patchwork",
  "officer", "flextable", "writexl", "svglite", "ragg"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install the following R packages before running the analysis: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(gamlss)
  library(gamlss.dist)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(officer)
  library(flextable)
})

options(stringsAsFactors = FALSE, width = 160)

data_file <- Sys.getenv("GAMLSS_DATA_FILE", unset = file.path("..", "analysis_data.csv"))
output_root <- Sys.getenv("GAMLSS_OUTPUT_DIR", unset = file.path("..", "results", "02_gamlss"))
residual_seed_base <- 20260803L
residual_replicates <- as.integer(Sys.getenv("GAMLSS_RESIDUAL_REPLICATES", "50"))
if (!is.finite(residual_replicates) || residual_replicates < 1L)
  stop("GAMLSS_RESIDUAL_REPLICATES must be a positive integer.")
maximum_score <- 21L

output_dirs <- list(
  root = output_root,
  models = file.path(output_root, "models"),
  figures = file.path(output_root, "figures"),
  tables = file.path(output_root, "tables"),
  source_data = file.path(output_root, "source_data"),
  logs = file.path(output_root, "logs")
)
for (directory in unique(unlist(output_dirs, use.names = FALSE))) {
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create output directory: ", directory)
  }
}

if (!file.exists(data_file)) {
  stop(
    "Analysis data not found: ", data_file, "\n",
    "Participant-level data are intentionally excluded from the public repository.\n",
    "Set GAMLSS_DATA_FILE before running if the CSV is stored elsewhere."
  )
}
analysis_data <- read.csv(data_file, check.names = FALSE, fileEncoding = "UTF-8")

support_items <- paste0("support", 1:9, "_item")
support_scores <- c("kinship_support", "offline_friend_support", "online_friend_support")
if (any(support_items %in% names(analysis_data))) {
  if (!all(c(support_items, support_scores) %in% names(analysis_data))) {
    stop("Support item and subscale columns are incomplete")
  }
  if (!all(vapply(analysis_data[support_items], is.numeric, logical(1)))) {
    stop("Support item columns must be numeric")
  }
  for (i in seq_along(support_scores)) {
    items <- support_items[((i - 1L) * 3L + 1L):(i * 3L)]
    complete <- stats::complete.cases(analysis_data[items])
    observed <- analysis_data[[support_scores[i]]]
    inconsistent <- (!complete & !is.na(observed)) |
      (complete & (is.na(observed) |
                   observed != rowSums(analysis_data[items], na.rm = TRUE)))
    if (any(inconsistent)) {
      stop(support_scores[i], " does not match the available item responses")
    }
  }
}

required_variables <- c(
  "survey_year", "sex", "addiction_group", "gad7_total",
  "kinship_support", "offline_friend_support", "online_friend_support"
)
missing_variables <- setdiff(required_variables, names(analysis_data))
if (length(missing_variables)) stop("Missing required variables: ", paste(missing_variables, collapse = ", "))

analysis_data <- analysis_data |>
  transmute(
    survey_year = factor(.data$survey_year),
    sex = factor(.data$sex, levels = c("Male", "Female")),
    addiction_group = factor(.data$addiction_group, levels = c("Non-addicted", "Addicted")),
    gad7_total = as.numeric(.data$gad7_total),
    kinship_support = as.numeric(.data$kinship_support),
    offline_friend_support = as.numeric(.data$offline_friend_support),
    online_friend_support = as.numeric(.data$online_friend_support),
    across(any_of(c(
      "grade_code", "major_code", "hukou_code", "father_education_code",
      "mother_education_code", "only_child_code", "ipaq_category_code"
    )), identity)
  )

complete_rows <- stats::complete.cases(analysis_data[required_variables])
if (any(!complete_rows)) {
  message(sum(!complete_rows), " rows with missing GAMLSS variables were excluded")
  analysis_data <- analysis_data[complete_rows, , drop = FALSE]
}
if (nrow(analysis_data) == 0L) stop("No complete observations remain for GAMLSS")
if (any(!analysis_data$gad7_total %in% 0:maximum_score)) stop("gad7_total must be integer-valued from 0 to 21")
for (variable in c("kinship_support", "offline_friend_support", "online_friend_support")) {
  if (any(analysis_data[[variable]] < 0 | analysis_data[[variable]] > 15)) {
    stop(variable, " must range from 0 to 15")
  }
}
if (anyNA(analysis_data$sex) || anyNA(analysis_data$addiction_group)) {
  stop("sex and addiction_group contain values outside the documented categories")
}
if (nlevels(droplevels(analysis_data$survey_year)) < 2L) {
  stop("Both survey years are required for year-adjusted models")
}
for (variable in c("kinship_support", "offline_friend_support", "online_friend_support")) {
  if (any(abs(analysis_data[[variable]] - round(analysis_data[[variable]])) > 1e-12)) {
    stop(variable, " must be integer-valued")
  }
}

stratum_counts <- with(analysis_data, table(addiction_group, sex, useNA = "ifany"))
if (any(stratum_counts == 0L)) {
  stop("All addiction-group-by-sex strata must contain observations")
}

fingerprint_file <- tempfile(pattern = "gamlss_data_", fileext = ".rds")
saveRDS(analysis_data, fingerprint_file, compress = FALSE)
data_fingerprint <- unname(tools::md5sum(fingerprint_file))
unlink(fingerprint_file)
model_cache_tag <- paste0(
  substr(data_fingerprint, 1L, 12L),
  "_g", gsub("\\.", "", as.character(utils::packageVersion("gamlss"))),
  "_d", gsub("\\.", "", as.character(utils::packageVersion("gamlss.dist")))
)

support_variables <- c(
  Kinship = "kinship_support",
  Offline_friend = "offline_friend_support",
  Online_friend = "online_friend_support"
)
support_labels <- c(
  Kinship = "Kinship support",
  Offline_friend = "Offline-friend support",
  Online_friend = "Online-friend support"
)

analysis_specs <- tribble(
  ~model_id, ~stratum, ~addiction_group, ~sex, ~model_set,
  "Overall_Nonaddicted", "Non-addicted overall", "Non-addicted", NA_character_, "overall",
  "Overall_Addicted", "Addicted overall", "Addicted", NA_character_, "overall",
  "Nonaddicted_Male", "Non-addicted men", "Non-addicted", "Male", "sex_stratified",
  "Addicted_Male", "Addicted men", "Addicted", "Male", "sex_stratified",
  "Nonaddicted_Female", "Non-addicted women", "Non-addicted", "Female", "sex_stratified",
  "Addicted_Female", "Addicted women", "Addicted", "Female", "sex_stratified"
)

model_data_for <- function(spec) {
  out <- analysis_data |> filter(.data$addiction_group == .env$spec$addiction_group)
  if (!is.na(spec$sex)) out <- out |> filter(.data$sex == .env$spec$sex)
  droplevels(out)
}

for (spec_index in seq_len(nrow(analysis_specs))) {
  spec <- analysis_specs[spec_index, ]
  stratum_data <- model_data_for(spec)
  if (!nrow(stratum_data)) stop("Empty analysis stratum: ", spec$stratum)
  if (!any(stratum_data$gad7_total == 0) || !any(stratum_data$gad7_total > 0)) {
    stop(spec$stratum, " must contain both zero and positive GAD-7 scores")
  }
  for (variable in c("kinship_support", "offline_friend_support", "online_friend_support")) {
    if (length(unique(stratum_data[[variable]])) < 4L) {
      stop(spec$stratum, " must contain at least four observed values of ", variable)
    }
  }
}

rhs_terms <- function(model_set) {
  terms <- c(
    "pb(kinship_support)", "pb(offline_friend_support)", "pb(online_friend_support)"
  )
  if (identical(model_set, "overall")) c(terms, "sex", "survey_year") else c(terms, "survey_year")
}

formula_text <- function(formula) paste(deparse(formula), collapse = " ")

gamlss_model_data_match <- function(model, model_data, allow_row_permutation = FALSE) {
  response <- if (is.matrix(model$y)) model$y[, 1L] else model$y
  if (length(response) != nrow(model_data)) return(FALSE)
  observed <- data.frame(response = as.numeric(response))
  expected <- data.frame(response = model_data$gad7_total)
  for (parameter in model$parameters) {
    formula <- model[[paste0(parameter, ".formula")]]
    intended <- stats::model.matrix(stats::delete.response(stats::terms(formula)), model_data)
    actual <- model[[paste0(parameter, ".x")]]
    if (is.null(actual) || !setequal(colnames(actual), colnames(intended))) return(FALSE)
    actual <- actual[, colnames(intended), drop = FALSE]
    for (column in seq_len(ncol(intended))) {
      key <- paste(parameter, colnames(intended)[column], sep = "__")
      observed[[key]] <- as.numeric(actual[, column])
      expected[[key]] <- as.numeric(intended[, column])
    }
  }
  canonical <- function(x) {
    if (allow_row_permutation) x <- x[do.call(order, c(x, list(na.last = TRUE))), , drop = FALSE]
    rownames(x) <- NULL
    x
  }
  isTRUE(all.equal(canonical(observed), canonical(expected), tolerance = 1e-12, check.attributes = FALSE))
}

gamlss_cache_identity <- function(model_data, formulas) {
  temporary <- tempfile(pattern = "gamlss_cache_data_", fileext = ".rds")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(model_data, temporary, compress = FALSE)
  list(
    input_fingerprint = data_fingerprint,
    model_data_fingerprint = unname(tools::md5sum(temporary)),
    formulas = vapply(formulas, formula_text, character(1)),
    package_versions = c(gamlss = as.character(packageVersion("gamlss")),
                         gamlss.dist = as.character(packageVersion("gamlss.dist")))
  )
}

verify_gamlss_cache_identity <- function(model, expected, cache_file) {
  if (!identical(attr(model, "cache_identity", exact = TRUE), expected)) {
    stop("Cached model data, formulas or software could not be verified: ", cache_file,
         ". Preserve the file and use a separate output directory for a new fit.")
  }
  invisible(model)
}

save_publication_figure <- function(plot, basename, width_mm = 183, height_mm = 116) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  write_raster <- function(device, extension, ...) {
    target <- paste0(basename, extension)
    temporary <- tempfile(pattern = "gamlss_figure_", fileext = extension)
    device(temporary, width = width_in, height = height_in, units = "in", res = 600, ...)
    print(plot)
    grDevices::dev.off()
    copied <- file.copy(temporary, target, overwrite = TRUE)
    unlink(temporary)
    if (!copied) stop("Unable to write figure: ", target)
  }
  write_raster(ragg::agg_png, ".png", background = "white")
  write_vector <- function(device, extension, ...) {
    target <- paste0(basename, extension)
    temporary <- tempfile(pattern = "gamlss_figure_", fileext = extension)
    device(temporary, ...)
    print(plot)
    grDevices::dev.off()
    copied <- file.copy(temporary, target, overwrite = TRUE)
    unlink(temporary)
    if (!copied) stop("Unable to write figure: ", target)
  }
  write_vector(grDevices::cairo_pdf, ".pdf", width = width_in, height = height_in, family = "Arial")
  write_vector(svglite::svglite, ".svg", width = width_in, height = height_in)
  write_raster(ragg::agg_tiff, ".tiff", background = "white", compression = "lzw")
}
