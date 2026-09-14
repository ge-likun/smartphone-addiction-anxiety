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
residual_replicates <- 50L
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

if (anyNA(analysis_data[required_variables])) stop("Required analysis variables contain missing values")
if (any(!analysis_data$gad7_total %in% 0:maximum_score)) stop("gad7_total must be integer-valued from 0 to 21")
for (variable in c("kinship_support", "offline_friend_support", "online_friend_support")) {
  if (any(analysis_data[[variable]] < 0 | analysis_data[[variable]] > 15)) {
    stop(variable, " must range from 0 to 15")
  }
}
if (anyNA(analysis_data$sex) || anyNA(analysis_data$addiction_group)) {
  stop("sex and addiction_group contain values outside the documented categories")
}
if (!all(as.character(analysis_data$survey_year) %in% c("2023", "2024"))) {
  stop("survey_year must contain 2023 and/or 2024")
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
model_cache_tag <- substr(data_fingerprint, 1L, 12L)

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
  if (nrow(stratum_data) < 20L) {
    stop(spec$stratum, " must contain at least 20 observations")
  }
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
  grDevices::cairo_pdf(paste0(basename, ".pdf"), width = width_in, height = height_in, family = "Arial")
  print(plot)
  grDevices::dev.off()
  svglite::svglite(paste0(basename, ".svg"), width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()
  write_raster(ragg::agg_tiff, ".tiff", background = "white", compression = "lzw")
}
