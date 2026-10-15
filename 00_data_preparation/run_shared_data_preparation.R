root <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("..", winslash = "/"))
source(file.path(root, "R", "project.R"), encoding = "UTF-8")
cfg <- read_project_config(root, Sys.getenv("PROJECT_CONFIG_FILE", "local_config.R"))
read_shared <- function(path) {
  path <- absolute_path(path, root)
  if (!file.exists(path)) stop("Input file not found: ", path)
  x <- utils::read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM",
                       na.strings = c("", "NA"), stringsAsFactors = FALSE)
  if (!nrow(x) || anyDuplicated(names(x))) stop("Empty data or duplicate column names: ", path)
  x
}
require_columns <- function(x, columns) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
}
same_values <- function(x, y) {
  identical(is.na(x), is.na(y)) && all(x[!is.na(x)] == y[!is.na(y)])
}
check_values <- function(x, legal, name) {
  if (any(!is.na(x) & !x %in% legal)) stop("Invalid values in ", name)
}
analysis <- read_shared(cfg$shared$cross_sectional)
network <- read_shared(cfg$shared$longitudinal)
support <- c("kinship_support", "offline_friend_support", "online_friend_support")
gad <- paste0("gad", 1:7, "_item")
mpiq <- paste0("mpiq", 1:8, "_item")
support_items <- paste0("support", 1:9, "_item")
demographic <- c("sex_code", "age_years", "bmi_kg_m2", "grade_code", "major_code",
                 "hukou_code", "father_education_code", "mother_education_code",
                 "only_child_code", "ipaq_category_code")
require_columns(analysis, c("record_id", "survey_year", demographic, gad, mpiq,
  support_items, support, "gad7_total", "mpiq_total", "social_support_total",
  "retest_status"))
require_columns(network, c("baseline_record_id", "baseline_year", "followup_year",
  as.vector(outer(c(support, gad, mpiq), c("t1", "t2"), paste, sep = "_")),
  "retest_interval_months", paste0("retest_", support, "_t2")))
for (x in list(analysis$record_id, network$baseline_record_id)) {
  if (anyNA(x) || anyDuplicated(x)) stop("Record identifiers must be present and unique.")
}
baseline_index <- match(network$baseline_record_id, analysis$record_id)
if (anyNA(baseline_index)) stop("Longitudinal baseline records are missing from the cross-sectional file.")
if (!same_values(network$baseline_year, analysis$survey_year[baseline_index])) {
  stop("Baseline year disagrees between files.")
}
if (any(network$followup_year < network$baseline_year, na.rm = TRUE) ||
    any(network$retest_interval_months <= 0, na.rm = TRUE)) stop("Invalid follow-up timing.")
for (name in gad) check_values(analysis[[name]], 0:3, name)
for (name in mpiq) check_values(analysis[[name]], 1:7, name)
for (name in support_items) check_values(analysis[[name]], 0:5, name)
check_values(analysis$sex_code, 1:2, "sex_code")
if (!same_values(analysis$gad7_total, rowSums(analysis[gad])) ||
    !same_values(analysis$mpiq_total, rowSums(analysis[mpiq]))) {
  stop("Questionnaire totals disagree with item coding.")
}
derived <- list(sex = c("Male", "Female")[analysis$sex_code],
  addiction_group_code = as.integer(analysis$mpiq_total > 32),
  addiction_group = ifelse(analysis$mpiq_total > 32, "Addicted", "Non-addicted"),
  gad7_symptomatic_code = as.integer(analysis$gad7_total >= 5),
  grade_manuscript = ifelse(analysis$grade_code >= 5, "Postgraduate",
    ifelse(analysis$grade_code %in% 1:4, paste0("Year ", analysis$grade_code), NA_character_)))
for (name in names(derived)) {
  if (name %in% names(analysis) && !same_values(analysis[[name]], derived[[name]])) {
    stop("Derived field disagrees with its source variables: ", name)
  }
  analysis[[name]] <- derived[[name]]
}
for (name in c("sex", "addiction_group_t1")) {
  expected <- analysis[[if (name == "sex") "sex" else "addiction_group"]][baseline_index]
  if (name %in% names(network) && !same_values(network[[name]], expected)) stop("Baseline grouping disagrees between files.")
  network[[name]] <- expected
}
for (i in seq_along(support)) {
  columns <- support_items[((i - 1L) * 3L + 1L):(i * 3L)]
  if (!same_values(analysis[[support[i]]], rowSums(analysis[columns]))) stop("Support score disagrees with item sum.")
}
if (!same_values(analysis$social_support_total, rowSums(analysis[support]))) stop("Total support disagrees with its components.")
for (name in c(support, gad, mpiq)) {
  if (!same_values(network[[paste0(name, "_t1")]], analysis[[name]][baseline_index])) {
    stop("Baseline values disagree between files: ", name)
  }
  legal <- if (name %in% support) 0:15 else if (name %in% gad) 0:3 else 1:7
  check_values(network[[paste0(name, "_t2")]], legal, paste0(name, "_t2"))
}
for (name in paste0("retest_", support, "_t2")) check_values(network[[name]], 0:15, name)
status <- ifelse(analysis$record_id %in% network$baseline_record_id, "Retested", "Not retested")
if (!same_values(analysis$retest_status, status)) stop("Retest status disagrees with linked records.")
analysis$participant_id <- analysis$record_id
analysis$record_id <- NULL
analysis$retest_pair_id <- ifelse(analysis$participant_id %in% network$baseline_record_id,
                                 analysis$participant_id, NA_character_)
network$retest_pair_id <- network$baseline_record_id
network$baseline_participant_id <- network$baseline_record_id
network$mpiq_total_t1 <- rowSums(network[paste0(mpiq, "_t1")])
network$mpiq_total_t2 <- rowSums(network[paste0(mpiq, "_t2")])
retest <- network[, c("retest_pair_id", "baseline_participant_id",
  "baseline_year", "followup_year", "sex", "addiction_group_t1", paste0(support, "_t1")), drop = FALSE]
for (name in support) retest[[paste0(name, "_t2")]] <- network[[paste0("retest_", name, "_t2")]]
retest$social_support_total_t1 <- rowSums(retest[paste0(support, "_t1")])
retest$social_support_total_t2 <- rowSums(retest[paste0(support, "_t2")])
retest$interval_months <- network$retest_interval_months
output <- Sys.getenv("PROJECT_DATA_DIR", absolute_path(cfg$paths$data, root))
report <- file.path(Sys.getenv("PROJECT_RESULTS_DIR", absolute_path(cfg$paths$results, root)), "00_data_preparation")
for (d in c(output, report)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
write_csv(analysis, file.path(output, "analysis_data.csv"))
write_csv(network, file.path(output, "network_data.csv"))
write_csv(retest, file.path(output, "retest_data.csv"))
write_csv(data.frame(dataset = c("analysis", "retest", "network"),
  n = c(nrow(analysis), nrow(retest), nrow(network))), file.path(report, "sample_sizes.csv"))
sources <- c(cfg$shared$cross_sectional, cfg$shared$longitudinal)
write_csv(data.frame(file = sources, md5 = unname(tools::md5sum(vapply(sources,
  absolute_path, character(1), root = root)))), file.path(report, "input_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(report, "sessionInfo.txt"))
message("Prepared analysis, retest and network inputs from controlled-access data.")
