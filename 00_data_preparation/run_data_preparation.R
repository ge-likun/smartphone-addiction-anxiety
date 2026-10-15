if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.UTF-8"))
}
root <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("..", winslash = "/"))
source(file.path(root, "R", "project.R"), encoding = "UTF-8")
cfg <- read_project_config(root, Sys.getenv("PROJECT_CONFIG_FILE", "local_config.R"))
raw_config <- cfg$raw
read_input <- function(file, sheet = 1) {
  file <- absolute_path(file, root)
  if (!file.exists(file)) stop("Input file not found: ", file)
  if (tolower(tools::file_ext(file)) == "csv") {
    x <- utils::read.csv(file, check.names = FALSE, colClasses = "character", fileEncoding = "UTF-8")
  } else {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl to read Excel files.")
    x <- as.data.frame(readxl::read_excel(file, sheet = sheet, col_types = "text",
                                        .name_repair = "minimal"), check.names = FALSE)
  }
  mapped <- unique(c(raw_config$id_column, raw_config$timestamp_column,
    raw_config$time_id_column, raw_config$time_column, unname(raw_config$columns),
    raw_config$gad_columns, raw_config$mpiq_columns, raw_config$support_columns,
    unname(raw_config$reported_totals)))
  duplicates <- unique(names(x)[duplicated(names(x))])
  if (any(duplicates %in% mapped)) stop("Ambiguous mapped input column: ",
    paste(intersect(duplicates, mapped), collapse = ", "))
  names(x) <- make.unique(names(x))
  x
}
normalise_id <- function(x) {
  x <- trimws(as.character(x))
  scientific <- !is.na(x) & grepl("^[+-]?[0-9]+(\\.[0-9]+)?[eE][+-]?[0-9]+$", x)
  if (any(scientific)) {
    numeric_id <- as.numeric(x[scientific])
    safe <- is.finite(numeric_id) & abs(numeric_id) <= 2^53 & numeric_id == trunc(numeric_id)
    if (any(!safe)) message(sum(!safe), " scientific-notation IDs retained unchanged; exact integer conversion unavailable.")
    x[which(scientific)[safe]] <- sprintf("%.0f", numeric_id[safe])
  }
  x <- sub("\\.0+$", "", x)
  x[!nzchar(x)] <- NA_character_
  x
}
numeric_column <- function(x, name) {
  out <- suppressWarnings(as.numeric(x))
  invalid <- !is.na(x) & nzchar(trimws(as.character(x))) & is.na(out)
  if (any(invalid)) stop("Non-numeric values in ", name, ": ", sum(invalid), " rows.")
  out
}
parse_time <- function(x) {
  as.POSIXct(x, tz = raw_config$timezone,
             tryFormats = c("%Y/%m/%d %H:%M:%S", "%Y-%m-%d %H:%M:%S", "%Y/%m/%d", "%Y-%m-%d"))
}
read_wave <- function(file, sheet, year, followup = FALSE) {
  x <- read_input(file, sheet)
  item_columns <- c(raw_config$gad_columns, raw_config$mpiq_columns, raw_config$support_columns)
  required <- c(raw_config$id_column, item_columns, if (followup) raw_config$timestamp_column else raw_config$columns[["sex_code"]])
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(basename(file), " is missing: ", paste(missing, collapse = ", "))
  d <- data.frame(raw_id = normalise_id(x[[raw_config$id_column]]), survey_year = year,
                  source_row = seq_len(nrow(x)), stringsAsFactors = FALSE)
  for (name in names(raw_config$columns)) {
    col <- raw_config$columns[[name]]
    d[[name]] <- if (!followup && col %in% names(x)) numeric_column(x[[col]], col) else NA_real_
  }
  groups <- list(gad = raw_config$gad_columns, mpiq = raw_config$mpiq_columns, support = raw_config$support_columns)
  for (domain in names(groups)) for (i in seq_along(groups[[domain]])) {
    col <- groups[[domain]][i]
    value <- numeric_column(x[[col]], col)
    value[value %in% raw_config$item_missing_codes] <- NA_real_
    if (domain == "gad") value <- value - raw_config$gad_input_offset
    if (domain == "support") value <- value - raw_config$support_input_offset
    legal <- switch(domain, gad = 0:3, mpiq = 1:7, support = 0:5)
    if (any(!value %in% legal & !is.na(value))) stop("Invalid questionnaire responses in ", col)
    d[[paste0(domain, i, "_item")]] <- value
  }
  d$gad7_total <- rowSums(d[paste0("gad", seq_along(groups$gad), "_item")])
  d$mpiq_total <- rowSums(d[paste0("mpiq", seq_along(groups$mpiq), "_item")])
  support_names <- c("kinship_support", "offline_friend_support", "online_friend_support")
  for (i in seq_along(support_names)) d[[support_names[i]]] <- rowSums(d[paste0("support", ((i - 1L) * 3L + 1L):(i * 3L), "_item")])
  totals_map <- if (followup && !is.null(raw_config$followup_reported_totals))
    raw_config$followup_reported_totals else raw_config$reported_totals
  for (name in names(totals_map)) {
    col <- totals_map[[name]]
    if (col %in% names(x)) {
      supplied <- numeric_column(x[[col]], col)
      both <- !is.na(supplied) & !is.na(d[[name]])
      if (any(supplied[both] != d[[name]][both])) stop("Reported score disagrees with item sum: ", col)
    }
  }
  d$social_support_total <- rowSums(d[support_names])
  d$addiction_group_code <- as.integer(d$mpiq_total > 32)
  d$addiction_group <- factor(d$addiction_group_code, 0:1, c("Non-addicted", "Addicted"))
  d$gad7_symptomatic_code <- as.integer(d$gad7_total >= 5)
  d$sex <- factor(d$sex_code, c(1, 2), c("Male", "Female"))
  if (!followup && any(!d$sex_code %in% c(1, 2) & !is.na(d$sex_code))) stop("Sex codes must be 1 or 2.")
  d$grade_manuscript <- ifelse(d$grade_code >= 5, "Postgraduate", ifelse(d$grade_code %in% 1:4, paste0("Year ", d$grade_code), NA_character_))
  if (followup) d$timestamp <- parse_time(x[[raw_config$timestamp_column]])
  d
}
surveys <- cfg$raw$surveys
waves <- lapply(seq_len(nrow(surveys)), function(i) read_wave(surveys$file[i], surveys$sheet[i], surveys$year[i]))
pooled <- do.call(rbind, waves)
baseline <- pooled[pooled$survey_year == raw_config$baseline_year & !is.na(pooled$raw_id), , drop = FALSE]
if (anyDuplicated(baseline$raw_id)) stop("Duplicate identifiers in the baseline survey.")
files <- cfg$raw$followups
followups <- do.call(rbind, lapply(seq_len(nrow(files)), function(i) {
  d <- read_wave(files$file[i], files$sheet[i], raw_config$followup_year, TRUE)
  d$batch <- i
  d
}))
eligible <- !is.na(followups$raw_id) & !is.na(followups$timestamp) &
  nchar(followups$raw_id) >= raw_config$minimum_id_length & followups$raw_id %in% baseline$raw_id
followups <- followups[eligible, , drop = FALSE]
if (!nrow(followups)) stop("No eligible baseline-follow-up links.")
pick_followup <- function(last = FALSE, records = followups) {
  ordered <- if (last) order(records$raw_id, -as.numeric(records$timestamp), records$batch, records$source_row) else
    order(records$raw_id, records$timestamp, records$batch, records$source_row)
  d <- records[ordered, , drop = FALSE]
  d[!duplicated(d$raw_id), , drop = FALSE]
}
earliest <- pick_followup()
batch_rule <- raw_config$network_batch_rule
if (is.null(batch_rule)) batch_rule <- "earliest"
if (!batch_rule %in% c("earliest", "all")) stop("network_batch_rule must be earliest or all.")
network_records <- followups
if (batch_rule == "earliest") {
  network_records <- network_records[order(network_records$batch, network_records$raw_id,
    network_records$timestamp, network_records$source_row), , drop = FALSE]
  network_records <- network_records[!duplicated(network_records[c("batch", "raw_id")]), , drop = FALSE]
}
latest <- pick_followup(TRUE, network_records)
ids <- sort(unique(earliest$raw_id))
pair_ids <- ids
if (!is.null(raw_config$retest_order_seed)) {
  set.seed(raw_config$retest_order_seed)
  pair_ids <- sample(ids)
}
pair_map <- data.frame(raw_id = pair_ids, retest_pair_id = sprintf("R%05d", seq_along(pair_ids)))
if (!is.null(raw_config$analysis_order_seed)) {
  set.seed(raw_config$analysis_order_seed)
  pooled <- pooled[sample.int(nrow(pooled)), , drop = FALSE]
}
pooled$participant_id <- sprintf("P%06d", seq_len(nrow(pooled)))
pooled$retest_pair_id <- pair_map$retest_pair_id[match(pooled$raw_id, pair_map$raw_id)]
pooled$retest_pair_id[pooled$survey_year != raw_config$baseline_year] <- NA_character_
pooled$retest_status <- ifelse(is.na(pooled$retest_pair_id), "Not retested", "Retested")
baseline <- pooled[pooled$survey_year == raw_config$baseline_year & pooled$raw_id %in% ids, , drop = FALSE]
time_spec <- cfg$raw$baseline_times
times <- if (nrow(time_spec)) do.call(rbind, lapply(seq_len(nrow(time_spec)), function(i) {
  x <- read_input(time_spec$file[i], time_spec$sheet[i])
  required <- c(raw_config$time_id_column, raw_config$time_column)
  if (!all(required %in% names(x))) stop("Missing baseline timestamp columns.")
  data.frame(raw_id = normalise_id(x[[required[1]]]), timestamp = parse_time(x[[required[2]]]))
})) else data.frame(raw_id = character(), timestamp = as.POSIXct(character()))
times <- times[order(times$raw_id, times$timestamp), , drop = FALSE]
times <- times[!duplicated(times$raw_id), , drop = FALSE]
core <- c(paste0("mpiq", 1:8, "_item"), "mpiq_total", "addiction_group_code",
          paste0("gad", 1:7, "_item"), "gad7_total", "gad7_symptomatic_code",
          paste0("support", 1:9, "_item"), "kinship_support", "offline_friend_support", "online_friend_support", "social_support_total")
build_pairs <- function(followup) {
  b <- baseline[match(ids, baseline$raw_id), , drop = FALSE]
  f <- followup[match(ids, followup$raw_id), , drop = FALSE]
  out <- data.frame(retest_pair_id = pair_map$retest_pair_id[match(ids, pair_map$raw_id)], baseline_participant_id = b$participant_id,
                     baseline_year = raw_config$baseline_year, followup_year = raw_config$followup_year,
                     followup_batch = f$batch, sex_code = b$sex_code, sex = b$sex,
                     addiction_group_t1 = b$addiction_group, addiction_group_t2 = f$addiction_group)
  for (name in core) {
    out[[paste0(name, "_t1")]] <- b[[name]]
    out[[paste0(name, "_t2")]] <- f[[name]]
  }
  for (name in names(raw_config$columns)) if (!name %in% names(out)) out[[name]] <- b[[name]]
  out$grade_manuscript <- b$grade_manuscript
  timestamp_t1 <- times$timestamp[match(ids, times$raw_id)]
  out$interval_days <- as.numeric(difftime(f$timestamp, timestamp_t1, units = "days"))
  out$interval_months <- out$interval_days / 30.4375
  if (any(out$interval_days <= 0, na.rm = TRUE)) stop("Follow-up must occur after baseline.")
  out
}
retest <- build_pairs(earliest)
retest <- retest[order(retest$retest_pair_id), , drop = FALSE]
network <- build_pairs(latest)
analysis <- pooled[, setdiff(names(pooled), c("raw_id", "source_row")), drop = FALSE]
output <- Sys.getenv("PROJECT_DATA_DIR", absolute_path(cfg$paths$data, root))
report <- file.path(Sys.getenv("PROJECT_RESULTS_DIR", absolute_path(cfg$paths$results, root)), "00_data_preparation")
for (d in c(output, report)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
write_csv(analysis, file.path(output, "analysis_data.csv"))
write_csv(retest, file.path(output, "retest_data.csv"))
write_csv(network, file.path(output, "network_data.csv"))
write_csv(data.frame(dataset = c("analysis", "retest", "network"), n = c(nrow(analysis), nrow(retest), nrow(network))), file.path(report, "sample_sizes.csv"))
sources <- unique(c(surveys$file, files$file, time_spec$file))
write_csv(data.frame(file = sources, md5 = unname(tools::md5sum(vapply(sources, absolute_path, character(1), root = root)))), file.path(report, "input_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(report, "sessionInfo.txt"))
message("Prepared analysis, retest and network datasets. Network within-batch rule: ", batch_rule,
        "; cross-batch rule: latest. Retest selects earliest eligible follow-up.")
