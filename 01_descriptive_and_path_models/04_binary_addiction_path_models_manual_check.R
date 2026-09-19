analysis_file <- Sys.getenv(
  "ANALYSIS_DATA_FILE",
  unset = file.path("..", "analysis_data.csv")
)
model_output_dir <- Sys.getenv(
  "BINARY_PATH_OUTPUT_DIR",
  unset = file.path("..", "results", "01_binary_addiction_path_models")
)

if (!file.exists(analysis_file)) stop("Analysis data not found: ", analysis_file)
dat <- read.csv(analysis_file, check.names = FALSE)
required <- c(
  "mpiq_total", "addiction_group", "gad7_total", "survey_year", "sex",
  "kinship_support", "offline_friend_support", "online_friend_support"
)
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
dat <- dat[stats::complete.cases(dat[, required]), , drop = FALSE]
if (nrow(dat) == 0L) stop("No complete observations remain.")

expected_group <- ifelse(dat$mpiq_total > 32, "Addicted", "Non-addicted")
if (any(as.character(dat$addiction_group) != expected_group)) {
  stop("addiction_group does not match the MPIQ >32 rule.")
}
if (any(!dat$sex %in% c("Female", "Male")) ||
    any(!dat$survey_year %in% c(2023, 2024))) {
  stop("Unexpected sex or survey_year value.")
}

dat$x <- as.numeric(dat$addiction_group == "Addicted")
dat$year2024 <- as.numeric(dat$survey_year == 2024)
supports <- c(
  kinship_support = "Kinship support",
  offline_friend_support = "Offline-friend support",
  online_friend_support = "Online-friend support"
)

rows <- list()
row_number <- 0L
for (variable in names(supports)) {
  dat$m <- dat[[variable]] - mean(dat[[variable]])
  dat$xm <- dat$x * dat$m
  for (sex_value in c("Female", "Male")) {
    subset <- dat[dat$sex == sex_value, , drop = FALSE]
    support_fit <- stats::lm(m ~ x + year2024, data = subset)
    gad7_fit <- stats::lm(
      gad7_total ~ x + m + xm + year2024,
      data = subset
    )
    coefficients <- c(
      a = unname(stats::coef(support_fit)["x"]),
      b1 = unname(stats::coef(gad7_fit)["m"]),
      b2 = unname(stats::coef(gad7_fit)["xm"]),
      cp = unname(stats::coef(gad7_fit)["x"])
    )
    for (path in names(coefficients)) {
      row_number <- row_number + 1L
      rows[[row_number]] <- data.frame(
        support = unname(supports[variable]),
        sex = sex_value,
        path = path,
        lm_estimate = unname(coefficients[path]),
        n = nrow(subset)
      )
    }
  }
}

check <- do.call(rbind, rows)
model_file <- file.path(
  model_output_dir,
  "binary_path_coefficient_summary.csv"
)
if (file.exists(model_file)) {
  published <- read.csv(model_file, check.names = FALSE)
  published$path <- sub("_[fm]$", "", published$path)
  match_key <- function(data) {
    paste(data$support, data$sex, data$path, sep = "|")
  }
  index <- match(match_key(check), match_key(published))
  if (anyNA(index)) stop("The fitted-model summary is missing path rows.")
  check$lavaan_estimate <- published$unstandardized_estimate[index]
  check$absolute_difference <- abs(
    check$lm_estimate - check$lavaan_estimate
  )
  check$agrees_with_lavaan <- check$absolute_difference < 1e-6
  if (!all(check$agrees_with_lavaan)) {
    stop("At least one point estimate differs between lm and lavaan.")
  }
}

print(check, row.names = FALSE)
if (dir.exists(model_output_dir)) {
  write.csv(
    check,
    file.path(model_output_dir, "binary_point_estimate_check.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}
