if (!requireNamespace("lavaan", quietly = TRUE)) {
  stop("Install the lavaan package before running this analysis.", call. = FALSE)
}
suppressPackageStartupMessages(library(lavaan))

data_file <- Sys.getenv(
  "ANALYSIS_DATA_FILE",
  unset = file.path("..", "data", "analysis_data.csv")
)
output_dir <- Sys.getenv(
  "CFA_OUTPUT_DIR",
  unset = file.path("..", "results", "04_measurement", "cfa")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_file)) {
  stop("Analysis data not found: ", data_file)
}

dat <- read.csv(data_file, check.names = FALSE)
item_names <- paste0("support", 1:9, "_item")
if (!all(item_names %in% names(dat))) {
  stop("The analysis file does not contain all nine support items.")
}

support_items <- dat[item_names]
names(support_items) <- paste0("SS", 1:9)

base_model <- paste(
  "kinship =~ SS1 + SS2 + SS3",
  "offline =~ SS4 + SS5 + SS6",
  "online =~ SS7 + SS8 + SS9",
  sep = "\n"
)

final_model <- paste(
  base_model,
  "SS5 ~~ SS8",
  "SS4 ~~ SS7",
  "SS2 ~~ SS5",
  "SS2 ~~ SS8",
  "SS1 ~~ SS7",
  sep = "\n"
)

fit_base_ml <- cfa(
  base_model,
  data = support_items,
  estimator = "ML",
  std.lv = TRUE,
  missing = "listwise"
)

fit_final_ml <- cfa(
  final_model,
  data = support_items,
  estimator = "ML",
  std.lv = TRUE,
  missing = "listwise"
)

fit_final_ml_fiml <- cfa(
  final_model,
  data = support_items,
  estimator = "ML",
  std.lv = TRUE,
  missing = "fiml"
)

fit_final_wlsmv <- cfa(
  final_model,
  data = support_items,
  estimator = "WLSMV",
  ordered = names(support_items),
  std.lv = TRUE
)

extract_fit <- function(fit, model, estimator, scaled = FALSE) {
  measures <- if (scaled) {
    c(
      chisq = "chisq.scaled",
      df = "df.scaled",
      cfi = "cfi.scaled",
      tli = "tli.scaled",
      rmsea = "rmsea.scaled",
      rmsea_ci_lower = "rmsea.ci.lower.scaled",
      rmsea_ci_upper = "rmsea.ci.upper.scaled",
      srmr = "srmr"
    )
  } else {
    c(
      chisq = "chisq",
      df = "df",
      cfi = "cfi",
      tli = "tli",
      rmsea = "rmsea",
      rmsea_ci_lower = "rmsea.ci.lower",
      rmsea_ci_upper = "rmsea.ci.upper",
      srmr = "srmr"
    )
  }
  values <- vapply(
    unname(measures),
    function(measure) as.numeric(fitMeasures(fit, measure)[1]),
    numeric(1)
  )
  names(values) <- names(measures)
  result <- data.frame(
    model = model,
    estimator = estimator,
    n = nobs(fit),
    converged = lavInspect(fit, "converged"),
    as.list(values),
    check.names = FALSE
  )
  result
}

fit_summary <- rbind(
  extract_fit(fit_base_ml, "Three-factor model", "ML, listwise"),
  extract_fit(
    fit_final_ml,
    "Three-factor model with five residual covariances",
    "ML, listwise"
  ),
  extract_fit(
    fit_final_ml_fiml,
    "Three-factor model with five residual covariances",
    "ML, FIML"
  ),
  extract_fit(
    fit_final_wlsmv,
    "Three-factor model with five residual covariances",
    "WLSMV",
    scaled = TRUE
  )
)

standardized_parameters <- parameterEstimates(
  fit_final_ml,
  standardized = TRUE,
  ci = TRUE
)
standardized_parameters <- standardized_parameters[
  standardized_parameters$op == "=~" |
    (
      standardized_parameters$op == "~~" &
        standardized_parameters$lhs %in% c("kinship", "offline", "online") &
        standardized_parameters$rhs %in% c("kinship", "offline", "online") &
        standardized_parameters$lhs != standardized_parameters$rhs
    ),
  c("lhs", "op", "rhs", "est", "se", "pvalue", "ci.lower", "ci.upper", "std.all")
]

raw_alpha <- function(x) {
  x <- x[complete.cases(x), , drop = FALSE]
  k <- ncol(x)
  k / (k - 1) * (1 - sum(vapply(x, var, numeric(1))) / var(rowSums(x)))
}

questionnaire_items <- list("GAD-7" = paste0("gad", 1:7, "_item"),
                            "MPIQ" = paste0("mpiq", 1:8, "_item"))
questionnaire_reliability <- do.call(rbind, lapply(names(questionnaire_items), function(scale) {
  x <- dat[questionnaire_items[[scale]]]
  data.frame(scale = scale, items = ncol(x), complete_records = sum(complete.cases(x)),
             cronbach_alpha = raw_alpha(x))
}))
write.csv(questionnaire_reliability,
          file.path(output_dir, "questionnaire_internal_consistency.csv"), row.names = FALSE)

loading_rows <- standardized_parameters[standardized_parameters$op == "=~", ]
reliability_from_loadings <- function(factor_name) {
  lambda <- loading_rows$std.all[loading_rows$lhs == factor_name]
  error_variance <- 1 - lambda^2
  c(
    composite_reliability = sum(lambda)^2 / (sum(lambda)^2 + sum(error_variance)),
    average_variance_extracted = sum(lambda^2) / (sum(lambda^2) + sum(error_variance))
  )
}

kinship_reliability <- reliability_from_loadings("kinship")
offline_reliability <- reliability_from_loadings("offline")
online_reliability <- reliability_from_loadings("online")

reliability <- data.frame(
  scale = c("Kinship support", "Offline-friend support", "Online-friend support", "Total support"),
  items = c(3L, 3L, 3L, 9L),
  cronbach_alpha = c(
    raw_alpha(support_items[paste0("SS", 1:3)]),
    raw_alpha(support_items[paste0("SS", 4:6)]),
    raw_alpha(support_items[paste0("SS", 7:9)]),
    raw_alpha(support_items)
  ),
  composite_reliability = c(
    kinship_reliability[["composite_reliability"]],
    offline_reliability[["composite_reliability"]],
    online_reliability[["composite_reliability"]],
    NA_real_
  ),
  average_variance_extracted = c(
    kinship_reliability[["average_variance_extracted"]],
    offline_reliability[["average_variance_extracted"]],
    online_reliability[["average_variance_extracted"]],
    NA_real_
  ),
  standardized_loading_minimum = c(
    min(loading_rows$std.all[loading_rows$lhs == "kinship"]),
    min(loading_rows$std.all[loading_rows$lhs == "offline"]),
    min(loading_rows$std.all[loading_rows$lhs == "online"]),
    NA_real_
  ),
  standardized_loading_maximum = c(
    max(loading_rows$std.all[loading_rows$lhs == "kinship"]),
    max(loading_rows$std.all[loading_rows$lhs == "offline"]),
    max(loading_rows$std.all[loading_rows$lhs == "online"]),
    NA_real_
  )
)

write.csv(fit_summary, file.path(output_dir, "cfa_fit_indices.csv"), row.names = FALSE)
write.csv(
  standardized_parameters,
  file.path(output_dir, "cfa_standardized_parameters.csv"),
  row.names = FALSE
)
write.csv(reliability, file.path(output_dir, "internal_consistency.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
