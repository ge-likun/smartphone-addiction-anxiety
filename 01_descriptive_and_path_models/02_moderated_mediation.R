required_packages(c("lavaan", "flextable", "officer"))

analysis_file <- Sys.getenv(
  "ANALYSIS_DATA_FILE",
  unset = file.path("..", "analysis_data.csv")
)
output_dir <- Sys.getenv(
  "MEDIATION_OUTPUT_DIR",
  unset = file.path("..", "results", "01_path_models")
)
bootstrap_reps <- as.integer(Sys.getenv("BOOTSTRAP_REPS", unset = "5000"))
save_models <- identical(tolower(Sys.getenv("SAVE_MODEL_OBJECTS", unset = "false")), "true")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(analysis_file)) stop("Analysis data not found: ", analysis_file)
if (is.na(bootstrap_reps) || bootstrap_reps < 1L) stop("BOOTSTRAP_REPS must be a positive integer.")

dat <- read.csv(analysis_file, check.names = FALSE)
required <- c(
  "mpiq_total", "gad7_total", "survey_year", "sex",
  "kinship_support", "offline_friend_support", "online_friend_support"
)
check_columns(dat, required, "analysis_data.csv")
check_support_scores(dat)
complete_rows <- stats::complete.cases(dat[, required])
if (any(!complete_rows)) {
  message(sum(!complete_rows), " rows with missing model variables were excluded.")
  dat <- dat[complete_rows, , drop = FALSE]
}
if (nrow(dat) == 0L) stop("No complete observations remain.")
if (any(!is.finite(dat$mpiq_total)) ||
    any(dat$mpiq_total != floor(dat$mpiq_total)) ||
    any(dat$mpiq_total < 8 | dat$mpiq_total > 56)) {
  stop("mpiq_total must contain integer scores from 8 to 56.")
}
wave_years <- sort(unique(dat$survey_year))
if (length(wave_years) != 2L) stop("The path model specification requires two survey waves.")

dat$sex_group <- factor(dat$sex, levels = c("Female", "Male"))
if (anyNA(dat$sex_group)) stop("sex must contain only Female and Male.")
dat$x <- dat$mpiq_total - mean(dat$mpiq_total)
dat$followup_indicator <- as.numeric(dat$survey_year == max(dat$survey_year))
x_sd <- stats::sd(dat$x)

sex_counts <- table(dat$sex_group)
if (any(sex_counts == 0L)) {
  stop("Each sex group must contain complete observations.")
}
if (length(unique(dat$survey_year)) < 2L) {
  stop("At least two survey years are required for year adjustment.")
}
for (sex_value in levels(dat$sex_group)) {
  group_data <- dat[dat$sex_group == sex_value, , drop = FALSE]
  for (variable in c("mpiq_total", "kinship_support", "offline_friend_support", "online_friend_support")) {
    if (length(unique(group_data[[variable]])) < 2L) {
      stop(variable, " has no variation in the ", sex_value, " group.")
    }
  }
}

supports <- c(
  kinship_support = "Kinship support",
  offline_friend_support = "Offline-friend support",
  online_friend_support = "Online-friend support"
)

model_syntax <- function(sd_x) {
  paste0(
    "m ~ c(a_f, a_m) * x + c(year_m_f, year_m_m) * followup_indicator\n",
    "gad7_total ~ c(cp_f, cp_m) * x + c(b1_f, b1_m) * m + ",
    "c(b2_f, b2_m) * xm + c(year_y_f, year_y_m) * followup_indicator\n",
    "index_f := a_f * b2_f\n",
    "index_m := a_m * b2_m\n",
    "index_diff := index_f - index_m\n",
    "cond_low_f := a_f * (b1_f + b2_f * ", format(-sd_x, scientific = FALSE), ")\n",
    "cond_mean_f := a_f * b1_f\n",
    "cond_high_f := a_f * (b1_f + b2_f * ", format(sd_x, scientific = FALSE), ")\n",
    "cond_low_m := a_m * (b1_m + b2_m * ", format(-sd_x, scientific = FALSE), ")\n",
    "cond_mean_m := a_m * b1_m\n",
    "cond_high_m := a_m * (b1_m + b2_m * ", format(sd_x, scientific = FALSE), ")\n",
    "cond_low_diff := cond_low_f - cond_low_m\n",
    "cond_mean_diff := cond_mean_f - cond_mean_m\n",
    "cond_high_diff := cond_high_f - cond_high_m\n",
    "diff_a := a_f - a_m\n",
    "diff_b1 := b1_f - b1_m\n",
    "diff_b2 := b2_f - b2_m\n",
    "diff_cp := cp_f - cp_m"
  )
}

all_estimates <- list()
run_status <- list()
model_objects <- list()

for (i in seq_along(supports)) {
  variable <- names(supports)[i]
  support_label <- unname(supports[i])
  model_data <- dat
  model_data$m <- model_data[[variable]] - mean(model_data[[variable]])
  model_data$xm <- model_data$x * model_data$m
  seed <- 20260908L + i
  set.seed(seed)
  started <- Sys.time()
  fit <- lavaan::sem(
    model_syntax(x_sd),
    data = model_data,
    group = "sex_group",
    group.label = c("Female", "Male"),
    fixed.x = TRUE,
    se = "bootstrap",
    bootstrap = bootstrap_reps
  )
  finished <- Sys.time()

  estimates <- lavaan::parameterEstimates(
    fit, standardized = TRUE, ci = TRUE, boot.ci.type = "bca.simple"
  )
  estimates$support_variable <- variable
  estimates$support <- support_label
  estimates$seed <- seed
  estimates$bootstrap_replicates <- bootstrap_reps
  estimates$converged <- lavaan::lavInspect(fit, "converged")
  estimates$x_sd <- x_sd
  all_estimates[[variable]] <- estimates
  model_objects[[variable]] <- fit

  bootstrap_coefficients <- fit@boot$coef
  failed_bootstraps <- if (is.null(bootstrap_coefficients)) NA_integer_ else sum(!stats::complete.cases(bootstrap_coefficients))
  run_status[[variable]] <- data.frame(
    support = support_label,
    n = lavaan::lavInspect(fit, "ntotal"),
    seed = seed,
    bootstrap_replicates = bootstrap_reps,
    converged = lavaan::lavInspect(fit, "converged"),
    failed_bootstraps = failed_bootstraps,
    elapsed_minutes = as.numeric(difftime(finished, started, units = "mins"))
  )

  write.csv(
    estimates,
    file.path(output_dir, paste0(variable, "_parameter_estimates.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

combined <- do.call(rbind, all_estimates)
status <- do.call(rbind, run_status)
write.csv(status, file.path(output_dir, "run_status.csv"), row.names = FALSE, fileEncoding = "UTF-8")
if (save_models) saveRDS(model_objects, file.path(output_dir, "fitted_models.rds"))

defined_names <- c(
  "index_f", "index_m", "index_diff",
  "cond_low_f", "cond_mean_f", "cond_high_f",
  "cond_low_m", "cond_mean_m", "cond_high_m",
  "cond_low_diff", "cond_mean_diff", "cond_high_diff",
  "diff_a", "diff_b1", "diff_b2", "diff_cp"
)
defined_summary <- combined[
  combined$op == ":=" & combined$lhs %in% defined_names,
  c("support", "lhs", "est", "se", "pvalue", "ci.lower", "ci.upper",
    "seed", "bootstrap_replicates", "converged", "x_sd")
]
path_summary <- combined[
  combined$op == "~" & combined$label %in% c(
    "a_f", "a_m", "b1_f", "b1_m", "b2_f", "b2_m", "cp_f", "cp_m"
  ),
  c("support", "group", "label", "lhs", "rhs", "est", "se", "pvalue",
    "ci.lower", "ci.upper", "std.all", "seed", "bootstrap_replicates", "converged")
]
path_summary$sex <- ifelse(path_summary$group == 1L, "Female", "Male")
path_export <- path_summary[, c(
  "support", "sex", "label", "lhs", "rhs", "est", "std.all", "se", "pvalue",
  "ci.lower", "ci.upper", "seed", "bootstrap_replicates", "converged"
)]
names(path_export) <- c(
  "support", "sex", "path", "outcome", "predictor", "unstandardized_estimate",
  "standardized_estimate", "standard_error", "normal_theory_p",
  "bootstrap_ci_lower", "bootstrap_ci_upper", "seed", "bootstrap_replicates", "converged"
)
defined_export <- defined_summary[, c(
  "support", "lhs", "est", "se", "pvalue", "ci.lower", "ci.upper",
  "seed", "bootstrap_replicates", "converged", "x_sd"
)]
names(defined_export) <- c(
  "support", "quantity", "estimate", "standard_error", "normal_theory_p",
  "bootstrap_ci_lower", "bootstrap_ci_upper", "seed", "bootstrap_replicates",
  "converged", "x_sd"
)
write.csv(path_export, file.path(output_dir, "path_coefficient_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(defined_export, file.path(output_dir, "conditional_indirect_association_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

get_path <- function(support, label) {
  row <- path_summary[path_summary$support == support & path_summary$label == label, ]
  if (nrow(row) != 1L) stop("Path coefficient not found: ", support, " / ", label)
  row
}
get_defined <- function(support, quantity) {
  row <- defined_summary[defined_summary$support == support & defined_summary$lhs == quantity, ]
  if (nrow(row) != 1L) stop("Defined parameter not found: ", support, " / ", quantity)
  row
}
format_estimate <- function(x, digits) sprintf(paste0("%.", digits, "f"), x)
format_ci <- function(lower, upper, digits) {
  paste0(format_estimate(lower, digits), " to ", format_estimate(upper, digits))
}

publication_rows <- list()
row_number <- 0L
add_section <- function(label) {
  row_number <<- row_number + 1L
  publication_rows[[row_number]] <<- data.frame(
    row_type = "header", Quantity = label, Description = "", `Sex or comparison` = "",
    `B or product` = "", `Standardized beta` = "", `95% BC bootstrap CI` = "",
    check.names = FALSE
  )
}
add_result <- function(quantity, description, comparison, estimate, std_beta, lower, upper, digits) {
  row_number <<- row_number + 1L
  publication_rows[[row_number]] <<- data.frame(
    row_type = "body", Quantity = quantity, Description = description,
    `Sex or comparison` = comparison,
    `B or product` = format_estimate(estimate, digits),
    `Standardized beta` = if (is.na(std_beta)) "" else format_estimate(std_beta, 3),
    `95% BC bootstrap CI` = format_ci(lower, upper, digits),
    check.names = FALSE
  )
}

path_definitions <- list(
  a = list(description = "MPIQ score -> support", female = "a_f", male = "a_m", difference = "diff_a"),
  b1 = list(description = "Support -> GAD-7", female = "b1_f", male = "b1_m", difference = "diff_b1"),
  b2 = list(description = "MPIQ score x support -> GAD-7", female = "b2_f", male = "b2_m", difference = "diff_b2"),
  cp = list(description = "MPIQ score -> GAD-7", female = "cp_f", male = "cp_m", difference = "diff_cp")
)

for (support in unname(supports)) {
  add_section(support)
  for (quantity in names(path_definitions)) {
    definition <- path_definitions[[quantity]]
    female <- get_path(support, definition$female)
    male <- get_path(support, definition$male)
    difference <- get_defined(support, definition$difference)
    display_quantity <- switch(quantity, b1 = "b1", b2 = "b2", cp = "c'", quantity)
    add_result(display_quantity, definition$description, "Female", female$est, female$std.all,
               female$ci.lower, female$ci.upper, 4)
    add_result(display_quantity, definition$description, "Male", male$est, male$std.all,
               male$ci.lower, male$ci.upper, 4)
    add_result(display_quantity, definition$description, "Female - male", difference$est, NA,
               difference$ci.lower, difference$ci.upper, 4)
  }

  for (sex_code in c("f", "m", "diff")) {
    row <- get_defined(support, paste0("index_", sex_code))
    comparison <- c(f = "Female", m = "Male", diff = "Female - male")[[sex_code]]
    add_result("a x b2", "Change in the indirect-association product", comparison,
               row$est, NA, row$ci.lower, row$ci.upper, 6)
  }
  condition_labels <- c(low = "At MPIQ mean - 1 s.d.", mean = "At MPIQ mean", high = "At MPIQ mean + 1 s.d.")
  for (condition in names(condition_labels)) {
    for (sex_code in c("f", "m", "diff")) {
      row <- get_defined(support, paste0("cond_", condition, "_", sex_code))
      comparison <- c(f = "Female", m = "Male", diff = "Female - male")[[sex_code]]
      add_result("a(b1 + b2X)", condition_labels[[condition]], comparison,
                 row$est, NA, row$ci.lower, row$ci.upper, 6)
    }
  }
}

publication_table <- do.call(rbind, publication_rows)
write.csv(publication_table[, -1], file.path(output_dir, "Supplementary_Table_5_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
save_table_docx(
  publication_table,
  "Supplementary Table 5 | Sex-grouped path coefficients and conditional indirect-association products",
  paste(
    paste0(
      "Female n = ", format(unname(sex_counts[["Female"]]), big.mark = ","),
      "; male n = ", format(unname(sex_counts[["Male"]]), big.mark = ","),
      ". MPIQ and support scores were mean-centred, and survey year was included in both equations."
    ),
    "B denotes an unstandardized coefficient. Product terms and female-minus-male differences are also unstandardized.",
    "Standardized beta is shown only for sex-specific path coefficients.",
    paste0("Intervals are bias-corrected bootstrap intervals based on ", format(bootstrap_reps, big.mark = ","), " resamples."),
    "The products summarize associations measured at the same assessment."
  ),
  file.path(output_dir, "Supplementary_Table_5_publication.docx")
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Moderated mediation results written to ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
