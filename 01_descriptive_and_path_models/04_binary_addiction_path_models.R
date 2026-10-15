required_packages(c("lavaan", "flextable", "officer"))

analysis_file <- Sys.getenv(
  "ANALYSIS_DATA_FILE",
  unset = file.path("..", "analysis_data.csv")
)
output_dir <- Sys.getenv(
  "BINARY_PATH_OUTPUT_DIR",
  unset = file.path("..", "results", "01_binary_addiction_path_models")
)
bootstrap_reps <- as.integer(Sys.getenv("BOOTSTRAP_REPS", unset = "5000"))
save_models <- identical(
  tolower(Sys.getenv("SAVE_MODEL_OBJECTS", unset = "false")),
  "true"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(analysis_file)) stop("Analysis data not found: ", analysis_file)
if (is.na(bootstrap_reps) || bootstrap_reps < 1L) {
  stop("BOOTSTRAP_REPS must be a positive integer.")
}

dat <- read.csv(analysis_file, check.names = FALSE)
required <- c(
  "mpiq_total", "addiction_group", "gad7_total", "survey_year", "sex",
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
expected_group <- ifelse(dat$mpiq_total > 32, "Addicted", "Non-addicted")
if (any(as.character(dat$addiction_group) != expected_group)) {
  stop("addiction_group must use MPIQ >32 for Addicted and MPIQ <=32 for Non-addicted.")
}

dat$sex_group <- factor(dat$sex, levels = c("Female", "Male"))
dat$addiction_group <- factor(
  dat$addiction_group,
  levels = c("Non-addicted", "Addicted")
)
if (anyNA(dat$sex_group)) stop("sex must contain only Female and Male.")
if (anyNA(dat$addiction_group)) {
  stop("addiction_group must contain only Non-addicted and Addicted.")
}

dat$x <- as.numeric(dat$addiction_group == "Addicted")
dat$followup_indicator <- as.numeric(dat$survey_year == max(dat$survey_year))

cell_counts <- table(dat$sex_group, dat$addiction_group)
if (any(cell_counts == 0L)) {
  stop("Each sex-by-addiction group must contain complete observations.")
}
for (sex_value in levels(dat$sex_group)) {
  group_data <- dat[dat$sex_group == sex_value, , drop = FALSE]
  if (length(unique(group_data$survey_year)) < 2L) {
    stop("Both survey years are required in the ", sex_value, " group.")
  }
  for (variable in c(
    "kinship_support", "offline_friend_support", "online_friend_support"
  )) {
    for (addiction_value in levels(dat$addiction_group)) {
      values <- group_data[
        group_data$addiction_group == addiction_value,
        variable
      ]
      if (length(unique(values)) < 2L) {
        stop(
          variable, " has no variation among ", sex_value, " ",
          tolower(addiction_value), " participants."
        )
      }
    }
  }
}

all_supports <- c(
  kinship_support = "Kinship support",
  offline_friend_support = "Offline-friend support",
  online_friend_support = "Online-friend support"
)
supports <- all_supports
selected_support <- Sys.getenv("SUPPORT_SOURCE", unset = "all")
if (selected_support != "all") {
  if (!selected_support %in% names(supports)) {
    stop("SUPPORT_SOURCE must be all or one of: ",
         paste(names(supports), collapse = ", "))
  }
  supports <- supports[selected_support]
}

model_syntax <- paste(
  "m ~ c(a_f, a_m) * x + c(year_m_f, year_m_m) * followup_indicator",
  paste0(
    "gad7_total ~ c(cp_f, cp_m) * x + c(b1_f, b1_m) * m + ",
    "c(b2_f, b2_m) * xm + c(year_y_f, year_y_m) * followup_indicator"
  ),
  "b_add_f := b1_f + b2_f",
  "b_add_m := b1_m + b2_m",
  "b_add_diff := b_add_f - b_add_m",
  "product_nonadd_f := a_f * b1_f",
  "product_nonadd_m := a_m * b1_m",
  "product_nonadd_diff := product_nonadd_f - product_nonadd_m",
  "product_add_f := a_f * (b1_f + b2_f)",
  "product_add_m := a_m * (b1_m + b2_m)",
  "product_add_diff := product_add_f - product_add_m",
  "index_f := a_f * b2_f",
  "index_m := a_m * b2_m",
  "index_diff := index_f - index_m",
  "diff_a := a_f - a_m",
  "diff_b1 := b1_f - b1_m",
  "diff_b2 := b2_f - b2_m",
  "diff_cp := cp_f - cp_m",
  sep = "\n"
)

all_estimates <- list()
run_status <- list()
model_objects <- list()

for (i in seq_along(supports)) {
  variable <- names(supports)[i]
  support_label <- unname(supports[i])
  model_data <- dat
  model_data$m <- model_data[[variable]] - mean(model_data[[variable]])
  model_data$xm <- model_data$x * model_data$m

  seed <- 20260920L + match(variable, names(all_supports))
  set.seed(seed)
  started <- Sys.time()
  fit <- lavaan::sem(
    model_syntax,
    data = model_data,
    group = "sex_group",
    group.label = c("Female", "Male"),
    fixed.x = TRUE,
    se = "bootstrap",
    bootstrap = bootstrap_reps
  )
  finished <- Sys.time()

  estimates <- lavaan::parameterEstimates(
    fit,
    standardized = TRUE,
    ci = TRUE,
    boot.ci.type = "bca.simple"
  )
  estimates$support_variable <- variable
  estimates$support <- support_label
  estimates$seed <- seed
  estimates$bootstrap_replicates <- bootstrap_reps
  estimates$converged <- lavaan::lavInspect(fit, "converged")
  all_estimates[[variable]] <- estimates
  model_objects[[variable]] <- fit

  bootstrap_coefficients <- fit@boot$coef
  failed_bootstraps <- if (is.null(bootstrap_coefficients)) {
    NA_integer_
  } else {
    sum(!stats::complete.cases(bootstrap_coefficients))
  }
  run_status[[variable]] <- data.frame(
    support = support_label,
    n = lavaan::lavInspect(fit, "ntotal"),
    female_n = unname(cell_counts["Female", "Non-addicted"] +
                        cell_counts["Female", "Addicted"]),
    male_n = unname(cell_counts["Male", "Non-addicted"] +
                      cell_counts["Male", "Addicted"]),
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
write.csv(
  status,
  file.path(output_dir, "run_status.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
if (save_models) {
  saveRDS(model_objects, file.path(output_dir, "fitted_models.rds"))
}

path_labels <- c(
  "a_f", "a_m", "b1_f", "b1_m", "b2_f", "b2_m", "cp_f", "cp_m"
)
path_summary <- combined[
  combined$op == "~" & combined$label %in% path_labels,
  c(
    "support", "group", "label", "lhs", "rhs", "est", "se", "pvalue",
    "ci.lower", "ci.upper", "std.all", "seed", "bootstrap_replicates",
    "converged"
  )
]
path_summary$sex <- ifelse(path_summary$group == 1L, "Female", "Male")
path_export <- path_summary[, c(
  "support", "sex", "label", "lhs", "rhs", "est", "std.all", "se",
  "pvalue", "ci.lower", "ci.upper", "seed", "bootstrap_replicates",
  "converged"
)]
names(path_export) <- c(
  "support", "sex", "path", "outcome", "predictor",
  "unstandardized_estimate", "standardized_estimate", "standard_error",
  "normal_theory_p", "bootstrap_ci_lower", "bootstrap_ci_upper", "seed",
  "bootstrap_replicates", "converged"
)

defined_names <- c(
  "b_add_f", "b_add_m", "b_add_diff",
  "product_nonadd_f", "product_nonadd_m", "product_nonadd_diff",
  "product_add_f", "product_add_m", "product_add_diff",
  "index_f", "index_m", "index_diff",
  "diff_a", "diff_b1", "diff_b2", "diff_cp"
)
defined_summary <- combined[
  combined$op == ":=" & combined$lhs %in% defined_names,
  c(
    "support", "lhs", "est", "se", "pvalue", "ci.lower", "ci.upper",
    "seed", "bootstrap_replicates", "converged"
  )
]
names(defined_summary) <- c(
  "support", "quantity", "estimate", "standard_error", "normal_theory_p",
  "bootstrap_ci_lower", "bootstrap_ci_upper", "seed",
  "bootstrap_replicates", "converged"
)

write.csv(
  path_export,
  file.path(output_dir, "binary_path_coefficient_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  defined_summary,
  file.path(output_dir, "binary_defined_parameter_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

get_path <- function(support, label) {
  row <- path_summary[
    path_summary$support == support & path_summary$label == label,
  ]
  if (nrow(row) != 1L) {
    stop("Path coefficient not found: ", support, " / ", label)
  }
  row
}

get_defined <- function(support, quantity) {
  row <- defined_summary[
    defined_summary$support == support & defined_summary$quantity == quantity,
  ]
  if (nrow(row) != 1L) {
    stop("Defined parameter not found: ", support, " / ", quantity)
  }
  row
}

format_estimate <- function(x, digits) sprintf(paste0("%.", digits, "f"), x)
format_ci <- function(lower, upper, digits) {
  paste0(
    format_estimate(lower, digits),
    " to ",
    format_estimate(upper, digits)
  )
}

publication_rows <- list()
row_number <- 0L

add_section <- function(label) {
  row_number <<- row_number + 1L
  publication_rows[[row_number]] <<- data.frame(
    row_type = "header",
    Quantity = label,
    Description = "",
    `Sex or comparison` = "",
    `B or product` = "",
    `Standardized beta` = "",
    `P value` = "",
    `95% BC bootstrap CI` = "",
    check.names = FALSE
  )
}

add_result <- function(
    quantity, description, comparison, estimate, std_beta,
    p_value, lower, upper, digits = 4L) {
  row_number <<- row_number + 1L
  publication_rows[[row_number]] <<- data.frame(
    row_type = "body",
    Quantity = quantity,
    Description = description,
    `Sex or comparison` = comparison,
    `B or product` = format_estimate(estimate, digits),
    `Standardized beta` = if (is.na(std_beta)) {
      ""
    } else {
      format_estimate(std_beta, 3L)
    },
    `P value` = format_p(p_value),
    `95% BC bootstrap CI` = format_ci(lower, upper, digits),
    check.names = FALSE
  )
}

path_definitions <- list(
  a = list(
    description = "Addicted versus non-addicted -> support",
    female = "a_f", male = "a_m", difference = "diff_a"
  ),
  b1 = list(
    description = "Support -> GAD-7, non-addicted",
    female = "b1_f", male = "b1_m", difference = "diff_b1"
  ),
  b2 = list(
    description = "Addicted minus non-addicted support slope",
    female = "b2_f", male = "b2_m", difference = "diff_b2"
  ),
  cp = list(
    description = "Addicted versus non-addicted -> GAD-7",
    female = "cp_f", male = "cp_m", difference = "diff_cp"
  )
)

for (support in unname(supports)) {
  add_section(support)

  for (quantity in names(path_definitions)) {
    definition <- path_definitions[[quantity]]
    female <- get_path(support, definition$female)
    male <- get_path(support, definition$male)
    difference <- get_defined(support, definition$difference)
    display_quantity <- switch(
      quantity,
      b1 = "b1",
      b2 = "b2",
      cp = "c'",
      quantity
    )
    add_result(
      display_quantity, definition$description, "Female",
      female$est, female$std.all, female$pvalue,
      female$ci.lower, female$ci.upper
    )
    add_result(
      display_quantity, definition$description, "Male",
      male$est, male$std.all, male$pvalue,
      male$ci.lower, male$ci.upper
    )
    add_result(
      display_quantity, definition$description, "Female - male",
      difference$estimate, NA_real_, difference$normal_theory_p,
      difference$bootstrap_ci_lower, difference$bootstrap_ci_upper
    )
  }

  for (sex_code in c("f", "m", "diff")) {
    row <- get_defined(support, paste0("b_add_", sex_code))
    comparison <- c(
      f = "Female", m = "Male", diff = "Female - male"
    )[[sex_code]]
    add_result(
      "b1 + b2", "Support -> GAD-7, addicted", comparison,
      row$estimate, NA_real_, row$normal_theory_p,
      row$bootstrap_ci_lower, row$bootstrap_ci_upper
    )
  }

  for (sex_code in c("f", "m", "diff")) {
    row <- get_defined(support, paste0("index_", sex_code))
    comparison <- c(
      f = "Female", m = "Male", diff = "Female - male"
    )[[sex_code]]
    add_result(
      "a x b2", "Moderated product", comparison,
      row$estimate, NA_real_, row$normal_theory_p,
      row$bootstrap_ci_lower, row$bootstrap_ci_upper,
      digits = 6L
    )
  }
}

publication_table <- do.call(rbind, publication_rows)
write.csv(
  publication_table[, -1],
  file.path(output_dir, "binary_sensitivity_publication_source_data.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

save_table_docx(
  publication_table,
  paste(
    "Binary addiction-status sensitivity analysis of sex-grouped",
    "moderated mediation path models"
  ),
  paste(
    paste0(
      "Female n = ", format(sum(dat$sex_group == "Female"), big.mark = ","),
      "; male n = ", format(sum(dat$sex_group == "Male"), big.mark = ","),
      ". Addiction status was coded 0 for MPIQ scores of 32 or lower and 1 for scores greater than 32."
    ),
    "Support scores were centred at their full-sample means, and survey year was included in both equations.",
    "Path b1 is the support coefficient in the non-addicted group; b1 + b2 is the corresponding coefficient in the addicted group.",
    "Path b2 is the addicted-minus-non-addicted difference in the support coefficient.",
    "B and product terms are unstandardized. P values use two-sided normal-approximation tests with bootstrap standard errors; intervals are bias-corrected bootstrap intervals.",
    "The products summarize associations measured at the same assessment."
  ),
  file.path(output_dir, "binary_addiction_status_sensitivity_publication.docx")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)
message(
  "Binary addiction-status path results written to ",
  normalizePath(output_dir, winslash = "/", mustWork = FALSE)
)
