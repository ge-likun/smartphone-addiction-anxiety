optional_covariates <- c(
  "grade_code", "major_code", "hukou_code", "father_education_code",
  "mother_education_code", "only_child_code", "ipaq_category_code"
)

if (!all(optional_covariates %in% names(analysis_data)) ||
    any(vapply(analysis_data[intersect(optional_covariates, names(analysis_data))], function(x) all(is.na(x)), logical(1)))) {
  writeLines(
    "Supplementary Figure 5 was not generated because one or more optional covariates were absent.",
    file.path(output_dirs$logs, "Supplementary_Figure_5_status.txt")
  )
} else {
  covariate_data <- analysis_data |>
    mutate(
      grade = factor(case_when(
        .data$grade_code %in% 1:4 ~ as.character(.data$grade_code),
        .data$grade_code %in% 5:6 ~ "5+",
        TRUE ~ NA_character_
      ), levels = c("1", "2", "3", "4", "5+")),
      major = factor(ifelse(.data$major_code %in% 1:3, .data$major_code, NA)),
      hukou = factor(ifelse(.data$hukou_code %in% 1:2, .data$hukou_code, NA)),
      father_education = factor(ifelse(.data$father_education_code %in% 1:4, .data$father_education_code, NA)),
      mother_education = factor(ifelse(.data$mother_education_code %in% 1:4, .data$mother_education_code, NA)),
      only_child = factor(ifelse(.data$only_child_code %in% 1:2, .data$only_child_code, NA)),
      physical_activity = factor(ifelse(.data$ipaq_category_code %in% 1:3, .data$ipaq_category_code, NA))
    )

  covariate_rhs <- c(
    "pb(kinship_support)", "pb(offline_friend_support)", "pb(online_friend_support)",
    "survey_year", "grade", "major", "hukou", "father_education",
    "mother_education", "only_child", "physical_activity"
  )
  primary_rhs <- c(
    "pb(kinship_support)", "pb(offline_friend_support)",
    "pb(online_friend_support)", "survey_year"
  )

  fit_adjusted_model <- function(sex_value) {
    model_data <- covariate_data |>
      filter(.data$addiction_group == "Addicted", .data$sex == .env$sex_value) |>
      filter(if_all(all_of(c("gad7_total", "kinship_support", "offline_friend_support", "online_friend_support", "survey_year", "grade", "major", "hukou", "father_education", "mother_education", "only_child", "physical_activity")), ~ !is.na(.x))) |>
      droplevels()
    mu_formula <- as.formula(paste("gad7_total ~", paste(covariate_rhs, collapse = " + ")))
    parameter_formula <- as.formula(paste("~", paste(covariate_rhs, collapse = " + ")))
    primary_parameter_formula <- as.formula(paste("~", paste(primary_rhs, collapse = " + ")))
    primary_model <- primary_models[[paste0("Addicted_", sex_value)]]
    fit_seed <- residual_seed_base + 200000L + match(sex_value, c("Male", "Female"))
    if (!nrow(model_data)) stop("No complete covariate-sensitivity rows for ", sex_value)
    if (nrow(model_data) != primary_model$N) {
      if (identical(Sys.getenv("GAMLSS_CACHE_ONLY", "false"), "true"))
        stop("Matched-sample primary fit is needed; fitting is disabled.")
      set.seed(fit_seed)
      primary_model <- gamlss(
        as.formula(paste("gad7_total ~", paste(primary_rhs, collapse = " + "))),
        sigma.formula = primary_parameter_formula, nu.formula = primary_parameter_formula,
        family = ZAGA(), data = model_data, control = fit_control())
      if (!isTRUE(primary_model$converged)) stop("Matched-sample primary model did not converge.")
    }
    matched_primary_models[[sex_value]] <<- primary_model
    cache_file <- file.path(
      output_dirs$models,
      paste0("Addicted_", sex_value, "__ZAGA_covariate_adjusted__", model_cache_tag, ".rds")
    )
    cache_identity <- gamlss_cache_identity(
      model_data, list(mu = mu_formula, sigma = parameter_formula, nu = parameter_formula))
    model <- if (file.exists(cache_file)) {
      cached <- readRDS(cache_file)
      verify_gamlss_cache_identity(cached, cache_identity, cache_file)
      cached
    } else {
      if (identical(Sys.getenv("GAMLSS_CACHE_ONLY", "false"), "true"))
        stop("Missing covariate-adjusted cache; fitting is disabled: ", cache_file)
      set.seed(fit_seed)
      fit_stage <- function(mu, sigma, nu, start) {
        stage_model <- gamlss(
          mu, sigma.formula = sigma, nu.formula = nu, family = ZAGA(),
          data = model_data, start.from = start,
          control = gamlss.control(n.cyc = 200, mu.step = 0.5, sigma.step = 0.25, nu.step = 0.5, trace = FALSE)
        )
        if (!isTRUE(stage_model$converged)) stop("Covariate-adjusted model did not converge for ", sex_value)
        stage_model
      }
      stage_mu <- fit_stage(mu_formula, primary_parameter_formula, primary_parameter_formula, primary_model)
      stage_sigma <- fit_stage(mu_formula, parameter_formula, primary_parameter_formula, stage_mu)
      fitted_model <- fit_stage(mu_formula, parameter_formula, parameter_formula, stage_sigma)
      attr(fitted_model, "cache_identity") <- cache_identity
      attr(fitted_model, "fit_seed") <- fit_seed
      saveRDS(fitted_model, cache_file)
      fitted_model
    }
    expected_mu <- gsub("\\s+", " ", formula_text(mu_formula))
    expected_rhs <- gsub("\\s+", " ", formula_text(parameter_formula))
    formulas_ok <- all(vapply(c("mu", "sigma", "nu"), function(parameter) {
      observed <- gsub("\\s+", " ", formula_text(model[[paste0(parameter, ".formula")]]))
      expected <- if (parameter == "mu") expected_mu else expected_rhs
      identical(observed, expected)
    }, logical(1)))
    if (!identical(as.character(model$family[[1L]]), "ZAGA") || !isTRUE(model$converged) || model$N != nrow(model_data) || !formulas_ok) {
      stop("Covariate-adjusted model verification failed for ", sex_value)
    }
    model
  }

  matched_primary_models <- list()
  adjusted_models <- list(Male = fit_adjusted_model("Male"), Female = fit_adjusted_model("Female"))
  sensitivity_curves <- bind_rows(lapply(c("Male", "Female"), function(sex_value) {
    bind_rows(
      full_model_curve(matched_primary_models[[sex_value]], "kinship_support", "mu") |> mutate(specification = "Primary"),
      full_model_curve(adjusted_models[[sex_value]], "kinship_support", "mu") |> mutate(specification = "Covariate-adjusted")
    ) |>
      mutate(sex = sex_value, .before = 1)
  })) |>
    mutate(
      sex = factor(.data$sex, levels = c("Male", "Female")),
      specification = factor(.data$specification, levels = c("Primary", "Covariate-adjusted"))
    )
  write.csv(sensitivity_curves, file.path(output_dirs$source_data, "Supplementary_Figure_5_covariate_sensitivity.csv"), row.names = FALSE)

  colours <- c(Primary = "#4D4D4D", `Covariate-adjusted` = "#009E73")
  if (max(abs(sensitivity_curves$ratio[sensitivity_curves$support_score == 0] - 1)) > 1e-10) {
    stop("Supplementary Figure 5 is not referenced to support score 0")
  }
  sensitivity_limits <- c(0, 1.08 * max(c(1, sensitivity_curves$ratio_high), na.rm = TRUE))
  supplementary_figure_5 <- ggplot(
    sensitivity_curves,
    aes(.data$support_score, .data$ratio, colour = .data$specification, fill = .data$specification, linetype = .data$specification, group = .data$specification)
  ) +
    geom_ribbon(aes(ymin = .data$ratio_low, ymax = .data$ratio_high), alpha = 0.09, linewidth = 0, show.legend = FALSE) +
    geom_line(linewidth = 0.72, lineend = "round") +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.34, colour = "grey30") +
    facet_wrap(~ sex, nrow = 1, labeller = as_labeller(c(Male = "Addicted men", Female = "Addicted women"))) +
    scale_colour_manual(values = colours, name = NULL) +
    scale_fill_manual(values = colours, name = NULL) +
    scale_linetype_manual(values = c(Primary = "solid", `Covariate-adjusted` = "longdash"), name = NULL) +
    scale_x_continuous(breaks = c(0, 5, 10, 15), expand = expansion(mult = 0)) +
    scale_y_continuous(limits = sensitivity_limits, expand = expansion(mult = 0)) +
    coord_cartesian(xlim = c(-0.5, 15.5), expand = FALSE, clip = "off") +
    labs(x = "Kinship-support score", y = parameter_labels[["mu"]]) +
    theme_manuscript() +
    theme(
      strip.background = element_rect(fill = "grey96", colour = "grey45"),
      strip.text = element_text(face = "plain"),
      legend.position = "bottom"
    )
  save_publication_figure(supplementary_figure_5, file.path(output_dirs$figures, "Supplementary_Figure_5_covariate_sensitivity_score0"), 183, 82)
  writeLines("Supplementary Figure 5 was generated.", file.path(output_dirs$logs, "Supplementary_Figure_5_status.txt"))
}
