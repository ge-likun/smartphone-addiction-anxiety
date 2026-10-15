format_number <- function(x, digits = 1) formatC(x, format = "f", digits = digits, big.mark = ",")

write_publication_docx <- function(document, target) {
  temporary <- tempfile(pattern = "gamlss_table_", fileext = ".docx")
  on.exit(unlink(temporary), add = TRUE)
  print(document, target = temporary)
  if (!file.copy(temporary, target, overwrite = TRUE)) {
    stop("Unable to write table: ", target)
  }
}

residual_columns <- residual_diagnostic_summary |>
  select(
    family, model_id,
    central_deviation = central_maximum_absolute_deviation_median,
    upper_tail_deviation = upper_tail_mean_deviation_q_gt_2_5_median
  )

table6_base <- information_criteria |>
  left_join(residual_columns, by = c("family", "model_id")) |>
  transmute(
    Stratum = .data$stratum,
    Family = .data$family_label,
    n = format(.data$n, big.mark = ",", scientific = FALSE, trim = TRUE),
    Converged = .data$converged,
    AIC = format_number(.data$AIC),
    BIC = format_number(.data$BIC),
    `Maximum absolute central Q–Q deviation` = round(.data$central_deviation, 3),
    `Mean upper-tail deviation` = round(.data$upper_tail_deviation, 3),
    comparison_set = .data$comparison_set
  )

table6a <- table6_base |> filter(.data$comparison_set == "Continuous/mixed") |> select(-comparison_set)
table6b <- table6_base |> filter(.data$comparison_set == "Discrete") |> select(-comparison_set)
table6c <- boundary_summary |>
  transmute(
    Stratum = .data$stratum,
    Family = .data$family_label,
    n = format(.data$n, big.mark = ",", scientific = FALSE, trim = TRUE),
    `Observed n at 21` = .data$observed_n_at_21,
    `Expected n at/near 21*` = round(.data$expected_n_in_score21_bin, 1),
    `Expected n above 21` = round(.data$expected_n_above_21, 1)
  ) |>
  arrange(match(.data$Stratum, analysis_specs$stratum),
          match(.data$Family, vapply(family_specs, function(x) x$label, character(1))))

write.csv(table6a, file.path(output_dirs$tables, "Supplementary_Table_6A_continuous_mixed.csv"), row.names = FALSE)
write.csv(table6b, file.path(output_dirs$tables, "Supplementary_Table_6B_discrete.csv"), row.names = FALSE)
write.csv(table6c, file.path(output_dirs$tables, "Supplementary_Table_6C_boundary.csv"), row.names = FALSE)

table7 <- family_fit_status |>
  filter(.data$family == "ZAGA") |>
  transmute(
    Stratum = .data$stratum,
    n = format(.data$n, big.mark = ",", scientific = FALSE, trim = TRUE),
    Adjustment = ifelse(grepl("overall$", .data$stratum), "Sex and survey year", "Survey year"),
    Converged = .data$converged,
    `Effective df` = round(.data$effective_df, 1),
    `Log likelihood` = format_number(.data$log_likelihood),
    AIC = format_number(.data$AIC),
    BIC = format_number(.data$BIC)
  )
write.csv(table7, file.path(output_dirs$tables, "Supplementary_Table_7_primary_ZAGA_summary.csv"), row.names = FALSE)

interaction_comparisons <- list(
  Addiction_status = list(label = "Addiction status", data = analysis_data, group = "addiction_group", adjustment = c("sex", "survey_year")),
  Sex_within_nonaddicted = list(label = "Sex within non-addicted", data = filter(analysis_data, .data$addiction_group == "Non-addicted") |> droplevels(), group = "sex", adjustment = "survey_year"),
  Sex_within_addicted = list(label = "Sex within addicted", data = filter(analysis_data, .data$addiction_group == "Addicted") |> droplevels(), group = "sex", adjustment = "survey_year")
)

interaction_formulas <- function(target_parameter, target_support, group, adjustment, alternative) {
  other_supports <- setdiff(unname(support_variables), target_support)
  common_target <- c(paste0("pb(", other_supports, ")"), target_support, paste0("I(", target_support, "^2)"), adjustment, group)
  if (alternative) {
    common_target <- c(common_target, paste0(target_support, ":", group), paste0("I(", target_support, "^2):", group))
  }
  standard <- c(paste0("pb(", unname(support_variables), ")"), adjustment, group)
  rhs <- list(mu = standard, sigma = standard, nu = standard)
  rhs[[target_parameter]] <- common_target
  list(
    mu = as.formula(paste("gad7_total ~", paste(rhs$mu, collapse = " + "))),
    sigma = as.formula(paste("~", paste(rhs$sigma, collapse = " + "))),
    nu = as.formula(paste("~", paste(rhs$nu, collapse = " + ")))
  )
}

fit_interaction_model <- function(comparison_id, target_parameter, support_key, alternative) {
  comparison <- interaction_comparisons[[comparison_id]]
  target_support <- support_variables[[support_key]]
  formulas <- interaction_formulas(target_parameter, target_support, comparison$group, comparison$adjustment, alternative)
  cache_identity <- gamlss_cache_identity(comparison$data, formulas)
  suffix <- if (alternative) "alternative" else "null"
  cache_file <- file.path(
    output_dirs$models,
    paste(
      "quadratic_interaction", comparison_id, target_parameter, support_key,
      suffix, model_cache_tag, sep = "__"
    ) |> paste0(".rds")
  )
  model <- if (file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    verify_gamlss_cache_identity(cached, cache_identity, cache_file)
    cached
  } else {
    if (identical(Sys.getenv("GAMLSS_CACHE_ONLY", "false"), "true"))
      stop("Missing comparison-model cache; fitting is disabled: ", cache_file)
    fit_seed <- residual_seed_base + match(comparison_id, names(interaction_comparisons)) * 1000L +
      match(target_parameter, c("mu", "nu", "sigma")) * 100L +
      match(support_key, names(support_variables)) * 10L + as.integer(alternative)
    set.seed(fit_seed)
    fitted <- gamlss(
      formulas$mu, sigma.formula = formulas$sigma, nu.formula = formulas$nu,
      family = ZAGA(), data = comparison$data,
      control = gamlss.control(n.cyc = 150, trace = FALSE)
    )
    attr(fitted, "cache_identity") <- cache_identity
    attr(fitted, "fit_seed") <- fit_seed
    fitted
  }
  expected_formulas <- vapply(formulas, function(value) gsub("\\s+", " ", formula_text(value)), character(1))
  observed_formulas <- vapply(c("mu", "sigma", "nu"), function(parameter) {
    gsub("\\s+", " ", formula_text(model[[paste0(parameter, ".formula")]]))
  }, character(1))
  expected_mu_rhs <- sub("^[^~]+~", "~", expected_formulas[["mu"]])
  observed_mu_rhs <- sub("^[^~]+~", "~", observed_formulas[["mu"]])
  formulas_ok <- identical(expected_mu_rhs, observed_mu_rhs) &&
    identical(expected_formulas[["sigma"]], observed_formulas[["sigma"]]) &&
    identical(expected_formulas[["nu"]], observed_formulas[["nu"]])
  if (
    !identical(as.character(model$family[[1L]]), "ZAGA") ||
    !isTRUE(model$converged) || model$N != nrow(comparison$data) || !formulas_ok
  ) {
    stop("Quadratic interaction model failed: ", comparison_id, "/", target_parameter, "/", support_key, "/", suffix)
  }
  if (!file.exists(cache_file)) saveRDS(model, cache_file)
  model
}

interaction_rows <- list()
for (comparison_id in names(interaction_comparisons)) {
  for (parameter in c("mu", "nu", "sigma")) {
    for (support_key in names(support_variables)) {
      null_model <- fit_interaction_model(comparison_id, parameter, support_key, FALSE)
      alternative_model <- fit_interaction_model(comparison_id, parameter, support_key, TRUE)
      deviance_difference <- null_model$G.deviance - alternative_model$G.deviance
      df_difference <- alternative_model$df.fit - null_model$df.fit
      interaction_rows[[paste(comparison_id, parameter, support_key)]] <- tibble(
        comparison_id = comparison_id,
        Comparison = interaction_comparisons[[comparison_id]]$label,
        parameter = parameter,
        Support = sub("_", " ", support_key),
        `Deviance difference` = deviance_difference,
        `Effective-df difference` = df_difference,
        p = pchisq(deviance_difference, df = df_difference, lower.tail = FALSE),
        `ΔAIC (common − group-specific)` = AIC(null_model) - AIC(alternative_model),
        `ΔBIC (common − group-specific)` = GAIC(null_model, k = log(null_model$N)) - GAIC(alternative_model, k = log(alternative_model$N))
      )
    }
  }
}

parameter_table_labels <- c(
  mu = "μ: positive-score mean",
  nu = "ν: zero-score probability",
  sigma = "σ: positive-score scale"
)
table8_numeric <- bind_rows(interaction_rows) |>
  group_by(.data$comparison_id) |>
  mutate(`BH q` = p.adjust(.data$p, method = "BH")) |>
  ungroup() |>
  mutate(Parameter = parameter_table_labels[.data$parameter]) |>
  select(all_of(c(
    "Comparison", "Parameter", "Support", "Deviance difference",
    "Effective-df difference", "p", "BH q",
    "ΔAIC (common − group-specific)", "ΔBIC (common − group-specific)"
  )))
write.csv(table8_numeric, file.path(output_dirs$tables, "Supplementary_Table_8_group_shape_tests_numeric.csv"), row.names = FALSE)

format_probability <- function(x) {
  ifelse(x < 1e-16, "<1 × 10⁻¹⁶", ifelse(x < 0.001, formatC(x, format = "e", digits = 2), formatC(x, format = "f", digits = 3)))
}
table8 <- table8_numeric |>
  mutate(
    `Deviance difference` = round(.data$`Deviance difference`, 2),
    `Effective-df difference` = round(.data$`Effective-df difference`, 2),
    p = format_probability(.data$p),
    `BH q` = format_probability(.data$`BH q`),
    `ΔAIC (common − group-specific)` = round(.data$`ΔAIC (common − group-specific)`, 1),
    `ΔBIC (common − group-specific)` = round(.data$`ΔBIC (common − group-specific)`, 1)
  )
write.csv(table8, file.path(output_dirs$tables, "Supplementary_Table_8_group_shape_tests.csv"), row.names = FALSE)

publication_flextable <- function(data, left_columns = 1L) {
  flextable(data) |>
    theme_booktabs() |>
    font(fontname = "Times New Roman", part = "all") |>
    fontsize(size = 10, part = "all") |>
    bold(part = "header") |>
    align(align = "center", part = "all") |>
    align(j = left_columns, align = "left", part = "body") |>
    valign(valign = "center", part = "all") |>
    padding(padding = 1.5, part = "all") |>
    autofit()
}

add_document_text <- function(document, text, size = 10, bold = FALSE, space_after = 4) {
  body_add_fpar(
    document,
    fpar(
      ftext(text, fp_text(font.family = "Times New Roman", font.size = size, bold = bold)),
      fp_p = fp_par(padding.bottom = space_after)
    )
  )
}

new_landscape_document <- function() {
  body_set_default_section(
    read_docx(),
    value = prop_section(page_size = page_size(orient = "landscape"))
  )
}

family_code <- function(value) sub("^.*\\((NO|ZAGA|ZAIG|ZANBI|ZABB)\\)$", "\\1", value)
table6a_word <- table6a |> mutate(Family = family_code(.data$Family))
table6b_word <- table6b |> mutate(Family = family_code(.data$Family))
table6c_word <- table6c |>
  mutate(
    Family = family_code(.data$Family),
    `Expected n at/near 21*` = ifelse(
      as.numeric(.data$`Expected n at/near 21*`) < 0.05 & as.numeric(.data$`Expected n at/near 21*`) > 0,
      "<0.1",
      formatC(as.numeric(.data$`Expected n at/near 21*`), format = "f", digits = 1)
    ),
    `Expected n above 21` = ifelse(
      as.numeric(.data$`Expected n above 21`) < 0.05 & as.numeric(.data$`Expected n above 21`) > 0,
      "<0.1",
      formatC(as.numeric(.data$`Expected n above 21`), format = "f", digits = 1)
    )
  ) |>
  rename(`Expected n at/near 21` = `Expected n at/near 21*`)

table6_note <- paste(
  paste0(
    "All ", nrow(family_fit_status),
    " candidate models converged and used the intended observations and parameter equations."
  ),
  "AIC and BIC comparisons within a stratum are limited to ZAGA versus ZAIG and ZANBI versus ZABB. NO is shown as a Gaussian baseline but is not ranked against the zero-adjusted models because its density at zero and their point probability at zero have different likelihood scales. Values are not compared across strata.",
  "NO, Gaussian; ZAGA, zero-adjusted Gamma; ZAIG, zero-adjusted inverse Gaussian; ZANBI, zero-adjusted negative binomial type I; ZABB, zero-adjusted beta-binomial.",
  paste0("Residual summaries are medians across ", residual_replicates,
         " normalized quantile-residual evaluations; randomized draws were used at probability masses, whereas Gaussian residuals were deterministic."),
  "The maximum absolute central Q-Q deviation is evaluated over theoretical unit-normal quantiles from -2.5 to 2.5, and the mean upper-tail deviation over quantiles greater than 2.5. These summaries and the percentages outside the pointwise envelopes are descriptive rather than formal goodness-of-fit tests.",
  "In section c, expected n at/near 21 is n x P(20.5 < Y <= 21.5) for NO, ZAGA and ZAIG and n x P(Y = 21) for ZANBI and ZABB; expected n above 21 is n x P(Y > 21). Only ZABB is bounded at 21. The ZABB specification is a bounded-discrete sensitivity analysis and does not imply that GAD-7 arose from 21 exchangeable Bernoulli trials.",
  "Information-criterion differences are descriptive and are not statistical significance tests."
)

table6_doc <- new_landscape_document()
table6_doc <- add_document_text(table6_doc, "Supplementary Table 6 | Distribution-family sensitivity analysis across the six GAMLSS strata", 11, TRUE)
table6_doc <- add_document_text(table6_doc, "a  Continuous and mixed-distribution candidates", 10, TRUE)
table6_doc <- body_add_flextable(table6_doc, publication_flextable(table6a_word, c(1, 2)))
table6_doc <- body_add_break(table6_doc)
table6_doc <- add_document_text(table6_doc, "Supplementary Table 6 | Distribution-family sensitivity analysis across the six GAMLSS strata (continued)", 11, TRUE)
table6_doc <- add_document_text(table6_doc, "b  Discrete candidates", 10, TRUE)
table6_doc <- body_add_flextable(table6_doc, publication_flextable(table6b_word, c(1, 2)))
table6_doc <- body_add_break(table6_doc)
table6_doc <- add_document_text(table6_doc, "Supplementary Table 6 | Distribution-family sensitivity analysis across the six GAMLSS strata (continued)", 11, TRUE)
table6_doc <- add_document_text(table6_doc, "c  Upper-boundary behaviour across all candidate families", 10, TRUE)
table6_doc <- body_add_flextable(table6_doc, publication_flextable(table6c_word, c(1, 2)))
table6_doc <- add_document_text(table6_doc, paste0("Note. ", table6_note), 10, FALSE, 0)
write_publication_docx(table6_doc, file.path(output_dirs$tables, "Supplementary_Table_6.docx"))

table7_note <- paste(
  "All models used the ZAGA family and converged.",
  "In the overall strata, the \u03bc, \u03c3 and \u03bd equations each included penalized B-splines for kinship, offline-friend and online-friend support, with adjustment for sex and survey year.",
  "In the sex-stratified models, each equation included the same support terms and survey year. edf, effective degrees of freedom. AIC and BIC are not compared across strata because the fitted samples differ."
)
table7_doc <- new_landscape_document()
table7_doc <- add_document_text(table7_doc, "Supplementary Table 7 | Fit statistics for the six primary ZAGA models", 11, TRUE)
table7_doc <- body_add_flextable(table7_doc, publication_flextable(table7, c(1, 3)))
table7_doc <- add_document_text(table7_doc, paste0("Note. ", table7_note), 10, FALSE, 0)
write_publication_docx(table7_doc, file.path(output_dirs$tables, "Supplementary_Table_7.docx"))

table8_word <- table8 |>
  mutate(
    Parameter = recode(
      .data$Parameter,
      `μ: positive-score mean` = "\u03bc",
      `ν: zero-score probability` = "\u03bd",
      `σ: positive-score scale` = "\u03c3"
    ),
    Support = recode(.data$Support, `Offline friend` = "Offline-friend", `Online friend` = "Online-friend")
  ) |>
  select(all_of(c(
    "Comparison", "Parameter", "Support",
    "Deviance difference", "Effective-df difference", "p", "BH q",
    "\u0394AIC (common \u2212 group-specific)",
    "\u0394BIC (common \u2212 group-specific)"
  )))
names(table8_word) <- c(
  "Comparison", "Parameter", "Support source", "\u0394 penalized deviance",
  "\u0394edf", "p", "BH q", "\u0394AIC", "\u0394BIC"
)
table8_ft <- publication_flextable(table8_word, c(1, 2, 3)) |>
  merge_v(j = c("Comparison", "Parameter")) |>
  valign(j = c("Comparison", "Parameter"), valign = "center", part = "body")
table8_note <- paste(
  "\u03bc denotes the mean of the positive-score component, \u03bd the probability of GAD-7 = 0 and \u03c3 the distribution-specific scale/dispersion parameter of the positive-score component.",
  "For each comparison, support source and distributional parameter, the common-association model was compared with a model containing group-specific linear and quadratic support terms. Other support sources remained penalized B-splines, and the intended covariates were retained in all three parameter equations.",
  "\u0394 values are calculated as common-association minus group-specific; positive \u0394AIC or \u0394BIC therefore favours the group-specific specification on the corresponding information criterion. Approximate p values use differences in penalized deviance and effective degrees of freedom (\u0394edf), which may be non-integer.",
  "Benjamini-Hochberg correction was applied to the nine tests within each comparison set. These analyses compare overall linear-plus-quadratic association shapes; they are not exact tests of the plotted penalized-spline curves and do not test a three-way support-by-addiction-status-by-sex interaction.",
  "Information-criterion differences are descriptive and are not statistical significance tests."
)
table8_doc <- new_landscape_document()
table8_doc <- add_document_text(table8_doc, "Supplementary Table 8 | Approximate joint comparisons of support-association shapes in pooled ZAGA models", 11, TRUE)
table8_doc <- body_add_flextable(table8_doc, table8_ft)
table8_doc <- add_document_text(table8_doc, paste0("Note. ", table8_note), 10, FALSE, 0)
write_publication_docx(table8_doc, file.path(output_dirs$tables, "Supplementary_Table_8.docx"))

writexl::write_xlsx(
  list(
    Table_6A = table6a, Table_6B = table6b, Table_6C = table6c,
    Table_7 = table7, Table_8 = table8_numeric
  ),
  file.path(output_dirs$tables, "Supplementary_Tables_6_to_8.xlsx")
)
