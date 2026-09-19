if (!exists("ZANBI_order_safe", mode = "function", inherits = TRUE)) {
  source("00_zanbi_order_safe.R", encoding = "UTF-8")
}

zanbi_preflight <- zanbi_order_safe_audit()
write.csv(
  zanbi_preflight,
  file.path(output_dirs$logs, "ZANBI_order_safe_preflight.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
if (!all(zanbi_preflight$passed)) {
  stop(
    "Order-safe ZANBI preflight failed: ",
    paste(zanbi_preflight$component[!zanbi_preflight$passed], collapse = ", ")
  )
}
message("Order-safe ZANBI preflight passed (", nrow(zanbi_preflight), " checks).")

family_specs <- list(
  NO = list(label = "Gaussian (NO)", set = "Continuous/mixed", ic_set = "Gaussian baseline", family = NO, response = "numeric", type = "Continuous", nu = FALSE),
  ZAGA = list(label = "Zero-adjusted Gamma (ZAGA)", set = "Continuous/mixed", ic_set = "Zero-adjusted continuous", family = ZAGA, response = "numeric", type = "Mixed", nu = TRUE),
  ZAIG = list(label = "Zero-adjusted inverse Gaussian (ZAIG)", set = "Continuous/mixed", ic_set = "Zero-adjusted continuous", family = ZAIG, response = "numeric", type = "Mixed", nu = TRUE),
  ZANBI = list(label = "Zero-adjusted negative binomial type I (ZANBI)", set = "Discrete", ic_set = "Discrete", family = ZANBI_order_safe, response = "count", type = "Discrete", nu = TRUE),
  ZABB = list(label = "Zero-adjusted beta-binomial (ZABB)", set = "Discrete", ic_set = "Discrete", family = ZABB, response = "bounded", type = "Discrete", nu = TRUE)
)

model_formulas <- function(spec, family_id) {
  rhs <- rhs_terms(spec$model_set)
  response <- if (identical(family_id, "ZABB")) {
    paste0("cbind(gad7_total, ", maximum_score, " - gad7_total)")
  } else {
    "gad7_total"
  }
  list(
    mu = as.formula(paste(response, "~", paste(rhs, collapse = " + "))),
    parameter = as.formula(paste("~", paste(rhs, collapse = " + ")))
  )
}

fit_control <- function() {
  gamlss.control(n.cyc = 150, trace = TRUE)
}

check_stage_fit <- function(model, spec, family_id, stage_label) {
  model_data <- model_data_for(spec)
  checks <- c(
    inherits_gamlss = inherits(model, "gamlss"),
    converged = isTRUE(model$converged),
    family = identical(as.character(model$family[[1L]]), family_id),
    sample_size = identical(as.integer(model$N), as.integer(nrow(model_data)))
  )
  if (identical(family_id, "ZABB")) {
    checks <- c(
      checks,
      bounded_at_21 = !is.null(model$bd) &&
        length(model$bd) == model$N &&
        all(is.finite(model$bd)) &&
        all(as.numeric(model$bd) == maximum_score)
    )
  }
  if (!all(checks)) {
    stop(
      "Staged fit failed verification: ", spec$model_id, "__", family_id,
      " | ", stage_label, " | failed checks: ",
      paste(names(checks)[!checks], collapse = ", ")
    )
  }
  message(
    "  Completed ", stage_label,
    " (iterations ", model$iter,
    "; Global Deviance ", round(model$G.deviance, 4), ")"
  )
  invisible(model)
}

fit_three_stage_model <- function(spec, family_id, model_data, formulas, family_spec) {
  key <- paste(spec$model_id, family_id, sep = "__")

  fit_stage <- function(sigma_formula, nu_formula, start_from = NULL, stage_label) {
    message("  ", key, " | ", stage_label)
    model <- if (is.null(start_from)) {
      gamlss(
        formula = formulas$mu,
        sigma.formula = sigma_formula,
        nu.formula = nu_formula,
        family = family_spec$family(),
        data = model_data,
        control = fit_control()
      )
    } else {
      gamlss(
        formula = formulas$mu,
        sigma.formula = sigma_formula,
        nu.formula = nu_formula,
        family = family_spec$family(),
        data = model_data,
        start.from = start_from,
        control = fit_control()
      )
    }
    if (identical(family_id, "ZANBI")) {
      attr(model, "zanbi_implementation") <- ZANBI_ORDER_SAFE_ID
      attr(model, "zanbi_preflight") <- zanbi_preflight
    }
    check_stage_fit(model, spec, family_id, stage_label)
    model
  }

  stage_1 <- fit_stage(
    sigma_formula = ~1,
    nu_formula = ~1,
    stage_label = "stage 1/3: full mu; intercept-only sigma and nu"
  )
  stage_2 <- fit_stage(
    sigma_formula = formulas$parameter,
    nu_formula = ~1,
    start_from = stage_1,
    stage_label = "stage 2/3: full mu and sigma; intercept-only nu"
  )
  stage_3 <- fit_stage(
    sigma_formula = formulas$parameter,
    nu_formula = formulas$parameter,
    start_from = stage_2,
    stage_label = "stage 3/3: full mu, sigma and nu"
  )

  attr(stage_3, "fit_strategy") <- "three-stage hierarchical warm start"
  attr(stage_3, "staged_start_audit") <- data.frame(
    stage = c("mu_full", "mu_sigma_full", "mu_sigma_nu_full"),
    iterations = as.integer(c(stage_1$iter, stage_2$iter, stage_3$iter)),
    converged = c(
      isTRUE(stage_1$converged),
      isTRUE(stage_2$converged),
      isTRUE(stage_3$converged)
    ),
    global_deviance = as.numeric(c(
      stage_1$G.deviance,
      stage_2$G.deviance,
      stage_3$G.deviance
    )),
    stringsAsFactors = FALSE
  )
  stage_3
}

fit_candidate_model <- function(spec, family_id) {
  model_data <- model_data_for(spec)
  formulas <- model_formulas(spec, family_id)
  family_spec <- family_specs[[family_id]]

  if (family_id %in% c("ZANBI", "ZABB")) {
    return(fit_three_stage_model(spec, family_id, model_data, formulas, family_spec))
  }

  arguments <- list(
    formula = formulas$mu,
    sigma.formula = formulas$parameter,
    family = family_spec$family(),
    data = model_data,
    control = fit_control()
  )
  if (family_spec$nu) arguments$nu.formula <- formulas$parameter
  model <- do.call(gamlss, arguments)
  attr(model, "fit_strategy") <- "direct final-formula fit"
  model
}

verify_model <- function(model, spec, family_id) {
  formulas <- model_formulas(spec, family_id)
  expected_mu <- gsub("\\s+", " ", formula_text(formulas$mu))
  expected_parameter <- formula_text(formulas$parameter)
  parameters <- if (family_specs[[family_id]]$nu) c("mu", "sigma", "nu") else c("mu", "sigma")
  observed_formulas <- setNames(
    vapply(parameters, function(parameter) {
      gsub("\\s+", " ", formula_text(model[[paste0(parameter, ".formula")]]))
    }, character(1)),
    parameters
  )
  mu_formula_ok <- identical(observed_formulas[["mu"]], expected_mu)
  parameter_formula_ok <- all(vapply(
    setdiff(parameters, "mu"),
    function(parameter) identical(observed_formulas[[parameter]], gsub("\\s+", " ", expected_parameter)),
    logical(1)
  ))
  formula_ok <- mu_formula_ok && parameter_formula_ok
  family_ok <- identical(as.character(model$family[[1L]]), family_id)
  zanbi_implementation <- attr(model, "zanbi_implementation", exact = TRUE)
  zanbi_implementation_ok <- if (identical(family_id, "ZANBI")) {
    identical(zanbi_implementation, ZANBI_ORDER_SAFE_ID)
  } else {
    TRUE
  }
  type_ok <- identical(as.character(model$type), family_specs[[family_id]]$type)
  boundary_ok <- if (identical(family_id, "ZABB")) {
    !is.null(model$bd) && all(is.finite(model$bd)) && all(as.numeric(model$bd) == maximum_score)
  } else {
    TRUE
  }
  model_data <- model_data_for(spec)
  calculated_aic <- as.numeric(model$G.deviance) + 2 * as.numeric(model$df.fit)
  calculated_bic <- as.numeric(model$G.deviance) + log(as.numeric(model$N)) * as.numeric(model$df.fit)
  package_aic <- as.numeric(AIC(model))
  package_bic <- as.numeric(GAIC(model, k = log(as.numeric(model$N))))
  aic_ok <- isTRUE(all.equal(calculated_aic, package_aic, tolerance = 1e-8))
  bic_ok <- isTRUE(all.equal(calculated_bic, package_bic, tolerance = 1e-8))
  fit_strategy <- attr(model, "fit_strategy", exact = TRUE)
  if (is.null(fit_strategy)) fit_strategy <- "cache without fit-strategy metadata"
  staged_start_audit <- attr(model, "staged_start_audit", exact = TRUE)
  staged_iterations <- if (is.null(staged_start_audit)) {
    NA_character_
  } else {
    paste(
      paste(staged_start_audit$stage, staged_start_audit$iterations, sep = "="),
      collapse = "; "
    )
  }
  staged_global_deviance <- if (is.null(staged_start_audit)) {
    NA_character_
  } else {
    paste(
      paste(
        staged_start_audit$stage,
        format(round(staged_start_audit$global_deviance, 4), trim = TRUE),
        sep = "="
      ),
      collapse = "; "
    )
  }
  tibble(
    model_id = spec$model_id,
    stratum = spec$stratum,
    data_fingerprint = data_fingerprint,
    family = family_id,
    family_label = family_specs[[family_id]]$label,
    comparison_set = family_specs[[family_id]]$set,
    likelihood_set = family_specs[[family_id]]$ic_set,
    n = as.integer(model$N),
    intended_n = nrow(model_data),
    converged = isTRUE(model$converged),
    family_ok = family_ok,
    zanbi_implementation = if (is.null(zanbi_implementation)) NA_character_ else zanbi_implementation,
    zanbi_implementation_ok = zanbi_implementation_ok,
    response_type = as.character(model$type),
    response_type_ok = type_ok,
    expected_mu_formula = expected_mu,
    observed_mu_formula = observed_formulas[["mu"]],
    mu_formula_ok = mu_formula_ok,
    expected_parameter_formula = gsub("\\s+", " ", expected_parameter),
    observed_sigma_formula = observed_formulas[["sigma"]],
    observed_nu_formula = if ("nu" %in% parameters) observed_formulas[["nu"]] else NA_character_,
    parameter_formula_ok = parameter_formula_ok,
    formula_ok = formula_ok,
    response_formula_ok = mu_formula_ok,
    fit_strategy = fit_strategy,
    staged_iterations = staged_iterations,
    staged_global_deviance = staged_global_deviance,
    boundary_ok = boundary_ok,
    sample_ok = as.integer(model$N) == nrow(model_data),
    iterations = as.integer(model$iter),
    effective_df = as.numeric(model$df.fit),
    global_deviance = as.numeric(model$G.deviance),
    log_likelihood = -as.numeric(model$G.deviance) / 2,
    AIC = calculated_aic,
    package_AIC = package_aic,
    AIC_ok = aic_ok,
    BIC = calculated_bic,
    package_BIC = package_bic,
    BIC_ok = bic_ok
  ) |>
    mutate(
      all_checks_pass = .data$converged & .data$family_ok & .data$zanbi_implementation_ok & .data$response_type_ok & .data$formula_ok &
        .data$boundary_ok & .data$sample_ok & .data$AIC_ok & .data$BIC_ok
    )
}

candidate_models <- list()
status_rows <- list()
for (spec_index in seq_len(nrow(analysis_specs))) {
  spec <- analysis_specs[spec_index, ]
  for (family_id in names(family_specs)) {
    key <- paste(spec$model_id, family_id, sep = "__")
    cache_file <- file.path(output_dirs$models, paste0(key, "__", model_cache_tag, ".rds"))
    fit_started <- Sys.time()
    cache_exists <- file.exists(cache_file)
    cache_reused <- FALSE
    model <- if (cache_exists) {
      message("Checking cached model: ", key)
      cached_model <- readRDS(cache_file)
      zanbi_cache_ok <- !identical(family_id, "ZANBI") ||
        identical(attr(cached_model, "zanbi_implementation", exact = TRUE), ZANBI_ORDER_SAFE_ID)
      if (zanbi_cache_ok) {
        cache_reused <- TRUE
        cached_model
      } else {
        message("Rejecting incompatible ZANBI cache and refitting: ", key)
        fit_candidate_model(spec, family_id)
      }
    } else {
      message("Fitting model: ", key)
      fit_candidate_model(spec, family_id)
    }
    status <- verify_model(model, spec, family_id)
    if (!status$all_checks_pass) stop("Model verification failed: ", key)
    if (!cache_reused) {
      if (cache_exists) {
        backup_file <- sub("\\.rds$", "__incompatible_backup.rds", cache_file)
        if (!file.exists(backup_file) && !file.copy(cache_file, backup_file, overwrite = FALSE)) {
          stop("Unable to preserve the incompatible cache: ", cache_file)
        }
      }
      saveRDS(model, cache_file)
    }
    message(
      "Verified model: ", key,
      " (iterations ", status$iterations,
      "; elapsed ", round(as.numeric(difftime(Sys.time(), fit_started, units = "mins")), 1), " min)"
    )
    candidate_models[[key]] <- model
    status_rows[[key]] <- status
  }
}

family_fit_status <- bind_rows(status_rows)
write.csv(family_fit_status, file.path(output_dirs$tables, "Family_fit_status.csv"), row.names = FALSE, fileEncoding = "UTF-8")

information_criteria <- family_fit_status |>
  group_by(.data$stratum, .data$likelihood_set) |>
  mutate(
    delta_AIC_within_set = if (n() > 1L) .data$AIC - min(.data$AIC) else NA_real_,
    delta_BIC_within_set = if (n() > 1L) .data$BIC - min(.data$BIC) else NA_real_
  ) |>
  ungroup()
write.csv(information_criteria, file.path(output_dirs$tables, "Family_information_criteria.csv"), row.names = FALSE, fileEncoding = "UTF-8")

primary_models <- setNames(
  lapply(analysis_specs$model_id, function(id) candidate_models[[paste(id, "ZAGA", sep = "__")]]),
  analysis_specs$model_id
)
