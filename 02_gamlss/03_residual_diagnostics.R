fitted_arguments <- function(model) {
  arguments <- lapply(model$parameters, function(parameter) fitted(model, parameter))
  names(arguments) <- model$parameters
  if (!is.null(model$bd)) arguments$bd <- model$bd
  arguments
}

quantile_residual_intervals <- function(model, family_id) {
  pfun <- distribution_function_for(family_id, "p")
  arguments <- fitted_arguments(model)
  cdf <- function(q) do.call(pfun, c(list(q = q), arguments))

  if (identical(model$type, "Continuous")) {
    probability <- cdf(model$y)
    return(list(type = "Continuous", lower = probability, upper = probability))
  }

  if (identical(model$type, "Mixed")) {
    probability <- cdf(model$y)
    return(list(
      type = "Mixed", lower = probability, upper = probability,
      is_zero = model$y == 0, zero_probability = fitted(model, "nu")
    ))
  }

  if (identical(model$type, "Discrete")) {
    lower <- cdf(ifelse(model$y == 0, 0, model$y - 1))
    lower[model$y == 0] <- 0
    upper <- cdf(model$y)
    invalid <- which(!is.finite(lower) | !is.finite(upper) | lower < 0 | upper > 1 | lower > upper)
    if (length(invalid)) {
      dfun <- distribution_function_for(family_id, "d")
      for (index in invalid) {
        response <- as.integer(model$y[index])
        scalar_arguments <- lapply(arguments, function(value) value[index])
        mass <- do.call(dfun, c(list(x = 0:response), scalar_arguments))
        lower[index] <- if (response == 0L) 0 else sum(mass[-length(mass)])
        upper[index] <- sum(mass)
      }
    }
    return(list(type = "Discrete", lower = lower, upper = upper, fallback_index = invalid))
  }

  stop("Unsupported GAMLSS response type: ", model$type)
}

draw_quantile_residuals <- function(intervals, seed) {
  set.seed(seed)
  probability <- if (identical(intervals$type, "Mixed")) {
    ifelse(
      intervals$is_zero,
      runif(length(intervals$upper), min = 0, max = intervals$zero_probability),
      intervals$upper
    )
  } else if (identical(intervals$type, "Discrete")) {
    runif(length(intervals$lower), min = intervals$lower, max = intervals$upper)
  } else {
    intervals$upper
  }
  probability <- pmin(pmax(probability, .Machine$double.eps), 1 - .Machine$double.eps)
  qnorm(probability)
}

worm_components <- function(residual_values) {
  finite <- is.finite(residual_values)
  residual <- sort(residual_values[finite])
  n <- length(residual)
  if (n < 20L) stop("Too few finite residuals for diagnostics")
  index <- seq_len(n)
  theoretical <- qnorm(ppoints(n))
  deviation <- residual - theoretical
  envelope_low <- qnorm(qbeta(0.025, index, n + 1L - index)) - theoretical
  envelope_high <- qnorm(qbeta(0.975, index, n + 1L - index)) - theoretical
  fitted_curve <- as.numeric(predict(lm(deviation ~ theoretical + I(theoretical^2) + I(theoretical^3))))
  central <- abs(theoretical) <= 2.5
  upper_tail <- theoretical > 2.5
  bins <- cut(theoretical[central], seq(-2.5, 2.5, length.out = 21L), include.lowest = TRUE)
  binned <- as.numeric(tapply(deviation[central], bins, mean))
  binned <- binned[is.finite(binned)]
  signs <- sign(binned)
  signs <- signs[signs != 0]
  list(
    data = tibble(
      theoretical = theoretical, residual = residual, deviation = deviation,
      envelope_low = envelope_low, envelope_high = envelope_high, fitted_curve = fitted_curve
    ),
    metrics = tibble(
      residual_n_total = length(residual_values),
      residual_n_finite = n,
      residual_n_nonfinite = sum(!finite),
      residual_mean = mean(residual),
      residual_sd = sd(residual),
      maximum_absolute_deviation = max(abs(deviation)),
      central_maximum_absolute_deviation = max(abs(deviation[central])),
      percentage_outside_pointwise_envelope = 100 * mean(deviation < envelope_low | deviation > envelope_high),
      upper_tail_mean_deviation_q_gt_2_5 = if (any(upper_tail)) mean(deviation[upper_tail]) else NA_real_,
      upper_tail_minimum_deviation_q_gt_2_5 = if (any(upper_tail)) min(deviation[upper_tail]) else NA_real_,
      central_wave_range = if (length(binned)) diff(range(binned)) else NA_real_,
      central_wave_sign_changes = if (length(signs) > 1L) sum(diff(signs) != 0) else 0L
    )
  )
}

residual_metric_rows <- list()
representative_worm_rows <- list()
fallback_rows <- list()
for (family_id in names(family_specs)) {
  for (group_index in seq_len(nrow(analysis_specs))) {
    spec <- analysis_specs[group_index, ]
    key <- paste(spec$model_id, family_id, sep = "__")
    intervals <- quantile_residual_intervals(candidate_models[[key]], family_id)
    if (length(intervals$fallback_index)) {
      fallback_rows[[key]] <- tibble(family = family_id, model_id = spec$model_id, observation_index = intervals$fallback_index)
    }
    for (replicate_index in seq_len(residual_replicates)) {
      seed <- residual_seed_base + group_index * 1000L + replicate_index
      worm <- worm_components(draw_quantile_residuals(intervals, seed))
      row_key <- paste(key, replicate_index, sep = "__")
      residual_metric_rows[[row_key]] <- worm$metrics |>
        mutate(
          family = family_id, family_label = family_specs[[family_id]]$label,
          comparison_set = family_specs[[family_id]]$set, model_id = spec$model_id,
          stratum = spec$stratum, replicate = replicate_index, seed = seed, .before = 1
        )
      if (replicate_index == 1L) {
        keep <- unique(round(seq(1, nrow(worm$data), length.out = min(1200L, nrow(worm$data)))))
        representative_worm_rows[[key]] <- worm$data[keep, ] |>
          mutate(family = family_id, model_id = spec$model_id, stratum = spec$stratum, .before = 1)
      }
    }
  }
}

residual_replicate_summary <- bind_rows(residual_metric_rows)
residual_diagnostic_summary <- residual_replicate_summary |>
  group_by(.data$family, .data$family_label, .data$comparison_set, .data$model_id, .data$stratum) |>
  summarise(
    residual_draws = n(),
    representative_seed = first(.data$seed),
    maximum_nonfinite_residuals = max(.data$residual_n_nonfinite),
    across(
      c(
        residual_mean, residual_sd, maximum_absolute_deviation,
        central_maximum_absolute_deviation, percentage_outside_pointwise_envelope,
        upper_tail_mean_deviation_q_gt_2_5, upper_tail_minimum_deviation_q_gt_2_5,
        central_wave_range, central_wave_sign_changes
      ),
      list(median = ~ median(.x, na.rm = TRUE), q025 = ~ quantile(.x, 0.025, na.rm = TRUE), q975 = ~ quantile(.x, 0.975, na.rm = TRUE))
    ),
    .groups = "drop"
  )
representative_worm <- bind_rows(representative_worm_rows)

write.csv(residual_replicate_summary, file.path(output_dirs$tables, "Family_residual_diagnostic_replicates_50.csv"), row.names = FALSE)
write.csv(residual_diagnostic_summary, file.path(output_dirs$tables, "Family_residual_diagnostic_summary.csv"), row.names = FALSE)
write.csv(bind_rows(fallback_rows), file.path(output_dirs$tables, "Residual_CDF_fallback_cases.csv"), row.names = FALSE)
saveRDS(representative_worm, file.path(output_dirs$source_data, "Family_representative_worm_data.rds"))

boundary_rows <- list()
for (family_id in names(family_specs)) {
  for (spec_index in seq_len(nrow(analysis_specs))) {
    spec <- analysis_specs[spec_index, ]
    model <- candidate_models[[paste(spec$model_id, family_id, sep = "__")]]
    arguments <- fitted_arguments(model)
    pfun <- distribution_function_for(family_id, "p")
    probability_above_21 <- pmax(
      0,
      1 - do.call(pfun, c(list(q = rep(maximum_score, model$N)), arguments))
    )
    if (identical(model$type, "Discrete")) {
      dfun <- distribution_function_for(family_id, "d")
      score21_probability <- do.call(
        dfun,
        c(list(x = rep(maximum_score, model$N)), arguments)
      )
      score21_definition <- "P(Y = 21)"
    } else {
      cdf_low <- do.call(
        pfun,
        c(list(q = rep(maximum_score - 0.5, model$N)), arguments)
      )
      cdf_high <- do.call(
        pfun,
        c(list(q = rep(maximum_score + 0.5, model$N)), arguments)
      )
      score21_probability <- pmax(0, cdf_high - cdf_low)
      score21_definition <- "P(20.5 < Y <= 21.5)"
    }
    boundary_rows[[paste(family_id, spec$model_id)]] <- tibble(
      family = family_id,
      family_label = family_specs[[family_id]]$label,
      comparison_set = family_specs[[family_id]]$set,
      model_id = spec$model_id,
      stratum = spec$stratum,
      n = model$N,
      observed_n_at_21 = sum(model$y == maximum_score),
      score21_bin_definition = score21_definition,
      expected_n_in_score21_bin = sum(score21_probability),
      expected_n_above_21 = sum(probability_above_21),
      support_bounded_at_21 = identical(family_id, "ZABB")
    )
  }
}
boundary_summary <- bind_rows(boundary_rows)
write.csv(
  boundary_summary,
  file.path(output_dirs$tables, "Boundary_score_21_all_families.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
