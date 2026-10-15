model_diagnostics <- function(fit) {
  draws <- posterior::as_draws_array(fit)
  rhat <- posterior::rhat
  ess_bulk <- posterior::ess_bulk
  ess_tail <- posterior::ess_tail
  diagnostics <- posterior::summarise_draws(
    draws,
    rhat,
    ess_bulk,
    ess_tail
  )
  sampled <- diagnostics[diagnostics$variable != "lp__", , drop = FALSE]
  regression <- diagnostics[grepl("^b_", diagnostics$variable), , drop = FALSE]
  residual_sd <- sampled[grepl("^sigma_", sampled$variable), , drop = FALSE]
  residual_correlations <- sampled[grepl("^rescor__", sampled$variable), , drop = FALSE]
  other <- sampled[!sampled$variable %in% c(regression$variable, residual_sd$variable, residual_correlations$variable), , drop = FALSE]

  sampler <- brms::nuts_params(fit)
  required_sampler_parameters <- c("divergent__", "treedepth__", "energy__")
  if (!all(required_sampler_parameters %in% sampler$Parameter)) {
    stop("The fitted model is missing NUTS sampler diagnostics.", call. = FALSE)
  }
  max_depth <- attr(fit, "analysis_spec")$max_treedepth
  if (is.null(max_depth)) max_depth <- fit$stan_args$control$max_treedepth
  if (is.null(max_depth) && inherits(fit$fit, "stanfit")) {
    max_depth <- unique(vapply(
      fit$fit@stan_args, function(x) {
        z <- x$control$max_treedepth
        if (is.null(z)) z <- x$max_depth
        if (is.null(z)) NA_real_ else as.numeric(z)
      }, numeric(1)
    ))
  }
  if (length(max_depth) != 1L || !is.finite(max_depth)) {
    stop("The fitted model does not record one maximum tree depth.", call. = FALSE)
  }

  divergent <- sampler[sampler$Parameter == "divergent__", ]
  depth <- sampler[sampler$Parameter == "treedepth__", ]
  energy <- sampler[sampler$Parameter == "energy__", ]
  retained <- nrow(divergent)
  if (nrow(depth) != retained || nrow(energy) != retained) {
    stop("NUTS diagnostic lengths differ.", call. = FALSE)
  }
  divergences <- sum(divergent$Value)
  depth_hits <- sum(depth$Value >= max_depth)
  energy <- energy[order(energy$Chain, energy$Iteration), ]
  energy_by_chain <- split(energy$Value, energy$Chain)
  ebfmi <- vapply(energy_by_chain, function(energy_values) {
    if (length(energy_values) < 3L || stats::var(energy_values) <= 0) {
      return(NA_real_)
    }
    energy <- energy_values
    mean(diff(energy)^2) / stats::var(energy)
  }, numeric(1))
  if (any(!is.finite(ebfmi))) stop("E-BFMI could not be calculated.", call. = FALSE)

  table <- data.frame(
    diagnostic = c(
      "Regression coefficients",
      "All sampled parameters excluding lp__",
      "Divergent transitions",
      paste0("Transitions at maximum tree depth ", max_depth),
      "Minimum E-BFMI across chains"
    ),
    parameters_or_transitions = c(nrow(regression), nrow(sampled), retained, retained, length(ebfmi)),
    result = c(
      sprintf(
        "max Rhat %.4f; min bulk ESS %.0f; min tail ESS %.0f",
        max(regression$rhat, na.rm = TRUE),
        min(regression$ess_bulk, na.rm = TRUE),
        min(regression$ess_tail, na.rm = TRUE)
      ),
      sprintf(
        "max Rhat %.4f; min bulk ESS %.0f; min tail ESS %.0f",
        max(sampled$rhat, na.rm = TRUE),
        min(sampled$ess_bulk, na.rm = TRUE),
        min(sampled$ess_tail, na.rm = TRUE)
      ),
      sprintf("%d (%.1f%%)", divergences, 100 * divergences / retained),
      sprintf("%d (%.1f%%)", depth_hits, 100 * depth_hits / retained),
      sprintf("%.3f", min(ebfmi))
    ),
    stringsAsFactors = FALSE
  )

  groups <- list("Residual SD parameters" = residual_sd, "Residual correlations" = residual_correlations,
                 "Other sampled parameters" = other)
  additional <- do.call(rbind, lapply(names(groups), function(label) {
    x <- groups[[label]]
    data.frame(diagnostic = label, parameters_or_transitions = nrow(x),
      result = if (nrow(x)) sprintf("max Rhat %.4f; min bulk ESS %.0f; min tail ESS %.0f",
                                   max(x$rhat), min(x$ess_bulk), min(x$ess_tail)) else "No parameters in this category")
  }))
  table <- rbind(table[1L, ], additional, table[-1L, ])
  chain_table <- data.frame(diagnostic = "E-BFMI by chain", parameters_or_transitions = length(ebfmi),
                            result = paste(sprintf("%.3f", ebfmi), collapse = ", "))
  table <- rbind(table, chain_table)

  list(full = diagnostics, table = table, ebfmi = ebfmi)
}
