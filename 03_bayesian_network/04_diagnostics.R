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

  sampler <- rstan::get_sampler_params(fit$fit, inc_warmup = FALSE)
  retained <- sum(vapply(sampler, nrow, integer(1)))
  max_depth <- unique(vapply(fit$fit@stan_args, function(x) as.numeric(x$max_depth), numeric(1)))
  if (length(max_depth) != 1L) stop("Chains used different maximum tree depths.", call. = FALSE)

  divergences <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1)))
  depth_hits <- sum(vapply(sampler, function(x) sum(x[, "treedepth__"] >= max_depth), numeric(1)))
  ebfmi <- vapply(sampler, function(x) {
    energy <- x[, "energy__"]
    mean(diff(energy)^2) / stats::var(energy)
  }, numeric(1))

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

  list(full = diagnostics, table = table, ebfmi = ebfmi)
}
