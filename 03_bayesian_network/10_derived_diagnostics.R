derived_diagnostics <- function(matrix, fit, type) {
  dims <- dim(posterior::as_draws_array(fit))
  if (nrow(matrix) != dims[1] * dims[2]) stop("Derived diagnostics require all draws in original chain order.")
  if (is.null(colnames(matrix))) colnames(matrix) <- paste0(type, "_", seq_len(ncol(matrix)))
  draws <- posterior::as_draws_array(array(as.numeric(matrix), c(dims[1], dims[2], ncol(matrix)),
                                           dimnames = list(NULL, NULL, colnames(matrix))))
  sm <- posterior::summarise_draws(draws, "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean",
    q025 = function(x) unname(stats::quantile(x, probs = .025)),
    q975 = function(x) unname(stats::quantile(x, probs = .975)),
    mcse_q025 = function(x) unname(posterior::mcse_quantile(x, probs = .025)),
    mcse_q975 = function(x) unname(posterior::mcse_quantile(x, probs = .975)))
  sm$type <- type
  as.data.frame(sm)
}
write_derived_diagnostics <- function(fit, arrays, metrics, directory) {
  edge <- do.call(cbind, lapply(names(arrays), function(g) {
    x <- matrix(arrays[[g]], nrow = dim(arrays[[g]])[1])
    colnames(x) <- paste(g, rep(network_nodes, times = length(network_nodes)),
                         rep(network_nodes, each = length(network_nodes)), sep = "|")
    x
  }))
  group_metrics <- matrix(metrics, nrow = dim(metrics)[1])
  colnames(group_metrics) <- paste(rep(dimnames(metrics)$group, each = dim(metrics)[2]),
                                   rep(dimnames(metrics)$metric, times = dim(metrics)[3]), sep = "|")
  female <- metrics[, , "addicted_female"] - metrics[, , "nonaddicted_female"]
  male <- metrics[, , "addicted_male"] - metrics[, , "nonaddicted_male"]
  contrasts <- cbind((female + male) / 2, female, male, male - female)
  colnames(contrasts) <- paste(rep(c("overall", "women", "men", "sex_difference"), each = ncol(female)),
                              rep(colnames(female), 4), sep = "|")
  diagnostics <- rbind(derived_diagnostics(edge, fit, "standardized_edges"),
                        derived_diagnostics(group_metrics, fit, "group_metrics"),
                        derived_diagnostics(contrasts, fit, "metric_contrasts"))
  write_csv(diagnostics, file.path(directory, "derived_mcmc_diagnostics.csv"))
  summaries <- lapply(split(diagnostics, diagnostics$type), function(x) data.frame(
    diagnostic = x$type[1], parameters_or_transitions = nrow(x),
    result = sprintf("max R-hat %.4f; min bulk ESS %.0f; min tail ESS %.0f", max(x$rhat), min(x$ess_bulk), min(x$ess_tail))))
  diagnostics_primary <- model_diagnostics(fit)
  diagnostics_primary$table <- rbind(diagnostics_primary$table, do.call(rbind, summaries))
  flags <- diagnostics[diagnostics$type == "metric_contrasts" &
    (abs(diagnostics$q025) <= 2 * diagnostics$mcse_q025 |
     abs(diagnostics$q975) <= 2 * diagnostics$mcse_q975), , drop = FALSE]
  write_csv(flags, file.path(directory, "interval_endpoint_mcse_flags.csv"))
  diagnostics_primary$table <- rbind(diagnostics_primary$table,
    data.frame(diagnostic = "Interval endpoints near zero", parameters_or_transitions = nrow(flags),
               result = "within two endpoint MCSE of zero; consult interval_endpoint_mcse_flags.csv"))
  write_diagnostic_table(diagnostics_primary, directory)
}
