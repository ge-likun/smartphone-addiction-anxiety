if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.UTF-8"))
}
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- if (length(args)) dirname(dirname(normalizePath(sub("^--file=", "", args[1]), winslash = "/"))) else getwd()
source(file.path(root, "03_bayesian_network", "06_metrics.R"))
network_nodes <- c("kin", "offline", "online", paste0("gad", 1:7), paste0("mpiq", 1:8))
network_communities <- list(Support = network_nodes[1:3], Anxiety = network_nodes[4:10], Smartphone_addiction = network_nodes[11:18])
impulse_steps <- 15L
W <- matrix(0, length(network_nodes), length(network_nodes), dimnames = list(network_nodes, network_nodes))
W["kin", "gad1"] <- -0.4
W["gad1", "gad1"] <- 0.5
m <- network_metrics(W)
stopifnot(abs(m[["Spectral radius: anxiety"]] - .5) < 1e-12)
stopifnot(abs(m[["Outgoing sum: kin to anxiety"]] + .4) < 1e-12)
stopifnot(abs(m[["Bridge strength: kin"]] - .4) < 1e-12)
stopifnot(abs(m[["Cumulative response: kin to anxiety"]] - (-.4 * sum(.5^(0:14)))) < 1e-12)
groups <- c("nonaddicted_female", "addicted_female", "nonaddicted_male", "addicted_male")
for (draws in c(9L, 17L)) {
  arrays <- setNames(lapply(seq_along(groups), function(i) {
    a <- array(rep(as.vector(W), each = draws), c(draws,nrow(W),ncol(W)), dimnames=list(NULL, network_nodes, network_nodes))
    a * i
  }), groups)
  md <- calculate_metric_draws(arrays)
  sm <- summarize_metric_contrasts(md)
  own <- md[, , "addicted_female"] - md[, , "nonaddicted_female"]
  male <- md[, , "addicted_male"] - md[, , "nonaddicted_male"]
  for (metric in dimnames(md)$metric) {
    row <- sm[sm$metric == metric & sm$contrast == "Addicted minus non-addicted, averaged across sexes", ]
    stopifnot(abs(row$posterior_mean - mean((own[,metric] + male[,metric])/2)) < 1e-12)
    stopifnot(row$n_draws == draws)
  }
}
p <- posterior_summary(c(rep(-1,95),rep(1,5)))
stopifnot(p[["probability_negative"]] == .95, p[["probability_of_direction"]] == .95)
cat("Metric orientation, propagation signs, contrasts and draw-count checks passed.\n")
