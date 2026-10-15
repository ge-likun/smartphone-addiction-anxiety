ordinal_context <- function(data = read_network_data()) {
  root <- Sys.getenv("PROJECT_ROOT", normalizePath("..", winslash = "/"))
  source(file.path(root, "R", "project.R"), encoding = "UTF-8")
  cfg <- read_project_config(root, Sys.getenv("PROJECT_CONFIG_FILE", "local_config.R"))$network
  if (nzchar(Sys.getenv("CMDSTAN"))) cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN"))
  nodes <- network_nodes
  means <- setNames(vapply(paste0(nodes, "_t1"), function(x) mean(data[[x]]), numeric(1)), nodes)
  sds <- setNames(vapply(paste0(nodes, "_t1"), function(x) sd(data[[x]]), numeric(1)), nodes)
  if (any(!is.finite(sds) | sds <= 0)) stop("Baseline nodes must have positive finite standard deviations.")
  raw_data <- data
  for (node in nodes) {
    data[[paste0("z_", node, "_t1")]] <- (data[[paste0(node, "_t1")]] - means[[node]]) / sds[[node]]
    values <- if (grepl("^gad", node)) 1:4 else if (grepl("^mpiq", node)) 1:7 else 0:15
    data[[paste0(node, "_t2")]] <- ordered(data[[paste0(node, "_t2")]], levels = values)
  }
  groups <- list(nonaddicted_female = which(raw_data$addiction == "Non-Addiction" & raw_data$gender == "Female"),
                 addicted_female = which(raw_data$addiction == "Addiction" & raw_data$gender == "Female"),
                 nonaddicted_male = which(raw_data$addiction == "Non-Addiction" & raw_data$gender == "Male"),
                 addicted_male = which(raw_data$addiction == "Addiction" & raw_data$gender == "Male"))
  directory <- file.path(output_dir, "ordinal")
  for (p in c(directory, file.path(directory, "models"))) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  list(data = data, raw = raw_data, nodes = nodes, groups = groups, means = means, sds = sds,
       scales = common_node_scales(raw_data), cfg = cfg, versions = software_versions(),
       hash = object_hash(raw_data), directory = directory)
}
ordinal_formula <- function(target, ctx) {
  as.formula(paste0(target, "_t2 ~ addiction * gender * (",
                    paste(paste0("z_", ctx$nodes, "_t1"), collapse = " + "), ")"))
}
ordinal_spec <- function(target, ctx) {
  c(list(data_hash = ctx$hash, target = target, n = nrow(ctx$data),
         formula = paste(deparse(ordinal_formula(target, ctx)), collapse = ""),
         means = ctx$means, sds = ctx$sds, versions = ctx$versions,
         family = "cumulative", link = "probit", coefficient_prior = "normal(0,0.5)",
         threshold_prior = "normal(0,3)", init = 0),
    list(chains = ctx$cfg$chains, iter = ctx$cfg$iter, warmup = ctx$cfg$warmup,
         seed = ctx$cfg$ordinal_seed + match(target, ctx$nodes),
         adapt_delta = ctx$cfg$ordinal_adapt_delta, max_treedepth = ctx$cfg$ordinal_max_treedepth))
}
sampling_summary <- function(fit, depth) {
  draws <- posterior::as_draws_array(fit)
  sm <- posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail")
  np <- brms::nuts_params(fit)
  chains <- lapply(sort(unique(np$Chain)), function(chain) {
    z <- np[np$Chain == chain, , drop = FALSE]
    energy <- z$Value[z$Parameter == "energy__"]
    data.frame(chain = chain, divergences = sum(z$Value[z$Parameter == "divergent__"]),
                depth_hits = sum(z$Value[z$Parameter == "treedepth__"] >= depth),
                E_BFMI = mean(diff(energy)^2) / var(energy))
  })
  ch <- do.call(rbind, chains)
  result <- data.frame(parameters = nrow(sm), draws = posterior::ndraws(draws),
                       max_rhat = max(sm$rhat, na.rm = TRUE), min_bulk = min(sm$ess_bulk, na.rm = TRUE),
                       min_tail = min(sm$ess_tail, na.rm = TRUE), divergences = sum(ch$divergences),
                       depth_hits = sum(ch$depth_hits), min_E_BFMI = min(ch$E_BFMI))
  result$pass <- all(is.finite(unlist(result))) && result$max_rhat <= 1.01 &&
    result$min_bulk >= 400 && result$min_tail >= 400 && result$divergences == 0 && result$min_E_BFMI >= 0.3
  list(summary = result, parameters = as.data.frame(sm), chains = ch)
}
verify_ordinal_bundle <- function(bundle, target, ctx) {
  spec <- ordinal_spec(target, ctx)
  if (!is.list(bundle) || !inherits(bundle$fit, "brmsfit") || !identical(bundle$spec, spec))
    stop("Ordinal cache specification mismatch: ", target, ". Use a separate results directory for a new specification.")
  fit <- bundle$fit
  columns <- c("addiction", "gender", paste0("z_", ctx$nodes, "_t1"), paste0(target, "_t2"))
  if (!identical(model_data_signature(fit$data, columns), model_data_signature(ctx$data, columns)) ||
      !identical(fit$family$family, "cumulative") || !identical(fit$family$link, "probit") ||
      !identical(attr(terms(fit$formula$formula), "term.labels"), attr(terms(ordinal_formula(target, ctx)), "term.labels")))
    stop("Cached ordinal model differs from intended data or formula: ", target)
  expected <- ctx$cfg$chains * (ctx$cfg$iter - ctx$cfg$warmup)
  if (posterior::ndraws(posterior::as_draws_array(fit)) != expected) stop("Incomplete posterior cache: ", target)
  prior_text <- gsub("\\s+", "", fit$prior$prior)
  if (!any(fit$prior$class == "b" & prior_text == "normal(0,0.5)") ||
      !any(fit$prior$class == "Intercept" & prior_text == "normal(0,3)"))
    stop("Ordinal cache priors differ: ", target)
  invisible(TRUE)
}
fit_ordinal_targets <- function(ctx, selector = "all", action = "fit") {
  selected <- switch(selector, all = ctx$nodes, support = c("kin", "offline", "online"),
                     anxiety = grep("^gad", ctx$nodes, value = TRUE), mpiq = grep("^mpiq", ctx$nodes, value = TRUE),
                     strsplit(selector, ",", fixed = TRUE)[[1]])
  if (!all(selected %in% ctx$nodes)) stop("Unknown ordinal targets: ", selector)
  for (target in selected) {
    path <- file.path(ctx$directory, "models", paste0(target, ".rds"))
    if (file.exists(path)) {
      bundle <- readRDS(path)
      verify_ordinal_bundle(bundle, target, ctx)
      sm <- sampling_summary(bundle$fit, ctx$cfg$ordinal_max_treedepth)
      if (!isTRUE(sm$summary$pass)) stop("Cached fit requires diagnostic review: ", target)
      message("Verified cache; no refit: ", target)
      next
    }
    if (action == "plan") { message("Not fitted: ", target); next }
    lock <- paste0(path, ".lock")
    if (!dir.create(lock, showWarnings = FALSE)) stop("Target is locked by another process: ", target)
    tryCatch({
      message("Fitting ordinal target: ", target)
      start <- Sys.time()
      fit <- brms::brm(ordinal_formula(target, ctx), data = ctx$data,
                        family = brms::cumulative(link = "probit"),
                        prior = c(brms::set_prior("normal(0,0.5)", class = "b"),
                                  brms::set_prior("normal(0,3)", class = "Intercept")),
                        chains = ctx$cfg$chains, cores = min(ctx$cfg$cores, ctx$cfg$chains),
                        iter = ctx$cfg$iter, warmup = ctx$cfg$warmup,
                        seed = ctx$cfg$ordinal_seed + match(target, ctx$nodes), init = 0,
                        backend = "cmdstanr", control = list(adapt_delta = ctx$cfg$ordinal_adapt_delta,
                                                              max_treedepth = ctx$cfg$ordinal_max_treedepth))
      sm <- sampling_summary(fit, ctx$cfg$ordinal_max_treedepth)
      sm$summary$target <- target
      sm$summary$elapsed_minutes <- as.numeric(difftime(Sys.time(), start, units = "mins"))
      write_csv(sm$summary, file.path(ctx$directory, paste0("diagnostics_", target, ".csv")))
      write_csv(sm$parameters, file.path(ctx$directory, paste0("parameters_", target, ".csv")))
      write_csv(sm$chains, file.path(ctx$directory, paste0("chains_", target, ".csv")))
      bundle <- list(fit = fit, spec = ordinal_spec(target, ctx))
      temp <- tempfile(tmpdir = dirname(path), fileext = ".rds")
      saveRDS(bundle, temp)
      if (!file.rename(temp, path)) stop("Could not finalize ordinal cache: ", target)
      if (!isTRUE(sm$summary$pass)) stop("Fit saved; sampling diagnostics require review: ", target)
    }, finally = unlink(lock, recursive = FALSE))
  }
  writeLines(capture.output(sessionInfo()), file.path(ctx$directory, paste0("sessionInfo_", gsub(",", "_", selector), ".txt")))
}
ordinal_group_slope <- function(draws, source, group) {
  z <- paste0("z_", source, "_t1")
  terms <- paste0("b_", z)
  if (group %in% c("addicted_female", "addicted_male")) terms <- c(terms, paste0("b_addictionAddiction:", z))
  if (group %in% c("nonaddicted_male", "addicted_male")) terms <- c(terms, paste0("b_genderMale:", z))
  if (group == "addicted_male") terms <- c(terms, paste0("b_addictionAddiction:genderMale:", z))
  if (!all(terms %in% colnames(draws))) stop("Missing ordinal group-slope terms.")
  rowSums(draws[, terms, drop = FALSE])
}
ordinal_jacobian <- function(fit, target, ctx, chunk = 250L, draw_ids = NULL) {
  draws <- as.matrix(fit)
  original_ids <- if (is.null(draw_ids)) seq_len(nrow(draws)) else as.integer(draw_ids)
  if (!length(original_ids) || anyNA(original_ids) || any(original_ids < 1L | original_ids > nrow(draws)))
    stop("Invalid posterior draw indices.")
  draws <- draws[original_ids, , drop = FALSE]
  values <- as.numeric(levels(fit$data[[paste0(target, "_t2")]]))
  design <- model.matrix(fit$formula$formula, fit$data)
  design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  slopes <- paste0("b_", colnames(design))
  cuts <- paste0("b_Intercept[", seq_len(length(values) - 1L), "]")
  if (!all(c(slopes, cuts) %in% colnames(draws))) stop("Missing ordinal parameters: ", target)
  result <- array(NA_real_, c(nrow(draws), length(ctx$nodes), length(ctx$groups)),
                  dimnames = list(NULL, ctx$nodes, names(ctx$groups)))
  for (begin in seq.int(1L, nrow(draws), by = chunk)) {
    ids <- begin:min(begin + chunk - 1L, nrow(draws))
    eta <- draws[ids, slopes, drop = FALSE] %*% t(design)
    phi <- matrix(0, length(ids), nrow(ctx$data))
    for (k in seq_along(cuts)) phi <- phi + diff(values)[k] *
      dnorm(matrix(draws[ids, cuts[k]], length(ids), nrow(ctx$data)) - eta)
    for (group in names(ctx$groups)) {
      average <- rowMeans(phi[, ctx$groups[[group]], drop = FALSE])
      for (source in ctx$nodes) result[ids, source, group] <-
        ordinal_group_slope(draws[ids, , drop = FALSE], source, group) / ctx$sds[[source]] *
        average * ctx$scales[[source]] / ctx$scales[[target]]
    }
  }
  # Verify analytic expected-score derivatives against brms finite differences.
  ids <- unique(round(seq(1, nrow(draws), length.out = min(3L, nrow(draws)))))
  rows <- vapply(ctx$groups, function(x) x[1], integer(1))
  base <- fit$data[rows, , drop = FALSE]
  expected_score <- function(prob) {
    out <- matrix(0, dim(prob)[1], dim(prob)[2])
    for (k in seq_along(values)) out <- out + values[k] * prob[, , k]
    out
  }
  eta <- brms::posterior_linpred(fit, newdata = base, draw_ids = original_ids[ids], transform = FALSE)
  analytic_score <- matrix(max(values), length(ids), nrow(base))
  phi <- matrix(0, length(ids), nrow(base))
  colnames(phi) <- names(ctx$groups)
  for (k in seq_along(cuts)) {
    u <- matrix(draws[ids, cuts[k]], length(ids), nrow(base)) - eta
    analytic_score <- analytic_score - diff(values)[k] * pnorm(u)
    phi <- phi + diff(values)[k] * dnorm(u)
  }
  score_error <- max(abs(analytic_score - expected_score(brms::posterior_epred(fit, newdata = base, draw_ids = original_ids[ids]))))
  perturb_data <- function(direction) {
    do.call(rbind, lapply(ctx$nodes, function(source) {
      shifted <- base
      column <- paste0("z_", source, "_t1")
      shifted[[column]] <- shifted[[column]] + direction * 0.001 * ctx$scales[[source]] / ctx$sds[[source]]
      shifted
    }))
  }
  plus_score <- expected_score(brms::posterior_epred(fit, newdata = perturb_data(1), draw_ids = original_ids[ids]))
  minus_score <- expected_score(brms::posterior_epred(fit, newdata = perturb_data(-1), draw_ids = original_ids[ids]))
  derivative_error <- 0
  for (source_index in seq_along(ctx$nodes)) {
    source <- ctx$nodes[source_index]
    h <- 0.001 * ctx$scales[[source]]
    columns <- (source_index - 1L) * nrow(base) + seq_len(nrow(base))
    finite <- (plus_score[, columns, drop = FALSE] - minus_score[, columns, drop = FALSE]) / (2 * h)
    analytic <- sapply(names(ctx$groups), function(g) ordinal_group_slope(draws[ids, , drop = FALSE], source, g) / ctx$sds[[source]] * phi[, g])
    derivative_error <- max(derivative_error, abs(finite - analytic))
  }
  if (score_error > 1e-8 || derivative_error > 1e-5) stop("Ordinal derivative verification failed: ", target)
  list(jacobian = result, qa = data.frame(target = target, expected_score_max_error = score_error,
                                         raw_derivative_max_error = derivative_error))
}
ordinal_prediction_engine <- function(target, cfg) {
  engine <- if (is.null(cfg$ordinal_ppc_engine)) "inverse_cdf" else cfg$ordinal_ppc_engine
  overrides <- cfg$ordinal_ppc_engine_overrides
  if (!is.null(overrides) && target %in% names(overrides)) engine <- unname(overrides[[target]])
  if (length(engine) != 1L || !engine %in% c("inverse_cdf", "brms"))
    stop("Ordinal prediction engine must be 'inverse_cdf' or 'brms'.")
  engine
}
ordinal_category_checks <- function(fit, target, ctx) {
  set.seed(ctx$cfg$ordinal_ppc_seed + match(target, ctx$nodes))
  available <- posterior::ndraws(posterior::as_draws_array(fit))
  ids <- sample.int(available, min(available, ctx$cfg$ordinal_ppc_draws))
  values <- as.numeric(levels(fit$data[[paste0(target, "_t2")]]))
  if (length(values) < 2L || any(!is.finite(values)) || any(diff(values) <= 0))
    stop("Ordinal response categories must be increasing finite scores.")
  engine <- ordinal_prediction_engine(target, ctx$cfg)
  if (engine == "brms") {
    indices <- brms::posterior_predict(fit, draw_ids = ids)
  } else {
    draws <- as.matrix(fit)
    design <- model.matrix(fit$formula$formula, fit$data)
    design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
    slopes <- paste0("b_", colnames(design))
    cuts <- paste0("b_Intercept[", seq_len(length(values) - 1L), "]")
    if (!all(c(slopes, cuts) %in% colnames(draws))) stop("Missing ordinal prediction parameters.")
    eta <- draws[ids, slopes, drop = FALSE] %*% t(design)
    u <- matrix(runif(length(ids) * nrow(fit$data)), length(ids), nrow(fit$data))
    indices <- matrix(1L, length(ids), nrow(fit$data))
    for (k in seq_along(cuts)) indices <- indices +
      (u > pnorm(matrix(draws[ids, cuts[k]], length(ids), nrow(fit$data)) - eta))
  }
  if (any(!indices %in% seq_along(values))) stop("Unexpected ordinal prediction category indices.")
  replicated <- matrix(values[indices], nrow(indices), ncol(indices))
  observed <- ctx$raw[[paste0(target, "_t2")]]
  do.call(rbind, lapply(names(ctx$groups), function(group) {
    rows <- ctx$groups[[group]]
    do.call(rbind, lapply(values, function(category) {
      obs <- mean(observed[rows] == category)
      reps <- rowMeans(replicated[, rows, drop = FALSE] == category)
      interval <- unname(quantile(reps, c(.025, .975)))
      data.frame(target = target, group = group, category = category, n = length(rows),
                  observed_fraction = obs, replicated_mean = mean(reps),
                  predictive_low = interval[1], predictive_high = interval[2],
                  outside_pointwise_interval = obs < interval[1] || obs > interval[2],
                  absolute_fraction_gap = abs(obs - mean(reps)), replicates = length(ids),
                  prediction_engine = if (engine == "brms") "brms posterior_predict category indices"
                    else "exact cumulative-probit inverse CDF")
    }))
  }))
}
summarize_ordinal <- function(ctx) {
  arrays <- NULL
  diagnostics <- ppc <- qa <- list()
  for (target in ctx$nodes) {
    path <- file.path(ctx$directory, "models", paste0(target, ".rds"))
    if (!file.exists(path)) stop("Missing ordinal fit: ", target, ". Fit missing targets before summarizing.")
    bundle <- readRDS(path)
    verify_ordinal_bundle(bundle, target, ctx)
    sm <- sampling_summary(bundle$fit, ctx$cfg$ordinal_max_treedepth)
    if (!isTRUE(sm$summary$pass)) stop("Sampling diagnostic review required: ", target)
    sm$summary$target <- target
    diagnostics[[target]] <- sm$summary
    checkpoint <- file.path(ctx$directory, paste0("jacobian_", target, ".rds"))
    spec <- list(fit_md5 = unname(tools::md5sum(path)), scales = ctx$scales,
                 ppc_draws = ctx$cfg$ordinal_ppc_draws, ppc_seed = ctx$cfg$ordinal_ppc_seed,
                 prediction_engine = ordinal_prediction_engine(target, ctx$cfg),
                 method = "expected-score derivative and categorical prediction v2")
    if (file.exists(checkpoint)) {
      one <- readRDS(checkpoint)
      if (!identical(one$spec, spec)) stop("Ordinal checkpoint differs from current specification: ", target)
    } else {
      one <- ordinal_jacobian(bundle$fit, target, ctx)
      one$ppc <- ordinal_category_checks(bundle$fit, target, ctx)
      one$spec <- spec
      saveRDS(one, checkpoint)
    }
    if (is.null(arrays)) arrays <- setNames(lapply(ctx$groups, function(x)
      array(NA_real_, c(dim(one$jacobian)[1], length(ctx$nodes), length(ctx$nodes)),
            dimnames = list(NULL, source = ctx$nodes, target = ctx$nodes))), names(ctx$groups))
    for (g in names(arrays)) arrays[[g]][, , target] <- one$jacobian[, , g]
    ppc[[target]] <- one$ppc
    qa[[target]] <- one$qa
    rm(bundle)
    gc(verbose = FALSE)
  }
  metrics <- calculate_metric_draws(arrays)
  contrasts <- summarize_metric_contrasts(metrics)
  group_summary <- dplyr::bind_rows(lapply(dimnames(metrics)$group, function(g)
    dplyr::bind_rows(lapply(dimnames(metrics)$metric, function(m)
      data.frame(group = g, metric = m, as.list(posterior_summary(metrics[,m,g])))))))
  write_csv(group_summary, file.path(ctx$directory, "group_network_metric_summaries.csv"))
  means <- posterior_mean_networks(arrays)
  for (g in names(means)) utils::write.csv(means[[g]], file.path(ctx$directory, paste0("mean_network_", g, ".csv")))
  write_csv(contrasts, file.path(ctx$directory, "network_metric_contrasts.csv"))
  saveRDS(metrics, file.path(ctx$directory, "network_metric_draws.rds"))
  saveRDS(arrays, file.path(ctx$directory, "expected_score_network_draws.rds"))
  write_csv(do.call(rbind, diagnostics), file.path(ctx$directory, "all_node_mcmc_diagnostics.csv"))
  write_csv(do.call(rbind, qa), file.path(ctx$directory, "expected_score_derivative_checks.csv"))
  write_csv(do.call(rbind, ppc), file.path(ctx$directory, "category_predictive_checks.csv"))
  write_csv(data.frame(node = ctx$nodes, baseline_mean = ctx$means, baseline_sd = ctx$sds,
                       crosswave_sd = ctx$scales), file.path(ctx$directory, "node_scales.csv"))
  write_csv(data.frame(group = names(ctx$groups), n = vapply(ctx$groups, length, integer(1))), file.path(ctx$directory, "group_sizes.csv"))
  write_distance_table(summarize_network_distances(network_distance_draws(arrays)), ctx$directory)
  # Reorder entire target columns, retaining dependence among sources and groups.
  pairing <- list(data.frame(pairing_seed = NA_integer_, contrasts))
  for (seed in ctx$cfg$pairing_seeds) {
    set.seed(seed)
    alternative <- arrays
    for (target in ctx$nodes) {
      order <- sample.int(dim(arrays[[1]])[1])
      for (g in names(arrays)) alternative[[g]][, , target] <- arrays[[g]][order, , target]
    }
    pairing[[length(pairing) + 1L]] <- data.frame(pairing_seed = seed,
      summarize_metric_contrasts(calculate_metric_draws(alternative)))
  }
  write_csv(do.call(rbind, pairing), file.path(ctx$directory, "draw_pairing_sensitivity.csv"))
  gaussian_path <- file.path(output_dir, "network_metric_contrasts.csv")
  if (file.exists(gaussian_path)) {
    gaussian <- read.csv(gaussian_path, check.names = FALSE)
    comparison <- merge(gaussian, contrasts, by = c("category", "metric", "contrast"), suffixes = c("_gaussian", "_ordinal"), all = TRUE)
    comparison$same_mean_direction <- sign(comparison$posterior_mean_gaussian) == sign(comparison$posterior_mean_ordinal)
    for (threshold in c(.90, .95)) comparison[[paste0("retained_pd_", threshold)]] <-
      comparison$probability_of_direction_ordinal > threshold & comparison$same_mean_direction
    write_csv(comparison, file.path(ctx$directory, "gaussian_ordinal_comparison.csv"))
  }
  gaussian_groups <- file.path(output_dir, "group_network_metric_summaries.csv")
  if (file.exists(gaussian_groups)) {
    write_group_support_table(read.csv(gaussian_groups, check.names = FALSE), group_summary,
      data.frame(group = names(ctx$groups), n = vapply(ctx$groups, length, integer(1))), ctx$directory)
  } else {
    message("Supplementary Table 15 requires Gaussian group summaries; run gaussian_results before ordinal_results.")
  }
  table <- contrasts |>
    dplyr::mutate(Category = vapply(.data$metric, metric_section_label, character(1)),
      Metric = metric_display_label(.data$metric),
      Comparison = dplyr::recode(.data$contrast,
        `Addicted minus non-addicted, averaged across sexes` = "Equal-sex average",
        `Addicted minus non-addicted, female` = "Women",
        `Addicted minus non-addicted, male` = "Men",
        `Male minus female addiction contrast` = "Men minus women"),
      Estimate = paste0(format_interval(.data$posterior_mean, .data$cri_lower, .data$cri_upper),
                         " (Pd = ", format_probability(.data$probability_of_direction, .data$n_draws), ")")) |>
    dplyr::select("Category", "Metric", "Comparison", "Estimate") |>
    tidyr::pivot_wider(names_from = "Comparison", values_from = "Estimate") |>
    sectioned_table()
  save_publication_table(table, "Supplementary Table 14 | Bayesian ordinal outcome sensitivity: posterior network-metric contrasts",
    paste0("Means, 95% credible intervals and posterior direction probabilities refer to group-averaged observed-score derivative matrices. Outcomes were fitted separately; residual correlations between outcomes were not estimated. The ", impulse_steps, "-step propagation holds the derivative matrix fixed."),
    file.path(ctx$directory, "Supplementary_Table_14.docx"), landscape = TRUE)
  diagnostic_table <- do.call(rbind, diagnostics) |>
    dplyr::transmute(Target = .data$target, `Maximum R-hat` = sprintf("%.4f", .data$max_rhat),
      `Minimum bulk ESS` = round(.data$min_bulk), `Minimum tail ESS` = round(.data$min_tail),
      `Minimum E-BFMI` = sprintf("%.3f", .data$min_E_BFMI),
      `Divergences / depth hits` = paste0(.data$divergences, " / ", .data$depth_hits))
  save_publication_table(diagnostic_table, "Supplementary Table 13 | MCMC diagnostics for Bayesian ordinal models",
    "Divergences and maximum-tree-depth hits are separate counts. R-hat and ESS assess computation; they do not measure predictive fit.",
    file.path(ctx$directory, "Supplementary_Table_13.docx"), landscape = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(ctx$directory, "sessionInfo_postprocessing.txt"))
}
