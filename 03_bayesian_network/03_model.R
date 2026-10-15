build_network_formula <- function(nodes = network_nodes) {
  outcomes <- paste0(nodes, "_t2")
  predictors <- paste0(nodes, "_t1")
  formula_text <- paste0(
    "mvbind(", paste(outcomes, collapse = ", "), ") ~ ",
    "addiction * gender * (", paste(predictors, collapse = " + "), ")"
  )
  brms::bf(stats::as.formula(formula_text)) + brms::set_rescor(TRUE)
}

fit_network_model <- function(
  data,
  output_file = model_file,
  chains = network_chains,
  iter = network_iterations,
  warmup = network_warmup,
  threads_per_chain = network_threads_per_chain,
  seed = network_seed,
  adapt_delta = network_adapt_delta,
  max_treedepth = network_max_treedepth
) {
  formula <- build_network_formula()
  priors <- brms::set_prior("normal(0, 0.5)", class = "b")
  available_cores <- parallel::detectCores(logical = FALSE)
  if (!is.finite(available_cores)) available_cores <- chains

  fit <- brms::brm(
    formula = formula,
    data = data,
    family = stats::gaussian(),
    prior = priors,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = min(chains, as.integer(Sys.getenv("NETWORK_CORES", available_cores))),
    threads = brms::threading(threads_per_chain),
    backend = "cmdstanr",
    seed = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    save_pars = brms::save_pars(all = TRUE)
  )

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  attr(fit, "analysis_spec") <- gaussian_spec(data)
  saveRDS(fit, output_file, compress = FALSE)
  fit
}

load_network_model <- function(path = model_file) {
  if (!file.exists(path)) {
    stop(
      "Fitted model not found: ", path,
      ". Set NETWORK_MODEL_FILE or run run_fit_model.R.",
      call. = FALSE
    )
  }

  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    fit <- readRDS(path)
  } else {
    objects <- new.env(parent = emptyenv())
    loaded <- load(path, envir = objects)
    brms_objects <- loaded[vapply(loaded, function(x) inherits(objects[[x]], "brmsfit"), logical(1))]
    if (length(brms_objects) != 1L) {
      stop("The RData file must contain one brmsfit object.", call. = FALSE)
    }
    fit <- objects[[brms_objects]]
  }

  if (!inherits(fit, "brmsfit")) stop("The model object is not a brmsfit.", call. = FALSE)
  fit
}

expected_edge_parameters <- function(nodes = network_nodes, prefix = "b_") {
  unlist(lapply(nodes, function(target) {
    response <- gsub("_", "", paste0(target, "_t2"))
    unlist(lapply(nodes, function(source) {
      predictor <- paste0(source, "_t1")
      paste0(
        prefix,
        response,
        c(
          paste0("_", predictor),
          paste0("_addictionAddiction:", predictor),
          paste0("_genderMale:", predictor),
          paste0("_addictionAddiction:genderMale:", predictor)
        )
      )
    }))
  }))
}

gaussian_spec <- function(data) {
  root <- Sys.getenv("PROJECT_ROOT", normalizePath("..", winslash = "/"))
  source(file.path(root, "R", "project.R"), encoding = "UTF-8")
  list(data_hash = object_hash(data), versions = software_versions(),
       formula = lapply(build_network_formula()$forms, function(x) paste(deparse(x$formula), collapse = "")),
       family = "gaussian", rescor = TRUE, prior = "normal(0,0.5)",
       chains = network_chains, iter = network_iterations, warmup = network_warmup,
       seed = network_seed, adapt_delta = network_adapt_delta,
       max_treedepth = network_max_treedepth,
       threads_per_chain = network_threads_per_chain)
}

validate_network_model <- function(fit, nodes = network_nodes, expected_data = NULL) {
  required_columns <- network_required_columns(nodes)
  absent_data <- setdiff(required_columns, names(fit$data))
  if (length(absent_data)) {
    stop("The fitted model data are missing: ", paste(absent_data, collapse = ", "), call. = FALSE)
  }

  parameters <- posterior::variables(posterior::as_draws_array(fit))
  missing_parameters <- setdiff(expected_edge_parameters(nodes), parameters)
  if (length(missing_parameters)) {
    stop(
      "The fitted model does not contain all required edge parameters: ",
      paste(missing_parameters, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(expected_data)) {
    cols <- network_required_columns(nodes)
    if (!identical(validate_network_data(fit$data, nodes), validate_network_data(expected_data, nodes)))
      stop("Gaussian cache data differ from current input.")
    if (!identical(attr(fit, "analysis_spec"), gaussian_spec(expected_data)))
      stop("Gaussian cache specification or software version differs. Use a separate results directory for a new fit.")
    expected_formulas <- build_network_formula()$forms
    observed_terms <- lapply(fit$formula$forms, function(x) attr(terms(x$formula), "term.labels"))
    expected_terms <- lapply(expected_formulas, function(x) attr(terms(x$formula), "term.labels"))
    if (!identical(observed_terms, expected_terms))
      stop("Gaussian model formula differs from the intended specification.")
    families <- if (!is.null(fit$family$family)) list(fit$family) else fit$family
    if (!all(vapply(families, function(x) identical(x$family, "gaussian"), logical(1))))
      stop("Expected Gaussian response families.")
    if (!isTRUE(fit$formula$rescor)) stop("Expected jointly estimated residual correlations.")
    b_priors <- fit$prior$prior[fit$prior$class == "b" & nzchar(fit$prior$prior)]
    if (!length(b_priors) || any(gsub("\\s+", "", b_priors) != "normal(0,0.5)"))
      stop("Gaussian slope priors differ from the intended specification.")
    if (posterior::ndraws(posterior::as_draws_array(fit)) != network_chains * (network_iterations - network_warmup))
      stop("Gaussian cache has an unexpected posterior draw count.")
  }

  invisible(TRUE)
}
