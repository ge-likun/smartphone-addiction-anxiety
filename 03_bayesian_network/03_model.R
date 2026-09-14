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
    cores = min(chains, available_cores),
    threads = brms::threading(threads_per_chain),
    backend = "cmdstanr",
    seed = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    save_pars = brms::save_pars(all = TRUE)
  )

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
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

validate_network_model <- function(fit, nodes = network_nodes) {
  expected_data <- network_required_columns(nodes)
  absent_data <- setdiff(expected_data, names(fit$data))
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

  invisible(TRUE)
}
