extract_edge_draws <- function(fit, nodes = network_nodes) {
  parameters <- expected_edge_parameters(nodes)
  draws <- posterior::as_draws_matrix(fit, variable = parameters)
  missing_parameters <- setdiff(parameters, colnames(draws))
  if (length(missing_parameters)) {
    stop("Missing posterior parameters: ", paste(missing_parameters, collapse = ", "), call. = FALSE)
  }
  draws[, parameters, drop = FALSE]
}

reconstruct_network_arrays <- function(edge_draws, nodes = network_nodes) {
  n_draws <- nrow(edge_draws)
  n_nodes <- length(nodes)
  dimensions <- c(n_draws, n_nodes, n_nodes)
  dim_names <- list(NULL, source = nodes, target = nodes)

  arrays <- list(
    nonaddicted_female = array(0, dimensions, dimnames = dim_names),
    addicted_female = array(0, dimensions, dimnames = dim_names),
    nonaddicted_male = array(0, dimensions, dimnames = dim_names),
    addicted_male = array(0, dimensions, dimnames = dim_names)
  )

  for (target in nodes) {
    response <- gsub("_", "", paste0(target, "_t2"))
    for (source in nodes) {
      predictor <- paste0(source, "_t1")
      names <- paste0(
        "b_", response,
        c(
          paste0("_", predictor),
          paste0("_addictionAddiction:", predictor),
          paste0("_genderMale:", predictor),
          paste0("_addictionAddiction:genderMale:", predictor)
        )
      )
      coefficients <- edge_draws[, names, drop = FALSE]

      arrays$nonaddicted_female[, source, target] <- coefficients[, 1]
      arrays$addicted_female[, source, target] <- coefficients[, 1] + coefficients[, 2]
      arrays$nonaddicted_male[, source, target] <- coefficients[, 1] + coefficients[, 3]
      arrays$addicted_male[, source, target] <- rowSums(coefficients)
    }
  }

  arrays
}

common_node_scales <- function(model_data, nodes = network_nodes) {
  scales <- vapply(nodes, function(node) {
    stats::sd(c(model_data[[paste0(node, "_t1")]], model_data[[paste0(node, "_t2")]]))
  }, numeric(1))
  if (any(!is.finite(scales) | scales <= 0)) {
    stop("Every node must have a positive finite standard deviation.", call. = FALSE)
  }
  scales
}

standardize_network_array <- function(array, scales) {
  multiplier <- outer(scales, 1 / scales)
  sweep(array, c(2, 3), multiplier, `*`)
}

standardize_network_arrays <- function(arrays, model_data, nodes = network_nodes) {
  scales <- common_node_scales(model_data, nodes)
  lapply(arrays, standardize_network_array, scales = scales)
}

posterior_mean_networks <- function(arrays) {
  lapply(arrays, function(x) apply(x, c(2, 3), mean))
}

