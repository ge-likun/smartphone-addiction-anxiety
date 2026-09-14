spectral_radius <- function(matrix) {
  max(Mod(eigen(matrix, only.values = TRUE)$values))
}

impulse_response <- function(matrix, source, target_nodes, steps = impulse_steps) {
  state <- setNames(rep(0, nrow(matrix)), rownames(matrix))
  state[source] <- 1
  cumulative <- 0
  for (step in seq_len(steps)) {
    state <- drop(t(matrix) %*% state[colnames(matrix)])
    names(state) <- colnames(matrix)
    cumulative <- cumulative + sum(state[target_nodes])
  }
  cumulative
}

network_metrics <- function(matrix, communities = network_communities) {
  support <- communities$Support
  anxiety <- communities$Anxiety
  addiction <- communities$Smartphone_addiction
  nodes <- rownames(matrix)

  stability <- c(
    `Spectral radius: support` = spectral_radius(matrix[support, support, drop = FALSE]),
    `Spectral radius: anxiety` = spectral_radius(matrix[anxiety, anxiety, drop = FALSE]),
    `Spectral radius: smartphone addiction` = spectral_radius(matrix[addiction, addiction, drop = FALSE])
  )

  bridge <- vapply(nodes, function(node) {
    own <- names(communities)[vapply(communities, function(x) node %in% x, logical(1))]
    other <- unlist(communities[setdiff(names(communities), own)], use.names = FALSE)
    sum(abs(matrix[node, other])) + sum(abs(matrix[other, node]))
  }, numeric(1))
  names(bridge) <- paste0("Bridge strength: ", nodes)

  outgoing <- unlist(lapply(support, function(node) {
    c(
      setNames(sum(matrix[node, anxiety]), paste0("Outgoing sum: ", node, " to anxiety")),
      setNames(sum(matrix[node, addiction]), paste0("Outgoing sum: ", node, " to smartphone addiction"))
    )
  }))

  propagated <- unlist(lapply(support, function(node) {
    c(
      setNames(
        impulse_response(matrix, node, anxiety),
        paste0("Cumulative response: ", node, " to anxiety")
      ),
      setNames(
        impulse_response(matrix, node, addiction),
        paste0("Cumulative response: ", node, " to smartphone addiction")
      )
    )
  }))

  c(stability, bridge, outgoing, propagated)
}

calculate_metric_draws <- function(arrays) {
  metric_names <- names(network_metrics(arrays[[1]][1, , ]))
  result <- array(
    NA_real_,
    dim = c(dim(arrays[[1]])[1], length(metric_names), length(arrays)),
    dimnames = list(NULL, metric = metric_names, group = names(arrays))
  )

  for (group in names(arrays)) {
    for (draw in seq_len(dim(arrays[[group]])[1])) {
      result[draw, , group] <- network_metrics(arrays[[group]][draw, , ])
    }
  }
  result
}

posterior_summary <- function(x) {
  direction_probability <- max(mean(x > 0), mean(x < 0))
  c(
    posterior_mean = mean(x),
    cri_lower = unname(stats::quantile(x, 0.025)),
    cri_upper = unname(stats::quantile(x, 0.975)),
    probability_of_direction = direction_probability
  )
}

metric_category <- function(metric) {
  if (grepl("^Spectral radius", metric)) return("Spectral radius")
  if (grepl("^Bridge strength", metric)) return("Bridge strength")
  if (grepl("^Outgoing sum", metric)) return("Outgoing support sum")
  "Impulse propagation"
}

metric_display_label <- function(metric) {
  label <- sub("^[^:]+: ", "", metric)
  label <- sub("^support$", "Social support", label)
  label <- sub("^anxiety$", "Anxiety", label)
  label <- sub("^smartphone addiction$", "Smartphone addiction", label)
  label <- sub("^gad([1-7])$", "GAD-7 item \\1", label)
  label <- sub("^mpiq([1-8])$", "MPIQ item \\1", label)
  label <- sub("^kin$", "Kinship support", label)
  label <- sub("^offline$", "Offline-friend support", label)
  label <- sub("^online$", "Online-friend support", label)
  label <- sub("^kin to ", "Kinship support to ", label)
  label <- sub("^offline to ", "Offline-friend support to ", label)
  label <- sub("^online to ", "Online-friend support to ", label)
  label
}

summarize_metric_contrasts <- function(metric_draws) {
  output <- lapply(dimnames(metric_draws)$metric, function(metric) {
    female <- metric_draws[, metric, "addicted_female"] -
      metric_draws[, metric, "nonaddicted_female"]
    male <- metric_draws[, metric, "addicted_male"] -
      metric_draws[, metric, "nonaddicted_male"]
    overall <- (female + male) / 2
    sex_difference <- male - female

    contrasts <- list(
      `Addicted minus non-addicted, averaged across sexes` = overall,
      `Addicted minus non-addicted, female` = female,
      `Addicted minus non-addicted, male` = male,
      `Male minus female addiction contrast` = sex_difference
    )

    do.call(rbind, lapply(names(contrasts), function(contrast) {
      data.frame(
        category = metric_category(metric),
        metric = metric,
        contrast = contrast,
        as.list(posterior_summary(contrasts[[contrast]])),
        check.names = FALSE
      )
    }))
  })

  dplyr::bind_rows(output)
}

network_distance_draws <- function(arrays) {
  pairs <- utils::combn(names(arrays), 2, simplify = FALSE)
  output <- lapply(pairs, function(pair) {
    difference <- arrays[[pair[1]]] - arrays[[pair[2]]]
    frobenius <- sqrt(apply(difference^2, 1, sum))
    maed <- apply(abs(difference), 1, mean)
    data.frame(
      draw = seq_along(frobenius),
      group_1 = pair[1],
      group_2 = pair[2],
      comparison = paste(pair, collapse = " versus "),
      frobenius = frobenius,
      maed = maed
    )
  })
  dplyr::bind_rows(output)
}

summarize_network_distances <- function(distance_draws) {
  pairwise <- distance_draws |>
    dplyr::group_by(.data$comparison) |>
    dplyr::summarise(
      frobenius_mean = mean(.data$frobenius),
      frobenius_lower = stats::quantile(.data$frobenius, 0.025),
      frobenius_upper = stats::quantile(.data$frobenius, 0.975),
      maed_mean = mean(.data$maed),
      maed_lower = stats::quantile(.data$maed, 0.025),
      maed_upper = stats::quantile(.data$maed, 0.975),
      .groups = "drop"
    )

  wide <- distance_draws |>
    dplyr::select("draw", "comparison", "frobenius", "maed") |>
    tidyr::pivot_wider(names_from = "comparison", values_from = c("frobenius", "maed"))

  male_pair <- "nonaddicted_male versus addicted_male"
  female_pair <- "nonaddicted_female versus addicted_female"
  addicted_pair <- "addicted_female versus addicted_male"
  nonaddicted_pair <- "nonaddicted_female versus nonaddicted_male"

  differences <- list(
    `Male minus female addiction-related distance` = list(male_pair, female_pair),
    `Addicted minus non-addicted sex-related distance` = list(addicted_pair, nonaddicted_pair)
  )

  distance_differences <- dplyr::bind_rows(lapply(names(differences), function(label) {
    pair <- differences[[label]]
    frobenius <- wide[[paste0("frobenius_", pair[[1]])]] - wide[[paste0("frobenius_", pair[[2]])]]
    maed <- wide[[paste0("maed_", pair[[1]])]] - wide[[paste0("maed_", pair[[2]])]]
    data.frame(
      comparison = label,
      frobenius_mean = mean(frobenius),
      frobenius_lower = stats::quantile(frobenius, 0.025),
      frobenius_upper = stats::quantile(frobenius, 0.975),
      frobenius_pr_positive = mean(frobenius > 0),
      maed_mean = mean(maed),
      maed_lower = stats::quantile(maed, 0.025),
      maed_upper = stats::quantile(maed, 0.975),
      maed_pr_positive = mean(maed > 0),
      check.names = FALSE
    )
  }))

  largest <- distance_draws |>
    dplyr::group_by(.data$draw) |>
    dplyr::mutate(
      largest_frobenius = .data$frobenius == max(.data$frobenius),
      largest_maed = .data$maed == max(.data$maed)
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$comparison) |>
    dplyr::summarise(
      probability_largest_frobenius = mean(.data$largest_frobenius),
      probability_largest_maed = mean(.data$largest_maed),
      .groups = "drop"
    )

  list(pairwise = pairwise, differences = distance_differences, largest = largest)
}
