format_interval <- function(mean, lower, upper, digits = 3) {
  sprintf(
    paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
    mean, lower, upper
  )
}

format_probability <- function(x, n = NULL, digits = 3L) {
  vapply(x, function(p) {
    if (!is.finite(p)) return("")
    if (p == 0) return(if (length(n)) paste0("0/", n[1]) else "0 (MC)")
    if (p == 1) return(if (length(n)) paste0(n[1], "/", n[1]) else "1 (MC)")
    boundary <- 10^-digits
    if (p < boundary) return(paste0("<", formatC(boundary, format = "f", digits = digits)))
    if (p > 1 - boundary) return(paste0(">", formatC(1 - boundary, format = "f", digits = digits)))
    formatC(p, format = "f", digits = digits)
  }, character(1))
}

sectioned_table <- function(data, section_column = "Category") {
  sections <- unique(data[[section_column]])
  columns <- setdiff(names(data), section_column)
  parts <- lapply(sections, function(section) {
    heading <- as.data.frame(setNames(rep(list(""), length(columns)), columns),
                             check.names = FALSE)
    heading[[1]] <- section
    rbind(heading, data[data[[section_column]] == section, columns, drop = FALSE])
  })
  result <- do.call(rbind, parts)
  attr(result, "section_rows") <- cumsum(c(1L, head(vapply(parts, nrow, integer(1)), -1L)))
  result
}

metric_section_label <- function(metric) {
  switch(metric_category(metric),
    `Spectral radius` = "Sub-network Stability (spectral radius)",
    `Bridge strength` = "Bridge Centrality (bridge strength)",
    `Outgoing support sum` = "Potential Influence (outgoing support sum)",
    `Impulse propagation` = paste0("Cumulative Propagated Response (", impulse_steps,
                                  "-step impulse propagation)"))
}

three_line_flextable <- function(data, widths = NULL) {
  border <- officer::fp_border(color = "black", width = 1)
  table <- flextable::flextable(data)
  table <- flextable::border_remove(table)
  table <- flextable::hline_top(table, border = border, part = "header")
  table <- flextable::hline_bottom(table, border = border, part = "header")
  table <- flextable::hline_bottom(table, border = border, part = "body")
  table <- flextable::font(table, fontname = "Times New Roman", part = "all")
  table <- flextable::fontsize(table, size = 10, part = "all")
  table <- flextable::bold(table, part = "header")
  table <- flextable::align(table, align = "center", part = "header")
  table <- flextable::align(table, align = "center", part = "body")
  table <- flextable::align(table, j = 1, align = "left", part = "body")
  for (i in attr(data, "section_rows")) {
    table <- flextable::merge_at(table, i = i, j = seq_len(ncol(data)), part = "body")
    table <- flextable::bold(table, i = i, part = "body")
    table <- flextable::align(table, i = i, align = "left", part = "body")
  }
  table <- flextable::valign(table, valign = "center", part = "all")
  table <- flextable::padding(table, padding = 3, part = "all")
  table <- flextable::autofit(table)
  if (!is.null(widths)) table <- flextable::width(table, width = widths)
  flextable::set_table_properties(table, layout = "autofit", width = 1)
}

save_publication_table <- function(data, title, note, file, landscape = FALSE) {
  document <- officer::read_docx()
  title_format <- officer::fp_text(
    font.family = "Times New Roman",
    font.size = 11,
    bold = TRUE,
    color = "000000"
  )
  note_format <- officer::fp_text(
    font.family = "Times New Roman",
    font.size = 10,
    color = "000000"
  )
  document <- officer::body_add_fpar(
    document,
    officer::fpar(officer::ftext(title, title_format))
  )
  document <- flextable::body_add_flextable(document, three_line_flextable(data))
  document <- officer::body_add_fpar(
    document,
    officer::fpar(officer::ftext(paste0("Note. ", note), note_format))
  )
  if (landscape) document <- officer::body_end_section_landscape(document)

  staged_file <- tempfile(pattern = "network-table-", fileext = ".docx")
  on.exit(unlink(staged_file), add = TRUE)
  print(document, target = staged_file)
  if (!file.copy(staged_file, file, overwrite = TRUE)) {
    stop("Could not save table: ", file)
  }
  invisible(file)
}

write_diagnostic_table <- function(diagnostics, directory) {
  table <- diagnostics$table
  names(table) <- c("Diagnostic", "Number assessed", "Result")
  note <- paste(
    "Rhat, bulk effective sample size (ESS) and tail ESS summarize MCMC mixing.",
    "E-BFMI is the energy Bayesian fraction of missing information.",
    "Frequent maximum-tree-depth hits indicate limited sampling efficiency. These diagnostics summarize computation and do not assess predictive fit."
  )
  utils::write.csv(diagnostics$full, file.path(directory, "Supplementary_Table_9_source_data.csv"), row.names = FALSE)
  save_publication_table(
    table,
    "Supplementary Table 9 | MCMC diagnostics for the Bayesian longitudinal network model",
    note,
    file.path(directory, "Supplementary_Table_9.docx")
  )
}

write_distance_table <- function(summaries, directory) {
  pairwise <- summaries$pairwise |>
    dplyr::transmute(
      Category = "Pairwise network distances",
      Comparison = .data$comparison,
      Frobenius = format_interval(.data$frobenius_mean, .data$frobenius_lower, .data$frobenius_upper),
      MAED = format_interval(.data$maed_mean, .data$maed_lower, .data$maed_upper)
    )
  differences <- summaries$differences |>
    dplyr::transmute(
      Category = "Distance differences",
      Comparison = .data$comparison,
      Frobenius = paste0(format_interval(.data$frobenius_mean, .data$frobenius_lower, .data$frobenius_upper),
        " (", format_probability(.data$frobenius_pr_positive, attr(summaries, "n_draws")), ")"),
      MAED = paste0(format_interval(.data$maed_mean, .data$maed_lower, .data$maed_upper),
        " (", format_probability(.data$maed_pr_positive, attr(summaries, "n_draws")), ")")
    )
  largest <- summaries$largest |>
    dplyr::transmute(
      Category = "Probability of largest distance",
      Comparison = .data$comparison,
      Frobenius = format_probability(.data$probability_largest_frobenius, attr(summaries, "n_draws")),
      MAED = format_probability(.data$probability_largest_maed, attr(summaries, "n_draws"))
    )
  table <- dplyr::bind_rows(pairwise, differences, largest)
  for (group in names(network_groups)) {
    table$Comparison <- gsub(group, network_groups[[group]], table$Comparison, fixed = TRUE)
  }
  table <- sectioned_table(table)
  note <- paste(
    "All distances were calculated from standardized directed adjacency matrices.",
    "MAED, mean absolute edge difference; CrI, credible interval.",
    "Pairwise distances are non-negative descriptive measures and are not tests against zero.",
    "Distance-difference cells include the probability of a positive difference in parentheses. Extreme Monte Carlo probabilities are displayed as qualifying draws divided by retained draws; zero qualifying draws do not establish an impossible event."
  )
  utils::write.csv(table, file.path(directory, "Supplementary_Table_10_source_data.csv"), row.names = FALSE)
  save_publication_table(
    table,
    "Supplementary Table 10 | Posterior network-distance summaries across addiction-by-sex groups",
    note,
    file.path(directory, "Supplementary_Table_10.docx"),
    landscape = FALSE
  )
}

prepare_metric_table <- function(summary, contrasts) {
  summary |>
    dplyr::filter(.data$contrast %in% contrasts) |>
    dplyr::transmute(
      Category = vapply(.data$metric, metric_section_label, character(1)),
      Metric = metric_display_label(.data$metric),
      Contrast = .data$contrast,
      `Posterior mean [95% CrI]` = format_interval(.data$posterior_mean, .data$cri_lower, .data$cri_upper),
      `Probability of direction` = format_probability(.data$probability_of_direction, .data$n_draws)
    )
}

write_metric_tables <- function(summary, directory) {
  main <- prepare_metric_table(
    summary,
    "Addicted minus non-addicted, averaged across sexes"
  ) |>
    dplyr::select(-"Contrast") |>
    sectioned_table()
  sex_long <- prepare_metric_table(
    summary,
    c(
      "Addicted minus non-addicted, female",
      "Addicted minus non-addicted, male",
      "Male minus female addiction contrast"
    )
  )
  sex <- sex_long |>
    dplyr::mutate(Contrast = dplyr::recode(.data$Contrast,
      `Addicted minus non-addicted, female` = "Women",
      `Addicted minus non-addicted, male` = "Men",
      `Male minus female addiction contrast` = "Men minus women"),
      Estimate = paste0(.data$`Posterior mean [95% CrI]`, " (Pd = ", .data$`Probability of direction`, ")")) |>
    dplyr::select("Category", "Metric", "Contrast", "Estimate") |>
    tidyr::pivot_wider(names_from = "Contrast", values_from = "Estimate") |>
    sectioned_table()

  common_note <- paste(
    "CrI, credible interval.",
    "Networks use standardized lagged coefficients; rows are source nodes and columns are target nodes.",
    "Bridge strength is the sum of absolute incoming and outgoing cross-community coefficients.",
    "Outgoing sums retain coefficient signs.",
    paste0("Cumulative Propagated Response applies a one-standard-deviation input for ", impulse_steps,
           " propagation steps. Signed responses are observational model summaries. Extreme direction probabilities are Monte Carlo estimates; zero qualifying draws do not imply impossibility.")
  )

  utils::write.csv(main, file.path(directory, "Supplementary_Table_11_source_data.csv"), row.names = FALSE)
  utils::write.csv(sex, file.path(directory, "Supplementary_Table_12_source_data.csv"), row.names = FALSE)
  save_publication_table(
    main,
    "Supplementary Table 11 | Addiction-status contrasts in posterior network metrics",
    paste(common_note, "The reported contrast is addicted minus non-addicted, averaged across the female and male contrasts."),
    file.path(directory, "Supplementary_Table_11.docx"),
    landscape = FALSE
  )
  save_publication_table(
    sex,
    "Supplementary Table 12 | Sex-specific addiction-status contrasts and their formal difference",
    paste(common_note, "The male-minus-female column is the formal comparison of the two addiction-status contrasts."),
    file.path(directory, "Supplementary_Table_12.docx"),
    landscape = TRUE
  )
}

write_group_support_table <- function(gaussian, ordinal, group_sizes, directory) {
  required <- c("group", "metric", "n_draws", "posterior_mean", "cri_lower", "cri_upper", "probability_of_direction")
  if (!all(required %in% names(gaussian)) || !all(required %in% names(ordinal)))
    stop("Group-specific network summaries are incomplete.")
  groups <- c(addicted_female = "Women", addicted_male = "Men")
  long <- dplyr::bind_rows(
    dplyr::mutate(gaussian, likelihood = "Gaussian"),
    dplyr::mutate(ordinal, likelihood = "Ordinal")) |>
    dplyr::filter(.data$group %in% names(groups), grepl("^(Outgoing sum|Cumulative response):", .data$metric)) |>
    dplyr::left_join(group_sizes, by = "group")
  if (anyNA(long$n) || anyDuplicated(long[c("group", "metric", "likelihood")]))
    stop("Missing sample sizes or duplicate group-specific network summaries.")
  if (!setequal(unique(long$likelihood), c("Gaussian", "Ordinal")))
    stop("Both likelihood specifications are required for Supplementary Table 15.")
  write.csv(long, file.path(directory, "Supplementary_Table_15_source_data.csv"), row.names = FALSE)
  table <- long |>
    dplyr::mutate(
      Category = vapply(.data$metric, metric_section_label, character(1)),
      Metric = metric_display_label(.data$metric),
      Comparison = paste0(groups[.data$group], ": ", .data$likelihood,
                          " (n = ", format(.data$n, big.mark = ",", trim = TRUE), ")"),
      Estimate = paste0(format_interval(.data$posterior_mean, .data$cri_lower, .data$cri_upper),
                        " (Pd = ", vapply(seq_len(nrow(long)), function(i)
                          format_probability(long$probability_of_direction[i], long$n_draws[i]), character(1)), ")")) |>
    dplyr::select("Category", "Metric", "Comparison", "Estimate") |>
    tidyr::pivot_wider(names_from = "Comparison", values_from = "Estimate") |>
    sectioned_table()
  save_publication_table(table,
    "Supplementary Table 15 | Group-specific support-related network indices in addicted women and men",
    paste0("Cells report posterior means, 95% credible intervals and probability of direction (Pd). Sample sizes and posterior draw counts are calculated from the current data and fits. Signed outgoing sums and ",
           impulse_steps, "-step cumulative responses are absolute within-group indices, distinct from addiction-status contrasts. Positive values indicate higher modelled symptom scores; these observational summaries do not establish causal intervention benefits."),
    file.path(directory, "Supplementary_Table_15.docx"), landscape = TRUE)
  invisible(long)
}
