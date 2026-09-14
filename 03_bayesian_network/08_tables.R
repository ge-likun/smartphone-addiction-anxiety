format_interval <- function(mean, lower, upper, digits = 3) {
  sprintf(
    paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
    mean, lower, upper
  )
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
  table <- flextable::align(table, j = 1:min(2, ncol(data)), align = "left", part = "body")
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

  print(document, target = file)
  invisible(file)
}

write_diagnostic_table <- function(diagnostics, directory) {
  table <- diagnostics$table
  names(table) <- c("Diagnostic", "Number assessed", "Result")
  note <- paste(
    "Rhat, bulk effective sample size (ESS) and tail ESS summarize MCMC mixing.",
    "E-BFMI is the energy Bayesian fraction of missing information.",
    "Transitions reaching the maximum tree depth are reported because frequent hits indicate incomplete Hamiltonian trajectories even when Rhat and ESS are satisfactory."
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
      Analysis = "Pairwise distance",
      Comparison = .data$comparison,
      Frobenius = format_interval(.data$frobenius_mean, .data$frobenius_lower, .data$frobenius_upper),
      MAED = format_interval(.data$maed_mean, .data$maed_lower, .data$maed_upper),
      `Pr(difference > 0)` = ""
    )
  differences <- summaries$differences |>
    dplyr::transmute(
      Analysis = "Distance difference",
      Comparison = .data$comparison,
      Frobenius = format_interval(.data$frobenius_mean, .data$frobenius_lower, .data$frobenius_upper),
      MAED = format_interval(.data$maed_mean, .data$maed_lower, .data$maed_upper),
      `Pr(difference > 0)` = sprintf("%.3f / %.3f", .data$frobenius_pr_positive, .data$maed_pr_positive)
    )
  largest <- summaries$largest |>
    dplyr::transmute(
      Analysis = "Probability of largest distance",
      Comparison = .data$comparison,
      Frobenius = sprintf("%.3f", .data$probability_largest_frobenius),
      MAED = sprintf("%.3f", .data$probability_largest_maed),
      `Pr(difference > 0)` = ""
    )
  table <- dplyr::bind_rows(pairwise, differences, largest)
  note <- paste(
    "All distances were calculated from standardized directed adjacency matrices.",
    "MAED, mean absolute edge difference; CrI, credible interval.",
    "Pairwise distances are non-negative descriptive measures and are not tests against zero.",
    "For distance-difference rows, the two probabilities refer to Frobenius distance and MAED, respectively."
  )
  utils::write.csv(table, file.path(directory, "Supplementary_Table_10_source_data.csv"), row.names = FALSE)
  save_publication_table(
    table,
    "Supplementary Table 10 | Posterior network-distance summaries across addiction-by-sex groups",
    note,
    file.path(directory, "Supplementary_Table_10.docx"),
    landscape = TRUE
  )
}

prepare_metric_table <- function(summary, contrasts) {
  summary |>
    dplyr::filter(.data$contrast %in% contrasts) |>
    dplyr::transmute(
      Category = .data$category,
      Metric = metric_display_label(.data$metric),
      Contrast = .data$contrast,
      `Posterior mean [95% CrI]` = format_interval(.data$posterior_mean, .data$cri_lower, .data$cri_upper),
      `Probability of direction` = sprintf("%.3f", .data$probability_of_direction)
    )
}

write_metric_tables <- function(summary, directory) {
  main <- prepare_metric_table(
    summary,
    "Addicted minus non-addicted, averaged across sexes"
  ) |>
    dplyr::select(-"Contrast")
  sex <- prepare_metric_table(
    summary,
    c(
      "Addicted minus non-addicted, female",
      "Addicted minus non-addicted, male",
      "Male minus female addiction contrast"
    )
  )

  common_note <- paste(
    "CrI, credible interval.",
    "Networks use standardized lagged coefficients; rows are source nodes and columns are target nodes.",
    "Bridge strength is the sum of absolute incoming and outgoing cross-community coefficients.",
    "Outgoing sums retain coefficient signs.",
    "Impulse-propagation summaries apply a one-standard-deviation input and repeatedly apply the fitted transition matrix for 15 steps; they are descriptive model-based quantities and do not estimate intervention effects."
  )

  utils::write.csv(main, file.path(directory, "Supplementary_Table_11_source_data.csv"), row.names = FALSE)
  utils::write.csv(sex, file.path(directory, "Supplementary_Table_12_source_data.csv"), row.names = FALSE)
  save_publication_table(
    main,
    "Supplementary Table 11 | Addiction-status contrasts in posterior network metrics",
    paste(common_note, "The reported contrast is addicted minus non-addicted, averaged across the female and male contrasts."),
    file.path(directory, "Supplementary_Table_11.docx"),
    landscape = TRUE
  )
  save_publication_table(
    sex,
    "Supplementary Table 12 | Sex-specific addiction-status contrasts and their formal difference",
    paste(common_note, "The male-minus-female row is the formal comparison of the two addiction-status contrasts."),
    file.path(directory, "Supplementary_Table_12.docx"),
    landscape = TRUE
  )
}
