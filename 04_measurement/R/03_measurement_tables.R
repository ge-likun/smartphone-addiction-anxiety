required_packages <- c("flextable", "officer")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the following packages before creating the measurement tables: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

format_number <- function(x, digits = 3L) {
  ifelse(is.na(x), "—", formatC(x, digits = digits, format = "f"))
}

publication_table <- function(data) {
  flextable::flextable(data) |>
    flextable::theme_booktabs() |>
    flextable::font(fontname = "Times New Roman", part = "all") |>
    flextable::fontsize(size = 10, part = "all") |>
    flextable::bold(part = "header") |>
    flextable::align(j = 1, align = "left", part = "all") |>
    flextable::align(j = 2:ncol(data), align = "center", part = "all") |>
    flextable::valign(valign = "center", part = "all") |>
    flextable::autofit()
}

write_table <- function(data, title, note, path) {
  document <- officer::read_docx()
  document <- officer::body_add_fpar(
    document,
    officer::fpar(
      officer::ftext(
        title,
        officer::fp_text(font.family = "Times New Roman", font.size = 10, bold = TRUE)
      )
    )
  )
  document <- flextable::body_add_flextable(document, publication_table(data))
  document <- officer::body_add_fpar(
    document,
    officer::fpar(
      officer::ftext(
        paste0("Note. ", note),
        officer::fp_text(font.family = "Times New Roman", font.size = 10)
      )
    )
  )
  print(document, target = path)
}

table_3 <- data.frame(
  Scale = reliability$scale,
  Items = reliability$items,
  `Cronbach's alpha` = format_number(reliability$cronbach_alpha),
  `Standardized loading range` = ifelse(
    is.na(reliability$standardized_loading_minimum),
    "—",
    paste0(
      format_number(reliability$standardized_loading_minimum, 3L),
      "–",
      format_number(reliability$standardized_loading_maximum, 3L)
    )
  ),
  CR = format_number(reliability$composite_reliability),
  AVE = format_number(reliability$average_variance_extracted),
  check.names = FALSE
)

table_4 <- data.frame(
  Model = fit_summary$model,
  Estimator = fit_summary$estimator,
  n = format(fit_summary$n, big.mark = ",", scientific = FALSE),
  Converged = fit_summary$converged,
  `Chi-square` = format_number(fit_summary$chisq, 2L),
  df = format_number(fit_summary$df, 0L),
  CFI = format_number(fit_summary$cfi),
  TLI = format_number(fit_summary$tli),
  RMSEA = format_number(fit_summary$rmsea),
  `RMSEA 90% CI` = paste0(
    format_number(fit_summary$rmsea_ci_lower),
    "–",
    format_number(fit_summary$rmsea_ci_upper)
  ),
  SRMR = format_number(fit_summary$srmr),
  check.names = FALSE
)

write.csv(table_3, file.path(output_dir, "Supplementary_Table_3_source_data.csv"), row.names = FALSE)
write.csv(table_4, file.path(output_dir, "Supplementary_Table_4_source_data.csv"), row.names = FALSE)

write_table(
  table_3,
  "Supplementary Table 3 | Reliability and factor-loading summary for source-specific support",
  paste(
    "Standardized loadings are from the respecified complete-case maximum-likelihood CFA.",
    "CR, composite reliability; AVE, average variance extracted.",
    "CR and AVE are reported for the three latent subscales."
  ),
  file.path(output_dir, "Supplementary_Table_3.docx")
)

write_table(
  table_4,
  "Supplementary Table 4 | Fit of the initial and respecified three-factor support models",
  paste(
    "The respecified model included residual covariances for SS5–SS8, SS4–SS7, SS2–SS5, SS2–SS8 and SS1–SS7.",
    "ML, maximum likelihood; FIML, full-information maximum likelihood; WLSMV, diagonally weighted least squares with mean- and variance-adjusted test statistics.",
    "CFI, comparative fit index; TLI, Tucker–Lewis index; RMSEA, root mean square error of approximation; SRMR, standardized root mean square residual.",
    if (all(fit_summary$converged)) "All models converged." else "Some models did not converge; consult the convergence column."
  ),
  file.path(output_dir, "Supplementary_Table_4.docx")
)
