required_packages <- c(
  "brms", "posterior", "rstan", "bayesplot", "qgraph", "ggplot2", "dplyr",
  "tidyr", "patchwork", "officer", "flextable"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop(
    "Install the required R packages before running the analysis: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}
