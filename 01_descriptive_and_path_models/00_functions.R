required_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install the following packages before running this script: ",
         paste(missing, collapse = ", "))
  }
}

check_columns <- function(data, columns, data_name = "data") {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(data_name, " is missing: ", paste(missing, collapse = ", "))
  }
}

normalise_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.0$", "", x)
  x[x == ""] <- NA_character_
  x
}

format_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

format_mean_sd <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("")
  sprintf("%.2f (%.2f)", mean(x), stats::sd(x))
}

format_n_pct <- function(x, level) {
  keep <- !is.na(x)
  denominator <- sum(keep)
  numerator <- sum(x[keep] == level)
  sprintf("%s (%.1f%%)", format(numerator, big.mark = ","), 100 * numerator / denominator)
}

cohens_d <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  nx <- length(x)
  ny <- length(y)
  pooled_sd <- sqrt(
    ((nx - 1) * stats::var(x) + (ny - 1) * stats::var(y)) /
      (nx + ny - 2)
  )
  (mean(x) - mean(y)) / pooled_sd
}

cramers_v <- function(x, group) {
  keep <- stats::complete.cases(x, group)
  x <- droplevels(as.factor(x[keep]))
  group <- droplevels(as.factor(group[keep]))
  tab <- table(x, group)
  test <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  n <- sum(tab)
  k <- min(nrow(tab) - 1, ncol(tab) - 1)
  sqrt(unname(test$statistic) / (n * k))
}

continuous_tost <- function(x, y, bound = 0.10, alpha = 0.05) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  nx <- length(x)
  ny <- length(y)
  pooled_sd <- sqrt(((nx - 1) * stats::var(x) + (ny - 1) * stats::var(y)) /
                    (nx + ny - 2))
  difference <- mean(x) - mean(y)
  standard_error <- pooled_sd * sqrt(1 / nx + 1 / ny)
  degrees_freedom <- nx + ny - 2
  lower_bound <- -bound * pooled_sd
  upper_bound <- bound * pooled_sd
  p_lower <- stats::pt((difference - lower_bound) / standard_error,
                       df = degrees_freedom, lower.tail = FALSE)
  p_upper <- stats::pt((difference - upper_bound) / standard_error,
                       df = degrees_freedom, lower.tail = TRUE)
  critical <- stats::qt(1 - alpha, df = degrees_freedom)
  ci_raw <- difference + c(-1, 1) * critical * standard_error
  data.frame(
    effect_measure = "Cohen's d",
    estimate = difference / pooled_sd,
    interval_lower = ci_raw[1] / pooled_sd,
    interval_upper = ci_raw[2] / pooled_sd,
    confidence_level = 0.90,
    interval_type = "Two-sided",
    equivalence_bound = bound,
    equivalence_p = max(p_lower, p_upper),
    equivalent = max(p_lower, p_upper) < alpha
  )
}

categorical_equivalence <- function(x, group, bound = 0.10, alpha = 0.05) {
  keep <- stats::complete.cases(x, group)
  x <- droplevels(as.factor(x[keep]))
  group <- droplevels(as.factor(group[keep]))
  tab <- table(x, group)
  test <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  chi_square <- unname(test$statistic)
  degrees_freedom <- unname(test$parameter)
  n <- sum(tab)
  k <- min(nrow(tab) - 1, ncol(tab) - 1)
  estimate <- sqrt(chi_square / (n * k))
  boundary_ncp <- n * k * bound^2
  p_equivalence <- stats::pchisq(
    chi_square, df = degrees_freedom, ncp = boundary_ncp, lower.tail = TRUE
  )
  objective <- function(ncp) {
    stats::pchisq(chi_square, df = degrees_freedom, ncp = ncp) - alpha
  }
  upper_ncp <- NA_real_
  if (objective(0) > 0) {
    upper_limit <- 1
    while (objective(upper_limit) > 0 && upper_limit < 1e9) {
      upper_limit <- upper_limit * 2
    }
    if (objective(upper_limit) <= 0) {
      upper_ncp <- stats::uniroot(objective, interval = c(0, upper_limit))$root
    }
  }
  data.frame(
    effect_measure = "Cram\u00e9r's V",
    estimate = estimate,
    interval_lower = 0,
    interval_upper = sqrt(upper_ncp / (n * k)),
    confidence_level = 0.95,
    interval_type = "One-sided upper",
    equivalence_bound = bound,
    equivalence_p = if (is.finite(upper_ncp)) p_equivalence else NA_real_,
    equivalent = is.finite(upper_ncp) && p_equivalence < alpha
  )
}

run_equivalence_tests <- function(data, group, group_levels, variables,
                                  bound_d = 0.10, bound_v = 0.10,
                                  alpha = 0.05) {
  group <- factor(group, levels = group_levels)
  rows <- lapply(variables, function(spec) {
    value <- data[[spec$variable]]
    if (spec$type == "continuous") {
      x <- value[group == group_levels[1]]
      y <- value[group == group_levels[2]]
      result <- continuous_tost(x, y, bound = bound_d, alpha = alpha)
      difference_p <- stats::wilcox.test(x, y, exact = FALSE)$p.value
    } else {
      result <- categorical_equivalence(value, group, bound = bound_v, alpha = alpha)
      keep <- stats::complete.cases(value, group)
      test_table <- table(
        droplevels(as.factor(value[keep])),
        droplevels(as.factor(group[keep]))
      )
      difference_p <- suppressWarnings(
        stats::chisq.test(test_table, correct = FALSE)$p.value
      )
    }
    data.frame(
      variable = spec$variable,
      Characteristic = spec$label,
      result,
      difference_p = difference_p,
      dagger = isTRUE(result$equivalent),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

add_equivalence_markers <- function(table_data, equivalence_results) {
  marked <- equivalence_results[equivalence_results$dagger, , drop = FALSE]
  for (i in seq_len(nrow(marked))) {
    row <- which(table_data$Characteristic == marked$Characteristic[i])
    if (length(row) == 1L && nzchar(table_data$`Effect size`[row])) {
      table_data$`Effect size`[row] <- paste0(table_data$`Effect size`[row], "\u2020")
    }
  }
  table_data
}

comparison_table <- function(data, group, group_levels, variables, overall_label) {
  group <- factor(group, levels = group_levels)
  output <- list()
  row_index <- 0L
  for (spec in variables) {
    value <- data[[spec$variable]]
    if (spec$type == "continuous") {
      x <- value[group == group_levels[1]]
      y <- value[group == group_levels[2]]
      row_index <- row_index + 1L
      output[[row_index]] <- data.frame(
        row_type = "body",
        Characteristic = spec$label,
        Overall = format_mean_sd(value),
        Group1 = format_mean_sd(x),
        Group2 = format_mean_sd(y),
        `Effect size` = sprintf("d = %.2f", cohens_d(x, y)),
        `P value` = format_p(stats::wilcox.test(x, y, exact = FALSE)$p.value),
        check.names = FALSE
      )
    } else {
      effect <- cramers_v(value, group)
      keep <- stats::complete.cases(value, group)
      test_table <- table(droplevels(as.factor(value[keep])), droplevels(as.factor(group[keep])))
      p_value <- suppressWarnings(stats::chisq.test(test_table, correct = FALSE)$p.value)
      if (!is.null(spec$single_level)) {
        row_index <- row_index + 1L
        output[[row_index]] <- data.frame(
          row_type = "body",
          Characteristic = spec$label,
          Overall = format_n_pct(value, spec$single_level),
          Group1 = format_n_pct(value[group == group_levels[1]], spec$single_level),
          Group2 = format_n_pct(value[group == group_levels[2]], spec$single_level),
          `Effect size` = sprintf("V = %.2f", effect),
          `P value` = format_p(p_value),
          check.names = FALSE
        )
      } else {
        row_index <- row_index + 1L
        output[[row_index]] <- data.frame(
          row_type = "header",
          Characteristic = spec$label,
          Overall = "", Group1 = "", Group2 = "",
          `Effect size` = sprintf("V = %.2f", effect),
          `P value` = format_p(p_value),
          check.names = FALSE
        )
        for (level in spec$levels) {
          row_index <- row_index + 1L
          level_label <- if (!is.null(spec$level_labels)) spec$level_labels[[as.character(level)]] else as.character(level)
          output[[row_index]] <- data.frame(
            row_type = "subrow",
            Characteristic = paste0("  ", level_label),
            Overall = format_n_pct(value, level),
            Group1 = format_n_pct(value[group == group_levels[1]], level),
            Group2 = format_n_pct(value[group == group_levels[2]], level),
            `Effect size` = "", `P value` = "",
            check.names = FALSE
          )
        }
      }
    }
  }
  result <- do.call(rbind, output)
  names(result)[names(result) == "Overall"] <- overall_label
  names(result)[names(result) == "Group1"] <- group_levels[1]
  names(result)[names(result) == "Group2"] <- group_levels[2]
  result
}

save_table_docx <- function(table_data, title, note, path) {
  required_packages(c("flextable", "officer"))
  display <- table_data[, setdiff(names(table_data), "row_type"), drop = FALSE]
  ft <- flextable::flextable(display)
  ft <- flextable::font(ft, fontname = "Times New Roman", part = "all")
  ft <- flextable::fontsize(ft, size = 10, part = "all")
  ft <- flextable::bold(ft, part = "header")
  header_rows <- which(table_data$row_type == "header")
  if (length(header_rows)) ft <- flextable::bold(ft, i = header_rows, bold = TRUE, part = "body")
  ft <- flextable::align(ft, j = 1, align = "left", part = "all")
  if (ncol(display) > 1) ft <- flextable::align(ft, j = 2:ncol(display), align = "center", part = "all")
  ft <- flextable::border_remove(ft)
  black <- officer::fp_border(color = "black", width = 1)
  ft <- flextable::hline_top(ft, border = black, part = "header")
  ft <- flextable::hline(ft, i = 1, border = black, part = "header")
  ft <- flextable::hline_bottom(ft, border = black, part = "body")
  ft <- flextable::autofit(ft)
  ft <- flextable::add_footer_lines(ft, values = note)
  ft <- flextable::font(ft, fontname = "Times New Roman", part = "footer")
  ft <- flextable::fontsize(ft, size = 10, part = "footer")
  document <- officer::read_docx()
  title_format <- officer::fp_text(
    font.family = "Times New Roman", font.size = 10, bold = TRUE
  )
  document <- officer::body_add_fpar(
    document,
    officer::fpar(officer::ftext(title, title_format))
  )
  document <- flextable::body_add_flextable(document, value = ft)
  staged_file <- tempfile(pattern = "publication-table-", fileext = ".docx")
  on.exit(unlink(staged_file), add = TRUE)
  print(document, target = staged_file)
  if (!file.copy(staged_file, path, overwrite = TRUE)) {
    stop("Could not save table: ", path)
  }
}
