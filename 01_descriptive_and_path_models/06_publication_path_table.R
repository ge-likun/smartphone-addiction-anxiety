binary_publication_path_table <- function(estimates) {
  check_columns(estimates, c("support", "lhs", "op", "label", "est", "std.all", "ci.lower", "ci.upper"),
                "binary parameter estimates")
  fields <- c("Quantity", "Description", "Sex or comparison", "B or product", "Standardized beta", "95% BC bootstrap CI")
  rows <- list()
  header <- function(label) {
    row <- as.data.frame(setNames(rep(list(""), length(fields)), fields), check.names = FALSE)
    row[[1]] <- label
    rows[[length(rows) + 1L]] <<- row
  }
  result <- function(quantity, description, sex, value, digits, standardized = FALSE) {
    row <- data.frame(Quantity = quantity, Description = description, `Sex or comparison` = sex,
      `B or product` = sprintf(paste0("%.", digits, "f"), value$est),
      `Standardized beta` = if (standardized) sprintf("%.3f", value$std.all) else "",
      `95% BC bootstrap CI` = paste(sprintf(paste0("%.", digits, "f"), value$ci.lower), "to",
                                   sprintf(paste0("%.", digits, "f"), value$ci.upper)), check.names = FALSE)
    rows[[length(rows) + 1L]] <<- row
  }
  supports <- unique(estimates$support)
  definitions <- list(
    a = c("Addicted versus non-addicted -> support", "a", "diff_a"),
    b1 = c("Support -> GAD-7, non-addicted", "b1", "diff_b1"),
    b2 = c("Addicted minus non-addicted support slope", "b2", "diff_b2"),
    "c'" = c("Addicted versus non-addicted -> GAD-7", "cp", "diff_cp"))
  products <- list(
    "b1 + b2" = c("Support -> GAD-7, addicted", "b_add", 4L),
    "a x b1" = c("Association product in the non-addicted group", "product_nonadd", 6L),
    "a x (b1 + b2)" = c("Association product in the addicted group", "product_add", 6L),
    "a x b2" = c("Change in the association product", "index", 6L))
  for (support in supports) {
    header(support)
    data <- estimates[estimates$support == support, ]
    lookup <- function(name, path = FALSE) {
      value <- if (path) data[data$op == "~" & data$label == name, ]
        else data[data$op == ":=" & data$lhs == name, ]
      if (nrow(value) != 1L) stop("Missing or duplicate binary quantity: ", support, "/", name)
      value
    }
    for (quantity in names(definitions)) {
      definition <- definitions[[quantity]]
      for (sex in c("Female", "Male")) {
        suffix <- if (sex == "Female") "f" else "m"
        result(quantity, definition[1], sex, lookup(paste0(definition[2], "_", suffix), TRUE), 4L, TRUE)
      }
      result(quantity, definition[1], "Female - male", lookup(definition[3]), 4L)
    }
    for (quantity in names(products)) {
      definition <- products[[quantity]]
      for (suffix in c("f", "m", "diff")) {
        sex <- c(f = "Female", m = "Male", diff = "Female - male")[[suffix]]
        result(quantity, definition[1], sex, lookup(paste0(definition[2], "_", suffix)), as.integer(definition[3]))
      }
    }
  }
  do.call(rbind, rows)
}
combine_publication_path_tables <- function(continuous, binary) {
  fields <- c("Quantity", "Description", "Sex or comparison", "B or product",
              "Standardized beta", "95% BC bootstrap CI")
  check_columns(continuous, fields, "continuous path table")
  check_columns(binary, fields, "binary path table")
  section <- function(label) {
    result <- as.data.frame(setNames(rep(list(""), length(fields)), fields), check.names = FALSE)
    result[[1]] <- label
    result
  }
  result <- rbind(section("Continuous MPIQ-score models"), continuous[fields],
                  section("Binary addiction-status sensitivity models"), binary[fields])
  symbols <- c(b1 = "b\u2081", b2 = "b\u2082", "c'" = "c\u2032",
               "a x b2" = "a \u00d7 b\u2082", "a(b1 + b2X)" = "a(b\u2081 + b\u2082X)",
               "b1 + b2" = "b\u2081 + b\u2082", "a x b1" = "a \u00d7 b\u2081",
               "a x (b1 + b2)" = "a \u00d7 (b\u2081 + b\u2082)")
  mapped <- match(result$Quantity, names(symbols))
  result$Quantity[!is.na(mapped)] <- unname(symbols[mapped[!is.na(mapped)]])
  result$Description <- gsub("->", "\u2192", result$Description, fixed = TRUE)
  result$Description <- gsub(" x ", " \u00d7 ", result$Description, fixed = TRUE)
  result$row_type <- ifelse(result[["B or product"]] == "", "header", "body")
  result[c("row_type", fields)]
}
write_publication_path_table <- function() {
  continuous_dir <- Sys.getenv("MEDIATION_OUTPUT_DIR", file.path("..", "results", "01_path_models"))
  binary_dir <- Sys.getenv("BINARY_PATH_OUTPUT_DIR", file.path("..", "results", "01_binary_path_models"))
  read_table <- function(directory, file) {
    path <- file.path(directory, file)
    if (!file.exists(path)) stop("Missing path results: ", path)
    read.csv(path, check.names = FALSE, colClasses = "character", na.strings = character())
  }
  continuous_status <- read_table(continuous_dir, "run_status.csv")
  binary_status <- read_table(binary_dir, "run_status.csv")
  if (!identical(sort(unique(continuous_status$n)), sort(unique(binary_status$n))) ||
      !identical(sort(unique(continuous_status$bootstrap_replicates)), sort(unique(binary_status$bootstrap_replicates))) ||
      any(continuous_status$converged != "TRUE") || any(binary_status$converged != "TRUE"))
    stop("Continuous and binary path runs have incompatible sample sizes, bootstrap settings or convergence status.")
  continuous <- read_table(continuous_dir, "Supplementary_Table_5_source_data.csv")
  estimate_files <- list.files(binary_dir, pattern = "_parameter_estimates[.]csv$", full.names = TRUE)
  if (!length(estimate_files)) stop("Missing binary parameter-estimate files.")
  binary_estimates <- do.call(rbind, lapply(estimate_files, read.csv, check.names = FALSE))
  binary <- binary_publication_path_table(binary_estimates)
  support_order <- continuous$Quantity[continuous[["B or product"]] == ""]
  binary <- do.call(rbind, lapply(support_order, function(support) {
    start <- which(binary$Quantity == support)
    if (length(start) != 1L) stop("Continuous and binary support sources differ: ", support)
    next_header <- which(binary[["B or product"]] == "" & seq_len(nrow(binary)) > start)
    end <- if (length(next_header)) min(next_header) - 1L else nrow(binary)
    binary[start:end, , drop = FALSE]
  }))
  result <- combine_publication_path_tables(continuous, binary)
  write.csv(result[-1], file.path(continuous_dir, "Supplementary_Table_5_combined_source_data.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  save_table_docx(result, "Supplementary Table 5. Sex-grouped path-model associations and binary addiction-status sensitivity",
    paste("Continuous MPIQ and support scores were mean-centred. Binary addiction status was coded 0 for scores of 32 or lower and 1 for higher scores.",
          "Survey year was included in both equations. B and product terms are unstandardized; standardized beta is reported for sex-specific paths.",
          paste0("Intervals are bias-corrected bootstrap intervals based on ", unique(continuous_status$bootstrap_replicates), " resamples."),
          "Products summarize associations measured at the same assessment."),
    file.path(continuous_dir, "Supplementary_Table_5_combined_publication.docx"))
  invisible(result)
}
