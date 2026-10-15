network_required_columns <- function(nodes = network_nodes) {
  c(
    "addiction", "gender",
    paste0(rep(nodes, each = 2), "_", rep(c("t1", "t2"), times = length(nodes)))
  )
}

map_integrated_retest_data <- function(data, nodes = network_nodes) {
  required <- c(
    "sex", "addiction_group_t1",
    "kinship_support_t1", "offline_friend_support_t1", "online_friend_support_t1",
    "kinship_support_t2", "offline_friend_support_t2", "online_friend_support_t2",
    paste0("gad", 1:7, "_item_t1"), paste0("gad", 1:7, "_item_t2"),
    paste0("mpiq", 1:8, "_item_t1"), paste0("mpiq", 1:8, "_item_t2")
  )
  absent <- setdiff(required, names(data))
  if (length(absent)) {
    stop("The retest file is missing: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  if (any(!data$addiction_group_t1 %in% c("Non-addicted", "Addicted"), na.rm = TRUE)) {
    stop("addiction_group_t1 must contain Non-addicted or Addicted.", call. = FALSE)
  }

  out <- data.frame(
    addiction = ifelse(data$addiction_group_t1 == "Addicted", "Addiction", "Non-Addiction"),
    gender = data$sex,
    stringsAsFactors = FALSE
  )

  out$kin_t1 <- data$kinship_support_t1
  out$offline_t1 <- data$offline_friend_support_t1
  out$online_t1 <- data$online_friend_support_t1
  out$kin_t2 <- data$kinship_support_t2
  out$offline_t2 <- data$offline_friend_support_t2
  out$online_t2 <- data$online_friend_support_t2

  for (i in 1:7) {
    out[[paste0("gad", i, "_t1")]] <- data[[paste0("gad", i, "_item_t1")]] + 1L
    out[[paste0("gad", i, "_t2")]] <- data[[paste0("gad", i, "_item_t2")]] + 1L
  }
  for (i in 1:8) {
    out[[paste0("mpiq", i, "_t1")]] <- data[[paste0("mpiq", i, "_item_t1")]]
    out[[paste0("mpiq", i, "_t2")]] <- data[[paste0("mpiq", i, "_item_t2")]]
  }

  if ("mpiq_total_t1" %in% names(data)) {
    expected_group <- ifelse(data$mpiq_total_t1 > 32, "Addiction", "Non-Addiction")
    if (!all(out$addiction == expected_group, na.rm = TRUE)) {
      stop("The baseline addiction group does not match MPIQ total score >32.", call. = FALSE)
    }
  }

  out
}

validate_network_data <- function(data, nodes = network_nodes) {
  data <- as.data.frame(data)
  required <- network_required_columns(nodes)
  absent <- setdiff(required, names(data))
  if (length(absent)) {
    stop("The network data are missing: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  if (any(!as.character(data$addiction) %in% c("Non-Addiction", "Addiction"), na.rm = TRUE) ||
      any(!as.character(data$gender) %in% c("Female", "Male"), na.rm = TRUE)) {
    stop("Use Non-Addiction/Addiction and Female/Male as the group labels.", call. = FALSE)
  }

  data <- data[, required, drop = FALSE]
  numeric_columns <- setdiff(required, c("addiction", "gender"))
  data[numeric_columns] <- lapply(data[numeric_columns], function(x) suppressWarnings(as.numeric(x)))

  data$addiction <- factor(
    data$addiction,
    levels = c("Non-Addiction", "Addiction")
  )
  data$gender <- factor(data$gender, levels = c("Female", "Male"))

  complete <- stats::complete.cases(data)
  if (!all(complete)) {
    message(sum(!complete), " rows with missing required values were excluded.")
    data <- data[complete, , drop = FALSE]
  }

  if (!nrow(data)) stop("No complete rows remain for the network model.", call. = FALSE)
  if (anyNA(data$addiction) || anyNA(data$gender)) {
    stop("Use Non-Addiction/Addiction and Female/Male as the group labels.", call. = FALSE)
  }

  support_columns <- paste0(c("kin", "offline", "online"), "_", rep(c("t1", "t2"), each = 3))
  gad_columns <- grep("^gad[1-7]_t[12]$", names(data), value = TRUE)
  mpiq_columns <- grep("^mpiq[1-8]_t[12]$", names(data), value = TRUE)

  if (any(vapply(data[support_columns], function(x) any(x < 0 | x > 15), logical(1)))) {
    stop("Support scores must range from 0 to 15.", call. = FALSE)
  }
  if (any(vapply(data[gad_columns], function(x) any(x < 1 | x > 4), logical(1)))) {
    stop("The model input uses GAD-7 item scores coded 1 to 4.", call. = FALSE)
  }
  if (any(vapply(data[mpiq_columns], function(x) any(x < 1 | x > 7), logical(1)))) {
    stop("MPIQ item scores must range from 1 to 7.", call. = FALSE)
  }

  group_counts <- table(data$addiction, data$gender)
  if (any(group_counts == 0L)) {
    stop(
      "Each addiction-by-sex group must be represented. Observed counts: ",
      paste(group_counts, collapse = ", "),
      call. = FALSE
    )
  }

  rownames(data) <- NULL
  data
}

read_network_data <- function(path = data_file) {
  if (!file.exists(path)) stop("Analysis data not found: ", path, call. = FALSE)
  data <- utils::read.csv(path, check.names = FALSE)

  if (all(network_required_columns() %in% names(data))) {
    ready <- data
  } else {
    ready <- map_integrated_retest_data(data)
  }

  validate_network_data(ready)
}
