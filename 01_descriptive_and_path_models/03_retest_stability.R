retest_file <- Sys.getenv(
  "RETEST_DATA_FILE",
  unset = file.path("..", "retest_data.csv")
)
output_dir <- Sys.getenv(
  "RETEST_OUTPUT_DIR",
  unset = file.path("..", "results", "01_retest_stability")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(retest_file)) stop("Retest data not found: ", retest_file)
dat <- read.csv(retest_file, check.names = FALSE)
required <- c(
  "retest_pair_id", "interval_months",
  "kinship_support_t1", "kinship_support_t2",
  "offline_friend_support_t1", "offline_friend_support_t2",
  "online_friend_support_t1", "online_friend_support_t2",
  "social_support_total_t1", "social_support_total_t2"
)
check_columns(dat, required, "retest_data.csv")
if (nrow(dat) == 0L || anyDuplicated(dat$retest_pair_id)) {
  stop("The retest dataset must contain unique linked participants.")
}

pairs <- list(
  kinship = c("kinship_support_t1", "kinship_support_t2"),
  offline = c("offline_friend_support_t1", "offline_friend_support_t2"),
  online = c("online_friend_support_t1", "online_friend_support_t2"),
  total_support = c("social_support_total_t1", "social_support_total_t2")
)

results <- do.call(rbind, lapply(names(pairs), function(measure) {
  variables <- pairs[[measure]]
  keep <- stats::complete.cases(dat[, variables])
  n <- sum(keep)
  r <- if (n > 1L) stats::cor(dat[[variables[1]]][keep], dat[[variables[2]]][keep]) else NA_real_
  z <- atanh(r)
  se <- if (n > 3L) 1 / sqrt(n - 3) else NA_real_
  data.frame(
    measure = measure,
    complete_pairs = n,
    pearson_r = r,
    ci_lower = tanh(z - stats::qnorm(0.975) * se),
    ci_upper = tanh(z + stats::qnorm(0.975) * se)
  )
}))

valid_intervals <- dat$interval_months[is.finite(dat$interval_months)]
interval_summary <- data.frame(
  n = length(valid_intervals),
  mean_months = if (length(valid_intervals)) mean(valid_intervals) else NA_real_,
  sd_months = stats::sd(valid_intervals),
  minimum_months = if (length(valid_intervals)) min(valid_intervals) else NA_real_,
  maximum_months = if (length(valid_intervals)) max(valid_intervals) else NA_real_
)

write.csv(results, file.path(output_dir, "retest_stability.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(interval_summary, file.path(output_dir, "retest_interval.csv"), row.names = FALSE, fileEncoding = "UTF-8")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Retest summaries written to ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
