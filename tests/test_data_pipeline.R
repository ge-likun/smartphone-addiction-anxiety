if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.UTF-8"))
}
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- if (length(args)) dirname(dirname(normalizePath(sub("^--file=", "", args[1]), winslash = "/"))) else getwd()
source(file.path(root, "R", "project.R"))
test_root <- file.path(root, "results", "tests", "data_pipeline")
dir.create(test_root, recursive = TRUE, showWarnings = FALSE)
base_config <- read_project_config(root, "config.R.example")
set.seed(4401L)
for (n in c(96L, 137L)) {
  case <- file.path(test_root, paste0("n", n))
  dir.create(case, recursive = TRUE, showWarnings = FALSE)
  d <- data.frame(id = sprintf("S%08d", seq_len(n)), sex = rep(1:2, length.out=n), age = sample(18:25,n,TRUE),
                  bmi = runif(n,18,28), grade = sample(1:5,n,TRUE), major=sample(1:3,n,TRUE),
                  hukou=sample(1:2,n,TRUE), father=sample(1:4,n,TRUE), mother=sample(1:4,n,TRUE),
                  only=sample(1:2,n,TRUE), activity=sample(1:3,n,TRUE), met=runif(n,0,3000))
  for (i in 1:7) d[[paste0("GAD",i)]] <- sample(1:4,n,TRUE)
  for (i in 1:8) d[[paste0("MPIQ",i)]] <- sample(1:7,n,TRUE)
  support_traits <- matrix(rnorm(n * 3L), n, 3L)
  for (i in 1:9) {
    latent <- support_traits[, ceiling(i / 3)] + rnorm(n, sd = .45)
    d[[paste0("SS",i)]] <- as.integer(cut(latent, breaks = c(-Inf,-1.1,-.6,0,.6,1.1,Inf)))
  }
  zero_scores <- seq_len(n) %% 7L == 0L
  for (i in 1:7) d[[paste0("GAD",i)]][zero_scores] <- 1L
  # Generate observations in each sex-by-addiction group.
  low <- seq_len(n) %% 4 %in% c(0,1)
  for (i in 1:8) d[[paste0("MPIQ",i)]][low] <- sample(1:3,sum(low),TRUE)
  for (i in 1:8) d[[paste0("MPIQ",i)]][!low] <- sample(5:7,sum(!low),TRUE)
  d$SS3[n] <- -2L
  follow <- d[seq_len(floor(n*.7)), , drop=FALSE]
  follow$time <- "2022-05-01 10:00:00"
  repeat_record <- follow[1, , drop=FALSE]
  repeat_record$time <- "2022-06-01 10:00:00"
  repeat_record$GAD1 <- if (repeat_record$GAD1 == 4) 1 else 4
  follow <- rbind(follow, repeat_record)
  times <- data.frame(id=d$id,time="2021-05-01 10:00:00")
  write_csv(d,file.path(case,"baseline.csv"))
  write_csv(d,file.path(case,"survey_next.csv"))
  write_csv(follow,file.path(case,"followup.csv"))
  write_csv(times,file.path(case,"times.csv"))
  cfg <- base_config
  cfg$input_mode <- "raw"
  cfg$raw$surveys <- data.frame(file=file.path(case,c("baseline.csv","survey_next.csv")),sheet=1,year=c(2021L,2022L))
  cfg$raw$followups <- data.frame(file=file.path(case,"followup.csv"),sheet=1)
  cfg$raw$baseline_times <- data.frame(file=file.path(case,"times.csv"),sheet=1)
  cfg$raw$baseline_year <- 2021L; cfg$raw$followup_year <- 2022L
  cfg$raw$network_batch_rule <- "all"
  cfg$raw$id_column <- "id"; cfg$raw$timestamp_column <- "time"
  cfg$raw$time_id_column <- "id"; cfg$raw$time_column <- "time"
  cfg$raw$columns <- c(sex_code="sex",age_years="age",bmi_kg_m2="bmi",grade_code="grade",major_code="major",
                       hukou_code="hukou",father_education_code="father",mother_education_code="mother",
                       only_child_code="only",ipaq_category_code="activity",ipaq_total_met_min_week="met")
  cfg$raw$support_columns <- paste0("SS",1:9)
  cfg$raw$reported_totals <- character()
  cfg$paths$data <- file.path(case,"data"); cfg$paths$results <- file.path(case,"output")
  config_path <- file.path(case,"test_config.R")
  lines <- c("config <-", capture.output(dput(cfg)))
  writeLines(lines, config_path)
  Sys.setenv(PROJECT_ROOT=root,PROJECT_CONFIG_FILE=config_path,PROJECT_DATA_DIR=cfg$paths$data,PROJECT_RESULTS_DIR=cfg$paths$results)
  run_module(root,"00_data_preparation","run_data_preparation.R")
  a <- read.csv(file.path(cfg$paths$data,"analysis_data.csv"))
  r <- read.csv(file.path(cfg$paths$data,"retest_data.csv"))
  w <- read.csv(file.path(cfg$paths$data,"network_data.csv"))
  stopifnot(nrow(a)==2*n,nrow(r)==floor(n*.7),nrow(w)==nrow(r),!"raw_id" %in% names(a))
  stopifnot(sum(is.na(a$kinship_support)) == 2L, sum(is.na(a$social_support_total)) == 2L)
  linked <- match(r$retest_pair_id, w$retest_pair_id)
  stopifnot(any(r$gad1_item_t2 != w$gad1_item_t2[linked]), all(r$interval_days>0))
  cfg$raw$network_batch_rule <- "earliest"
  cfg$paths$data <- file.path(case, "data_earliest")
  cfg$paths$results <- file.path(case, "output_earliest")
  writeLines(c("config <-", capture.output(dput(cfg))), config_path)
  Sys.setenv(PROJECT_DATA_DIR = cfg$paths$data, PROJECT_RESULTS_DIR = cfg$paths$results)
  run_module(root, "00_data_preparation", "run_data_preparation.R")
  earliest <- read.csv(file.path(cfg$paths$data, "network_data.csv"))
  linked_earliest <- match(r$retest_pair_id, earliest$retest_pair_id)
  stopifnot(identical(r$gad1_item_t2, earliest$gad1_item_t2[linked_earliest]))
  cfg$paths$data <- file.path(case, "data")
  old <- getwd(); setwd(file.path(root,"03_bayesian_network"))
  project_root <- getwd(); source("01_config.R"); source("02_data.R")
  checked <- read_network_data(file.path(cfg$paths$data,"network_data.csv"))
  stopifnot(nrow(checked)==nrow(w))
  setwd(old)
  shared_cross <- a[, setdiff(names(a), c("participant_id", "retest_pair_id", "ipaq_total_met_min_week",
    "sex", "grade_manuscript", "addiction_group_code", "addiction_group", "gad7_symptomatic_code")), drop = FALSE]
  shared_cross$record_id <- a$participant_id
  paired_columns <- c("kinship_support", "offline_friend_support", "online_friend_support",
                      paste0("gad", 1:7, "_item"), paste0("mpiq", 1:8, "_item"))
  shared_long <- w[, c("baseline_year", "followup_year",
    as.vector(outer(paired_columns, c("t1", "t2"), paste, sep = "_"))), drop = FALSE]
  shared_long$baseline_record_id <- w$baseline_participant_id
  retest_index <- match(w$retest_pair_id, r$retest_pair_id)
  shared_long$retest_interval_months <- r$interval_months[retest_index]
  for (name in c("kinship_support", "offline_friend_support", "online_friend_support")) {
    shared_long[[paste0("retest_", name, "_t2")]] <- r[[paste0(name, "_t2")]][retest_index]
  }
  write_csv(shared_cross, file.path(case, "cross_sectional.csv"))
  write_csv(shared_long, file.path(case, "longitudinal.csv"))
  cfg$input_mode <- "shared"
  cfg$shared <- list(cross_sectional = file.path(case, "cross_sectional.csv"),
                     longitudinal = file.path(case, "longitudinal.csv"))
  cfg$paths$data <- file.path(case, "data_shared")
  cfg$paths$results <- file.path(case, "output_shared")
  writeLines(c("config <-", capture.output(dput(cfg))), config_path)
  Sys.setenv(PROJECT_DATA_DIR = cfg$paths$data, PROJECT_RESULTS_DIR = cfg$paths$results)
  run_module(root, "00_data_preparation", "run_shared_data_preparation.R")
  shared_checked <- read_network_data(file.path(cfg$paths$data, "network_data.csv"))
  shared_retest <- read.csv(file.path(cfg$paths$data, "retest_data.csv"))
  retest_values <- c("interval_months", paste0(c("kinship_support", "offline_friend_support",
                      "online_friend_support", "social_support_total"), "_t1"),
                     paste0(c("kinship_support", "offline_friend_support", "online_friend_support", "social_support_total"), "_t2"))
  stopifnot(identical(checked, shared_checked),
            isTRUE(all.equal(r[retest_index, retest_values], shared_retest[retest_values], tolerance = 0,
                            check.attributes = FALSE)))
}
cat("Raw and shared-data scoring/linkage checks passed for different sample sizes and survey years.\n")
