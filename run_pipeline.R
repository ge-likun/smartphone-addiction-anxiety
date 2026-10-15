if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.UTF-8"))
}
pipeline_root <- local({
  paths <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(paths)) dirname(normalizePath(paths[[1]], winslash = "/")) else {
    a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash = "/")) else getwd()
  }
})
source(file.path(pipeline_root, "R", "project.R"), encoding = "UTF-8")
run_pipeline <- function(stage = "plan", config_file = "local_config.R", ordinal_targets = "all") {
  cfg <- read_project_config(pipeline_root, config_file)
  input_mode <- if (is.null(cfg$input_mode)) "raw" else cfg$input_mode
  if (!input_mode %in% c("shared", "raw")) stop("input_mode must be shared or raw.")
  stages <- c("prepare", "descriptive", "measurement", "binary", "path_tables", "gamlss", "gaussian_fit",
              "gaussian_results", "gaussian_ppc", "ordinal_fit", "ordinal_results")
  if (!stage %in% c("plan", "all", stages)) stop("Unknown stage: ", stage)
  data_dir <- absolute_path(cfg$paths$data, pipeline_root)
  results <- absolute_path(cfg$paths$results, pipeline_root)
  env <- c(PROJECT_ROOT = pipeline_root,
           PROJECT_CONFIG_FILE = absolute_path(config_file, pipeline_root),
           PROJECT_DATA_DIR = data_dir, PROJECT_RESULTS_DIR = results,
           ANALYSIS_DATA_FILE = file.path(data_dir, "analysis_data.csv"),
           RETEST_DATA_FILE = file.path(data_dir, "retest_data.csv"),
           GAMLSS_DATA_FILE = file.path(data_dir, "analysis_data.csv"),
           NETWORK_DATA_FILE = file.path(data_dir, "network_data.csv"),
           NETWORK_MODEL_FILE = file.path(results, "03_bayesian_network", "models", "gaussian.rds"),
           DESCRIPTIVE_OUTPUT_DIR = file.path(results, "01_descriptive_tables"),
           MEDIATION_OUTPUT_DIR = file.path(results, "01_path_models"),
           RETEST_OUTPUT_DIR = file.path(results, "01_retest_stability"),
           BINARY_PATH_OUTPUT_DIR = file.path(results, "01_binary_path_models"),
           GAMLSS_OUTPUT_DIR = file.path(results, "02_gamlss"),
           NETWORK_OUTPUT_DIR = file.path(results, "03_bayesian_network"),
           CFA_OUTPUT_DIR = file.path(results, "04_measurement", "cfa"),
           BOOTSTRAP_REPS = cfg$bootstrap_reps,
           NETWORK_CHAINS = cfg$network$chains, NETWORK_ITERATIONS = cfg$network$iter,
           NETWORK_WARMUP = cfg$network$warmup, NETWORK_CORES = cfg$network$cores,
           NETWORK_THREADS_PER_CHAIN = if (is.null(cfg$network$threads_per_chain)) 1L else cfg$network$threads_per_chain,
           NETWORK_ADAPT_DELTA = cfg$network$adapt_delta,
           NETWORK_MAX_TREEDEPTH = cfg$network$max_treedepth,
           NETWORK_IMPULSE_STEPS = cfg$network$impulse_steps,
           NETWORK_DISPLAY_THRESHOLD = cfg$network$display_threshold,
           NETWORK_SEED = cfg$network$seed,
           NETWORK_DISPLAY_SCALE = if (is.null(cfg$network$display_scale)) "raw" else cfg$network$display_scale,
           BASELINE_YEAR = cfg$raw$baseline_year, FOLLOWUP_YEAR = cfg$raw$followup_year,
           GAMLSS_RESIDUAL_REPLICATES = cfg$gamlss_residual_replicates,
           SAVE_NETWORK_POSTERIOR = "true")
  previous <- Sys.getenv(names(env), unset = NA_character_)
  on.exit({
    Sys.unsetenv(names(previous)[is.na(previous)])
    restored <- previous[!is.na(previous)]
    if (length(restored)) do.call(Sys.setenv, as.list(restored))
  })
  do.call(Sys.setenv, as.list(env))
  input_files <- if (input_mode == "shared") c(cfg$shared$cross_sectional, cfg$shared$longitudinal) else
    c(cfg$raw$surveys$file, cfg$raw$followups$file, cfg$raw$baseline_times$file)
  if (stage == "plan") {
    cat("Pipeline stages:\n", paste(stages, collapse = " -> "), "\n\n", sep = "")
    print(data.frame(file = input_files, exists = file.exists(vapply(input_files, absolute_path, character(1), root = pipeline_root))), row.names = FALSE)
    cat("\nResults directory:", results, "\nNo models were fitted.\n")
    return(invisible(stages))
  }
  dir.create(results, recursive = TRUE, showWarnings = FALSE)
  saveRDS(cfg, file.path(results, "run_configuration.rds"))
  writeLines(capture.output(sessionInfo()), file.path(results, "sessionInfo_pipeline.txt"))
  for (s in if (stage == "all") stages else stage) {
    message("Running stage: ", s)
    switch(s,
      prepare = run_module(pipeline_root, "00_data_preparation", if (input_mode == "shared")
        "run_shared_data_preparation.R" else "run_data_preparation.R"),
      descriptive = run_module(pipeline_root, "01_descriptive_and_path_models", "run_all.R"),
      binary = run_module(pipeline_root, "01_descriptive_and_path_models", "run_binary_sensitivity.R"),
      path_tables = run_module(pipeline_root, "01_descriptive_and_path_models", "run_publication_path_table.R"),
      measurement = run_module(pipeline_root, "04_measurement", "run_all.R"),
      gamlss = run_module(pipeline_root, "02_gamlss", "run_all.R"),
      gaussian_fit = run_module(pipeline_root, "03_bayesian_network", "run_fit_model.R"),
      gaussian_results = run_module(pipeline_root, "03_bayesian_network", "run_postprocessing.R"),
      gaussian_ppc = run_module(pipeline_root, "03_bayesian_network", "run_predictive_checks.R"),
      ordinal_fit = run_module(pipeline_root, "03_bayesian_network", "run_ordinal.R", c("fit", ordinal_targets)),
      ordinal_results = run_module(pipeline_root, "03_bayesian_network", "run_ordinal.R", "summarize"))
  }
  invisible(results)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  run_pipeline(if (length(args)) args[1] else "plan",
               if (length(args) >= 2L) args[2] else "local_config.R",
               if (length(args) >= 3L) args[3] else "all")
}
