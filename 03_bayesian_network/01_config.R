network_nodes <- c(
  "kin", "offline", "online",
  paste0("gad", 1:7),
  paste0("mpiq", 1:8)
)

network_node_labels <- c(
  "Kin", "Offline", "Online",
  paste0("G", 1:7),
  paste0("MP", 1:8)
)

network_communities <- list(
  Support = c("kin", "offline", "online"),
  Anxiety = paste0("gad", 1:7),
  Smartphone_addiction = paste0("mpiq", 1:8)
)

network_groups <- c(
  nonaddicted_female = "Non-addicted women",
  addicted_female = "Addicted women",
  nonaddicted_male = "Non-addicted men",
  addicted_male = "Addicted men"
)

network_seed <- 123L
layout_seed <- 20260913L
edge_display_threshold <- 0.10
impulse_steps <- 15L
minimum_group_n <- as.integer(Sys.getenv("NETWORK_MINIMUM_GROUP_N", unset = "20"))
network_chains <- as.integer(Sys.getenv("NETWORK_CHAINS", unset = "4"))
network_iterations <- as.integer(Sys.getenv("NETWORK_ITERATIONS", unset = "4000"))
network_warmup <- as.integer(Sys.getenv("NETWORK_WARMUP", unset = "2000"))
network_threads_per_chain <- as.integer(Sys.getenv("NETWORK_THREADS_PER_CHAIN", unset = "1"))
network_adapt_delta <- as.numeric(Sys.getenv("NETWORK_ADAPT_DELTA", unset = "0.90"))
network_max_treedepth <- as.integer(Sys.getenv("NETWORK_MAX_TREEDEPTH", unset = "10"))

if (!is.finite(minimum_group_n) || minimum_group_n < 1L) stop("NETWORK_MINIMUM_GROUP_N must be a positive integer.")
if (!is.finite(network_chains) || network_chains < 1L) stop("NETWORK_CHAINS must be a positive integer.")
if (!is.finite(network_iterations) || network_iterations < 2L) stop("NETWORK_ITERATIONS must be at least 2.")
if (!is.finite(network_warmup) || network_warmup < 1L || network_warmup >= network_iterations) {
  stop("NETWORK_WARMUP must be positive and smaller than NETWORK_ITERATIONS.")
}
if (!is.finite(network_threads_per_chain) || network_threads_per_chain < 1L) {
  stop("NETWORK_THREADS_PER_CHAIN must be a positive integer.")
}
if (!is.finite(network_adapt_delta) || network_adapt_delta <= 0 || network_adapt_delta >= 1) {
  stop("NETWORK_ADAPT_DELTA must be between 0 and 1.")
}
if (!is.finite(network_max_treedepth) || network_max_treedepth < 1L) {
  stop("NETWORK_MAX_TREEDEPTH must be a positive integer.")
}

data_file <- chartr("\\", "/", Sys.getenv(
  "NETWORK_DATA_FILE",
  unset = file.path(project_root, "..", "network_data.csv")
))

model_file <- chartr("\\", "/", Sys.getenv(
  "NETWORK_MODEL_FILE",
  unset = file.path(project_root, "brms_model.rds")
))

output_dir <- chartr("\\", "/", Sys.getenv(
  "NETWORK_OUTPUT_DIR",
  unset = file.path(project_root, "..", "results", "03_bayesian_network")
))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
