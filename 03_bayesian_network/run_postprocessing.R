arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments, value = TRUE)
frame_files <- vapply(
  sys.frames(),
  function(frame) if (is.null(frame$ofile)) NA_character_ else frame$ofile,
  character(1)
)
script_file <- if (length(file_argument)) {
  sub("^--file=", "", file_argument[1])
} else if (any(!is.na(frame_files))) {
  tail(frame_files[!is.na(frame_files)], 1)
} else {
  NA_character_
}
project_root <- if (is.na(script_file)) "." else dirname(script_file)

source(file.path(project_root, "00_packages.R"), encoding = "UTF-8")
source(file.path(project_root, "01_config.R"), encoding = "UTF-8")
source(file.path(project_root, "02_data.R"), encoding = "UTF-8")
source(file.path(project_root, "03_model.R"), encoding = "UTF-8")
source(file.path(project_root, "04_diagnostics.R"), encoding = "UTF-8")
source(file.path(project_root, "05_networks.R"), encoding = "UTF-8")
source(file.path(project_root, "06_metrics.R"), encoding = "UTF-8")
source(file.path(project_root, "07_figures.R"), encoding = "UTF-8")
source(file.path(project_root, "08_tables.R"), encoding = "UTF-8")

fit <- load_network_model(model_file)
validate_network_model(fit)

diagnostics <- model_diagnostics(fit)
write_diagnostic_table(diagnostics, output_dir)

edge_draws <- extract_edge_draws(fit)
maximum_draws <- as.integer(Sys.getenv("NETWORK_MAX_DRAWS", unset = as.character(nrow(edge_draws))))
if (!is.finite(maximum_draws) || maximum_draws < 1L) stop("NETWORK_MAX_DRAWS must be a positive integer.", call. = FALSE)
edge_draws <- edge_draws[seq_len(min(maximum_draws, nrow(edge_draws))), , drop = FALSE]

raw_arrays <- reconstruct_network_arrays(edge_draws)
standardized_arrays <- standardize_network_arrays(raw_arrays, fit$data)
raw_means <- posterior_mean_networks(raw_arrays)
standardized_means <- posterior_mean_networks(standardized_arrays)

for (group in names(raw_means)) {
  utils::write.csv(
    raw_means[[group]],
    file.path(output_dir, paste0("mean_network_raw_", group, ".csv"))
  )
  utils::write.csv(
    standardized_means[[group]],
    file.path(output_dir, paste0("mean_network_standardized_", group, ".csv"))
  )
}

draw_network_figure(
  standardized_means,
  file.path(output_dir, "Figure_4_temporal_networks.pdf"),
  "pdf"
)
draw_network_figure(
  standardized_means,
  file.path(output_dir, "Figure_4_temporal_networks.tiff"),
  "tiff"
)

metric_draws <- calculate_metric_draws(standardized_arrays)
metric_summary <- summarize_metric_contrasts(metric_draws)
utils::write.csv(metric_summary, file.path(output_dir, "network_metric_contrasts.csv"), row.names = FALSE)
write_metric_tables(metric_summary, output_dir)
plot_metric_contrasts(metric_summary, file.path(output_dir, "Figure_5_network_metric_contrasts.pdf"))
plot_metric_contrasts(metric_summary, file.path(output_dir, "Figure_5_network_metric_contrasts.tiff"))

distance_draws <- network_distance_draws(standardized_arrays)
distance_summaries <- summarize_network_distances(distance_draws)
write_distance_table(distance_summaries, output_dir)

if (identical(tolower(Sys.getenv("SAVE_NETWORK_POSTERIOR", unset = "false")), "true")) {
  saveRDS(metric_draws, file.path(output_dir, "network_metric_draws.rds"), compress = FALSE)
  saveRDS(distance_draws, file.path(output_dir, "network_distance_draws.rds"), compress = FALSE)
}

figure_legends <- c(
  "Figure 4 | Standardized lagged networks by baseline smartphone-addiction group and sex.",
  "Panels show posterior-mean standardized lagged coefficients for non-addicted women (a), addicted women (b), non-addicted men (c) and addicted men (d). Rows of each adjacency matrix are source nodes at T1 and columns are target nodes at T2. Blue edges are positive and orange edges are negative. Edge width represents absolute coefficient magnitude. For legibility, coefficients with an absolute posterior mean below 0.10 are not displayed. This display threshold is descriptive and does not constitute an uncertainty-based edge selection rule. All networks use the same node positions and edge-width scale.",
  "",
  "Figure 5 | Posterior contrasts in standardized network summaries.",
  "Points show posterior means and horizontal lines show 95% credible intervals. Addiction-status contrasts are calculated as addicted minus non-addicted. The male-minus-female addiction contrast directly compares the two sex-specific addiction-status contrasts. Spectral radii summarize domain-specific transition matrices. Bridge strength sums absolute incoming and outgoing cross-community coefficients. Outgoing support associations retain coefficient signs. Impulse-propagation summaries apply a one-standard-deviation input to a support node and repeatedly apply the fitted transition matrix for 15 steps. The latter summaries describe the fitted model and do not estimate intervention effects."
)
writeLines(figure_legends, file.path(output_dir, "Figure_legends.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_postprocessing.txt"))
