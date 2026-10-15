if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.UTF-8"))
}
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
source(file.path(project_root, "10_derived_diagnostics.R"), encoding = "UTF-8")
source(file.path(Sys.getenv("PROJECT_ROOT", normalizePath("..", winslash = "/")), "R", "project.R"), encoding = "UTF-8")

fit <- load_network_model(model_file)
validate_network_model(fit, expected_data = read_network_data())

diagnostics <- model_diagnostics(fit)
write_diagnostic_table(diagnostics, output_dir)

edge_draws <- extract_edge_draws(fit)

raw_arrays <- reconstruct_network_arrays(edge_draws)
standardized_arrays <- standardize_network_arrays(raw_arrays, fit$data)
raw_means <- posterior_mean_networks(raw_arrays)
standardized_means <- posterior_mean_networks(standardized_arrays)
display_scale <- Sys.getenv("NETWORK_DISPLAY_SCALE", "raw")
if (!display_scale %in% c("raw", "standardized")) stop("NETWORK_DISPLAY_SCALE must be raw or standardized.")
display_means <- if (display_scale == "raw") raw_means else standardized_means

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
  display_means,
  file.path(output_dir, "Figure_4_temporal_networks.pdf"),
  "pdf"
)
draw_network_figure(
  display_means,
  file.path(output_dir, "Figure_4_temporal_networks.tiff"),
  "tiff"
)

metric_draws <- calculate_metric_draws(standardized_arrays)
metric_summary <- summarize_metric_contrasts(metric_draws)
write_derived_diagnostics(fit, standardized_arrays, metric_draws, output_dir)
group_summaries <- dplyr::bind_rows(lapply(dimnames(metric_draws)$group, function(g) {
  dplyr::bind_rows(lapply(dimnames(metric_draws)$metric, function(m)
    data.frame(group = g, metric = m, as.list(posterior_summary(metric_draws[, m, g])))))
}))
write_csv(group_summaries, file.path(output_dir, "group_network_metric_summaries.csv"))
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
  paste0("Figure 4 | ", if (display_scale == "raw") "Unstandardized" else "Standardized", " lagged networks by baseline smartphone-addiction group and sex."),
  paste0("Panels show non-addicted women (A), addicted women (B), non-addicted men (C) and addicted men (D). Edges are posterior-mean lagged coefficients; baseline sources are in rows and follow-up targets in columns. Blue edges are positive and red edges are negative. Widths share a magnitude scale, with no upper display cap. Absolute means below ", edge_display_threshold, " are omitted from the plots only. This display threshold is descriptive. All networks use the same node positions. Source data are the exported mean matrices on the indicated coefficient scale."),
  "",
  "Figure 5 | Posterior contrasts in standardized network summaries.",
  paste0("Points show posterior means and lines show 95% credible intervals. The left column averages addicted-minus-non-addicted contrasts equally across sexes; the right column shows within-sex contrasts. Orange denotes women and blue denotes men. Filled circles identify intervals excluding zero, triangles identify remaining Pd > 0.90 estimates and open circles identify the rest. Pd is the larger sign probability; symbols are descriptive and not multiplicity-adjusted tests. Sub-network Stability is the within-domain spectral radius, Bridge Centrality sums absolute cross-domain weights, Potential Influence retains outgoing coefficient signs, and Cumulative Propagated Response sums signed responses over ", impulse_steps, " steps after a +1 SD input. All metrics use standardized matrices. These observational propagation summaries do not estimate causal intervention benefits. Source data: network_metric_contrasts.csv.")
)
writeLines(figure_legends, file.path(output_dir, "Figure_legends.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_postprocessing.txt"))
