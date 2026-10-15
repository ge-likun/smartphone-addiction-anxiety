parameter_labels <- c(
  mu = "Relative mean of the positive-score component (\u03bc)",
  nu = "Relative odds of GAD-7 = 0 (\u03bd)",
  sigma = "Relative positive-score dispersion (\u03c3)"
)
grid_length <- 180L
loess_span <- 1.0

support_score <- function(model, variable) {
  column <- paste0("pb(", variable, ")")
  if (is.null(model$mu.x) || !column %in% colnames(model$mu.x)) stop("Support score not found in model matrix: ", column)
  as.numeric(model$mu.x[, column])
}

full_model_curve <- function(model, variable, parameter) {
  score <- support_score(model, variable)
  fitted_link <- as.numeric(predict(model, what = parameter, type = "link"))
  grouped <- tibble(support_score = score, fitted_link = fitted_link) |>
    filter(is.finite(.data$support_score), is.finite(.data$fitted_link)) |>
    group_by(.data$support_score) |>
    summarise(fitted_link = mean(.data$fitted_link), frequency = n(), .groups = "drop") |>
    arrange(.data$support_score)
  smoother <- loess(
    fitted_link ~ support_score, data = grouped, weights = grouped$frequency,
    span = loess_span, degree = 2, family = "gaussian",
    control = loess.control(surface = "direct")
  )
  grid <- seq(min(grouped$support_score), max(grouped$support_score), length.out = grid_length)
  prediction <- predict(smoother, newdata = data.frame(support_score = grid), se = TRUE)
  fit <- as.numeric(prediction$fit)
  se <- as.numeric(prediction$se.fit)
  if (any(!is.finite(fit))) fit <- approx(grid[is.finite(fit)], fit[is.finite(fit)], grid, rule = 2)$y
  if (any(!is.finite(se))) se <- approx(grid[is.finite(se)], se[is.finite(se)], grid, rule = 2)$y
  reference <- as.numeric(predict(smoother, newdata = data.frame(support_score = 0)))
  if (!is.finite(reference)) stop("The score-0 reference cannot be estimated for ", variable,
                                  "; preserve the specified reference and review the observed score range.")
  tibble(
    support_score = grid,
    ratio = exp(fit - reference),
    ratio_low = exp(fit - 1.96 * se - reference),
    ratio_high = exp(fit + 1.96 * se - reference)
  )
}

curve_rows <- list()
histogram_rows <- list()
for (spec_index in seq_len(nrow(analysis_specs))) {
  spec <- analysis_specs[spec_index, ]
  model <- primary_models[[spec$model_id]]
  for (support_key in names(support_variables)) {
    variable <- support_variables[[support_key]]
    metadata <- list(
      model_id = spec$model_id, model_set = spec$model_set, stratum = spec$stratum,
      addiction_group = spec$addiction_group, sex = spec$sex,
      support_key = support_key, support_label = support_labels[[support_key]]
    )
    histogram_rows[[paste(spec$model_id, support_key)]] <- tibble(support_score = support_score(model, variable)) |>
      mutate(!!!metadata, .before = 1)
    for (parameter in names(parameter_labels)) {
      curve_rows[[paste(spec$model_id, support_key, parameter)]] <- full_model_curve(model, variable, parameter) |>
        mutate(!!!c(metadata, list(parameter = parameter)), .before = 1)
    }
  }
}
curve_data <- bind_rows(curve_rows)
histogram_observations <- bind_rows(histogram_rows)
histogram_density <- histogram_observations |>
  mutate(support_score = as.integer(round(.data$support_score))) |>
  count(
    .data$model_id, .data$model_set, .data$stratum, .data$addiction_group,
    .data$sex, .data$support_key, .data$support_label, .data$support_score,
    name = "frequency"
  ) |>
  group_by(.data$model_id, .data$support_key) |>
  mutate(density = .data$frequency / sum(.data$frequency)) |>
  ungroup()
data_limits <- function(model_set, parameter) {
  rows <- curve_data$model_set == model_set & curve_data$parameter == parameter
  upper <- max(c(1, curve_data$ratio_high[rows]), na.rm = TRUE)
  c(0, upper * 1.08)
}
overall_limits <- setNames(lapply(names(parameter_labels), function(p) data_limits("overall", p)), names(parameter_labels))
sex_limits <- setNames(lapply(names(parameter_labels), function(p) data_limits("sex_stratified", p)), names(parameter_labels))
write.csv(curve_data, file.path(output_dirs$source_data, "GAMLSS_association_curve_data.csv"), row.names = FALSE)
write.csv(histogram_density, file.path(output_dirs$source_data, "GAMLSS_support_histogram_data.csv"), row.names = FALSE)

theme_manuscript <- function() {
  theme_bw(base_size = 7.1, base_family = "Arial") +
    theme(
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.20),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey45", fill = NA, linewidth = 0.38),
      plot.title = element_text(hjust = 0.5, face = "plain", size = 7.6),
      axis.title = element_text(size = 6.9), axis.text = element_text(size = 6.3),
      legend.text = element_text(size = 6.7), legend.key = element_blank(),
      legend.key.width = grid::unit(13, "pt"), plot.margin = margin(5, 6, 5, 7)
    )
}

density_reference_max <- max(.50, max(histogram_density$density, na.rm = TRUE) * 1.3)
overall_density_breaks <- if (density_reference_max == .50) c(0, .20, .40) else pretty(c(0, density_reference_max), n = 3)
sex_density_breaks <- if (density_reference_max == .50) c(0, .15, .30) else overall_density_breaks
overall_histogram_display_multiplier <- 1.00
sex_histogram_display_multiplier <- 1.12

overall_colours <- c(`Non-addicted` = "#009E73", Addicted = "#CC79A7")
overall_linetypes <- c(`Non-addicted` = "solid", Addicted = "22")
overall_histogram_colours <- c(`Non-addicted` = "#E3F3EE", Addicted = "#F6E8F0")
sex_colours <- c(Female = "#D55E00", Male = "#0072B2")
sex_ribbon_colours <- c(Female = "#D55E00", Male = "#56B4E9")
sex_histogram_colours <- c(Female = "#F8EADF", Male = "#E7F2F8")

histogram_rectangles <- function(data, grouping, y_limits, display_multiplier) {
  offsets <- if (identical(grouping, "sex")) c(Female = -0.23, Male = 0.23) else c(`Non-addicted` = -0.23, Addicted = 0.23)
  density_scale <- y_limits[2] / density_reference_max * display_multiplier
  data |>
    mutate(
      group_value = as.character(.data[[grouping]]),
      xmin = .data$support_score + offsets[.data$group_value] - 0.22,
      xmax = .data$support_score + offsets[.data$group_value] + 0.22,
      ymin = 0, ymax = .data$density * density_scale
    )
}

check_curve_limits <- function(curves, y_limits, panel_name) {
  bounds <- c(curves$ratio_low, curves$ratio_high)
  if (any(!is.finite(bounds)) || any(bounds < y_limits[1] | bounds > y_limits[2])) {
    stop("Curve interval exceeds the plotting range in ", panel_name)
  }
}

figure2_row_tag <- function(label) {
  ggplot() +
    annotate(
      "text", x = 0.5, y = 0.98, label = label,
      hjust = 0.5, vjust = 1, fontface = "bold", size = 3.0
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(5, 0, 0, 0))
}

overall_panel <- function(parameter, support_key) {
  y_limits <- overall_limits[[parameter]]
  curves <- curve_data |>
    filter(.data$model_set == "overall", .data$parameter == .env$parameter, .data$support_key == .env$support_key)
  check_curve_limits(curves, y_limits, paste("overall", parameter, support_key))
  hist <- histogram_density |>
    filter(.data$model_set == "overall", .data$support_key == .env$support_key) |>
    histogram_rectangles("addiction_group", y_limits, overall_histogram_display_multiplier)
  right <- support_key == "Online_friend"
  secondary_axis <- if (right) sec_axis(~ . / (y_limits[2] / density_reference_max * overall_histogram_display_multiplier), name = "Density", breaks = overall_density_breaks) else waiver()
  ggplot() +
    geom_rect(data = hist, aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0, ymax = .data$ymax, fill = .data$group_value), alpha = 0.78, colour = NA) +
    geom_ribbon(data = filter(curves, .data$addiction_group == "Non-addicted"), aes(.data$support_score, ymin = .data$ratio_low, ymax = .data$ratio_high), fill = overall_colours[["Non-addicted"]], alpha = 0.09, colour = NA) +
    geom_ribbon(data = filter(curves, .data$addiction_group == "Addicted"), aes(.data$support_score, ymin = .data$ratio_low, ymax = .data$ratio_high), fill = overall_colours[["Addicted"]], alpha = 0.09, colour = NA) +
    geom_line(data = curves, aes(.data$support_score, .data$ratio, colour = .data$addiction_group, linetype = .data$addiction_group), linewidth = 0.72, lineend = "round") +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.34, colour = "grey30") +
    scale_colour_manual(values = overall_colours, breaks = c("Addicted", "Non-addicted"), name = NULL) +
    scale_fill_manual(values = overall_histogram_colours, guide = "none") +
    scale_linetype_manual(values = overall_linetypes, breaks = c("Addicted", "Non-addicted"), name = NULL) +
    guides(
      fill = "none", linetype = "none",
      colour = guide_legend(override.aes = list(linetype = c("22", "solid"), linewidth = 0.85))
    ) +
    scale_x_continuous(breaks = c(0, 5, 10, 15), expand = expansion(mult = 0)) +
    scale_y_continuous(name = if (support_key == "Kinship") parameter_labels[[parameter]] else NULL, limits = y_limits, sec.axis = secondary_axis, expand = expansion(mult = 0)) +
    coord_cartesian(xlim = c(-0.5, 15.5), clip = "off", expand = FALSE) +
    labs(title = if (parameter == "mu") support_labels[[support_key]] else NULL, x = if (parameter == "sigma") "Support score" else NULL) +
    theme_manuscript()
}

overall_panels <- unlist(
  lapply(names(parameter_labels), function(parameter) {
    lapply(names(support_variables), function(key) overall_panel(parameter, key))
  }),
  recursive = FALSE
)
figure2_design <- "
ABCD
EFGH
IJKL
"
figure2 <- wrap_plots(
  A = figure2_row_tag("a"), B = overall_panels[[1]], C = overall_panels[[2]], D = overall_panels[[3]],
  E = figure2_row_tag("b"), F = overall_panels[[4]], G = overall_panels[[5]], H = overall_panels[[6]],
  I = figure2_row_tag("c"), J = overall_panels[[7]], K = overall_panels[[8]], L = overall_panels[[9]],
  design = figure2_design
) +
  plot_layout(widths = c(0.035, 1, 1, 1), heights = c(1, 1, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")
save_publication_figure(figure2, file.path(output_dirs$figures, "Figure_2_full_model_score0_common_axes"), 183, 176)

row_style <- list(
  background = c(
    `Non-addicted` = grDevices::adjustcolor("#E8EEF5", 0.80),
    Addicted = grDevices::adjustcolor("#F4E8E4", 0.80)
  ),
  text = c(`Non-addicted` = "#30485A", Addicted = "#8A452C")
)

sex_panel <- function(parameter, addiction_group, support_key) {
  y_limits <- sex_limits[[parameter]]
  curves <- curve_data |>
    filter(
      .data$model_set == "sex_stratified", .data$parameter == .env$parameter,
      .data$addiction_group == .env$addiction_group, .data$support_key == .env$support_key
    )
  check_curve_limits(curves, y_limits, paste("sex-stratified", addiction_group, parameter, support_key))
  hist <- histogram_density |>
    filter(
      .data$model_set == "sex_stratified", .data$addiction_group == .env$addiction_group,
      .data$support_key == .env$support_key
    ) |>
    histogram_rectangles("sex", y_limits, sex_histogram_display_multiplier)
  top <- addiction_group == "Non-addicted"
  right <- support_key == "Online_friend"
  secondary_axis <- if (right) sec_axis(~ . / (y_limits[2] / density_reference_max * sex_histogram_display_multiplier), name = "Density", breaks = sex_density_breaks) else waiver()
  ggplot() +
    geom_rect(data = hist, aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0, ymax = .data$ymax, fill = .data$sex), alpha = 0.78, colour = NA) +
    geom_ribbon(data = filter(curves, .data$sex == "Female"), aes(.data$support_score, ymin = .data$ratio_low, ymax = .data$ratio_high), fill = sex_ribbon_colours[["Female"]], alpha = 0.09, colour = NA) +
    geom_ribbon(data = filter(curves, .data$sex == "Male"), aes(.data$support_score, ymin = .data$ratio_low, ymax = .data$ratio_high), fill = sex_ribbon_colours[["Male"]], alpha = 0.09, colour = NA) +
    geom_line(data = curves, aes(.data$support_score, .data$ratio, colour = .data$sex, group = .data$sex), linewidth = 0.72, lineend = "round") +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.34, colour = "grey30") +
    scale_colour_manual(values = sex_colours, breaks = c("Female", "Male"), name = NULL) +
    scale_fill_manual(values = sex_histogram_colours, guide = "none") +
    scale_x_continuous(breaks = c(0, 5, 10, 15), expand = expansion(mult = 0)) +
    scale_y_continuous(name = NULL, limits = y_limits, sec.axis = secondary_axis, expand = expansion(mult = 0)) +
    coord_cartesian(xlim = c(-0.5, 15.5), clip = "off", expand = FALSE) +
    labs(title = if (top) support_labels[[support_key]] else NULL, x = if (!top) "Support score" else NULL) +
    theme_manuscript()
}

row_label <- function(label) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = row_style$background[[label]], colour = NA) +
    annotate("text", x = 0.5, y = 0.5, label = label, angle = 90, colour = row_style$text[[label]], size = 2.8) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + theme_void()
}

axis_label <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, angle = 90, size = 2.8) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
}

divider <- ggplot() + annotate("segment", x = 0, xend = 1, y = 0.5, yend = 0.5, colour = "#D8D8D8", linewidth = 0.35) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + theme_void()

sex_figure <- function(parameter) {
  design <- "\nJKABCD\nJKEEEE\nJKFGHI\n"
  wrap_plots(
    J = axis_label(parameter_labels[[parameter]]),
    K = plot_spacer(),
    A = row_label("Non-addicted"),
    B = sex_panel(parameter, "Non-addicted", "Kinship"),
    C = sex_panel(parameter, "Non-addicted", "Offline_friend"),
    D = sex_panel(parameter, "Non-addicted", "Online_friend"),
    E = divider,
    F = row_label("Addicted"),
    G = sex_panel(parameter, "Addicted", "Kinship"),
    H = sex_panel(parameter, "Addicted", "Offline_friend"),
    I = sex_panel(parameter, "Addicted", "Online_friend"),
    design = design
  ) +
    plot_layout(widths = c(0.05, 0.03, 0.07, 1, 1, 1), heights = c(1, 0.018, 1), guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal")
}

save_publication_figure(sex_figure("mu"), file.path(output_dirs$figures, "Figure_3_mu_full_model_score0_common_axes"), 183, 116)
save_publication_figure(sex_figure("nu"), file.path(output_dirs$figures, "Supplementary_Figure_3_nu_full_model_score0_common_axes"), 183, 116)
save_publication_figure(sex_figure("sigma"), file.path(output_dirs$figures, "Supplementary_Figure_4_sigma_full_model_score0_common_axes"), 183, 116)

zaga_worm <- representative_worm |>
  filter(.data$family == "ZAGA") |>
  left_join(analysis_specs |> select(model_id, addiction_group, sex, model_set), by = "model_id") |>
  mutate(
    row = factor(ifelse(.data$model_set == "overall", "Overall", ifelse(.data$sex == "Male", "Men", "Women")), levels = c("Overall", "Men", "Women")),
    addiction_group = factor(.data$addiction_group, levels = c("Non-addicted", "Addicted"))
  )
write.csv(zaga_worm, file.path(output_dirs$source_data, "Supplementary_Figure_2_ZAGA_residuals.csv"), row.names = FALSE)

supplementary_figure_4 <- ggplot(zaga_worm, aes(.data$theoretical, .data$deviation)) +
  geom_ribbon(aes(ymin = .data$envelope_low, ymax = .data$envelope_high), fill = "grey85", colour = NA) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.35) +
  geom_line(aes(y = .data$fitted_curve), colour = "#0072B2", linewidth = 0.65) +
  facet_grid(rows = vars(row), cols = vars(addiction_group)) +
  labs(x = "Unit-normal quantile", y = "Residual deviation") +
  theme_bw(base_size = 8, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"))
save_publication_figure(supplementary_figure_4, file.path(output_dirs$figures, "Supplementary_Figure_2_ZAGA_residuals"), 150, 170)
