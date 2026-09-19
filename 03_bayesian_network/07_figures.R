network_plot_settings <- function(mean_networks) {
  average <- Reduce(`+`, mean_networks) / length(mean_networks)
  set.seed(layout_seed)
  layout <- qgraph::qgraph(
    average,
    layout = "spring",
    DoNotPlot = TRUE,
    repulsion = 1.1
  )$layout
  maximum <- max(abs(unlist(mean_networks)))
  list(layout = layout, maximum = maximum)
}

draw_network_panel <- function(matrix, title, settings) {
  group_indices <- lapply(network_communities, function(x) match(x, rownames(matrix)))
  qgraph::qgraph(
    matrix,
    layout = settings$layout,
    groups = group_indices,
    color = c("#E69F00", "#56B4E9", "#009E73"),
    labels = network_node_labels,
    shape = "circle",
    vsize = 7.5,
    label.cex = 1.05,
    label.color = "white",
    label.font = 2,
    borders = TRUE,
    border.color = "white",
    border.width = 1.5,
    directed = TRUE,
    threshold = edge_display_threshold,
    maximum = settings$maximum,
    esize = 9,
    asize = 3,
    posCol = "#0072B2",
    negCol = "#D55E00",
    fade = TRUE,
    diag = TRUE,
    curveAll = TRUE,
    curve = 0.22,
    legend = FALSE,
    mar = c(3, 3, 5, 3)
  )
  title(main = title, line = 2, cex.main = 1.35)
}

draw_network_figure <- function(mean_networks, file, device = c("pdf", "tiff")) {
  device <- match.arg(device)
  settings <- network_plot_settings(mean_networks)
  staged_file <- tempfile(
    pattern = "network-figure-",
    fileext = if (device == "pdf") ".pdf" else ".tiff"
  )
  if (device == "pdf") {
    grDevices::pdf(staged_file, width = 14, height = 11, family = "Times", useDingbats = FALSE)
  } else {
    grDevices::tiff(staged_file, width = 14, height = 11, units = "in", res = 600, compression = "lzw")
  }
  device_open <- TRUE
  on.exit({
    if (device_open) grDevices::dev.off()
    unlink(staged_file)
  }, add = TRUE)

  graphics::par(mfrow = c(2, 2), oma = c(2, 2, 3, 9), family = "serif")
  panel_order <- c(
    "nonaddicted_female", "addicted_female",
    "nonaddicted_male", "addicted_male"
  )
  panel_labels <- c("a", "b", "c", "d")
  for (i in seq_along(panel_order)) {
    group <- panel_order[i]
    draw_network_panel(mean_networks[[group]], network_groups[[group]], settings)
    graphics::mtext(panel_labels[i], side = 3, line = 1, adj = 0, font = 2, cex = 1.2)
  }
  graphics::par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend(
    "right",
    inset = c(0.015, 0),
    legend = c("Social support", "Anxiety", "Smartphone addiction"),
    pch = 21,
    pt.bg = c("#E69F00", "#56B4E9", "#009E73"),
    pt.cex = 2.1,
    bty = "n",
    title = "Node community",
    cex = 0.95
  )
  grDevices::dev.off()
  device_open <- FALSE
  if (!file.copy(staged_file, file, overwrite = TRUE)) {
    stop("Could not save figure: ", file)
  }
}

plot_metric_contrasts <- function(summary, output_file) {
  display_contrasts <- c(
    "Addicted minus non-addicted, averaged across sexes",
    "Addicted minus non-addicted, female",
    "Addicted minus non-addicted, male",
    "Male minus female addiction contrast"
  )
  plot_data <- summary |>
    dplyr::filter(.data$contrast %in% display_contrasts) |>
    dplyr::mutate(
      contrast = factor(.data$contrast, levels = display_contrasts),
      category = factor(
        .data$category,
        levels = c(
          "Spectral radius", "Bridge strength",
          "Outgoing support sum", "Impulse propagation"
        )
      ),
      metric_label = metric_display_label(.data$metric)
    )

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$posterior_mean, y = .data$metric_label)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$cri_lower, xmax = .data$cri_upper),
      orientation = "y",
      width = 0,
      linewidth = 0.45
    ) +
    ggplot2::geom_point(size = 1.7) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$category),
      cols = ggplot2::vars(.data$contrast),
      scales = "free",
      space = "free_y"
    ) +
    ggplot2::labs(x = "Posterior contrast", y = NULL) +
    ggplot2::theme_bw(base_size = 9, base_family = "serif") +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 8),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.7, "lines")
    )

  extension <- tolower(tools::file_ext(output_file))
  staged_file <- tempfile(pattern = "network-metrics-", fileext = paste0(".", extension))
  on.exit(unlink(staged_file), add = TRUE)
  if (extension %in% c("tif", "tiff")) {
    ggplot2::ggsave(
      staged_file,
      plot,
      width = 15,
      height = 13,
      device = "tiff",
      dpi = 600,
      compression = "lzw"
    )
  } else {
    ggplot2::ggsave(
      staged_file,
      plot,
      width = 15,
      height = 13,
      device = grDevices::cairo_pdf
    )
  }
  if (!file.copy(staged_file, output_file, overwrite = TRUE)) {
    stop("Could not save figure: ", output_file)
  }
  plot
}
