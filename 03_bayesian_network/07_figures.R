network_plot_settings <- function(mean_networks) {
  set.seed(layout_seed)
  layout <- qgraph::qgraph(Reduce(`+`, mean_networks) / length(mean_networks),
    layout = "spring", repulsion = 1.1, DoNotPlot = TRUE)$layout
  maximum <- max(abs(unlist(mean_networks, use.names = FALSE)))
  list(layout = layout, maximum = maximum)
}

draw_network_figure <- function(mean_networks, file, device = c("pdf", "tiff")) {
  device <- match.arg(device)
  settings <- network_plot_settings(mean_networks)
  panels <- character(length(mean_networks))
  on.exit(unlink(panels), add = TRUE)
  colours <- c(rep("#E19C29", 3), rep("#5D8CA8", 7), rep("#6FAB5E", 8))
  labs <- c("Kin", "Off", "On", paste0("G", 1:7), paste0("MP", 1:8))
  titles <- c(nonaddicted_female = "Non-Addicted Female", addicted_female = "Addicted Female",
              nonaddicted_male = "Non-Addicted Male", addicted_male = "Addicted Male")
  for (i in seq_along(network_groups)) {
    group <- names(network_groups)[i]
    panels[i] <- tempfile(fileext = ".png")
    grDevices::png(panels[i], width = 2000, height = 2000, res = 300, bg = "white")
    graphics::par(xpd = NA, family = "sans")
    qgraph::qgraph(mean_networks[[group]], layout = settings$layout, color = colours,
      labels = labs, legend = FALSE, directed = TRUE, minimum = edge_display_threshold, maximum = settings$maximum,
      cut = .15, vsize = 8.5, label.cex = 1.2, label.font = 2, esize = 10, asize = 3.5,
      fade = TRUE, diag = TRUE, curveAll = TRUE, curve = .25,
      label.color = "white", borders = TRUE, border.color = "white", border.width = 2,
      posCol = "#003049", negCol = "#D62828", mar = c(4, 4, 6, 4))
    graphics::text(0, 1.4, titles[[group]], cex = 2.5, font = 2)
    grDevices::dev.off()
  }
  panel_plots <- lapply(panels, function(p) cowplot::ggdraw() + cowplot::draw_image(p))
  networks <- cowplot::plot_grid(plotlist = panel_plots, ncol = 2, labels = LETTERS[1:4],
    label_size = 35, label_x = .02, label_y = .95)
  legend_data <- data.frame(community = factor(c("Social Support", "Anxiety", "Phone Addiction"),
    levels = c("Social Support", "Anxiety", "Phone Addiction")), x = 1, y = 1)
  legend <- cowplot::get_legend(ggplot2::ggplot(legend_data,
      ggplot2::aes(x, y, colour = community)) + ggplot2::geom_point(size = 10) +
    ggplot2::scale_colour_manual(values = c("#E19C29", "#5D8CA8", "#6FAB5E"),
      name = "Node Communities") + ggplot2::theme_void() + ggplot2::theme(
      legend.position = "right", legend.text = ggplot2::element_text(size = 24),
      legend.title = ggplot2::element_text(size = 28, face = "bold",
        margin = ggplot2::margin(b = 20)), legend.key.height = grid::unit(1.5, "cm"),
      legend.spacing.y = grid::unit(1, "cm")))
  figure <- cowplot::plot_grid(networks, legend, nrow = 1, rel_widths = c(1, .25))
  if (device == "pdf") ggplot2::ggsave(file, figure, width = 22, height = 16,
    device = grDevices::cairo_pdf, bg = "white")
  else ggplot2::ggsave(file, figure, width = 22, height = 16, dpi = 300,
    device = ragg::agg_tiff, compression = "lzw", bg = "white")
  invisible(figure)
}
plot_metric_contrasts <- function(summary, output_file) {
  suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
d <- summary
cats <- c("Spectral radius","Bridge strength","Outgoing support sum","Impulse propagation")
d$cat <- factor(d$category,levels=cats)
d$sex <- factor(ifelse(grepl("female$",d$contrast),"Female","Male"),levels=c("Female","Male"))
d$evidence <- factor(ifelse(d$cri_lower>0|d$cri_upper<0,"95% CrI excludes 0",
                            ifelse(d$probability_of_direction>.90,"CrI includes 0\nPd > 0.90","Remaining estimates")),
                     levels=c("95% CrI excludes 0","CrI includes 0\nPd > 0.90","Remaining estimates"))
label <- function(category,m) {
  if (category==cats[1]) {
    x <- sub("^Spectral radius: ","",m)
    return(paste(c("support"="Support","anxiety"="Anxiety","smartphone addiction"="Addiction")[[x]],"Stability"))
  }
  if (category==cats[2]) return(paste0("Bridge_",sub("^Bridge strength: ","",m)))
  x <- sub("^(Outgoing sum|Cumulative response): ","",m)
  b <- strsplit(x," to ",fixed=TRUE)[[1]]
  src <- c(kin="Kin",offline="Offline",online="Online")[[b[1]]]
  tar <- if(b[2]=="anxiety") "Anxiety" else "Addiction"
  paste0(if(category==cats[3]) "Potential" else "Simulated"," ",src," (vs. ",tar,")")
}
d$lab <- mapply(label,d$category,d$metric,USE.NAMES=FALSE)
orders <- list(c("Addiction Stability","Anxiety Stability","Support Stability"),
               c(paste0("Bridge_gad",1:7),"Bridge_kin",paste0("Bridge_mpiq",1:8),"Bridge_offline","Bridge_online"),
               c("Potential Kin (vs. Addiction)","Potential Kin (vs. Anxiety)","Potential Offline (vs. Addiction)","Potential Offline (vs. Anxiety)","Potential Online (vs. Addiction)","Potential Online (vs. Anxiety)"),
               c("Simulated Kin (vs. Addiction)","Simulated Kin (vs. Anxiety)","Simulated Offline (vs. Addiction)","Simulated Offline (vs. Anxiety)","Simulated Online (vs. Addiction)","Simulated Online (vs. Anxiety)"))
cs <- scale_color_manual(name="Sex",values=c(Female="#D55E00",Male="#0072B2"))
ss <- scale_shape_manual(name="Posterior summary",values=c("95% CrI excludes 0"=16,"CrI includes 0\nPd > 0.90"=17,"Remaining estimates"=1))
th <- theme_bw(base_size=10) + theme(panel.grid.major=element_line(colour="#E7E7E7",linewidth=.28),
        panel.grid.minor=element_line(colour="#F3F3F3",linewidth=.2),
        plot.title=element_text(face="bold",hjust=.5,size=10,margin=margin(b=4)),
        axis.text.y=element_text(size=7.1,colour="#555555"),axis.text.x=element_text(size=7.8),
        axis.title.y=element_text(size=8.5),legend.position="none",plot.margin=margin(3,4,3,4))
titles <- c("Sub-network Stability","Bridge Centrality","Potential Influence (Static)","Cumulative Propagated Response")
panel <- function(i,sex_specific) {
  x <- d[d$category==cats[i] & d$contrast %in% if(sex_specific) c("Addicted minus non-addicted, female","Addicted minus non-addicted, male") else "Addicted minus non-addicted, averaged across sexes",]
  x$lab <- factor(x$lab,levels=rev(orders[[i]]))
  if (sex_specific) {
    p <- ggplot(x,aes(posterior_mean,lab,colour=sex))+
      geom_vline(xintercept=0,linetype="dashed",colour="#777777",linewidth=.35)+
      geom_errorbarh(aes(xmin=cri_lower,xmax=cri_upper,group=interaction(lab,sex)),height=0,linewidth=.55,position=position_dodge(width=.58))+
      geom_point(aes(shape=evidence,group=interaction(lab,sex)),size=2,position=position_dodge(width=.58))+cs
  } else {
    p <- ggplot(x,aes(posterior_mean,lab))+
      geom_vline(xintercept=0,linetype="dashed",colour="#777777",linewidth=.35)+
      geom_errorbarh(aes(xmin=cri_lower,xmax=cri_upper),height=0,linewidth=.55,colour="#1A1A1A")+
      geom_point(aes(shape=evidence),size=2,colour="#1A1A1A")
  }
  p+ss+labs(title=paste0(LETTERS[i+if(sex_specific)4 else 0],": ",titles[i]),x=NULL,y=NULL)+th
}
headplot <- function(main,sub,fill) ggplot()+
  annotate("rect",xmin=0,xmax=1,ymin=0,ymax=1,fill=fill)+
  annotate("text",x=.5,y=.66,label=main,fontface="bold",size=3.8)+
  annotate("text",x=.5,y=.25,label=sub,fontface="italic",size=2.8)+
  xlim(0,1)+ylim(0,1)+theme_void()+theme(plot.margin=margin(0,4,0,4))
left <- wrap_plots(c(list(headplot("Overall addiction-group contrast","(Addicted - Non-addicted)","#F1F1F1")),lapply(1:4,panel,sex_specific=FALSE)),ncol=1,heights=c(.06,.15,.49,.17,.17))
right <- wrap_plots(c(list(headplot("Sex-specific addiction-group contrast","(Addicted - Non-addicted within sex)","#EDF3F7")),lapply(1:4,panel,sex_specific=TRUE)),ncol=1,heights=c(.06,.15,.49,.17,.17))
dummy <- data.frame(sex=factor(c("Female","Male","Female"),levels=c("Female","Male")),
                    evidence=factor(c("95% CrI excludes 0","CrI includes 0\nPd > 0.90","Remaining estimates"),levels=levels(d$evidence)))
leg <- cowplot::get_legend(ggplot(dummy,aes(1,1,colour=sex,shape=evidence))+geom_point(size=2)+cs+ss+theme_void()+theme(legend.position="right",legend.box="vertical",legend.text=element_text(size=7.2),legend.title=element_text(size=8)))
fig <- cowplot::plot_grid(left,right,leg,nrow=1,rel_widths=c(1,1,.28),align="h")

  extension <- tolower(tools::file_ext(output_file))
  temporary <- tempfile(fileext=paste0(".",extension))
  on.exit(unlink(temporary))
  if (extension %in% c("tif","tiff")) {
    ragg::agg_tiff(temporary,width=12.8,height=10.2,units="in",res=600,compression="lzw")
  } else grDevices::cairo_pdf(temporary,width=12.8,height=10.2,family="Arial")
  print(fig)
  grDevices::dev.off()
  if (!file.copy(temporary,output_file,overwrite=TRUE)) stop("Could not write metric figure.")
  invisible(fig)
}
