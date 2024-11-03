# function to draw the gene reads mapping


draw_isoforms_mapped <- function(
    gene_name = NULL,
    gene_reads = NULL,
    
    x = NULL,
    y = NULL,
    color = NULL,
    
    xlabel = "Nucleotide position",
    ylabel = "Abundance (RPM)",
    title = NULL,
    
    gene_color = 'Set1',
    plot_type = "line", # line or bar, area
    
    line_width = 1,
    grid_width = 0.2,
    fill_alpha = 1,
    
    xmin = NULL,
    xmax = NULL,
    xbreaks = 5,
    
    ymin = NULL,
    ymax = NULL,
    ybreaks = 2,
    
    sqrt = 'none',
    
    facet = TRUE,
    show_legend = FALSE,

    font_size = 12) {
  
  # generate the color palette
  gene_color_num <- length(unique(gene_reads[[color]]))
  
  if (gene_color_num >= 8) {
  	color_list <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = gene_color))(gene_color_num)
  } else if (gene_color_num < 8 && gene_color_num >= 3) {
  	color_list <- RColorBrewer::brewer.pal(n = gene_color_num, name = gene_color)
  } else if (gene_color_num < 3) {
    color_list <- RColorBrewer::brewer.pal(n = 3, name = gene_color)
    color_list <- color_list[1:gene_color_num]
  }
  
  # create the gene reads plot
  gene_reads_plot <- ggplot(data = gene_reads, aes(x = !!sym(x), y = !!sym(y)))
    
  # draw the gene reads
  if (plot_type == "line") {
    gene_reads_plot <- gene_reads_plot +
      geom_line(aes(color = !!sym(color)), alpha = fill_alpha, linewidth = line_width) +
      scale_color_manual(values = color_list)
    
  } else if (plot_type == "area") {
    gene_reads_plot <- gene_reads_plot +
      geom_area(aes(fill = !!sym(color)), alpha = fill_alpha) +
      scale_fill_manual(values = color_list)
    
  } else if (plot_type == "bar") {
    gene_reads_plot <- gene_reads_plot +
      geom_col(aes(fill = !!sym(color), color = !!sym(color)), alpha = fill_alpha) +
      scale_color_manual(values = color_list)
  }

  # set the theme
  gene_reads_plot <- gene_reads_plot +
    theme_bw() +
    theme(
      plot.title = element_text(size = font_size + 2, vjust = 0.5, hjust = 0.5),
      axis.title = element_text(size = font_size + 1),
      axis.text = element_text(size = font_size),
      strip.text = element_text(size = font_size - 1),
      panel.grid.major = element_line(linewidth = grid_width),
      panel.grid.minor = element_blank(),
      strip.background = element_blank()) +
    labs(title = title, x = xlabel, y = ylabel) + 
    theme(legend.position = show_legend,
          legend.text = element_text(size = font_size - 2))
  
  # reset the x-axis limits
  if (!is.na(xmin) & !is.na(xmax)) {
    gene_reads_plot <- gene_reads_plot +
      coord_cartesian(xlim = c(xmin, xmax), clip = "on") + 
      scale_x_continuous(limits = c(xmin, xmax), n.breaks = xbreaks)
  }
  
  # reset the y-axis limits
  if (!is.na(ymin) & !is.na(ymax)) {
    gene_reads_plot <- gene_reads_plot +
      scale_y_continuous(limits = c(ymin, ymax), n.breaks = ybreaks)
  } else {
    ymin <- min(gene_reads[[y]], na.rm = TRUE)
    ymax <- max(gene_reads[[y]], na.rm = TRUE)
    
    gene_reads_plot <- gene_reads_plot +
      scale_y_continuous(limits = c(ymin, ymax), n.breaks = ybreaks)
  }
  
  # set the sqrt
  if (sqrt == 'sqrt') {
    gene_reads_plot <- gene_reads_plot +
      scale_y_sqrt()
  }
  
  # facet the plot
  if (facet) {
    facet_formula <- as.formula(paste0(color, "~."))
    
    gene_reads_plot <- gene_reads_plot +
      facet_grid(facet_formula, scales = "free_x")
  }
  
  return(gene_reads_plot)
}

