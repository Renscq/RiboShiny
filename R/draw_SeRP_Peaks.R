# function to draw enrichment peaks plot


library(ggplot2)
library(RColorBrewer)
library(tidyverse)


draw_serp_peaks <- function(
    gene_name = NULL,
    gene_peaks = NULL,
    
    x = NULL,
    enrich = NULL,
    edge = NULL,
    bound = NULL,

    xlabel = "Distance from TIS [AA.]",
    ylabel = "Enrichment [a.u.]",
    title = NULL,
    
    line_color = c("#0072bd", "#edb120", "#d95319"),
    
    line_width = 1,
    fill_alpha = 0.8,
    
    xmin = NULL,
    xmax = NULL,
    xbreaks = 5,
    
    ymin = NULL,
    ymax = NULL,
    ybreaks = 2,
    
    grid_width = 0.2,
    grid_color = "grey85",
    sqrt = FALSE,
    
    show_legend = FALSE,
    
    font_size = 12) {

  # get the max value of the gene peaks
  max_value <- max(gene_peaks[enrich])
  
  # browser()

  # set the colors
  line_color <- alpha(line_color, fill_alpha)
  names(line_color) <- c(enrich, edge, bound)
  
  # draw the line plot
  plots <- ggplot(data = gene_peaks) +
    geom_line(aes(x = !!sym(x), y = !!sym(enrich), color = "enrich"),
              linewidth = line_width) + 
    geom_line(aes(x = !!sym(x), y = !!sym(edge), color = "edge"),
              linewidth = line_width) + 
    geom_line(aes(x = !!sym(x), y = !!sym(bound), color = "bound"),
              linewidth = line_width) + 
    scale_color_manual(values = line_color) +
    theme_bw() +
    theme(legend.text = element_text(size = font_size - 3),
          legend.title = element_text(size = font_size - 2),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = grid_width, color = grid_color),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, size = font_size + 2),
          axis.title = element_text(color = "black", hjust = 0.5, vjust = 0.5, size = font_size + 1),
          axis.text = element_text(color = "black", hjust = 0.5, vjust = 0.5, size = font_size)) +
    # geom_vline(aes(xintercept = 2), size = 1, colour = "#BB0000", linetype = "dashed") + 
    # xlim(xmin, xmax) +
    # scale_y_sqrt() +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(gene_name)
  
  # set the limits
  if (!is.na(xmin) && !is.na(xmax) && !is.na(xbreaks)) {
    plots <- plots + 
      scale_x_continuous(limits = c(xmin, xmax), n.breaks = xbreaks)
  }
  
  # set the limits
  if (!is.na(ymin) && !is.na(ymax) && !is.na(ybreaks)) {
    plots <- plots + 
      scale_y_continuous(limits = c(ymin, ymax), n.breaks = ybreaks)
  }

  # sqrt the plot
  if (isTRUE(sqrt)) {
    plots <- plots +  scale_y_sqrt()
  }
  
  # ggplot2::ggsave(plot = plots, filename = './images/reads_heat_length.png', width = 9, height = 6, dpi = 300)
  
  return(plots)
}