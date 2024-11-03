# function to draw coverage plot

library(ggplot2)
library(RColorBrewer)
library(MetBrewer)
library(tidyverse)


draw_coverage <- function(
    coverage = NULL,
    x = NULL,
    y = NULL,
    xstart = NULL,
    xend = NULL,
    font_size = 12,
    line_color = "Blues", # only for the gene plot
    line_width = 0.8, # only for the line plot
    group = NULL, # group for fill
    fill_color = "Blues", # Blues, Reds, Greens ...
    fill_alpha = 1,
    facet = FALSE,
    wrap_group = NULL, # group for facet_wrap
    scales = "free_y",
    ncol = NULL,
    nrow = NULL,
    plot_type = "Line" # line plot, bar plot or heat map
) {
  
  # browser()
  
  # draw the coverage plot
  if (plot_type == "Line") {
    
    # retrieve the colors
    color_num <- length(unique(coverage[[group]]))
    line_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = line_color))(color_num + 2)
    # line_colors <- line_color_set[2:c(color_num + 1)]
    line_colors <- line_color_set[1:color_num]
    line_colors <- alpha(line_colors, fill_alpha)
    
    coverage_plot <- ggplot(data = coverage) +
      geom_line(aes(x = !!sym(x), 
                    y = !!sym(y), 
                    color = !!sym(group)),
                linewidth = line_width) +
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.text = element_text(color = 'black', size = font_size),
            axis.text = element_text(color = 'black', size = font_size),
            axis.title = element_text(color = 'black', size = font_size + 1, vjust = 0.5, hjust = 0.5),
            plot.title = element_text(color = 'black', size = font_size + 2, vjust = 0.5, hjust = 0.5),
            legend.text = element_text(color = 'black', size = font_size - 3),
            legend.title = element_text(color = 'black', size = font_size - 2),
            panel.grid = element_blank(),
            legend.position = "right") +
      scale_color_manual(values = line_colors) +
      # xlim(xstart, xend) +
      xlab(x) +
      ylab(y)
    
    # warp the coverage plot
    if (facet) {
      coverage_plot <- coverage_plot +
        facet_wrap(`group`, scales = "free_x")
    }
    
    # ggplot2::ggsave(plot = coverage_plot, filename = './output/coverage_line_plot.png', width = 12, height = 12, dpi = 300)
  } else if(plot_type == "Heat") {
    # retrieve the colors
    fill_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = fill_color))(9)
    
    # browser()
    
    coverage_plot <- ggplot(data = coverage) +
      geom_tile(aes(x = !!sym(x),
                    y = !!sym(y),
                    fill = !!sym(group)),
                alpha = fill_alpha) +
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.text = element_text(color = 'black', size = font_size),
            axis.text.x = element_text(color = 'black', size = font_size),
            axis.text.y = element_text(color = 'black', size = font_size, angle = 0, hjust = 0.5, vjust = 0.5),
            axis.title = element_text(color = 'black', size = font_size + 1, vjust = 0.5, hjust = 0.5),
            plot.title = element_text(color = 'black', size = font_size + 2, vjust = 0.5, hjust = 0.5),
            legend.text = element_text(color = 'black', size = font_size - 3),
            legend.title = element_text(color = 'black', size = font_size - 2),
            panel.grid = element_blank()) +
      scale_fill_gradient(low = fill_color_set[2], high = fill_color_set[8]) +
      # facet_wrap(~Region, scales = "free_x") +
      # xlim(xstart, xend) +
      xlab(x) +
      ylab(y)
    
    # warp the coverage plot
    if (facet) {
      coverage_plot <- coverage_plot +
        facet_wrap(`group`, scales = "free_x")
    }
  }

  # ggplot2::ggsave(plot = coverage_plot, filename = './output/coverage_heat_plot.png', width = 9, height = 9, dpi = 300)
  
  return(coverage_plot)
  
}