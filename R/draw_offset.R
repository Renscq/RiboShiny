# function to draw reads offset

library(ggplot2)
library(RColorBrewer)
library(tidyverse)


draw_offset <- function(
    offset = NULL,
    x = NULL,
    y = NULL,
    xstart = NULL,
    xend = NULL,
    wrap_group = NULL, # group for facet_wrap
    facet = FALSE,
    scales = "free_y",
    ncol = NULL,
    nrow = NULL,
    fill_group = NULL, # group for fill
    fill_color = "Blues", # Blues, Reds, Greens ...
    fill_alpha = 1,
    font_size = 12,
    bar_width = 0.8, # only for the bar plot
    plot_type = "bar" # bar plot or heat map
) {
  
  # browser()
  
  # retrieve the colors
  fill_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = fill_color))(9)
  
  # draw the offset
  if (plot_type == "bar") {

    offset_plot <- ggplot(data = offset) +
      geom_col(aes(x = !!sym(x), 
                   y = !!sym(y), 
                   fill = !!sym(fill_group),
                   group = !!sym(fill_group)),
               alpha = fill_alpha,
               position = position_dodge(width = bar_width),
               width = bar_width) +
      theme_bw() +
      theme(strip.background = element_blank(),
            panel.grid = element_blank(),
            legend.text = element_text(size = font_size - 3),
            legend.title = element_text(size = font_size - 2),
            axis.text = element_text(color = "black", size = font_size),
            axis.title = element_text(color = "black", size = font_size + 1),
            plot.title = element_text(size = font_size + 2, hjust = 0.5, vjust = 0.5)) +
      scale_fill_manual(values = rev(fill_color_set)[c(2, 4, 6)]) +
      # xlim(xstart, xend) +x
      xlab(x) +
      ylab(y)
    
    # warp the offset plot
    if (facet) {
      offset_plot <- offset_plot + 
        facet_wrap(`wrap_group`) +
        theme(strip.background = element_blank(),
              strip.text = element_text(size = font_size),
              panel.grid = element_blank())
    } 
    
    # ggplot2::ggsave(plot = offset_plot, filename = './output/frame_offset_bar_plot.png', width = 9, height = 5, dpi = 300)
    
  } else if(plot_type == "heat") {
    offset_plot <- ggplot(data = offset) +
      geom_tile(aes(x = !!sym(fill_group),
                    y = !!sym(x),
                    fill = !!sym(y)),
                alpha = fill_alpha) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.text = element_text(size = font_size - 3),
            legend.title = element_text(size = font_size - 2),
            axis.text = element_text(color = "black", size = font_size, hjust = 0.5, vjust = 0.5),
            axis.title = element_text(color = "black", size = font_size + 1, hjust = 0.5, vjust = 0.5),
            plot.title = element_text(color = "black", size = font_size + 2, hjust = 0.5, vjust = 0.5)) +
      scale_fill_gradient(low = fill_color_set[2], high = fill_color_set[8]) +
      # xlim(xstart, xend) +
      xlab(fill_group) +
      ylab(x)
    
    # warp the offset plot
    if (facet) {
      offset_plot <- offset_plot + 
        facet_grid(`wrap_group`) +
        theme(strip.background = element_blank(),
              strip.text = element_text(size = font_size),
              panel.grid = element_blank())
    }
  }
  
  # ggplot2::ggsave(plot = offset_plot, filename = './output/frame_offset_heat_plot.png', width = 6, height = 6, dpi = 300)
  
  return(offset_plot)
  
}