# function to draw meta plot

library(ggplot2)
library(RColorBrewer)
library(tidyverse)


draw_metaplot <- function(
    meta = NULL,
    x = NULL,
    y = NULL,
    xstart = NULL,
    xend = NULL,
    
    ystart = NULL,
    yend = NULL,
    
    font_size = 12,
    line_color = "Blues", # only for the gene plot
    line_width = 0.8, # only for the line plot
    group = NULL, # group for fill
    fill_color = "Blues", # Blues, Reds, Greens ...
    fill_alpha = 1,
    line_alpha = 1,
    
    legend_position = "bottom",
    
    facet = NULL, 
    wrap_group = NULL, # group for facet_wrap
    scales = "free_y",
    ncol = NULL,
    nrow = NULL,
    plot_type = "Line" # line plot, bar plot or heat map
) {
  
  # browser()

  # draw the meta
  if (plot_type == "Line") {
    
    # retrieve the colors
    line_number <- length(unique(meta[[group]]))
    edge_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = line_color))(line_number)
    line_colors <- edge_color_set[1:line_number]
    
    meta_plot <- ggplot(data = meta) +
      geom_line(aes(x = !!sym(x), 
                    y = !!sym(y), 
                    color = !!sym(group)),
                alpha = line_alpha,
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
            legend.position = legend_position) +
      scale_color_manual(values = line_colors) +
      facet_wrap(~Meta, scales = "free_x") +
      scale_y_continuous(limits = c(ystart, yend)) +
      xlab(x) +
      ylab(y)
    
    # warp the meta plot
    if (facet) {
      meta_plot <- meta_plot +
        facet_grid(as.formula(paste0(wrap_group, " ~ Meta")),
                   scales = "free_x") +
        theme(strip.background = element_blank(),
              strip.text = element_text(color = 'black', size = font_size),
              panel.grid = element_blank()) 
    }
    
    # ggplot2::ggsave(plot = meta_plot, filename = './output/meta_line_plot.png', width = 9, height = 12, dpi = 300)
  } else if (plot_type == "Bar") {
    
    # retrieve the colors
    fill_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = fill_color))(9)
    
    meta_plot <- ggplot(data = meta) +
      geom_bar(aes(x = !!sym(x), 
                   y = !!sym(y), 
                   fill = !!sym(group),
                   group = !!sym(group)),
               alpha = fill_alpha,
               stat = "identity",
               position = "dodge") +
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.text = element_text(color = 'black', size = font_size),
            axis.text = element_text(color = 'black', size = font_size),
            axis.title = element_text(color = 'black', size = font_size + 1, vjust = 0.5, hjust = 0.5),
            plot.title = element_text(color = 'black', size = font_size + 2, vjust = 0.5, hjust = 0.5),
            legend.text = element_text(color = 'black', size = font_size - 3),
            legend.title = element_text(color = 'black', size = font_size - 2),
            panel.grid = element_blank(),
            legend.position = legend_position) +
      scale_fill_manual(values = rev(fill_color_set)[c(2, 5, 8)]) +
      facet_wrap(~Meta, scales = "free_x") +
      scale_y_continuous(limits = c(ystart, yend)) +
      xlab(x) +
      ylab(y)
    
    # warp the meta plot
    if (facet) {
      meta_plot <- meta_plot +
        facet_grid(as.formula(paste0(wrap_group, " ~ Meta")),
                   scales = "free_x") +
        theme(strip.background = element_blank(),
              strip.text = element_text(color = 'black', size = font_size),
              panel.grid = element_blank()) 
    }
    
    # ggplot2::ggsave(plot = meta_plot, filename = './output/meta_bar_plot.png', width = 9, height = 9, dpi = 300)
    
    
  } else if(plot_type == "Heat") {
    # retrieve the colors
    fill_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = fill_color))(8)
    
    # browser()
    
    meta_plot <- ggplot(data = meta) +
      geom_tile(aes(x = !!sym(x),
                    y = !!sym(y),
                    fill = !!sym(group)),
                alpha = fill_alpha) +
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.text = element_text(color = 'black', size = font_size),
            axis.text = element_text(color = 'black', size = font_size),
            axis.title = element_text(color = 'black', size = font_size + 1, vjust = 0.5, hjust = 0.5),
            plot.title = element_text(color = 'black', size = font_size + 2, vjust = 0.5, hjust = 0.5),
            legend.text = element_text(color = 'black', size = font_size - 3),
            legend.title = element_text(color = 'black', size = font_size - 2),
            legend.position = legend_position,
            panel.grid = element_blank()) +
      scale_fill_gradient(low = fill_color_set[2], high = fill_color_set[7]) +
      facet_wrap(~Meta, scales = "free_x") +
      # xlim(xstart, xend) +
      xlab(x) +
      ylab(group)
    
    # warp the meta plot
    if (facet) {
      meta_plot <- meta_plot +
        facet_grid(as.formula(paste0(wrap_group, " ~ Meta")),
                   scales = "free") +
        theme(strip.background = element_blank(),
              strip.text = element_text(color = 'black', size = font_size),
              panel.grid = element_blank()) 
    }
  }
  
  # ggplot2::ggsave(plot = meta_plot, filename = './output/meta_heat_plot.png', width = 9, height = 9, dpi = 300)
  
  return(meta_plot)
  
}