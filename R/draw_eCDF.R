####################################################
# draw_eCDF
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

draw_ecdf <- function(
    in_data = NULL, 
    x = NULL,
    group = NULL, 
    
    title = "cumulative",
    xlabel = "x",
    ylabel = "eCDF",
    
    line_color = "Blues",
    line_alpha = 0.8,
    line_width = 0.9,
    
    x_min = NA,
    x_max = NA,
    x_breaks = 5,
    
    legend_position = "right", # c("none", "top", "bottom", "left", "right")
    
    log2_transform = FALSE,
    facet = FALSE,
    font_size = 12
    ) {

  ######################################################
  # retrieve the colors
  color_num <- length(unique(in_data[[group]]))
  line_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = line_color))(color_num)
  
  # set the color alpha = 0.8
  line_color_set <- alpha(line_color_set, line_alpha)
  
  # browser()
  
  if (is.na(x_min)) {
    x_min <- min(in_data[[x]])
  }
  
  if (is.na(x_max)) {
    x_max <- max(in_data[[x]])
  }
  
  # draw the line plot
  plots <- ggplot(data = in_data, 
                  mapping = aes(x = !!sym(x),
                                group = !!sym(group), 
                                color = !!sym(group))) +
    stat_ecdf(linewidth = line_width) +
    theme_bw() +
    theme(axis.text.x = element_text(color = 'black', size = font_size),
          axis.text.y = element_text(color = 'black', size = font_size, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(color = 'black', size = font_size + 1, vjust = 0.5, hjust = 0.5),
          plot.title = element_text(color = 'black', size = font_size + 2, vjust = 0.5, hjust = 0.5),
          legend.position = legend_position,
          legend.text = element_text(color = 'black', size = font_size - 3),
          legend.title = element_text(color = 'black', size = font_size - 2)) +
    scale_color_manual(values = line_color_set) +
    scale_x_continuous(limits = c(x_min, x_max), n.breaks = x_breaks) +
    labs(x = xlabel, y = ylabel) +
    ggtitle(title)
  
  # ggplot2::ggsave(plot = plots, filename = './output/cumulative.png', width = 9, height = 6, dpi = 300)
  
  # wrap the plots
  if (facet) {
    plots <- plots + 
      facet_wrap(`group`) +
      theme(strip.background = element_blank(),
            strip.text = element_text(color = 'black', size = font_size),
            legend.position = "none")
  }
  
  #######################
  return(plots)
}