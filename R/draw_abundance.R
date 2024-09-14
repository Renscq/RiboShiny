# normalise the abundance data and draw the abundance heatmap

library(ggplot2)
library(RColorBrewer)
library(tidyverse)


draw_abundance <- function(
    abundance = NULL,
    x = NULL,
    y = NULL,
    xstart = NULL,
    xend = NULL,
    line_color = "Blues", # only for the gene plot
    line_width = 0.8, # only for the line plot
    group = NULL, # group for fill
    fill_color = "Blues", # Blues, Reds, Greens ...
    fill_alpha = 1,
    facet = NULL, 
    wrap_group = NULL, # group for facet_wrap
    scales = "free_y",
    ncol = NULL,
    nrow = NULL,
    plot_type = "Line" # line plot, box plot
) {
  
  # browser()
  if (plot_type == "Line"){
    
    
  } else if (plot_type == "Box") {
    
  }
  
  
}


