# function to draw meta enrich

library(ggplot2)
library(RColorBrewer)
library(tidyverse)


draw_meta_enrich <- function(
    meta = NULL,
    x = NULL,
    y = NULL,
    xstart = NULL,
    xend = NULL,
    
    font_size = 12,
    line_color = "Blues", # only for the gene plot
    line_width = 0.8, # only for the line plot
    cl_alpha = 0.5, # only for the confidence interval
    
    group = NULL, # group for fill
    facet = NULL, 
    wrap_group = NULL, # group for facet_wrap
    scales = "free_y",
    
    ncol = NULL,
    nrow = NULL,
    
    xlabel = "Codon",
    ylabel = "Enrichment [A.U.]"
) {
  
  # browser()
  
  # draw the meta
  # retrieve the colors
  line_number <- length(unique(meta[[group]]))
  edge_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = line_color))(line_number)
  line_colors <- edge_color_set[1:line_number]
  
  meta_plot <- ggplot(meta,
                      aes(x = !!sym(x), 
                          y = !!sym(y), 
                          group = !!sym(group), 
                          color = !!sym(group), 
                          fill = !!sym(group))) +
    stat_summary(fun.data = "mean_cl_boot", geom = "smooth", se = TRUE, alpha = cl_alpha) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_text(color = 'black', size = font_size),
          axis.text = element_text(color = 'black', size = font_size),
          axis.title = element_text(color = 'black', size = font_size + 1, vjust = 0.5, hjust = 0.5),
          plot.title = element_text(color = 'black', size = font_size + 2, vjust = 0.5, hjust = 0.5),
          legend.text = element_text(color = 'black', size = font_size - 3),
          legend.title = element_text(color = 'black', size = font_size - 2),
          panel.grid = element_blank(),
          legend.position = "bottom") +
    scale_color_manual(values = line_colors) +
    scale_fill_manual(values = line_colors) +
    facet_wrap(~Meta, scales = "free_x") +
    xlab(xlabel) +
    ylab(ylabel)
  
  # warp the meta plot
  if (facet) {
    meta_plot <- meta_plot +
      facet_grid(as.formula(paste0(wrap_group, " ~ Meta")),
                 scales = "free_x") +
      theme(strip.background = element_blank(),
            strip.text = element_text(color = 'black', size = font_size),
            panel.grid = element_blank()) 
  }
  
  # ggplot2::ggsave(plot = meta_plot, filename = './output/meta_heat_plot.png', width = 9, height = 9, dpi = 300)
  
  return(meta_plot)
  
}