# function to draw reads offset

library(ggplot2)
library(ggplotify)
library(RColorBrewer)
library(tidyverse)
# library(pheatmap)

draw_offset_end <- function(
    offset = NULL,
    x = NULL,
    y = NULL,
    xstart = NULL,
    xend = NULL,
    ystart = NULL,
    yend = NULL,
    
    wrap_split = FALSE,
    wrap_group = NULL,  # group for facet_wrap
    ncol = NULL,
    nrow = NULL,
    
    font_size = 12,
    edge_color = "grey85",
    edge_width = 0.5,
    
    fill_group = NULL, # group for fill
    fill_color = "Blues", # Blues, Reds, Greens ...
    fill_alpha = 1,
    
    xbreaks = 3,
    ybreaks = 3,
    
    xlabel = 'Offset',
    ylabel = 'Length',
    title = 'TIS'
) {
  
  # retrieve the colors
  fill_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = fill_color))(9)
  
  # browser()
  
  # draw the offset
  end_offset_plot <- ggplot(offset) + 
    # geom_tile(aes(x = Offset, y = Length, fill = Scaled)) + 
    geom_tile(aes(x = !!sym(x), 
                  y = !!sym(y), 
                  fill = !!sym(fill_group)),
              alpha = fill_alpha,
              color = edge_color,
              size = edge_width) + 
    theme_bw() +
    theme(legend.text = element_text(size = font_size - 3),
          legend.title = element_text(size = font_size - 2),
          axis.text = element_text(color = "black", size = font_size, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(color = "black", size = font_size + 1, hjust = 0.5, vjust = 0.5),
          plot.title = element_text(size = font_size + 2, hjust = 0.5, vjust = 0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = font_size - 1),
          panel.grid = element_blank()
          ) +
    scale_fill_gradient(low = fill_color_set[2], high = fill_color_set[8]) +
    # scale_x_continuous(labels = c(xstart, xend), n.breaks = xbreaks) +
    # scale_y_continuous(limits = c(ystart, yend), n.breaks = ybreaks) +
    # xlim(xstart, xend) +
    facet_grid(~End, scales = 'free') +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(title)
  
  # facet_wrap(wrap_group, ncol = ncol, nrow = nrow) +
  if (isTRUE(wrap_split)) {
    end_offset_plot <- end_offset_plot + 
      facet_grid(as.formula(paste0(wrap_group, " ~ End")), scales = "free_x") +
      theme(strip.background = element_blank(),
            strip.text = element_text(size = font_size - 1))
  }
  
  # ggplot2::ggsave(plot = end_offset_plot, filename = './output/end_of_offset_heat_plot.png', width = 12, height = 6, dpi = 300)
  
  return(end_offset_plot)
  
}