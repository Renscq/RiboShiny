# function to draw reads period figure

library(ggplot2)
library(RColorBrewer)
library(tidyverse)

draw_period <- function(
    period = NULL,
    x = 'Sample',
    y = 'Count',
    title = "3nt periodicity",
    edge_color = 'grey85',
    edge_width = 0.5,
    fill_group = "Frame",
    fill_color = 'Set1',
    fill_alpha = 0.8,
    bar_width = 0.8,
    font_size = 12,
    label_size = 3,
    trans = FALSE,
    text_label = FALSE) {
  
  # browser()
  
  # retrieve the colors
  fill_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = fill_color))(9)
  
  # draw the figure
  p <- ggplot(period) +
    geom_bar(aes(x = !!sym(x),
                 y = !!sym(y),
                 fill = !!sym(fill_group)), 
             stat = 'identity', 
             position = "stack", 
             color = edge_color,
             linewidth = edge_width,
             alpha = fill_alpha, 
             width = bar_width) +
    theme_bw() + 
    theme(axis.text.y = element_text(color = 'black', size = font_size),
          axis.text.x = element_text(color = 'black', angle = 90, vjust = 0.5, hjust = 1, size = font_size),
          axis.title = element_text(size = font_size + 1),
          plot.title = element_text(vjust = 0.5, hjust = 0.5, size = font_size + 2),
          legend.text = element_text(size = font_size - 3),
          legend.title = element_text(size = font_size - 2),
          panel.grid = element_blank(),
    ) +
    scale_fill_manual(values = fill_color_set[c(2, 5, 8)]) +
    ggtitle(title) +
    labs(y = y, x = x)
  
  # ggplot2::ggsave(plot = p, filename = './output/period.png', width = 9, height = 6, dpi = 300)
  
  # label the data
  if (text_label) {
    
    if (y == "Count") {
      p <- p + geom_text(data = period, 
                         aes(x = !!sym(x), 
                             y = !!sym(y), 
                             label = Label), 
                         position = position_stack(vjust = 0.5), size = label_size)
    } else if(y == "Ratio") {
      p <- p + geom_text(data = period, 
                         aes(x = !!sym(x), 
                             y = !!sym(y), 
                             label = Ratio), 
                         position = position_stack(vjust = 0.5), size = label_size)
    }
  }
  # ggplot2::ggsave(plot = p, filename = './output/period_label.png', width = 9, height = 6, dpi = 300)
  
  
  if (trans) {
    p <- p + coord_flip()
  }
  
  # ggplot2::ggsave(plot = p, filename = './output/period_flip.png', width = 9, height = 6, dpi = 300)
  
  return(p)
}