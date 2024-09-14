# function to draw reads alignment figure

library(ggplot2)
library(RColorBrewer)
library(tidyverse)

draw_alignment <- function(
    align = NULL,
    x = 'Sample',
    y = 'Count',
    group = 'Database',
    title = "Reads alignment",
    edge_color = 'grey85',
    fill_color = 'Set1',
    fill_alpha = 0.8,
    bar_width = 0.8,
    trans = FALSE,
    text_label = FALSE,
    font_size = 12,
    label_size = 3
) {
  
  # browser()
  
  # draw the figure
  p <- ggplot(align) +
    geom_bar(aes(x = !!sym(x),
                 y = !!sym(y),
                 fill = !!sym(group)), 
             stat = 'identity', position = "stack", color = edge_color,
             alpha = fill_alpha, width = bar_width) +
    theme_bw() + 
    theme(axis.text.x = element_text(size = font_size, color = 'black', angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = font_size, color = 'black'),
          axis.title = element_text(vjust = 0.5, hjust = 0.5, size = font_size + 1),
          plot.title = element_text(vjust = 0.5, hjust = 0.5, size = font_size + 2),
          legend.text = element_text(size = font_size - 2),
          legend.title = element_text(size = font_size - 1),
          panel.grid = element_blank()) +
    scale_fill_brewer(palette = fill_color) +
    ggtitle(title) +
    labs(y = y, x = x)
  
  # ggplot2::ggsave(plot = p, filename = './images/alignment.png', width = 9, height = 6, dpi = 300)
  
  if (text_label) {
    data = align %>% 
      dplyr::mutate(Database = factor(Database, levels = rev(unique(Database)))) %>% 
      dplyr::arrange(!!sym(x), Database)
    if (y == "Count"){
      p <- p + geom_text(data = data, 
                         aes(x = !!sym(x), 
                             y = !!sym(y), 
                             label = Label), 
                         size = label_size, 
                         position = position_stack(vjust = 0.5))
    } else if(y == "Ratio"){
      p <- p + geom_text(data = data, 
                         aes(x = !!sym(x), 
                             y = !!sym(y), 
                             label = Ratio), 
                         size = label_size, 
                         position = position_stack(vjust = 0.5))
    }
  }
  
  # ggplot2::ggsave(plot = p, filename = './images/alignment_label.png', width = 9, height = 6, dpi = 300)
  
  
  if (trans) {
    p <- p + coord_flip()
  }
  
  # ggplot2::ggsave(plot = p, filename = './images/alignment_flip.png', width = 9, height = 6, dpi = 300)
  
  return(p)
}