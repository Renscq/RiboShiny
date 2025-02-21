library(tidyverse)
library(ggplot2)
library(RColorBrewer)

########################################################
# draw the RLE
draw_rle <- function(
    rle_mat = NULL,
    x = "Sample",
    y = "RLE",
    ymin = -4,
    ymax = 4,
    group = "Group",
    fill_color = "Blues",
    fill_alpha = 0.9,
    font_size = 12, 
    title = "Expression level"
) {
  
  color_num <- length(unique(rle_mat[[group]]))
  fill_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = fill_color))(color_num + 2)
  fill_colors <- fill_color_set[2:c(color_num + 1)]
  
  
  rle_plot <- ggplot(rle_mat) + 
    geom_boxplot(aes(x = !!sym(x),
                     y = !!sym(y),
                     fill = !!sym(group)),
                 alpha = fill_alpha,
                 outlier.shape = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(color = 'black', size = font_size, angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(color = 'black', size = font_size, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(color = 'black', size = font_size + 1, vjust = 0.5, hjust = 0.5),
          plot.title = element_text(color = 'black', size = font_size + 2, vjust = 0.5, hjust = 0.5),
          legend.text = element_text(color = 'black', size = font_size - 3),
          legend.title = element_text(color = 'black', size = font_size - 2),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    scale_fill_manual(values = fill_colors) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    ggtitle(title)
  
  # ggplot2::ggsave(plot = rle_plot, filename = './output/rle_plot.png', width = 9, height = 6, dpi = 300)
  
  return(rle_plot)
}