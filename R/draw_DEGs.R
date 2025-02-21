##########################################################
# draw_DEGs.R
# 
# this script is used to draw the DEGs bar plot
#

require(ggplot2)


draw_degs_bar <- function(
  degs = NULL,
  x = NULL,
  y = NULL,
  fill_degs = NULL,
  
  up_color = "red",
  down_color = "blue",
  alpha = 1,
  edge_color = 'white',
  bar_width = 1,
  
  bar_degs = NULL,
  bar_group = NULL,
  
  facet_groups = NULL,
  font_size = 12,
  
  label = FALSE,
  label_size = 3,
  label_color = "black",
  
  sqrt_y = FALSE,
  
  title = "",
  xlab = "",
  ylab = ""
) {
  
  # browser()
  # create a ggplot object
  degs_bar <- ggplot(degs, aes(x = !!sym(x), y = !!sym(y))) + 
    geom_col(aes(fill = !!sym(fill_degs)), width = bar_width, alpha = alpha) + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = font_size + 2, vjust = 0.5, hjust = 0.5, color = 'black'),
          plot.subtitle = element_text(size = font_size - 2, vjust = 0.5, hjust = 0.5, color = 'black'),
          axis.title = element_text(size = font_size + 1, color = "black"),
          axis.text.x = element_text(size = font_size, color = 'black', angle = 90, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = font_size, color = 'black', angle = 0, hjust = 1),
          panel.background = element_rect(color = 'white', fill = 'transparent'),
          strip.background = element_rect(fill = 'transparent'),
          strip.text = element_text(size = font_size - 2, color = 'black'),
          legend.title = element_blank(),
          legend.key = element_rect(fill = 'transparent'),
          legend.text = element_text(size = font_size - 2),
          legend.background = element_rect(fill = 'transparent')) +
    scale_fill_manual(values = c(down_color, up_color), guide = FALSE) +
    ggtitle("DEGs") +
    labs(x = xlab,
         y = ylab) + 
    facet_grid(as.formula(paste0(bar_degs, "~", bar_group)),
               scales = "free")

  # add the facet groups
  if (!is.null(label)) {
    degs_bar <- degs_bar +
      geom_text(aes(label = !!sym(y)), 
                position = position_stack(vjust = 0.5), 
                size = label_size)
  }
  
  # scale the y axis
  if (sqrt_y) {
    degs_bar <- degs_bar + scale_y_sqrt()
  }
  
  # return the ggplot object
  return(degs_bar)
  
}
