# function to draw the length distribution of the sequences

draw_cds_len_distr <- function(
    distr = NULL,
    label = 'length',
    title = "CDS length distribution",
    edge_color = '#1155FF',
    fill_color = '#1155FF',
    fill_alpha = 0.9,
    xmin = 0,
    xmax = 10000,
    binwidth = 10,
    sqrty = FALSE,
    font_size = 12) {

  # draw the figure
  p <- ggplot(distr) +
    geom_histogram(aes_string(x = label), 
                   binwidth = binwidth, 
                   color = edge_color, 
                   fill = fill_color,
                   alpha = fill_alpha) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = font_size, color = 'black'),
          axis.text.y = element_text(size = font_size, color = 'black'),
          axis.title = element_text(size = font_size + 1, color = 'black'),
          plot.title = element_text(hjust = 0.5, size = font_size + 2, color = 'black'),
          legend.text = element_text(size = font_size - 3),
          legend.title = element_text(size = font_size - 2),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = "bottom") +
    xlim(xmin, xmax) +
    labs(x = "Length", y = "Count", title = title)

  # grid the figure
  if (sqrty) {
    p <- p + scale_y_sqrt()
  }
  
  # ggplot2::ggsave(plot = p, filename = './images/codon_usage.png', width = 9, height = 6, dpi = 300)
  
  return(p)
}