# function to draw the codon usage line plot
# draw_codon
# Codon table contains the following columns:
# \itemize{} 
# \item Codon: the codon sequence
# \item AntiCodon: the anticodon sequence
# \item AA: the amino acid
# \item Abbreviation: the abbreviation of the amino acid
# \item Name: the name of the amino acid
# \item Frequency: the frequency of the codon
# \item RSCU: the relative synonymous codon usage
# \item CAI: the codon adaptation index

# @param codon A \code{Dataframe} object.
# @param x A \code{character} object. The x-axis variable.
# @param y A \code{character} object. The y-axis variable.
# @param xlab A \code{character} object. The x-axis label.
# @param ylab A \code{character} object. The y-axis label.
# @param title A \code{character} object. The title of the plot.
# @param dot_color A \code{character} object. The color of the dot.
# @param dot_size A \code{numeric} object. The size of the dot.
# @param codon_usage A \code{Dataframe} object.
# @param cu_class A \code{character} object. cu_class can be "Frequency", "RSCU", "CAI".
# @param grid A \code{logical} object. Whether to split the plot into grids.
# @param grid_width A \code{numeric} object. The width of the grid line.
# @param grid_color A \code{character} object. The color of the grid line.
# @param wrap A \code{logical} object. Whether to wrap the plot.
# @param wrap_group A \code{character} object. The group to wrap the plot.
# @param wrap_row A \code{numeric} object. The number of rows to wrap the plot.
#
# @return Returns a \code{ggplot} object.
#
# @examples
#
# codon_freq <- draw_codon(codon = codon_table, cu_class = "Frequency")
# codon_rscu <- draw_codon(codon = codon_table, cu_class = "RSCU")
# codon_cai <- draw_codon(codon = codon_table, cu_class = "CAI")
#
# @importFrom ggplot2 ggplot geom_point theme_bw theme element_text labs facet_wrap
# @importFrom dplyr filter mutate
# 
# @export
# 

require(RColorBrewer)
require(tidyverse)
require(ggplot2)


draw_codon <- function(
    codon = NULL,
    x = "Codon",
    y = "PauseScore",
    
    xlab = "Codon",
    ylab = "PauseScore",
    title = "Codon pausing score",
    
    fill_color = "AA",
    color_palette = "Set1",
    
    dot_size = 2,
    line_size = 0.5,
    color_alpha = 0.8,
    
    stop_codon = FALSE,
    grid_width = 0.2,
    grid_color = 'grey85',
    
    font_size = 12,
    font_color = 'black',
    
    legend_row = 3,
    
    wrap = TRUE,
    wrap_group = "AA",
    wrap_row = 3,
    wrap_scales = "free") {
  
  # browser()
  
  # set the color palette
  sample_num <- length(unique(codon[[fill_color]]))
  
  color_list <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = color_palette))(sample_num)
  
  # filter the stop codon
  if (stop_codon) {
    codon <- codon %>% 
      dplyr::filter(Codon != "TAA" & Codon != "TAG" & Codon != "TGA")
  }
  
  # arrange the data
  codon <- codon %>% 
    # dplyr::arrange(AA, !!sym(y)) %>%
    dplyr::arrange(AA, Codon) %>% 
    dplyr::mutate(Codon = factor(Codon, levels = unique(Codon)))
  
  # draw the figure
  p <- ggplot(codon, 
              aes(x = !!sym(x), 
                  y = !!sym(y), 
                  color = !!sym(fill_color),
                  group = !!sym(fill_color))) +
    geom_point(size = dot_size, alpha = color_alpha) +
    geom_line(linewidth = line_size, alpha = color_alpha) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = font_size, color = font_color),
          axis.text.y = element_text(size = font_size, color = font_color),
          axis.title = element_text(size = font_size + 1, color = font_color),
          plot.title = element_text(hjust = 0.5, size = font_size + 2, color = font_color),
          legend.text = element_text(size = font_size - 3),
          legend.title = element_text(size = font_size - 2),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = grid_width, color = grid_color),
          legend.position = "bottom") +
    labs(x = xlab, y = ylab, title = title) +
    scale_color_manual(values = color_list) +
    guides(color = guide_legend(nrow = legend_row))
  
  # grid the figure
  if (wrap) {
    p <- p + 
      facet_wrap(formula(paste("~", wrap_group)), nrow = wrap_row, scales = wrap_scales) +
      # facet_grid(formula(paste("~", wrap_group)), scales = wrap_scales, space = "free") +
      theme(strip.background = element_blank(),
            strip.text = element_text(size = font_size))
  }
  
  # ggplot2::ggsave(plot = p, filename = './images/codon_usage.png', width = 9, height = 6, dpi = 300)
  
  return(p)
}