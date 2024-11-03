# function to draw the codon usage line plot
# draw_codon_usage
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
# 
# @param codon_usage A \code{Dataframe} object.
# @param cu_class A \code{character} object. cu_class can be "Frequency", "RSCU", "CAI".
# @param grid A \code{logical} object. Whether to split the plot into grids.
# @param grid_width A \code{numeric} object. The width of the grid line.
# @param grid_color A \code{character} object. The color of the grid line.
# @param wrap A \code{logical} object. Whether to wrap the plot.
# @param wrap_row A \code{numeric} object. The number of rows to wrap the plot.
#
# @return Returns a \code{ggplot} object.
#
# @examples
#
# codon_freq <- draw_codon_usage(codon_usage = codon_table, cu_class = "Frequency")
# codon_rscu <- draw_codon_usage(codon_usage = codon_table, cu_class = "RSCU")
# codon_cai <- draw_codon_usage(codon_usage = codon_table, cu_class = "CAI")
#
# @importFrom ggplot2 ggplot geom_point theme_bw theme element_text labs facet_wrap
# @importFrom dplyr filter mutate
# 
# @export
# 

draw_codon_usage <- function(
    codon_usage = NULL,
    cu_class = "Frequency",
    stop_codon = FALSE,
    
    grid_width = 0.2,
    grid_color = 'grey85',
    font_size = 12,
    font_color = 'black',
    
    dot_size = 2,
    color_map = "Paired",
    color_alpha = 0.8,
    
    line_width = 0.8,
    legend_row = 3,
    
    wrap = TRUE,
    wrap_row = 3,
    wrap_scales = "free",
    
    type = "dot" # Dot, Circle
    
    ) {

  # filter the data
  if (!stop_codon) {
    codon_usage <- codon_usage %>% 
      dplyr::filter(AA != "TER")
  }
  
  codon_usage <- codon_usage %>% 
    arrange(AA, Codon) %>% 
    mutate(Codon = factor(Codon, levels = unique(Codon)))

  # set the color palette
  # browser()
  
  color_num <- codon_usage %>% distinct(AA) %>% nrow()
  color_list <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = color_map))(color_num)
  
  if (type == "dot") {
    # draw the codon usage dot plot
    cu_plot <- ggplot(codon_usage, aes(x = Codon, 
                                           y = !!sym(cu_class), 
                                           color = AA, 
                                           group = AA)) +
      geom_line(linewidth = line_width, alpha = color_alpha) +
      geom_point(size = dot_size, alpha = color_alpha) +
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
      scale_color_manual(values = color_list) +
      labs(x = "Codon", y = cu_class, title = "Codon usage") +
      guides(color = guide_legend(nrow = legend_row))
  
  } else if (type == "circle") {
    
    # browser()
    
    # draw the codon usage circle plot
    codon_usage$row <- c(1:nrow(codon_usage))
    codon_usage$angle <- 90 - codon_usage$row * (360 / nrow(codon_usage))
    
    cu_plot <- ggplot(codon_usage, aes(x = Codon, 
                                       y = !!sym(cu_class), 
                                       fill = AA, 
                                       group = AA)) +
      geom_bar(stat = "identity", width = 0.8, alpha = color_alpha) +
      coord_polar(start = 0) +
      ylim(-max(codon_usage[cu_class]), max(codon_usage[cu_class])) +
      geom_text(aes(x = Codon, y = max(codon_usage[cu_class]), label = Codon, angle = angle), 
                size = 4, hjust = 0) +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(size = font_size, color = font_color),
            axis.title = element_text(size = font_size + 1, color = font_color),
            plot.title = element_text(hjust = 0.5, size = font_size + 2, color = font_color),
            legend.text = element_text(size = font_size - 3),
            legend.title = element_text(size = font_size - 2),
            panel.border = element_blank(),
            plot.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(linewidth = grid_width, color = grid_color),
            legend.position = "bottom") +
      scale_fill_manual(values = color_list) +
      labs(x = "", y = "", title = cu_class) +
      guides(fill = guide_legend(nrow = legend_row))
  }

  # browser()
  
  # grid the figure
  if (wrap) {
    cu_plot <- cu_plot + 
      facet_wrap(~AA, nrow = wrap_row, scales = wrap_scales) +
      # facet_grid(~AA, scales = wrap_scales, space = "free") +
      theme(strip.background = element_blank(),
            strip.text = element_text(size = font_size))
  }
  
  # ggplot2::ggsave(plot = cu_plot, filename = './images/codon_usage.png', width = 9, height = 6, dpi = 300)
  
  return(cu_plot)
}