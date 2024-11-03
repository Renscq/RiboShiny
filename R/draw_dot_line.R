# function to draw reads length distribution

library(ggplot2)
library(RColorBrewer)
library(tidyverse)

draw_dot_line <- function(
    distr = NULL,
    x = 'Length',
    y = 'Count',
    titles = 'Length distribution',
    group = 'Sample',
    
    facet = TRUE,
    wrap_num = 2,
    
    fill_color = 'Blues',
    fill_alpha = 0.9,
    
    line_color = 'Paired',
    line_width = 1,
    
    dot_size = 1,
    font_size = 12,
    
    xstart = 20,
    xend = 40,
    xbreaks = 5,
    
    ystart = NA,
    yend = NA,
    ybreaks = 3,
    
    xpeak = NA,
    xpeak_color = 'Red',
    
    ypeak = NA,
    ypeak_color = 'Red',
    
    grid_width = 0.2,
    grid_color = 'grey80',
    
    type = 'line'
    ) {

  # get the max value
  max_value <- max(distr[y])

  # browser()
  # draw the plot
  if (type == "line") {
    # summary the number of class
    color_num <- nrow(unique(distr[group]))
    
    # retrieve the colors
    line_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = line_color))(color_num)
    
    line_color_set <- line_color_set[1:color_num]
    line_color_set <- alpha(line_color_set, fill_alpha)
    
    # draw the line plot
    plots <- ggplot(data = distr, 
                    mapping = aes(x = !!sym(x),
                                  y = !!sym(y), 
                                  group = !!sym(group), 
                                  color = !!sym(group))) +
      geom_line(linewidth = line_width, alpha = fill_alpha) +
      geom_point(size = dot_size, alpha = fill_alpha) +
      theme_bw() +
      theme(legend.text = element_text(size = font_size - 3),
            legend.title = element_text(size = font_size - 2),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(size = grid_width, color = grid_color),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, size = font_size + 2),
            axis.title = element_text(color = "black", hjust = 0.5, vjust = 0.5, size = font_size + 1),
            axis.text = element_text(color = "black", hjust = 0.5, vjust = 0.5, size = font_size)) +
      scale_color_manual(values = line_color_set) +
      # geom_vline(aes(xintercept = peak), size = 1, colour = "#BB0000", linetype = "dashed") + 
      # xlim(xstart, xend) +
      # scale_y_sqrt() +
      xlab(x) +
      ylab(y) +
      ggtitle(titles)
    
    if (!is.na(xstart) && !is.na(xend) && !is.na(xbreaks)) {
      plots <- plots + 
        scale_x_continuous(limits = c(xstart, xend), n.breaks = xbreaks)
    }
    
    if (!is.na(ystart) && !is.na(yend) && !is.na(ybreaks)) {
      plots <- plots + 
        scale_y_continuous(limits = c(ystart, yend), n.breaks = ybreaks)
    }
    
    # ggplot2::ggsave(plot = plots, filename = './images/reads_line_length.png', width = 9, height = 6, dpi = 300)
    
  } else if (type == "heat") {
    # retrieve the colors
    fill_color_set <- RColorBrewer::brewer.pal(n = 8, name = fill_color)
    
    # draw the heat plot
    plots <- ggplot(data = distr) +
      geom_tile(aes(x = !!sym(x),
                    y = !!sym(group),
                    fill = !!sym(y)),
                alpha = fill_alpha) +
      theme_bw() +
      theme(legend.text = element_text(size = font_size - 3),
            legend.title = element_text(size = font_size - 2),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, size = font_size + 2),
            axis.title = element_text(color = "black", size = font_size + 1, hjust = 0.5, vjust = 0.5),
            axis.text = element_text(color = "black", size = font_size, hjust = 0.5, vjust = 0.5),
            panel.grid = element_blank()) +
      scale_fill_gradient(low = fill_color_set[2], high = fill_color_set[7]) +
      xlab(x) +
      ylab(group) +
      ggtitle(titles)
    
    # ggplot2::ggsave(plot = plots, filename = './images/reads_heat_length.png', width = 9, height = 6, dpi = 300)
  }
  
  # label the peak
  if (!is.na(xpeak)) {
    plots <- plots + 
      geom_vline(aes(xintercept = xpeak), linewidth = 1, colour = xpeak_color, linetype = "dashed") +
      geom_text(x = xpeak, y = max_value, label = xpeak, color = xpeak_color, size = 4, hjust = 1, vjust = 0)
  }
  
  if (!is.na(ypeak)) {
    # browser()
    plots <- plots + 
      geom_hline(aes(yintercept = ypeak), linewidth = 1, colour = xpeak_color, linetype = "dashed") +
      geom_text(x = (xstart + xend) / 2, y = ypeak, label = ypeak, color = xpeak_color, size = 4, hjust = 1, vjust = 0)
  }
  
  # wrap the plot
  if (facet) {
    facet_formula <- as.formula(paste("~", group))
    
    plots <- plots + 
      facet_wrap(facet_formula, scales = "free_x", ncol = wrap_num) +
      theme(strip.background = element_blank(),
            strip.text = element_text(size = font_size - 1))
  }
  
  # ggplot2::ggsave(plot = plots, filename = './images/reads_heat_length.png', width = 9, height = 6, dpi = 300)
  
  return(plots)
}