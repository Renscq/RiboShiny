# function to draw coefficient of variation line plot

library(minpack.lm)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

draw_cov_line <- function(
    cov_table = NULL,
    x = 'log2Mean',
    y = 'log2CoV',
    
    group = 'Sample',
    
    facet = TRUE,
    wrap_num = 2,
    
    dot_color = 'Paired',
    dot_alpha = 0.2,
    dot_size = 1,
    
    line_color = 'Paired',
    line_alpha = 0.9,
    line_width = 1,
    
    font_size = 12,
    
    xstart = NA,
    xend = NA,
    xbreaks = 5,
    
    ystart = NA,
    yend = NA,
    ybreaks = 3,

    type = 'fitted', # scatter, fitted, merge
    
    xlabel = "log2(Mean)",
    ylabel = "log2(CoV)",
    title = "Coefficient of Variation"
    
    ) {

  # browser()
  
  # summary the number of class
  color_num <- nrow(unique(cov_table[group]))
  dot_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = dot_color))(color_num)
  line_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = line_color))(color_num)
  
  dot_color_set <- dot_color_set[1:color_num]
  line_color_set <- line_color_set[1:color_num]
  
  # set the nonlinear model
  nonlinear_model <- function(mu, alpha, beta) {
    return(0.5 * log2((beta / mu) + alpha))
  }
  
  # create the ggplot object
  plots <- ggplot(data = cov_table,
         mapping = aes(x = !!sym(x), y = !!sym(y), color = !!sym(group), fill = !!sym(group)))
  
  # draw the plot
  # draw the scatter plot
  if (type == "scatter") {
    plots <- plots + 
      geom_point(shape = 21, alpha = dot_alpha, colour = "white")

  # draw the fitted line plot
  } else if (type == "fitted") {
    plots <- plots +
      geom_smooth(method = "nlsLM", 
                  formula = y ~ nonlinear_model(2^(x), alpha, beta),
                  method.args = list(start = c(alpha = 1, beta = 1)), 
                  se = F, 
                  alpha = line_alpha,
                  linewidth = line_width)
    
  # draw the scatter and fitted line plot
  } else if (type == "merge") {
    plots <- plots + 
      geom_point(shape = 21, alpha = dot_alpha, colour = "white") +
      geom_smooth(method = "nlsLM", 
                  formula = y ~ nonlinear_model(2^(x), alpha, beta),
                  method.args = list(start = c(alpha = 1, beta = 1)), 
                  se = F, 
                  alpha = line_alpha,
                  linewidth = line_width)

  }

  # set the theme
  plots <- plots + 
    theme_bw() +
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = font_size - 1),
          plot.title = element_text(size = font_size + 2, vjust = 0.5, hjust = 0.5),
          axis.title = element_text(size = font_size),
          axis.text = element_text(size = font_size - 1))
  
  # set the dot and line color
  if (!is.na(dot_color)) {
    plots <- plots + scale_fill_manual(values = dot_color_set)
  }
  
  if (!is.na(line_color)) {
    plots <- plots + scale_color_manual(values = line_color_set)
  }
  
  # set the axis
  if (!is.na(xstart) && !is.na(xend) && !is.na(xbreaks)) {
    plots <- plots + 
      scale_x_continuous(limits = c(xstart, xend), n.breaks = xbreaks)
  }
  
  if (!is.na(ystart) && !is.na(yend) && !is.na(ybreaks)) {
    plots <- plots + 
      scale_y_continuous(limits = c(ystart, yend), n.breaks = ybreaks)
  }
  
  # wrap the plot
  if (facet) {
    facet_formula <- as.formula(paste("~", group))
    
    plots <- plots + 
      facet_wrap(facet_formula, ncol = wrap_num) +
      theme(strip.background = element_blank(),
            strip.text = element_text(size = font_size - 1))
  }
  
  # ggplot2::ggsave(plot = plots, filename = './images/cov_fitted.png', width = 9, height = 6, dpi = 300)
  
  return(plots)
}
