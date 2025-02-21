################################################
# draw codon usage boxplot plot
# 
# degs: the DEGs data
# log2fc_column: the log2FC column name
# pvalue_column: the pvalue column name
# log2fc: the log2FC cutoff
# pvalue: the pvalue cutoff

# class: the DEGs class
# up_color: the up DEGs color
# down_color: the down DEGs color
# ns_color: the non-significant DEGs color
# dot_size: the dot size
# font_size: the font size

# x_max: the x-axis max value
# x_min: the x-axis min value

# fill_0: fill the 0 Padj to 0.1 * min(Padj)

# vol_sqrt: scale the y-axis to sqrt
# vol_log2: scale the y-axis to log2

# vol_title: the title of the volcano plot
# vol_xlab: the x-axis label
# vol_ylab: the y-axis label

# label_col: the gene name column
# label_color: the gene name color
# label_size: the gene name size
# label_list: the gene name list
# label_num: the gene name number

require(ggplot2)
require(tidyverse)
require(RColorBrewer)


draw_cu_box <- function(
    codon_usage = NULL,
    
    x = "Codon",
    y = "Usage",

    box_color = "Paired",
    box_alpha = 0.8,
    box_width = 0.5,
    
    line_width = 0.2,
    line_color = "black",
    
    outlier_show = FALSE,
    outlier_color = "black",
    outlier_size = 3,

    font_size = 12,
    
    y_min = NA,
    y_max = NA,
    y_breaks = 5,

    y_sqrt = "none", # c("none", "log10", "sqrt")
    
    x_angle = 90,
    
    legend_position = "top", # c("none", "top", "bottom", "left", "right")
    
    cu_title = "Codon usage",
    cu_xlab = "CLass",
    cu_ylab = "Frquency"

) {
  

  # browser()
  ######################################################
  # retrieve the colors
  color_num <- length(unique(codon_usage[[x]]))
  box_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = box_color))(color_num)
  
  ######################################################
  ## set the ylim
  
  if (is.na(y_max)) {
    y_max <- max(codon_usage$Value)
  }
  
  if (is.na(y_min)) {
    y_min <- min(codon_usage$Value)
  }
  
  codon_usage <- codon_usage %>% 
    dplyr::filter(Value >= y_min, Value <= y_max)

  ######################################################
  # set the label type
  
  # browser()
  
  cu_boxplot <- ggplot(codon_usage, 
                       aes(x = !!sym(x), 
                           y = !!sym(y), 
                           fill = !!sym(x))) +
    geom_boxplot(
      width = box_width,
      alpha = box_alpha,
      size = line_width,
      color = line_color,
      outliers = outlier_show,
      outlier.shape = 19,
      outlier.size = outlier_size,
      outlier.color = outlier_color
    ) +
    scale_fill_manual(values = box_color_set) +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = font_size + 2, vjust = 0.5, hjust = 0.5, color = 'black'),
          plot.subtitle = element_text(size = font_size - 2, vjust = 0.5, hjust = 0.5, color = 'black'),
          axis.title = element_text(size = font_size + 1, color = "black"),
          axis.text.x = element_text(size = font_size, color = 'black', angle = x_angle, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = font_size, color = 'black'),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          legend.position = legend_position,
          legend.title = element_blank(),
          legend.key = element_rect(fill = 'transparent', colour = 'white'),
          legend.text = element_text(size = font_size - 2),
          legend.background = element_rect(fill = 'transparent', colour = 'white')) +
    scale_y_continuous(limits = c(y_min, y_max), n.breaks = y_breaks) +
    ggtitle(cu_title) +
    labs(x = cu_xlab, y = cu_ylab)
  
  ######################################################
  ## scale x/y sqrt
  # browser()
  if (y_sqrt == "sqrt") {
    cu_boxplot <- cu_boxplot + 
      scale_y_sqrt()
  } else if (y_sqrt == "log10") {
    cu_boxplot <- cu_boxplot + 
      scale_y_log10()
  }
  
  
  return(cu_boxplot)
  
}