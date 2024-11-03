###########################################
# draw_delta.R
#
# This script draws a quadratic function
# using ggplot2
#
# degs: the data frame of the DEGs
# x: the x-axis of the plot

require(ggplot2)
require(ggrepel)
require(tidyverse)
require(RColorBrewer)


draw_delta <- function(
    degs = NULL,
    x = NULL,
    y = NULL,
    x_logfc = NULL,
    y_logfc = NULL,
    Delta = "Delta",
    
    remove_others = FALSE, 
    remove_ns = FALSE, 
    
    dot_size = 1.5,
    color = "Set2",
    font_size = 12,
    
    x_max = NULL,
    x_min = NULL,
    y_max = NULL,
    y_min = NULL,
    
    title = "Delta Plot",
    xlabel = "RNA_log2FC",
    ylabel = "Ribo_log2FC",
    
    label_column = "Gene",
    label_color = 'black',
    label_size = 3,
    label_list = NULL,
    overlaps_num = 10
    
) {
  
  ######################################################
  # set the color list
  color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = color))(8)
  colors <- color_set[1:4]
  
  ######################################################
  # Remove NA and NS
  # browser()
  
  degs <- degs %>% 
    tidyr::drop_na(DEGs) %>%
    dplyr::arrange(desc(Delta)) %>%
    dplyr::mutate(Delta = as_factor(Delta))
    
  if (remove_others) {
    degs <- degs %>%
      dplyr::filter(Delta != "others")
  }
  
  if (remove_ns) {
    degs <- degs %>%
      dplyr::filter(DEGs != "NS")
  }
  
  ######################################################
  # browser()
  
  # set the x and y axis range
  if (!is.numeric(x_min)) {
    x_min <- min(degs[[x]] - degs[[x]] * 0.1)
  }
  
  if (!is.numeric(x_max)) {
    x_max <- max(degs[[x]] + degs[[x]] * 0.1)
  }
  
  if (!is.numeric(y_min)) {
    y_min <- min(degs[[y]] - degs[[y]] * 0.1)
  }
  
  if (!is.numeric(y_max)) {
    y_max <- max(degs[[y]] + degs[[y]] * 0.1)
  }
  
  ######################################################
  # calculate the correlation
  # browser()
  corr <- cor.test(degs[[x]], degs[[y]])
  R_value <- corr$estimate
  P_value <- format(corr$p.value, scientific = TRUE, digits = 3)
  
  subtitle <- paste("N = ", nrow(degs), ", ",
                    "R = ", round(R_value, 3), ", ", 
                    "P = ", P_value)
  
  ######################################################
  # browser()
  
  # set the color for the plot
  if (length(unique(degs[[Delta]])) > 4) {
    delta_color <- c('grey90', color_set[1:4])
  } else {
    delta_color <- c(color_set[1:4])
  }
  
  # Draw the plot
  delta_plot <- ggplot(degs, aes(x = !!sym(x), y = !!sym(y), color = !!sym(Delta))) +
    geom_point(size = dot_size) +
    scale_color_manual(values = delta_color) +
    theme_bw() +
    theme(axis.text.x = element_text(color = 'black', size = font_size),
          axis.text.y = element_text(color = 'black', size = font_size, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(color = 'black', size = font_size + 1, vjust = 0.5, hjust = 0.5),
          plot.title = element_text(color = 'black', size = font_size + 2, vjust = 0.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size = font_size - 2, vjust = 0.5, hjust = 0.5),
          legend.text = element_text(color = 'black', size = font_size - 3),
          legend.title = element_text(color = 'black', size = font_size - 2),
          panel.grid.minor = element_blank()) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggtitle(title, subtitle = subtitle) +
    labs(x = xlabel, y = ylabel) +
    xlim(x_min, x_max) +
    ylim(y_min, y_max)
  
  ######################################################
  # Add labels
  # browser()
  
  if (length(label_list) > 0 & label_column != "") {
    label_col_gene <- degs %>%
      dplyr::filter(!!sym(label_column) %in% label_list)
    
    ## add the gene name
    delta_plot <- delta_plot +
      geom_text_repel(data = label_col_gene,
                      aes(x = !!sym(x), y = !!sym(y), label = !!sym(label_column)),
                      size = label_size, color = label_color, max.overlaps = overlaps_num)
    
  }
  
  return(delta_plot)
}