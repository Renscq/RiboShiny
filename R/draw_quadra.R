###########################################
# draw_quadra.R
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


draw_quadra <- function(
    degs = NULL,
    x = NULL,
    y = NULL,
    x_logfc = NULL,
    y_logfc = NULL,
    
    gene_class = "Class",
    
    remove_ns = TRUE, 
    
    dot_size = 3,
    color = "Set2",
    font_size = 12,
    
    x_max = NULL,
    x_min = NULL,
    y_max = NULL,
    y_min = NULL,
    
    title = "Quadrant Plot",
    xlabel = "RNA_log2FC",
    ylabel = "Ribo_log2FC",
    
    label_column = 'Gene',
    label_color = 'black',
    label_size = 3,
    label_list = NULL,
    overlaps_num = 10
    
) {
  
  ######################################################
  # set the color list
  color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = color))(9)
  colors <- color_set[1:9]
  
  ######################################################
  # Remove NA and NS
  # browser()

  degs <- degs %>% 
    tidyr::drop_na(DEGs)
  
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
  
  # Draw the plot
  quadra_plot <- ggplot(degs, aes(x = !!sym(x), y = !!sym(y), color = !!sym(gene_class))) +
    geom_point(size = dot_size) +
    scale_color_manual(values = colors) +
    theme_bw() +
    theme(axis.text.x = element_text(color = 'black', size = font_size),
          axis.text.y = element_text(color = 'black', size = font_size, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(color = 'black', size = font_size + 1, vjust = 0.5, hjust = 0.5),
          plot.title = element_text(color = 'black', size = font_size + 2, vjust = 0.5, hjust = 0.5),
          plot.subtitle = element_text(color = 'black', size = font_size - 2, vjust = 0.5, hjust = 0.5),
          legend.text = element_text(color = 'black', size = font_size - 3),
          legend.title = element_text(color = 'black', size = font_size - 2),
          panel.grid.minor = element_blank()) +
    geom_vline(xintercept = c(-x_logfc, x_logfc), linetype = "dashed", color = "black") +
    geom_hline(yintercept = c(-y_logfc, y_logfc), linetype = "dashed", color = "black") +
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
    quadra_plot <- quadra_plot +
      geom_text_repel(data = label_col_gene,
                      aes(x = !!sym(x), y = !!sym(y), label = !!sym(label_column)),
                      size = label_size, color = label_color, max.overlaps = overlaps_num)

  }
  
  return(quadra_plot)
}