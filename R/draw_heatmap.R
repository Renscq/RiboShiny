# ####################################################
# # draw the heatmap of gene expression
library(pheatmap)
# library(ComplexHeatmap)
# library(ggplotify)
library(RColorBrewer)


draw_heat <- function(
    in_dat = NULL,
    color = "RdYlBu",
    alpha = 0.8,
    border_color = NA,
    
    scale = "none",
    
    gaps_row = 0,
    gaps_col = 0,
    
    cluster_cols = T,
    cluster_rows = T,
    
    xlab = "",
    ylab = "",
    
    show_rownames = T,
    show_colnames = T,
    angle_col = 90,
    
    display_numbers = F,
    number_format = "%.2f",
    fontsize_number = 12,
    
    cluster_distance = "euclidean",
    cluster_method = "ward.D2",
    
    silent = T,
    fontsize = 12,
    title = "heatmap"
    ) {
  
  # browser()
  
  # set the colors
  my_color <- rev(colorRampPalette(RColorBrewer::brewer.pal(8, color))(100))
  my_color <- alpha(my_color, alpha)
  
  # draw the heat map
  heat_plot <- pheatmap::pheatmap(mat = in_dat %>% as.matrix(),
                                  clustering_distance_rows = cluster_distance,
                                  clustering_distance_cols = cluster_distance,
                                  
                                  clustering_method = cluster_method,
                                  
                                  cluster_cols = cluster_cols,
                                  cluster_rows = cluster_rows,
                                  
                                  scale = scale,
                                  show_rownames = show_rownames,
                                  show_colnames = show_colnames,
                                  
                                  angle_col = angle_col,
                                  
                                  gaps_row = gaps_row,
                                  gaps_col = gaps_col,
                                  
                                  fontsize = fontsize,
                                  color = my_color,
                                  alpha = alpha,
                                  border_color = border_color,
                                  na_col = "grey",
                                  
                                  display_numbers = display_numbers,
                                  number_format = number_format,
                                  
                                  xlab = xlab,
                                  ylab = ylab,
                                  main = title,
                                  silent = silent)
  
  # draw the heat map
  # heat_plot <- ComplexHeatmap::Heatmap(
  #   mat = in_dat,
  #   
  #   clustering_method_rows = "single",
  #   clustering_distance_rows = cluster_distance,
  #   
  #   clustering_method_columns = "single",
  #   clustering_distance_columns = cluster_distance,
  # 
  #   cluster_rows = cluster_rows,
  #   show_row_dend = T,
  #   
  #   cluster_columns = cluster_cols,
  #   show_column_dend = T,
  #   
  #   scale = scale,
  #   show_rownames = show_rownames,
  #   show_colnames = show_colnames,
  #   angle_col = angle_col,
  # 
  #   row_km = gaps_row,
  #   column_km = gaps_col,
  # 
  #   fontsize = fontsize,
  #   color = my_color,
  #   alpha = alpha,
  #   border = border_color,
  #   na_col = "grey",
  #   xlab = xlab,
  #   ylab = ylab,
  #   name = title,
  #   silent = silent) %>%
  #   as.ggplot()
  
  
  
  # ggsave(filename = "output/heat_plot.pdf", plot = heat_plot, width = 6, height = 8)
  
  return(heat_plot)

}