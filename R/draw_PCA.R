####################################################
# draw the PCA plot
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

require(FactoMineR)
require(factoextra)
require(ggrepel)
require(cowplot)


draw_pca <- function(
    in_data = NULL, 
    
    phenodata = NULL, 
    group = NULL, 
    
    fill_color = "Blues",
    fill_alpha = 0.9,
    dot_size = 5,
    
    x_min = NULL,
    x_max = NULL,
    
    y_min = NULL,
    y_max = NULL,
    
    gene_num = 50,
    
    label_sample = T,
    label_size = 3,
    font_size = 12,
    title = "PCA"
    
    ) {
  
  # browser()
  
  ######################################################
  # set the color list
  color_num <- length(unique(phenodata[[group]]))
  fill_color_set <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = fill_color))(color_num + 2)
  fill_colors <- fill_color_set[2:c(color_num + 1)]
  
  ######################################################
  # main function
  k = rep(colnames(in_data), 1)
  res.pca <- FactoMineR::PCA(in_data, graph = FALSE, scale.unit = T)
  res.f.pca <- FactoMineR::PCA(in_data, graph = FALSE, scale.unit = F)
  
  ######################################################
  # Eigenvalues / Variances
  eig_plot <- fviz_eig(res.f.pca, addlabels = TRUE) +
    theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = font_size + 2),
          axis.title = element_text(hjust = 0.5, vjust = 0.5, size = font_size + 1),
          axis.text = element_text(hjust = 0.5, vjust = 0.5, size = font_size))
  
  # var_plot <- fviz_pca_var(res.f.pca, col.var = "cos2") +
  #   theme(plot.title = element_text(vjust = 0.5, hjust = 0.5),
  #         axis.title = element_text(hjust = 0.5, vjust = 0.5),
  #         axis.text = element_text(hjust = 0.5, vjust = 0.5))
  # eig_var_plot <- cowplot::plot_grid(eig_plot, var_plot, nrow = 1, align = "hv")
  
  # ggsave(filename = "output/eig_var_plot.pdf", plot = eig_var_plot, width = 8, height = 4)
  
  pca_res <- list(
    eigenvalue = res.f.pca$eig %>% as.data.frame(),
    var.coord = res.f.pca$var$coord %>% as.data.frame(),
    var.cor = res.f.pca$var$cor %>% as.data.frame(),
    var.cos2 = res.f.pca$var$cos2 %>% as.data.frame(),
    var.contrib = res.f.pca$var$contrib %>% as.data.frame(),
    ind.coord = res.f.pca$ind$coord %>% as.data.frame(),
    ind.cos2 = res.f.pca$ind$cos2 %>% as.data.frame(),
    ind.contrib = res.f.pca$ind$contrib %>% as.data.frame()
  )
  
  ######################################################
  # set the x and y axis range
  if (!is.numeric(x_min)) {
    x_min <- min(res.f.pca$var$coord[, 1]) - abs(min(res.f.pca$var$coord[, 1]) * 0.1)
  }
  
  if (!is.numeric(x_max)) {
    x_max <- max(res.f.pca$var$coord[, 1]) + abs(max(res.f.pca$var$coord[, 1]) * 0.1)
  }
  
  if (!is.numeric(y_min)) {
    y_min <- min(res.f.pca$var$coord[, 2]) - abs(min(res.f.pca$var$coord[, 2]) * 0.1)
  }
  
  if (!is.numeric(y_max)) {
    y_max <- max(res.f.pca$var$coord[, 2]) + abs(max(res.f.pca$var$coord[, 2]) * 0.1)
  }

  # x_min <- min(res.f.pca$var$coord[, 1]) - abs(min(res.f.pca$var$coord[, 1]) * 0.1)
  # x_max <- max(res.f.pca$var$coord[, 1]) + abs(max(res.f.pca$var$coord[, 1]) * 0.1)
  # y_min <- min(res.f.pca$var$coord[, 2]) - abs(min(res.f.pca$var$coord[, 2]) * 0.1)
  # y_max <- max(res.f.pca$var$coord[, 2]) + abs(max(res.f.pca$var$coord[, 2]) * 0.1)

  ######################################################
  # link the sample group to the PCA result
  pca_res$var.coord <- pca_res$var.coord %>% 
    rownames_to_column(var = "Sample") %>%
    left_join(phenodata, by = "Sample")
  
  # browser()
  ######################################################
  # draw the PCA plot
  pca_scatter <- ggplot(pca_res$var.coord) +
    geom_point(aes(x = Dim.1, 
                   y = Dim.2, 
                   color = !!sym(group)),
               alpha = fill_alpha,
               size = dot_size) +
    # scale_color_brewer(palette = 'Set2') +
    scale_color_manual(values = fill_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(color = 'black', size = font_size),
          axis.text.y = element_text(color = 'black', size = font_size, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(color = 'black', size = font_size + 1, vjust = 0.5, hjust = 0.5),
          plot.title = element_text(color = 'black', size = font_size + 2, vjust = 0.5, hjust = 0.5),
          legend.text = element_text(color = 'black', size = font_size - 3),
          legend.title = element_text(color = 'black', size = font_size - 2),
          panel.grid.minor = element_blank()) +
    xlab(paste0("Dim1 (", round(pca_res$eigenvalue[1, 2], 2), "%)")) +
    ylab(paste0("Dim2 (", round(pca_res$eigenvalue[2, 2], 2 ), "%)")) +
    ggtitle(title) +
    xlim(x_min, x_max) +
    ylim(y_min, y_max)
  
  # browser()
  
  # label the sample name
  if (label_sample) {
    pca_scatter <- pca_scatter +
      ggrepel::geom_text_repel(data = pca_res$var.coord,
                aes(x = Dim.1,
                    y = Dim.2,
                    label = Sample),
                size = label_size,
                alpha = fill_alpha,
                max.overlaps = 20)
    
  }
  
  # ggsave(filename = "output/PCA_scatter.pdf", plot = pca_scatter, width = 8, height = 6)
  
  #######################
  # draw the contribution plot
  pca_contrib_var <- fviz_contrib(res.f.pca, choice = "var", axes = 1:2) +
    coord_flip() +
    theme_light() +
    theme(axis.text = element_text(size = font_size),
          axis.title = element_text(vjust = 0.5, hjust = 0.5, size = font_size + 1),
          plot.title = element_text(vjust = 0.5, hjust = 0.5, size = font_size + 2))
  
  pca_contrib_ind <- fviz_contrib(res.f.pca, choice = "ind", axes = 1:2, top = gene_num) +
    coord_flip() +
    theme_light() +
    theme(axis.text = element_text(size = font_size),
          axis.title = element_text(vjust = 0.5, hjust = 0.5, size = font_size + 1),
          plot.title = element_text(vjust = 0.5, hjust = 0.5, size = font_size + 2))
  
  #######################
  return(list(eig_plot = eig_plot,
              pca_scatter = pca_scatter, 
              pca_contrib_var = pca_contrib_var,
              pca_contrib_ind = pca_contrib_ind,
              pca_res = pca_res))
}


