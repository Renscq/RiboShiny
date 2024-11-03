###########################################
# draw_kegg_Terms.R
#
# This script enrich the KEGG terms for the DEGs
#
# kegg_results: the data frame of the KEGG terms
# x: the x-axis of the KEGG terms
# color: the color of the KEGG terms
# showCategory: the number of the KEGG terms to show
# dotsize: the size of the dots
# chr_width: the width of the KEGG terms
# dot_alpha: the alpha of the dots
# low_color: the low color of the KEGG terms
# high_color: the high color of the KEGG terms
# facet_group: whether to facet the KEGG terms
# gridwidth: the width of the grid
# font_size: the font size of the KEGG terms
# xmin: the minimum x-axis of the KEGG terms
# xmax: the maximum x-axis of the KEGG terms
# x_breaks: the breaks of the x-axis of the KEGG terms
# clustered: whether to cluster the KEGG terms
# title: the title of the KEGG terms
# xlab: the x-axis label of the KEGG terms
# ylab: the y-axis label of the KEGG terms
#
# usage: 


library(clusterProfiler)
library(enrichplot)
library(stringr)


draw_kegg_terms <- function(
    kegg_results = NULL,
    x = "GeneRatio",
    color = "p.adjust",
    
    showCategory = 10,
    
    dotsize = c(3, 9),
    chr_width = 80,
    
    dot_alpha = 0.9,
    low_color = "blue",
    high_color = "red",

    facet_group1 = "Cluster",
    facet_group2 = NULL,
    
    gridwidth = 0.4,
    font_size = 12,
    
    xmin = NA,
    xmax = NA,
    x_breaks = NA,
    
    clustered = FALSE,
    
    title = "KEGG Enrichment Analysis",
    xlab = "GeneRatio",
    ylab = "KEGG Pathway"
    
) {

  ######################################################
  # draw the GO terms
  kegg_plots <- enrichplot::dotplot(kegg_results, showCategory = showCategory)
  
  ######################################################
  # facet the GO terms
  # browser()
  
  if ("FALSE" %in% facet_group1 | length(facet_group1) == 0) {facet_group1 = FALSE} 
  
  if (!isFALSE(facet_group1) && clustered == TRUE) {
    
    if (inherits(kegg_results, 'gseaResult')) {
      if (facet_group1 == "vertical") {
        kegg_plots <- kegg_plots + 
          facet_grid(.~.sign, 
                     scales = "free",
                     space = "free")
        
      } else if (facet_group1 == "horizontal") {
        kegg_plots <- kegg_plots + 
          facet_grid(.sign~., 
                     scales = "free",
                     space = "free")
      }
      
    } else {
      if (facet_group1 == "vertical") {
        kegg_plots <- kegg_plots + 
          facet_grid(.~Cluster, 
                     scales = "free",
                     space = "free")
        
      } else if (facet_group1 == "horizontal") {
        kegg_plots <- kegg_plots + 
          facet_grid(Cluster~., 
                     scales = "free",
                     space = "free")
      }
    }

    
  } else if (!isFALSE(facet_group1) && clustered == FALSE) {
    
    equal_group <- all.equal(sort(facet_group1), sort(facet_group2))
    
    if (!inherits(equal_group, 'logical')) {
      
      overlap <- intersect(facet_group1, facet_group2)
      facet_group2 <- setdiff(facet_group2, overlap)
      
      # browser()
      
      if ("NULL" %in% facet_group1 | length(facet_group1) == 0) {facet_group1 = NULL}
      if ("NULL" %in% facet_group2 | length(facet_group2) == 0) {facet_group2 = NULL}
      
      if (!is.null(facet_group1) && is.null(facet_group2)) {
        facet_group1 = paste0(facet_group1, collapse = "+")
        facet_formula <- as.formula(paste0(facet_group1, '~ .'))
        
        kegg_plots <- kegg_plots + 
          facet_grid(facet_formula, 
                     scales = "free",
                     space = "free")
        
      } else if (is.null(facet_group1) && !is.null(facet_group2)) {
        facet_group2 = paste0(facet_group2, collapse = "+")
        facet_formula <- as.formula(paste0('. ~', facet_group2))
        
        kegg_plots <- kegg_plots + 
          facet_grid(facet_formula, 
                     scales = "free",
                     space = "free")
        
      } else if (!is.null(facet_group1) && !is.null(facet_group2)) {
        facet_group1 = paste0(facet_group1, collapse = "+")
        facet_group2 = paste0(facet_group2, collapse = "+")
        
        facet_formula <- as.formula(paste0(facet_group1, '~', facet_group2))
        
        kegg_plots <- kegg_plots + 
          facet_grid(facet_formula, 
                     scales = "free",
                     space = "free")
        
      }
      
    }

  } 
  
  ######################################################
  # set the colors
  kegg_plots <- kegg_plots +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = chr_width)) +
    scale_fill_gradient(low = alpha(low_color, alpha = dot_alpha), high = alpha(high_color, alpha = dot_alpha))
  
  ######################################################
  # set the theme
  kegg_plots <- kegg_plots +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 0.5, size = font_size),
          axis.text.y = element_text(color = "black", size = font_size),
          axis.title = element_text(color = "black", size = font_size + 1),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, size = font_size + 2),
          panel.grid = element_line(linewidth = gridwidth),
          strip.text = element_text(color = "black", size = font_size),
          strip.background = element_rect(fill = "white", color = "grey"),
          legend.text = element_text(size = font_size - 2)
    ) +
    labs(y = ylab, x = xlab, title = title)
  
  ######################################################
  # set the dot size
  kegg_plots <- kegg_plots +
    scale_size(range = dotsize)
  
  ######################################################
  # set the x-axis
  # browser()
  
  if (!is.na(xmin) && !is.na(xmax) && clustered == FALSE) {
    kegg_plots <- kegg_plots +
      scale_x_continuous(limits = c(xmin, xmax), n.breaks = x_breaks)
  }
  
  return(kegg_plots)
}

