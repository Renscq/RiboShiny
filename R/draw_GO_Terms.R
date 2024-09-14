###########################################
# draw_GO_Terms.R
#
# This script enrich the GO terms for the DEGs
#
# go_results: the data frame of the GO terms
# x: the x-axis of the GO terms
# color: the color of the GO terms
# showCategory: the number of the GO terms to show
# dotsize: the size of the dots
# chr_width: the width of the GO terms
# dot_alpha: the alpha of the dots
# low_color: the low color of the GO terms
# high_color: the high color of the GO terms
# facet_group: whether to facet the GO terms
# gridwidth: the width of the grid
# font_size: the font size of the GO terms
# xmin: the minimum x-axis of the GO terms
# xmax: the maximum x-axis of the GO terms
# x_breaks: the breaks of the x-axis of the GO terms
# clustered: whether to cluster the GO terms
# title: the title of the GO terms
# xlab: the x-axis label of the GO terms
# ylab: the y-axis label of the GO terms
#
# usage: run_go_terms(degs, orgdb, key_type, pvalueCutoff, qvalueCutoff, ontology, pAdjustMethod, readable)


library(clusterProfiler)
library(enrichplot)
library(stringr)


draw_go_terms <- function(
    go_results = NULL,
    x = "GeneRatio",
    color = "p.adjust",
    
    showCategory = 10,
    
    dotsize = c(3, 9),
    chr_width = 80,
    
    dot_alpha = 0.9,
    low_color = "blue",
    high_color = "red",

    ontology = "ALL",
    
    facet_group1 = ".",
    facet_group2 = "ONTOLOGY",
    
    gridwidth = 0.4,
    font_size = 12,
    
    xmin = NULL,
    xmax = NULL,
    x_breaks = NULL,
    
    clustered = FALSE,
    
    title = "GO Enrichment Analysis",
    xlab = "GeneRatio",
    ylab = "Go Terms"
    
) {
  
  ######################################################
  # draw the GO terms
  # browser()
  
  if ("FALSE" %in% ontology | length(ontology) == 0) {ontology = FALSE}
  
  if (!isTRUE(ontology)) {
    if ("ALL" %in% ontology) {
      go_plots <- enrichplot::dotplot(go_results, 
                                      split = "ONTOLOGY",
                                      showCategory = showCategory)
    } else {
      go_plots <- enrichplot::dotplot(go_results %>% 
                                        dplyr::filter(ONTOLOGY %in% ontology), 
                                      split = "ONTOLOGY",
                                      showCategory = showCategory)
    }

  } else {
    go_plots <- enrichplot::dotplot(go_results, 
                                    showCategory = showCategory)
  }

  ######################################################
  # facet the GO terms
  # browser()
  
  if (inherits(go_results, 'compareClusterResult')) {
    
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
        
        go_plots <- go_plots + 
          facet_grid(facet_formula, 
                     scales = "free",
                     space = "free")
        
      } else if (is.null(facet_group1) && !is.null(facet_group2)) {
        facet_group2 = paste0(facet_group2, collapse = "+")
        facet_formula <- as.formula(paste0('. ~', facet_group2))
        
        go_plots <- go_plots + 
          facet_grid(facet_formula, 
                     scales = "free",
                     space = "free")
        
      } else if (!is.null(facet_group1) && !is.null(facet_group2)) {
        facet_group1 = paste0(facet_group1, collapse = "+")
        facet_group2 = paste0(facet_group2, collapse = "+")
        
        facet_formula <- as.formula(paste0(facet_group1, '~', facet_group2))
        
        go_plots <- go_plots + 
          facet_grid(facet_formula, 
                     scales = "free",
                     space = "free")
        
      }

    }
  } else if (inherits(go_results, 'gseaResult')) {
    
    if ("FALSE" %in% facet_group1 | length(facet_group1) == 0) {facet_group1 = FALSE} 

    if (!isFALSE(facet_group1)) {
      if (facet_group1 == "vertical") {
        facet_formula <- as.formula(paste0(".~", facet_group2))
        
        go_plots <- go_plots + 
          facet_grid(facet_formula, 
                     scales = "free",
                     space = "free")
        
      } else if (facet_group1 == "horizontal") {
        facet_formula <- as.formula(paste0(facet_group2, "~."))
        
        go_plots <- go_plots + 
          facet_grid(facet_formula, 
                     scales = "free",
                     space = "free")
      }
    
    }
  }
  
  ######################################################
  # set the colors
  # browser()
  go_plots <- go_plots +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = chr_width)) +
    scale_fill_gradient(low = alpha(low_color, alpha = dot_alpha), high = alpha(high_color, alpha = dot_alpha))
  
  ######################################################
  # set the theme
  go_plots <- go_plots +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 1, size = font_size),
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
  go_plots <- go_plots +
    scale_size(range = dotsize)
  
  ######################################################
  # set the x-axis
  # browser()
  
  if (!is.null(xmin) && !is.null(xmax) && clustered == FALSE) {
    go_plots <- go_plots +
      scale_x_continuous(limits = c(xmin, xmax), n.breaks = x_breaks)
  }
  
  return(go_plots)
}

