###########################################
# draw_GO_GSEA_Terms.R
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
# usage: run_go_gsea_terms(degs, orgdb, key_type, pvalueCutoff, qvalueCutoff, ontology, pAdjustMethod, readable)


library(clusterProfiler)
library(enrichplot)


draw_gsea_terms <- function(
    gsea_results = NULL,
    geneSetID = NULL,
    title = "GO GSEA",
    base_size = 12,
    rel_heights = c(1.5, 0.5, 1),
    subplots = 1:3,
    color = NULL,
    pvalue_table = FALSE,
    ES_geom = "line" # line or dot
) {
  
  ######################################################
  # draw the GO GSEA terms
  # browser()
  
  gene_set_num <- length(geneSetID)
  color_list <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = color))(gene_set_num)
  
  gsea_plot <- gseaplot2(x = gsea_results,
                         geneSetID = geneSetID,
                         title = title,
                         color = color_list,
                         base_size = base_size,
                         rel_heights = rel_heights,
                         subplots = subplots,
                         pvalue_table = pvalue_table,
                         ES_geom = ES_geom
                         )

  return(gsea_plot)
}

