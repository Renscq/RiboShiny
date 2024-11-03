# Normalize count data using RUVseq
library(RUVSeq)
library(EDASeq)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)


source("R/draw_RLE.R")
source("R/calc_RLE.R")
source("R/draw_PCA.R")

########################################################
# normalize the count data using RUVSeq
ruvseq_norm <- function(
    count = NULL,
    sample_design = NULL,
    group = NULL,
    method = "full", # upper, full
    fill_color = "Blues",
    fill_alpha = 0.9,
    profile = "RNA" # RNA, RIBO
) {
  
  # browser()
  ######################################################
  # get the sample information
  phenoData <- sample_design %>%
    column_to_rownames(var = "Sample") %>% 
    dplyr::filter(SeqType == profile) %>%
    dplyr::mutate(Group = factor(!!sym(group), levels = unique(!!sym(group)))) %>% 
    dplyr::mutate(Batch = Time)
  
  ######################################################
  ## make the expression table
  
  # browser()
  count_set <- newSeqExpressionSet(as.matrix(count) %>% round(), phenoData = phenoData)
  # phenoData <- pData(count_set)
  raw_rle_mat <- calculate_rle(count_set@assayData$counts)
  raw_rle_longer <- make_rle_mat(raw_rle_mat, phenoData)
  raw_rle_plot <- draw_rle(rle_mat = raw_rle_longer, x = "Sample", y = "RLE", title = paste(profile, "(Raw count)"),
                           ymin = -4, ymax = 4, group = group,
                           fill_color = fill_color, fill_alpha = fill_alpha)
  raw_pca_res <- draw_pca(in_data = count_set@assayData$counts, 
                          phenodata = phenoData %>% tibble::rownames_to_column(var = "Sample"), 
                          group = group, title = paste(profile, "(Raw count)"),
                          fill_color = fill_color, fill_alpha = fill_alpha)
  # browser()
  ######################################################
  # normalize the count data
  if (method == "upper") {
    ## normalize by upper
    count_set2 <- betweenLaneNormalization(count_set, which = "upper")

  } else if (method == "full") {
    ## normalize by full
    count_set2 <- betweenLaneNormalization(count_set, which = "full")
  }
  
  norm_rle_mat <- calculate_rle(count_set2@assayData$normalizedCounts)
  norm_rle_longer <- make_rle_mat(norm_rle_mat, phenoData)
  norm_rle_plot <- draw_rle(rle_mat = norm_rle_longer, x = "Sample", y = "RLE", title = paste(profile, "(Normalized)"),
                            ymin = -4, ymax = 4, group = group,
                            fill_color = fill_color, fill_alpha = fill_alpha)
  norm_pca_res <- draw_pca(in_data = count_set2@assayData$normalizedCounts, 
                           phenodata = phenoData %>% tibble::rownames_to_column(var = "Sample"), 
                           group = group, title = paste(profile, "(Normalized)"),
                           fill_color = fill_color, fill_alpha = fill_alpha)
  
  ######################################################
  # format the output file
  count_set2_norm <- count_set2@assayData$normalizedCounts %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Gene")
  
  count_set2_rpm <- count_set2@assayData$normalizedCounts %>%
    as.data.frame() %>%
    dplyr::mutate_all(function(x) x / sum(x) * 1e6) %>%
    tibble::rownames_to_column(var = "Gene")
  
  # format the figure
  rle_plot <- cowplot::plot_grid(raw_rle_plot, norm_rle_plot, ncol = 1)
  pca_plot <- cowplot::plot_grid(raw_pca_res$pca_scatter, norm_pca_res$pca_scatter, ncol = 2)
  
  ######################################################
  return(list(
    count_set2_norm = count_set2_norm,
    count_set2_rpm = count_set2_rpm,
    rle_plot = rle_plot,
    pca_plot = pca_plot))
  
}






