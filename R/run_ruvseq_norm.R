# Normalize count data using RUVseq
library(RUVSeq)
library(EDASeq)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)


########################################################
# normalize the count data using RUVSeq
run_ruvseq_norm <- function(
    count = NULL,
    phenodata = NULL,
    group = NULL,
    method = "full", # upper, full
    fill_color = "Blues",
    fill_alpha = 0.9,
    profile = "RNA" # RNA, RIBO
) {
  
  # browser()

  ######################################################
  ## make the expression table
  
  # browser()
  count_set <- newSeqExpressionSet(as.matrix(count) %>% round(), phenodata = phenodata)
  # phenodata <- pData(count_set)
  raw_rle_mat <- calculate_rle(count_set@assayData$counts)

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

  ######################################################
  # format the output file
  count_set2_norm <- count_set2@assayData$normalizedCounts %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Gene")
  
  count_set2_rpm <- count_set2@assayData$normalizedCounts %>%
    as.data.frame() %>%
    dplyr::mutate_all(function(x) x / sum(x) * 1e6) %>%
    tibble::rownames_to_column(var = "Gene")

  ######################################################
  return(list(count_set2_norm = count_set2_norm,
              count_set2_rpm = count_set2_rpm,
              phenodata = pData(count_set)))
  
}






