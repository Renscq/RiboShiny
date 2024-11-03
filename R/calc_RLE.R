library(tidyverse)

########################################################
# calculate the RLE
calculate_rle <- function(
    data = NULL) {
  medians <- apply(data, 1, median)
  rle_mat <- log2(data / medians)
  rle_mat <- rle_mat %>% 
    as.data.frame() %>% 
    dplyr::filter_all(all_vars(!is.na(.) & !is.infinite(.)))
  
  return(rle_mat)
}

########################################################
# make_rle_mat
make_rle_mat <- function(
    rle_mat = NULL, 
    phenodata = NULL) {

  rle_longer <- rle_mat %>% 
    tibble::rownames_to_column(var = "Gene") %>%
    tidyr::pivot_longer(cols = -Gene, names_to = "Sample", values_to = "RLE") %>%
    dplyr::left_join(phenodata, by = "Sample")
  
  return(rle_longer)
}