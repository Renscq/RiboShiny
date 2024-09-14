#################################
# fill the missing values with the mean of the column
# input: data frame
# output: data frame with missing values filled
# example: df <- fill_empty(df)

library(tidyverse)

fill_empty <- function(
    expr_df = NULL, 
    fill = "mean" # "none", mean", "median", or min value (> 0)
    ) {

  if (fill == "none"){
    return(expr_df)
    
  } else if (fill == "min") {
    # filter the min value
    fill_num <- sapply(expr_df, function(x) min(x[x > 0], na.rm = TRUE))
    
  } else if (fill == "mean") {
    # filter the min value
    fill_num <- colMeans(expr_df, na.rm = TRUE)
    
  } else if (fill == "median") {
    # filter the min value
    fill_num <- sapply(expr_df, function(x) median(x[x > 0], na.rm = TRUE))
  }
  
  # replace the min value with fill_num
  for (i in 1:ncol(expr_df)) {
    expr_df[expr_df[, i] == 0, i] <- fill_num[i] / 2
  }

  return(expr_df)
  
}