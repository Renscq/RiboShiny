# function convert pfm to ppm and pwm


calc_pfm_ppm_pwm <- function(
    seq_pfm = NULL,
    seq_freq = 0.25,
    method = "custom"
    ) {

  # convert pfm to ppm and pwm
  seq_ppm <- seq_pfm / colSums(seq_pfm)
  seq_pwm <- log2(seq_ppm / seq_freq)
  
  # browser()
  
  if (method == "foldchange") {
    seq_mat <- log2(seq_pfm / rowMeans(seq_pfm))
  
  } else if (method == "delta") {
    seq_mat <- seq_pfm - rowMeans(seq_pfm)
  
  } else if (method == "enrichment") {
    col_mean <- mean(colMeans(seq_pfm))
    
    seq_mat <- (seq_pfm - rowMeans(seq_pfm)) / mean(colMeans(seq_pfm))
    
  } else {
    seq_mat <- seq_pfm
  }
  
  return(list(seq_mat = seq_mat, seq_ppm = seq_ppm, seq_pwm = seq_pwm))
  
}




