# function to draw reads digestion

library(ggplot2)
library(ggseqlogo)
library(RColorBrewer)
library(tidyverse)

draw_seqlogo <- function(
    pwm = NULL, 
    xstart = NULL,
    xend = NULL,
    x_breaks = 1,
    
    facet = "wrap",
    scales = "free_x", 
    ncol = NULL,
    nrow = NULL,
    method = "probability", # bits, information, probability, custom
    col_scheme = "auto", # auto chemistry chemistry2 hydrophobicity nucleotide nucleotide2 base_pairing clustalx taylor
    col_alpha = 1,
    seq_type = "dna", # dna, rna, aa
    stack_width = 0.8,
    font_family = "helvetica_regular",
    # helvetica_regular helvetica_bold helvetica_light roboto_medium roboto_bold roboto_regular akrobat_bold 
    # akrobat_regular roboto_slab_bold roboto_slab_regular roboto_slab_light xkcd_regular
    font_size = 10,
    
    xlabel = NULL,
    ylabel = NULL,
    title = NULL
    ) {

  
  # browser()
  
  # chekc the breaks
  if (x_breaks < 1) {
    x_breaks = 1
  } else if (x_breaks > (xend - xstart)) {
    x_breaks = round((xend - xstart) / 4)
  }
  
  # draw the seqlogo
  sequence_logo <- ggseqlogo(data = pwm, 
                              method = method,
                              facet = facet,
                              scales = scales,
                              ncol = ncol, nrow = nrow,
                              seq_type = seq_type,
                              col_scheme = col_scheme,
                              alpha = col_alpha,
                              stack_width = stack_width,
                              font = font_family) +
    scale_x_continuous(# limits = c(xstart, xend),
                       breaks = seq(1, (xend - xstart + 1), x_breaks), 
                       # breaks = seq(xstart, xend, 2),
                       labels = seq(xstart, xend, x_breaks)) +
    # annotate("rect", xmin = 5.5, xmax = 6.5, 
    #          ymin = -end5_limits, ymax = end5_limits, 
    #          alpha = 0.1, col = "black", fill = "yellow") +
    xlab('Site') +
    ylab(method) +
    theme_bw() +
    theme(strip.background = element_blank(),
          panel.grid = element_blank(),
          strip.text = element_text(size = font_size - 1),
          axis.text = element_text(size = font_size),
          axis.title = element_text(size = font_size + 1),
          plot.title = element_text(size = font_size + 2, vjust = 0.5, hjust = 0.5))

  if (!is.null(title)) {
    sequence_logo <- sequence_logo + 
      ggtitle(title)
  }
  
  if (!is.null(xlabel)) {
    sequence_logo <- sequence_logo + 
      xlab(xlabel)
  }
  
  if (!is.null(ylabel)) {
    sequence_logo <- sequence_logo + 
      ylab(ylabel)
  }
  
  # ggplot2::ggsave(plot = sequence_logo, filename = './output/reads_digest_seqlogo.png', width = 9, height = 4, dpi = 300)
  
  return(sequence_logo)
  
}




