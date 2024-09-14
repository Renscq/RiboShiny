# function to draw the gene schematic


draw_isoforms_schematic <- function(
    gene_name = NULL,
    
    utr5_length = NULL,
    cds_length = NULL,
    utr3_length = NULL,
    
    gene_color = c("#75aadb", "#264abd", "#75aadb"),
    gene_width = c(2, 5, 2),
    fill_alpha = 1,
    
    xmin = NA,
    xmax = NA,
    xbreaks = 5,
    
    arrow_length = 0.1,
    arrow_width = 0.8,
    arrow_color = "white",
    
    font_size = 12) {

  # generate the gene region table
  gene_length <- utr5_length + cds_length + utr3_length
  
  isoforms_region <- data.frame(Name = gene_name,
                                Region = c("UTR5", "CDS", "UTR3"),
                                Start = c(0, utr5_length, utr5_length + cds_length),
                                End = c(utr5_length, utr5_length + cds_length, gene_length)) %>% 
    dplyr::mutate(Region = factor(Region, levels = c("UTR5", "CDS", "UTR3")))
  
  # generate the arrow position by arrow_num
  start_pos <- c(0.2, 0.5, 0.8)
  arrow_start <- cds_length * start_pos + utr5_length
  arrow_end <- arrow_start + 1
  
  # draw the gene schematic
  isoforms_plot <- ggplot(isoforms_region, 
                          aes(x = Start, xend = End, y = Name, yend = Name, color = Region)) +
    geom_segment(linewidth = c("UTR5" = gene_width[1], "CDS" = gene_width[2], "UTR3" = gene_width[3]),
                 alpha = fill_alpha) +
    scale_color_manual(values = c("UTR5" = gene_color[1], "CDS" = gene_color[2], "UTR3" = gene_color[3])) +
    geom_segment(aes(x = arrow_start, xend = arrow_end, y = Name, yend = Name),
                 arrow = arrow(type = "open", ends = "last", length = unit(arrow_length, "inches")), 
                 color = arrow_color, linewidth = arrow_width) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(size = font_size, angle = 0, vjust = 0.5, hjust = 0.5),
          axis.text.y = element_text(size = font_size, angle = 90, vjust = 1, hjust = 0.5),
          legend.position = "none") +
    labs(x = "", y = "Gene")
  
  # reset the x-axis limits
  if (!is.na(xmin) & !is.na(xmax)) {
    isoforms_plot <- isoforms_plot +
      coord_cartesian(xlim = c(xmin, xmax), clip = "on")
      # scale_x_continuous(limits = c(xmin, xmax), n.breaks = xbreaks)
  }
  
  
  return(isoforms_plot)
}

