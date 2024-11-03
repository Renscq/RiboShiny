################################################
# draw MA plot
# 
# degs: the DEGs data
# log2fc_column: the log2FC column name
# pvalue_column: the pvalue column name
# log2fc: the log2FC cutoff
# pvalue: the pvalue cutoff

# class: the DEGs class
# up_color: the up DEGs color
# down_color: the down DEGs color
# ns_color: the non-significant DEGs color
# dot_size: the dot size
# font_size: the font size

# x_max: the x-axis max value
# x_min: the x-axis min value

# fill_0: fill the 0 Padj to 0.1 * min(Padj)

# vol_sqrt: scale the y-axis to sqrt
# vol_log2: scale the y-axis to log2

# vol_title: the title of the volcano plot
# vol_xlab: the x-axis label
# vol_ylab: the y-axis label

# label_col: the gene name column
# label_color: the gene name color
# label_size: the gene name size
# label_list: the gene name list
# label_num: the gene name number

require(ggplot2)
require(tidyverse)
require(ggrepel)
require(RColorBrewer)


draw_MAplot <- function(
    degs = NULL,
    
    log2fc_column = "log2FoldChange",
    basemean_column = "baseMean",
    pvalue_column = "padj",
    log2fc = 1, 
    pvalue = 0.05,
    
    class = "DEGs",
    up_color = "red",
    down_color = "blue",
    ns_color = "grey",
    
    dot_size = 3,
    font_size = 12,
    
    x_max = NA,
    x_min = NA,
    
    fill_0 = FALSE,
    
    vol_sqrt = FALSE,
    vol_log2 = FALSE,
    
    vol_title = "Volcano Plot",
    vol_xlab = "log2FC",
    vol_ylab = "-log10(Padj)",
    
    label_type = "text", # text, dot
    
    label_col = "Gene",
    label_color = "black",
    label_size = 3,
    label_list = NULL,
    label_num = 0,
    overlaps_num = 10
    
    ) {
  
  ######################################################
  ## filter the sig gene
  degs_data <- degs %>% 
    dplyr::select(!!sym(label_col), 
                  !!sym(log2fc_column), 
                  !!sym(basemean_column), 
                  !!sym(pvalue_column), 
                  !!sym(class)) %>% 
    rename_with(~ c('Gene', 'log2FC', 'baseMean', 'Padj', 'DEGs'), everything()) %>% 
    arrange(Padj, desc(abs(log2FC)), desc(baseMean))
  
  ######################################################
  ## drop NA and INF
  degs_data <- degs_data %>% 
    tidyr::drop_na(Padj) %>%
    tidyr::drop_na(log2FC) %>%
    dplyr::filter(if_any(everything(), ~ !is.infinite(.)))
  
  # browser()
  ######################################################
  ## set the xlim and ylim
  if (is.na(x_max)) {
    x_max <- max(degs_data$log2FC)
  }
  
  if (is.na(x_min)) {
    x_min <- min(degs_data$log2FC)
  }
  
  degs_data <- degs_data %>% 
    dplyr::filter(log2FC >= x_min, log2FC <= x_max)
  
  ######################################################
  ## convert Padj = 0 to Padj * 0.1
  # browser()
  if (fill_0) {
    ylim_min <- degs_data %>% dplyr::select(Padj) %>% dplyr::filter(. != 0) %>% min()
    degs_data <- degs_data %>%
      mutate(Padj = ifelse(Padj == 0, ylim_min / 10, Padj))
  }
  
  ######################################################
  # summary the DEGs
  # browser()
  degs_num <- degs_data %>% 
    dplyr::reframe(UP = sum(DEGs == "UP"),
                   DOWN = sum(DEGs == "DOWN"),
                   NS = sum(DEGs == "NS"))
  
  subtitle <- paste0("DOWN: ", degs_num$DOWN, ", NS: ", degs_num$NS, ", UP: ", degs_num$UP)
  
  ######################################################
  # set the label type
  if (label_type == "text") {
    ######################################################
    ## draw the volcano
    ma_plot <- ggplot(degs_data, aes(x = baseMean, y = log2FC)) +
      geom_point(aes(color = DEGs), size = dot_size, alpha = 1) +
      scale_color_manual(limits = c("DOWN", "NS", "UP"),
                         values = c(down_color, ns_color, up_color)) +
      theme_bw() + 
      theme(panel.grid = element_blank(),
            plot.title = element_text(size = font_size + 2, vjust = 0.5, hjust = 0.5, color = 'black'),
            plot.subtitle = element_text(size = font_size - 2, vjust = 0.5, hjust = 0.5, color = 'black'),
            axis.title = element_text(size = font_size + 1, color = "black"),
            axis.text = element_text(size = font_size, color = 'black'),
            panel.background = element_rect(color = 'black', fill = 'transparent'),
            legend.title = element_blank(),
            legend.key = element_rect(fill = 'transparent'),
            legend.text = element_text(size = font_size - 2),
            legend.background = element_rect(fill = 'transparent')) +
      # plot the line to mark log2FC and padj
      geom_hline(yintercept = c(-log2fc, log2fc), color = 'black', linewidth = 0.4, linetype = 2) +
      # geom_vline(xintercept = -log(pvalue, 10), color = 'black', linewidth = 0.4, linetype = 2) +
      scale_y_continuous(limits = c(x_min, x_max)) +
      ggtitle(vol_title, subtitle = subtitle) +
      labs(x = vol_xlab, y = vol_ylab)
    
    ######################################################
    ## scale x/y sqrt
    # browser()
    if (vol_sqrt) {
      ma_plot <- ma_plot + 
        scale_x_sqrt()
    } else if (vol_log2) {
      ma_plot <- ma_plot + 
        scale_x_log10()
    }
    
    ######################################################
    ## filter the gene label
    label_num_gene <- data.frame()
    label_list_gene <- data.frame()
    
    if (label_num > 0) {
      ## filter the gene name
      label_num_gene <- degs_data %>%
        head(label_num)
    }
    
    if (length(label_list) > 0) {
      ## retrieve the gene number
      
      label_list_gene <- degs_data %>%
        dplyr::filter(Gene %in% label_list) 
    }
    
    label_col_gene <- rbind(label_num_gene, label_list_gene)
    
    if (nrow(label_col_gene) > 0) {
      ## add the gene name
      ma_plot <- ma_plot +
        geom_text_repel(data = label_col_gene %>% dplyr::distinct(),
                        aes(x = baseMean, y = log2FC, label = Gene),
                        size = label_size, color = label_color, fontface = "bold",
                        max.overlaps = overlaps_num)
    }

    
  } else if (label_type == "dot") {
    ######################################################
    ## draw the volcano
    ma_plot <- ggplot(degs_data, aes(x = baseMean, y = log2FC)) +
      geom_point(color = 'grey80', size = dot_size, alpha = 1) +
      theme_bw() + 
      theme(panel.grid = element_blank(),
            plot.title = element_text(size = font_size + 2, vjust = 0.5, hjust = 0.5, color = 'black'),
            plot.subtitle = element_text(size = font_size - 2, vjust = 0.5, hjust = 0.5, color = 'black'),
            axis.title = element_text(size = font_size + 1, color = "black"),
            axis.text = element_text(size = font_size, color = 'black'),
            panel.background = element_rect(color = 'black', fill = 'transparent'),
            legend.title = element_blank(),
            legend.key = element_rect(fill = 'transparent'),
            legend.text = element_text(size = font_size - 2),
            legend.background = element_rect(fill = 'transparent')) +
      # plot the line to mark log2FC and padj
      geom_hline(yintercept = c(-log2fc, log2fc), color = 'black', linewidth = 0.4, linetype = 2) +
      # geom_vline(xintercept = -log(pvalue, 10), color = 'black', linewidth = 0.4, linetype = 2) +
      scale_y_continuous(limits = c(x_min, x_max)) +
      ggtitle(vol_title, subtitle = subtitle) +
      labs(x = vol_xlab, y = vol_ylab)
    
    ######################################################
    ## scale x/y sqrt
    # browser()
    if (vol_sqrt) {
      ma_plot <- ma_plot + 
        scale_x_sqrt()
    } else if (vol_log2) {
      ma_plot <- ma_plot + 
        scale_x_log10()
    }
    
    ######################################################
    ## filter the gene label
    label_num_gene <- data.frame()
    label_list_gene <- data.frame()
    
    if (label_num > 0) {
      ## filter the gene name
      label_num_gene <- degs_data %>%
        head(label_num)
    }
    
    if (length(label_list) > 0) {
      ## retrieve the gene number
      
      label_list_gene <- degs_data %>%
        dplyr::filter(Gene %in% label_list) 
    }
    
    label_col_gene <- rbind(label_num_gene, label_list_gene)
      
    if (nrow(label_col_gene) > 0) {
    ## add the gene name
    ma_plot <- ma_plot +
      geom_point(data = label_col_gene %>% dplyr::distinct(),
                 aes(x = baseMean, y = log2FC), 
                 color = label_color, size = dot_size, alpha = 1)
    }
  }
  
  return(ma_plot)
  
}