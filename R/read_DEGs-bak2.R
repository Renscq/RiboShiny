#########################################
# read_DEGs.R
# 
# This script reads in the quadra data and
# creates a data frame with the data.
# 
# data_file: the file path of the first data
# data_sheet: the sheet name of the first data
# category: the category of the data, e.g. RNA-Ribo, RNA-TE, Ribo-TE
# gene_column: the column name of the gene
# logfc_column: the column name of the log2FC
# padj_column: the column name of the adjusted p-value
# logfc_cutoff: the cutoff of log2FC
# pvalue_cutoff: the cutoff of adjusted p-value

require(tidyverse)
require(openxlsx)

anno_sections <- function(
    degs = NULL,
    logfc = "log2foldchange",
    pvalue = "padjust",
    fc_cutoff = 1,
    p_cutoff = 0.05
) {
  degs <- degs %>%
    dplyr::mutate(DEGs = if_else(!!sym(pvalue) < p_cutoff & abs(!!sym(logfc)) >= fc_cutoff, "UP",
                         if_else(!!sym(pvalue) < p_cutoff & abs(!!sym(logfc)) >= fc_cutoff, "DOWN", "NS")))
  
  return(degs)
}


read_degs <- function(
    data_file = NULL,
    data_sheet = NULL,
    gene_column = c("Gene", "baseMean", "log2FoldChange", "pvalue", "padj"), # or column names for edger: c("Gene", "LogFC", "logCPM", "PValue", "FDR")
    
    logfc_cutoff = 1,
    pvalue_cutoff = 0.05
    
) {
  
  ##############################################
  # import the raw data table
  # browser()
  # try to read the data file and catch the error return null
  tryCatch({
    degs <- read.xlsx(data_file, sheet = data_sheet, colNames = T, rowNames = F)
    
  }, error = function(e) {
    message("Error: ", e)
    degs <- NULL
    return(degs)
  })
  
  ##############################################
  # check the title
  if ('log2foldchange' %in% colnames(degs)) {
    gene_column <- c("Gene", "baseMean", "log2FoldChange", "pvalue", "padj")
    logfc_column <- "log2FoldChange"
    padj_column <- "padj"
    
  } else if ('FDR' %in% colnames(degs)) {
    gene_column <- c("Gene", "LogFC", "logCPM", "PValue", "FDR")
    logfc_column <- "LogFC"
    padj_column <- "FDR"
    
  } else {
    return(NULL)
  }
  
  ##############################################
  # filter the raw data table
  if (is.list(data_file)) {
    rna_name <- unlist(str_split(data_file, '\n|,'))
    file_num <- length(rna_name)
    
    for (i in 1:file_num) {
      degs <- read.xlsx(input$in_step15_merge_rna$datapath[i],
                        sheet = input$in_step15_merge_rna_sheet,
                        colNames = T, rowNames = F) %>%
        dplyr::select(selected_column) %>%
        dplyr::mutate(Groups = "rna_DEGs",
                      Files = rna_name[i])
      
      if (is.null(rna_degs)) {
        rna_degs <- degs
      } else {
        rna_degs <- rbind(rna_degs, degs)
      }
    }
    
    degs <- degs %>% 
      # select the columns contain the any one match gene_column
      dplyr::select(gene_column)
    
    degs <- anno_sections(degs,
                          fc_cutoff = logfc_cutoff, p_cutoff = pvalue_cutoff,
                          logfc = logfc_column, pvalue = padj_column)
    
  } else {
    degs <- degs %>% 
      # select the columns contain the any one match gene_column
      dplyr::select(gene_column)
    
    degs <- anno_sections(degs,
                          fc_cutoff = logfc_cutoff, p_cutoff = pvalue_cutoff,
                          logfc = logfc_column, pvalue = padj_column)
  }
  
  # return the degs data
  return(degs)
}