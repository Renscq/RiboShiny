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
    logfc = "RNA_log2FC",
    pvalue = "RNA_Padj",
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
    category = "RNA", # RNA, Ribo, TE
    gene_column = "Gene",

    logfc_column = "log2FC",
    padj_column = "padj",
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
  # filter the raw data table
  if (category == "RNA") {
    degs <- degs %>% 
      dplyr::select(!!sym(gene_column), !!sym(logfc_column), !!sym(padj_column)) %>%
      dplyr::rename_all(~c("Gene", "RNA_log2FC", "RNA_Padj"))
    
    degs <- anno_sections(degs,
                          fc_cutoff = logfc_cutoff, p_cutoff = pvalue_cutoff,
                          logfc = "RNA_log2FC", pvalue = "RNA_Padj")
    
    degs <- degs %>% dplyr::mutate(Type = "RNA")
    
  } else if (category == "Ribo") {
    degs <- degs %>%
      dplyr::select(!!sym(gene_column), !!sym(logfc_column), !!sym(padj_column)) %>%
      dplyr::rename_all(~c("Gene", "Ribo_log2FC", "Ribo_Padj"))
    
    degs <- anno_sections(degs,
                          fc_cutoff = logfc_cutoff, p_cutoff = pvalue_cutoff,
                          logfc = "Ribo_log2FC", pvalue = "Ribo_Padj")
    
    degs <- degs %>% dplyr::mutate(Type = "Ribo")
    
  } else if (category == "TE") {
    degs <- degs %>%
      dplyr::select(!!sym(gene_column), !!sym(logfc_column), !!sym(padj_column)) %>%
      dplyr::rename_all(~c("Gene", "TE_log2FC", "TE_Padj"))
    
    degs <- anno_sections(degs,
                          fc_cutoff = logfc_cutoff, p_cutoff = pvalue_cutoff,
                          logfc = "TE_log2FC", pvalue = "TE_Padj")
    
    degs <- degs %>% dplyr::mutate(Type = "TE")
  }
  
  # return the degs data
  return(degs)
}