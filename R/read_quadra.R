#########################################
# read_quadra.R
# 
# This script reads in the quadra data and
# creates a data frame with the data.
# 
# data1_file: the file path of the first data
# data1_sheet: the sheet name of the first data
# data2_file: the file path of the second data
# data2_sheet: the sheet name of the second data
# category: the category of the data, e.g. RNA-Ribo, RNA-TE, Ribo-TE
# delta_column: the column name of the delta value, Buffered Forward exclusive intensified
# gene_column: the column name of the gene
# logfc_column: the column name of the log2FC
# padj_column: the column name of the adjusted p-value
# logfc_cutoff: the cutoff of log2FC
# pvalue_cutoff: the cutoff of adjusted p-value

require(tidyverse)
require(openxlsx)

anno_sections <- function(
    quadra_degs = NULL,
    x_fc = "RNA_log2FC",
    y_fc = "Ribo_log2FC",
    x_p = "RNA_Padj",
    y_p = "Ribo_Padj",
    sig1 = "RNA", 
    sig2 = "Ribo",
    fc_cutoff = 1,
    p_cutoff = 0.05
    
    ) {
  quadra_degs <- quadra_degs %>%
    dplyr::mutate(DEGs = if_else(!!sym(x_p) < p_cutoff & abs(!!sym(x_fc)) >= fc_cutoff & !!sym(y_p) < p_cutoff & abs(!!sym(y_fc)) >= fc_cutoff, "Both",
                         if_else(!!sym(x_p) < p_cutoff & abs(!!sym(x_fc)) >= fc_cutoff, sig1,
                         if_else(!!sym(y_p) < p_cutoff & abs(!!sym(y_fc)) >= fc_cutoff, sig2, "NS")))) %>%
    dplyr::mutate(Sections = if_else(!!sym(x_fc) <= -fc_cutoff & !!sym(y_fc) >= fc_cutoff, "S1",
                             if_else(!!sym(x_fc) > -fc_cutoff & !!sym(x_fc) < fc_cutoff & !!sym(y_fc) >= fc_cutoff, "S2", 
                             if_else(!!sym(x_fc) >= fc_cutoff & !!sym(y_fc) >= fc_cutoff, "S3",
                             if_else(!!sym(x_fc) <= -fc_cutoff & !!sym(y_fc) > -fc_cutoff & !!sym(y_fc) < fc_cutoff, "S4",
                             if_else(!!sym(x_fc) > -fc_cutoff & !!sym(x_fc) < fc_cutoff & !!sym(y_fc) > -fc_cutoff & !!sym(y_fc) < fc_cutoff, "S5",
                             if_else(!!sym(x_fc) >= fc_cutoff & !!sym(y_fc) > -fc_cutoff & !!sym(y_fc) < fc_cutoff, "S6",
                             if_else(!!sym(x_fc) <= -fc_cutoff & !!sym(y_fc) <= -fc_cutoff, "S7",
                             if_else(!!sym(x_fc) > -fc_cutoff & !!sym(x_fc) < fc_cutoff & !!sym(y_fc) <= -fc_cutoff, "S8",
                             if_else(!!sym(x_fc) >= fc_cutoff & !!sym(y_fc) <= -fc_cutoff, "S9", "Others"))))))))))
  
  return(quadra_degs)
  
}
  
  
read_quadra <- function(
    data1_file = NULL,
    data1_sheet = NULL,
    
    data2_file = NULL,
    data2_sheet = NULL,
    
    category = "RNA-Ribo",
    
    delta_column = "",
    gene_column = "Gene",
    
    label_column = "",
    logfc_column = "log2FC",
    padj_column = "padj",
    logfc_cutoff = 1,
    pvalue_cutoff = 0.05
    
    ) {
  

  ##############################################
  # import the raw data table
  # category <- unlist(strsplit(category, "-"))
  # browser()
  
  if (grepl("RNA-Ribo", category)) {
    tryCatch({
      rna_degs <- read.xlsx(data1_file, sheet = data1_sheet, colNames = T, rowNames = F)
      ribo_degs <- read.xlsx(data2_file, sheet = data2_sheet, colNames = T, rowNames = F)
    }, error = function(e) {
      message("Error: ", e)
      degs <- NULL
      return(degs)
    })
    
    if (label_column != "") {
      anno_degs <- rna_degs %>% 
        dplyr::select(!!sym(gene_column), !!sym(label_column))
    }
    if (delta_column != "" & delta_column != gene_column) {
      delta_degs <- rna_degs %>% 
        dplyr::select(!!sym(gene_column), !!sym(delta_column)) %>%
        dplyr::rename_all(~c("Gene", "Delta"))
    }
  }
  
  if (grepl("RNA-TE", category)) {
    tryCatch({
      rna_degs <- read.xlsx(data1_file, sheet = data1_sheet, colNames = T, rowNames = F)
      te_degs <- read.xlsx(data2_file, sheet = data2_sheet, colNames = T, rowNames = F)
    }, error = function(e) {
      message("Error: ", e)
      degs <- NULL
      return(degs)
    })
    
    if (label_column != "") {
      anno_degs <- rna_degs %>% 
        dplyr::select(!!sym(gene_column), !!sym(label_column))
    }
    if (delta_column != "" & delta_column != gene_column) {
      delta_degs <- rna_degs %>% 
        dplyr::select(!!sym(gene_column), !!sym(delta_column)) %>%
        dplyr::rename_all(~c("Gene", "Delta"))
    }
  }

  if (grepl("Ribo-TE", category)) {
    tryCatch({
      ribo_degs <- read.xlsx(data1_file, sheet = data1_sheet, colNames = T, rowNames = F)
      te_degs <- read.xlsx(data2_file, sheet = data2_sheet, colNames = T, rowNames = F)
    }, error = function(e) {
      message("Error: ", e)
      degs <- NULL
      return(degs)
    })
    
    if (label_column != "") {
      anno_degs <- ribo_degs %>% 
        dplyr::select(!!sym(gene_column), !!sym(label_column))
    }
    if (delta_column != "" & delta_column != gene_column) {
      delta_degs <- ribo_degs %>% 
        dplyr::select(!!sym(gene_column), !!sym(delta_column)) %>%
        dplyr::rename_all(~c("Gene", "Delta"))
    }
  }

  ##############################################
  # filter the raw data table
  if (grepl("RNA", category)) {
    rna_degs <- rna_degs %>% 
      dplyr::select(!!sym(gene_column), !!sym(logfc_column), !!sym(padj_column)) %>%
      dplyr::rename_all(~c("Gene", "RNA_log2FC", "RNA_Padj"))
  }
  
  if (grepl("Ribo", category)) {
    ribo_degs <- ribo_degs %>%
      dplyr::select(!!sym(gene_column), !!sym(logfc_column), !!sym(padj_column)) %>%
      dplyr::rename_all(~c("Gene", "Ribo_log2FC", "Ribo_Padj"))
  }
  
  if (grepl("TE", category)) {
    te_degs <- te_degs %>%
      dplyr::select(!!sym(gene_column), !!sym(logfc_column), !!sym(padj_column)) %>%
      dplyr::rename_all(~c("Gene", "TE_log2FC", "TE_Padj"))
  }

  ##############################################
  # merge the rna and ribo DEGs
  if (category == "RNA-Ribo") {
    quadra_degs <- rna_degs %>% 
      dplyr::inner_join(ribo_degs, by = "Gene")
    
    # merge the annotation data
    if (label_column != "") {
      quadra_degs <- quadra_degs %>% 
        dplyr::inner_join(anno_degs, by = "Gene")
    }
    
    # merge the delta column
    if (delta_column != "" & delta_column != gene_column) {
      quadra_degs <- quadra_degs %>% 
        dplyr::inner_join(delta_degs, by = "Gene")
    }
    
    quadra_degs <- anno_sections(quadra_degs, sig1 = "RNA", sig2 = "Ribo",
                                 fc_cutoff = logfc_cutoff, p_cutoff = pvalue_cutoff,
                                 x_fc = "RNA_log2FC", y_fc = "Ribo_log2FC", 
                                 x_p = "RNA_Padj", y_p = "Ribo_Padj")
  } else if (category == "RNA-TE") {
    quadra_degs <- rna_degs %>% 
      dplyr::inner_join(te_degs, by = "Gene")
    
    quadra_degs <- anno_sections(quadra_degs, sig1 = "RNA", sig2 = "TE",
                                 fc_cutoff = logfc_cutoff, p_cutoff = pvalue_cutoff,
                                 x_fc = "RNA_log2FC", y_fc = "TE_log2FC", 
                                 x_p = "RNA_Padj", y_p = "TE_Padj")
  } else if (category == "Ribo-TE") {
    quadra_degs <- ribo_degs %>% 
      dplyr::inner_join(te_degs, by = "Gene")
    
    quadra_degs <- anno_sections(quadra_degs, sig1 = "Ribo", sig2 = "TE",
                                 fc_cutoff = logfc_cutoff, p_cutoff = pvalue_cutoff,
                                 x_fc = "Ribo_log2FC", y_fc = "TE_log2FC", 
                                 x_p = "Ribo_Padj", y_p = "TE_Padj")
  } 
  
  # return the quadra data
  return(quadra_degs)
}