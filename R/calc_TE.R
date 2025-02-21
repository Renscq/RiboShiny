###############################################
# calc_TE.R
# 
# phenodata: data frame contain the sample names, seqtype, group, and rank
# example like this:
# Sample  Seqtype Group   Rank
# m_SRR1    RNA     A       1
# m_SRR2    RNA     A       2
# m_SRR3    RNA     B       3
# m_SRR4    RNA     B       4
# R_SRR1    RIBO    A       1
# R_SRR2    RIBO    A       2
# R_SRR3    RIBO    B       3
# R_SRR4    RIBO    B       4
# T_SRR1    TE      A       1
# T_SRR2    TE      A       2
# T_SRR3    TE      B       3
# T_SRR4    TE      B       4


calc_te <- function(
    rna = NULL, # RNA-seq data frame contain the gene rownames
    ribo = NULL, # Ribo-seq data frame contain the gene rownames
    phenodata = NULL, # phenodata data frame contain the sample names and group
    linked = FALSE, # merge the RNA-seq and Ribo-seq data frame
    clean = TRUE
    ) {
  
  # browser()
  ##############################################
  # check the overlap sample of RNA-seq and Ribo-seq
  phenodata <- phenodata %>% 
    dplyr::arrange(SeqType, Rank)
  
  rna_sample <- phenodata %>% 
    dplyr::filter(Sample %in% colnames(rna))
  ribo_sample <- phenodata %>% 
    dplyr::filter(Sample %in% colnames(ribo))
  
  both_sample <- intersect(rna_sample$Rank, ribo_sample$Rank)
  
  phenodata <- phenodata %>% 
    dplyr::filter(Rank %in% both_sample)
  
  # filter the sample of RNA-seq and Ribo-seq
  rna_flt <- rna %>% 
    dplyr::select(any_of(phenodata$Sample))
  
  ribo_flt <- ribo %>%
    dplyr::select(any_of(phenodata$Sample))
  
  
  ##############################################
  # check the overlap gene of RNA-seq and Ribo-seq
  # browser()
  
  rna_gene <- rownames(rna_flt)
  ribo_gene <- rownames(ribo_flt)
  
  both_gene <- intersect(rna_gene, ribo_gene)
  
  # filter the gene of RNA-seq and Ribo-seq
  rna_flt <- rna_flt %>% 
    dplyr::filter(rownames(.) %in% both_gene) %>% 
    dplyr::arrange(rownames(.))
  
  ribo_flt <- ribo_flt %>%
    dplyr::filter(rownames(.) %in% both_gene) %>% 
    dplyr::arrange(rownames(.))
  
  ##############################################
  # calculate the TE
  # browser()
  
  te <- ribo_flt / rna_flt
  
  # clean the data
  if (clean) {
    # remove the NA and inf
    te <- te %>% 
      dplyr::filter(if_all(where(is.numeric), ~ !is.na(.))) %>% 
      dplyr::filter(if_all(where(is.numeric), ~ !is.infinite(.)))
  }
  
  # rename the colnames
  te_sample <- phenodata %>% 
    dplyr::filter(SeqType == "TE")
  
  te <- te %>% 
    dplyr::rename_all(~te_sample$Sample) %>% 
    tibble::rownames_to_column(var = "Gene") %>% 
    as.data.frame()
  
  ##############################################
  # merge all data
  if (isTRUE(linked)) {
    rna_flt <- rna_flt %>% tibble::rownames_to_column(var = "Gene")
    ribo_flt <- ribo_flt %>% tibble::rownames_to_column(var = "Gene")

    te <- te %>% 
      dplyr::left_join(rna_flt, by = "Gene") %>% 
      dplyr::left_join(ribo_flt, by = "Gene") %>% 
      as.data.frame()
    
  }
  
  return(list(te = te, phenodata = phenodata))
  
}