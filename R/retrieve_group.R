############################################################
# retrieve_group.R
# retrieve group message from uploaded file
# 
# exprs: data frame, expression matrix
# design: data frame, design matrix
# group1: character, group name
# group2: character, group name
# rowsum: numeric, row sum cutoff
# stdev: numeric, stdev cutoff

require(tidyverse)



retrieve_group <- function(
    exprs = NULL,
    design_table = NULL,
    sample_flt = "Sample",
    group_flt = NULL,
    group1 = NULL,
    group2 = NULL,
    rowsum = 0,
    stdev = NA
    ) {
  
  # select the design
  group1_samples <- design_table %>%
    dplyr::filter(!!sym(group_flt) == group1)
  
  group2_samples <- design_table %>%
    dplyr::filter(!!sym(group_flt) == group2)
  
  # select the samples
  group1_exprs <- exprs %>%
    dplyr::select(any_of(group1_samples[[sample_flt]]))
  group2_exprs <- exprs %>%
    dplyr::select(any_of(group2_samples[[sample_flt]]))
  
  # filter the gene average expression
  if (rowsum > 0) {
    group1_exprs <- group1_exprs %>% 
      dplyr::filter(rowMeans(.) >= rowsum)
    group2_exprs <- group2_exprs %>%
      dplyr::filter(rowMeans(.) >= rowsum)
  }
  
  # filter the variance
  if (!is.na(stdev)) {
    group1_genes <- log2(group1_exprs + 1) %>% 
      dplyr::filter(apply(., 1, sd) <= stdev) %>% 
      rownames()
    group1_exprs <- group1_exprs %>% 
      dplyr::filter(rownames(.) %in% group1_genes)
    
    group2_genes <- log2(group2_exprs + 1) %>% 
      dplyr::filter(apply(., 1, sd) <= stdev) %>% 
      rownames()
    group2_exprs <- group2_exprs %>%
      dplyr::filter(rownames(.) %in% group2_genes)
  }
  
  # join the two groups
  gene_exprs <- group1_exprs %>% 
    tibble::rownames_to_column("Gene") %>%
    dplyr::inner_join(group2_exprs %>% tibble::rownames_to_column("Gene"), 
                      by = "Gene") %>% 
    dplyr::mutate(across(where(is.numeric), as.double))
  
  # make the flt design
  flt_design <- rbind(group1_samples, group2_samples) %>% 
    dplyr::group_by(!!sym(group_flt)) %>%
    dplyr::mutate(Condition = cur_group_id())
  
  return(list(gene_exprs = gene_exprs,
              flt_design = flt_design))
  
}