####################################################
# function for edgeR analysis
# run_edgeR.R
# 2024-02-28
# exprs: ribo counts dataframe
# design: design infomation dataframe
# padj: padj for filter the DEGs
# log2fc: log2fc for filter the DEGs
# anno: gene annotation dataframe
# join_flag: join the anno table by which column


require(tidyverse)
require(edgeR)


run_edger <- function(
    exprs = NULL,
    design = NULL,
    method = "LRT", # LRT or glm, exact
    bcv = 0.01, # bcv for exact test
    padj_cutoff = 0.05, # padj for filter the DEGs
    log2fc_cutoff = 1, # log2fc for filter the DEGs
    flt_sig = F, # filter the sig diff genes
    shinkage = F, # shinkage the log2fc
    anno = NULL,
    join_flag = "Gene"
) {
  
  # browser()
  
  ####################################################
  ## 1. check the expression table and design table
  exprs <- exprs %>% 
    dplyr::rename(Gene = colnames(.)[1]) %>% 
    column_to_rownames(var = "Gene")
  
  sample_info <- design %>%
    mutate(Condition = as.factor(Condition)) %>%
    mutate(SeqType = as.factor(SeqType))
  
  # browser()
  ####################################################
  ## 4. run the edgeR
  if (method == "LRT") {
    group <- sample_info$Condition
    design <- model.matrix(~group)
    y <- DGEList(counts = exprs, group = group)
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTrendedDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef = 2)
    
    # merge results
    fit_values <- lrt$fitted.values %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'Gene')
    
    res_norm <- topTags(lrt, 
                         n = nrow(exprs),
                         adjust.method = "BH", 
                         sort.by = "PValue", 
                         p.value = 1) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'Gene') %>% 
      left_join(fit_values, by = 'Gene')
    
  } else if (method == "QLF") {
    group <- sample_info$Condition
    design <- model.matrix(~group)
    y <- DGEList(counts = exprs, group = group)
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTrendedDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    qlfit <- glmQLFit(y, design)
    qlf <- glmQLFTest(qlfit, coef = 2)
    
    # merge results
    fit_values <- qlf$fitted.values %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'Gene')
    
    res_norm <- topTags(qlf, 
                         n = nrow(exprs),
                         adjust.method = "BH", 
                         sort.by = "PValue", 
                         p.value = 1) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'Gene') %>% 
      left_join(fit_values, by = 'Gene')
    
  } else if (method == "exact") {
    y <- DGEList(counts = exprs, group = gl(2,1))
    y <- calcNormFactors(y)
    bcv <- bcv
    et <- exactTest(y, dispersion = bcv^2)
    
    # merge results
    fit_values <- et$fitted.values %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'Gene')
    
    res_norm <- topTags(et, 
                         n = nrow(exprs),
                         adjust.method = "BH", 
                         sort.by = "PValue", 
                         p.value = 1) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = 'Gene') %>% 
      left_join(fit_values, by = 'Gene')
    
  }
  
  ####################################################
  ## 5. annotate the gene sig table
  res_norm <- res_norm %>% 
    mutate(DEGs = if_else(FDR < padj_cutoff & logFC >= log2fc_cutoff, 'UP',
                          if_else(FDR < padj_cutoff & logFC <= -log2fc_cutoff, 'DOWN', 'NS')))
  
  ####################################################
  ## 5. annotate the gene message
  if (!is.null(anno)) {
    res_norm <- res_norm %>% 
      left_join(anno, by = c("Gene" = join_flag))
  }
  
  ####################################################
  ## 6. summary the table
  
  # browser()
  
  summary_table <- rbind(
    table(res_norm$DEGs) %>% 
      as.data.frame() %>% 
      rename_all(~c('Class', 'Count')) %>% 
      mutate(Condition = paste("log2FC:", log2fc_cutoff, "padj:", padj_cutoff))
  )
  
  ####################################################
  ## 11. output the merge results
  if (flt_sig) {
    edger_res <- list(
      DEGs = res_norm %>% filter(FDR < padj_cutoff, abs(logFC) >= log2fc_cutoff),
      summary_table = summary_table)
    
  } else {
    edger_res <- list(
      DEGs = res_norm,
      summary_table = summary_table)
  }
  
  return(edger_res)
  
}