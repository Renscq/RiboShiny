####################################################
# function for DESeq2 analysis
# run_DESeq2.R
# 2024-02-28
# exprs: ribo counts dataframe
# design: design infomation dataframe
# padj: padj for filter the DEGs
# log2fc: log2fc for filter the DEGs
# anno: gene annotation dataframe
# join_flag: join the anno table by which column

require(tidyverse)
require(DESeq2)


run_deseq2 <- function(
    exprs = NULL,
    design = NULL,
    method = "Wald",
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
    column_to_rownames(var = "Gene") %>% 
    round()
  
  sample_info <- design %>%
    mutate(Condition = as.factor(Condition)) %>%
    mutate(SeqType = as.factor(SeqType))
  
  ####################################################
  ## 2. Create DESeq2 object
  # browser()
  ddsMat <- DESeqDataSetFromMatrix(
    countData = exprs,
    colData = sample_info,
    design =~ Condition)

  ####################################################
  ## 3. Run mRNA DESeq2

  ddsMat <- DESeq(ddsMat, test = method)
  res_rna <- results(ddsMat)
  if (shinkage) {
    res_rna <- lfcShrink(ddsMat, res = res_rna, coef = 2,  type = "apeglm")
  }
  
  res_norm <- res_rna %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'Gene') %>% 
    left_join(counts(ddsMat, normalized = TRUE) %>% 
                as.data.frame() %>% 
                rownames_to_column('Gene'), 
              by = 'Gene')
  
  ####################################################
  ## 4. annotate the gene sig table
  res_norm <- res_norm %>% 
    mutate(DEGs = if_else(padj < padj_cutoff & log2FoldChange >= log2fc_cutoff, 'UP',
                          if_else(padj < padj_cutoff & log2FoldChange <= -log2fc_cutoff, 'DOWN', 'NS')))

  ####################################################
  ## 5. annotate the gene message
  # browser()
  
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
  ## 7. output the merge results
  if (flt_sig) {
    deseq2_res <- list(
      DEGs = res_norm %>% filter(padj < padj_cutoff, abs(log2FoldChange) >= log2fc_cutoff),
      summary_table = summary_table)

  } else {
    deseq2_res <- list(
      DEGs = res_norm,
      summary_table = summary_table)
  }
  
  return(deseq2_res)
  
}