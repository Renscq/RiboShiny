####################################################
# function to calculate the TE
# run_deltaTE.R
# 2023-11-20
# dt_ribo: ribo counts dataframe
# dt_mrna: mrna counts dataframe
# dt_design: design infomation dataframe
# padj: padj for filter the DEGs
# log2fc: log2fc for filter the DEGs
# anno: gene annotation dataframe
# join_flag: join the anno table by which column


require(DESeq2)


run_delta_te <- function(
    dt_ribo = NULL,
    dt_mrna = NULL,
    dt_design = NULL,
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
  ribo_counts <- dt_ribo %>% 
    dplyr::rename(Gene = colnames(.)[1])
  rna_counts <- dt_mrna %>% 
    dplyr::rename(Gene = colnames(.)[1])
  
  both_gene <- intersect(ribo_counts$Gene, rna_counts$Gene)
  
  ribo_counts <- ribo_counts %>% 
    dplyr::filter(Gene %in% both_gene) %>% 
    dplyr::arrange(Gene) %>% 
    tibble::column_to_rownames('Gene') %>% 
    round()
  
  rna_counts <- rna_counts %>% 
    dplyr::filter(Gene %in% both_gene) %>% 
    dplyr::arrange(Gene) %>% 
    tibble::column_to_rownames('Gene') %>%
    round()
  
  sample_info <- dt_design %>%
    mutate(Condition = as.factor(Condition)) %>%
    mutate(SeqType = as.factor(SeqType))
  
  # sample_info <- dt_design %>% 
  #   mutate(Condition = ave(as.character(Group), SeqType, FUN = function(x) as.integer(factor(x)))) %>% 
  #   mutate(Condition = as.factor(Condition)) %>% 
  #   mutate(SeqType = as.factor(SeqType))
  
  # browser()
  
  ####################################################
  ## 2. Create DESeq2 object
  ddsMat <- DESeqDataSetFromMatrix(
    countData = cbind(ribo_counts, rna_counts),
    colData = sample_info,
    design =~ Condition + SeqType + Condition:SeqType)
  
  ####################################################
  ## 3. Run all DESeq2
  ddsMat <- DESeq(ddsMat)
  resultsNames(ddsMat)
  
  # res <- results(ddsMat, name = 'Condition2.SeqTypeRNA') 
  res <- results(ddsMat) 
  res_norm <- res %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'Gene') %>% 
    left_join(counts(ddsMat, normalized = TRUE) %>% 
                as.data.frame() %>% 
                rownames_to_column('Gene'), by = 'Gene')
  
  ####################################################
  ## 4. Run mRNA DESeq2
  ddsMat_rna <- DESeqDataSetFromMatrix(
    countData = rna_counts,
    colData = sample_info[which(sample_info$SeqType == "RNA"),],
    design =~ Condition)
  
  ddsMat_rna <- DESeq(ddsMat_rna, test = "Wald")
  # res_rna <- results(ddsMat_rna, name = 'Condition_2_vs_1')
  res_rna <- results(ddsMat_rna)
  if (shinkage) {
    res_rna <- lfcShrink(ddsMat_rna, res = res_rna, coef = 2,  type = "apeglm")
  }
  
  res_rna_norm <- res_rna %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'Gene') %>% 
    left_join(counts(ddsMat_rna, normalized = TRUE) %>% 
                as.data.frame() %>% 
                rownames_to_column('Gene'), by = 'Gene')
  
  ####################################################
  ## 5. Run Ribo DESeq2
  ddsMat_ribo <- DESeqDataSetFromMatrix(
    countData = ribo_counts,
    colData = sample_info[which(sample_info$SeqType == "RIBO"),],
    design =~ Condition)
  
  ddsMat_ribo <- DESeq(ddsMat_ribo)
  # res_ribo <- results(ddsMat_ribo, name = "Condition_2_vs_1")
  res_ribo <- results(ddsMat_ribo)
  if (shinkage) {
    res_ribo <- lfcShrink(ddsMat_ribo, res = res_ribo, coef = 2,  type = "apeglm")
  }
  
  res_ribo_norm <- res_ribo %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'Gene') %>% 
    left_join(counts(ddsMat_ribo, normalized = TRUE) %>% 
                as.data.frame() %>% 
                rownames_to_column('Gene'), by = 'Gene')

  ####################################################
  ## 6. filter the forward genes
  forwarded <- rownames(res)[which(res$padj > 0.05 & res_ribo$padj < 0.05 & res_rna$padj < 0.05)]
  exclusive <- rownames(res)[which(res$padj < 0.05 & res_ribo$padj < 0.05 & res_rna$padj > 0.05)]
  
  both <- which(res$padj < 0.05 & res_ribo$padj < 0.05 & res_rna$padj < 0.05)
  intensified <- rownames(res)[both[which(res[both, 2] * res_rna[both, 2] > 0)]]
  
  buffered <- rownames(res)[both[which(res[both, 2] * res_rna[both, 2] < 0)]]
  buffered <- c(rownames(res)[which(res$padj < 0.05 & res_ribo$padj > 0.05 & res_rna$padj < 0.05)], buffered)
  
  ####################################################
  ## 7. annotate the gene delta table
  res_norm <- res_norm %>% 
    mutate(Delta = if_else(Gene %in% forwarded, 'forwarded',
                           if_else(Gene %in% exclusive, 'exclusive',
                                   if_else(Gene %in% intensified, 'intensified',
                                           if_else(Gene %in% buffered, 'buffered', 'others')))))
  res_rna_norm <- res_rna_norm %>% 
    mutate(Delta = if_else(Gene %in% forwarded, 'forwarded',
                           if_else(Gene %in% exclusive, 'exclusive',
                                   if_else(Gene %in% intensified, 'intensified',
                                           if_else(Gene %in% buffered, 'buffered', 'others')))))
  res_ribo_norm <- res_ribo_norm %>% 
    mutate(Delta = if_else(Gene %in% forwarded, 'forwarded',
                           if_else(Gene %in% exclusive, 'exclusive',
                                   if_else(Gene %in% intensified, 'intensified',
                                           if_else(Gene %in% buffered, 'buffered', 'others')))))
  
  ####################################################
  ## 8. annotate the gene sig table
  res_norm <- res_norm %>% 
    mutate(DEGs = if_else(padj < padj_cutoff & log2FoldChange >= log2fc_cutoff, 'UP',
                           if_else(padj < padj_cutoff & log2FoldChange <= -log2fc_cutoff, 'DOWN', 'NS')))
  
  res_rna_norm <- res_rna_norm %>% 
    mutate(DEGs = if_else(padj < padj_cutoff & log2FoldChange >= log2fc_cutoff, 'UP',
                           if_else(padj < padj_cutoff & log2FoldChange <= -log2fc_cutoff, 'DOWN', 'NS')))
  
  res_ribo_norm <- res_ribo_norm %>% 
    mutate(DEGs = if_else(padj < padj_cutoff & log2FoldChange >= log2fc_cutoff, 'UP',
                           if_else(padj < padj_cutoff & log2FoldChange <= -log2fc_cutoff, 'DOWN', 'NS')))
  
  
  ####################################################
  ## 9. join the annotation to TE / mRNA / Ribo DESeq2 results
  # browser()
  
  if (class(anno) == "data.frame") {
    res_norm <- res_norm %>% 
      left_join(., anno, by = c('Gene' = join_flag))
    
    res_rna_norm <- res_rna_norm %>% 
      left_join(., anno, by = c('Gene' = join_flag))
    
    res_ribo_norm <- res_ribo_norm %>% 
      left_join(., anno, by = c('Gene' = join_flag))
  }
  
  ####################################################
  ## 10. summary the table
  
  # browser()
  
  summary_table <- rbind(
    table(res_norm$DEGs) %>% 
      as.data.frame() %>% 
      rename_all(~c('Class', 'Count')) %>% 
      mutate(Group = "TE", Condition = paste("log2FC:", log2fc_cutoff, "padj:", padj_cutoff)),
    table(res_rna_norm$DEGs) %>% 
      as.data.frame() %>% 
      rename_all(~c('Class', 'Count')) %>% 
      mutate(Group = "RNA", Condition = paste("log2FC:", log2fc_cutoff, "padj:", padj_cutoff)),
    table(res_ribo_norm$DEGs) %>% 
      as.data.frame() %>% 
      rename_all(~c('Class', 'Count')) %>% 
      mutate(Group = "RIBO", Condition = paste("log2FC:", log2fc_cutoff, "padj:", padj_cutoff))
  )
  
  ####################################################
  ## 11. output the merge results
  if (flt_sig) {
    delta_te_res <- list(
      DTEGs = res_norm %>% filter(padj < padj_cutoff, abs(log2FoldChange) >= log2fc_cutoff),
      rna_DEGs = res_rna_norm %>% filter(padj < padj_cutoff, abs(log2FoldChange) >= log2fc_cutoff),
      ribo_DEGs = res_ribo_norm %>% filter(padj < padj_cutoff, abs(log2FoldChange) >= log2fc_cutoff),
      forwarded = forwarded,
      exclusive = exclusive,
      intensified = intensified,
      buffered = buffered,
      summary_table = summary_table)
    
    # write.xlsx(x = delta_te_res, file = paste0(delte_te_path, outname, '.xlsx'))
  
  } else {
    delta_te_res <- list(
      DTEGs = res_norm,
      rna_DEGs = res_rna_norm,
      ribo_DEGs = res_ribo_norm,
      forwarded = forwarded,
      exclusive = exclusive,
      intensified = intensified,
      buffered = buffered,
      summary_table = summary_table)
    
    # write.xlsx(x = delta_te_sig_res, file = paste0(delte_te_path, outname, '_sig.xlsx'))
  }
  
  return(delta_te_res)
  
}