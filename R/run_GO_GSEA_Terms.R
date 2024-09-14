###########################################
# run_GO_GSEA_Terms.R
#
# This script enrich the GO terms for the DEGs
#
# degs: the character vector of list of the DEGs
# degs_group: whether the DEGs are grouped
# orgdb: the organism database
# key_type: the key type of the DEGs
# pvalueCutoff: the p-value cutoff
# qvalueCutoff: the q-value cutoff
# ontology: the ontology of the GO terms
# pAdjustMethod: the method of p-value adjustment
# readable: whether to return the readable GO terms
#
# usage: run_go_gsea_terms(degs, orgdb, key_type, pvalueCutoff, qvalueCutoff, ontology, pAdjustMethod, readable)


library(clusterProfiler)
library(gson)


run_go_gsea_terms <- function(
    degs = NULL,
    
    database = "gson",
    gson_file = NULL,
    orgdb = NULL,
    ontology = "ALL", # One of "BP", "MF", and "CC", or "ALL" 
    
    keytype = NULL, # ENTREZID or GID
    
    exponent = 1,
    eps = 1e-10,
    
    minGSSize = 10,
    maxGSSize = 500,
    
    pAdjustMethod = FALSE, # one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
    pvalueCutoff = NULL
) {
  
  ######################################################
  # import the organism database
  if (database == "gson") {
    gson_go <- gson::read.gson(gson_file)

    # browser()
    go_gsea_terms <- GSEA(geneList = degs,
                          exponent = exponent,
                          eps = eps,
                          
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          
                          gson = gson_go,
                          seed = FALSE,
                          by = "fgsea")
    
    go_terms_ont <- go2ont(go_gsea_terms@result$ID) %>% 
      dplyr::rename_all(~c('ID', 'ONTOLOGY'))
    
    go_gsea_terms@result <- go_gsea_terms@result %>% 
      dplyr::left_join(go_terms_ont, by = "ID") %>% 
      dplyr::select("ONTOLOGY", "ID", "Description", "setSize", "enrichmentScore", "NES",
                    "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment")
    
    rownames(go_gsea_terms@result) <- go_gsea_terms@result$ID
    
  } else if (database == "orgdb") {
    library(orgdb, character.only = TRUE)
    
    # browser()
    
    tryCatch({
      go_gsea_terms <- gseGO(geneList = degs,
                             ont = "ALL",
                             OrgDb = orgdb,
                             keyType = "ENTREZID",
                             exponent = exponent,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             eps = eps,
                             pvalueCutoff = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             verbose = TRUE,
                             seed = FALSE,
                             by = "fgsea")

    }, error = function(e) {
      
      gson_go <- gson_GO(OrgDb = orgdb, keytype = keytype, ont = 'ALL')
      
      go_gsea_terms <- GSEA(geneList = degs,
                            exponent = exponent,
                            minGSSize = minGSSize,
                            maxGSSize = maxGSSize,
                            eps = eps,
                            pvalueCutoff = pvalueCutoff,
                            pAdjustMethod = pAdjustMethod,
                            gson = gson_go,
                            seed = FALSE,
                            by = "fgsea")
      
      # browser()
      
      go_terms_ont <- go2ont(go_gsea_terms@result$ID) %>% 
        dplyr::rename_all(~c('ID', 'ONTOLOGY'))
      
      go_gsea_terms@result <- go_gsea_terms@result %>% 
        dplyr::left_join(go_terms_ont, by = "ID")
      
    })
    
    # browser()
    
    go_gsea_terms@result <- go_gsea_terms@result %>% 
      dplyr::select("ONTOLOGY", "ID", "Description", "setSize", "enrichmentScore", "NES",
                    "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment")
    
    rownames(go_gsea_terms@result) <- go_gsea_terms@result$ID
    
  }

  return(go_gsea_terms)
}

