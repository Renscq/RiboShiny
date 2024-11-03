###########################################
# run_GO_Terms.R
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
# usage: run_go_terms(degs, orgdb, key_type, pvalueCutoff, qvalueCutoff, ontology, pAdjustMethod, readable)


library(clusterProfiler)
library(gson)


run_go_terms <- function(
    degs = NULL,
    degs_group = FALSE,
    
    database = "gson",
    gson_file = NULL,
    orgdb = NULL,
    
    keytype = NULL, # ENTREZID or GID
    
    pvalueCutoff = NULL,
    qvalueCutoff = NULL,
    
    ontology = FALSE, # One of "BP", "MF", and "CC", or "ALL" 
    pAdjustMethod = FALSE, # one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
    
    readable = FALSE
) {
  
  ######################################################
  # import the organism database
  if (database == "gson") {
    gson_go <- gson::read.gson(gson_file)
    
    ######################################################
    # Enrich the GO terms
    # browser()
    
    if (degs_group == TRUE) {
      go_terms <- compareCluster(geneCluster = degs,
                                 fun = enricher,
                                 gson = gson_go,
                                 pAdjustMethod = pAdjustMethod,
                                 pvalueCutoff = pvalueCutoff,
                                 qvalueCutoff = qvalueCutoff)
      
      go_terms_ont <- go2ont(go_terms@compareClusterResult$ID) %>% 
        dplyr::rename_all(~c('ID', 'ONTOLOGY'))
      
      go_terms@compareClusterResult <- go_terms@compareClusterResult %>% 
        dplyr::left_join(go_terms_ont, by = "ID") %>% 
        dplyr::select("Cluster", "ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio",
                      "pvalue", "p.adjust", "qvalue", "geneID", "Count")
      
      rownames(go_terms@compareClusterResult) <- go_terms@compareClusterResult$ID
      
    } else {
      go_terms <- enricher(degs,
                           gson = gson_go,
                           pAdjustMethod = pAdjustMethod,
                           pvalueCutoff = pvalueCutoff,
                           qvalueCutoff = qvalueCutoff)
      
      go_terms_ont <- go2ont(go_terms@result$ID) %>% 
        dplyr::rename_all(~c('ID', 'ONTOLOGY'))
      
      go_terms@result <- go_terms@result %>% 
        dplyr::left_join(go_terms_ont, by = "ID")
      go_terms@result$Cluster <- 'DEGs'
      go_terms@result <- go_terms@result %>% 
        dplyr::select("Cluster", "ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio",
                      "pvalue", "p.adjust", "qvalue", "geneID", "Count")
      
      rownames(go_terms@result) <- go_terms@result$ID
      
    }
    
  } else if (database == "orgdb") {
    library(orgdb, character.only = TRUE)
    
    ######################################################
    # Enrich the GO terms
    # browser()
    
    if (degs_group == TRUE) {
      go_terms <- compareCluster(geneCluster = degs, 
                                 fun = enrichGO, 
                                 OrgDb = orgdb, 
                                 keyType = keytype, 
                                 ont = ontology, 
                                 pAdjustMethod = pAdjustMethod,
                                 pvalueCutoff = pvalueCutoff,
                                 qvalueCutoff = qvalueCutoff)
      
    } else {
      go_terms <- enrichGO(degs,
                           OrgDb = orgdb,
                           keyType = keytype,
                           ont = ontology,
                           pAdjustMethod = pAdjustMethod,
                           pvalueCutoff = pvalueCutoff,
                           qvalueCutoff = qvalueCutoff,
                           readable = readable)
      
      go_terms@result$Cluster <- 'DEGs'
    }
  }

  return(go_terms)
}

