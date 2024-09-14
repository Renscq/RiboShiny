###########################################
# run_Group_GO_Terms.R
#
# This script enrich the GO terms for the DEGs
#
# degs: the data frame of the DEGs
# degs_id: the ID of the DEGs
# degs_group: whether the DEGs are grouped
# database: the database of the GO, gson or orgdb
# orgdb: the organism database
# gson_file: the gson file
# key_type: the key type of the DEGs
# pvalueCutoff: the p-value cutoff
# qvalueCutoff: the q-value cutoff
# ontology: the ontology of the GO terms
# pAdjustMethod: the method of p-value adjustment
# readable: whether to return the readable GO terms
#
# usage: run_group_go_terms(degs, degs_group, orgdb, keytype, pvalueCutoff, qvalueCutoff, ontology, pAdjustMethod, readable)


library(clusterProfiler)
library(gson)


run_group_go_terms <- function(
    degs = NULL,
    degs_id = NULL, # GID or ENTREZID or ENSEMBL
    degs_group = NULL, # DEGs, Groups, Files
    
    database = "gson",
    orgdb = NULL,
    gson_file = NULL,
    
    keytype = NULL, # ENTREZID or GID
    
    pvalueCutoff = NULL,
    qvalueCutoff = NULL,
    
    ontology = FALSE, # One of "BP", "MF", and "CC", or "ALL" 
    pAdjustMethod = FALSE, # one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
    
    readable = FALSE
) {
  
  ######################################################
  # set the groups formula
  # browser()
  degs_group = paste0(degs_group, collapse = '+')
  group_formula <- as.formula(paste0(degs_id, "~", degs_group))
  
  ######################################################
  # import the organism database
  if (database == "gson") {
    gson_file <- gson::read.gson(orgdb)
    
    ######################################################
    # Enrich the GO terms
    # browser()
    
    go_terms <- compareCluster(data = degs,
                               geneClusters = group_formula,
                               fun = enricher,
                               gson = gson_file,
                               pAdjustMethod = pAdjustMethod,
                               pvalueCutoff = pvalueCutoff,
                               qvalueCutoff = qvalueCutoff)
    
    go_terms_ont <- go2ont(go_terms@compareClusterResult$ID) %>% 
      dplyr::rename_all(~c('ID', 'ONTOLOGY'))
    
    go_terms@compareClusterResult <- go_terms@compareClusterResult %>% 
      dplyr::left_join(go_terms_ont, by = "ID")

  } else if (database == "orgdb") {
    library(orgdb, character.only = TRUE)
    
    ######################################################
    # Enrich the GO terms
    # browser()
    
    go_terms <- compareCluster(data = degs,
                               geneClusters = group_formula,
                               fun = enrichGO, 
                               OrgDb = orgdb, 
                               keyType = keytype, 
                               ont = ontology, 
                               pvalueCutoff = pvalueCutoff,
                               qvalueCutoff = qvalueCutoff)
  }
  
  return(go_terms)
}

