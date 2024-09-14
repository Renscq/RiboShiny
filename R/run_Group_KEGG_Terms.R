###########################################
# run_Group_KEGG_Terms.R
#
# This script enrich the KEGG terms for the DEGs
#
# degs: the character vector of list of the DEGs
# degs_group: whether the DEGs are grouped
# orgdb: the organism database
# key_type: the key type of the DEGs
# pvalueCutoff: the p-value cutoff
# qvalueCutoff: the q-value cutoff
# ontology: the ontology of the KEGG terms
# pAdjustMethod: the method of p-value adjustment
# readable: whether to return the readable KEGG terms
#
# usage: run_KEGG_terms(degs, orgdb, key_type, pvalueCutoff, qvalueCutoff, ontology, pAdjustMethod, readable)


library(clusterProfiler)
library(gson)


run_group_kegg_terms <- function(
    degs = NULL,
    degs_id = NULL, # GID or ENTREZID or ENSEMBL
    degs_group = NULL, # DEGs, Groups, Files
    
    database = "gson",
    species = NULL,
    gson_file = NULL,
    
    keytype = 'kegg', #  one of "kegg", 'ncbi-geneid', 'ncbi-proteinid' and 'uniprot'
    
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.2,
    
    pAdjustMethod = "BH", # one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
    
    use_internal_data = FALSE
) {
  
  ######################################################
  # set the groups formula
  # browser()
  degs_group = paste0(degs_group, collapse = '+')
  group_formula <- as.formula(paste(degs_id, "~", degs_group))

  ######################################################
  # import the organism database
  # browser()
  
  if (database == "gson") {
    gson_kegg <- gson::read.gson(gson_file)
    
    ####################################################
    # enrich the kegg with down and up genes
    KEGG_terms <- compareCluster(data = degs,
                                 geneClusters = group_formula,
                                 fun = enricher,
                                 gson = gson_kegg,
                                 pAdjustMethod = pAdjustMethod,
                                 pvalueCutoff = pvalueCutoff,
                                 qvalueCutoff = qvalueCutoff)
    
    
  } else if (database == "kegg") {
    species <- species
    
    ####################################################
    # enrich the kegg with down and up genes
    KEGG_terms <- compareCluster(data = degs,
                                 geneClusters = group_formula,
                                 fun = enrichKEGG,
                                 organism = species,
                                 keyType = keytype,
                                 pvalueCutoff = pvalueCutoff,
                                 qvalueCutoff = qvalueCutoff)
    
  }
  
  return(KEGG_terms)
}

