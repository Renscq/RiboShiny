###########################################
# run_KEGG_Terms.R
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
# library(AnnotationHub)
library(gson)


run_kegg_terms <- function(
    degs = NULL,
    degs_group = FALSE,
    
    database = "gson",
    species = NULL,
    kegg_type = "KEGG", # MKEGG, KEGG
    
    gson_file = NULL,
    
    keytype = 'kegg', #  one of "kegg", 'ncbi-geneid', 'ncbi-proteinid' and 'uniprot'
    
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.2,
    
    pAdjustMethod = "BH", # one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
    
    use_internal_data = FALSE
) {
  
  ######################################################
  # import the organism database
  # browser()
  
  if (database == "gson") {
    gson_kegg <- gson::read.gson(gson_file)
    
    ####################################################
    # enrich the kegg with down and up genes
    if (degs_group == TRUE) {
      KEGG_terms <- compareCluster(geneCluster = degs,
                                   fun = enricher,
                                   gson = gson_kegg,
                                   pAdjustMethod = pAdjustMethod,
                                   pvalueCutoff = pvalueCutoff,
                                   qvalueCutoff = qvalueCutoff)
      
    } else {
      KEGG_terms <- enricher(degs,
                             gson = gson_kegg,
                             pAdjustMethod = pAdjustMethod,
                             pvalueCutoff = pvalueCutoff,
                             qvalueCutoff = qvalueCutoff)
      
    }
    
  } else if (database == "kegg") {
    species <- species
    
    ####################################################
    # enrich the kegg with down and up genes
    if (degs_group == TRUE) {
      if (kegg_type == "MKEGG") {
        KEGG_terms <- compareCluster(geneCluster = degs,
                                     fun = enrichMKEGG,
                                     organism = species,
                                     keyType = keytype,
                                     pvalueCutoff = pvalueCutoff,
                                     qvalueCutoff = qvalueCutoff)
      } else if (kegg_type == "KEGG") {
        KEGG_terms <- compareCluster(geneCluster = degs,
                                     fun = enrichKEGG,
                                     organism = species,
                                     keyType = keytype,
                                     pvalueCutoff = pvalueCutoff,
                                     qvalueCutoff = qvalueCutoff)
      }
      
    } else {
      if (kegg_type == "MKEGG") {
        KEGG_terms <- enrichMKEGG(degs,
                                 organism = species,
                                 keyType = keytype,
                                 pAdjustMethod = pAdjustMethod,
                                 pvalueCutoff = pvalueCutoff,
                                 qvalueCutoff = qvalueCutoff,
                                 use_internal_data = use_internal_data)
      } else if (kegg_type == "KEGG") {
        KEGG_terms <- enrichKEGG(degs,
                                 organism = species,
                                 keyType = keytype,
                                 pAdjustMethod = pAdjustMethod,
                                 pvalueCutoff = pvalueCutoff,
                                 qvalueCutoff = qvalueCutoff,
                                 use_internal_data = use_internal_data)
      }
    }
    
  }
  
  return(KEGG_terms)
}

