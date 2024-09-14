###########################################
# run_KEGG_GSEA_Terms.R
#
# This script enrich the GO terms for the DEGs
#
# degs: the character vector of list of the DEGs
# database: the organism database
# key_type: the key type of the DEGs
# pvalueCutoff: the p-value cutoff
# qvalueCutoff: the q-value cutoff
# ontology: the ontology of the KEGG terms
#
# usage: run_kegg_gsea_terms(degs, orgdb, key_type, pvalueCutoff, qvalueCutoff, ontology, pAdjustMethod, readable)


library(clusterProfiler)
library(gson)


run_kegg_gsea_terms <- function(
    degs = NULL,
    
    database = "gson",
    gson_file = NULL,
    
    species = NULL,
    kegg_type = "KEGG", # MKEGG, KEGG
    
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
    gson_kegg <- gson::read.gson(gson_file)

    # browser()
    kegg_gsea_terms <- GSEA(geneList = degs,
                          exponent = exponent,
                          eps = eps,
                          
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          
                          gson = gson_kegg,
                          seed = FALSE,
                          by = "fgsea")

  } else if (database == "kegg") {

    # gson_kegg <- gson_KEGG(species = species, KEGG_Type = kegg_type, keyType = keytype)

    # browser()
    if (kegg_type == "MKEGG") {
      kegg_gsea_terms <- gseMKEGG(geneList = degs,
                                  organism = species,
                                  keyType = keytype,
                                  exponent = exponent,
                                  eps = eps,
                                  minGSSize = minGSSize,
                                  maxGSSize = maxGSSize,
                                  pvalueCutoff = pvalueCutoff,
                                  pAdjustMethod = pAdjustMethod,
                                  seed = FALSE,
                                  by = "fgsea")
      
    } else if (kegg_type == "KEGG") {
      kegg_gsea_terms <- gseKEGG(geneList = degs,
                                 organism = species,
                                 keyType = keytype,
                                 exponent = exponent,
                                 eps = eps,
                                 minGSSize = minGSSize,
                                 maxGSSize = maxGSSize,
                                 pvalueCutoff = pvalueCutoff,
                                 pAdjustMethod = pAdjustMethod,
                                 seed = FALSE,
                                 by = "fgsea")
    }

  }

  return(kegg_gsea_terms)
}

