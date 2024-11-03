###############################################
# calc_Enrichment.R
# 

calc_enrichment <- function(
    meta_table = NULL, # metagene data frame contain the region and density
    group = NULL, # group data frame contain the group information
    group1 = "IP", # IP, MOCK, WT
    group2 = "WT", # IP, MOCK, WT
    category = "IP/WT", # IP/WT, IP/MOCK, WT/MOCK
    method = "mix" # mix, individual
) {
  
  # browser()
  meta_enrich <- NULL

  ##############################################
  # retrieve individual group1 and group2 samples for enrichment
  # browser()
  
  if (method == "individual") {

    if (length(group1) != length(group2)) {
      return(NULL)
    }
    
    for (now_sample in 1:length(group1)) {
      now_g1 = group1[now_sample]
      now_g2 = group2[now_sample]
      
      now_g1_dat = meta_table %>% dplyr::filter(!!sym(group) == now_g1)
      now_g2_dat = meta_table %>% dplyr::filter(!!sym(group) == now_g2)
      
      Enrichment = now_g1_dat$Density / now_g2_dat$Density
      now_g1_g2 = cbind(now_g1_dat[, c("Meta", "Codon", "Nucleotide")], Enrichment)
  
      now_g1_g2$Samples1 = now_g1
      now_g1_g2$Samples2 = now_g2
  
      now_g1_g2$Density1 = now_g1_dat$Density
      now_g1_g2$Density2 = now_g2_dat$Density
  
      now_g1_g2$Category = category
      now_g1_g2$Paired = now_sample
  
      # merge the data
      if (is.null(meta_enrich)) {
        meta_enrich = now_g1_g2
      } else {
        meta_enrich = rbind(meta_enrich, now_g1_g2)
      }
  
    }
    
  } else if (method == "mix") {
    ##############################################
    # retrieve all the group1 and group2 samples for enrichment
    # browser()
    now_sample = 1
    for (now_g1 in group1) {
      now_g1_dat = meta_table %>% dplyr::filter(!!sym(group) == now_g1)
      
      for (now_g2 in group2) {
        now_g2_dat = meta_table %>% dplyr::filter(!!sym(group) == now_g2)
        
        Enrichment = now_g1_dat$Density / now_g2_dat$Density
        now_g1_g2 = cbind(now_g1_dat[, c("Meta", "Codon", "Nucleotide")], Enrichment)
        
        now_g1_g2$Samples1 = now_g1
        now_g1_g2$Samples2 = now_g2
        
        now_g1_g2$Density1 = now_g1_dat$Density
        now_g1_g2$Density2 = now_g2_dat$Density
        
        now_g1_g2$Category = category
        now_g1_g2$Paired = now_sample
        now_sample = now_sample + 1
        
        # merge the data
        if (is.null(meta_enrich)) {
          meta_enrich = now_g1_g2
        } else {
          meta_enrich = rbind(meta_enrich, now_g1_g2)
        }
        
      }
      
    }
    
  }

  
  return(meta_enrich)
  
}