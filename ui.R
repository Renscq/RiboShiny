#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
############################ packages ############################
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager");
if (!require("devtools")) 
  install.packages("devtools");


# shiny and dashboard packages
if (!require("shiny")) 
  install.packages("shiny");
if (!require("shinyFiles")) 
  install.packages("shinyFiles");
if (!require("shinycssloaders"))
  install.packages("shinycssloaders");
if (!require("shinydashboard")) 
  install.packages("shinydashboard");
if (!require("dashboardthemes")) 
  install.packages("dashboardthemes");
# if (!require("shinySignals")) 
#   devtools::install_github("hadley/shinySignals");


# import and export data
if (!require("openxlsx")) 
  install.packages("openxlsx");
if (!require("DT")) 
  install.packages("DT");
if (!require("data.table")) 
  install.packages("data.table");
if (!require("showtext")) 
  install.packages("showtext");

# data analysis
if (!require("tidyverse")) 
  install.packages("tidyverse");
if (!require("zoo"))
  install.packages("zoo");
if (!require("FactoMineR"))
  install.packages("FactoMineR");
if (!require("factoextra"))
  install.packages("factoextra");


# data visualization
if (!require("ggplot2")) 
  install.packages("ggplot2");
# if (!require("ggpubr")) 
#   install.packages("ggpubr");
if (!require("MetBrewer")) 
  install.packages("MetBrewer")
if (!require("RColorBrewer")) 
  install.packages("RColorBrewer");
if (!require("pheatmap"))
  install.packages("pheatmap");
if (!require("minpack.lm")) 
  install.packages("minpack.lm");
# if (!require("ComplexHeatmap")) 
#   install.packages("ComplexHeatmap");
if (!require("export"))
  install.packages("export");
if (!require("aplot"))
  install.packages("aplot");
# if (!require("cowplot"))
#   install.packages("cowplot");
if (!require("ggplotify"))
  install.packages("ggplotify");
# if (!require("ggrepel"))
#   install.packages("ggrepel");
if (!require("ggnewscale"))
  install.packages("ggnewscale");
if (!require("ggseqlogo"))
  install.packages("ggseqlogo");
if (!require("ggVennDiagram"))
  install.packages("ggVennDiagram");
if (!require("UpSetR"))
  install.packages("UpSetR");

# BiocManager
if (!require("Biostrings"))
  BiocManager::install("Biostrings");
if (!require("coRdon")) 
  BiocManager::install("coRdon");
# if (!require("RSamtools"))
#   BiocManager::install("RSamtools");
if (!require("ShortRead")) 
  BiocManager::install("ShortRead");
if (!require("DESeq2"))
  BiocManager::install("DESeq2");
if (!require("edgeR"))
  BiocManager::install("edgeR");
if (!require("RUVSeq"))
  BiocManager::install("RUVSeq");
# if (!require("GOSE"))
#   BiocManager::install("GOSE");
if (!require("ReactomePA"))
  BiocManager::install("ReactomePA");
if (!require("clusterProfiler"))
  BiocManager::install("clusterProfiler");
if (!require("AnnotationDbi"))
  BiocManager::install("AnnotationDbi");
if (!require("AnnotationHub"))
  BiocManager::install("AnnotationHub");
if (!require("AnnotationForge"))
  BiocManager::install("AnnotationForge");


############################ packages ############################
suppressMessages(require("shiny"))
suppressMessages(require("shinyFiles"))
suppressMessages(require("shinycssloaders"))
suppressMessages(require("shinydashboard"))
suppressMessages(require("dashboardthemes"))
# suppressMessages(require("shinySignals"))
suppressMessages(require("openxlsx"))
suppressMessages(require("DT"))
suppressMessages(require("data.table"))
suppressMessages(require("zoo"))
suppressMessages(require("showtext"))
suppressMessages(require("tidyverse"))
suppressMessages(require("Biostrings"))
suppressMessages(require("coRdon"))
suppressMessages(require("FactoMineR"))
suppressMessages(require("factoextra"))
suppressMessages(require("ggplot2"))
suppressMessages(require("ggplotify"))
suppressMessages(require("ggpubr"))
suppressMessages(require("RColorBrewer"))
suppressMessages(require("pheatmap"))
suppressMessages(require("minpack.lm"))
suppressMessages(require("export"))
suppressMessages(require("aplot"))
suppressMessages(require("ggrepel"))
suppressMessages(require("ggnewscale"))
suppressMessages(require("ggVennDiagram"))
suppressMessages(require("VennDiagram"))
suppressMessages(require("UpSetR"))
suppressMessages(require("AnnotationHub"))
suppressMessages(require("AnnotationDbi"))
suppressMessages(require("AnnotationForge"))




############################ set the default options ############################
# set the color list
color_list <- c("Set1", "Set2", "Set3", "Dark2", "Paired", "Accent", "Pastel1", "Pastel2", "Spectral", 
                "Blues", "Greens", "Reds", "Oranges", "Purples", "Greys",
                'Purples', 'YlGnBu', 'Spectral', 'RdYlGn', 'RdYlBu', 'RdBu', 'PRGn', 'PiYG')

heat_color_list <- c('Purples', 'Blues', 'YlGnBu', 'Spectral', 'RdYlGn', 'RdYlBu', 'RdBu', 'PRGn', 'PiYG')

gradient_colors <- c("Blues", "Greens", "Reds", "Oranges", "Purples", "Greys")

seqlogo_colors <- c("auto", "chemistry", "chemistry2", "hydrophobicity", 
                    "nucleotide", "nucleotide2", "base_pairing", "clustalx", "taylor")

seqlogo_font_family <- c("helvetica_regular", "helvetica_bold", "helvetica_light", "roboto_medium",
                         "roboto_bold", "roboto_regular", "akrobat_bold", "akrobat_regular", 
                         "roboto_slab_bold", "roboto_slab_regular", "roboto_slab_light", "xkcd_regular")

genetic_code_index <- c("SGC0", "SGC1", "SGC2", "SGC3", "SGC4", "SGC5", "SGC8", "SGC9")
genetic_code_anno <- c("Standard", "Vertebrate Mitochondrial", "Yeast Mitochondrial", 
                       "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma",
                       "Invertebrate Mitochondrial", "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
                       "Echinoderm Mitochondrial; Flatworm Mitochondrial", "Euplotid Nuclear")
genetic_code_list <- paste(genetic_code_index, genetic_code_anno, sep = " - ")

default_codon_list <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", 
                        "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", 
                        "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", 
                        "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", 
                        "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", 
                        "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", 
                        "TTA", "TTC", "TTG", "TTT")


dashboardPage(
  
  ############################ dashboardHeader ############################
  
  dashboardHeader(title = "RiboShiny (v0.1.0)",
                  titleWidth = 300,
                  ## set the message menu ############################
                  dropdownMenu(type = "messages",
                               icon = icon("user"),
                               messageItem(
                                 from = "Ren Shuchao",
                                 message = "PhD in bioinformatics."
                               ),
                               messageItem(
                                 from = "E-mail",
                                 message = "rensc0718@163.com",
                                 icon = icon("envelope")
                               ),
                               messageItem(
                                 from = "CopyRight",
                                 message = "All copyrights reserved.",
                                 # time = "2024-06-30",
                                 icon = icon("copyright")
                               )
                  ),
                  
                  dropdownMenu(type = "notifications",
                               icon = icon("quote-left"),
                               messageItem(
                                 message  = "RiboParer and RiboShiny: A comprehensive framework for RNA-seq and Ribo-seq data mining and visualization",
                                 href = "https://www.bing.com",
                                 icon("link")
                               )
                  )
  ),
  
  
  
  ############################ dashboardSidebar ############################
  
  dashboardSidebar(
    sidebarMenu(
      ## step0 OverView ############################
      menuItem(text = "OverView", tabName = "OverView", icon = icon("table")
               # menuSubItem(text = "GeneLevel", tabName = "GeneLevel"),
               # menuSubItem(text = "CodonLevel", tabName = "CodonLevel")
      ),
      
      ## step1_Database ############################
      menuItem(text = "step1_Database", tabName = "step1_Database", icon = icon("table"),
               # menuSubItem(text = "ip", tabName = "step1_ip"),
               menuSubItem(text = "Options", tabName = "step1_options"),
               menuSubItem(text = "Annotation", tabName = "step1_annotation"),
               menuSubItem(text = "OrgDb", tabName = "step1_orgdb"),
               menuSubItem(text = "Gson_KEGG", tabName = "step1_kegg"),
               menuSubItem(text = "Sequence", tabName = "step1_sequence")
      ),
      
      ## step2_Mapping ############################
      menuItem(text = "step2_Mapping", tabName = "step2_Mapping", icon = icon("chart-bar"),
               # menuSubItem(text = "RNAseq", tabName = "step2_rnaseq"),
               menuSubItem(text = "RiboSeq", tabName = "step2_riboseq")
      ),
      
      ## step3_Distribution ############################
      menuItem(text = "step3_Distribution", tabName = "step3_Distribution", icon = icon("chart-bar"),
               menuSubItem(text = "RNAseq", tabName = "step3_rnaseq"),
               menuSubItem(text = "RiboSeq", tabName = "step3_riboseq")
      ),
      
      ## step4_Saturation ############################
      menuItem(text = "step4_Saturation", tabName = "step4_Saturation", icon = icon("chart-bar"),
               menuSubItem(text = "RNAseq", tabName = "step4_rnaseq"),
               menuSubItem(text = "RiboSeq", tabName = "step4_riboseq")
      ),
      
      ## step5_Digestion ############################
      menuItem(text = "step5_Digestion", tabName = "step5_Digestion", icon = icon("chart-bar"),
               # menuSubItem(text = "RNAseq", tabName = "step5_rnaseq"),
               menuSubItem(text = "RiboSeq", tabName = "step5_riboseq")
      ),
      
      ## step6_Offset ############################
      menuItem(text = "step6_Offset", tabName = "step6_Offset", icon = icon("chart-bar"),
               menuSubItem(text = "RSBM", tabName = "step6_frame"),
               menuSubItem(text = "SSCBM", tabName = "step6_tis"),
               menuSubItem(text = "Detail", tabName = "step6_end")
      ),
      
      ## step7_Period ############################
      menuItem(text = "step7_Period", tabName = "step7_Period", icon = icon("chart-bar"),
               menuSubItem(text = "RNAseq", tabName = "step7_rnaseq"),
               menuSubItem(text = "RiboSeq", tabName = "step7_riboseq")
      ),
      
      ## step8_Metaplot ############################
      menuItem(text = "step8_Metaplot", tabName = "step8_Metaplot", icon = icon("chart-bar"),
               menuSubItem(text = "RNAseq", tabName = "step8_rnaseq"),
               menuSubItem(text = "RiboSeq", tabName = "step8_riboseq")
      ),
      
      ## step9_Coverage ############################
      menuItem(text = "step9_Coverage", tabName = "step9_Coverage", icon = icon("chart-line"),
               menuSubItem(text = "RNAseq", tabName = "step9_rnaseq"),
               menuSubItem(text = "RiboSeq", tabName = "step9_riboseq")
      ),
      
      ## step10_Expression ############################
      menuItem(text = "step10_Expression", tabName = "step10_Expression", icon = icon("table"),
               menuSubItem(text = "RNAseq", tabName = "step10_rnaseq"),
               menuSubItem(text = "RiboSeq", tabName = "step10_riboseq"),
               menuSubItem(text = "TE", tabName = "step10_te")
               # only show the raw and normalized table here, show the boxplot and cdf plot
      ),
      
      ## step11_Heatmap ############################
      menuItem(text = "step11_Heatmap", tabName = "step11_Heatmap", icon = icon("chart-bar"),
               menuSubItem(text = "RNAseq", tabName = "step11_rnaseq"),
               menuSubItem(text = "RiboSeq", tabName = "step11_riboseq"),
               menuSubItem(text = "TE", tabName = "step11_te")
               # show the heatmap and correlation plot
      ),
      
      ## step12_PCA ############################
      menuItem(text = "step12_PCA", tabName = "step12_PCA", icon = icon("chart-bar"),
               menuSubItem(text = "RNAseq", tabName = "step12_rnaseq"),
               menuSubItem(text = "RiboSeq", tabName = "step12_riboseq"),
               menuSubItem(text = "TE", tabName = "step12_te")
      ),
      
      ## step13_deltaTE ############################
      menuItem(text = "step13_deltaTE", tabName = "step13_deltaTE", icon = icon("chart-bar"),
               menuSubItem(text = "edgeR", tabName = "step13_edger"),
               menuSubItem(text = "DESeq2", tabName = "step13_deseq2"),
               menuSubItem(text = "DeltaTE", tabName = "step13_delta")
      ),
      
      ## step14_volcano ############################
      menuItem(text = "step14_volcano", tabName = "step14_volcano", icon = icon("chart-bar"),
               menuSubItem(text = "Volcano", tabName = "step14_rna_ribo"),
               menuSubItem(text = "Quadrant", tabName = "step14_quadrant"),
               menuSubItem(text = "DeltaTE", tabName = "step14_deltate")
      ),
      
      ## step15_DEGs ############################
      menuItem(text = "step15_DEGs", tabName = "step15_DEGs", icon = icon("chart-bar"),
               menuSubItem(text = "Venn_diagram", tabName = "step15_venn"),
               menuSubItem(text = "DEGs_Venn", tabName = "step15_degs_venn"),
               menuSubItem(text = "DEGs_Merge", tabName = "step15_merge"),
               menuSubItem(text = "DEGs_Bar", tabName = "step15_barplot")
      ),
      
      ## step16_GO ############################
      menuItem(text = "step16_GO", tabName = "step16_GO", icon = icon("chart-bar"),
               menuSubItem(text = "EnrichGO", tabName = "step16_enrich_go"),
               menuSubItem(text = "CompareCluster", tabName = "step16_compare_cluster")
      ),
      
      ## step17_GO_GSEA ############################
      menuItem(text = "step17_GO_GSEA", tabName = "step17_GO_GSEA", icon = icon("chart-bar"),
               menuSubItem(text = "GSEA", tabName = "step17_go_gsea")
      ),
      
      ## step18_KEGG ############################
      menuItem(text = "step18_KEGG", tabName = "step18_KEGG", icon = icon("chart-bar"),
               menuSubItem(text = "EnrichKEGG", tabName = "step18_enrich_kegg"),
               menuSubItem(text = "CompareCluster", tabName = "step18_compare_cluster")
      ),
      
      ## step19_KEGG_GSEA ############################
      menuItem(text = "step19_KEGG_GSEA", tabName = "step19_KEGG_GSEA", icon = icon("chart-bar"),
               menuSubItem(text = "GSEA", tabName = "step19_kegg_gsea")
      ),
      
      ## step20_Codon_Enrich ############################
      menuItem(text = "step20_Codon_Enrich", tabName = "step20_Codon_Enrich", icon = icon("chart-line"),
               menuSubItem(text = "Codon_Usage", tabName = "step20_codon_usage"),
               menuSubItem(text = "GO_Enrich", tabName = "step20_go_gsea_enrich"),
               menuSubItem(text = "KEGG_Enrich", tabName = "step20_kegg_gsea_enrich")
      ),
      
      ## step21_Gene_plot ############################
      menuItem(text = "step21_Gene_plot", tabName = "step21_Gene_plot", icon = icon("chart-bar"),
               menuSubItem(text = "Isoforms", tabName = "step21_isoforms"),
               menuSubItem(text = "Gene", tabName = "step21_gene")
      ),
      
      ## step22_Pausing ############################
      menuItem(text = "step22_Pausing", tabName = "step22_Pausing", icon = icon("chart-bar"),
               menuSubItem(text = "Pausing", tabName = "step22_pausing"),
               menuSubItem(text = "Diff_pausing", tabName = "step22_diff_pausing")
      ),
      
      ## step23_Occupancy ############################
      menuItem(text = "step23_Occupancy", tabName = "step23_Occupancy", icon = icon("chart-bar"),
               menuSubItem(text = "Occupancy", tabName = "step23_occupancy"),
               menuSubItem(text = "Diff_occupancy", tabName = "step23_diff_occupancy")
      ),
      
      ## step24_CDT ############################
      menuItem(text = "step24_CDT", tabName = "step24_CDT", icon = icon("chart-bar"),
               menuSubItem(text = "CDT", tabName = "step24_cdt"),
               menuSubItem(text = "Diff_CDT", tabName = "step24_diff_cdt")
      ),
      
      ## step25_CST ############################
      menuItem(text = "step25_CST", tabName = "step25_CST", icon = icon("chart-bar"),
               menuSubItem(text = "CST", tabName = "step25_cst"),
               menuSubItem(text = "Diff_CST", tabName = "step25_diff_cst"),
               menuSubItem(text = "Iterative", tabName = "step25_iterative_cst")
      ),
      
      ## step26_OddRatio ############################
      menuItem(text = "step26_OddRatio", tabName = "step26_OddRatio", icon = icon("chart-bar"),
               menuSubItem(text = "OddRatio", tabName = "step26_odd_ratio")
      ),
      
      ## step27_CoV ############################
      menuItem(text = "step27_CoV", tabName = "step27_CoV", icon = icon("chart-bar"),
               menuSubItem(text = "CoV", tabName = "step27_cov"),
               menuSubItem(text = "CoV_eCDF", tabName = "step27_cov_cdf")
      ),
      
      ## step28_MetaCodon ############################
      menuItem(text = "step28_MetaCodon", tabName = "step28_MetaCodon", icon = icon("chart-bar"),
               menuSubItem(text = "MetaPlot", tabName = "step28_meta_codon_plot"),
               menuSubItem(text = "SeqLogo", tabName = "step28_meta_codon_seqlogo")
      ),
      
      ## step29_SeRP enrichment ############################
      menuItem(text = "step29_SeRP_Enrich", tabName = "step29_SeRP_Enrich", icon = icon("chart-bar"),
               # menuSubItem(text = "SeRP_Meta", tabName = "step29_serp_meta"),
               menuSubItem(text = "SeRP_Enrich", tabName = "step29_serp_enrich")
      ),
      
      ## step30_SeRP Peaks ############################
      menuItem(text = "step30_SeRP_Peaks", tabName = "step30_SeRP_Peaks", icon = icon("chart-bar"),
               menuSubItem(text = "SeRP_Peaks", tabName = "step30_serp_peaks")
      ),
      
      ## step31_SeRP motif ############################
      menuItem(text = "step31_SeRP_Motif", tabName = "step31_SeRP_Motif", icon = icon("chart-bar"),
               menuSubItem(text = "SeRP_Motif", tabName = "step31_serp_motif")
      )
      
      # ## step29_smORF ############################
      # menuItem(text = "step32_smORF", tabName = "step32_smORF", icon = icon("chart-bar"),
      #          menuSubItem(text = "Length", tabName = "step32_smORF_length"),
      #          menuSubItem(text = "Codon_Usage", tabName = "step32_smORF_codon_usage"),
      #          menuSubItem(text = "Start_Codon", tabName = "step32_smORF_start_codon")
      # ),
      # 
      # ## step30_smORF_Expr ############################
      # menuItem(text = "step33_smORF_Expr", tabName = "step33_smORF_Expr", icon = icon("chart-bar"),
      #          menuSubItem(text = "Correlation", tabName = "step33_smORF_Expr_correlation"),
      #          menuSubItem(text = "Cluster", tabName = "step33_smORF_Expr_cluster"),
      #          menuSubItem(text = "PCA", tabName = "step33_smORF_Expr_pca")
      # ),
      # 
      # ## step31_smORF_mORF ############################
      # menuItem(text = "step34_smORF_mORF", tabName = "step34_smORF_mORF", icon = icon("chart-bar"),
      #          menuSubItem(text = "Association", tabName = "step34_smORF_mORF_association")
      # )
    )
  ),
  
  
  ############################ dashboardBody ############################
  
  dashboardBody(
    
    tabItems(
      
      ############################################################################################################
      ## figure overview ###############################
      ############################################################################################################
      
      tabItem(tabName = "OverView",
              h2("A comprehensive framework for RNA-seq and Ribo-seq data mining and visualization."),
              
              column(10,
                     tabBox(id = "tabset1", height = "800px", width = NULL, title = "Figure Overview",
                            tabPanel(title = "Frame Work", icon = icon("chart-bar"),
                                     h4("Show the flowchart of RiboParser and RiboShiny."),
                                     imageOutput("Plot_fig_framework")),
                            
                            tabPanel(title = "Offset detection", icon = icon("chart-bar"),
                                     h4("Offset detection model."),
                                     imageOutput("Plot_fig_offset_detection")),
                            
                            tabPanel(title = "Quality Control", icon = icon("chart-bar"),
                                     h4("Quality control and visualization."),
                                     imageOutput("Plot_fig_quality_control")),
                            
                            tabPanel(title = "Codon Level", icon = icon("chart-bar"),
                                     h4("Codon level analysis and visualization."),
                                     imageOutput("Plot_fig_codon_level")),
                            
                            tabPanel(title = "Gene Level", icon = icon("chart-bar"),
                                     h4("Gene level analysis and visualization."),
                                     imageOutput("Plot_fig_gene_level")),
                            
                            tabPanel(title = "SeRP", icon = icon("chart-bar"),
                                     h4("SeRP analysis and visualization."),
                                     imageOutput("Plot_fig_serp")),
                            
                            tabPanel(title = "smORF", icon = icon("chart-bar"),
                                     h4("To be continue."),
                                     imageOutput("Plot_fig_smorf"))
                     )
              )
              
              # box(title = "Plot_fig_framework", collapsible = TRUE, collapsed = FALSE,
              #     width = NULL, status = "primary", solidHeader = TRUE,
              #     imageOutput("Plot_fig_framework")
              # ),
              # 
              # box(title = "Plot_fig_offset_detection", collapsible = TRUE, collapsed = FALSE,
              #     width = NULL, status = "primary", solidHeader = TRUE,
              #     imageOutput("Plot_fig_offset_detection")
              # ),
              # 
              # box(title = "Plot_fig_genelevel", collapsible = TRUE, collapsed = FALSE,
              #     width = NULL, status = "primary", solidHeader = TRUE,
              #     imageOutput("Plot_fig_genelevel")
              # ),
              # 
              # box(title = "Plot_fig_serp", collapsible = TRUE, collapsed = FALSE,
              #     width = NULL, status = "primary", solidHeader = TRUE,
              #     imageOutput("Plot_fig_serp")
              # ),
              # 
              # box(title = "Plot_fig_smorf", collapsible = TRUE, collapsed = FALSE,
              #     width = NULL, status = "primary", solidHeader = TRUE,
              #     imageOutput("Plot_fig_smorf")
              # )
              
              # fluidRow(
              #   # imageOutput("Plot_fig_framework"),
              #   # imageOutput("Plot_fig_offset_detection"),
              #   # imageOutput("Plot_fig_genelevel"),
              #   imageOutput("Plot_fig_serp"),
              #   imageOutput("Plot_fig_smorf")
              # )
      ),
      
      
      ############################################################################################################
      # tabItem(tabName = "step1_ip",
      #         h2("A comprehensive framework for RNA-seq and Ribo-seq data mining and visualization."),
      #         
      #         column(7,
      #                tabsetPanel(
      #                  tabPanel(title = "ip_address",
      #                           h4("Show the ip_address."),
      #                           textOutput("ip_address"))
      #                )
      #         )
      # ),

      
      ############################################################################################################
      ## step1 set the design and database ###############################
      
      ### step1 options ###############################
      tabItem(tabName = "step1_options",
              
              h2("Set the options."),
              fluidRow(
                
                column(5,
                       box(title = "Set the design", collapsible = TRUE, collapsed = FALSE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           
                           # set the input panel
                           fluidRow(
                             column(8, fileInput("in_step1_design", label = "experiment design excel file", 
                                                 multiple = FALSE, accept = c(".xlsx", ".XLSX"))),
                             column(4, numericInput("in_step1_design_sheet", label = "sheet number", value = 1, min = 1))
                           ),
                           p("The design file should be a excel file contains four columns: ['Sample', 'SeqType', 'Rank', 'Group']."),
                           
                           ## click to import the design file
                           fluidRow(
                             column(6, actionButton("act_step1_import_design", label = "import the design", icon = icon("file"))),
                             column(6, actionButton("act_step1_show_example", label = "show the example", icon = icon("file")))
                           ),
                           
                           tags$style(HTML("#act_step1_import_design {background-color: #75aadb; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#act_step1_show_example {background-color: #75aadb; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Set the output", collapsible = TRUE, collapsed = FALSE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           
                           ## click to check and create the output folder
                           textInput("output_folder", label = "create output folder", value = "C:/Users/Default/Desktop/RiboParser/test/output"),
                           
                           actionButton("act_create_folder", label = "create the folder", icon = icon("folder")),
                           tags$style(HTML("#act_create_folder {background-color: #75aadb; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                ## set the output panel
                column(7,
                       tabsetPanel(
                         tabPanel(title = "Design",
                                  h4("Show the input design file."),
                                  DT::dataTableOutput("out_step1_design")),
                         
                         tabPanel(title = "Example",
                                  h4("Show the example of design file."),
                                  DT::dataTableOutput("out_step1_example")),
                         
                         tabPanel(title = "Folder",
                                  h4("Show the output folder file."),
                                  tableOutput("out_step1_folder"))
                         
                       )
                )
              )
      ),
      
      ### step1 annotation ###############################
      tabItem(tabName = "step1_annotation",
              h2("Check the gene annotation."),
              fluidRow(
                tags$style(HTML(".download-btn {margin-top: 300px;}")),
                
                column(5,
                       
                       box(title = "Import gene annotation.", collapsed = FALSE,
                           collapsible = TRUE, width = NULL, status = "primary", solidHeader = TRUE,
                           
                           # set the input panel
                           fluidRow(
                             column(8, fileInput(inputId = "in_step1_anno", label = "annotation excel file", accept = c(".xlsx", ".XLSX"))),
                             column(4, numericInput("in_step1_anno_sheet", label = "sheet number", value = 1, min = 1))
                           ),
                           
                           ## click to import the annotation
                           actionButton("act_step1_import_anno", label = "import the annotation", icon = icon("file")),
                           tags$style(HTML("#act_step1_import_anno {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Import gene annotation.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           fluidRow(
                             column(4, checkboxInput("in_step1_anno_longest", label = "longest", value = TRUE)),
                             column(4, checkboxInput("in_step1_sqrty", label = "sqrt_y", value = FALSE)),
                             column(4, selectizeInput("in_step1_anno_type", label = "type", selected = "cds_length", 
                                                      choices = c("cds_length", "utr5_length", "utr3_length")))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step1_edge_color", label = "edge color", value = "#2244aa")),
                             column(6, textInput("in_step1_fill_color", label = "fill color", value = "#2244aa"))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step1_font_size", label = "Font size", value = 15, min = 5, max = 25,  step = 0.5)),
                             column(6, sliderInput("in_step1_fill_alpha", label = "Alpha", min = 0, max = 1, value = 0.9, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step1_length_xlim", label = "xlim", min = 0, max = 10000, value = c(0, 10000))),
                             column(6, sliderInput("in_step1_bin_range", label = "bin width", min = 1, max = 500, value = 10))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step1_hist_width", label = "Figure width", value = 5, min = 1)),
                             column(6, numericInput("out_step1_hist_height", label = "Figure height", value = 4, min = 1))
                           ),

                           fluidRow(
                             column(6, actionButton("act_step1_draw_histplot", label = "draw the histplot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step1_hist", label = "Save hist plot",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step1_draw_histplot {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step1_hist {background-color: #e78c45; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Preview",
                                  h4("Show the annotation file."),
                                  conditionalPanel(condition = ("input.act_step1_import_anno > 0"),
                                                   DT::dataTableOutput("out_step1_annotation") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         ## set the design panel
                         tabPanel(title = "length",
                                  h4("Show the hist plot file."),
                                  plotOutput("out_step1_hist_plot"))
                         
                       )
                )
              )
      ),
      
      ### step1 OrgDb #####################################
      tabItem(tabName = "step1_orgdb",
              h2("Create the GO orgdb"),
              fluidRow(
                tags$style(HTML(".download-btn {margin-top: 300px;}")),
                
                column(5,
                       # h4("The orgdb is used to convert the gene id to gene symbol."),
                       box(title = "Retrieve the orgdb.", collapsible = TRUE, collapsed = FALSE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           
                           ## click to import the go orgdb
                           checkboxInput("in_step1_orgdb_local", label = "Local", value = FALSE),
                           actionButton("act_step1_retrieve_orgdb", label = "Retrieve OrgDb", icon = icon("file")),
                           tags$style(HTML("#act_step1_retrieve_orgdb {background-color: #75aadb; color: black;}")),
                           
                           # set the input panel
                           fluidRow(
                             column(4, textInput(inputId = "in_step1_orgdb_genus", label = "Genus", value = "", placeholder = "e.g. homo")),
                             column(4, textInput(inputId = "in_step1_orgdb_species", label = "Species", value = "", placeholder = "e.g. sapiens")),
                             column(4, textInput(inputId = "in_step1_orgdb_taxonomy", label = "Taxonomy", value = "", placeholder = "e.g. 9606"))
                           ),
                           
                           actionButton("act_step1_filter_orgdb", label = "Filter OrgDb", icon = icon("file")),
                           tags$style(HTML("#act_step1_filter_orgdb {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Create the orgdb.", collapsed = FALSE, collapsible = TRUE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           # create the orgdb
                           textInput("out_step1_orgdb_dir", label = "Output directory", value = "D:/database/orgdb/", placeholder = "D:/database/orgdb/"),
                           
                           textInput("in_step1_orgdb_records", label = "Records", value = "", placeholder = "e.g. AH114084"),
                           
                           selectizeInput("in_step1_orgdb_column", label = "Annotation column", 
                                          selected = c("GID", "SYMBOL", "GENENAME", "GO", "ONTOLOGY", "EVIDENCE"),
                                          choices = c("CHR", "GID", "ENTREZID", "ALIAS", "ACCNUM", "SYMBOL", "REFSEQ", 
                                                      "GENENAME", "GO", "ONTOLOGY", "EVIDENCE", "PMID"),
                                          multiple = TRUE),
                           
                           actionButton("act_step1_create_orgdb", label = "Create OrgDb", icon = icon("file")),
                           tags$style(HTML("#act_step1_create_orgdb {background-color: #75aadb; color: black;}")),
                           
                           fluidRow(
                             column(4, downloadButton("save_step1_gson_go", label = "Save gson GO",
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step1_gene_mess", label = "Save gene mess",
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step1_orgdb", label = "Save OrgDb",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step1_gson_go {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step1_gene_mess {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step1_orgdb {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                       box(title = "Install the orgdb.", collapsed = FALSE, collapsible = TRUE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           # install the orgdb
                           fileInput("in_step1_orgdb_file", label = "OrgDb", accept = ".gz", multiple = FALSE),
                           
                           actionButton("act_step1_install_orgdb", label = "Install OrgDb", icon = icon("file")),
                           tags$style(HTML("#act_step1_install_orgdb {background-color: #75aadb; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "AnnotationHub",
                                  h4("Show the annotationhub file."),
                                  conditionalPanel(condition = ("input.act_step1_retrieve_orgdb > 0"),
                                                   DT::dataTableOutput("out_step1_annotationhub") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Preview",
                                  h4("Show the orgdb file."),
                                  DT::dataTableOutput("out_step1_orgdb_preview")),
                         
                         ## set the design panel
                         tabPanel(title = "OrgDb",
                                  h4("Show the OrgDb gene message."),
                                  conditionalPanel(condition = ("input.act_step1_create_orgdb > 0"),
                                                   DT::dataTableOutput("out_step1_orgdb_mess") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Install",
                                  h4("Install the OrgDb."),
                                  conditionalPanel(condition = ("input.act_step1_install_orgdb > 0"),
                                                   verbatimTextOutput("out_step1_orgdb_install") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                         
                       )
                )
              )
      ),
      
      ### step1 KEGG ###################################
      tabItem(tabName = "step1_kegg",
              h2("Create the gson KEGG"),
              fluidRow(
                tags$style(HTML(".download-btn {margin-top: 300px;}")),
                
                column(5,
                       box(title = "Set the species.", collapsible = TRUE,
                           collapsed = FALSE, width = NULL, status = "primary", solidHeader = TRUE,
                           
                           # set the input panel
                           textInput(inputId = "in_step1_kegg_species", label = "Species", value = "hsa", placeholder = "e.g. hsa"),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step1_kegg_type", label = "KEGG type", selected = "KEGG",
                                                      choices = c("KEGG", "MKEGG"))),
                             column(6, selectizeInput("in_step1_kegg_key_type", label = "KeyType", selected = "kegg",
                                                      choices = c("kegg", "ncbi-geneid", "ncbi-proteinid", "uniprot")))
                           ),
                           
                           textInput("out_step1_kegg_name", label = "Output", value = "", placeholder = "e.g. hsa"),
                           
                           fluidRow(
                             column(6, actionButton("act_step1_retrieve_kegg", label = "Retrieve KEGG", icon = icon("file"))),
                             column(6, downloadButton("save_step1_gson_kegg", label = "Save gson kegg", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step1_retrieve_kegg {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step1_gson_kegg {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                column(5,
                       tabsetPanel(
                         
                         ## set the design panel
                         tabPanel(title = "KEGG",
                                  h4("Show the kegg gene message."),
                                  conditionalPanel(condition = ("input.act_step1_retrieve_kegg > 0"),
                                                   DT::dataTableOutput("out_step1_kegg_mess") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                         
                       )
                )
              )
      ),
      
      ### step1 sequence ###############################
      tabItem(tabName = "step1_sequence",
              h2("Codon usage."),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center;}")),
                column(5,
                       box(title = "Import gene sequence.", collapsible = TRUE, 
                           collapsed = FALSE, width = NULL, status = "primary", solidHeader = TRUE,
                           
                           ## click to import the sequence
                           fileInput(inputId = "in_step1_seq", label = "mRNA sequence", accept = c(".fa", ".fasta", ".FA", ".FASTA")),
                           
                           actionButton("act_step1_import_seq", label = "import the mRNA sequence", icon = icon("file")),
                           tags$style(HTML("#act_step1_import_seq {background-color: #75aadb; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Calculate the codon usage.",
                           collapsible = TRUE, width = NULL, status = "primary", solidHeader = TRUE,
                           ## click to calculate the codon usage
                           p("Recommend using the most representative transcript in each gene to calculate the codon usage."),
                           
                           fluidRow(
                             column(6, checkboxInput("in_step1_stop", label = "Stop codon", value = TRUE)),
                             column(6, numericInput("in_step1_float", label = "decimal precision", value = 4, min = 0))
                           ),
                           
                           textInput("out_step1_codon_count_name", label = "Codon table name", value = "CDS"),
                           
                           fluidRow(
                             column(6, actionButton("act_step1_calc_cu", label = "calculate codon usage", icon = icon("calculator"))),
                             column(6, downloadButton("save_step1_codon_table", label = "save codon table",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step1_calc_cu {background-color: #75aadb; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step1_codon_table {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                       box(title = "Draw codon usage.", collapsible = TRUE, 
                           collapsed = FALSE, width = NULL, status = "primary", solidHeader = TRUE,
                           ## click to draw the codon usage
                           tags$hr(),
                           div(style = "border-bottom: 2px solid #0099ff; margin-bottom: 5px;"),
                           fluidRow(
                             column(6, selectizeInput("in_step1_cu_plot_type", label = "Type", selected = "dot",
                                                      choices = c("dot", "circle"))),
                             column(6, checkboxInput("in_step1_cu_wrap", label = "Wrap plot", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step1_cu_wrap_row", label = "Warp rows", value = 3, min = 0)),
                             column(6, numericInput("in_step1_cu_legend_row", label = "Legend rows", value = 3, min = 0))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step1_cu_grid_width", label = "Grid width", value = 0.2, min = 0, step = 0.05)),
                             column(6, textInput("in_step1_cu_grid_color", label = "Grid color", value = "grey85"))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step1_cu_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(6, numericInput("in_step1_cu_dot_size", label = "Dot size", value = 3, min = 0, max = 10, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step1_cu_color_map", label = "Color", selected = "Paired", 
                                                      choices = color_list, multiple = FALSE)),
                             column(4, sliderInput("in_step1_cu_color_alpha", label = "Color alpha", value = 0.8, min = 0, max = 1, step = 0.05)),
                             column(4, sliderInput("in_step1_cu_line_width", label = "Line width", value = 0.8, min = 0, max = 5, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step1_cu_fig_name", label = "Figure name", value = "Codon")),
                             column(4, numericInput("out_step1_cu_fig_width", label = "Figure width", value = 9, min = 1)),
                             column(4, numericInput("out_step1_cu_fig_height", label = "Figure height", value = 6, min = 1))
                           ),
                           
                           actionButton("act_step1_draw_cu", label = "draw codon usage", icon = icon("chart-line")),
                           tags$style(HTML("#act_step1_draw_cu {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           
                           fluidRow(
                             column(4, downloadButton("save_step1_cu_freq", label = "Save freq", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step1_cu_rscu", label = "Save rscu", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step1_cu_cai", label = "Save cai", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step1_cu_freq {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step1_cu_rscu {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step1_cu_cai {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Preview",
                                  h4("Show the annotation file."),
                                  verbatimTextOutput("out_step1_sequence")),
                         
                         ## set the codon table
                         tabPanel(title = "Codon table",
                                  h4("Show the codon table."),
                                  conditionalPanel(condition = ("input.act_step1_calc_cu > 0"),
                                                   DT::dataTableOutput("out_step1_codon_table") %>% 
                                                     withSpinner(color = "#3c8dbc", type = 5, size = 1))),
                         
                         ## set the codon usage panel
                         tabPanel(title = "Codon usage",
                                  h4("Show the codon usage table."),
                                  conditionalPanel(condition = ("input.act_step1_calc_cu > 0"),
                                                   DT::dataTableOutput("out_step1_cu_table") %>% 
                                                     withSpinner(color = "#3c8dbc", type = 5, size = 1))),
                         
                         ## set the codon usage panel
                         tabPanel(title = "Frequency",
                                  h4("Show the codon frequency plot."),
                                  plotOutput("out_step1_freq_plot")),
                         
                         tabPanel(title = "RSCU",
                                  h4("Show the codon RSCU plot."),
                                  plotOutput("out_step1_rscu_plot")),
                         
                         tabPanel(title = "CAI",
                                  h4("Show the codon CAI plot."),
                                  plotOutput("out_step1_cai_plot"))
                         
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step2 mapping condition ###############################
      
      ### step2 mapping table ###############################
      tabItem(tabName = "step2_riboseq",
              # h2("draw the RNA-seq mapping plot"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center;}")),
                # set the input panel
                column(5,
                       box(title = "Import the alignment.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           ## click to import the sequence
                           fileInput(inputId = "in_step2_align_file", label = "alignments", 
                                     multiple = FALSE, accept = c(".txt", ".TXT")),
                           
                           actionButton("act_step2_import_align", label = "import alignment", icon = icon("file")),
                           tags$style(HTML("#act_step2_import_align {background-color: #75aadb; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Draw the RNA-seq alignment plot.", collapsed = FALSE, collapsible = TRUE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           # click to draw the alignment plot
                           textInput("in_step2_x", label = "X-axis", value = "Sample", placeholder = "need design file"),
                           
                           fluidRow(
                             column(6, checkboxInput("in_step2_coord", label = "Coord flip", value = FALSE)),
                             column(6, checkboxInput("in_step2_label", label = "Label", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step2_edge_color", label = "Edge", value = "NA", placeholder = "e.g. grey85")),
                             column(6, selectizeInput("in_step2_fill_color", label = "Fill", choices = color_list, selected = "Set1"))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step2_fill_alpha", label = "Color alpha", value = 0.9, min = 0, max = 1, step = 0.1)),
                             column(6, sliderInput("in_step2_bar_width", label = "Bar width", value = 0.8, min = 0, max = 1, step = 0.1)),
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step2_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(6, sliderInput("in_step2_label_size", label = "Label size", value = 3, min = 0, max = 10, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step2_fig_width", label = "Figure width", value = 6, min = 1)),
                             column(6, numericInput("out_step2_fig_height", label = "Figure height", value = 6, min = 1))
                           ),
                           
                           textInput("out_step2_fig_name", label = "Output", value = "Alignment", placeholder = "file name prefix"),
                           
                           fluidRow(
                             column(4, actionButton("act_step2_draw_align", label = "draw alignment", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step2_count_plot", label = "Save count", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step2_ratio_plot", label = "Save ratio", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step2_draw_align {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step2_count_plot {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step2_ratio_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Alignment",
                                  h4("Show the alignment file"),
                                  DT::dataTableOutput("out_step2_align_table")),
                         
                         ## set the alignment table
                         tabPanel(title = "Count",
                                  h4("Show the alignment plot"),
                                  plotOutput("out_step2_count_plot")),
                         
                         tabPanel(title = "Ratio",
                                  h4("Show the alignment plot"),
                                  plotOutput("out_step2_ratio_plot"))
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step3 length distribution ###############################
      
      ### step3 RNA-seq distribution ###############################
      tabItem(tabName = "step3_rnaseq",
              h2("draw the RNA-seq length distribution"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center;}")),
                # set the input panel
                column(5,
                       box(title = "Import the RNA-seq length.", collapsed = FALSE, collapsible = TRUE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           ## click to import the reads length distribution
                           fileInput(inputId = "in_step3_rnaseq", label = "length table", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           textInput("in_step3_rna_x", label = "X-axis", value = "Sample", placeholder = "need design file"),
                           
                           actionButton("act_step3_import_rna_length", label = "import the distribution", icon = icon("file")),
                           tags$style(HTML("#act_step3_import_rna_length {background-color: #75aadb; color: black; margin-top: 0px;}")),
                       ),
                       
                       box(title = "Draw the RNA-seq distribution plot.", collapsed = FALSE, collapsible = TRUE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           fluidRow(
                             column(6, selectizeInput("in_step3_rna_strand", label = "Strand", choices = c("Plus", "Minus"), selected = "Plus")),
                             column(6, selectizeInput("in_step3_rna_type", label = "Type", choices = c("Count", "Ratio"), selected = "Count"))
                           ),
                           
                           fluidRow(
                             column(6, checkboxInput("in_step3_rna_wrap", label = "Wrap", value = FALSE)),
                             column(6, sliderInput("in_step3_rna_xlim", label = "Xlim", value = c(20, 151), min = 1, max = 151, dragRange = TRUE))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step3_rna_line_color", label = "Line color", choices = color_list, selected = "Paired")),
                             column(6, sliderInput("in_step3_rna_line_width", label = "Line width", value = 1, min = 0, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step3_rna_fill_color", label = "Fill color", choices = color_list, selected = "Blues")),
                             column(6, sliderInput("in_step3_rna_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step3_rna_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(6, sliderInput("in_step3_rna_dot_size", label = "Dot size", value = 1, min = 0, max = 10, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step3_rna_width", label = "Figure width", value = 9, min = 1)),
                             column(6, numericInput("out_step3_rna_height", label = "Figure height", value = 6, min = 1))
                           ),
                           
                           textInput("out_step3_rna_fig_name", label = "Output", value = "RNA-seq", placeholder = "file name prefix"),
                           
                           fluidRow(
                             column(4, actionButton("act_step3_draw_rna_length", label = "draw the distribution", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step3_rna_line", label = "Save lineplot", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step3_rna_heat", label = "Save heatmap", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step3_draw_rna_length {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step3_rna_line {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step3_rna_heat {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Distribution",
                                  h4("Show the distribution table"),
                                  DT::dataTableOutput("out_step3_rna_distr")),
                         
                         ## set the length distribution plot
                         tabPanel(title = "Line",
                                  h4("Show the RNA-seq line plot"),
                                  plotOutput("out_step3_rna_line")),
                         
                         tabPanel(title = "Heat",
                                  h4("Show the RNA-seq heat plot"),
                                  plotOutput("out_step3_rna_heat"))
                       )
                )
              )
      ),
      
      ### step3 Ribo-seq distribution ###############################
      tabItem(tabName = "step3_riboseq",
              h2("draw the Ribo-seq length distribution"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center;}")),
                
                ## set the input panel
                column(5,
                       box(title = "Import the Ribo-seq length.", collapsed = FALSE, collapsible = TRUE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           ## click to import the reads length distribution
                           fileInput(inputId = "in_step3_riboseq", label = "length table", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           textInput("in_step3_ribo_x", label = "X-axis", value = "Sample", placeholder = "need design file"),
                           
                           actionButton("act_step3_import_ribo_length", label = "import the distribution", icon = icon("file")),
                           tags$style(HTML("#act_step3_import_ribo_length {background-color: #75aadb; color: black; margin-top: 0px;}")),
                       ),
                       
                       box(title = "Import the Ribo-seq distribution plot.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           fluidRow(
                             column(6, selectizeInput("in_step3_ribo_strand", label = "Strand", choices = c("Plus", "Minus"), selected = "Plus")),
                             column(6, selectizeInput("in_step3_ribo_type", label = "Type", choices = c("Count", "Ratio"), selected = "Count"))
                           ),
                           
                           fluidRow(
                             column(6, checkboxInput("in_step3_ribo_wrap", label = "Wrap", value = FALSE)),
                             column(6, sliderInput("in_step3_ribo_xlim", label = "Xlim", value = c(20, 40), min = 1, max = 60, dragRange = TRUE))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step3_ribo_line_color", label = "Line color", choices = color_list, selected = "Paired")),
                             column(6, sliderInput("in_step3_ribo_line_width", label = "Line width", value = 1, min = 0, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step3_ribo_fill_color", label = "Fill color", choices = color_list, selected = "Blues")),
                             column(6, sliderInput("in_step3_ribo_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step3_ribo_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(6, sliderInput("in_step3_ribo_dot_size", label = "Dot size", value = 1, min = 0, max = 10, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step3_ribo_width", label = "Figure width", value = 7, min = 1)),
                             column(6, numericInput("out_step3_ribo_height", label = "Figure height", value = 5, min = 1))
                           ),
                           
                           textInput("out_step3_ribo_fig_name", label = "Output", value = "Ribo-seq", placeholder = "file name prefix"),
                           
                           fluidRow(
                             column(4, actionButton("act_step3_draw_ribo_length", label = "draw the distribution", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step3_ribo_line", label = "Save lineplot", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step3_ribo_heat", label = "Save heatmap", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step3_draw_ribo_length {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step3_ribo_line {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step3_ribo_heat {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Distribution",
                                  h4("Show the distribution table"),
                                  DT::dataTableOutput("out_step3_ribo_distr")),
                         
                         ## set the length distribution plot
                         tabPanel(title = "Line",
                                  h4("Show the RNA-seq line plot"),
                                  plotOutput("out_step3_ribo_line")),
                         
                         tabPanel(title = "Heat",
                                  h4("Show the RNA-seq heat plot"),
                                  plotOutput("out_step3_ribo_heat"))
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step4 gene saturation ###############################
      
      ### step4 RNA-seq saturation ###############################
      tabItem(tabName = "step4_rnaseq",
              h2("draw the RNA-seq saturation"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center;}")),
                ## set the input panel
                column(5,
                       box(title = "Import the RNA-seq saturation.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           ## click to import the reads length distribution
                           fileInput(inputId = "in_step4_rnaseq", label = "Satutation", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           actionButton("act_step4_import_rna_saturation", label = "import the saturation", icon = icon("file")),
                           tags$style(HTML("#act_step4_import_rna_saturation {background-color: #75aadb; color: black; margin-top: 0px;}")),
                       ),
                       
                       box(title = "Draw the RNA-seq saturation.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           fluidRow(
                             column(6, textInput("in_step4_rna_x", label = "X-axis", value = "Sample", placeholder = "need design file")),
                             column(6, checkboxInput("in_step4_rna_label", label = "Total", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step4_rna_dot_size", label = "Dot size", value = 3, min = 1, max = 6, step = 0.1)),
                             column(6, sliderInput("in_step4_rna_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step4_rna_line_color", label = "Line color", choices = color_list, selected = "Paired")),
                             column(6, sliderInput("in_step4_rna_line_width", label = "Line width", value = 1, min = 0, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step4_rna_fill_color", label = "Fill color", choices = color_list, selected = "Blues")),
                             column(6, sliderInput("in_step4_rna_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step4_rna_width", label = "Figure width", value = 10, min = 1)),
                             column(6, numericInput("out_step4_rna_height", label = "Figure height", value = 4.5, min = 1))
                           ),
                           
                           textInput("out_step4_rna_fig_name", label = "Output", value = "RNA-seq", placeholder = "file name prefix"),
                           
                           fluidRow(
                             column(4, actionButton("act_step4_draw_rna_saturation", label = "draw the saturation", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step4_rna_line", label = "Save lineplot", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step4_rna_heat", label = "Save heatmap", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step4_draw_rna_saturation {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step4_rna_line {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step4_rna_heat {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Saturation",
                                  h4("Show the saturation table"),
                                  DT::dataTableOutput("out_step4_rna_saturation")),
                         
                         ## set the length saturation plot
                         tabPanel(title = "Line",
                                  h4("Show the RNA-seq saturation line plot"),
                                  plotOutput("out_step4_rna_line")),
                         
                         tabPanel(title = "Heat",
                                  h4("Show the RNA-seq saturation heat plot"),
                                  plotOutput("out_step4_rna_heat"))
                       )
                )
              )
      ),
      
      ### step4 Ribo-seq saturation ###############################
      tabItem(tabName = "step4_riboseq",
              h2("draw the Ribo-seq saturation"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center;}")),
                
                ## set the input panel
                column(5,
                       box(title = "Import the RNA-seq saturation.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           ## click to import the reads length distribution
                           fileInput(inputId = "in_step4_riboseq", label = "Satutation", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           actionButton("act_step4_import_ribo_saturation", label = "import the saturation", icon = icon("file")),
                           tags$style(HTML("#act_step4_import_ribo_saturation {background-color: #75aadb; color: black; margin-top: 0px;}")),
                       ),
                       
                       box(title = "Import the RNA-seq saturation.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           fluidRow(
                             column(6, textInput("in_step4_ribo_x", label = "X-axis", value = "Sample", placeholder = "need design file")),
                             column(6, checkboxInput("in_step4_ribo_label", label = "Total", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step4_ribo_dot", label = "Dot size", value = 3, min = 1, max = 6, step = 0.1)),
                             column(6, sliderInput("in_step4_ribo_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step4_ribo_line_color", label = "Line color", choices = color_list, selected = "Paired")),
                             column(6, sliderInput("in_step4_ribo_line_width", label = "Line width", value = 1, min = 0, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step4_ribo_fill_color", label = "Fill color", choices = color_list, selected = "Blues")),
                             column(6, sliderInput("in_step4_ribo_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step4_ribo_width", label = "Figure width", value = 10, min = 1)),
                             column(6, numericInput("out_step4_ribo_height", label = "Figure height", value = 5, min = 1))
                           ),
                           
                           textInput("out_step4_ribo_fig_name", label = "Output", value = "Ribo-seq", placeholder = "file name prefix"),
                           
                           fluidRow(
                             column(4, actionButton("act_step4_draw_ribo_saturation", label = "draw the saturation", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step4_ribo_line", label = "Save line", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step4_ribo_heat", label = "Save heatmap", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step4_draw_ribo_saturation {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step4_ribo_line {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step4_ribo_heat {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Saturation",
                                  h4("Show the saturation table"),
                                  DT::dataTableOutput("out_step4_ribo_saturation")),
                         
                         ## set the length saturation plot
                         tabPanel(title = "Line",
                                  h4("Show the Ribo-seq saturation line plot"),
                                  plotOutput("out_step4_ribo_line")),
                         
                         tabPanel(title = "Heat",
                                  h4("Show the Ribo-seq saturation heat plot"),
                                  plotOutput("out_step4_ribo_heat"))
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step5 reads digestion ###############################
      
      ### step5 Ribo-seq digestion ###############################
      tabItem(tabName = "step5_riboseq",
              h2("draw the Ribo-seq digestion"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                ## set the input panel
                column(5,
                       box(title = "Import the digestion.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           ## click to import the digestion
                           fileInput(inputId = "in_step5_riboseq", label = "Digestion", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           actionButton("act_step5_import_ribo_digestion", label = "import the digestion", icon = icon("file")),
                           tags$style(HTML("#act_step5_import_ribo_digestion {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Draw the digestion.", collapsed = FALSE, collapsible = TRUE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           
                           fluidRow(
                             column(6, textInput("in_step5_ribo_x", label = "X-axis", value = "Sample", placeholder = "need design file")),
                             column(6, selectizeInput("in_step5_ribo_method", label = "Method", 
                                                      choices = c('bits', 'probability', 'custom'), selected = "probability"))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step5_ribo_fill_color", label = "Fill color", selected = "auto", choices = seqlogo_colors)),
                             column(6, sliderInput("in_step5_ribo_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step5_ribo_font_family", label = "Font family", 
                                                      choices = seqlogo_font_family, selected = "helvetica_regular")),
                             column(4, sliderInput("in_step5_ribo_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(4, sliderInput("in_step5_ribo_font_stack", label = "Stack width", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step5_ribo_fig_name", label = "Output", value = "Ribo-seq", placeholder = "file name prefix")),
                             column(4, numericInput("out_step5_ribo_width", label = "Figure width", value = 10, min = 1)),
                             column(4, numericInput("out_step5_ribo_height", label = "Figure height", value = 8, min = 1))
                           ),
                           
                           
                           
                           fluidRow(
                             column(4, actionButton("act_step5_draw_ribo_digestion", label = "draw the digestion", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step5_ribo_seqlogo_end_5p", label = "Save end5p", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step5_ribo_seqlogo_end_3p", label = "Save end3p", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step5_draw_ribo_digestion {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step5_ribo_seqlogo_end_5p {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step5_ribo_seqlogo_end_3p {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Digestion",
                                  h4("Show the digestion table"),
                                  DT::dataTableOutput("out_step5_ribo_digestion")),
                         
                         ## set the length digestion plot
                         tabPanel(title = "SeqLogo_End_5p",
                                  h4("Show the Ribo-seq digestion line plot"),
                                  plotOutput("out_step5_ribo_seqlogo_end_5p")),
                         
                         tabPanel(title = "SeqLogo_End_3p",
                                  h4("Show the Ribo-seq digestion line plot"),
                                  plotOutput("out_step5_ribo_seqlogo_end_3p"))
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step6 Ribo-seq offset ###############################
      
      ### step6 Ribo-seq RSBM offset ###############################
      tabItem(tabName = "step6_frame",
              h2("draw the offset of RSBM"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                ## set the input panel
                column(5,
                       box(title = "Import the RSBM offset.", collapsed = FALSE, collapsible = TRUE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           ## click to import the RSBM
                           fileInput(inputId = "in_step6_frame", label = "Offset", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           actionButton("act_step6_import_frame_offset", label = "import the offset", icon = icon("file")),
                           tags$style(HTML("#act_step6_import_frame_offset {background-color: #75aadb; color: black; margin-top: 0px;}")),
                       ),
                       
                       box(title = "Draw the frame offset.", collapsed = FALSE, collapsible = TRUE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           ## set the input panel
                           textInput("in_step6_frame_x", label = "X-axis", value = "Sample", placeholder = "need design file"),
                           fluidRow(
                             column(6, selectizeInput("in_step6_frame_value", label = "Value", choices = c("Count", "Ratio"), selected = "Count")),
                             column(6, checkboxInput("in_step6_frame_warp", label = "Warp", value = TRUE))
                           ),
                           
                           sliderInput("in_step6_frame_xlim", label = "Xlim", value = c(27, 33), min = 20, max = 40, step = 1),
                           
                           fluidRow(
                             column(6, sliderInput("in_step6_frame_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(6, sliderInput("in_step6_frame_bar_width", label = "Bar width", value = 0.8, min = 0.1, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step6_frame_fill_color", label = "Fill color", selected = "Blues", choices = color_list)),
                             column(6, sliderInput("in_step6_frame_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step6_frame_width", label = "Figure width", value = 10, min = 1)),
                             column(6, numericInput("out_step6_frame_height", label = "Figure height", value = 5, min = 1))
                           ),
                           
                           textInput("out_step6_frame_fig_name", label = "Output", value = "Ribo-seq", placeholder = "file name prefix"),
                           
                           fluidRow(
                             column(4, actionButton("act_step6_draw_frame_offset", label = "draw the offset", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step6_frame_bar", label = "Save barplot", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step6_frame_heat", label = "Save heatmap", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step6_draw_frame_offset {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step6_frame_bar {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step6_frame_heat {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Offset",
                                  h4("Show the RSBM offset table"),
                                  DT::dataTableOutput("out_step6_frame_offset")),
                         
                         ## set the length frame offset plot
                         tabPanel(title = "Bar",
                                  h4("Show the Ribo-seq RSBM offset bar plot"),
                                  plotOutput("out_step6_frame_bar")),
                         
                         tabPanel(title = "Heat",
                                  h4("Show the Ribo-seq RSBM offset heat plot"),
                                  plotOutput("out_step6_frame_heat"))
                       )
                )
              )
      ),
      
      ### step6 Ribo-seq SSCBM offset ###############################
      tabItem(tabName = "step6_tis",
              h2("draw the offset of SSCBM"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                ## set the input panel
                column(5,
                       box(title = "Import the SSCBM offset.", collapsed = FALSE, collapsible = TRUE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           ## click to import the SSCBM
                           fileInput(inputId = "in_step6_tis", label = "Offset", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           actionButton("act_step6_import_tis_offset", label = "import the offset", icon = icon("file")),
                           tags$style(HTML("#act_step6_import_tis_offset {background-color: #75aadb; color: black;}")),
                       ),
                       
                       box(title = "Draw the tis offset.", collapsed = FALSE, collapsible = TRUE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           textInput("in_step6_tis_x", label = "X-axis", value = "Sample", placeholder = "need design file"),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step6_tis_value", label = "Value", choices = c("Count", "Ratio"), selected = "Count")),
                             column(6, checkboxInput("in_step6_tis_warp", label = "Warp", value = TRUE))
                           ),
                           
                           sliderInput("in_step6_tis_xlim", label = "Xlim", value = c(27, 33), min = 20, max = 40, step = 1),
                           
                           fluidRow(
                             column(6, sliderInput("in_step6_tis_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(6, sliderInput("in_step6_tis_bar_width", label = "Bar width", value = 0.8, min = 0.1, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step6_tis_fill_color", label = "Fill color", selected = "Blues", choices = color_list)),
                             column(6, sliderInput("in_step6_tis_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step6_tis_width", label = "Figure width", value = 10, min = 1)),
                             column(6, numericInput("out_step6_tis_height", label = "Figure height", value = 5, min = 1))
                           ),
                           
                           textInput("out_step6_tis_fig_name", label = "Output", value = "Ribo-seq", placeholder = "file name prefix"),
                           
                           fluidRow(
                             column(4, actionButton("act_step6_draw_tis_offset", label = "draw the offset", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step6_tis_bar", label = "Save barplot", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step6_tis_heat", label = "Save heatmap", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step6_draw_tis_offset {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step6_tis_bar {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step6_tis_heat {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Offset",
                                  h4("Show the SSCBM offset table"),
                                  DT::dataTableOutput("out_step6_tis_offset")),
                         
                         ## set the tis offset plot
                         tabPanel(title = "Bar",
                                  h4("Show the Ribo-seq SSCBM offset bar plot"),
                                  plotOutput("out_step6_tis_bar")),
                         
                         tabPanel(title = "Heat",
                                  h4("Show the Ribo-seq SSCBM offset heat plot"),
                                  plotOutput("out_step6_tis_heat"))
                       )
                )
              )
      ),
      
      
      ### step6 Ribo-seq tis offset end ###############################
      tabItem(tabName = "step6_end",
              h2("draw the end of offset of SSCBM"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                ## set the input panel
                column(5,
                       box(title = "Import the offset end.", collapsed = FALSE,
                           collapsible = TRUE, width = NULL, status = "primary", solidHeader = TRUE,
                           ## click to import the end of offset
                           fileInput(inputId = "in_step6_end", label = "Offset", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           actionButton("act_step6_import_end_offset", label = "import the offset", icon = icon("file")),
                           tags$style(HTML("#act_step6_import_end_offset {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Import the offset end.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           fluidRow(
                             column(4, textInput("in_step6_end_group", label = "Group", value = "Sample", placeholder = "need design file")),
                             column(4, checkboxInput("in_step6_end_wrap", label = "Warp", value = FALSE)),
                             column(4, selectizeInput("in_step6_end_value", label = "Value", choices = c("Count", "Log2", "Scaled"), selected = "Count"))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step6_end_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(4, textInput("in_step6_end_edge_color", label = "Edge color", value = "white")),
                             column(4, sliderInput("in_step6_end_edge_width", label = "Edge width", value = 0.3, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step6_end_x_min", label = "Xlim min", value = -30, min = -40, max = 40, step = 1)),
                             column(4, numericInput("in_step6_end_x_max", label = "Xlim max", value = 30, min = -40, max = 40, step = 1))
                             # column(4, numericInput("in_step6_end_x_breaks", label = "X breaks", value = 3, min = 1, max = 10, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step6_end_y_min", label = "Ylim min", value = 25, min = 20, max = 40, step = 1)),
                             column(4, numericInput("in_step6_end_y_max", label = "Ylim max", value = 33, min = 20, max = 40, step = 1))
                             # column(4, numericInput("in_step6_end_y_breaks", label = "Y breaks", value = 3, min = 1, max = 10, step = 1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step6_end_fill_color", label = "Fill color", selected = "Blues", choices = color_list)),
                             column(6, sliderInput("in_step6_end_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step6_end_xlabel", label = "Xlabel", value = "Offset")),
                             column(4, textInput("in_step6_end_ylabel", label = "Ylabel", value = "RPFs length")),
                             column(4, textInput("in_step6_end_title", label = "Title", value = "End of RPFs"))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step6_end_width", label = "Figure width", value = 10, min = 1)),
                             column(6, numericInput("out_step6_end_height", label = "Figure height", value = 5, min = 1))
                           ),
                           
                           textInput("out_step6_end_fig_name", label = "Output", value = "Ribo-seq", placeholder = "file name prefix"),
                           
                           fluidRow(
                             column(4, actionButton("act_step6_draw_end_offset", label = "draw the offset", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step6_tis_end_heat", label = "Save tis heat", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step6_tts_end_heat", label = "Save tts heat", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step6_draw_end_offset {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step6_tis_end_heat {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step6_tts_end_heat {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Offset",
                                  h4("Show the end of TIS/TTS offset table"),
                                  DT::dataTableOutput("out_step6_end_offset")),
                         
                         ## set the TIS offset end plot
                         tabPanel(title = "TIS_Heat",
                                  h4("Show the Ribo-seq end TIS offset heat plot"),
                                  plotOutput("out_step6_tis_end_heat")),
                         
                         tabPanel(title = "TTS_Heat",
                                  h4("Show the Ribo-seq end TTS offset heat plot"),
                                  plotOutput("out_step6_tts_end_heat"))
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step7 3nt periodicity ###############################
      
      ### step7 RNA-seq period ###############################
      tabItem(tabName = "step7_rnaseq",
              h2("draw the period of RNA-seq"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                ## set the input panel
                column(5,
                       box(title = "Import the periodicity.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           ## click to import the 3-nt periodicity
                           fileInput(inputId = "in_step7_rnaseq", label = "Period", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           actionButton("act_step7_import_rna_period", label = "import the period", icon = icon("file")),
                           tags$style(HTML("#act_step7_import_rna_period {background-color: #75aadb; color: black; margin-top: 0px;}")),
                       ),
                       
                       box(title = "Draw the periodicity.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           textInput("in_step7_rna_x", label = "X-axis", value = "Sample", placeholder = "need design file"),
                           
                           fluidRow(
                             column(6, checkboxInput("in_step7_rna_text", label = "Label", value = FALSE)),
                             column(6, checkboxInput("in_step7_rna_trans", label = "Flip", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step7_rna_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(6, sliderInput("in_step7_rna_label_size", label = "Label size", value = 3, min = 0, max = 10, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step7_rna_edge_color", label = "Edge color", value = "white")),
                             column(6, selectizeInput("in_step7_rna_fill_color", label = "Fill color", selected = "Blues", choices = color_list))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step7_rna_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1)),
                             column(6, sliderInput("in_step7_rna_edge_width", label = "Edge width", value = 0.3, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step7_rna_width", label = "Figure width", value = 10, min = 1)),
                             column(6, numericInput("out_step7_rna_height", label = "Figure height", value = 5, min = 1))
                           ),
                           
                           textInput("out_step7_rna_fig_name", label = "Output", value = "RNA-seq", placeholder = "file name prefix"),
                           
                           fluidRow(
                             column(4, actionButton("act_step7_draw_rna_period", label = "draw the period", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step7_rna_count", label = "Save count", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step7_rna_ratio", label = "Save ratio", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step7_draw_rna_period {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step7_rna_count {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step7_rna_ratio {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Periodicity",
                                  h4("Show the RNA-seq period table"),
                                  DT::dataTableOutput("out_step7_rna_period")),
                         
                         ## set the 3-nt periodicity plot
                         tabPanel(title = "Count",
                                  h4("Show the RNA-seq period bar plot"),
                                  plotOutput("out_step7_rna_count")),
                         
                         tabPanel(title = "Proportion",
                                  h4("Show the RNA-seq period bar plot"),
                                  plotOutput("out_step7_rna_ratio"))
                       )
                )
              )
      ),
      
      ### step7 Ribo-seq period ###############################
      tabItem(tabName = "step7_riboseq",
              h2("draw the period of Ribo-seq"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the periodicity.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           
                           ## click to import the 3-nt periodicity
                           fileInput(inputId = "in_step7_riboseq", label = "Period", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           actionButton("act_step7_import_ribo_period", label = "import the period", icon = icon("file")),
                           tags$style(HTML("#act_step7_import_ribo_period {background-color: #75aadb; color: black; margin-top: 0px;}")),
                       ),
                       
                       box(title = "Draw the periodicity.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           textInput("in_step7_ribo_x", label = "X-axis", value = "Sample", placeholder = "need design file"),
                           
                           fluidRow(
                             column(3, checkboxInput("in_step7_ribo_text", label = "Label", value = FALSE)),
                             column(3, checkboxInput("in_step7_ribo_trans", label = "Flip", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step7_ribo_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(6, sliderInput("in_step7_ribo_label_size", label = "Label size", value = 3, min = 0, max = 10, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step7_ribo_edge_color", label = "Edge color", value = "white")),
                             column(6, selectizeInput("in_step7_ribo_fill_color", label = "Fill color", selected = "Blues", choices = color_list))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step7_ribo_edge_width", label = "Edge width", value = 0.3, min = 0, max = 1, step = 0.1)),
                             column(6, sliderInput("in_step7_ribo_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step7_ribo_width", label = "Figure width", value = 10, min = 1)),
                             column(6, numericInput("out_step7_ribo_height", label = "Figure height", value = 5, min = 1))
                           ),
                           
                           textInput("out_step7_ribo_fig_name", label = "Output", value = "Ribo-seq", placeholder = "file name prefix"),
                           
                           fluidRow(
                             column(4, actionButton("act_step7_draw_ribo_period", label = "draw the period", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step7_ribo_count", label = "Save count", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step7_ribo_ratio", label = "Save ratio", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step7_draw_ribo_period {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step7_ribo_count {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step7_ribo_ratio {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Periodicity",
                                  h4("Show the Ribo-seq period table"),
                                  DT::dataTableOutput("out_step7_ribo_period")),
                         
                         ## set the 3-nt periodicity plot
                         tabPanel(title = "Count",
                                  h4("Show the Ribo-seq period bar plot"),
                                  plotOutput("out_step7_ribo_count")),
                         
                         tabPanel(title = "Proportion",
                                  h4("Show the Ribo-seq period bar plot"),
                                  plotOutput("out_step7_ribo_ratio"))
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step8 metaplot ###############################
      
      ### step8 RNA-seq metaplot ###############################
      tabItem(tabName = "step8_rnaseq",
              h2("draw the metaplot of RNA-seq"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the RNA-seq meta gene.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           
                           ## click to import the metaplot
                           fileInput(inputId = "in_step8_rnaseq", label = "Metagene", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           actionButton("act_step8_import_rna_meta", label = "import the metagene", icon = icon("file")),
                           tags$style(HTML("#act_step8_import_rna_meta {background-color: #75aadb; color: black; margin-top: 0px;}")),
                       ),
                       
                       box(title = "Draw the RNA-seq metaplot.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           
                           fluidRow(
                             column(4, textInput("in_step8_rna_x", label = "Group", value = "Sample", placeholder = "Group")),
                             column(4, selectizeInput("in_step8_rna_x_type", label = "Type", selected = "Nucleotide",
                                                      choices = c("Codon", "Nucleotide"))),
                             conditionalPanel(condition = ("input.in_step8_rna_x_type == 'Codon'"),
                                              column(4, numericInput("in_step8_rna_smooth", label = "smooth (k)", value = 0, step = 1)))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step8_rna_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(6, selectizeInput("in_step8_rna_legend", label = "Legend", selected = "top",
                                                      choices = c("top", "bottom", "left", "right")))
                           ),
                           
                           fluidRow(
                             column(6, checkboxInput("in_step8_rna_wrap", label = "Wrap", value = FALSE)),
                             column(6, sliderInput("in_step8_rna_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step8_rna_tis_x_min", label = "TIS min", value = -10, min = -30, max = 90, step = 1)),
                             column(4, numericInput("in_step8_rna_tis_x_max", label = "TIS max", value = 90, min = 90, max = 300, step = 1)),
                             column(4, checkboxInput("in_step8_rna_scale", label = "Scale", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step8_rna_tts_x_min", label = "TTS min", value = -90, min = -30, max = 90, step = 1)),
                             column(4, numericInput("in_step8_rna_tts_x_max", label = "TTS max", value = 10, min = 90, max = 300, step = 1)),
                             column(4, checkboxInput("in_step8_rna_wrap", label = "Wrap", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step8_rna_fill_color", label = "Fill color", choices = color_list, selected = "Blues")),
                             column(6, sliderInput("in_step8_rna_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step8_rna_line_color", label = "Line color", choices = color_list, selected = "Blues")),
                             column(6, sliderInput("in_step8_rna_line_width", label = "Line width", value = 1, min = 0.1, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step8_rna_width", label = "Figure width", value = 10, min = 1)),
                             column(6, numericInput("out_step8_rna_height", label = "Figure height", value = 7, min = 1))
                           ),
                           
                           textInput("out_step8_rna_fig_name", label = "Output", value = "RNA-seq", placeholder = "file name prefix"),
                           
                           actionButton("act_step8_draw_rna_meta", label = "draw the metagene", icon = icon("chart-line")),
                           tags$style(HTML("#act_step8_draw_rna_meta {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           
                           fluidRow(
                             column(4, actionButton("save_step8_rna_bar", label = "Save barplot", 
                                                    class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step8_rna_line", label = "Save lineplot", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step8_rna_heat", label = "Save heatmap", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step8_rna_bar {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step8_rna_line {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step8_rna_heat {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Metagene",
                                  h4("Show the RNA-seq metaplot table"),
                                  DT::dataTableOutput("out_step8_rna_meta")),
                         
                         ## set the metaplot plot
                         tabPanel(title = "Bar",
                                  h4("Show the RNA-seq metaplot bar plot"),
                                  plotOutput("out_step8_rna_bar")),
                         
                         tabPanel(title = "Line",
                                  h4("Show the RNA-seq metaplot line plot"),
                                  plotOutput("out_step8_rna_line")),
                         
                         tabPanel(title = "Heat",
                                  h4("Show the RNA-seq metaplot bar plot"),
                                  plotOutput("out_step8_rna_heat"))
                       )
                )
              )
      ),
      
      ### step8 Ribo-seq metaplot ###############################
      tabItem(tabName = "step8_riboseq",
              h2("draw the metaplot of Ribo-seq"),
              fluidRow(
                tags$style(HTML(".download-btn {margin-top: 250px;}",
                                ".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the Ribo-seq meta gene.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           
                           ## click to import the metaplot
                           fileInput(inputId = "in_step8_riboseq", label = "Metagene", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           actionButton("act_step8_import_ribo_meta", label = "import the metagene", icon = icon("file")),
                           tags$style(HTML("#act_step8_import_ribo_meta {background-color: #75aadb; color: black; margin-top: 0px;}")),
                       ),

                       box(title = "Draw the Ribo-seq metaplot.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           
                           fluidRow(
                             column(4, textInput("in_step8_ribo_x", label = "Group", value = "Sample", placeholder = "Group")),
                             column(4, selectizeInput("in_step8_ribo_x_type", label = "Type", selected = "Nucleotide",
                                                      choices = c("Codon", "Nucleotide"))),
                             conditionalPanel(condition = ("input.in_step8_ribo_x_type == 'Codon'"),
                                              column(4, numericInput("in_step8_ribo_smooth", label = "smooth (k)", value = 0, step = 1)))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step8_ribo_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(4, selectizeInput("in_step8_ribo_legend", label = "Legend", selected = "bottom",
                                                      choices = c("top", "bottom", "left", "right")))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step8_ribo_tis_x_min", label = "TIS min", value = -10, min = -30, max = 90, step = 1)),
                             column(4, numericInput("in_step8_ribo_tis_x_max", label = "TIS max", value = 90, min = 90, max = 300, step = 1)),
                             column(4, checkboxInput("in_step8_ribo_scale", label = "Scale", value = FALSE))
                           ),
                           
                           # sliderInput("in_step8_ribo_tis", label = "TIS xlim", value = c(-10, 30), min = -30, max = 90, step = 1),
                           # sliderInput("in_step8_ribo_tts", label = "TTS xlim", value = c(-30, 10), min = -90, max = 30, step = 1),
                           
                           fluidRow(
                             column(4, numericInput("in_step8_ribo_tts_x_min", label = "TTS min", value = -90, min = -300, max = 30, step = 1)),
                             column(4, numericInput("in_step8_ribo_tts_x_max", label = "TTS max", value = 10, min = 30, max = 90, step = 1)),
                             column(4, checkboxInput("in_step8_ribo_wrap", label = "Wrap", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step8_ribo_y_min", label = "Ylim min", value = 0, step = 1)),
                             column(4, numericInput("in_step8_ribo_y_max", label = "Ylim max", value = 10, step = 1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step8_ribo_fill_color", label = "Fill color", choices = color_list, selected = "Blues")),
                             column(6, sliderInput("in_step8_ribo_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step8_ribo_line_color", label = "Line color", choices = color_list, selected = "Blues")),
                             column(6, sliderInput("in_step8_ribo_line_width", label = "Line width", value = 1, min = 0.1, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step8_ribo_width", label = "Figure width", value = 10, min = 1)),
                             column(6, numericInput("out_step8_ribo_height", label = "Figure height", value = 7, min = 1))
                           ),
                           
                           textInput("out_step8_ribo_fig_name", label = "Output", value = "Ribo-seq", placeholder = "file name prefix"),
                           
                           actionButton("act_step8_draw_ribo_meta", label = "draw the metagene", icon = icon("chart-line")),
                           tags$style(HTML("#act_step8_draw_ribo_meta {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           
                           fluidRow(
                             column(4, downloadButton("save_step8_ribo_bar", label = "Save barplot", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step8_ribo_line", label = "Save lineplot", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step8_ribo_heat", label = "Save heatmap", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step8_ribo_bar {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step8_ribo_line {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step8_ribo_heat {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Metagene",
                                  h4("Show the Ribo-seq metaplot table"),
                                  DT::dataTableOutput("out_step8_ribo_meta")),

                         ## set the metaplot plot
                         tabPanel(title = "Bar",
                                  h4("Show the RNA-seq metaplot bar plot"),
                                  plotOutput("out_step8_ribo_bar")),
                         
                         tabPanel(title = "Line",
                                  h4("Show the RNA-seq metaplot line plot"),
                                  plotOutput("out_step8_ribo_line")),
                         
                         tabPanel(title = "Heat",
                                  h4("Show the Ribo-seq metaplot bar plot"),
                                  plotOutput("out_step8_ribo_heat"))
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step9 coverage ###############################
      
      ### step9 RNA-seq coverage ###############################
      tabItem(tabName = "step9_rnaseq",
              h2("draw the coverage of RNA-seq"),
              fluidRow(
                tags$style(HTML(".download-btn {margin-top: 250px;}",
                                ".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the RNA-seq coverage.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           ## click to import the coverage
                           fileInput(inputId = "in_step9_rnaseq", label = "Coverage", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           actionButton("act_step9_import_rna_coverage", label = "import the coverage", icon = icon("file")),
                           tags$style(HTML("#act_step9_import_rna_coverage {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Draw the RNA-seq coverage.", collapsed = FALSE, collapsible = TRUE,
                           width = NULL, status = "primary", solidHeader = TRUE,
                           textInput("in_step9_rna_x", label = "X-axis", value = "Sample", placeholder = "need design file"), 
                           
                           fluidRow(
                             column(3, checkboxInput("in_step9_rna_wrap", label = "Wrap", value = FALSE)),
                             column(3, checkboxInput("in_step9_rna_scale", label = "Scale", value = FALSE)),
                             column(6, sliderInput("in_step9_rna_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5))
                           ), 
                           
                           fluidRow(
                             column(6, selectizeInput("in_step9_rna_line_color", label = "Line color", choices = color_list, selected = "Paired")),
                             column(6, sliderInput("in_step9_rna_line_width", label = "Line width", value = 1, min = 0.1, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step9_rna_fill_color", label = "Fill color", choices = color_list, selected = "Blues")),
                             column(6, sliderInput("in_step9_rna_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step9_rna_width", label = "Figure width", value = 10, min = 1)),
                             column(6, numericInput("out_step9_rna_height", label = "Figure height", value = 7, min = 1))
                           ),
                           
                           textInput("out_step9_rna_fig_name", label = "Output", value = "RNA-seq", placeholder = "file name prefix"),
                           
                           fluidRow(
                             column(4, actionButton("act_step9_draw_rna_coverage", label = "draw the coverage", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step9_rna_line", label = "Save lineplot", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step9_rna_heat", label = "Save heatmap", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step9_draw_rna_coverage {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step9_rna_line {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step9_rna_heat {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Coverage",
                                  h4("Show the RNA-seq coverage table"),
                                  DT::dataTableOutput("out_step9_rna_coverage")),
                         
                         ## set the output panel
                         tabPanel(title = "Asymmetry score",
                                  h4("Show the RNA-seq asymmetry score table"),
                                  DT::dataTableOutput("out_step9_rna_as_score")),
                         
                         ## set the coverage plot
                         tabPanel(title = "Line",
                                  h4("Show the RNA-seq coverage line plot"),
                                  plotOutput("out_step9_rna_line")),
                         
                         tabPanel(title = "Heat",
                                  h4("Show the RNA-seq coverage bar plot"),
                                  plotOutput("out_step9_rna_heat"))
                       )
                )
              )
      ),
      
      ### step9 Ribo-seq coverage ###############################
      tabItem(tabName = "step9_riboseq",
              h2("draw the coverage of Ribo-seq"),
              fluidRow(
                tags$style(HTML(".download-btn {margin-top: 250px;}",
                                ".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the RNA-seq coverage.", collapsed = FALSE,
                           collapsible = TRUE, width = NULL, status = "primary", solidHeader = TRUE,
                           
                           ## click to import the coverage
                           fileInput(inputId = "in_step9_riboseq", label = "Coverage", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           actionButton("act_step9_import_ribo_coverage", label = "import the coverage", icon = icon("file")),
                           tags$style(HTML("#act_step9_import_ribo_coverage {background-color: #75aadb; color: black;}")),
                       ),
                       
                       box(title = "Import the RNA-seq coverage.", collapsed = FALSE, collapsible = TRUE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           textInput("in_step9_ribo_x", label = "X-axis", value = "Sample", placeholder = "need design file"),
                           
                           fluidRow(
                             column(3, checkboxInput("in_step9_ribo_wrap", label = "Wrap", value = FALSE)),
                             column(3, checkboxInput("in_step9_ribo_scale", label = "Scale", value = FALSE)),
                             column(6, sliderInput("in_step9_ribo_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step9_ribo_line_color", label = "Line color", choices = color_list, selected = "Paired")),
                             column(6, sliderInput("in_step9_ribo_line_width", label = "Line width", value = 1, min = 0.1, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step9_ribo_fill_color", label = "Fill color", choices = color_list, selected = "Blues")),
                             column(6, sliderInput("in_step9_ribo_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step9_ribo_width", label = "Figure width", value = 8, min = 1)),
                             column(6, numericInput("out_step9_ribo_height", label = "Figure height", value = 5, min = 1))
                           ),
                           
                           textInput("out_step9_ribo_fig_name", label = "Output", value = "Ribo-seq", placeholder = "file name prefix"),
                           
                           fluidRow(
                             column(4, actionButton("act_step9_draw_ribo_coverage", label = "draw the coverage", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step9_ribo_line", label = "Save lineplot", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step9_ribo_heat", label = "Save heatmap", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step9_draw_ribo_coverage {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step9_ribo_line {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step9_ribo_heat {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Coverage",
                                  h4("Show the Ribo-seq coverage table"),
                                  DT::dataTableOutput("out_step9_ribo_coverage")),
                         
                         ## set the output panel
                         tabPanel(title = "Asymmetry score",
                                  h4("Show the Ribo-seq asymmetry score table"),
                                  DT::dataTableOutput("out_step9_ribo_as_score")),
                         
                         ## set the coverage plot
                         tabPanel(title = "Line",
                                  h4("Show the RNA-seq coverage line plot"),
                                  plotOutput("out_step9_ribo_line")),
                         
                         tabPanel(title = "Heat",
                                  h4("Show the Ribo-seq coverage bar plot"),
                                  plotOutput("out_step9_ribo_heat"))
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step10 quantification ###############################
      
      ### step10 RNA-seq quantification ###############################
      tabItem(tabName = "step10_rnaseq",
              h2("show the abundance of RNA-seq"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the RNA-seq data.", status = "primary", solidHeader = TRUE, 
                           collapsible = TRUE, collapsed = FALSE, width = NULL,
                           
                           ## click to import the coverage
                           fileInput(inputId = "in_step10_rnaseq", label = "Count", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           actionButton("act_step10_import_rna_count", label = "import the RNA", icon = icon("file")),
                           tags$style(HTML("#act_step10_import_rna_count {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Remove the batch effect.", status = "primary", solidHeader = TRUE, 
                           collapsible = TRUE, collapsed = FALSE, width = NULL,
                           
                           fluidRow(
                             column(6, textInput("in_step10_rna_group", label = "Group", value = "Group", placeholder = "need design file")),
                             column(6, selectizeInput("in_step10_rna_seq", label = "SeqType", 
                                                      choices = c("RNA", "RIBO"), selected = "RNA"))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step10_rna_filter", label = "Filter rowSums", value = 1, min = 0)),
                             column(6, selectizeInput("in_step10_rna_method", label = "Normalise method", 
                                                      choices = c("upper", "full"), selected = "upper")),
                           ),
                           
                           fluidRow(
                             column(4, actionButton("act_step10_norm_rna_count", label = "normalize", icon = icon("calculator"))),
                             column(4, downloadButton("save_step10_rna_norm", label = "Save norm", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step10_rna_rpm", label = "Save rpm", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step10_norm_rna_count {background-color: #75aadb; color: black;}")),
                           tags$style(HTML("#save_step10_rna_norm {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step10_rna_rpm {background-color: #grey; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Show the RNA-seq expression.", status = "primary", solidHeader = TRUE, 
                           collapsible = TRUE, collapsed = FALSE, width = NULL,
                           fluidRow(
                             column(6, numericInput("in_step10_rna_ylim_min", label = "Ylim min", value = -4, step = 0.5)),
                             column(6, numericInput("in_step10_rna_ylim_max", label = "Ylim max", value = 4, step = 0.5))
                           ),
                           p("Ylim used for boxplot."),
                           
                           fluidRow(
                             column(6, numericInput("in_step10_rna_xlim_min", label = "Xlim min", value = -1, step = 0.5)),
                             column(6, numericInput("in_step10_rna_xlim_max", label = "Xlim max", value = 14, step = 0.5))
                           ),
                           p("Xlim used for eCDF."),
                           
                           # fluidRow(
                           #   column(6, sliderInput("in_step10_rna_ylim", label = "Ylim", value = c(-4, 4), min = -10, max = 10, step = 0.5)),
                           #   column(6, sliderInput("in_step10_rna_xlim", label = "Xlim", value = c(-1, 14), min = -20, max = 20, step = 0.5))
                           # ),
                           # p("Ylim used for boxplot. Xlim used for eCDF."),
                           
                           fluidRow(
                             column(6, sliderInput("in_step10_rna_font_size", label = "Font size", value = 12, min = 0, max = 20, step = 0.5)),
                             column(3, checkboxInput("in_step10_rna_log2", label = "Log2", value = TRUE)),
                             column(3, checkboxInput("in_step10_rna_wrap", label = "Wrap", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step10_rna_fill_color", label = "Fill color", selected = "Paired", choices = color_list)),
                             column(6, sliderInput("in_step10_rna_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step10_rna_line_color", label = "Line color", selected = "Paired", choices = color_list)),
                             column(6, sliderInput("in_step10_rna_line_width", label = "Line width", value = 1, min = 0.1, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step10_rna_width", label = "Figure width", value = 9, min = 1)),
                             column(6, numericInput("out_step10_rna_height", label = "Figure height", value = 4, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step10_rna_name", label = "Output", value = "RNA-seq", placeholder = "file name prefix")),
                             column(6, selectizeInput("out_step10_rna_format", label = "Format", choices = c(".txt", ".xlsx"), selected = "txt"))
                           ),
                           
                           actionButton("act_step10_draw_rna_expr", label = "draw the expression", icon = icon("chart-line")),
                           tags$style(HTML("#act_step10_draw_rna_expr {background-color: #e78c45; color: black;}")),
                           
                           fluidRow(
                             column(4, downloadButton("save_step10_rna_box", label = "Save boxplot", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step10_rna_pca", label = "Save pca", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step10_rna_ecdf", label = "Save ecdf", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step10_rna_box {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step10_rna_pca {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step10_rna_ecdf {background-color: #grey; color: black; margin-top: 5px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Abundance",
                                  h4("Show the RNA-seq raw count table"),
                                  DT::dataTableOutput("out_step10_rna_count")),
                         
                         tabPanel(title = "Normalized",
                                  h4("Show the RNA-seq normalized table"),
                                  conditionalPanel(condition = ("input.act_step10_norm_rna_count > 0"),
                                                   DT::dataTableOutput("out_step10_rna_norm") %>% 
                                                     withSpinner(color = "#3c8dbc", type = 5, size = 1)),
                                  verbatimTextOutput("out_step10_rna_norm_info"),
                                  tags$style(HTML("#out_step10_rna_norm_info {color: red; font-size: 26px;}"))),
                         
                         tabPanel(title = "RPM",
                                  h4("Show the RNA-seq RPM table"),
                                  DT::dataTableOutput("out_step10_rna_rpm")),
                         
                         ## set the abundance plot
                         tabPanel(title = "Boxplot",
                                  h4("Show the RNA-seq box plot"),
                                  plotOutput("out_step10_rna_box")),
                         
                         tabPanel(title = "PCA",
                                  h4("Show the RNA-seq PCA plot"),
                                  plotOutput("out_step10_rna_pca")),
                         
                         tabPanel(title = "eCDF",
                                  h4("Show the RNA-seq cumulative plot"),
                                  plotOutput("out_step10_rna_ecdf"))
                         
                       )
                )
              )
      ),
      
      ### step10 Ribo-seq quantification ###############################
      tabItem(tabName = "step10_riboseq",
              h2("show the abundance of Ribo-seq"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the Ribo-seq data.", status = "primary", 
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, width = NULL,
                           
                           ## click to import the coverage
                           fileInput(inputId = "in_step10_riboseq", label = "Count", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           actionButton("act_step10_import_ribo_count", label = "import the Ribo", icon = icon("file")),
                           tags$style(HTML("#act_step10_import_ribo_count {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Remove the batch effect.", status = "primary", solidHeader = TRUE, 
                           collapsible = TRUE, collapsed = FALSE, width = NULL,
                           
                           fluidRow(
                             column(6, textInput("in_step10_ribo_group", label = "Group", value = "Group", placeholder = "need design file")),
                             column(6, selectizeInput("in_step10_ribo_seq", label = "SeqType", 
                                                      choices = c("RNA", "RIBO"), selected = "RIBO"))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step10_ribo_filter", label = "Filter rowSums", value = 1, min = 0)),
                             column(6, selectizeInput("in_step10_ribo_method", label = "Normalise method", 
                                                      choices = c("upper", "full"), selected = "upper")),
                           ),
                           
                           fluidRow(
                             column(4, actionButton("act_step10_norm_ribo_count", label = "normalize", icon = icon("calculator"))),
                             column(4, downloadButton("save_step10_ribo_norm", label = "Save norm", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step10_ribo_rpm", label = "Save rpm", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step10_norm_ribo_count {background-color: #75aadb; color: black;}")),
                           tags$style(HTML("#save_step10_ribo_norm {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step10_ribo_rpm {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                       box(title = "Show the Ribo-seq expression.", status = "primary", solidHeader = TRUE, 
                           collapsible = TRUE, collapsed = FALSE, width = NULL,
                           
                           fluidRow(
                             column(6, numericInput("in_step10_ribo_ylim_min", label = "Ylim min", value = -4, step = 0.5)),
                             column(6, numericInput("in_step10_ribo_ylim_max", label = "Ylim max", value = 4, step = 0.5))
                           ),
                           p("Ylim used for boxplot."),
                           
                           fluidRow(
                             column(6, numericInput("in_step10_ribo_xlim_min", label = "Xlim min", value = -1, step = 0.5)),
                             column(6, numericInput("in_step10_ribo_xlim_max", label = "Xlim max", value = 14, step = 0.5))
                           ),
                           p("Xlim used for eCDF."),
                           
                           fluidRow(
                             column(6, sliderInput("in_step10_ribo_font_size", label = "Font size", value = 12, min = 0, max = 20, step = 0.5)),
                             column(3, checkboxInput("in_step10_ribo_log2", label = "Log2", value = TRUE)),
                             column(3, checkboxInput("in_step10_ribo_wrap", label = "Wrap", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step10_ribo_fill_color", label = "Fill color", selected = "Paired", choices = color_list)),
                             column(6, sliderInput("in_step10_ribo_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step10_ribo_line_color", label = "Line color", selected = "Paired", choices = color_list)),
                             column(6, sliderInput("in_step10_ribo_line_width", label = "Line width", value = 1, min = 0.1, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step10_ribo_width", label = "Figure width", value = 9, min = 1)),
                             column(6, numericInput("out_step10_ribo_height", label = "Figure height", value = 4, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step10_ribo_name", label = "Output", value = "Ribo-seq", placeholder = "file name prefix")),
                             column(6, selectizeInput("out_step10_ribo_format", label = "Format", choices = c(".txt", ".xlsx"), selected = "txt"))
                           ),
                           
                           actionButton("act_step10_draw_ribo_expr", label = "draw the expression", icon = icon("chart-line")),
                           tags$style(HTML("#act_step10_draw_ribo_expr {background-color: #e78c45; color: black;}")),
                           
                           fluidRow(
                             column(4, downloadButton("save_step10_ribo_box", label = "Save boxplot", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step10_ribo_pca", label = "Save pca", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step10_ribo_ecdf", label = "Save ecdf", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step10_ribo_box {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step10_ribo_pca {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step10_ribo_ecdf {background-color: #grey; color: black; margin-top: 5px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Abundance",
                                  h4("Show the RNA-seq raw count table"),
                                  DT::dataTableOutput("out_step10_ribo_count")),
                         
                         tabPanel(title = "Normalized",
                                  h4("Show the RNA-seq normalized table"),
                                  conditionalPanel(condition = ("input.act_step10_norm_ribo_count > 0"),
                                                   DT::dataTableOutput("out_step10_ribo_norm") %>% 
                                                     withSpinner(color = "#3c8dbc", type = 5, size = 1)),
                                  verbatimTextOutput("out_step10_ribo_norm_info"),
                                  tags$style(HTML("#out_step10_ribo_norm_info {color: red; font-size: 26px;}"))),
                         
                         tabPanel(title = "RPM",
                                  h4("Show the RNA-seq RPM table"),
                                  DT::dataTableOutput("out_step10_ribo_rpm")),
                         
                         ## set the abundance plot
                         tabPanel(title = "Boxplot",
                                  h4("Show the Ribo-seq box plot"),
                                  plotOutput("out_step10_ribo_box")),
                         
                         tabPanel(title = "PCA",
                                  h4("Show the Ribo-seq PCA plot"),
                                  plotOutput("out_step10_ribo_pca")),
                         
                         tabPanel(title = "eCDF",
                                  h4("Show the Ribo-seq cumulative plot"),
                                  plotOutput("out_step10_ribo_ecdf"))
                         
                       )
                )
              )
      ),
      
      ### step10 TE quantification ###############################
      tabItem(tabName = "step10_te",
              h2("draw the abundance of TE"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the RNA-seq and RIbo-seq data", collapsible = TRUE, collapsed = FALSE,
                           status = "primary", solidHeader = TRUE, width = NULL,
                           
                           ## click to import the normalized RNA and Ribo expression table
                           fileInput(inputId = "in_step10_te_rna", label = "Norm RNA", accept = c(".txt", ".TXT"), multiple = FALSE),
                           fileInput(inputId = "in_step10_te_ribo", label = "Norm Ribo", accept = c(".txt", ".TXT"), multiple = FALSE),
                           
                           fluidRow(
                             column(6, numericInput("in_step10_te_filter", label = "Filter rowSums", value = 1, min = 0)),
                             column(6, checkboxInput("in_step10_te_linked", label = "Linked", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step10_te_fill", label = "Fill 0", 
                                                      choices = c("none", "min", "mean", "median"), selected = "none")),
                             column(6, checkboxInput("in_step10_te_clean", label = "NA / INF", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step10_import_te_calc", label = "Calculate the TE", icon = icon("calculator"))),
                             column(6, downloadButton("save_step10_te_norm", label = "Save data", class = "download-btn", icon = icon("download")))
                           ), 
                           
                           tags$style(HTML("#act_step10_import_te_calc {background-color: #75aadb; color: black;}")),
                           tags$style(HTML("#save_step10_te_norm {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                       box(title = "Show the TE abundance.", status = "primary", collapsed = FALSE, collapsible = TRUE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, textInput("in_step10_te_group", label = "Group", value = "Group", placeholder = "need design file")),
                             column(4, checkboxInput("in_step10_te_wrap", label = "Wrap", value = FALSE)),
                             column(4, checkboxInput("in_step10_te_log2", label = "Log2", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step10_te_xlim", label = "Xlim", value = c(-10, 10), min = -20, max = 20, step = 0.5)),
                             column(6, sliderInput("in_step10_te_font_size", label = "Font size", value = 12, min = 0, max = 20, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step10_te_fill_color", label = "Fill color", selected = "Paired", choices = color_list)),
                             column(6, sliderInput("in_step10_te_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step10_te_line_color", label = "Line color", selected = "Paired", choices = color_list)),
                             column(6, sliderInput("in_step10_te_line_width", label = "Line width", value = 1, min = 0.1, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step10_te_name", label = "Output", value = "TE", placeholder = "file name prefix")),
                             column(3, numericInput("out_step10_te_width", label = "Figure width", value = 5, min = 1)),
                             column(3, numericInput("out_step10_te_height", label = "Figure height", value = 4, min = 1))
                           ),
                           
                           fluidRow(
                             column(4, actionButton("act_step10_draw_te_expr", label = "draw the TE", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step10_te_pca", label = "Save pca", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step10_te_ecdf", label = "Save ecdf", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step10_draw_te_expr {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step10_te_pca {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step10_te_ecdf {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "TE",
                                  h4("Show the TE calculated with the RNA-seq and Ribo-seq"),
                                  DT::dataTableOutput("out_step10_te_norm")),
                         
                         ## set the abundance plot
                         tabPanel(title = "PCA",
                                  h4("Show the TE PCA plot"),
                                  conditionalPanel(condition = ("input.act_step10_draw_te_expr > 0"),
                                                   plotOutput("out_step10_te_pca") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "eCDF",
                                  h4("Show the TE cumulative plot"),
                                  conditionalPanel(condition = ("input.act_step10_draw_te_expr > 0"),
                                                   plotOutput("out_step10_te_ecdf") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                         
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step11 heatmap ###############################
      
      ### step11 RNA-seq heatmap ###############################
      tabItem(tabName = "step11_rnaseq",
              h2("show the abundance of RNA-seq"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the RNA-seq data.", status = "primary", solidHeader = TRUE, collapsed = FALSE,
                           collapsible = TRUE, width = NULL,
                           
                           ## click to import the expression table
                           fileInput(inputId = "in_step11_rnaseq", label = "Norm", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           fluidRow(
                             column(6, textInput("in_step11_rna_group", label = "Group", value = "Sample", placeholder = "need design file")),
                             column(6, numericInput("in_step11_rna_top", label = "Top number", value = 0))
                           ),
                           p("Only the top highly expressed genes will be shown in the heatmap. default: 0, use all genes.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           fluidRow(
                             column(6, numericInput("in_step11_rna_rowsum", label = "Filter rowSums", value = 1, min = 0)),
                             column(6, numericInput("in_step11_rna_colsum", label = "Filter sample ratio", value = 30, min = 0, max = 100))
                           ),
                           p("The rowSums and colSums are used to filter the low expression genes and samples.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           p("rowSums is the sum of each gene, and colSums is the ratio of effective samples.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           fluidRow(
                             column(6, actionButton("act_step11_import_rna_norm", label = "import the RNA", icon = icon("file"))),
                             column(6, downloadButton("save_step11_rna_corr", label = "save corr", class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step11_import_rna_norm {background-color: #75aadb; color: black;}")),
                           tags$style(HTML("#save_step11_rna_corr {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                       box(title = "Draw the correlation.", status = "primary", solidHeader = TRUE, collapsed = FALSE,
                           collapsible = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, checkboxInput("in_step11_rna_corr_number", label = "Number", value = TRUE)),
                             column(4, numericInput("in_step11_rna_corr_num_format", label = "Decimal", value = 2, min = 0, max = 10)),
                             column(4, textInput("in_step11_rna_corr_num_color", label = "Number color", value = "grey30"))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step11_rna_corr_number_size", label = "Number size", value = 10, min = 5, max = 20)),
                             column(4, numericInput("in_step11_rna_corr_font_size", label = "Font size", value = 15, min = 5, max = 20)),
                             column(4, numericInput("in_step11_rna_corr_angle", label = "Angle", value = 90, min = 0, max = 360))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step11_rna_corr_min", label = "Colormap min", value = 0, min = -1, max = 1, step = 0.1)),
                             column(4, numericInput("in_step11_rna_corr_max", label = "Colormap max", value = 1, min = -1, max = 1, step = 0.1)),
                             column(4, textInput("in_step11_rna_corr_border_color", label = "Border color", value = "white", placeholder = "eg: grey30"))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step11_rna_corr_fill_color", label = "Fill color", selected = 'RdBu', choices = color_list)),
                             column(6, sliderInput("in_step11_rna_corr_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step11_rna_corr_row_gap", label = "Row gaps", value = 0, placeholder = "eg: 2,4")),
                             column(6, textInput("in_step11_rna_corr_col_gap", label = "Col gaps", value = 0, placeholder = "eg: 2,4"))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step11_rna_corr_name", label = "Output", value = "RNA-seq", placeholder = "file name prefix")),
                             column(3, numericInput("out_step11_rna_corr_width", label = "Figure width", value = 8, min = 1)),
                             column(3, numericInput("out_step11_rna_corr_height", label = "Figure height", value = 8, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step11_draw_rna_corr", label = "draw the correlation", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step11_rna_corr_plot", label = "Save corr plot",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step11_draw_rna_corr {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step11_rna_corr_plot {background-color: #grey; color: black; margin-top: 5px;}"))
                           
                       ),
                       
                       box(title = "Draw the heatmap.", status = "primary", solidHeader = TRUE, collapsed = FALSE,
                           collapsible = TRUE, width = NULL,
                           
                           textAreaInput("in_step11_rna_gene_name", label = "Gene name", value = '', placeholder = "eg: gene1, gene2, gene3"),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step11_rna_method", label = "Cluster method", selected = "complete",
                                                      choices = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))),
                             column(4, selectizeInput("in_step11_rna_distance", label = "Cluster distance", selected = "correlation",
                                                      choices = c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))),
                             column(4, numericInput("in_step11_rna_font_size", label = "Font size", value = 15, min = 5, max = 20))
                           ),
                           
                           fluidRow(
                             column(3, checkboxInput("in_step11_rna_row_clust", label = "Row cluster", value = TRUE)),
                             column(3, checkboxInput("in_step11_rna_col_clust", label = "Col cluster", value = FALSE)),
                             column(3, checkboxInput("in_step11_rna_rowname", label = "Row name", value = FALSE)),
                             column(3, checkboxInput("in_step11_rna_colname", label = "Col name", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step11_rna_row_gap", label = "Row gaps", value = 0, placeholder = "eg: 2,4")),
                             column(4, textInput("in_step11_rna_col_gap", label = "Col gaps", value = 0, placeholder = "eg: 2,4")),
                             column(3, checkboxInput("in_step11_rna_log2", label = "Log2", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step11_rna_fill_color", label = "Fill color", selected = 'RdBu', choices = color_list)),
                             column(4, sliderInput("in_step11_rna_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1)),
                             column(4, checkboxInput("in_step11_rna_border_color", label = "Border color", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step11_rna_name", label = "Output", value = "RNA-seq", placeholder = "file name prefix")),
                             column(3, numericInput("out_step11_rna_width", label = "Figure width", value = 8, min = 1)),
                             column(3, numericInput("out_step11_rna_height", label = "Figure height", value = 8, min = 1))
                           ),
                           
                           fluidRow(
                             column(4, actionButton("act_step11_draw_rna_expr", label = "draw the heatmap", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step11_rna_scale", label = "Save scaled",
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step11_rna_unscale", label = "Save unscaled",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step11_draw_rna_expr {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step11_rna_scale {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step11_rna_unscale {background-color: #grey; color: black; margin-top: 5px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Abundance",
                                  h4("Show the RNA-seq norm table"),
                                  conditionalPanel(condition = ("input.act_step11_import_rna_norm > 0"),
                                                   DT::dataTableOutput("out_step11_rna_norm") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Correlation",
                                  h4("Show the RNA-seq correlation table"),
                                  DT::dataTableOutput("out_step11_rna_corr_tab")),
                         
                         tabPanel(title = "Heat-corr",
                                  h4("Show the RNA-seq correlation plot"),
                                  conditionalPanel(condition = ("input.act_step11_draw_rna_corr > 0"),
                                                   plotOutput("out_step11_rna_corr_plot") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         ## set the abundance plot
                         tabPanel(title = "Heat-scaled",
                                  h4("Show the RNA-seq scaled heat map"),
                                  conditionalPanel(condition = ("input.act_step11_draw_rna_expr > 0"),
                                                   plotOutput("out_step11_rna_scale") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Heat-unscaled",
                                  h4("Show the RNA-seq un-scaled heat map"),
                                  conditionalPanel(condition = ("input.act_step11_draw_rna_expr > 0"),
                                                   plotOutput("out_step11_rna_unscale") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                         
                       )
                )
              )
      ),
      
      ### step11 Ribo-seq heatmap ###############################
      tabItem(tabName = "step11_riboseq",
              h2("show the abundance of ribo-seq"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the Ribo-seq data", status = "primary", solidHeader = TRUE, collapsed = FALSE,
                           collapsible = TRUE, width = NULL,
                           
                           ## click to import the expression table
                           fileInput(inputId = "in_step11_riboseq", label = "Norm", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           fluidRow(
                             column(6, textInput("in_step11_ribo_group", label = "Group", value = "Sample", placeholder = "need design file")),
                             column(6, numericInput("in_step11_ribo_top", label = "Top number", value = 0))
                           ),
                           p("Only the top highly expressed genes will be shown in the heatmap. default: 0, use all genes.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           fluidRow(
                             column(6, numericInput("in_step11_ribo_rowsum", label = "Filter rowSums", value = 1, min = 0)),
                             column(6, numericInput("in_step11_ribo_colsum", label = "Filter sample ratio", value = 30, min = 0, max = 100))
                           ),
                           p("The rowSums and colSums are used to filter the low expression genes and samples.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           p("rowSums is the sum of each gene, and colSums is the ratio of effective samples.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           fluidRow(
                             column(6, actionButton("act_step11_import_ribo_norm", label = "import the ribo", icon = icon("file"))),
                             column(6, downloadButton("save_step11_ribo_corr", label = "save corr", class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step11_import_ribo_norm {background-color: #75aadb; color: black;}")),
                           tags$style(HTML("#save_step11_ribo_corr {background-color: #grey; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Draw the correlation.", status = "primary", solidHeader = TRUE, collapsed = FALSE,
                           collapsible = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, checkboxInput("in_step11_ribo_corr_number", label = "Number", value = TRUE)),
                             column(4, numericInput("in_step11_ribo_corr_num_format", label = "Decimal", value = 2, min = 0, max = 10)),
                             column(4, textInput("in_step11_ribo_corr_num_color", label = "Number color", value = "grey30"))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step11_ribo_corr_number_size", label = "Number size", value = 10, min = 5, max = 20)),
                             column(4, numericInput("in_step11_ribo_corr_font_size", label = "Font size", value = 15, min = 5, max = 20)),
                             column(4, numericInput("in_step11_ribo_corr_angle", label = "Angle", value = 90, min = 0, max = 360))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step11_ribo_corr_min", label = "Colormap min", value = 0, min = -1, max = 1, step = 0.1)),
                             column(4, numericInput("in_step11_ribo_corr_max", label = "Colormap max", value = 1, min = -1, max = 1, step = 0.1)),
                             column(4, textInput("in_step11_ribo_corr_border_color", label = "Border color", value = "white", placeholder = "eg: grey30"))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step11_ribo_corr_fill_color", label = "Fill color", choices = color_list, selected = 'RdBu')),
                             column(6, sliderInput("in_step11_ribo_corr_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step11_ribo_corr_row_gap", label = "Row gaps", value = 0, placeholder = "eg: 2,4")),
                             column(6, textInput("in_step11_ribo_corr_col_gap", label = "Col gaps", value = 0, placeholder = "eg: 2,4"))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step11_ribo_corr_name", label = "Output", value = "Ribo-seq", placeholder = "file name prefix")),
                             column(3, numericInput("out_step11_ribo_corr_width", label = "Figure width", value = 8, min = 1)),
                             column(3, numericInput("out_step11_ribo_corr_height", label = "Figure height", value = 8, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step11_draw_ribo_corr", label = "draw the correlation", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step11_ribo_corr_plot", label = "Save corr plot",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step11_draw_ribo_corr {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step11_ribo_corr_plot {background-color: #grey; color: black; margin-top: 5px;}"))
                           
                       ),
                       
                       box(title = "Draw the heatmap", status = "primary", collapsed = FALSE, collapsible = TRUE,
                           solidHeader = TRUE, width = NULL,
                           
                           textAreaInput("in_step11_ribo_gene_name", label = "Gene name", value = '', placeholder = "eg: gene1, gene2, gene3"),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step11_ribo_method", label = "Cluster method", selected = "complete",
                                                      choices = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))),
                             column(4, selectizeInput("in_step11_ribo_distance", label = "Cluster distance", selected = "correlation",
                                                      choices = c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))),
                             column(4, numericInput("in_step11_ribo_font_size", label = "Font size", value = 15, min = 5, max = 20))
                           ),
                           
                           fluidRow(
                             column(3, checkboxInput("in_step11_ribo_row_clust", label = "Row cluster", value = TRUE)),
                             column(3, checkboxInput("in_step11_ribo_col_clust", label = "Col cluster", value = FALSE)),
                             column(3, checkboxInput("in_step11_ribo_rowname", label = "Row name", value = FALSE)),
                             column(3, checkboxInput("in_step11_ribo_colname", label = "Col name", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step11_ribo_row_gap", label = "Row gaps", value = 0, placeholder = "eg: 2,4")),
                             column(4, textInput("in_step11_ribo_col_gap", label = "Col gaps", value = 0, placeholder = "eg: 2,4")),
                             column(3, checkboxInput("in_step11_ribo_log2", label = "Log2", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step11_ribo_fill_color", label = "Fill color", choices = color_list, selected = 'RdBu')),
                             column(4, sliderInput("in_step11_ribo_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1)),
                             column(4, checkboxInput("in_step11_ribo_border_color", label = "Border color", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step11_ribo_name", label = "Output", value = "Ribo-seq", placeholder = "file name prefix")),
                             column(3, numericInput("out_step11_ribo_width", label = "Figure width", value = 8, min = 1)),
                             column(3, numericInput("out_step11_ribo_height", label = "Figure height", value = 8, min = 1))
                           ),
                           
                           fluidRow(
                             column(4, actionButton("act_step11_draw_ribo_expr", label = "draw the expression", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step11_ribo_scale", label = "Save scaled",
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step11_ribo_unscale", label = "Save unscaled",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step11_draw_ribo_expr {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step11_ribo_scale {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step11_ribo_unscale {background-color: #grey; color: black; margin-top: 5px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Abundance",
                                  h4("Show the ribo-seq norm table"),
                                  conditionalPanel(condition = ("input.act_step11_import_ribo_norm > 0"),
                                                   DT::dataTableOutput("out_step11_ribo_norm") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Correlation",
                                  h4("Show the ribo-seq correlation table"),
                                  DT::dataTableOutput("out_step11_ribo_corr_tab")),
                         
                         tabPanel(title = "Heat-corr",
                                  h4("Show the ribo-seq correlation plot"),
                                  conditionalPanel(condition = ("input.act_step11_draw_ribo_corr > 0"),
                                                   plotOutput("out_step11_ribo_corr_plot") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         ## set the abundance plot
                         tabPanel(title = "Heat-scaled",
                                  h4("Show the ribo-seq scaled heat map"),
                                  conditionalPanel(condition = ("input.act_step11_draw_ribo_expr > 0"),
                                                   plotOutput("out_step11_ribo_scale") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Heat-unscaled",
                                  h4("Show the ribo-seq un-scaled heat map"),
                                  conditionalPanel(condition = ("input.act_step11_draw_ribo_expr > 0"),
                                                   plotOutput("out_step11_ribo_unscale") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                         
                       )
                )
              )
      ),
      
      ### step11 TE heatmap ###############################
      tabItem(tabName = "step11_te",
              h2("show the abundance of TE"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the TE data", status = "primary", collapsed = FALSE, collapsible = TRUE,
                           solidHeader = TRUE, width = NULL,
                           
                           ## click to import the expression table
                           fileInput(inputId = "in_step11_te", label = "Norm", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           fluidRow(
                             column(6, textInput("in_step11_te_group", label = "Group", value = "Sample", placeholder = "need design file")),
                             column(6, numericInput("in_step11_te_top", label = "Top number", value = 0))
                           ),
                           p("Only the top highly expressed genes will be shown in the heatmap. default: 0, use all genes.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           fluidRow(
                             column(6, numericInput("in_step11_te_rowsum", label = "Filter rowSums", value = 0, min = 0)),
                             column(6, numericInput("in_step11_te_colsum", label = "Filter sample ratio", value = 30, min = 0, max = 100))
                           ),
                           p("The rowSums and colSums are used to filter the low expression genes and samples.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           p("rowSums is the sum of each gene, and colSums is the ratio of effective samples.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           fluidRow(
                             column(6, actionButton("act_step11_import_te_norm", label = "import the TE", icon = icon("file"))),
                             column(6, downloadButton("save_step11_te_corr", label = "Save the TE", class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step11_import_te_norm {background-color: #75aadb; color: black;}")),
                           tags$style(HTML("#save_step11_te_corr {background-color: #grey; color: black; margin-top: 5px;}"))
                       ),
                       
                       box(title = "Draw the correlation.", status = "primary", solidHeader = TRUE, collapsed = FALSE,
                           collapsible = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, checkboxInput("in_step11_te_corr_number", label = "Number", value = TRUE)),
                             column(4, numericInput("in_step11_te_corr_num_format", label = "Decimal", value = 2, min = 0, max = 10)),
                             column(4, textInput("in_step11_te_corr_num_color", label = "Number color", value = "grey30"))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step11_te_corr_number_size", label = "Number size", value = 10, min = 5, max = 20)),
                             column(4, numericInput("in_step11_te_corr_font_size", label = "Font size", value = 15, min = 5, max = 20)),
                             column(4, numericInput("in_step11_te_corr_angle", label = "Angle", value = 90, min = 0, max = 360))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step11_te_corr_min", label = "Colormap min", value = 0, min = -1, max = 1, step = 0.1)),
                             column(4, numericInput("in_step11_te_corr_max", label = "Colormap max", value = 1, min = -1, max = 1, step = 0.1)),
                             column(4, textInput("in_step11_te_corr_border_color", label = "Border color", value = "white", placeholder = "eg: grey30"))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step11_te_corr_fill_color", label = "Fill color", choices = color_list, selected = 'RdBu')),
                             column(6, sliderInput("in_step11_te_corr_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step11_te_corr_row_gap", label = "Row gaps", value = 0, placeholder = "eg: 2,4")),
                             column(6, textInput("in_step11_te_corr_col_gap", label = "Col gaps", value = 0, placeholder = "eg: 2,4"))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step11_te_corr_name", label = "Output", value = "TE", placeholder = "file name prefix")),
                             column(3, numericInput("out_step11_te_corr_width", label = "Figure width", value = 8, min = 1)),
                             column(3, numericInput("out_step11_te_corr_height", label = "Figure height", value = 8, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step11_draw_te_corr", label = "draw the correlation", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step11_te_corr_plot", label = "Save corr plot",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step11_draw_te_corr {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step11_te_corr_plot {background-color: #grey; color: black; margin-top: 5px;}"))
                           
                       ),
                       
                       box(title = "Draw the heatmap", status = "primary", collapsed = FALSE, collapsible = TRUE,
                           solidHeader = TRUE, width = NULL,
                           
                           textAreaInput("in_step11_te_gene_name", label = "Gene name", value = '', placeholder = "eg: gene1,gene2,gene3"),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step11_te_method", label = "Cluster method", selected = "complete",
                                                      choices = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))),
                             column(4, selectizeInput("in_step11_te_distance", label = "Cluster distance", selected = "correlation",
                                                      choices = c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))),
                             column(4, numericInput("in_step11_te_font_size", label = "Font size", value = 15, min = 5, max = 20))
                           ),
                           
                           fluidRow(
                             column(3, checkboxInput("in_step11_te_row_clust", label = "Row cluster", value = TRUE)),
                             column(3, checkboxInput("in_step11_te_col_clust", label = "Col cluster", value = FALSE)),
                             column(3, checkboxInput("in_step11_te_rowname", label = "Row name", value = FALSE)),
                             column(3, checkboxInput("in_step11_te_colname", label = "Col name", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step11_te_row_gap", label = "Row gaps", value = 0, placeholder = "eg: 2,4")),
                             column(4, textInput("in_step11_te_col_gap", label = "Col gaps", value = 0, placeholder = "eg: 2,4")),
                             column(3, checkboxInput("in_step11_te_log2", label = "Log2", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step11_te_fill_color", label = "Fill color", choices = color_list, selected = 'RdBu')),
                             column(4, sliderInput("in_step11_te_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.1)),
                             column(4, checkboxInput("in_step11_te_border_color", label = "Border color", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step11_te_width", label = "Figure width", value = 8, min = 1)),
                             column(6, numericInput("out_step11_te_height", label = "Figure height", value = 8, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step11_te_name", label = "Output", value = "TE", placeholder = "file name prefix")),
                             column(6, selectizeInput("out_step11_te_format", label = "Format", choices = c(".pdf", ".png"), selected = "pdf"))
                           ),
                           
                           fluidRow(
                             column(4, actionButton("act_step11_draw_te_expr", label = "draw the expression", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step11_te_scale", label = "Save scaled",
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step11_te_unscale", label = "Save unscaled",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step11_draw_te_expr {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step11_te_scale {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step11_te_unscale {background-color: #grey; color: black; margin-top: 5px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Abundance",
                                  h4("Show the TE table"),
                                  conditionalPanel(condition = ("input.act_step11_import_te_norm > 0"),
                                                   DT::dataTableOutput("out_step11_te_norm") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Correlation",
                                  h4("Show the TE correlation table"),
                                  DT::dataTableOutput("out_step11_te_corr_tab")),
                         
                         tabPanel(title = "Heat-corr",
                                  h4("Show the TE correlation plot"),
                                  conditionalPanel(condition = ("input.act_step11_draw_te_corr > 0"),
                                                   plotOutput("out_step11_te_corr_plot") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         ## set the abundance plot
                         tabPanel(title = "Heat-scaled",
                                  h4("Show the TE scaled heat map"),
                                  conditionalPanel(condition = ("input.act_step11_draw_te_expr > 0"),
                                                   plotOutput("out_step11_te_scale") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Heat-unscaled",
                                  h4("Show the TE un-scaled heat map"),
                                  conditionalPanel(condition = ("input.act_step11_draw_te_expr > 0"),
                                                   plotOutput("out_step11_te_unscale") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                         
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step12 PCA ###############################
      
      ### step12 RNA-seq PCA ###############################
      tabItem(tabName = "step12_rnaseq",
              h2("show the PCA of RNA-seq"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the RNA-seq data", status = "primary", collapsed = FALSE, collapsible = TRUE, 
                           solidHeader = TRUE, width = NULL,
                           
                           ## click to import the expression table
                           fileInput(inputId = "in_step12_rnaseq", label = "Norm", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           fluidRow(
                             column(7, textInput("in_step12_rna_group", label = "Group", value = "Group", placeholder = "need design file")),
                             column(5, numericInput("in_step12_rna_top", label = "Top number", value = 0))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step12_rna_rowsum", label = "Filter rowSums", value = 1, min = 0)),
                             column(4, numericInput("in_step12_rna_colsum", label = "Filter sample ratio", value = 30, min = 0, max = 100)),
                             column(4, checkboxInput("in_step12_rna_log2", label = "Log2", value = TRUE))
                           ),
                           p("The rowSums and colSums are used to filter the low expression genes and samples.",
                             "rowSums is the sum of each gene, and colSums is the ratio of effective samples.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           textInput("in_step12_rna_gene_name", label = "Gene name", value = '', placeholder = "eg: gene1,gene2,gene3"),
                           
                           actionButton("act_step12_import_rna_norm", label = "import the RNA", icon = icon("file")),
                           tags$style(HTML("#act_step12_import_rna_norm {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Draw the PCA", status = "primary", collapsed = FALSE, collapsible = TRUE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(3, numericInput("in_step12_rna_xmin", label = "Xlim down", value = NULL, width = 100)),
                             column(3, numericInput("in_step12_rna_xmax", label = "Xlim up", value = NULL)),
                             column(3, numericInput("in_step12_rna_ymin", label = "Ylim down", value = NULL)),
                             column(3, numericInput("in_step12_rna_ymax", label = "Ylim up", value = NULL))
                           ),
                           
                           fluidRow(
                             column(4, checkboxInput("in_step12_rna_label_sample", label = "Label sample", value = TRUE)),
                             column(4, sliderInput("in_step12_rna_label_size", label = "Label size", value = 3, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step12_rna_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step12_rna_dot_color", label = "Dot color", choices = color_list, selected = "Paired")),
                             column(4, sliderInput("in_step12_rna_dot_size", label = "Dot size", value = 3, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step12_rna_dot_alpha", label = "Alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step12_rna_width", label = "Figure width", value = 7, min = 1)),
                             column(6, numericInput("out_step12_rna_height", label = "Figure height", value = 7, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step12_rna_name", label = "Output", value = "RNA-seq", placeholder = "file name prefix")),
                             column(6, selectizeInput("out_step12_rna_format", label = "Format", choices = c(".pdf", ".png"), selected = "pdf"))
                           ),
                           
                           fluidRow(
                             column(4, actionButton("act_step12_draw_rna_pca", label = "run PCA", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step12_rna_eigen", label = "Save results", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step12_rna_eigen_plot", label = "Save eigen plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step12_draw_rna_pca {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step12_rna_eigen {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step12_rna_eigen_plot {background-color: #grey; color: black; margin-top: 0px;}")),
                           
                           fluidRow(
                             column(4, downloadButton("save_step12_rna_scatter", label = "Save scatter", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step12_rna_variances", label = "Save variances", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step12_rna_individual", label = "Save individual", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step12_rna_scatter {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step12_rna_variances {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step12_rna_individual {background-color: #grey; color: black; margin-top: 5px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Abundance",
                                  h4("Show the RNA-seq norm table"),
                                  DT::dataTableOutput("out_step12_rna_norm")),
                         
                         tabPanel(title = "Eigenvalue",
                                  h4("Show the RNA-seq eigenvalue table"),
                                  DT::dataTableOutput("out_step12_rna_eigen")),
                         
                         tabPanel(title = "Var.contrib",
                                  h4("Show the RNA-seq contribution table"),
                                  DT::dataTableOutput("out_step12_rna_var_contrib")),
                         
                         tabPanel(title = "Ind.contrib",
                                  DT::dataTableOutput("out_step12_rna_ind_contrib"),
                                  h4("Show the RNA-seq contribution table")),
                         
                         ## set the PCA plot
                         tabPanel(title = "EigenPlot",
                                  h4("Show the RNA-seq eigen plot"),
                                  plotOutput("out_step12_rna_eigen_plot")),
                         
                         tabPanel(title = "Scatter",
                                  h4("Show the RNA-seq scatter plot"),
                                  plotOutput("out_step12_rna_scatter")),
                         
                         tabPanel(title = "Variances",
                                  h4("Show the RNA-seq variances plot"),
                                  plotOutput("out_step12_rna_variances")),
                         
                         tabPanel(title = "Individual",
                                  h4("Show the RNA-seq individual plot"),
                                  plotOutput("out_step12_rna_individual"))
                         
                       )
                )
              )
      ),
      
      ### step11 Ribo-seq PCA ###############################
      tabItem(tabName = "step12_riboseq",
              h2("show the PCA of Ribo-seq"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the Ribo-seq data", status = "primary", collapsed = FALSE, collapsible = TRUE, 
                           solidHeader = TRUE, width = NULL,
                           
                           ## click to import the expression table
                           fileInput(inputId = "in_step12_riboseq", label = "Norm", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           fluidRow(
                             column(7, textInput("in_step12_ribo_group", label = "Group", value = "Group", placeholder = "need design file")),
                             column(5, numericInput("in_step12_ribo_top", label = "Top number", value = 0))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step12_ribo_rowsum", label = "Filter rowSums", value = 1, min = 0)),
                             column(4, numericInput("in_step12_ribo_colsum", label = "Filter sample ratio", value = 30, min = 0, max = 100)),
                             column(4, checkboxInput("in_step12_ribo_log2", label = "Log2", value = TRUE))
                           ),
                           p("The rowSums and colSums are used to filter the low expression genes and samples.",
                             "rowSums is the sum of each gene, and colSums is the ratio of effective samples.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           textInput("in_step12_ribo_gene_name", label = "Gene name", value = '', placeholder = "eg: gene1,gene2,gene3"),
                           
                           actionButton("act_step12_import_ribo_norm", label = "import the Ribo", icon = icon("file")),
                           tags$style(HTML("#act_step12_import_ribo_norm {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Draw the PCA", status = "primary", collapsed = FALSE, collapsible = TRUE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(3, numericInput("in_step12_ribo_xmin", label = "Xlim down", value = NULL, width = 100)),
                             column(3, numericInput("in_step12_ribo_xmax", label = "Xlim up", value = NULL)),
                             column(3, numericInput("in_step12_ribo_ymin", label = "Ylim down", value = NULL)),
                             column(3, numericInput("in_step12_ribo_ymax", label = "Ylim up", value = NULL))
                           ),
                           
                           fluidRow(
                             column(4, checkboxInput("in_step12_ribo_label_sample", label = "Label sample", value = TRUE)),
                             column(4, sliderInput("in_step12_ribo_label_size", label = "Label size", value = 3, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step12_ribo_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step12_ribo_dot_color", label = "Dot color", choices = color_list, selected = "Paired")),
                             column(4, sliderInput("in_step12_ribo_dot_size", label = "Dot size", value = 3, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step12_ribo_dot_alpha", label = "Alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step12_ribo_width", label = "Figure width", value = 7, min = 1)),
                             column(6, numericInput("out_step12_ribo_height", label = "Figure height", value = 7, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step12_ribo_name", label = "Output", value = "Ribo-seq", placeholder = "file name prefix")),
                             column(6, selectizeInput("out_step12_ribo_format", label = "Format", choices = c(".pdf", ".png"), selected = "pdf"))
                           ),
                           
                           fluidRow(
                             column(4, actionButton("act_step12_draw_ribo_pca", label = "run PCA", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step12_ribo_eigen", label = "Save results", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step12_ribo_eigen_plot", label = "Save eigen plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step12_draw_ribo_pca {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step12_ribo_eigen {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step12_ribo_eigen_plot {background-color: #grey; color: black; margin-top: 0px;}")),
                           
                           fluidRow(
                             column(4, downloadButton("save_step12_ribo_scatter", label = "Save scatter", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step12_ribo_variances", label = "Save variances", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step12_ribo_individual", label = "Save individual", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step12_ribo_scatter {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step12_ribo_variances {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step12_ribo_individual {background-color: #grey; color: black; margin-top: 5px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Abundance",
                                  h4("Show the ribo-seq norm table"),
                                  DT::dataTableOutput("out_step12_ribo_norm")),
                         
                         tabPanel(title = "Eigenvalue",
                                  h4("Show the ribo-seq eigenvalue table"),
                                  DT::dataTableOutput("out_step12_ribo_eigen")),
                         
                         tabPanel(title = "Var.contrib",
                                  h4("Show the ribo-seq contribution table"),
                                  DT::dataTableOutput("out_step12_ribo_var_contrib")),
                         
                         tabPanel(title = "Ind.contrib",
                                  h4("Show the ribo-seq contribution table"),
                                  DT::dataTableOutput("out_step12_ribo_ind_contrib")),
                         
                         ## set the PCA plot
                         tabPanel(title = "EigenPlot",
                                  h4("Show the ribo-seq eigen plot"),
                                  plotOutput("out_step12_ribo_eigen_plot")),
                         
                         tabPanel(title = "Scatter",
                                  h4("Show the ribo-seq scatter plot"),
                                  plotOutput("out_step12_ribo_scatter")),
                         
                         tabPanel(title = "Variances",
                                  h4("Show the ribo-seq variances plot"),
                                  plotOutput("out_step12_ribo_variances")),
                         
                         tabPanel(title = "Individual",
                                  h4("Show the ribo-seq individual plot"),
                                  plotOutput("out_step12_ribo_individual"))
                         
                       )
                )
              )
      ),
      
      ### step11 TE PCA ###############################
      tabItem(tabName = "step12_te",
              h2("show the PCA of TE"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the TE data", status = "primary", collapsed = FALSE, collapsible = TRUE, 
                           solidHeader = TRUE, width = NULL,
                           
                           ## click to import the expression table
                           fileInput(inputId = "in_step12_te", label = "Norm", accept = c(".txt", ".TXT"), multiple = TRUE),
                           
                           fluidRow(
                             column(7, textInput("in_step12_te_group", label = "Group", value = "Group", placeholder = "need design file")),
                             column(5, numericInput("in_step12_te_top", label = "Top number", value = 0))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step12_te_rowsum", label = "Filter rowMeans", value = 0, min = 0)),
                             column(4, numericInput("in_step12_te_colsum", label = "Filter sample ratio", value = 30, min = 0, max = 100)),
                             column(4, checkboxInput("in_step12_te_log2", label = "Log2", value = TRUE))
                           ),
                           p("The rowSums and colSums are used to filter the low expression genes and samples.",
                             "rowSums is the sum of each gene, and colSums is the ratio of effective samples.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           textInput("in_step12_te_gene_name", label = "Gene name", value = '', placeholder = "eg: gene1,gene2,gene3"),
                           
                           actionButton("act_step12_import_te_norm", label = "import the TE", icon = icon("file")),
                           tags$style(HTML("#act_step12_import_te_norm {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Draw the PCA", status = "primary", collapsed = FALSE, collapsible = TRUE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(3, numericInput("in_step12_te_xmin", label = "Xlim down", value = NULL, width = 100)),
                             column(3, numericInput("in_step12_te_xmax", label = "Xlim up", value = NULL)),
                             column(3, numericInput("in_step12_te_ymin", label = "Ylim down", value = NULL)),
                             column(3, numericInput("in_step12_te_ymax", label = "Ylim up", value = NULL))
                           ),
                           
                           fluidRow(
                             column(4, checkboxInput("in_step12_te_label_sample", label = "Label sample", value = TRUE)),
                             column(4, sliderInput("in_step12_te_label_size", label = "Label size", value = 3, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step12_te_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step12_te_dot_color", label = "Dot color", choices = color_list, selected = "Paired")),
                             column(4, sliderInput("in_step12_te_dot_size", label = "Dot size", value = 3, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step12_te_dot_alpha", label = "Alpha", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("out_step12_te_width", label = "Figure width", value = 7, min = 1)),
                             column(6, numericInput("out_step12_te_height", label = "Figure height", value = 7, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step12_te_name", label = "Output", value = "TE", placeholder = "file name prefix")),
                             column(6, selectizeInput("out_step12_te_format", label = "Format", choices = c(".pdf", ".png"), selected = "pdf"))
                           ),
                           
                           fluidRow(
                             column(4, actionButton("act_step12_draw_te_pca", label = "run PCA", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step12_te_eigen", label = "Save results", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step12_te_eigen_plot", label = "Save eigen plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step12_draw_te_pca {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step12_te_eigen {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step12_te_eigen_plot {background-color: #grey; color: black; margin-top: 0px;}")),
                           
                           fluidRow(
                             column(4, downloadButton("save_step12_te_scatter", label = "Save scatter", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step12_te_variances", label = "Save variances", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step12_te_individual", label = "Save individual", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step12_te_scatter {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step12_te_variances {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step12_te_individual {background-color: #grey; color: black; margin-top: 5px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Abundance",
                                  h4("Show the TE norm table"),
                                  DT::dataTableOutput("out_step12_te_norm")),
                         
                         tabPanel(title = "Eigenvalue",
                                  h4("Show the TE eigenvalue table"),
                                  DT::dataTableOutput("out_step12_te_eigen")),
                         
                         tabPanel(title = "Var.contrib",
                                  h4("Show the TE contribution table"),
                                  DT::dataTableOutput("out_step12_te_var_contrib")),
                         
                         tabPanel(title = "Ind.contrib",
                                  h4("Show the TE contribution table"),
                                  DT::dataTableOutput("out_step12_te_ind_contrib")),
                         
                         ## set the PCA plot
                         tabPanel(title = "EigenPlot",
                                  h4("Show the TE eigen plot"),
                                  plotOutput("out_step12_te_eigen_plot")),
                         
                         tabPanel(title = "Scatter",
                                  h4("Show the TE scatter plot"),
                                  plotOutput("out_step12_te_scatter")),
                         
                         tabPanel(title = "Variances",
                                  h4("Show the TE variances plot"),
                                  plotOutput("out_step12_te_variances")),
                         
                         tabPanel(title = "Individual",
                                  h4("Show the TE individual plot"),
                                  plotOutput("out_step12_te_individual"))
                         
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step13 deltaTE ###############################
      
      ### step13 edgeR differential analysis ###############################
      tabItem(tabName = "step13_edger",
              h2("run the edgeR"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the expression", status = "primary",
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, width = NULL,
                           
                           ## click to import the expression table
                           fileInput(inputId = "in_step13_edger", label = "Expression", accept = c(".txt", ".TXT"),
                                     multiple = FALSE, placeholder = "expression data"),
                           
                           fluidRow(
                             column(6, selectizeInput(inputId = "in_step13_edger_type", label = "Type", choices = c("RNA", "RIBO"),
                                                      selected = "RNA")),
                             column(6, selectizeInput("in_step13_edger_group", label = "Group name", choices = NULL, selected = NULL))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step13_edger_group1", label = "Group 1", choices = NULL, selected = NULL)),
                             column(6, selectizeInput("in_step13_edger_group2", label = "Group 2", choices = NULL, selected = NULL))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step13_edger_anno_lab", label = "Gene annotation label",
                                                      choices = NULL, selected = NULL, multiple = FALSE)),
                             column(8, selectizeInput("in_step13_edger_anno_col", label = "Gene annotation column",
                                                      choices = NULL, selected = NULL, multiple = TRUE))
                           ),
                           actionButton("act_step13_import_edger", label = "import the exprs", icon = icon("file")),
                           tags$style(HTML("#act_step13_import_edger {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "edgeR differential analysis", status = "primary",
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, width = NULL,
                           
                           fluidRow(
                             column(6, numericInput("in_step13_edger_rowsum", label = "Filter rowSums", value = 5, min = 0)),
                             column(6, numericInput("in_step13_edger_stdev", label = "Filter st.dev.", value = NULL, min = 0, max = 100))
                           ),
                           p("The rowSums is used to filter the low expression genes and samples, defalut 5.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           p("The st.dev. is Standard deviation used to filter outliers genes in same group, default NULL.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           fluidRow(
                             column(6, numericInput("in_step13_edger_padj", label = "Filter padj", value = 0.05, min = 0, max = 1)),
                             column(6, numericInput("in_step13_edger_log2fc", label = "Filter log2FC", value = 1, min = 0, max = 10))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step13_edger_test", label = "Method", choices = c("LRT", "QLF", "exact"), selected = "LRT")),
                             column(6, numericInput("in_step13_edger_bcv", label = "BCV", value = 0.02, min = 0, max = 1)),
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step13_edger_name", label = "Output xlsx", value = "RNA-seq", placeholder = "file name prefix")),
                             column(6, textInput("out_step13_edger_sheet", label = "Sheet name", value = "rna_DEGs", placeholder = "sheet name of excel file"))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step13_run_edger", label = "run edgeR", icon = icon("calculator"))),
                             column(6, downloadButton("save_step13_edger_res", label = "Save results",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step13_run_edger {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step13_edger_res {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Abundance",
                                  h4("Show the edgeR exprs table"),
                                  DT::dataTableOutput("out_step13_edger_exprs"),
                                  textOutput("out_step13_edger_exprs_info"),
                                  tags$style(HTML("#out_step13_edger_exprs_info {color: red; font-size: 26px;}"))),
                         
                         tabPanel(title = "Design",
                                  h4("Show the DESeq2 design table"),
                                  DT::dataTableOutput("out_step13_edger_design")),
                         
                         tabPanel(title = "Results",
                                  h4("Show the edgeR results table"),
                                  conditionalPanel(condition = ("input.act_step13_run_edger > 0"),
                                                   DT::dataTableOutput("out_step13_edger_res") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         # tabPanel(title = "Down",
                         #          conditionalPanel(condition = ("input.act_step13_run_edger > 0"),
                         #                           DT::dataTableOutput("out_step13_edger_down") %>%
                         #                             withSpinner(color="black", type = 5, size = 0.5)),
                         #          h4("Show the edgeR down-regulated genes")),
                         #
                         # tabPanel(title = "Up",
                         #          conditionalPanel(condition = ("input.act_step13_run_edger > 0"),
                         #                           DT::dataTableOutput("out_step13_edger_up") %>%
                         #                             withSpinner(color="black", type = 5, size = 0.5)),
                         #          h4("Show the edgeR up-regulated genes")),
                         
                         tabPanel(title = "Summary",
                                  h4("Show the edgeR results summary"),
                                  DT::dataTableOutput("out_step13_edger_summary"))
                         
                       )
                )
              )
      ),
      
      ### step13 DESeq2 differential analysis ###############################
      tabItem(tabName = "step13_deseq2",
              h2("run the DESeq2"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the expression", status = "primary",
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, width = NULL,
                           
                           ## click to import the expression table
                           fileInput(inputId = "in_step13_deseq2", label = "Expression", accept = c(".txt", ".TXT"),
                                     multiple = FALSE, placeholder = "expression data"),
                           
                           fluidRow(
                             column(6, selectizeInput(inputId = "in_step13_deseq2_type", label = "Type", choices = c("RNA", "RIBO"),
                                                      selected = "RNA")),
                             column(6, selectizeInput("in_step13_deseq2_group", label = "Group name", choices = NULL, selected = NULL))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step13_deseq2_group1", label = "Group 1", choices = NULL, selected = NULL)),
                             column(6, selectizeInput("in_step13_deseq2_group2", label = "Group 2", choices = NULL, selected = NULL))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step13_deseq2_anno_lab", label = "Gene annotation label",
                                                      choices = NULL, selected = NULL, multiple = FALSE)),
                             column(8, selectizeInput("in_step13_deseq2_anno_col", label = "Gene annotation column",
                                                      choices = NULL, selected = NULL, multiple = TRUE))
                           ),
                           
                           actionButton("act_step13_import_deseq2", label = "import the exprs", icon = icon("file")),
                           tags$style(HTML("#act_step13_import_deseq2 {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "DESeq2 differential analysis", status = "primary",
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, width = NULL,
                           
                           fluidRow(
                             column(6, numericInput("in_step13_deseq2_rowsum", label = "Filter rowSums", value = 5, min = 0)),
                             column(6, numericInput("in_step13_deseq2_stdev", label = "Filter st.dev.", value = NULL, min = 0, max = 100))
                           ),
                           p("The rowSums is used to filter the low expression genes and samples, defalut 5.",
                             "The st.dev. is Standard deviation used to filter outliers genes in same group, default NULL.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           fluidRow(
                             column(4, numericInput("in_step13_deseq2_padj", label = "Filter padj", value = 0.05, min = 0, max = 1)),
                             column(4, numericInput("in_step13_deseq2_log2fc", label = "Filter log2FC", value = 1, min = 0, max = 10)),
                             column(4, selectizeInput("in_step13_deseq2_test", label = "Method", choices = c('Wald', 'LRT'), selected = "Wald")),
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step13_deseq2_name", label = "Output xlsx", value = "RNA-seq", placeholder = "file name prefix")),
                             column(6, textInput("out_step13_deseq2_sheet", label = "Sheet name", value = "rna_DEGs", placeholder = "sheet name of excel file"))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step13_run_deseq2", label = "run DESeq2", icon = icon("calculator"))),
                             column(6, downloadButton("save_step13_deseq2_res", label = "Save results",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step13_run_deseq2 {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step13_deseq2_res {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Abundance",
                                  h4("Show the DESeq2 exprs table"),
                                  DT::dataTableOutput("out_step13_deseq2_exprs"),
                                  textOutput("out_step13_deseq2_exprs_info"),
                                  tags$style(HTML("#out_step13_deseq2_exprs_info {color: red; font-size: 26px;}"))),
                         
                         tabPanel(title = "Design",
                                  h4("Show the DESeq2 design table"),
                                  DT::dataTableOutput("out_step13_deseq2_design")),
                         
                         tabPanel(title = "Results",
                                  h4("Show the DESeq2 results table"),
                                  conditionalPanel(condition = ("input.act_step13_run_deseq2 > 0"),
                                                   DT::dataTableOutput("out_step13_deseq2_res") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         # tabPanel(title = "Down",
                         #          conditionalPanel(condition = ("input.act_step13_run_deseq2 > 0"),
                         #                           DT::dataTableOutput("out_step13_deseq2_down") %>%
                         #                             withSpinner(color="#3c8dbc", type = 5, size = 1)),
                         #          h4("Show the DESeq2 down-regulated genes")),
                         #
                         # tabPanel(title = "Up",
                         #          conditionalPanel(condition = ("input.act_step13_run_deseq2 > 0"),
                         #                           DT::dataTableOutput("out_step13_deseq2_up") %>%
                         #                             withSpinner(color="#3c8dbc", type = 5, size = 1)),
                         #          h4("Show the DESeq2 up-regulated genes")),
                         
                         tabPanel(title = "Summary",
                                  h4("Show the DESeq2 results summary"),
                                  DT::dataTableOutput("out_step13_deseq2_summary"))
                         
                       )
                )
              )
      ),
      
      ### step13 deltaTE differential analysis ###############################
      tabItem(tabName = "step13_delta",
              h2("run the deltaTE"),
              
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       box(title = "Import the RNA-seq and Ribo-seq data", status = "primary", solidHeader = TRUE,
                           collapsed = FALSE, collapsible = TRUE, width = NULL,
                           
                           ## click to import the expression table
                           fileInput(inputId = "in_step13_deltate_rna", label = "RNA", accept = c(".txt", ".TXT"),
                                     multiple = F, placeholder = "RNA exprs"),
                           fileInput(inputId = "in_step13_deltate_ribo", label = "RIBO", accept = c(".txt", ".TXT"),
                                     multiple = F, placeholder = "Ribo exprs"),
                           
                           fluidRow(
                             column(3, selectizeInput("in_step13_deltate_rna_group", label = "RNA group", choices = NULL, selected = NULL)),
                             column(4, selectizeInput("in_step13_deltate_group1", label = "RNA group1", choices = NULL, selected = NULL)),
                             column(4, selectizeInput("in_step13_deltate_group2", label = "RNA group2", choices = NULL, selected = NULL))
                           ),
                           
                           fluidRow(
                             column(3, selectizeInput("in_step13_deltate_ribo_group", label = "Ribo group", choices = NULL, selected = NULL)),
                             column(4, selectizeInput("in_step13_deltate_group3", label = "Ribo group1", choices = NULL, selected = NULL)),
                             column(4, selectizeInput("in_step13_deltate_group4", label = "Ribo group2", choices = NULL, selected = NULL))
                           ),
                           
                           actionButton("act_step13_import_deltate", label = "import the exprs", icon = icon("file")),
                           tags$style(HTML("#act_step13_import_deltate {background-color: #75aadb; color: black;}"))
                           
                       ),
                       
                       box(title = "Run deltaTE", status = "primary", collapsed = FALSE, collapsible = TRUE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step13_deltate_anno_lab", label = "Gene annotation label",
                                                      choices = NULL, selected = NULL)),
                             column(8, selectizeInput("in_step13_deltate_anno_col", label = "Gene annotation column",
                                                      choices = NULL, selected = NULL, multiple = TRUE))
                           ),
                           p("The gene annotation is used to annotate the gene message.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           p("The annotation label is the column name to join the annotation and deltaTE results",
                             "The annotation column is the column name to show in the deltaTE results.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           fluidRow(
                             column(6, numericInput("in_step13_deltate_rowsum", label = "Filter rowSums", value = 5, min = 0)),
                             column(6, numericInput("in_step13_deltate_stdev", label = "Filter st.dev.", value = NULL, min = 0, max = 1e+5))
                           ),
                           p("The rowSums is used to filter the low expression genes and samples, defalut 5.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           p("The variance is used to filter outliers in the same group of samples, default NULL.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           fluidRow(
                             column(6, numericInput("in_step13_deltate_padj", label = "Filter padj", value = 0.05, min = 0, max = 1)),
                             column(6, numericInput("in_step13_deltate_log2fc", label = "Filter log2FC", value = 1, min = 0, max = 10))
                           ),
                           
                           textInput("out_step13_deltate_name", label = "Output xlsx", value = "deltaTE", placeholder = "file name prefix"),
                           
                           fluidRow(
                             column(6, actionButton("act_step13_run_deltate", label = "run deltaTE", icon = icon("calculator"))),
                             column(6, downloadButton("save_step13_deltate_res", label = "Save deltaTE results",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step13_run_deltate {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step13_deltate_res {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "RNA",
                                  h4("Show the RNA-seq exprs table"),
                                  DT::dataTableOutput("out_step13_deltate_rna_exprs"),
                                  textOutput("out_step13_deltate_rna_info"),
                                  tags$style(HTML("#out_step13_deltate_rna_info {color: red; font-size: 26px;}"))),
                         
                         tabPanel(title = "Ribo",
                                  h4("Show the Ribo-seq exprs table"),
                                  DT::dataTableOutput("out_step13_deltate_ribo_exprs"),
                                  textOutput("out_step13_deltate_ribo_info"),
                                  tags$style(HTML("#out_step13_deltate_ribo_info {color: red; font-size: 26px;}"))),
                         
                         tabPanel(title = "Design",
                                  h4("Show the design table"),
                                  DT::dataTableOutput("out_step13_deltate_design")),
                         
                         tabPanel(title = "DTEGs",
                                  h4("Show the deltaTE results table"),
                                  conditionalPanel(condition = ("input.act_step13_run_deltate > 0"),
                                                   DT::dataTableOutput("out_step13_deltate_res") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "RNA-DEGs",
                                  h4("Show the deltaTE down-regulated genes"),
                                  conditionalPanel(condition = ("input.act_step13_run_deltate > 0"),
                                                   DT::dataTableOutput("out_step13_deltate_rna_degs") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Ribo-DEGs",
                                  h4("Show the deltaTE up-regulated genes"),
                                  conditionalPanel(condition = ("input.act_step13_run_deltate > 0"),
                                                   DT::dataTableOutput("out_step13_deltate_ribo_degs") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Summary",
                                  h4("Show the deltaTE results summary"),
                                  conditionalPanel(condition = ("input.act_step13_run_deltate > 0"),
                                                   DT::dataTableOutput("out_step13_deltate_summary") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                         
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step14 volcano ###############################
      
      ### step14 draw the volcano of RNA-seq / Ribo-seq / TE ###############################
      tabItem(tabName = "step14_rna_ribo",
              h2("draw the volcano plot of RNA-seq / Ribo-seq / TE"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the expression table
                       box(title = "Import the DEGs table", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step14_vol", label = "DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = FALSE, placeholder = "DEGs table")),
                             column(4, selectizeInput(inputId = "in_step14_vol_sheet", label = "Sheet", choices = c("DTEGs", "rna_DEGs", "ribo_DEGs"),
                                                      selected = "rna_DEGs"))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput(inputId = "in_step14_vol_type", label = "Method", choices = c("DESeq2", "edgeR", "deltaTE"),
                                                      selected = "DESeq2")),
                             conditionalPanel(
                               condition = "input.in_step14_vol_type == 'DESeq2'",
                               column(6, textInput("in_step14_vol_title", label = "Title", value = "DESeq2-volcano", placeholder = "figure title"))
                             ),
                             conditionalPanel(
                               condition = "input.in_step14_vol_type == 'edgeR'",
                               column(6, textInput("in_step14_vol_title", label = "Title", value = "edgeR-volcano", placeholder = "figure title"))
                             ),
                             conditionalPanel(
                               condition = "input.in_step14_vol_type == 'deltaTE'",
                               column(6, textInput("in_step14_vol_title", label = "Title", value = "deltaTE-volcano", placeholder = "figure title"))
                             )
                           ),
                           
                           fluidRow(
                             conditionalPanel(
                               condition = "input.in_step14_vol_type == 'DESeq2'",
                               column(4, textInput("in_step14_vol_x", label = "X-axis", value = "log2FoldChange", placeholder = "x-axis data")),
                               column(4, textInput("in_step14_vol_y", label = "Y-axis", value = "padj", placeholder = "y-axis data")),
                               column(4, textInput("in_step14_vol_basemean", label = "BaseMean", value = "baseMean", placeholder = "x-axis data"))
                             ),
                             
                             conditionalPanel(
                               condition = "input.in_step14_vol_type == 'edgeR'",
                               column(4, textInput("in_step14_vol_x", label = "X-axis", value = "logFC", placeholder = "x-axis data")),
                               column(4, textInput("in_step14_vol_y", label = "Y-axis", value = "FDR", placeholder = "y-axis data")),
                               column(4, textInput("in_step14_vol_basemean", label = "BaseMean", value = "logCPM", placeholder = "x-axis data"))
                             ),
                             
                             conditionalPanel(
                               condition = "input.in_step14_vol_type == 'deltaTE'",
                               column(4, textInput("in_step14_vol_x", label = "X-axis", value = "log2FoldChange", placeholder = "x-axis data")),
                               column(4, textInput("in_step14_vol_y", label = "Y-axis", value = "padj", placeholder = "y-axis data")),
                               column(4, textInput("in_step14_vol_basemean", label = "BaseMean", value = "baseMean", placeholder = "y-axis data"))
                             )
                           ),
                           p("The DEGs table should contain the gene ID, log2FC, padj, and other information.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           p("The terms 'Padj' and 'LogFC' in this context correspond to the data used for plotting in the table..",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           actionButton("act_step14_import_vol_degs", label = "import the DEGs", icon = icon("file")),
                           tags$style(HTML("#act_step14_import_vol_degs {background-color: #75aadb; color: black;}"))
                           
                       ),
                       
                       box(title = "Draw the volcano", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(6, textInput("in_step14_vol_xlab", label = "Volcano-Xlab", value = "log2FC", placeholder = "label for x-axis")),
                             column(6, textInput("in_step14_vol_ylab", label = "Volcano-Ylab", value = "-log10(Padj)", placeholder = "label for y-axis"))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step14_vol_ma_xlab", label = "MAplot-Xlab", value = "baseMean", placeholder = "label for x-axis")),
                             column(6, textInput("in_step14_vol_ma_ylab", label = "MAplot-Ylab", value = "log2FC", placeholder = "label for y-axis"))
                           ),
                           
                           fluidRow(
                             column(3, numericInput("in_step14_vol_padj", label = "Padj-min", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(3, numericInput("in_step14_vol_logfc", label = "LogFC-min", value = 1, min = 0, max = 20, step = 0.1)),
                             column(3, numericInput("in_step14_vol_xlim_min", label = "X-min", value = NULL, min = -20, max = 20, step = 0.5)),
                             column(3, numericInput("in_step14_vol_xlim_max", label = "X-max", value = NULL, min = -20, max = 20, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step14_vol_down_color", label = "DOWN color", value = "#2b78e4", placeholder = "eg. #2b78e4 or blue")),
                             column(4, textInput("in_step14_vol_ns_color", label = "NS color", value = "grey", placeholder = "eg. #bebebe or grey")),
                             column(4, textInput("in_step14_vol_up_color", label = "UP color", value = "#ea6e3c", placeholder = "eg. #ea6e3c or red"))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step14_vol_dot_size", label = "Dot size", value = 1.5, min = 0, max = 10, step = 0.5)),
                             column(6, sliderInput("in_step14_vol_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(4, checkboxInput("in_step14_vol_sqrt", label = "Sqrt", value = FALSE)),
                             column(4, checkboxInput("in_step14_vol_log2", label = "Log", value = FALSE)),
                             column(4, checkboxInput("in_step14_vol_pvalue0", label = "Pvalue 0", value = TRUE))
                           ),
                           p("The Sqrt and Log2 are used to change the y-axis limits. Pvalue 0 is used to convert the 0 value to the smallest value.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           fluidRow(
                             column(4, textInput("in_step14_vol_label_col", label = "Label column", value = "Gene", placeholder = "column for gene label")),
                             column(4, selectizeInput("in_step14_vol_label_type", label = "Label type", selected = "text",
                                                      choices = c("text", "dot"))),
                             column(4, textInput("in_step14_vol_label_color", label = "Label color", value = "black"))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step14_vol_label_num", label = "Label number", value = 0, min = 0)),
                             column(4, numericInput("in_step14_vol_overlaps_num", label = "Max overlaps", value = 10, min = 0)),
                             column(4, numericInput("in_step14_vol_label_size", label = "Label size", value = 3, min = 1))
                           ),
                           
                           textAreaInput("in_step14_vol_label_degs", label = "Label genes", value = NULL, placeholder = "eg. gene1, gene2, gene3"),
                           
                           fluidRow(
                             column(6, textInput("out_step14_vol_out_name", label = "Output plot", value = "RNA-seq", placeholder = "file name prefix")),
                             column(3, numericInput("out_step14_vol_width", label = "Figure width", value = 6, min = 1)),
                             column(3, numericInput("out_step14_vol_height", label = "Figure height", value = 6, min = 1))
                           ),
                           
                           fluidRow(
                             column(4, actionButton("act_step14_draw_volcano", label = "draw volcano", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step14_vol_plot", label = "Save volcano", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step14_vol_maplot", label = "Save MAplot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step14_draw_volcano {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step14_vol_plot {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step14_vol_maplot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "DEGs",
                                  h4("Show the DEGs table"),
                                  DT::dataTableOutput("out_step14_vol_degs")),
                         
                         tabPanel(title = "Volcano",
                                  h4("Show the DEGs volcano."),
                                  plotOutput("out_step14_vol_plot")),
                         
                         tabPanel(title = "MAplot",
                                  h4("Show the DEGs MAplot."),
                                  plotOutput("out_step14_vol_maplot"))
                         
                       )
                )
              )
      ),
      
      ### step14 draw the quadrant plot of RNA-seq / Ribo-seq / TE ###############################
      tabItem(tabName = "step14_quadrant",
              h2("draw the quadrant plot of RNA-seq / Ribo-seq / TE"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the expression table
                       box(title = "Import the DEGs table", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step14_quadra_degs_x", label = "RNA/Ribo DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = FALSE, placeholder = "DEGs table")),
                             column(4, selectizeInput(inputId = "in_step14_quadra_degs_x_sheet", label = "Sheet", selected = "rna_DEGs",
                                                      choices = c("rna_DEGs", "ribo_DEGs")))
                           ),
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step14_quadra_degs_y", label = "Ribo/TE DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = FALSE, placeholder = "DEGs table")),
                             column(4, selectizeInput(inputId = "in_step14_quadra_degs_y_sheet", label = "Sheet", selected = "ribo_DEGs",
                                                      choices = c("ribo_DEGs", "DTEGs")))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput(inputId = "in_step14_quadra_type", label = "Method", choices = c("DESeq2", "edgeR", "deltaTE"),
                                                      selected = "DESeq2")),
                             
                             conditionalPanel(
                               condition = "input.in_step14_quadra_type == 'DESeq2'",
                               column(4, textInput("in_step14_quadra_y", label = "pvalue name", value = "padj", placeholder = "column for pvalue")),
                               column(4, textInput("in_step14_quadra_x", label = "log2FC name", value = "log2FoldChange", placeholder = "column for log2FC"))
                             ),
                             
                             conditionalPanel(
                               condition = "input.in_step14_quadra_type == 'edgeR'",
                               column(4, textInput("in_step14_quadra_y", label = "pvalue name", value = "FDR", placeholder = "column for pvalue")),
                               column(4, textInput("in_step14_quadra_x", label = "log2FC name", value = "logFC", placeholder = "column for log2FC"))
                             ),
                             
                             conditionalPanel(
                               condition = "input.in_step14_quadra_type == 'deltaTE'",
                               column(4, textInput("in_step14_quadra_y", label = "pvalue name", value = "padj", placeholder = "column for pvalue")),
                               column(4, textInput("in_step14_quadra_x", label = "log2FC name", value = "log2FoldChange", placeholder = "column for log2FC"))
                             )
                           ),
                           p("The DEGs table should contain the gene ID, log2FC, padj, and other information.",
                             style = "color: #7F7F7F; font-size: 13px; font-style: Italic"),
                           
                           fluidRow(
                             column(6, actionButton("act_step14_import_quadra_degs", label = "Import the DEGs", icon = icon("file"))),
                             column(6, downloadButton("save_step14_quadra_degs", label = "Save the DEGs", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step14_import_quadra_degs {background-color: #75aadb; color: black;}")),
                           tags$style(HTML("#save_step14_quadra_degs {background-color: #75aadb; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Draw nine-quadrant plot", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, numericInput("in_step14_quadra_padj", label = "Padj", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step14_quadra_logfc", label = "LogFC (abs)", value = 1, min = 0, max = 20, step = 0.1)),
                             column(4, checkboxInput("in_step14_quadra_ns", label = "Remove NS", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(3, numericInput("in_step14_quadra_xlim_min", label = "X-min", value = NULL, min = -20, max = 20, step = 0.5)),
                             column(3, numericInput("in_step14_quadra_xlim_max", label = "X-max", value = NULL, min = -20, max = 20, step = 0.5)),
                             column(3, numericInput("in_step14_quadra_ylim_min", label = "Y-min", value = NULL, min = -20, max = 20, step = 0.5)),
                             column(3, numericInput("in_step14_quadra_ylim_max", label = "Y-max", value = NULL, min = -20, max = 20, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step14_quadra_dot_size", label = "Dot size", value = 1.5, min = 0, max = 10, step = 0.5)),
                             column(4, selectizeInput("in_step14_quadra_color", label = "Dot color", choices = color_list, selected = "Set2")),
                             column(4, numericInput("in_step14_quadra_font_size", label = "Font size", value = 12, min = 0, max = 20, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step14_quadra_label_col", label = "Label column", value = NULL, placeholder = "column for gene label")),
                             column(4, textInput("in_step14_quadra_label_color", label = "Label color", value = "black"))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step14_quadra_label_size", label = "Label size", value = 3, min = 0, max = 10, step = 0.5)),
                             column(4, numericInput("in_step14_quadra_overlaps_num", label = "Max overlaps", value = 10, min = 0, step = 1))
                           ),
                           
                           textAreaInput("in_step14_quadra_label_degs", label = "Label genes", value = NULL, placeholder = "eg. gene1, gene2, gene3"),
                           
                           fluidRow(
                             column(4, textInput("in_step14_quadra_title", label = "Title", value = "", placeholder = "DEGs")),
                             column(4, textInput("in_step14_quadra_xlabel", label = "Xlabel", value = "RNA_log2FC", placeholder = "RNA_log2FC")),
                             column(4, textInput("in_step14_quadra_ylabel", label = "Ylabel", value = "Ribo_log2FC", placeholder = "Ribo_log2FC"))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step14_quadra_out_name", label = "Output plot", value = "quadrants", placeholder = "file name prefix")),
                             column(3, numericInput("out_step14_quadra_width", label = "Figure width", value = 6, min = 1)),
                             column(3, numericInput("out_step14_quadra_height", label = "Figure height", value = 6, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step14_draw_quadra_plot", label = "draw quadrant", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step14_quadra_plot", label = "Save quadrant",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step14_draw_quadra_plot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step14_quadra_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "DEGs",
                                  h4("Show the DEGs table"),
                                  DT::dataTableOutput("out_step14_quadra_degs")),
                         
                         tabPanel(title = "Quadrant",
                                  h4("Show the DEGs quadrant plot"),
                                  plotOutput("out_step14_quadra_plot"))
                         
                       )
                )
              )
      ),
      
      ### step14 draw the delta quadrant plot ###############################
      tabItem(tabName = "step14_deltate",
              h2("draw the deltaTE plot"),
              fluidRow(
                tags$style(HTML(".download-btn {margin-top: 250px;}",
                                ".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the degs table
                       box(title = "Import the DEGs table", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step14_delta_rna", label = "RNA DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = FALSE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step14_delta_rna_sheet", label = "Sheet", value = "rna_DEGs"))
                           ),
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step14_delta_ribo", label = "Ribo DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = FALSE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step14_delta_ribo_sheet", label = "Sheet", value = "ribo_DEGs"))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step14_import_delta_degs", label = "Import the DEGs", icon = icon("file"))),
                             column(6, downloadButton("save_step14_delta_degs", label = "Save the DEGs", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step14_import_delta_degs {background-color: #75aadb; color: black;}")),
                           tags$style(HTML("#save_step14_delta_degs {background-color: #75aadb; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Draw nine-quadrant plot", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(6, textInput("in_step14_delta_column", label = "Delta column", value = "Delta", placeholder = "column for deltaTE")),
                             column(6, textInput("in_step14_delta_gene_column", label = "Gene column", value = "Gene", placeholder = "column for GeneID"))
                           ),
                           
                           fluidRow(
                             column(6, checkboxInput("in_step14_delta_others", label = "remove others", value = FALSE)),
                             column(6, checkboxInput("in_step14_delta_ns", label = "remove ns", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step14_delta_logfc", label = "log2FC", value = 1, min = 0, step = 0.5)),
                             column(6, numericInput("in_step14_delta_padj", label = "Padj", value = 0.05, min = 0, max = 1, step = 0.01))
                           ),
                           
                           fluidRow(
                             column(3, numericInput("in_step14_delta_xlim_min", label = "Xlim min", value = NULL, min = -20, max = 20, step = 0.5)),
                             column(3, numericInput("in_step14_delta_xlim_max", label = "Xlim max", value = NULL, min = -20, max = 20, step = 0.5)),
                             column(3, numericInput("in_step14_delta_ylim_min", label = "Ylim min", value = NULL, min = -20, max = 20, step = 0.5)),
                             column(3, numericInput("in_step14_delta_ylim_max", label = "Ylim max", value = NULL, min = -20, max = 20, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step14_delta_dot_size", label = "Dot size", value = 1.5, min = 0, max = 10, step = 0.5)),
                             column(4, selectizeInput("in_step14_delta_color", label = "Dot color", choices = color_list, selected = "Set2")),
                             column(4, numericInput("in_step14_delta_font_size", label = "Font size", value = 12, min = 0, max = 20, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step14_delta_label_col", label = "Label column", value = NULL, placeholder = "column for gene label")),
                             column(4, textInput("in_step14_delta_label_color", label = "Label color", value = "black"))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step14_delta_label_size", label = "Label size", value = 3, min = 0, max = 10, step = 0.5)),
                             column(4, numericInput("in_step14_delta_overlaps_num", label = "Max overlaps", value = 10, min = 0, step = 1))
                           ),
                           
                           textAreaInput("in_step14_delta_label_degs", label = "Label genes", value = NULL, placeholder = "eg. gene1, gene2, gene3"),
                           
                           fluidRow(
                             column(4, textInput("in_step14_delta_title", label = "Title", value = "", placeholder = "DEGs")),
                             column(4, textInput("in_step14_delta_xlabel", label = "Xlabel", value = "RNA_log2FC", placeholder = "RNA_log2FC")),
                             column(4, textInput("in_step14_delta_ylabel", label = "Ylabel", value = "Ribo_log2FC", placeholder = "Ribo_log2FC"))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step14_delta_out_name", label = "Output plot", value = "deltaTE", placeholder = "file name prefix")),
                             column(3, numericInput("out_step14_delta_width", label = "Figure width", value = 6, min = 1)),
                             column(3, numericInput("out_step14_delta_height", label = "Figure height", value = 6, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step14_draw_delta_plot", label = "draw deltaTE", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step14_delta_plot", label = "Save deltaTE",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step14_draw_delta_plot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step14_delta_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "DEGs",
                                  h4("Show the DEGs table"),
                                  DT::dataTableOutput("out_step14_delta_degs")),
                         
                         tabPanel(title = "deltaTE",
                                  h4("Show the deltaTE plot."),
                                  plotOutput("out_step14_delta_plot"))
                       )
                )
              )
      ),
      
      
      ############################################################################################################
      ## step15 DEGs ###############################
      
      ### step15 venn ###############################
      tabItem(tabName = "step15_venn",
              h2("draw the venn diagram"),
              fluidRow(
                tags$style(HTML(".download-btn {margin-top: 250px;}",
                                ".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                tags$style(HTML("#in_step15_venn_set_1_title {color: red; background-color: #f0f0f0; border-color: red;}")),
                tags$style(HTML("#in_step15_venn_set_2_title {color: blue; background-color: #f0f0f0; border-color: blue;}")),
                tags$style(HTML("#in_step15_venn_set_3_title {color: green; background-color: #f0f0f0; border-color: green;}")),
                tags$style(HTML("#in_step15_venn_set_4_title {color: purple; background-color: #f0f0f0; border-color: purple;}")),
                
                # set the input panel
                column(5,
                       ## click to import the set
                       
                       box(title = "Import the data sets", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(6, textInput(inputId = "in_step15_venn_set_1_title", label = "Set 1 title", value = "Set1", 
                                                 placeholder = 'Set1')),
                             column(6, textInput(inputId = "in_step15_venn_set_2_title", label = "Set 2 title", value = "Set2",
                                                 placeholder = 'Set2'))
                           ),
                           
                           fluidRow(
                             column(6, textAreaInput(inputId = "in_step15_venn_set_1", label = "Set 1", value = "", height = "100px", 
                                           placeholder = "gene1\ngene2\ngene3")),
                             column(6, textAreaInput(inputId = "in_step15_venn_set_2", label = "Set 2", value = "", height = "100px", 
                                           placeholder = "gene1\ngene2\ngene3"))
                           ),
                           
                           fluidRow(
                             column(6, textInput(inputId = "in_step15_venn_set_3_title", label = "Set 3 title", value = "Set3", 
                                                 placeholder = 'Set3')),
                             column(6, textInput(inputId = "in_step15_venn_set_4_title", label = "Set 4 title", value = "Set4", 
                                                 placeholder = 'Set4'))
                           ),
                           
                           fluidRow(
                             column(6, textAreaInput(inputId = "in_step15_venn_set_3", label = "Set 3", value = "", height = "100px", 
                                           placeholder = "gene1\ngene2\ngene3")),
                             column(6, textAreaInput(inputId = "in_step15_venn_set_4", label = "Set 4", value = "", height = "100px", 
                                           placeholder = "gene1\ngene2\ngene3"))
                           ),
                           
                           fluidRow(
                             column(4, actionButton("act_step15_venn_calc", label = "Overlap sets", icon = icon("calculator"))),
                             column(4, downloadButton("save_step15_venn_set", label = "Save venn sets",
                                                   class = "download-btn", icon = icon("download"))),
                           ),
                           
                           fluidRow(
                             column(4, tags$style(HTML("#act_step15_venn_calc {background-color: #e78c45; color: black; margin-top: 5px;}"))),
                             column(4, tags$style(HTML("#save_step15_venn_set {background-color: #grey; color: black; margin-top: 5px;}")))
                           ),
                       ),
                       
                       
                       box(title = "Set the options for venn diagram", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(3, textInput("in_step15_venn_set_1_color", label = "Set1 color", value = "#0022cc")),
                             column(3, textInput("in_step15_venn_set_2_color", label = "Set2 color", value = "#22cc00")),
                             column(3, textInput("in_step15_venn_set_3_color", label = "Set3 color", value = "#cc2200")),
                             column(3, textInput("in_step15_venn_set_4_color", label = "Set4 color", value = "#aa22dd"))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step15_venn_font_size", label = "Font size", value = 3, min = 0, max = 10, step = 0.5)),
                             column(6, numericInput("in_step15_venn_number_size", label = "Number size", value = 2, min = 0, max = 10, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step15_venn_fontface", label = "Bold", choices = c("plain", "bold"), selected = "plain")),
                             column(6, selectizeInput("in_step15_venn_fontfamily", label = "Font family", choices = c("sans", "serif", "mono"), selected = "sans"))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step15_venn_line_width", label = "Line width", value = 1, step = 0.1)),
                             column(6, numericInput("in_step15_venn_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step15_venn_label_dist", label = "Label dist", value = 0.1, min = 0, max = 1, step = 0.1)),
                             column(6, checkboxInput("in_step15_venn_scale", label = "Scaled", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step15_venn_width", label = "Width", value = 6, min = 0, max = 20, step = 1)),
                             column(6, numericInput("in_step15_venn_height", label = "Length", value = 6, min = 0, max = 20, step = 1))
                           ),
                           
                           textInput("out_step15_venn_set_title", label = "Output", value = "venn-diagram", placeholder = "venn-diagram"),
                           
                           fluidRow(
                             column(4, actionButton("act_step15_venn_diagram", label = "Draw venn diagram", 
                                                      icon = icon("calculator"))),
                             column(4, downloadButton("save_step15_venn_diagram", label = "Save venn diagram",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step15_venn_diagram {background-color: #e78c45; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step15_venn_diagram {background-color: #grey; color: black; margin-top: 5px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "DataSets",
                                  h4("Show the data sets"),
                                  DT::dataTableOutput("out_step15_venn_data_sets")),
                         
                         tabPanel(title = "Venn diagram",
                                  h4("Show the venn diagram"),
                                  plotOutput("out_step15_venn_diagram"))
                       )
                )
              )
      ),
      
      
      ### step15 degs venn ###############################
      tabItem(tabName = "step15_degs_venn",
              h2("draw the DEGs venn plot"),
              fluidRow(
                tags$style(HTML(".download-btn {margin-top: 250px;}",
                                ".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the degs table
                       
                       box(title = "Import the DEGs table", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step15_degs_venn_rna", label = "RNA DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = FALSE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step15_degs_venn_rna_sheet", label = "Sheet", value = "rna_DEGs"))
                           ),
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step15_degs_venn_ribo", label = "Ribo DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = FALSE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step15_degs_venn_ribo_sheet", label = "Sheet", value = "ribo_DEGs"))
                           ),
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step15_degs_venn_te", label = "TE DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = FALSE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step15_degs_venn_te_sheet", label = "Sheet", value = "DTEGs"))
                           ),
                           
                           actionButton("act_step15_degs_venn_import", label = "import the DEGs", icon = icon("file")),
                           tags$style(HTML("#act_step15_degs_venn_import {background-color: #75aadb; color: black;}"))
                           
                       ),
                       
                       box(title = "Set the options for venn plot", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(6, textInput("in_step15_degs_venn_gene", label = "Gene column", value = "Gene")),
                             column(6, textInput("in_step15_degs_venn_degs", label = "DEGs column", value = "DEGs"))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step15_degs_venn_set_1", label = "Set1", value = "RNA", placeholder = "name of the set")),
                             column(4, textInput("in_step15_degs_venn_set_2", label = "Set2", value = "Ribo", placeholder = "name of the set")),
                             column(4, textInput("in_step15_degs_venn_set_3", label = "Set3", value = "TE", placeholder = "name of the set"))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step15_degs_venn_set_size", label = "Set size", value = 6, min = 0, max = 20, step = 0.5)),
                             column(6, textInput("in_step15_degs_venn_set_color", label = "Set color", value = "black"))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step15_degs_venn_low_color", label = "Low color", value = "white", placeholder = "fill the plot")),
                             column(6, textInput("in_step15_degs_venn_high_color", label = "High color", value = "red", placeholder = "fill the plot"))
                           ),
                           
                           fluidRow(
                             column(3, selectizeInput("in_step15_degs_venn_label", label = "Label", choices = c("both", "count", "percent", "none"), selected = "both")),
                             column(3, selectizeInput("in_step15_degs_venn_label_geom", label = "Label geom", choices = c("label", "text"), selected = "text")),
                             column(3, numericInput("in_step15_degs_venn_label_alpha", label = "Label alpha", value = 1, min = 0, max = 1, step = 0.1)),
                             column(3, textInput("in_step15_degs_venn_label_color", label = "Label color", value = "black"))
                           ),
                           
                           fluidRow(
                             column(3, numericInput("in_step15_degs_venn_label_size", label = "Label size", value = 8, min = 0, max = 20, step = 0.5)),
                             column(3, numericInput("in_step15_degs_venn_label_width", label = "Label width", value = 40, min = 0, max = 100, step = 1)),
                             column(3, selectizeInput("in_step15_degs_venn_edge_lty", label = "Edge lty", choices = c("dashed", "dotted", "solid"), selected = "solid")),
                             column(3, numericInput("in_step15_degs_venn_edge_size", label = "Edge size", value = 0.8, min = 0, max = 5, step = 0.1))
                           )
                       ),
                       
                       box(title = "Set the options for upset plot", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, numericInput("in_step15_degs_venn_set_num", label = "Set number", value = 40, min = 0, max = 40, step = 1)),
                             column(4, sliderInput("in_step15_degs_venn_text_scale", label = "Text scale", value = 2, min = 1, max = 6, step = 0.1)),
                             column(4, sliderInput("in_step15_degs_venn_mb_ratio", label = "Mb ratio", value = 0.6, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step15_degs_venn_up_color", label = "Top bar color", value = '#008f00')),
                             column(6, selectizeInput("in_step15_degs_venn_left_color", label = "Left bar color", choices = color_list, selected = "Paired"))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step15_degs_venn_dot_color", label = "Dot color", value = '#00008f')),
                             column(4, numericInput("in_step15_degs_venn_dot_size", label = "Dot size", value = 4, min = 0, max = 10, step = 0.1)),
                             column(4, numericInput("in_step15_degs_venn_line_size", label = "line size", value = 0.7, min = 0, max = 5, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step15_degs_venn_out_name", label = "Output plot", value = "venn", placeholder = "file name prefix")),
                             column(3, numericInput("out_step15_degs_venn_width", label = "Figure width", value = 6, min = 1)),
                             column(3, numericInput("out_step15_degs_venn_height", label = "Figure height", value = 6, min = 1))
                           ),
                           
                           actionButton("act_step15_draw_venn_plot", label = "draw venn", icon = icon("chart-line")),
                           tags$style(HTML("#act_step15_draw_venn_plot {background-color: #e78c45; color: black;}")),
                           
                           fluidRow(
                             column(4, downloadButton("save_step15_down_venn_plot", label = "Save down venn",
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step15_up_venn_plot", label = "Save up venn",
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step15_degs_venn_plot", label = "Save UpSet",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step15_down_venn_plot {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step15_up_venn_plot {background-color: #grey; color: black; margin-top: 5px;}")),
                           tags$style(HTML("#save_step15_degs_venn_plot {background-color: #grey; color: black; margin-top: 5px;}"))
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "DEGs",
                                  h4("Show the DEGs table"),
                                  DT::dataTableOutput("out_step15_degs_venn_degs")),
                         
                         tabPanel(title = "DOWN",
                                  h4("Show the down venn plot."),
                                  plotOutput("out_step15_down_venn_plot")),
                         
                         tabPanel(title = "UP",
                                  h4("Show the up venn plot."),
                                  plotOutput("out_step15_up_venn_plot")),
                         
                         tabPanel(title = "UpSet",
                                  h4("Show the venn plot."),
                                  plotOutput("out_step15_degs_venn_plot"))
                       )
                )
              )
      ),
      
      ### step15 merge ###############################
      tabItem(tabName = "step15_merge",
              h2("Merge all DEGs"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the degs table
                       
                       box(title = "Set the DEGs table", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step15_merge_rna", label = "RNA DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = TRUE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step15_merge_rna_sheet", label = "Sheet", value = "rna_DEGs"))
                           ),
                           
                           fluidRow(
                             column(6, textAreaInput("in_step15_merge_rna_file", label = "File name", value = "", height = "100px")),
                             column(6, textAreaInput("in_step15_merge_rna_name", label = "Rename", value = "", height = "100px"))
                           ),
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step15_merge_ribo", label = "Ribo DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = TRUE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step15_merge_ribo_sheet", label = "Sheet", value = "ribo_DEGs"))
                           ),
                           
                           fluidRow(
                             column(6, textAreaInput("in_step15_merge_ribo_file", label = "File name", value = "", height = "100px")),
                             column(6, textAreaInput("in_step15_merge_ribo_name", label = "Rename", value = "", height = "100px"))
                           ),
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step15_merge_te", label = "TE DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = TRUE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step15_merge_te_sheet", label = "Sheet", value = "DTEGs"))
                           ),
                           
                           fluidRow(
                             column(6, textAreaInput("in_step15_merge_te_file", label = "File name", value = "", height = "100px")),
                             column(6, textAreaInput("in_step15_merge_te_name", label = "Rename", value = "", height = "100px"))
                           )
                       ),
                       
                       box(title = "Merge all DEGs", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(6, textAreaInput("in_step15_merge_degs_column", label = "DEGs column",
                                                     value = "Gene\nbaseMean\nlog2FoldChange\npvalue\npadj\nDEGs",
                                                     width = "100%", height = "150px")),
                             column(6, textAreaInput("in_step15_merge_anno_column", label = "Annotation column",
                                                     value = "GeneID\nSymbol\nEnsembl\nName",
                                                     width = "100%", height = "150px"))
                           ),
                           p("Select the column names of the DEGs table to merge."),
                           
                           fluidRow(
                             column(8, textInput("out_step15_merge_out_name", label = "Output", value = "merge", placeholder = "file name prefix")),
                             column(4, checkboxInput("in_step15_merge_degs_sig", label = "Significant", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step15_merge_degs", label = "merge DEGs", icon = icon("file"))),
                             column(6, downloadButton("save_step15_merge_degs", label = "Save DEGs",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step15_merge_degs {background-color: #75aadb; color: black;}")),
                           tags$style(HTML("#save_step15_merge_degs {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Merge DEGs",
                                  h4("Show the merged DEGs table"),
                                  conditionalPanel(condition = ("input.act_step15_merge_degs > 0"),
                                                   DT::dataTableOutput("out_step15_merge_degs") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      ### step15 degs bar plot ###############################
      tabItem(tabName = "step15_barplot",
              h2("Darw DEGs barplot"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the degs table
                       
                       box(title = "Import the DEGs table", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step15_merge_degs", label = "Merged DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = TRUE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step15_merge_degs_sheet", label = "Sheet", value = "DEGs"))
                           ),
                           
                           actionButton("act_step15_import_merge_degs", label = "Import DEGs", icon = icon("file")),
                           tags$style(HTML("#act_step15_import_merge_degs {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Draw DEGs", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(6, textInput("in_step15_bar_x", label = "x column", value = "Files")),
                             column(6, textInput("in_step15_bar_y", label = "y column", value = "Count"))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step15_bar_degs", label = "DEGs column", value = "DEGs")),
                             column(6, textInput("in_step15_bar_group", label = "Group column", value = "Groups"))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step15_bar_down_color", label = "Down color", value = "#2b78e4", placeholder = "color")),
                             column(6, textInput("in_step15_bar_up_color", label = "Up color", value = "#ea6e3c", placeholder = "color"))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step15_bar_width", label = "Bar width", value = 0.8, min = 0, max = 1, step = 0.1)),
                             column(4, numericInput("in_step15_bar_alpha", label = "Alpha", value = 1, min = 0, max = 1, step = 0.1)),
                             column(4, checkboxInput("in_step15_bar_label", label = "Label", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step15_bar_font_size", label = "Font size", value = 12, min = 0, max = 20, step = 1)),
                             column(6, sliderInput("in_step15_bar_label_size", label = "Label size", value = 3, min = 0, max = 10, step = 1)),
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step15_bar_out_name", label = "Output", value = "merge", placeholder = "file name prefix")),
                             column(4, numericInput("out_step15_bar_width", label = "Figure width", value = 6, min = 1)),
                             column(4, numericInput("out_step15_bar_height", label = "Figure height", value = 6, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step15_draw_degs", label = "draw DEGs", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step15_degs_bar", label = "Save bar plot",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step15_draw_degs {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step15_degs_bar {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "DEGs",
                                  h4("Show the merged DEGs table"),
                                  DT::dataTableOutput("out_step15_import_degs")),
                         
                         tabPanel(title = "Barplot",
                                  h4("Show the DEGs barplot"),
                                  plotOutput("out_step15_degs_bar"))
                       )
                )
              )
      ),
      
      
      ############################################################################################################
      ## step16 GO ###############################
      ### step 16 enrich go ###############################
      tabItem(tabName = "step16_enrich_go",
              h2("Enrich the GO terms"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the degs table
                       
                       box(title = "Import the DEGs table", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step16_go_degs", label = "DEGs", accept = c(".xlsx", ".XLSX"), 
                                                 multiple = TRUE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step16_go_degs_sheet", label = "Sheet", value = "rna_DEGs"))
                           ),
                           
                           actionButton("act_step16_import_go_degs", label = "import the DEGs", icon = icon("file")),
                           tags$style(HTML("#act_step16_import_go_degs {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Enrich GO", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectInput("in_step16_go_database", label = "Database", choices = c("gson", "orgdb"), selected = "orgdb")),
                             conditionalPanel(condition = "input.in_step16_go_database == 'gson'",
                                              column(8, fileInput("in_step16_go_group_gson", label = "Gson file", accept = c(".gson", ".GSON"), multiple = FALSE))
                             ),
                             conditionalPanel(condition = "input.in_step16_go_database == 'orgdb'",
                                              column(8, selectizeInput("in_step16_go_orgdb", label = "Orgdb name", selected = "", choices = ""))
                             )
                           ),
                           p("The orgdb and gson are the name of the database for enrichment. 
                      Select one database of them is satisfied."),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step16_go_key_type", label = "Key Type", selected = "GID", 
                                                      choices = c("GID", "ENTREZID", "ENSEMBL"))),
                             column(4, textInput("in_step16_go_gene_column", label = "Gene column", value = "GeneID", placeholder = "gene column from enrich"))
                           ),
                           p("The Key Type is the IDs from OrgDb for GO enrichment."),
                           p("The Gene column is the IDs from DEGs table corresponding to the Key Type."),
                           
                           fluidRow(
                             column(4, textInput("in_step16_go_degs_column", label = "DEGs column", value = "DEGs", placeholder = "DEGs column from DEGs")),
                             column(4, checkboxInput("in_step16_go_degs_group", label = "Group DEGs", value = TRUE))
                           ),
                           p("The DEGs groups means the DEGs are divided into up and down groups. 
                      DEGs column is the column contains the DOWN/UP/NS groups."),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step16_go_pAdjustMethod", label = "pAdjustMethod", selected = "BH", 
                                                      choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))),
                             column(4, numericInput("in_step16_go_pvalue", label = "pvalueCutoff", value = 1, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step16_go_qvalue", label = "qvalueCutoff", value = 1, min = 0, max = 1, step = 0.01)),
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step16_go_enrich", label = "Enrich GO", icon = icon("calculator"))),
                             column(6, downloadButton("save_step16_go_results", label = "Save go results", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step16_go_enrich {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step16_go_results {background-color: #grey; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Simplify results", status = "primary", collapsible = TRUE, collapsed = TRUE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(6, numericInput("in_step16_go_simplify_cutoff", label = "Cutoff", value = 0.8, min = 0, max = 1, step = 0.01)),
                             column(6, selectizeInput("in_step16_go_simplify_by", label = "Simplify by", selected = "p.adjust", 
                                                      choices = c("p.adjust", "pvalue", "qvalue")))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step16_simplify_go_enrich", label = "Simplify GO", icon = icon("calculator"))),
                             column(6, downloadButton("save_step16_simplify_go_enrich_results", label = "Save simplify results", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step16_simplify_go_enrich {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step16_simplify_go_enrich_results {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                       box(title = "Draw dotplot", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step16_go_ontology", label = "Split Ontology", selected = "ALL",
                                                      choices = c("BP", "CC", "MF", "ALL", "FALSE"), multiple = TRUE)),
                             column(4, selectizeInput("in_step16_go_facet_group1", label = "Facet group x", selected = "", 
                                                      choices = c("ONTOLOGY", "Cluster"), multiple = TRUE)),
                             column(4, selectizeInput("in_step16_go_facet_group2", label = "Facet group y", selected = "",
                                                      choices = c("ONTOLOGY", "Cluster"), multiple = TRUE))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step16_go_pvalue_flt", label = "pvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step16_go_qvalue_flt", label = "qvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step16_go_count_flt", label = "Gene count", value = 2, min = 0, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step16_go_item_num", label = "Item number", value = 15, min = 1, max = 30, step = 1)),
                             column(4, sliderInput("in_step16_go_item_width", label = "Item width", value = 80, min = 0, max = 160, step = 1)),
                             column(4, sliderInput("in_step16_go_grid_width", label = "Grid width", value = 0.4, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step16_go_xmin", label = "X-min", value = NULL, step = 0.01)),
                             column(4, numericInput("in_step16_go_xmax", label = "X-max", value = NULL, step = 0.01)),
                             column(4, numericInput("in_step16_go_break", label = "Breaks", value = 3, min = 2, max = 10, step = 1))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step16_go_dot_size", label = "Dot size", value = c(2, 8), min = 0, max = 15, step = 0.5)),
                             column(6, sliderInput("in_step16_go_font_size", label = "Font size", value = 12, min = 1, max = 20, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step16_go_low_color", label = "Low color", value = 'red')),
                             column(4, textInput("in_step16_go_high_color", label = "High color", value = 'blue')),
                             column(4, numericInput("in_step16_go_dot_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step16_go_title", label = "Title", value = 'GO Enrichment Analysis')),
                             column(4, textInput("in_step16_go_xlab", label = "Xlab", value = 'GeneRatio')),
                             column(4, textInput("in_step16_go_ylab", label = "Ylab", value = 'Go Terms'))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step16_go_out_name", label = "Output plot", value = "enrich", placeholder = "file name prefix")),
                             column(3, numericInput("out_step16_go_width", label = "Figure width", value = 6, min = 1)),
                             column(3, numericInput("out_step16_go_height", label = "Figure height", value = 8, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step16_draw_go_enrich_plot", label = "draw enrich plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step16_go_dot_plot", label = "Save go plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step16_draw_go_enrich_plot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step16_go_dot_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "DEGs",
                                  h4("Show the DEGs table"),
                                  conditionalPanel(condition = ("input.act_step16_import_go_degs > 0"),
                                                   DT::dataTableOutput("out_step16_go_degs") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "GO Terms",
                                  h4("Show the go enrich results"),
                                  conditionalPanel(condition = ("input.act_step16_go_enrich > 0"),
                                                   DT::dataTableOutput("out_step16_go_results") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "GO simplify",
                                  h4("Show the GO simplify results"),
                                  conditionalPanel(condition = ("input.act_step16_simplify_go_enrich > 0"),
                                                   DT::dataTableOutput("out_step16_simplify_go_results") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "GO plot",
                                  h4("Show the go enrich dot plot."),
                                  plotOutput("out_step16_go_dot_plot"))
                       )
                )
              )
      ),
      
      
      ### step 16 compare cluster ###############################
      tabItem(tabName = "step16_compare_cluster",
              h2("Enrich the group GO terms"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the degs table
                       box(title = "Import the merge DEGs table", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step16_group_go_degs", label = "merged DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = FALSE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step16_group_go_degs_sheet", label = "Sheet", value = "DEGs"))
                           ),
                           
                           actionButton("act_step16_import_group_go_degs", label = "import the DEGs", icon = icon("file")),
                           tags$style(HTML("#act_step16_import_group_go_degs {background-color: #75aadb; color: black;}"))
                           
                       ),
                       
                       box(title = "Group enrich GO", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectInput("in_step16_group_go_database", label = "Database", choices = c("gson", "orgdb"), selected = "orgdb")),
                             conditionalPanel(condition = "input.in_step16_group_database == 'gson'",
                                              column(8, fileInput("in_step16_group_go_gson", label = "Gson file", accept = c(".gson", ".GSON"), multiple = FALSE))
                             ),
                             conditionalPanel(condition = "input.in_step16_group_go_database == 'orgdb'",
                                              column(8, selectizeInput("in_step16_group_go_orgdb", label = "Orgdb name", selected = "", choices = ""))
                             )
                           ),
                           p("The orgdb and gson are the name of the database for enrichment.
                      Select one database of them is satisfied."),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step16_group_go_key_type", label = "Key Type", selected = "GID", 
                                                      choices = c("GID", "ENTREZID", "ENSEMBL"))),
                             column(4, textInput("in_step16_group_go_gene_column", label = "Gene column", value = "GeneID", placeholder = "gene column for enrich")),
                             column(4, selectizeInput("in_step16_group_go_group_column", label = "Group column", selected = "",  choices = "", multiple = TRUE))
                           ),
                           p("The Key Type is the IDs from OrgDb for GO enrichment."),
                           p("The Gene column is the IDs from DEGs table corresponding to the Key Type."),
                           
                           fluidRow(
                             column(4, numericInput("in_step16_group_go_pvalue", label = "pvalueCutoff", value = 1, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step16_group_go_qvalue", label = "qvalueCutoff", value = 1, min = 0, max = 1, step = 0.01)),
                             column(4, selectizeInput("in_step16_group_go_pAdjustMethod", label = "pAdjustMethod", selected = "BH",
                                                      choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step16_group_go_enrich", label = "Enrich GO", icon = icon("calculator"))),
                             column(6, downloadButton("save_step16_group_go_enrich_results", label = "Save enrich results", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step16_group_go_enrich {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step16_group_go_enrich_results {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                       box(title = "Simplify results", status = "primary", collapsible = TRUE, collapsed = TRUE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(6, numericInput("in_step16_group_go_simplify_cutoff", label = "Cutoff", value = 0.8, min = 0, max = 1, step = 0.01)),
                             column(6, selectizeInput("in_step16_group_go_simplify_by", label = "Simplify by", selected = "p.adjust", 
                                                      choices = c("p.adjust", "pvalue", "qvalue")))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step16_simplify_group_go_enrich", label = "Simplify GO", icon = icon("calculator"))),
                             column(6, downloadButton("save_step16_simplify_group_go_enrich_results", label = "Save simplify results", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step16_simplify_group_go_enrich {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step16_simplify_group_go_enrich_results {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                       box(title = "Draw dotplot", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, numericInput("in_step16_group_go_pvalue_flt", label = "pvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step16_group_go_qvalue_flt", label = "qvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step16_group_go_count_flt", label = "Gene count", value = 2, min = 0, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step16_group_go_item_num", label = "Item number", value = 5, min = 0, max = 30, step = 1)),
                             column(4, sliderInput("in_step16_group_go_item_width", label = "Item width", value = 80, min = 0, max = 160, step = 1)),
                             column(4, sliderInput("in_step16_group_go_grid_width", label = "Grid width", value = 0.4, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step16_group_go_ontology", label = "Split Ontology", selected = "ALL",
                                                      choices = c("BP", "CC", "MF", "ALL", "FALSE"), multiple = TRUE)),
                             column(4, selectizeInput("in_step16_group_go_facet_group1", label = "Facet group x", selected = "", 
                                                      choices = "", multiple = TRUE)),
                             column(4, selectizeInput("in_step16_group_go_facet_group2", label = "Facet group y", selected = "",
                                                      choices = "", multiple = TRUE))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step16_group_go_dot_size", label = "Dot size", value = c(2, 8), min = 0, max = 15, step = 0.5)),
                             column(6, sliderInput("in_step16_group_go_font_size", label = "Font size", value = 12, min = 1, max = 20, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step16_group_go_low_color", label = "Low color", value = 'red')),
                             column(4, textInput("in_step16_group_go_high_color", label = "High color", value = 'blue')),
                             column(4, sliderInput("in_step16_group_go_dot_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step16_group_go_title", label = "Title", value = 'GO Enrichment Analysis')),
                             column(4, textInput("in_step16_group_go_xlab", label = "Xlab", value = 'Groups')),
                             column(4, textInput("in_step16_group_go_ylab", label = "Ylab", value = 'Go Terms'))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step16_group_go_out_name", label = "Output plot", value = "enrich", placeholder = "file name prefix")),
                             column(3, numericInput("out_step16_group_go_width", label = "Figure width", value = 10, min = 1)),
                             column(3, numericInput("out_step16_group_go_height", label = "Figure height", value = 10, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step16_draw_group_go_enrich_plot", label = "draw enrich plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step16_group_go_enrich_plot", label = "Save enrich plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step16_draw_group_go_enrich_plot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step16_group_go_enrich_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "DEGs",
                                  h4("Show the group DEGs table"),
                                  conditionalPanel(condition = ("input.act_step16_import_group_go_degs > 0"),
                                                   DT::dataTableOutput("out_step16_group_go_degs") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "GO results",
                                  h4("Show the GO enrich results"),
                                  conditionalPanel(condition = ("input.act_step16_group_go_enrich > 0"),
                                                   DT::dataTableOutput("out_step16_group_go_results") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "GO simplify",
                                  h4("Show the GO simplify results"),
                                  conditionalPanel(condition = ("input.act_step16_simplify_group_go_enrich > 0"),
                                                   DT::dataTableOutput("out_step16_simplify_group_go_results") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "GO plot",
                                  h4("Show the GO group enrich plot."),
                                  plotOutput("out_step16_group_go_plot"))
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step17 GO_GSEA ###############################
      
      ### step 17 enrich GO GSEA ###############################
      tabItem(tabName = "step17_go_gsea",
              h2("Enrich the GO GSEA"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the degs table
                       
                       box(title = "Import the DEGs table", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step17_go_gsea_degs", label = "DEGs", accept = c(".xlsx", ".XLSX"), 
                                                 multiple = TRUE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step17_go_gsea_degs_sheet", label = "Sheet", value = "rna_DEGs"))
                           ),
                           p("The DEGs table should contain at least two columns, gene-name and log2FC."),
                           
                           fluidRow(
                             column(6, textInput("in_step17_go_gsea_gene_column", label = "Gene column", value = "GeneID", 
                                                 placeholder = "gene column from enrich")),
                             column(6, textInput("in_step17_go_gsea_log2fc_column", label = "Log2FC column", value = "log2FoldChange", 
                                                 placeholder = "DEGs column from DEGs"))
                           ),
                           p("The Gene column is the IDs from DEGs table corresponding to the Key Type."),
                           
                           fluidRow(
                             column(6, actionButton("act_step17_import_go_gsea_degs", label = "import DEGs table", icon = icon("file"))),
                             column(6, actionButton("act_step17_create_go_gsea_list", label = "create gene list", icon = icon("calculator")))
                           ),
                           
                           tags$style(HTML("#act_step17_import_go_gsea_degs {background-color: #75aadb; color: black;}")),
                           tags$style(HTML("#act_step17_create_go_gsea_list {background-color: #e78c45; color: black;}"))
                       ),
                       
                       box(title = "Enrich GO GSEA", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectInput("in_step17_go_gsea_database", label = "Database", choices = c("gson", "orgdb"), selected = "orgdb")),
                             conditionalPanel(condition = "input.in_step17_go_gsea_database == 'gson'",
                                              column(8, fileInput("in_step17_go_group_gson", label = "Gson file", accept = c(".gson", ".GSON"), multiple = FALSE))
                             ),
                             conditionalPanel(condition = "input.in_step17_go_gsea_database == 'orgdb'",
                                              column(8, selectizeInput("in_step17_go_gsea_orgdb", label = "Orgdb name", selected = "", choices = ""))
                             )
                           ),
                           p("The orgdb and gson are the name of the database for enrichment. 
                             Select one database of them is satisfied."),
                           
                           fluidRow(
                             column(4, numericInput("in_step17_go_gsea_eps", label = "EPS", value = 1e-10)),
                             column(4, numericInput("in_step17_go_gsea_minsize", label = "minGSSize", value = 10)),
                             column(4, numericInput("in_step17_go_gsea_maxsize", label = "maxGSSize", value = 500))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step17_go_gsea_key_type", label = "Key Type", selected = "GID", 
                                                      choices = c("GID", "ENTREZID", "ENSEMBL"))),
                             column(4, selectizeInput("in_step17_go_gsea_pAdjustMethod", label = "pAdjustMethod", selected = "BH", 
                                                      choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))),
                             column(4, numericInput("in_step17_go_gsea_pvalue", label = "pvalueCutoff", value = 1, min = 0, max = 1, step = 0.01))
                           ),
                           p("The Key Type is the IDs from OrgDb for GO enrichment."),
                           
                           fluidRow(
                             column(6, actionButton("act_step17_go_gsea_enrich", label = "Enrich GO GSEA", icon = icon("calculator"))),
                             column(6, downloadButton("save_step17_go_gsea_results", label = "Save GO GSEA results", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step17_go_gsea_enrich {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step17_go_gsea_results {background-color: #grey; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Draw dotplot", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step17_go_gsea_ontology", label = "Split Ontology", selected = "ALL",
                                                      choices = c("BP", "CC", "MF", "ALL", "FALSE"), multiple = TRUE)),
                             column(4, selectizeInput("in_step17_go_gsea_facet_group1", label = "Facet orientation", selected = "horizontal",
                                                      choices = c("horizontal", "vertical", "FALSE"))),
                             column(4, selectizeInput("in_step17_go_gsea_facet_group2", label = "Facet group", selected = "ONTOLOGY",
                                                      choices = c("ONTOLOGY", ".sign")))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step17_go_gsea_pvalue_flt", label = "pvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(6, numericInput("in_step17_go_gsea_qvalue_flt", label = "qvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step17_go_gsea_item_num", label = "Item number", value = 15, min = 1, max = 30, step = 1)),
                             column(4, sliderInput("in_step17_go_gsea_item_width", label = "Item width", value = 80, min = 0, max = 160, step = 1)),
                             column(4, sliderInput("in_step17_go_gsea_grid_width", label = "Grid width", value = 0.4, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step17_go_gsea_xmin", label = "X-min", value = NULL, step = 0.01)),
                             column(4, numericInput("in_step17_go_gsea_xmax", label = "X-max", value = NULL, step = 0.01)),
                             column(4, numericInput("in_step17_go_gsea_break", label = "Breaks", value = 3, min = 2, max = 10, step = 1))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step17_go_gsea_dot_size", label = "Dot size", value = c(2, 8), min = 0, max = 15, step = 0.5)),
                             column(6, sliderInput("in_step17_go_gsea_font_size", label = "Font size", value = 12, min = 1, max = 20, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step17_go_gsea_low_color", label = "Low color", value = 'red')),
                             column(4, textInput("in_step17_go_gsea_high_color", label = "High color", value = 'blue')),
                             column(4, numericInput("in_step17_go_gsea_dot_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step17_go_gsea_title", label = "Title", value = 'GO GSEA')),
                             column(4, textInput("in_step17_go_gsea_xlab", label = "Xlab", value = 'GeneRatio')),
                             column(4, textInput("in_step17_go_gsea_ylab", label = "Ylab", value = 'Go Terms'))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step17_go_gsea_out_name", label = "Output plot", value = "enrich", 
                                                 placeholder = "file name prefix")),
                             column(3, numericInput("out_step17_go_gsea_width", label = "Dotplot width", value = 8, min = 1)),
                             column(3, numericInput("out_step17_go_gsea_height", label = "Dotplot height", value = 8, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step17_draw_go_gsea_dot_plot", label = "draw enrich plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step17_go_gsea_dot_plot", label = "Save GO GSEA plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step17_draw_go_gsea_dot_plot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step17_go_gsea_dot_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                       box(title = "Draw gseaplot", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           textAreaInput("in_step17_go_gsea_line_item_name", label = "Terms IDs", value = "", height = "100px"),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step17_go_gsea_line_color", label = "Terms color", selected = "Paired",
                                                      choices = color_list)),
                             column(4, numericInput("in_step17_go_gsea_line_fontsize", label = "Font size", value = 12, min = 1, max = 20, step = 1)),
                             column(4, selectizeInput("in_step17_go_gsea_line_es_geom", label = "ES geom", selected = "", 
                                                      choices = c("line", "dot")))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step17_go_gsea_line_title", label = "Title", value = 'GO GSEA')),
                             column(4, selectizeInput("in_step17_go_gsea_line_subplot", label = "Subplot", selected = c(1, 2, 3),
                                                      choices = c(1, 2 ,3), multiple = TRUE)),
                             column(4, checkboxInput("in_step17_go_gsea_line_pvalue", label = "Pvalue table", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step17_go_gsea_line_height1", label = "Height 1", value = 1.5, min = 0, max = 5, step = 0.1)),
                             column(4, numericInput("in_step17_go_gsea_line_height2", label = "Height 2", value = 0.5, min = 0, max = 5, step = 0.1)),
                             column(4, numericInput("in_step17_go_gsea_line_height3", label = "Height 3", value = 1, min = 0, max = 5, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step17_go_gsea_line_out_name", label = "Output plot", value = "enrich", 
                                                 placeholder = "file name prefix")),
                             column(3, numericInput("out_step17_go_gsea_line_width", label = "GSEAplot width", value = 9, min = 1)),
                             column(3, numericInput("out_step17_go_gsea_line_height", label = "GSEAplot height", value = 6, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step17_draw_go_gsea_lineplot", label = "draw gsea plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step17_go_gsea_lineplot", label = "Save gsea plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step17_draw_go_gsea_lineplot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step17_go_gsea_lineplot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "DEGs",
                                  h4("Show the DEGs table"),
                                  conditionalPanel(condition = ("input.act_step17_import_go_gsea_degs > 0"),
                                                   DT::dataTableOutput("out_step17_go_gsea_degs") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Gene List",
                                  h4("Show the gene list"),
                                  conditionalPanel(condition = ("input.act_step17_create_go_gsea_list > 0"),
                                                   verbatimTextOutput("out_step17_go_gsea_gene_list") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "GO Terms",
                                  h4("Show the GO enrich results"),
                                  conditionalPanel(condition = ("input.act_step17_go_gsea_enrich > 0"),
                                                   DT::dataTableOutput("out_step17_go_gsea_results") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "GO dotplot",
                                  h4("Show the GO enrich dot plot."),
                                  plotOutput("out_step17_go_gsea_dot_plot")),
                         
                         tabPanel(title = "GO gseaplot",
                                  h4("Show the GO enrich GSEA plot."),
                                  plotOutput("out_step17_go_gsea_line_plot"))
                       )
                )
              )
      ),
      
      
      
      ############################################################################################################
      ## step18 KEGG ###############################
      
      ### step 18 enrich kegg ###############################
      tabItem(tabName = "step18_enrich_kegg",
              h2("Enrich the KEGG terms"),
              fluidRow(
                tags$style(HTML(".download-btn {margin-top: 250px;}",
                                ".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the degs table
                       box(title = "Import the DEGs table", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step18_kegg_degs", label = "DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = TRUE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step18_kegg_degs_sheet", label = "Sheet", value = "rna_DEGs"))
                           ),
                           
                           actionButton("act_step18_import_kegg_degs", label = "import the DEGs", icon = icon("file")),
                           tags$style(HTML("#act_step18_import_kegg_degs {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Enrich KEGG", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step18_kegg_database", label = "Database", selected = "kegg",
                                                      choices = c("gson", "kegg"), multiple = FALSE)),
                             conditionalPanel(condition = "input.in_step18_kegg_database == 'gson'",
                                              column(8, fileInput("in_step18_kegg_gson", label = "Gson file", accept = c(".gson", ".GSON"), multiple = FALSE))
                             ),
                             conditionalPanel(condition = "input.in_step18_kegg_database == 'kegg'",
                                              column(8, textInput("in_step18_kegg_species", label = "Species", value = "", placeholder = "eg. hsa"))
                             ),
                           ),
                           p("The kegg and gson file are the name of the database for enrichment. 
                        Select one database of them is satisfied."),
                           
                           fluidRow(
                             column(6, textInput("in_step18_kegg_gene_column", label = "Gene column", value = "GeneID", placeholder = "gene class from DEGs")),
                             column(6, textInput("in_step18_kegg_degs_column", label = "DEGs column", value = "DEGs", placeholder = "DEGs column from DEGs"))
                           ),
                           p("Gene column is the gene id from DEGs table used for GO enrichment. 
                        DEGs column is the column contains the DOWN/UP/NS class."),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step18_kegg_key_type", label = "Key type", selected = "kegg", 
                                                      choices = c("kegg", "ncbi-geneid", "ncbi-proteinid", "uniprot"))),
                             column(4, checkboxInput("in_step18_kegg_degs_group", label = "Group DEGs", value = TRUE)),
                             column(4, checkboxInput("in_step18_kegg_internal", label = "Internal", value = FALSE))
                           ),
                           p("The (key type) is the gene id from KEGG.db used for the enrichment, sample with (gene column). 
                        The DEGs groups means the DEGs are divided into up and down groups. 
                        The internal means use KEGG.db or latest online KEGG data."),
                           
                           fluidRow(
                             column(4, numericInput("in_step18_kegg_pvalue", label = "pvalueCutoff", value = 1, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step18_kegg_qvalue", label = "qvalueCutoff", value = 1, min = 0, max = 1, step = 0.01)),
                             column(4, selectizeInput("in_step18_kegg_pAdjustMethod", label = "pAdjustMethod", selected = "BH",
                                                      choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step18_kegg_enrich", label = "Enrich KEGG", icon = icon("calculator"))),
                             column(6, downloadButton("save_step18_kegg_results", label = "Save kegg enrich", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step18_kegg_enrich {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step18_kegg_results {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                       box(title = "Enrich dotplot", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, numericInput("in_step18_kegg_pvalue_flt", label = "pvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step18_kegg_qvalue_flt", label = "qvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step18_kegg_count_flt", label = "Gene count", value = 2, min = 0, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step18_kegg_item_num", label = "Item number", value = 15, min = 1, max = 30, step = 1)),
                             column(4, sliderInput("in_step18_kegg_item_width", label = "Item width", value = 80, min = 0, max = 160, step = 1)),
                             column(4, sliderInput("in_step18_kegg_grid_width", label = "Grid width", value = 0.4, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step18_kegg_xmin", label = "X-min", value = NULL, step = 0.01)),
                             column(4, numericInput("in_step18_kegg_xmax", label = "X-max", value = NULL, step = 0.01)),
                             column(4, numericInput("in_step18_kegg_break", label = "Breaks", value = 3, min = 2, max = 10, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step18_kegg_font_size", label = "Font size", value = 12, min = 1, max = 20, step = 1)),
                             column(4, sliderInput("in_step18_kegg_dot_size", label = "Dot size", value = c(2, 8), min = 1, max = 15, step = 0.5)),
                             column(4, selectizeInput("in_step18_kegg_facet_group1", label = "Facet orientation", selected = "horizontal",
                                                      choices = c("horizontal", "vertical", "FALSE")))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step18_kegg_low_color", label = "Low color", value = 'red')),
                             column(4, textInput("in_step18_kegg_high_color", label = "High color", value = 'blue')),
                             column(4, sliderInput("in_step18_kegg_dot_alpha", label = "Alpha", value = 0.9, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step18_kegg_title", label = "Title", value = 'KEGG Enrichment Analysis')),
                             column(4, textInput("in_step18_kegg_xlab", label = "Xlab", value = 'GeneRatio')),
                             column(4, textInput("in_step18_kegg_ylab", label = "Ylab", value = 'KEGG Pathways'))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step18_kegg_out_name", label = "Output plot", value = "enrich", placeholder = "file name prefix")),
                             column(3, numericInput("out_step18_kegg_width", label = "Figure width", value = 6, min = 1)),
                             column(3, numericInput("out_step18_kegg_height", label = "Figure height", value = 6, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step18_draw_enrich_plot", label = "draw enrich plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step18_kegg_plot", label = "Save kegg plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step18_draw_enrich_plot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step18_kegg_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "DEGs",
                                  h4("Show the DEGs table"),
                                  conditionalPanel(condition = ("input.act_step18_import_kegg_degs > 0"),
                                                   DT::dataTableOutput("out_step18_kegg_degs") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "KEGG Pathway",
                                  h4("Show the kegg enrich results"),
                                  conditionalPanel(condition = ("input.act_step18_kegg_enrich > 0"),
                                                   DT::dataTableOutput("out_step18_kegg_results") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "KEGG plot",
                                  h4("Show the kegg enrich dot plot."),
                                  plotOutput("out_step18_kegg_dot_plot"))
                       )
                )
              )
      ),
      
      ### step 18 compare cluster ###############################
      tabItem(tabName = "step18_compare_cluster",
              h2("Enrich the group KEGG terms"),
              fluidRow(
                tags$style(HTML(".download-btn {margin-top: 250px;}",
                                ".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the degs table
                       box(title = "Import the DEGs table", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step18_group_kegg_degs", label = "DEGs", accept = c(".xlsx", ".XLSX"),
                                                 multiple = TRUE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step18_group_kegg_degs_sheet", label = "Sheet", value = "DEGs"))
                           ),
                           
                           actionButton("act_step18_import_group_kegg_degs", label = "import the DEGs", icon = icon("file")),
                           tags$style(HTML("#act_step18_import_group_kegg_degs {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Enrich KEGG", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step18_group_kegg_database", label = "Database", selected = "kegg",
                                                      choices = c("gson", "kegg"), multiple = FALSE)),
                             conditionalPanel(condition = "input.in_step18_group_kegg_database == 'gson'",
                                              column(8, fileInput("in_step18_group_kegg_gson", label = "Gson file", accept = c(".gson", ".GSON"), multiple = FALSE))
                             ),
                             conditionalPanel(condition = "input.in_step18_group_kegg_database == 'kegg'",
                                              column(8, textInput("in_step18_group_kegg_species", label = "Species", value = "", placeholder = "eg. hsa"))
                             ),
                           ),
                           p("The kegg and gson are the name of the database for enrichment.
                      Select one database of them is satisfied."),
                           
                           fluidRow(
                             column(6, textInput("in_step18_group_kegg_gene_column", label = "Gene column", value = "GeneID", placeholder = "gene column from enrich")),
                             column(6, selectizeInput("in_step18_group_kegg_group_column", label = "Group column", selected = "",  choices = "", multiple = TRUE))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step18_group_kegg_key_type", label = "Key Type", selected = "kegg", 
                                                      choices = c("kegg", "ncbi-geneid", "ncbi-proteinid", "uniprot"))),
                             column(6, checkboxInput("in_step18_kegg_internal", label = "Internal", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step18_group_kegg_pvalue", label = "pvalueCutoff", value = 1, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step18_group_kegg_qvalue", label = "qvalueCutoff", value = 1, min = 0, max = 1, step = 0.01)),
                             column(4, selectizeInput("in_step18_group_kegg_pAdjustMethod", label = "pAdjustMethod", selected = "BH",
                                                      choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")))
                           ),
                           
                           fluidRow(
                             column(4, actionButton("act_step18_group_kegg_enrich", label = "Enrich KEGG", icon = icon("calculator"))),
                             column(4, downloadButton("save_step18_group_kegg_enrich_results", label = "Save down results",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step18_group_kegg_enrich {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step18_group_kegg_enrich_results {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                       box(title = "Draw dotplot", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, numericInput("in_step18_group_kegg_pvalue_flt", label = "pvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step18_group_kegg_qvalue_flt", label = "qvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step18_group_kegg_count_flt", label = "Gene count", value = 2, min = 0, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step18_group_kegg_item_num", label = "Item number", value = 5, min = 0, max = 30, step = 1)),
                             column(4, sliderInput("in_step18_group_kegg_item_width", label = "Item width", value = 80, min = 0, max = 160, step = 1)),
                             column(4, sliderInput("in_step18_group_kegg_grid_width", label = "Grid width", value = 0.4, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step18_group_kegg_font_size", label = "Font size", value = 12, min = 0, max = 20, step = 1)),
                             column(6, sliderInput("in_step18_group_kegg_dot_size", label = "Dot size", value = c(2, 8), min = 0, max = 15, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step18_group_kegg_facet_group1", label = "Facet Group x", selected = "", choices = "", multiple = TRUE)),
                             column(6, selectizeInput("in_step18_group_kegg_facet_group2", label = "Facet Group y", selected = "", choices = "", multiple = TRUE))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step18_group_kegg_low_color", label = "Low color", value = 'red')),
                             column(4, textInput("in_step18_group_kegg_high_color", label = "High color", value = 'blue')),
                             column(4, sliderInput("in_step18_group_kegg_dot_alpha", label = "Alpha", value = 0.9, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step18_group_kegg_title", label = "Title", value = 'KEGG Enrichment Analysis')),
                             column(4, textInput("in_step18_group_kegg_xlab", label = "Xlab", value = 'Groups')),
                             column(4, textInput("in_step18_group_kegg_ylab", label = "Ylab", value = 'KEGG Pathways'))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step18_group_kegg_out_name", label = "Output plot", value = "enrich", placeholder = "file name prefix")),
                             column(3, numericInput("out_step18_group_kegg_width", label = "Figure width", value = 6, min = 1)),
                             column(3, numericInput("out_step18_group_kegg_height", label = "Figure height", value = 6, min = 1))
                           ),
                           
                           fluidRow(
                             column(4, actionButton("act_step18_draw_group_kegg_enrich_plot", label = "draw enrich plot", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step18_group_kegg_enrich_plot", label = "Save enrich plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step18_draw_group_kegg_enrich_plot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step18_group_kegg_enrich_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "DEGs",
                                  h4("Show the group DEGs table"),
                                  conditionalPanel(condition = ("input.act_step18_import_group_kegg_degs > 0"),
                                                   DT::dataTableOutput("out_step18_group_kegg_degs") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "KEGG results",
                                  h4("Show the down enrich results"),
                                  conditionalPanel(condition = ("input.act_step18_group_kegg_enrich > 0"),
                                                   DT::dataTableOutput("out_step18_group_kegg_results") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "KEGG plot",
                                  h4("Show the up group enrich plot."),
                                  plotOutput("out_step18_group_kegg_plot"))
                       )
                )
              )
      ),
      
      
      
      ############################################################################################################
      ## step19 KEGG_GSEA ###############################
      
      ### step19 KEGG GSEA ############################
      tabItem(tabName = "step19_kegg_gsea",
              h2("Enrich the KEGG GSEA"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the degs table
                       
                       box(title = "Import the DEGs table", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step19_kegg_gsea_degs", label = "DEGs", accept = c(".xlsx", ".XLSX"), 
                                                 multiple = TRUE, placeholder = "DEGs table")),
                             column(4, textInput(inputId = "in_step19_kegg_gsea_degs_sheet", label = "Sheet", value = "rna_DEGs"))
                           ),
                           p("The DEGs table should contain at least two columns, gene-name and log2FC."),
                           
                           fluidRow(
                             column(6, textInput("in_step19_kegg_gsea_gene_column", label = "Gene column", value = "GeneID", 
                                                 placeholder = "gene column from enrich")),
                             column(6, textInput("in_step19_kegg_gsea_log2fc_column", label = "Log2FC column", value = "log2FoldChange", 
                                                 placeholder = "DEGs column from DEGs"))
                           ),
                           p("The Gene column is the IDs from DEGs table corresponding to the Key Type."),
                           
                           fluidRow(
                             column(6, actionButton("act_step19_import_kegg_gsea_degs", label = "import DEGs table", icon = icon("file"))),
                             column(6, actionButton("act_step19_create_kegg_gsea_list", label = "create gene list", icon = icon("calculator")))
                           ),
                           
                           tags$style(HTML("#act_step19_import_kegg_gsea_degs {background-color: #75aadb; color: black;}")),
                           tags$style(HTML("#act_step19_create_kegg_gsea_list {background-color: #e78c45; color: black;}"))
                       ),
                       
                       box(title = "Enrich KEGG GSEA", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step19_kegg_gsea_database", label = "Database", selected = "kegg",
                                                      choices = c("gson", "kegg"), multiple = FALSE)),
                             conditionalPanel(condition = "input.in_step19_kegg_gsea_database == 'gson'",
                                              column(8, fileInput("in_step19_kegg_gsea_gson", label = "Gson file", accept = c(".gson", ".GSON"), multiple = FALSE))
                             ),
                             conditionalPanel(condition = "input.in_step19_kegg_gsea_database == 'kegg'",
                                              column(4, textInput("in_step19_kegg_gsea_species", label = "Species", value = "", placeholder = "eg. hsa")),
                                              column(4, selectizeInput("in_step19_kegg_gsea_type", label = "KEGG type", selected = "KEGG", 
                                                                       choices = c("KEGG", "MKEGG")))
                             ),
                           ),
                           p("The kegg and gson are the name of the database for enrichment.
                             Select one database of them is satisfied."),
                           
                           fluidRow(
                             column(4, numericInput("in_step19_kegg_gsea_eps", label = "EPS", value = 1e-10)),
                             column(4, numericInput("in_step19_kegg_gsea_minsize", label = "minGSSize", value = 10)),
                             column(4, numericInput("in_step19_kegg_gsea_maxsize", label = "maxGSSize", value = 500))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step19_kegg_gsea_key_type", label = "Key Type", selected = "kegg", 
                                                      choices = c("kegg", "ncbi-geneid", "ncbi-proteinid", "uniprot"))),
                             column(4, selectizeInput("in_step19_kegg_gsea_pAdjustMethod", label = "pAdjustMethod", selected = "BH", 
                                                      choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))),
                             column(4, numericInput("in_step19_kegg_gsea_pvalue", label = "pvalueCutoff", value = 1, min = 0, max = 1, step = 0.01))
                           ),
                           p("The Key Type is the IDs for KEGG enrichment."),
                           
                           fluidRow(
                             column(6, actionButton("act_step19_kegg_gsea_enrich", label = "Enrich KEGG GSEA", icon = icon("calculator"))),
                             column(6, downloadButton("save_step19_kegg_gsea_results", label = "Save KEGG GSEA results", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step19_kegg_gsea_enrich {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step19_kegg_gsea_results {background-color: #grey; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Draw dotplot", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, numericInput("in_step19_kegg_gsea_pvalue_flt", label = "pvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step19_kegg_gsea_qvalue_flt", label = "qvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(4, selectizeInput("in_step19_kegg_gsea_facet_group1", label = "Facet orientation", selected = "horizontal",
                                                      choices = c("horizontal", "vertical", "FALSE")))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step19_kegg_gsea_item_num", label = "Item number", value = 15, min = 1, max = 30, step = 1)),
                             column(4, sliderInput("in_step19_kegg_gsea_item_width", label = "Item width", value = 80, min = 0, max = 160, step = 1)),
                             column(4, sliderInput("in_step19_kegg_gsea_grid_width", label = "Grid width", value = 0.4, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step19_kegg_gsea_xmin", label = "X-min", value = NULL, step = 0.01)),
                             column(4, numericInput("in_step19_kegg_gsea_xmax", label = "X-max", value = NULL, step = 0.01)),
                             column(4, numericInput("in_step19_kegg_gsea_break", label = "Breaks", value = 3, min = 2, max = 10, step = 1))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step19_kegg_gsea_dot_size", label = "Dot size", value = c(2, 8), min = 0, max = 15, step = 0.5)),
                             column(6, sliderInput("in_step19_kegg_gsea_font_size", label = "Font size", value = 12, min = 1, max = 20, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step19_kegg_gsea_low_color", label = "Low color", value = 'red')),
                             column(4, textInput("in_step19_kegg_gsea_high_color", label = "High color", value = 'blue')),
                             column(4, numericInput("in_step19_kegg_gsea_dot_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step19_kegg_gsea_title", label = "Title", value = 'KEGG GSEA')),
                             column(4, textInput("in_step19_kegg_gsea_xlab", label = "Xlab", value = 'GeneRatio')),
                             column(4, textInput("in_step19_kegg_gsea_ylab", label = "Ylab", value = 'KEGG Terms'))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step19_kegg_gsea_out_name", label = "Output plot", value = "enrich", 
                                                 placeholder = "file name prefix")),
                             column(3, numericInput("out_step19_kegg_gsea_width", label = "Dotplot width", value = 8, min = 1)),
                             column(3, numericInput("out_step19_kegg_gsea_height", label = "Dotplot height", value = 8, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step19_draw_kegg_gsea_dot_plot", label = "draw enrich plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step19_kegg_gsea_dot_plot", label = "Save KEGG GSEA plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step19_draw_kegg_gsea_dot_plot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step19_kegg_gsea_dot_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                       box(title = "Draw gseaplot", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           textAreaInput("in_step19_kegg_gsea_line_item_name", label = "Terms IDs", value = "", height = "100px"),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step19_kegg_gsea_line_color", label = "Terms color", selected = "Paired",
                                                      choices = color_list)),
                             column(4, numericInput("in_step19_kegg_gsea_line_fontsize", label = "Font size", value = 12, min = 1, max = 20, step = 1)),
                             column(4, selectizeInput("in_step19_kegg_gsea_line_es_geom", label = "ES geom", selected = "", 
                                                      choices = c("line", "dot")))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step19_kegg_gsea_line_title", label = "Title", value = 'KEGG GSEA')),
                             column(4, selectizeInput("in_step19_kegg_gsea_line_subplot", label = "Subplot", selected = c(1, 2, 3),
                                                      choices = c(1, 2 ,3), multiple = TRUE)),
                             column(4, checkboxInput("in_step19_kegg_gsea_line_pvalue", label = "Pvalue table", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step19_kegg_gsea_line_height1", label = "Height 1", value = 1.5, min = 0, max = 5, step = 0.1)),
                             column(4, numericInput("in_step19_kegg_gsea_line_height2", label = "Height 2", value = 0.5, min = 0, max = 5, step = 0.1)),
                             column(4, numericInput("in_step19_kegg_gsea_line_height3", label = "Height 3", value = 1, min = 0, max = 5, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step19_kegg_gsea_line_out_name", label = "Output plot", value = "enrich", 
                                                 placeholder = "file name prefix")),
                             column(3, numericInput("out_step19_kegg_gsea_line_width", label = "GSEAplot width", value = 9, min = 1)),
                             column(3, numericInput("out_step19_kegg_gsea_line_height", label = "GSEAplot height", value = 6, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step19_draw_kegg_gsea_lineplot", label = "draw gsea plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step19_kegg_gsea_lineplot", label = "Save gsea plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step19_draw_kegg_gsea_lineplot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step19_kegg_gsea_lineplot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "DEGs",
                                  h4("Show the DEGs table"),
                                  conditionalPanel(condition = ("input.act_step19_import_kegg_gsea_degs > 0"),
                                                   DT::dataTableOutput("out_step19_kegg_gsea_degs") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Gene List",
                                  h4("Show the gene list"),
                                  conditionalPanel(condition = ("input.act_step19_create_kegg_gsea_list > 0"),
                                                   verbatimTextOutput("out_step19_kegg_gsea_gene_list") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "KEGG Terms",
                                  h4("Show the KEGG enrich results"),
                                  conditionalPanel(condition = ("input.act_step19_kegg_gsea_enrich > 0"),
                                                   DT::dataTableOutput("out_step19_kegg_gsea_results") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "KEGG dotplot",
                                  h4("Show the KEGG enrich dot plot."),
                                  plotOutput("out_step19_kegg_gsea_dot_plot")),
                         
                         tabPanel(title = "KEGG gseaplot",
                                  h4("Show the KEGG enrich GSEA plot."),
                                  plotOutput("out_step19_kegg_gsea_line_plot"))
                       )
                )
              )
      ),
      
      
      
      ############################################################################################################
      ## step20 Codon Enrich ###############################
      
      ### step20.1 Codon usage ###############################
      tabItem(tabName = "step20_codon_usage",
              h2("Codon usage."),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the codon table
                       
                       box(title = "Import the codon usage table", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step20_codon_table", label = "Codon count", accept = c(".xlsx", ".XLSX"), 
                                                 multiple = F, placeholder = "codon table")),
                             column(4, selectizeInput("in_step20_codon_table_sheet", label = "Sheet", selected = "frequency",
                                                      choices = c("Frequency", "RSCU", "CAI"), multiple = FALSE))
                           ),
                           
                           actionButton("act_step20_codon_import_table", label = "import codon table", icon = icon("file")),
                           
                           tags$style(HTML("#act_step20_codon_import_table {background-color: #75aadb; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Filter the gene table", status = "primary", collapsible = TRUE, collapsed = TRUE, 
                           solidHeader = TRUE, width = NULL,
                           p("The file needs to contain a list of gene names and a list of classification message."),
                           fluidRow(
                             column(8, fileInput(inputId = "in_step20_codon_gene_table", label = "Gene table", accept = c(".xlsx", ".XLSX"), 
                                                 multiple = F, placeholder = "gene table with classification")),
                             column(4, numericInput("in_step20_codon_gene_table_sheet", label = "Sheet number", value = 1))
                           ),

                           actionButton("act_step20_codon_import_gene_table", label = "import gene table", icon = icon("file")),
                           
                           tags$style(HTML("#act_step20_codon_import_gene_table {background-color: #75aadb; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Filter codon usage table", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           selectizeInput("in_step20_codon_list", label = "Codon", selected = "AAA",
                                          choices = default_codon_list, multiple = TRUE),
                           
                           fluidRow(
                             column(6, textInput("in_step20_codon_table_column", label = "Codon table column", value = "Gene",
                                                 placeholder = "Gene")),
                             column(6, textInput("in_step20_codon_gene_table_column", label = "Gene table column", value = "GeneID",
                                                 placeholder = "Gene"))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step20_codon_gene_table_class", label = "Gene class column", value = "Sections",
                                                 placeholder = "Class")),
                             column(6, numericInput("in_step20_codon_value", label = "Min codon value", value = 0, step = 0.1))
                           ),
                           
                           textInput("out_step20_codon_table_name", label = "Output table name", value = "filtered-codon-usage", 
                                     placeholder = "eg.: filtered-codon-usage"),
                           
                           fluidRow(
                             column(6, actionButton("act_step20_codon_filter_table", label = "Filter codon table", icon = icon("list"))),
                             column(6, downloadButton("act_step20_codon_save_table", label = "Save codon table", class = "download-btn",
                                                      icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step20_codon_filter_table {background-color: #e78c45; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#act_step20_codon_save_table {background-color: #e78c45; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Draw codon usage boxplot", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,

                           fluidRow(
                             column(4, numericInput("in_step20_codon_box_y_min", label = "Ylim min", value = 0)),
                             column(4, numericInput("in_step20_codon_box_y_max", label = "Ylim max", value = NA)),
                             column(4, numericInput("in_step20_codon_box_y_breaks", label = "Ylim breaks", value = 5))
                           ),
                          
                          fluidRow(
                             column(4, selectizeInput("in_step20_codon_box_y_sqrt", label = "Y sqrt", selected = "none",
                                                      choices = c("none", "log10", "sqrt"))),
                             column(4, selectizeInput("in_step20_codon_box_legend", label = "Legend", selected = "none",
                                                      choices = c("none", "top", "bottom", "left", "right"))),
                             column(4, numericInput("in_step20_codon_box_x_angle", label = "Xticks angle", value = 90))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step20_codon_box_color", label = "Box color", selected = "Paired",
                                                      choices = color_list)),
                             column(4, numericInput("in_step20_codon_box_alpha", label = "Box alpha", value = 0.8, min = 0, max = 1, step = 0.05)),
                             column(4, numericInput("in_step20_codon_box_width", label = "Box width", value = 0.8, min = 0.1, max = 2, step = 0.05))
                           ),

                           fluidRow(
                             column(4, numericInput("in_step20_codon_box_line", label = "Line width", value = 0.2, min = 0, max = 1, step = 0.05)),
                             column(4, textInput("in_step20_codon_box_line_color", label = "Line color", value = "black")),
                             column(4, numericInput("in_step20_codon_box_fontsize", label = "Font size", value = 15, min = 1, max = 20, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, checkboxInput("in_step20_codon_box_outlier", label = "Outlier", value = TRUE)),
                             column(4, textInput("in_step20_codon_box_outlier_color", label = "Outlier color", value = "grey30")),
                             column(4, numericInput("in_step20_codon_box_outlier_size", label = "Outlier size", value = 1, min = 0, max = 5, step = 0.1))
                            ),
                           
                           fluidRow(
                             column(4, textInput("out_step20_codon_box_title", label = "Title", value = "Codon usage")),
                             column(4, textInput("out_step20_codon_box_xlabel", label = "Xlabel", value = "Class")),
                             column(4, textInput("out_step20_codon_box_ylabel", label = "Ylabel", value = "Frequency"))
                           ),
                           
                          fluidRow(
                            column(4, textInput("out_step20_codon_box_name", label = "Output", value = "codon-usage")),
                            column(4, numericInput("out_step20_codon_box_width", label = "Width", value = 5)),
                            column(4, numericInput("out_step20_codon_box_height", label = "Height", value = 6))
                          ),
                          
                           fluidRow(
                             column(6, actionButton("act_step20_codon_draw_box", label = "Draw box plot", icon = icon("chart-bar"))),
                             column(6, downloadButton("save_step20_codon_save_box", label = "Save box plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step20_codon_draw_box {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step20_codon_save_box {background-color: #grey; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Draw cudon usage cumulative", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, numericInput("in_step20_codon_cdf_x_min", label = "Xlim min", value = NA)),
                             column(4, numericInput("in_step20_codon_cdf_x_max", label = "Xlim max", value = NA)),
                             column(4, numericInput("in_step20_codon_cdf_x_breaks", label = "Xlim breaks", value = 5))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step20_codon_cdf_color", label = "Line color", selected = "Paired",
                                                      choices = color_list)),
                             column(4, numericInput("in_step20_codon_cdf_alpha", label = "Line alpha", value = 0.9, min = 0, max = 1, step = 0.05)),
                             column(4, numericInput("in_step20_codon_cdf_width", label = "Line width", value = 0.8, min = 0, max = 5, step = 0.05)),
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step20_codon_cdf_fontsize", label = "Font size", value = 15, min = 1, max = 20, step = 1)),
                             column(4, selectizeInput("in_step20_codon_cdf_legend", label = "Legend", selected = "none",
                                                      choices = c("none", "top", "bottom", "left", "right"))),
                             column(4, checkboxInput("in_step20_codon_cdf_facet", label = "Facet", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step20_codon_cdf_title", label = "Title", value = "Codon usage")),
                             column(4, textInput("in_step20_codon_cdf_xlabel", label = "Xlabel", value = "Codon")),
                             column(4, textInput("in_step20_codon_cdf_ylabel", label = "Ylabel", value = "eCDF (Frequency)"))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step20_codon_cdf_name", label = "Output", value = "codon-usage")),
                             column(4, numericInput("out_step20_codon_cdf_width", label = "Width", value = 5)),
                             column(4, numericInput("out_step20_codon_cdf_height", label = "Height", value = 6))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step20_codon_draw_cdf", label = "Draw cdf plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step20_codon_save_cdf", label = "Save cdf plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step20_codon_draw_cdf {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step20_codon_save_cdf {background-color: #grey; color: black; margin-top: 0px;}"))
                       ),
                       
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Codon table",
                                  h4("Show the codon table"),
                                  conditionalPanel(condition = ("input.act_step20_codon_import_table > 0"),
                                                   DT::dataTableOutput("out_step20_codon_table") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Gene table",
                                  h4("Show the gene table"),
                                  conditionalPanel(condition = ("input.act_step20_codon_import_gene_table > 0"),
                                                   DT::dataTableOutput("out_step20_codon_gene_table") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Filtered codon table",
                                  h4("Show the filtered codon table"),
                                  conditionalPanel(condition = ("input.act_step20_codon_filter_table > 0"),
                                                   DT::dataTableOutput("out_step20_codon_filter_table") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Box plot",
                                  h4("Show the box plot"),
                                  conditionalPanel(condition = ("input.act_step20_codon_draw_box > 0"),
                                                   plotOutput("out_step20_codon_save_box") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "CDF plot",
                                  h4("Show the CDF plot"),
                                  conditionalPanel(condition = ("input.act_step20_codon_draw_cdf > 0"),
                                                   plotOutput("out_step20_codon_save_cdf") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      
      ### step20.2 codon GO GSEA enrichment ###############################
      tabItem(tabName = "step20_go_gsea_enrich",
              h2("GO GSEA enrichment of the gene codon."),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the codon table
                       
                       box(title = "Import the codon count table", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step20_go_codon_table", label = "Codon count", accept = c(".xlsx", ".XLSX"), 
                                                 multiple = F, placeholder = "codon table")),
                             column(4, selectizeInput("in_step20_go_codon_table_sheet", label = "Sheet", selected = "frequency",
                                                      choices = c("Frequency", "RSCU", "CAI"), multiple = FALSE))
                           ),
                           
                           actionButton("act_step20_go_import_codon_table", label = "import codon table", icon = icon("file")),
                           
                           tags$style(HTML("#act_step20_go_import_codon_table {background-color: #75aadb; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Create codon table vector", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step20_go_codon_table_codon", label = "Codon", selected = "AAA",
                                                      choices = default_codon_list, multiple = TRUE)),
                             column(4, numericInput("in_step20_go_codon_table_number", label = "Number", value = 1)),
                             column(4, textInput("in_step20_go_codon_gene_column", label = "Gene column", value = 'Gene'))
                           ),
                           
                           actionButton("act_step20_go_create_codon_table_vector", label = "Create codon vector", icon = icon("list")),
                           
                           tags$style(HTML("#act_step20_go_create_codon_table_vector {background-color: #e78c45; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Enrich GO GSEA", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectInput("in_step20_go_codon_gsea_database", label = "Database", choices = c("gson", "orgdb"), selected = "orgdb")),
                             conditionalPanel(condition = "input.in_step20_go_codon_gsea_database == 'gson'",
                                              column(8, fileInput("in_step20_go_codon_group_gson", label = "Gson file", accept = c(".gson", ".GSON"), multiple = FALSE))
                             ),
                             conditionalPanel(condition = "input.in_step20_go_codon_gsea_database == 'orgdb'",
                                              column(8, selectizeInput("in_step20_go_codon_gsea_orgdb", label = "Orgdb name", selected = "", choices = ""))
                             )
                           ),
                           p("The orgdb and gson are the name of the database for enrichment. 
                             Select one database of them is satisfied."),
                           
                           fluidRow(
                             column(4, numericInput("in_step20_go_codon_gsea_eps", label = "EPS", value = 1e-10)),
                             column(4, numericInput("in_step20_go_codon_gsea_minsize", label = "minGSSize", value = 10)),
                             column(4, numericInput("in_step20_go_codon_gsea_maxsize", label = "maxGSSize", value = 500))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step20_go_codon_gsea_key_type", label = "Key Type", selected = "GID", 
                                                      choices = c("GID", "ENTREZID", "ENSEMBL"))),
                             column(4, selectizeInput("in_step20_go_codon_gsea_pAdjustMethod", label = "pAdjustMethod", selected = "BH", 
                                                      choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))),
                             column(4, numericInput("in_step20_go_codon_gsea_pvalue", label = "pvalueCutoff", value = 1, min = 0, max = 1, step = 0.01))
                           ),
                           p("The Key Type is the IDs from OrgDb for GO enrichment."),
                           
                           fluidRow(
                             column(6, actionButton("act_step20_go_codon_gsea_enrich", label = "Enrich GO GSEA", icon = icon("calculator"))),
                             column(6, downloadButton("save_step20_go_codon_gsea_results", label = "Save GO GSEA results", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step20_go_codon_gsea_enrich {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step20_go_codon_gsea_results {background-color: #grey; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Draw dotplot", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step20_go_codon_gsea_ontology", label = "Split Ontology", selected = "ALL",
                                                      choices = c("BP", "CC", "MF", "ALL", "FALSE"), multiple = TRUE)),
                             column(4, selectizeInput("in_step20_go_codon_gsea_facet_group1", label = "Facet orientation", selected = "horizontal",
                                                      choices = c("horizontal", "vertical", "FALSE"))),
                             column(4, selectizeInput("in_step20_go_codon_gsea_facet_group2", label = "Facet group", selected = "ONTOLOGY",
                                                      choices = c("ONTOLOGY", ".sign")))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step20_go_codon_gsea_pvalue_flt", label = "pvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(6, numericInput("in_step20_go_codon_gsea_qvalue_flt", label = "qvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step20_go_codon_gsea_item_num", label = "Item number", value = 15, min = 1, max = 30, step = 1)),
                             column(4, sliderInput("in_step20_go_codon_gsea_item_width", label = "Item width", value = 80, min = 0, max = 160, step = 1)),
                             column(4, sliderInput("in_step20_go_codon_gsea_grid_width", label = "Grid width", value = 0.4, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step20_go_codon_gsea_xmin", label = "X-min", value = NULL, step = 0.01)),
                             column(4, numericInput("in_step20_go_codon_gsea_xmax", label = "X-max", value = NULL, step = 0.01)),
                             column(4, numericInput("in_step20_go_codon_gsea_break", label = "Breaks", value = 3, min = 2, max = 10, step = 1))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step20_go_codon_gsea_dot_size", label = "Dot size", value = c(2, 8), min = 0, max = 15, step = 0.5)),
                             column(6, sliderInput("in_step20_go_codon_gsea_font_size", label = "Font size", value = 12, min = 1, max = 20, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step20_go_codon_gsea_low_color", label = "Low color", value = 'red')),
                             column(4, textInput("in_step20_go_codon_gsea_high_color", label = "High color", value = 'blue')),
                             column(4, numericInput("in_step20_go_codon_gsea_dot_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step20_go_codon_gsea_title", label = "Title", value = 'GO GSEA')),
                             column(4, textInput("in_step20_go_codon_gsea_xlab", label = "Xlab", value = 'GeneRatio')),
                             column(4, textInput("in_step20_go_codon_gsea_ylab", label = "Ylab", value = 'Go Terms'))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step20_go_codon_gsea_out_name", label = "Output plot", value = "enrich", 
                                                 placeholder = "file name prefix")),
                             column(3, numericInput("out_step20_go_codon_gsea_width", label = "Dotplot width", value = 8, min = 1)),
                             column(3, numericInput("out_step20_go_codon_gsea_height", label = "Dotplot height", value = 8, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step20_draw_go_codon_gsea_dot_plot", label = "draw enrich plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step20_go_codon_gsea_dot_plot", label = "Save GO plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step20_draw_go_codon_gsea_dot_plot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step20_go_codon_gsea_dot_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                       box(title = "Draw gseaplot", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           textAreaInput("in_step20_go_codon_gsea_line_item_name", label = "Terms IDs", value = "", height = "100px"),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step20_go_codon_gsea_line_color", label = "Terms color", selected = "Paired",
                                                      choices = color_list)),
                             column(4, numericInput("in_step20_go_codon_gsea_line_fontsize", label = "Font size", value = 12, min = 1, max = 20, step = 1)),
                             column(4, selectizeInput("in_step20_go_codon_gsea_line_es_geom", label = "ES geom", selected = "", 
                                                      choices = c("line", "dot")))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step20_go_codon_gsea_line_title", label = "Title", value = 'GO GSEA')),
                             column(4, selectizeInput("in_step20_go_codon_gsea_line_subplot", label = "Subplot", selected = c(1, 2, 3),
                                                      choices = c(1, 2 ,3), multiple = TRUE)),
                             column(4, checkboxInput("in_step20_go_codon_gsea_line_pvalue", label = "Pvalue table", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step20_go_codon_gsea_line_height1", label = "Height 1", value = 1.5, min = 0, max = 5, step = 0.1)),
                             column(4, numericInput("in_step20_go_codon_gsea_line_height2", label = "Height 2", value = 0.5, min = 0, max = 5, step = 0.1)),
                             column(4, numericInput("in_step20_go_codon_gsea_line_height3", label = "Height 3", value = 1, min = 0, max = 5, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step20_go_codon_gsea_line_out_name", label = "Output plot", value = "enrich", 
                                                 placeholder = "file name prefix")),
                             column(3, numericInput("out_step20_go_codon_gsea_line_width", label = "GSEAplot width", value = 9, min = 1)),
                             column(3, numericInput("out_step20_go_codon_gsea_line_height", label = "GSEAplot height", value = 6, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step20_draw_go_codon_gsea_lineplot", label = "draw gsea plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step20_go_codon_gsea_lineplot", label = "Save gsea plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step20_draw_go_codon_gsea_lineplot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step20_go_codon_gsea_lineplot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                           
                       )
                       
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Codon table",
                                  h4("Show the codon table"),
                                  conditionalPanel(condition = ("input.act_step20_go_import_codon_table > 0"),
                                                   DT::dataTableOutput("out_step20_go_codon_table") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Codon vector",
                                  h4("Show the codon table vector"),
                                  conditionalPanel(condition = ("input.act_step20_go_create_codon_table_vector > 0"),
                                                   verbatimTextOutput("out_step20_go_codon_table_vector") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "GO Terms",
                                  h4("Show the GO enrich results"),
                                  conditionalPanel(condition = ("input.act_step20_go_codon_gsea_enrich > 0"),
                                                   DT::dataTableOutput("out_step20_go_codon_gsea_results") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "GO dotplot",
                                  h4("Show the GO enrich dot plot."),
                                  plotOutput("out_step20_go_codon_gsea_dot_plot")),
                         
                         tabPanel(title = "GO gseaplot",
                                  h4("Show the GO enrich GSEA plot."),
                                  plotOutput("out_step20_go_codon_gsea_line_plot"))
                       )
                )
              )
      ),
      
      ### step20.3 codon KEGG GSEA enrichment ###############################
      tabItem(tabName = "step20_kegg_gsea_enrich",
              h2("Enrich the gene codon."),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the codon table
                       
                       box(title = "Import the codon count table", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(8, fileInput(inputId = "in_step20_kegg_codon_table", label = "Codon count", accept = c(".xlsx", ".XLSX"), 
                                                 multiple = F, placeholder = "codon table")),
                             column(4, selectizeInput("in_step20_kegg_codon_table_sheet", label = "Sheet", selected = "frequency",
                                                      choices = c("Frequency", "RSCU", "CAI"), multiple = FALSE))
                           ),
                           
                           actionButton("act_step20_kegg_import_codon_table", label = "import codon table", icon = icon("file")),
                           
                           tags$style(HTML("#act_step20_kegg_import_codon_table {background-color: #75aadb; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Create codon table vector", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step20_kegg_codon_table_codon", label = "Codon", selected = "AAA",
                                                      choices = default_codon_list, multiple = TRUE)),
                             column(4, numericInput("in_step20_kegg_codon_table_number", label = "Number", value = 1)),
                             column(4, textInput("in_step20_kegg_codon_gene_column", label = "Gene column", value = 'Gene'))
                           ),
                           
                           actionButton("act_step20_kegg_create_codon_table_vector", label = "Create codon vector", icon = icon("list")),
                           
                           tags$style(HTML("#act_step20_kegg_create_codon_table_vector {background-color: #e78c45; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Enrich codon KEGG GSEA", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step20_kegg_codon_gsea_database", label = "Database", selected = "kegg",
                                                      choices = c("gson", "kegg"), multiple = FALSE)),
                             conditionalPanel(condition = "input.in_step20_kegg_codon_gsea_database == 'gson'",
                                              column(8, fileInput("in_step20_kegg_codon_gsea_gson", label = "Gson file", accept = c(".gson", ".GSON"), multiple = FALSE))
                             ),
                             conditionalPanel(condition = "input.in_step20_kegg_codon_gsea_database == 'kegg'",
                                              column(4, textInput("in_step20_kegg_codon_gsea_species", label = "Species", value = "", placeholder = "eg. hsa")),
                                              column(4, selectizeInput("in_step20_kegg_codon_gsea_type", label = "KEGG type", selected = "KEGG", 
                                                                       choices = c("KEGG", "MKEGG")))
                             ),
                           ),
                           p("The kegg and gson are the name of the database for enrichment.
                             Select one database of them is satisfied."),
                           
                           fluidRow(
                             column(4, numericInput("in_step20_kegg_codon_gsea_eps", label = "EPS", value = 1e-10)),
                             column(4, numericInput("in_step20_kegg_codon_gsea_minsize", label = "minGSSize", value = 10)),
                             column(4, numericInput("in_step20_kegg_codon_gsea_maxsize", label = "maxGSSize", value = 500))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step20_kegg_codon_gsea_key_type", label = "Key Type", selected = "kegg", 
                                                      choices = c("kegg", "ncbi-geneid", "ncbi-proteinid", "uniprot"))),
                             column(4, selectizeInput("in_step20_kegg_codon_gsea_pAdjustMethod", label = "pAdjustMethod", selected = "BH", 
                                                      choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))),
                             column(4, numericInput("in_step20_kegg_codon_gsea_pvalue", label = "pvalueCutoff", value = 1, min = 0, max = 1, step = 0.01))
                           ),
                           p("The Key Type is the IDs for KEGG enrichment."),
                           
                           fluidRow(
                             column(6, actionButton("act_step20_kegg_codon_gsea_enrich", label = "Enrich KEGG GSEA", icon = icon("calculator"))),
                             column(6, downloadButton("save_step20_kegg_codon_gsea_results", label = "Save KEGG GSEA results", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step20_kegg_codon_gsea_enrich {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step20_kegg_codon_gsea_results {background-color: #grey; color: black; margin-top: 0px;}"))
                       ),
                       
                       box(title = "Draw dotplot", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, numericInput("in_step20_kegg_codon_gsea_pvalue_flt", label = "pvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(4, numericInput("in_step20_kegg_codon_gsea_qvalue_flt", label = "qvalueCutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                             column(4, selectizeInput("in_step20_kegg_codon_gsea_facet_group1", label = "Facet orientation", selected = "horizontal",
                                                      choices = c("horizontal", "vertical", "FALSE")))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step20_kegg_codon_gsea_item_num", label = "Item number", value = 15, min = 1, max = 30, step = 1)),
                             column(4, sliderInput("in_step20_kegg_codon_gsea_item_width", label = "Item width", value = 80, min = 0, max = 160, step = 1)),
                             column(4, sliderInput("in_step20_kegg_codon_gsea_grid_width", label = "Grid width", value = 0.4, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step20_kegg_codon_gsea_xmin", label = "X-min", value = NULL, step = 0.01)),
                             column(4, numericInput("in_step20_kegg_codon_gsea_xmax", label = "X-max", value = NULL, step = 0.01)),
                             column(4, numericInput("in_step20_kegg_codon_gsea_break", label = "Breaks", value = 3, min = 2, max = 10, step = 1))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step20_kegg_codon_gsea_dot_size", label = "Dot size", value = c(2, 8), min = 0, max = 15, step = 0.5)),
                             column(6, sliderInput("in_step20_kegg_codon_gsea_font_size", label = "Font size", value = 12, min = 1, max = 20, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step20_kegg_codon_gsea_low_color", label = "Low color", value = 'red')),
                             column(4, textInput("in_step20_kegg_codon_gsea_high_color", label = "High color", value = 'blue')),
                             column(4, numericInput("in_step20_kegg_codon_gsea_dot_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step20_kegg_codon_gsea_title", label = "Title", value = 'KEGG GSEA')),
                             column(4, textInput("in_step20_kegg_codon_gsea_xlab", label = "Xlab", value = 'GeneRatio')),
                             column(4, textInput("in_step20_kegg_codon_gsea_ylab", label = "Ylab", value = 'KEGG Terms'))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step20_kegg_codon_gsea_out_name", label = "Output plot", value = "enrich", 
                                                 placeholder = "file name prefix")),
                             column(3, numericInput("out_step20_kegg_codon_gsea_width", label = "Dotplot width", value = 8, min = 1)),
                             column(3, numericInput("out_step20_kegg_codon_gsea_height", label = "Dotplot height", value = 8, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step20_draw_kegg_codon_gsea_dot_plot", label = "draw enrich plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step20_kegg_codon_gsea_dot_plot", label = "Save KEGG plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step20_draw_kegg_codon_gsea_dot_plot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step20_kegg_codon_gsea_dot_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                       box(title = "Draw gseaplot", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           textAreaInput("in_step20_kegg_codon_gsea_line_item_name", label = "Terms IDs", value = "", height = "100px"),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step20_kegg_codon_gsea_line_color", label = "Terms color", selected = "Paired",
                                                      choices = color_list)),
                             column(4, numericInput("in_step20_kegg_codon_gsea_line_fontsize", label = "Font size", value = 12, min = 1, max = 20, step = 1)),
                             column(4, selectizeInput("in_step20_kegg_codon_gsea_line_es_geom", label = "ES geom", selected = "", 
                                                      choices = c("line", "dot")))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step20_kegg_codon_gsea_line_title", label = "Title", value = 'KEGG GSEA')),
                             column(4, selectizeInput("in_step20_kegg_codon_gsea_line_subplot", label = "Subplot", selected = c(1, 2, 3),
                                                      choices = c(1, 2 ,3), multiple = TRUE)),
                             column(4, checkboxInput("in_step20_kegg_codon_gsea_line_pvalue", label = "Pvalue table", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step20_kegg_codon_gsea_line_height1", label = "Height 1", value = 1.5, min = 0, max = 5, step = 0.1)),
                             column(4, numericInput("in_step20_kegg_codon_gsea_line_height2", label = "Height 2", value = 0.5, min = 0, max = 5, step = 0.1)),
                             column(4, numericInput("in_step20_kegg_codon_gsea_line_height3", label = "Height 3", value = 1, min = 0, max = 5, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("out_step20_kegg_codon_gsea_line_out_name", label = "Output plot", value = "enrich", 
                                                 placeholder = "file name prefix")),
                             column(3, numericInput("out_step20_kegg_codon_gsea_line_width", label = "GSEAplot width", value = 9, min = 1)),
                             column(3, numericInput("out_step20_kegg_codon_gsea_line_height", label = "GSEAplot height", value = 6, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step20_draw_kegg_codon_gsea_lineplot", label = "draw gsea plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step20_kegg_codon_gsea_lineplot", label = "Save gsea plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step20_draw_kegg_codon_gsea_lineplot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step20_kegg_codon_gsea_lineplot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Codon table",
                                  h4("Show the codon table"),
                                  conditionalPanel(condition = ("input.act_step20_kegg_import_codon_table > 0"),
                                                   DT::dataTableOutput("out_step20_kegg_codon_table") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Codon vector",
                                  h4("Show the codon table vector"),
                                  conditionalPanel(condition = ("input.act_step20_kegg_create_codon_table_vector > 0"),
                                                   verbatimTextOutput("out_step20_kegg_codon_table_vector") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "KEGG Terms",
                                  h4("Show the KEGG enrich results"),
                                  conditionalPanel(condition = ("input.act_step20_kegg_codon_gsea_enrich > 0"),
                                                   DT::dataTableOutput("out_step20_kegg_codon_gsea_results") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "KEGG dotplot",
                                  h4("Show the KEGG enrich dot plot."),
                                  plotOutput("out_step20_kegg_codon_gsea_dot_plot")),
                         
                         tabPanel(title = "KEGG gseaplot",
                                  h4("Show the KEGG enrich GSEA plot."),
                                  plotOutput("out_step20_kegg_codon_gsea_line_plot"))
                       )
                )
              )
      ),
      
      
      
      ############################################################################################################
      ## step21 Gene plot ###############################
      ### step21 draw the gene plot  ###############################
      tabItem(tabName = "step21_isoforms",
              h2("Draw the isoforms density."),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the reads density table
                       box(title = "Import the reads density", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step21_isoforms_rna_seq", label = "RNA-seq", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "RNA-seq density file"),
                           
                           fileInput(inputId = "in_step21_isoforms_ribo_seq", label = "Ribo-seq", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "Ribo-seq density file"),
                           
                           fluidRow(
                             column(6, actionButton("act_step21_isoforms_import_rna_seq", label = "import the RNA-seq", icon = icon("file"))),
                             column(6, actionButton("act_step21_isoforms_import_ribo_seq", label = "import the Ribo-seq", icon = icon("file")))
                           ),
                           
                           tags$style(HTML("#act_step21_isoforms_import_rna_seq {background-color: #75aadb; color: black;}")),
                           tags$style(HTML("#act_step21_isoforms_import_ribo_seq {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Set the groups", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, width = 14,
                           fluidRow(
                             column(4, selectizeInput("in_step21_isoforms_rna_column", label = "RNA column", choices = NULL, selected = NULL, multiple = FALSE)),
                             column(8, selectizeInput("in_step21_isoforms_rna_groups", label = "RNA groups", choices = NULL, selected = NULL, multiple = TRUE))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step21_isoforms_ribo_column", label = "Ribo column", choices = NULL, selected = NULL, multiple = FALSE)),
                             column(8, selectizeInput("in_step21_isoforms_ribo_groups", label = "Ribo groups", choices = NULL, selected = NULL, multiple = TRUE))
                           )
                       ),
                       
                       box(title = "Draw reads density", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step21_isoforms_gene", label = "Gene name", selected = "", 
                                                      choices = NULL, multiple = FALSE)),
                             column(4, selectizeInput("in_step21_isoforms_type", label = "Figure type", selected = "line", 
                                                      choices = c("line", "bar", "area"))),
                             column(4, selectizeInput("in_step21_isoforms_legend", label = "Legend", selected = "none", 
                                                      choices = c("none", "right", "bottom", "top")))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step21_isoforms_rna_x_min", label = "RNA X min", value = 1, min = 1, step = 1)),
                             column(4, numericInput("in_step21_isoforms_rna_x_max", label = "RNA X max", value = NA, min = 1, step = 1)),
                             column(4, numericInput("in_step21_isoforms_rna_x_break", label = "RNA X breaks", value = 5, min = 3, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step21_isoforms_rna_y_min", label = "RNA Y min", value = 0, min = 0, step = 1)),
                             column(4, numericInput("in_step21_isoforms_rna_y_max", label = "RNA Y max", value = NA, min = 0, step = 1)),
                             column(4, numericInput("in_step21_isoforms_rna_y_break", label = "RNA Y breaks", value = 2, min = 2, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step21_isoforms_ribo_x_min", label = "Ribo X min", value = 1, min = 1, step = 1)),
                             column(4, numericInput("in_step21_isoforms_ribo_x_max", label = "Ribo X max", value = NA, min = 1, step = 1)),
                             column(4, numericInput("in_step21_isoforms_ribo_x_break", label = "Ribo X breaks", value = 5, min = 3, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step21_isoforms_ribo_y_min", label = "Ribo Y min", value = 0, min = 0, step = 1)),
                             column(4, numericInput("in_step21_isoforms_ribo_y_max", label = "Ribo Y max", value = NA, min = 0, step = 1)),
                             column(4, numericInput("in_step21_isoforms_ribo_y_break", label = "Ribo Y breaks", value = 2, min = 2, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step21_isoforms_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(4, sliderInput("in_step21_isoforms_grid_width", label = "Grid width", value = 0.1, min = 0, max = 1, step = 0.05)),
                             column(4, sliderInput("in_step21_isoforms_line_width", label = "Line width", value = 1, min = 0, max = 5, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step21_isoforms_color", label = "Color", selected = "Paired", choices = color_list)),
                             column(4, sliderInput("in_step21_isoforms_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.05)),
                             column(4, selectizeInput("in_step21_isoforms_y_scale", label = "Y scale", selected = "none", 
                                                      choices = c("none", "sqrt", "log")))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step21_isoforms_top_height", label = "Top height", value = 3)),
                             column(4, numericInput("in_step21_isoforms_bottom_height", label = "Bottom height", value = 1)),
                             column(4, checkboxInput("in_step21_isoforms_facet", label = "Facet", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step21_isoforms_xlabel", label = "Xlabel", value = "Nucleotide position")),
                             column(6, textInput("in_step21_isoforms_ylabel", label = "Ylabel", value = "Abundance (RPM)"))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step21_isoforms_out_name", label = "Output", value = "Isoforms-reads", placeholder = "output prefix")),
                             column(4, numericInput("out_step21_isoforms_fig_width", label = "Figure width", value = 8, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step21_isoforms_fig_height", label = "Figure height", value = 8, min = 1, max = 20, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step21_isoforms_draw_rna_seq", label = "draw rna-seq", icon = icon("chart-line"))),
                             column(6, actionButton("act_step21_isoforms_draw_ribo_seq", label = "draw ribo-seq", icon = icon("chart-line")))
                           ),
                           
                           tags$style(HTML("#act_step21_isoforms_draw_rna_seq {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#act_step21_isoforms_draw_ribo_seq {background-color: #e78c45; color: black;}")),
                           
                           fluidRow(
                             column(6, downloadButton("save_step21_isoforms_draw_rna_seq", label = "Save RNA plot",
                                                      class = "download-btn", icon = icon("download"))),
                             column(6, downloadButton("save_step21_isoforms_draw_ribo_seq", label = "Save Ribo plot",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step21_isoforms_draw_rna_seq {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step21_isoforms_draw_ribo_seq {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "RNA table",
                                  h4("Show the RNA-seq reads table"),
                                  DT::dataTableOutput("out_step21_isoforms_rna_seq")),
                         
                         tabPanel(title = "Ribo table",
                                  h4("Show the Ribo-seq reads table"),
                                  DT::dataTableOutput("out_step21_isoforms_ribo_seq")),
                         
                         tabPanel(title = "RNA plot",
                                  h4("Show the RNA-seq reads plot."),
                                  conditionalPanel(condition = ("input.act_step21_isoforms_draw_rna_seq > 0"),
                                                   plotOutput("out_step21_isoforms_rna_seq_plot") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Ribo plot",
                                  h4("Show the Ribo-seq reads plot."),
                                  conditionalPanel(condition = ("input.act_step21_isoforms_draw_ribo_seq > 0"),
                                                   plotOutput("out_step21_isoforms_ribo_seq_plot") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      
      ############################################################################################################
      ## step22 Pausing ###############################
      ### step22 draw the pausing score  ###############################
      tabItem(tabName = "step22_pausing",
              h2("Draw the pausing score."),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the pausing score table
                       box(title = "Import the pausing score", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step22_pausing", label = "Pausing score", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "pausing score table"),
                           
                           actionButton("act_step22_import_pausing", label = "Import the pausing score", icon = icon("file")),
                           tags$style(HTML("#act_step22_import_pausing {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Set the groups", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, width = 14,
                           
                           selectizeInput("in_step22_column", label = "Design groups", choices = NULL, selected = NULL, multiple = FALSE),
                           selectizeInput("in_step22_groups", label = "Groups", choices = NULL, selected = NULL, multiple = TRUE)
                           # selectizeInput("in_step22_samples", label = "Samples", choices = NULL, selected = NULL, multiple = TRUE)
                       ),
                       
                       box(title = "Draw pausing score", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step22_valid", label = "Pausing class", selected = 'Total',
                                                      choices = c("Total", "Valid"),  multiple = FALSE)),
                             column(4, checkboxInput("in_step22_wrap", label = "Wrap plot", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step22_wrap_row", label = "Warp rows", value = 3, min = 0)),
                             column(4, selectizeInput("in_step22_wrap_scales", label = "Scales", selected = "free",
                                                      choices = c("fixed", "free", "free_x", "free_y"))),
                             column(4, numericInput("in_step22_legend_row", label = "Legend rows", value = 3, min = 1))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step22_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(4, numericInput("in_step22_grid_width", label = "Grid width", value = 0.2, min = 0, step = 0.05)),
                             column(4, textInput("in_step22_grid_color", label = "Grid color", value = "grey85"))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step22_dot_color", label = "Dot color", selected = "Paired", 
                                                      choices = color_list)),
                             column(4, sliderInput("in_step22_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step22_dot_size", label = "Dot size", value = 4, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step22_line_size", label = "Line size", value = 0.8, min = 0, max = 3, step = 0.1))
                           ),
                           # 
                           # fluidRow(
                           #   column(4, textInput("out_step22_xlabel", label = "Xlab", value = "Codon")),
                           #   column(4, textInput("out_step22_ylabel", label = "Ylab", value = "Pausing score"))
                           # ),
                           
                           fluidRow(
                             column(4, textInput("out_step22_out_name", label = "Output", value = "codon", placeholder = "output prefix")),
                             column(4, numericInput("out_step22_fig_width", label = "Figure width", value = 10, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step22_fig_height", label = "Figure height", value = 8, min = 1, max = 20, step = 0.5))
                           ),
                           
                           actionButton("act_step22_draw_pausing_score", label = "draw pausing score", icon = icon("chart-line")),
                           tags$style(HTML("#act_step22_draw_pausing_score {background-color: #e78c45; color: black;}")),
                           
                           fluidRow(
                             column(6, downloadButton("save_step22_abs_pausing_plot", label = "Save absolute pausing", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(6, downloadButton("save_step22_rel_pausing_plot", label = "Save relative pausing", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step22_abs_pausing_plot {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step22_rel_pausing_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       ),
                       
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Absolute pausing",
                                  h4("Show the absolute total pausing score table"),
                                  DT::dataTableOutput("out_step22_absolute_pausing")),
                         
                         tabPanel(title = "Relative pausing",
                                  h4("Show the relative valid pausing score table"),
                                  DT::dataTableOutput("out_step22_relative_pausing")),
                         
                         tabPanel(title = "Absolute dotplot",
                                  h4("Show the absolute pausing score plot."),
                                  conditionalPanel(condition = ("input.act_step22_draw_pausing_score > 0"),
                                                   plotOutput("out_step22_abs_pausing_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Relative dotplot",
                                  h4("Show the relative pausing score plot."),
                                  conditionalPanel(condition = ("input.act_step22_draw_pausing_score > 0"),
                                                   plotOutput("out_step22_rel_pausing_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      ### step22 draw the differential pausing score  ###############################
      tabItem(tabName = "step22_diff_pausing",
              h2("Draw the differential codon pausing score."),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the pausing score table
                       box(title = "Import the pausing score", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step22_diff_pausing", label = "Pausing score", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "pausing score table"),
                           
                           actionButton("act_step22_diff_import_pausing", label = "import the pausing score", icon = icon("file")),
                           tags$style(HTML("#act_step22_diff_import_pausing {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Calaulate the differential codon pausing score", status = "primary",
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, width = 14,
                           
                           fluidRow(
                             column(6, selectizeInput("in_step22_diff_column", label = "Design groups", 
                                                      choices = NULL, selected = NULL, multiple = FALSE)),
                             column(6, selectizeInput("in_step22_diff_method", label = "Method", selected = 'Delta',
                                                      choices = c("Delta", "FoldChange"),  multiple = FALSE)),
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step22_diff_groups1", label = "Groups 1", choices = NULL, selected = NULL, multiple = FALSE)),
                             column(6, selectizeInput("in_step22_diff_groups2", label = "Groups 2", choices = NULL, selected = NULL, multiple = FALSE))
                           ),
                           
                           actionButton("act_step22_diff_calculate", label = "calc. diff. pausing score", icon = icon("file")),
                           tags$style(HTML("#act_step22_diff_calculate {background-color: #75aadb; color: black;}"))
                           
                       ),
                       
                       box(title = "Draw the differential pausing score", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step22_diff_valid", label = "Pausing class", selected = 'Total',
                                                      choices = c("Total", "Valid"),  multiple = FALSE)),
                             column(4, checkboxInput("in_step22_diff_wrap", label = "Wrap plot", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step22_diff_wrap_row", label = "Warp rows", value = 3, min = 0)),
                             column(4, selectizeInput("in_step22_diff_wrap_scales", label = "Scales", selected = "free",
                                                      choices = c("free", "free_x", "free_y"))),
                             column(4, numericInput("in_step22_diff_legend_row", label = "Legend rows", value = 3, min = 1))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step22_diff_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(4, numericInput("in_step22_diff_grid_width", label = "Grid width", value = 0.2, min = 0, step = 0.05)),
                             column(4, textInput("in_step22_diff_grid_color", label = "Grid color", value = "grey85"))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step22_diff_dot_color", label = "Dot color", selected = "Paired", 
                                                      choices = color_list)),
                             column(4, sliderInput("in_step22_diff_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step22_diff_dot_size", label = "Dot size", value = 4, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step22_diff_line_size", label = "Line size", value = 0.8, min = 0, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step22_diff_out_name", label = "Output", value = "codon", placeholder = "output prefix")),
                             column(4, numericInput("out_step22_diff_fig_width", label = "Figure width", value = 10, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step22_diff_fig_height", label = "Figure height", value = 8, min = 1, max = 20, step = 0.5))
                           ),
                           
                           actionButton("act_step22_diff_draw_pausing_score", label = "draw diff pausing score", icon = icon("chart-line")),
                           tags$style(HTML("#act_step22_diff_draw_pausing_score {background-color: #e78c45; color: black;}")),
                           
                           fluidRow(
                             column(6, downloadButton("save_step22_diff_abs_pausing_plot", label = "Save diff absolute pausing", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(6, downloadButton("save_step22_diff_rel_pausing_plot", label = "Save diff relative pausing", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step22_diff_abs_pausing_plot {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step22_diff_rel_pausing_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Raw pausing",
                                  h4("Show the pausing score table"),
                                  DT::dataTableOutput("out_step22_raw_pausing")),
                         
                         tabPanel(title = "Differential absolute pausing",
                                  h4("Show the absolute total pausing score table"),
                                  DT::dataTableOutput("out_step22_diff_absolute_pausing")),
                         
                         tabPanel(title = "Differential relative pausing",
                                  h4("Show the relative valid pausing score table"),
                                  DT::dataTableOutput("out_step22_diff_relative_pausing")),
                         
                         tabPanel(title = "Absolute dotplot",
                                  h4("Show the differential absolute pausing score plot."),
                                  conditionalPanel(condition = ("input.act_step22_diff_draw_pausing_score > 0"),
                                                   plotOutput("out_step22_diff_abs_pausing_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Relative dotplot",
                                  h4("Show the differential relative pausing score plot."),
                                  conditionalPanel(condition = ("input.act_step22_diff_draw_pausing_score > 0"),
                                                   plotOutput("out_step22_diff_rel_pausing_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      
      ############################################################################################################
      ## step23 Occupancy ###############################
      ### step23 draw the Occupancy  ###############################
      tabItem(tabName = "step23_occupancy",
              h2("Draw the occupancy"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the occupancy table
                       box(title = "Import the occupancy", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step23_occupancy", label = "occupancy", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "occupancy table"),
                           
                           actionButton("act_step23_import_occupancy", label = "import the occupancy", icon = icon("file")),
                           tags$style(HTML("#act_step23_import_occupancy {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Set the groups", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           selectizeInput("in_step23_column", label = "Design groups", choices = NULL, selected = NULL, multiple = FALSE),
                           selectizeInput("in_step23_groups", label = "Groups", choices = NULL, selected = NULL, multiple = TRUE)
                           # selectizeInput("in_step22_samples", label = "Samples", choices = NULL, selected = NULL, multiple = TRUE)
                       ),
                       
                       box(title = "Draw occupancy", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, checkboxInput("in_step23_wrap", label = "Wrap plot", value = TRUE)),
                             column(4, numericInput("in_step23_wrap_row", label = "Warp rows", value = 3, min = 0)),
                             column(4, selectizeInput("in_step23_wrap_scales", label = "Scales", selected = "free",
                                                      choices = c("fixed", "free", "free_x", "free_y")))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step23_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(4, numericInput("in_step23_legend_row", label = "Legend rows", value = 3, min = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step23_grid_width", label = "Grid width", value = 0.2, min = 0, step = 0.05)),
                             column(4, textInput("in_step23_grid_color", label = "Grid color", value = "grey85"))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step23_dot_color", label = "Dot color", selected = "Paired", choices = color_list)),
                             column(4, sliderInput("in_step23_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step23_dot_size", label = "Dot size", value = 4, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step23_line_size", label = "Line size", value = 0.8, min = 0, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step23_out_name", label = "Output", value = "codon", placeholder = "output prefix")),
                             column(4, numericInput("out_step23_fig_width", label = "Figure width", value = 10, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step23_fig_height", label = "Figure height", value = 8, min = 1, max = 20, step = 0.5))
                           ),
                           
                           actionButton("act_step23_draw_occupancy", label = "draw occupancy", icon = icon("chart-line")),
                           tags$style(HTML("#act_step23_draw_occupancy {background-color: #e78c45; color: black;}")),
                           
                           fluidRow(
                             column(6, downloadButton("save_step23_absolute_occupancy_plot", label = "Save absolute occupancy", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(6, downloadButton("save_step23_relative_occupancy_plot", label = "Save relative occupancy", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step23_absolute_occupancy_plot {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step23_relative_occupancy_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Absolute occupancy",
                                  h4("Show the occupancy table"),
                                  DT::dataTableOutput("out_step23_absolute_occupancy")),
                         
                         tabPanel(title = "Relative occupancy",
                                  h4("Show the relative occupancy table"),
                                  DT::dataTableOutput("out_step23_relative_occupancy")),
                         
                         tabPanel(title = "Absolute dotplot",
                                  h4("Show the occupancy plot."),
                                  conditionalPanel(condition = ("input.act_step23_draw_occupancy > 0"),
                                                   plotOutput("out_step23_absolute_occupancy_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Relative dotplot",
                                  h4("Show the relative occupancy plot."),
                                  conditionalPanel(condition = ("input.act_step23_draw_occupancy > 0"),
                                                   plotOutput("out_step23_relative_occupancy_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      ### step23 draw the differential Occupancy  ###############################
      tabItem(tabName = "step23_diff_occupancy",
              h2("Draw the differential occupancy."),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the codon occupancy table
                       box(title = "Import the occupancy", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step23_diff_occupancy", label = "Occupancy", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "occupancy table"),
                           
                           actionButton("act_step23_diff_import_occupancy", label = "import the occupancy", icon = icon("file")),
                           tags$style(HTML("#act_step23_diff_import_occupancy {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Calaulate the differential occupancy", status = "primary",
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, width = 14,
                           
                           fluidRow(
                             column(6, selectizeInput("in_step23_diff_column", label = "Design groups", 
                                                      choices = NULL, selected = NULL, multiple = FALSE)),
                             column(6, selectizeInput("in_step23_diff_method", label = "Method", selected = 'delta',
                                                      choices = c("Delta", "FoldChange"),  multiple = FALSE)),
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step23_diff_groups1", label = "Groups 1", choices = NULL, selected = NULL, multiple = FALSE)),
                             column(6, selectizeInput("in_step23_diff_groups2", label = "Groups 2", choices = NULL, selected = NULL, multiple = FALSE))
                           ),
                           
                           actionButton("act_step23_diff_calculate", label = "calc. diff. occupancy", icon = icon("file")),
                           tags$style(HTML("#act_step23_diff_calculate {background-color: #75aadb; color: black;}"))
                           
                       ),
                       
                       box(title = "Draw the differential occupancy", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, checkboxInput("in_step23_diff_wrap", label = "Wrap plot", value = TRUE)),
                             column(4, numericInput("in_step23_diff_wrap_row", label = "Warp rows", value = 3, min = 0)),
                             column(4, selectizeInput("in_step23_diff_wrap_scales", label = "Scales", selected = "free",
                                                      choices = c("free", "free_x", "free_y")))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step23_diff_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(4, numericInput("in_step23_diff_legend_row", label = "Legend rows", value = 3, min = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step23_diff_grid_width", label = "Grid width", value = 0.2, min = 0, step = 0.05)),
                             column(4, textInput("in_step23_diff_grid_color", label = "Grid color", value = "grey85"))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step23_diff_dot_color", label = "Dot color", selected = "Paired", 
                                                      choices = color_list)),
                             column(4, sliderInput("in_step23_diff_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step23_diff_dot_size", label = "Dot size", value = 4, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step23_diff_line_size", label = "Line size", value = 0.8, min = 0, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step23_diff_out_name", label = "Output", value = "codon", placeholder = "output prefix")),
                             column(4, numericInput("out_step23_diff_fig_width", label = "Figure width", value = 10, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step23_diff_fig_height", label = "Figure height", value = 8, min = 1, max = 20, step = 0.5))
                           ),
                           
                           actionButton("act_step23_diff_draw_occupancy", label = "draw diff occupancy", icon = icon("chart-line")),
                           tags$style(HTML("#act_step23_diff_draw_occupancy {background-color: #e78c45; color: black;}")),
                           
                           fluidRow(
                             column(6, downloadButton("save_step23_diff_abs_occupancy_plot", label = "Save diff absolute occupancy", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(6, downloadButton("save_step23_diff_rel_occupancy_plot", label = "Save diff relative occupancy", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step23_diff_abs_occupancy_plot {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step23_diff_rel_occupancy_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Raw occupancy",
                                  h4("Show the occupancy table"),
                                  DT::dataTableOutput("out_step23_raw_occupancy")),
                         
                         tabPanel(title = "Differential absolute occupancy",
                                  h4("Show the absolute occupancy table"),
                                  DT::dataTableOutput("out_step23_diff_absolute_occupancy")),
                         
                         tabPanel(title = "Differential relative occupancy",
                                  h4("Show the relative occupancy table"),
                                  DT::dataTableOutput("out_step23_diff_relative_occupancy")),
                         
                         tabPanel(title = "Absolute dotplot",
                                  h4("Show the differential absolute occupancy plot."),
                                  conditionalPanel(condition = ("input.act_step23_diff_draw_occupancy > 0"),
                                                   plotOutput("out_step23_diff_abs_occupancy_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Relative dotplot",
                                  h4("Show the differential relative occupancy plot."),
                                  conditionalPanel(condition = ("input.act_step23_diff_draw_occupancy > 0"),
                                                   plotOutput("out_step23_diff_rel_occupancy_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      
      ############################################################################################################
      ## step24 CDT ###############################
      ### step24 draw the CDT  ###############################
      tabItem(tabName = "step24_cdt",
              h2("Draw the codon decoding time"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the cdt table
                       box(title = "Import the cdt", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step24_cdt", label = "cdt", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "cdt table"),
                           
                           actionButton("act_step24_import_cdt", label = "import the cdt", icon = icon("file")),
                           tags$style(HTML("#act_step24_import_cdt {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Set the groups", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           selectizeInput("in_step24_column", label = "Design groups", choices = NULL, selected = NULL, multiple = FALSE),
                           selectizeInput("in_step24_groups", label = "Groups", choices = NULL, selected = NULL, multiple = TRUE)
                           # selectizeInput("in_step22_samples", label = "Samples", choices = NULL, selected = NULL, multiple = TRUE)
                       ),
                       
                       box(title = "Draw cdt", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step24_norm", label = "CDT class", choices = c("Raw", "Normlized"), 
                                                      selected = 'Normlized', multiple = FALSE)),
                             column(4, checkboxInput("in_step24_wrap", label = "Wrap plot", value = TRUE)),
                             column(4, numericInput("in_step24_legend_row", label = "Legend rows", value = 3, min = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step24_wrap_row", label = "Warp rows", value = 3, min = 0)),
                             column(4, selectizeInput("in_step24_wrap_scales", label = "Scales", selected = "free",
                                                      choices = c("fixed", "free", "free_x", "free_y")))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step24_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(4, numericInput("in_step24_grid_width", label = "Grid width", value = 0.2, min = 0, step = 0.05)),
                             column(4, textInput("in_step24_grid_color", label = "Grid color", value = "grey85"))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step24_dot_color", label = "Dot color", selected = "Paired", 
                                                      choices = color_list)),
                             column(4, sliderInput("in_step24_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step24_dot_size", label = "Dot size", value = 4, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step24_line_size", label = "Line size", value = 0.8, min = 0, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step24_out_name", label = "Output", value = "codon", placeholder = "output prefix")),
                             column(4, numericInput("out_step24_fig_width", label = "Figure width", value = 10, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step24_fig_height", label = "Figure height", value = 8, min = 1, max = 20, step = 0.5))
                           ),
                           
                           actionButton("act_step24_draw_cdt", label = "draw cdt", icon = icon("chart-line")),
                           tags$style(HTML("#act_step24_draw_cdt {background-color: #e78c45; color: black;}")),
                           
                           fluidRow(
                             column(6, downloadButton("save_step24_absolute_cdt_plot", label = "Save absolute cdt", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(6, downloadButton("save_step24_relative_cdt_plot", label = "Save relative cdt", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step24_absolute_cdt_plot {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step24_relative_cdt_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Absolute cdt",
                                  h4("Show the cdt table"),
                                  DT::dataTableOutput("out_step24_absolute_cdt")),
                         
                         tabPanel(title = "Relative cdt",
                                  h4("Show the relative cdt table"),
                                  DT::dataTableOutput("out_step24_relative_cdt")),
                         
                         tabPanel(title = "Absolute dotplot",
                                  h4("Show the cdt plot."),
                                  conditionalPanel(condition = ("input.act_step24_draw_cdt > 0"),
                                                   plotOutput("out_step24_absolute_cdt_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Relative dotplot",
                                  h4("Show the relative cdt plot."),
                                  conditionalPanel(condition = ("input.act_step24_draw_cdt > 0"),
                                                   plotOutput("out_step24_relative_cdt_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      
      ### step24 draw the differential CDT  ###############################
      tabItem(tabName = "step24_diff_cdt",
              h2("Draw the differential cdt."),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the cdt table
                       box(title = "Import the cdt", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step24_diff_cdt", label = "cdt", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "cdt table"),
                           
                           actionButton("act_step24_diff_import_cdt", label = "import the cdt", icon = icon("file")),
                           tags$style(HTML("#act_step24_diff_import_cdt {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Calaulate the differential cdt", status = "primary",
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, width = 14,
                           
                           fluidRow(
                             column(6, selectizeInput("in_step24_diff_column", label = "Design groups", 
                                                      choices = NULL, selected = NULL, multiple = FALSE)),
                             column(6, selectizeInput("in_step24_diff_method", label = "Method", selected = 'delta',
                                                      choices = c("Delta", "FoldChange"),  multiple = FALSE)),
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step24_diff_groups1", label = "Groups 1", choices = NULL, selected = NULL, multiple = FALSE)),
                             column(6, selectizeInput("in_step24_diff_groups2", label = "Groups 2", choices = NULL, selected = NULL, multiple = FALSE))
                           ),
                           
                           actionButton("act_step24_diff_calculate", label = "calc. diff. cdt", icon = icon("file")),
                           tags$style(HTML("#act_step24_diff_calculate {background-color: #75aadb; color: black;}"))
                           
                       ),
                       
                       box(title = "Draw the differential cdt", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, checkboxInput("in_step24_diff_wrap", label = "Wrap plot", value = TRUE)),
                             column(4, numericInput("in_step24_diff_wrap_row", label = "Warp rows", value = 3, min = 0)),
                             column(4, selectizeInput("in_step24_diff_wrap_scales", label = "Scales", selected = "free",
                                                      choices = c("free", "free_x", "free_y")))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step24_diff_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(4, numericInput("in_step24_diff_legend_row", label = "Legend rows", value = 3, min = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step24_diff_grid_width", label = "Grid width", value = 0.2, min = 0, step = 0.05)),
                             column(4, textInput("in_step24_diff_grid_color", label = "Grid color", value = "grey85"))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step24_diff_dot_color", label = "Dot color", selected = "Paired", 
                                                      choices = color_list)),
                             column(4, sliderInput("in_step24_diff_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step24_diff_dot_size", label = "Dot size", value = 4, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step24_diff_line_size", label = "Line size", value = 0.8, min = 0, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step24_diff_out_name", label = "Output", value = "codon", placeholder = "output prefix")),
                             column(4, numericInput("out_step24_diff_fig_width", label = "Figure width", value = 10, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step24_diff_fig_height", label = "Figure height", value = 8, min = 1, max = 20, step = 0.5))
                           ),
                           
                           actionButton("act_step24_diff_draw_cdt", label = "draw diff cdt", icon = icon("chart-line")),
                           tags$style(HTML("#act_step24_diff_draw_cdt {background-color: #e78c45; color: black;}")),
                           
                           fluidRow(
                             column(6, downloadButton("save_step24_diff_abs_cdt_plot", label = "Save diff absolute cdt", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(6, downloadButton("save_step24_diff_rel_cdt_plot", label = "Save diff relative cdt", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step24_diff_abs_cdt_plot {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step24_diff_rel_cdt_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Raw cdt",
                                  h4("Show the cdt table"),
                                  DT::dataTableOutput("out_step24_raw_cdt")),
                         
                         tabPanel(title = "Differential absolute cdt",
                                  h4("Show the absolute total cdt table"),
                                  DT::dataTableOutput("out_step24_diff_absolute_cdt")),
                         
                         tabPanel(title = "Differential relative cdt",
                                  h4("Show the relative valid cdt table"),
                                  DT::dataTableOutput("out_step24_diff_relative_cdt")),
                         
                         tabPanel(title = "Absolute dotplot",
                                  h4("Show the differential absolute cdt plot."),
                                  conditionalPanel(condition = ("input.act_step24_diff_draw_cdt > 0"),
                                                   plotOutput("out_step24_diff_abs_cdt_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Relative dotplot",
                                  h4("Show the differential relative cdt plot."),
                                  conditionalPanel(condition = ("input.act_step24_diff_draw_cdt > 0"),
                                                   plotOutput("out_step24_diff_rel_cdt_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      
      ############################################################################################################
      ## step25 CST ###############################
      ### step25 draw the codon selection time ###############################
      tabItem(tabName = "step25_cst",
              h2("Draw the codon selection time"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the cst table
                       box(title = "Import the cst", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step25_cst", label = "cst", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "cst table"),
                           
                           actionButton("act_step25_import_cst", label = "import the cst", icon = icon("file")),
                           tags$style(HTML("#act_step25_import_cst {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Set the groups", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, width = 14,
                           
                           selectizeInput("in_step25_column", label = "Design groups", choices = NULL, selected = NULL, multiple = FALSE),
                           selectizeInput("in_step25_groups", label = "Groups", choices = NULL, selected = NULL, multiple = TRUE)
                           # selectizeInput("in_step22_samples", label = "Samples", choices = NULL, selected = NULL, multiple = TRUE)
                       ),
                       
                       box(title = "Draw cst", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, checkboxInput("in_step25_wrap", label = "Wrap plot", value = TRUE)),
                             column(4, numericInput("in_step25_wrap_row", label = "Warp rows", value = 3, min = 0))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step25_wrap_scales", label = "Scales", selected = "free", 
                                                      choices = c("fixed", "free", "free_x", "free_y"))),
                             column(4, numericInput("in_step25_legend_row", label = "Legend rows", value = 3, min = 1))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step25_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(4, numericInput("in_step25_grid_width", label = "Grid width", value = 0.2, min = 0, step = 0.05)),
                             column(4, textInput("in_step25_grid_color", label = "Grid color", value = "grey85"))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step25_dot_color", label = "Dot color", selected = "Paired", 
                                                      choices = color_list)),
                             column(4, sliderInput("in_step25_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.05)),
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step25_dot_size", label = "Dot size", value = 4, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step25_line_size", label = "Line size", value = 0.8, min = 0, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step25_out_name", label = "Output", value = "codon", placeholder = "output prefix")),
                             column(4, numericInput("out_step25_fig_width", label = "Figure width", value = 10, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step25_fig_height", label = "Figure height", value = 8, min = 1, max = 20, step = 0.5))
                           ),
                           
                           actionButton("act_step25_draw_cst", label = "draw cst", icon = icon("chart-line")),
                           tags$style(HTML("#act_step25_draw_cst {background-color: #e78c45; color: black;}")),
                           
                           fluidRow(
                             column(6, downloadButton("save_step25_absolute_cst_plot", label = "Save absolute cst", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(6, downloadButton("save_step25_relative_cst_plot", label = "Save relative cst", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step25_absolute_cst_plot {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step25_relative_cst_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Absolute CST",
                                  h4("Show the cst table"),
                                  DT::dataTableOutput("out_step25_absolute_cst")),
                         
                         tabPanel(title = "Relative CST",
                                  h4("Show the relative cst table"),
                                  DT::dataTableOutput("out_step25_relative_cst")),
                         
                         tabPanel(title = "Absolute dotplot",
                                  h4("Show the cst plot."),
                                  conditionalPanel(condition = ("input.act_step25_draw_cst > 0"),
                                                   plotOutput("out_step25_absolute_cst_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Relative dotplot",
                                  h4("Show the relative cst plot."),
                                  conditionalPanel(condition = ("input.act_step25_draw_cst > 0"),
                                                   plotOutput("out_step25_relative_cst_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      
      ### step25 draw the differential cst  ###############################
      tabItem(tabName = "step25_diff_cst",
              h2("Draw the differential cst."),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the cst table
                       box(title = "Import the cst", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step25_diff_cst", label = "cst", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "cst table"),
                           
                           actionButton("act_step25_diff_import_cst", label = "import the cst", icon = icon("file")),
                           tags$style(HTML("#act_step25_diff_import_cst {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Calaulate the differential cst", status = "primary",
                           solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, width = 14,
                           
                           fluidRow(
                             column(6, selectizeInput("in_step25_diff_column", label = "Design groups", 
                                                      choices = NULL, selected = NULL, multiple = FALSE)),
                             column(6, selectizeInput("in_step25_diff_method", label = "Method", selected = 'delta',
                                                      choices = c("Delta", "FoldChange"),  multiple = FALSE)),
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step25_diff_groups1", label = "Groups 1", choices = NULL, selected = NULL, multiple = FALSE)),
                             column(6, selectizeInput("in_step25_diff_groups2", label = "Groups 2", choices = NULL, selected = NULL, multiple = FALSE))
                           ),
                           
                           actionButton("act_step25_diff_calculate", label = "calc. diff. cst", icon = icon("file")),
                           tags$style(HTML("#act_step25_diff_calculate {background-color: #75aadb; color: black;}"))
                           
                       ),
                       
                       box(title = "Draw the differential cst", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, checkboxInput("in_step25_diff_wrap", label = "Wrap plot", value = TRUE)),
                             column(4, numericInput("in_step25_diff_wrap_row", label = "Warp rows", value = 3, min = 0)),
                             column(4, selectizeInput("in_step25_diff_wrap_scales", label = "Scales", selected = "free",
                                                      choices = c("free", "free_x", "free_y")))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step25_diff_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(4, numericInput("in_step25_diff_legend_row", label = "Legend rows", value = 3, min = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step25_diff_grid_width", label = "Grid width", value = 0.2, min = 0, step = 0.05)),
                             column(4, textInput("in_step25_diff_grid_color", label = "Grid color", value = "grey85"))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step25_diff_dot_color", label = "Dot color", selected = "Paired", 
                                                      choices = color_list)),
                             column(4, sliderInput("in_step25_diff_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step25_diff_dot_size", label = "Dot size", value = 4, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step25_diff_line_size", label = "Line size", value = 0.8, min = 0, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step25_diff_out_name", label = "Output", value = "codon", placeholder = "output prefix")),
                             column(4, numericInput("out_step25_diff_fig_width", label = "Figure width", value = 10, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step25_diff_fig_height", label = "Figure height", value = 8, min = 1, max = 20, step = 0.5))
                           ),
                           
                           actionButton("act_step25_diff_draw_cst", label = "draw diff cst", icon = icon("chart-line")),
                           tags$style(HTML("#act_step25_diff_draw_cst {background-color: #e78c45; color: black;}")),
                           
                           fluidRow(
                             column(6, downloadButton("save_step25_diff_abs_cst_plot", label = "Save diff absolute cst", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(6, downloadButton("save_step25_diff_rel_cst_plot", label = "Save diff relative cst", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step25_diff_abs_cst_plot {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step25_diff_rel_cst_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Raw cst",
                                  h4("Show the cst table"),
                                  DT::dataTableOutput("out_step25_raw_cst")),
                         
                         tabPanel(title = "Differential absolute cst",
                                  h4("Show the absolute cst table"),
                                  DT::dataTableOutput("out_step25_diff_absolute_cst")),
                         
                         tabPanel(title = "Differential relative cst",
                                  h4("Show the relative cst table"),
                                  DT::dataTableOutput("out_step25_diff_relative_cst")),
                         
                         tabPanel(title = "Absolute dotplot",
                                  h4("Show the differential absolute cst plot."),
                                  conditionalPanel(condition = ("input.act_step25_diff_draw_cst > 0"),
                                                   plotOutput("out_step25_diff_abs_cst_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Relative dotplot",
                                  h4("Show the differential relative cst plot."),
                                  conditionalPanel(condition = ("input.act_step25_diff_draw_cst > 0"),
                                                   plotOutput("out_step25_diff_rel_cst_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      # 
      # ### step25 draw the iterative codon selection time ###############################
      # tabItem(tabName = "step25_iterative_cst",
      #         h2("Draw the iterative codon selection time"),
      #         fluidRow(
      #           tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
      #           
      #           # set the input panel
      #           column(5,
      #                  ## click to import the cst table
      #                  box(title = "Import the cst", status = "primary", collapsible = TRUE, collapsed = FALSE, 
      #                      solidHeader = TRUE, width = NULL,
      #                      
      #                      fileInput(inputId = "in_step25_iterative_cst", label = "cst", accept = c(".txt", ".TXT"),
      #                                multiple = TRUE, placeholder = "cst table"),
      #                      
      #                      actionButton("act_step25_import_iterative_cst", label = "import the cst", icon = icon("file")),
      #                      tags$style(HTML("#act_step25_import_iterative_cst {background-color: #75aadb; color: black;}"))
      #                  ),
      #                  
      #                  box(title = "Set the groups", status = "primary", collapsible = TRUE, collapsed = FALSE, 
      #                      solidHeader = TRUE, width = NULL,
      #                      
      #                      selectizeInput("in_step25_iterative_column", label = "Design groups", choices = NULL, selected = NULL, multiple = FALSE),
      #                      selectizeInput("in_step25_iterative_groups", label = "Groups", choices = NULL, selected = NULL, multiple = TRUE)
      #                      # selectizeInput("in_step22_samples", label = "Samples", choices = NULL, selected = NULL, multiple = TRUE)
      #                  ),
      #                  
      #                  box(title = "Draw cst", status = "primary", collapsible = TRUE, collapsed = FALSE, 
      #                      solidHeader = TRUE, width = NULL,
      #                      
      #                      fluidRow(
      #                        column(4, checkboxInput("in_step25_iterative_wrap", label = "Wrap plot", value = TRUE)),
      #                        column(4, numericInput("in_step25_iterative_wrap_row", label = "Warp rows", value = 3, min = 0))
      #                      ),
      #                      
      #                      fluidRow(
      #                        column(4, sliderInput("in_step25_iterative_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
      #                        column(4, numericInput("in_step25_iterative_grid_width", label = "Grid width", value = 0.2, min = 0, step = 0.05)),
      #                        column(4, textInput("in_step25_iterative_grid_color", label = "Grid color", value = "grey85"))
      #                      ),
      #                      
      #                      fluidRow(
      #                        column(4, selectizeInput("in_step25_iterative_dot_color", label = "Dot color", selected = "Paired", 
      #                                                 choices = color_list)),
      #                        column(4, sliderInput("in_step25_iterative_dot_size", label = "Dot size", value = 4, min = 0, max = 10, step = 0.5)),
      #                        column(4, sliderInput("in_step25_iterative_line_size", label = "Line size", value = 0.8, min = 0, max = 3, step = 0.1))
      #                      ),
      #                      
      #                      fluidRow(
      #                        column(4, textInput("out_step25_iterative_out_name", label = "Output", value = "codon", placeholder = "output prefix")),
      #                        column(4, numericInput("out_step25_iterative_fig_width", label = "Figure width", value = 10, min = 1, max = 20, step = 0.5)),
      #                        column(4, numericInput("out_step25_iterative_fig_height", label = "Figure height", value = 8, min = 1, max = 20, step = 0.5))
      #                      ),
      #                      
      #                      actionButton("act_step25_iterative_draw_cst", label = "draw iterative cst", icon = icon("chart-line")),
      #                      tags$style(HTML("#act_step25_iterative_draw_cst {background-color: #e78c45; color: black;}")),
      #                      
      #                      fluidRow(
      #                        column(6, downloadButton("save_step25_iterative_absolute_cst_plot", label = "Save absolute cst", 
      #                                                 class = "download-btn", icon = icon("download"))),
      #                        column(6, downloadButton("save_step25_iterative_relative_cst_plot", label = "Save relative cst", 
      #                                                 class = "download-btn", icon = icon("download")))
      #                      ),
      #                      
      #                      tags$style(HTML("#save_step25_iterative_absolute_cst_plot {background-color: #grey; color: black; margin-top: 0px;}")),
      #                      tags$style(HTML("#save_step25_iterative_relative_cst_plot {background-color: #grey; color: black; margin-top: 0px;}"))
      #                      
      #                  )
      #           ),
      #           
      #           column(7,
      #                  tabsetPanel(
      #                    ## set the output panel
      #                    tabPanel(title = "CST",
      #                             h4("Show the iterative cst table"),
      #                             DT::dataTableOutput("out_step25_iterative_absolute_cst")),
      #                    
      #                    tabPanel(title = "Relative CST",
      #                             h4("Show the iterative relative cst table"),
      #                             DT::dataTableOutput("out_step25_iterative_relative_cst")),
      #                    
      #                    tabPanel(title = "Dotplot",
      #                             h4("Show the iterative cst plot."),
      #                             conditionalPanel(condition = ("input.act_step25_iterative_draw_cst > 0"),
      #                                              plotOutput("out_step25_iterative_absolute_cst_plot") %>% 
      #                                                withSpinner(color="#3c8dbc", type = 5, size = 1))),
      #                    
      #                    tabPanel(title = "Relative dotplot",
      #                             h4("Show the iterative relative cst plot."),
      #                             conditionalPanel(condition = ("input.act_step25_iterative_draw_cst > 0"),
      #                                              plotOutput("out_step25_iterative_relative_cst_plot") %>% 
      #                                                withSpinner(color="#3c8dbc", type = 5, size = 1)))
      #                  )
      #           )
      #         )
      # ),
      # 
      
      ############################################################################################################
      ## step26 Odd Ratio ###############################
      ### step26 draw the oddratio ###############################
      tabItem(tabName = "step26_odd_ratio",
              h2("Draw the oddratio"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the oddratio table
                       box(title = "Import the oddratio", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step26_odd_ratio", label = "oddratio", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "oddratio table"),
                           
                           actionButton("act_step26_import_odd_ratio", label = "import the oddratio", icon = icon("file")),
                           tags$style(HTML("#act_step26_import_odd_ratio {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Draw oddratio", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           selectizeInput("in_step26_control", label = "Control", choices = NULL, selected = NULL, multiple = TRUE),
                           selectizeInput("in_step26_treat", label = "Treated", choices = NULL, selected = NULL, multiple = TRUE),

                           fluidRow(
                             column(6, selectizeInput("in_step26_class", label = "Odd Ratio class", choices = c("number", "proportion"), 
                                                      selected = 'proportion', multiple = FALSE)),
                             column(4, checkboxInput("in_step26_wrap", label = "Wrap plot", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step26_wrap_row", label = "Warp rows", value = 3, min = 0)),
                             column(4, selectizeInput("in_step26_wrap_scales", label = "Scales", selected = "free",
                                                      choices = c("fixed", "free", "free_x", "free_y"))),
                             column(4, numericInput("in_step26_legend_row", label = "Legend rows", value = 3, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step26_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(6, sliderInput("in_step26_grid_width", label = "Grid width", value = 0.2, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step26_grid_color", label = "Grid color", value = "grey85")),
                             column(4, selectizeInput("in_step26_dot_color", label = "Dot color", selected = "Paired", choices = color_list)),
                             column(4, numericInput("in_step26_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             
                             column(4, sliderInput("in_step26_dot_size", label = "Dot size", value = 3, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step26_line_size", label = "Line size", value = 0.8, min = 0, max = 5, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step26_out_name", label = "Output", value = "odd_ratio", placeholder = "output prefix")),
                             column(4, numericInput("out_step26_fig_width", label = "Figure width", value = 10, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step26_fig_height", label = "Figure height", value = 8, min = 1, max = 20, step = 0.5))
                           ),
                           
                           actionButton("act_step26_draw_odd_ratio", label = "draw odd ratio", icon = icon("chart-line")),
                           tags$style(HTML("#act_step26_draw_odd_ratio {background-color: #e78c45; color: black;}")),
                           
                           fluidRow(
                             column(6, downloadButton("save_step26_absolute_odd_ratio_plot", label = "Save absolute oddratio", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(6, downloadButton("save_step26_relative_odd_ratio_plot", label = "Save relative oddratio", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#save_step26_absolute_odd_ratio_plot {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step26_relative_odd_ratio_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "OddRatio",
                                  h4("Show the Odd Ratio table"),
                                  DT::dataTableOutput("out_step26_absolute_odd_ratio")),
                         
                         tabPanel(title = "Relative OddRatio",
                                  h4("Show the relative Odd Ratio table"),
                                  DT::dataTableOutput("out_step26_relative_odd_ratio")),
                         
                         tabPanel(title = "Dotplot",
                                  h4("Show the Odd Ratio plot."),
                                  conditionalPanel(condition = ("input.act_step26_draw_odd_ratio > 0"),
                                                   plotOutput("out_step26_absolute_odd_ratio_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Relative dotplot",
                                  h4("Show the relative Odd Ratio plot."),
                                  conditionalPanel(condition = ("input.act_step26_draw_odd_ratio > 0"),
                                                   plotOutput("out_step26_relative_odd_ratio_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      
      ############################################################################################################
      ## step27 CoV ###############################
      ### step27 draw the coefficient of variation (CoV) #########################################################
      tabItem(tabName = "step27_cov",
              h2("Draw the coefficient of variation (CoV)"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the cov table
                       box(title = "Import the CoV", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step27_cov", label = "CoV", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "CoV table"),
                           
                           actionButton("act_step27_import_cov", label = "import the CoV", icon = icon("file")),
                           tags$style(HTML("#act_step27_import_cov {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Fitted the CoV", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step27_cov_column", label = "Design groups", choices = NULL, 
                                                      selected = NULL, multiple = FALSE)),
                             column(8, selectizeInput("in_step27_cov_groups", label = "Groups", choices = NULL, 
                                                      selected = NULL, multiple = TRUE))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step27_cov_sum", label = "Sum", value = 5, step = 1)),
                             column(6, numericInput("in_step27_cov_mean", label = "Mean", value = 0, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step27_cov_filter", label = "filter the CoV", icon = icon("calculator"))),
                             column(6, actionButton("act_step27_cov_fitted", label = "fitted the CoV", icon = icon("calculator")))
                           ),
                           
                           tags$style(HTML("#act_step27_cov_filter {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#act_step27_cov_fitted {background-color: #e78c45; color: black;}"))
                       ),
                       
                       box(title = "Draw cov", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(6, checkboxInput("in_step27_cov_facet", label = "Wrap", value = FALSE)),
                             column(6, numericInput("in_step27_cov_wrap_num", label = "Wrap number", value = 3))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step27_cov_figure_type", label = "Figure type", selected = "scatter",
                                                      choices = c("scatter", "fitted", "merge"))),
                             column(6, sliderInput("in_step27_cov_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step27_cov_dot_color", label = "Dot color", selected = "paired", choices = color_list)),
                             column(4, sliderInput("in_step27_cov_dot_size", label = "Dot size", value = 1, min = 0, max = 5, step = 0.5)),
                             column(4, sliderInput("in_step27_cov_dot_alpha", label = "Dot alpha", value = 0.3, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step27_cov_line_color", label = "Line color", selected = "paired", choices = color_list)),
                             column(4, sliderInput("in_step27_cov_line_width", label = "Line width", value = 1, min = 0, max = 5, step = 0.1)),
                             column(4, sliderInput("in_step27_cov_line_alpha", label = "Line alpha", value = 0.9, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step27_cov_x_min", label = "X min", value = NULL, step = 1)),
                             column(4, numericInput("in_step27_cov_x_max", label = "X max", value = NULL, step = 1)),
                             column(4, numericInput("in_step27_cov_x_breaks", label = "X breaks", value = 5, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step27_cov_y_min", label = "Y min", value = NULL, step = 1)),
                             column(4, numericInput("in_step27_cov_y_max", label = "Y max", value = NULL, step = 1)),
                             column(4, numericInput("in_step27_cov_y_breaks", label = "Y breaks", value = 5, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step27_cov_title", label = "Title", value = "Coefficient of Variation")),
                             column(4, textInput("in_step27_cov_xlabel", label = "Xlabel", value = "log2(Mean)")),
                             column(4, textInput("in_step27_cov_ylabel", label = "Ylabel", value = "log2(CoV)"))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step27_cov_out_name", label = "Output", value = "CoV", placeholder = "output prefix")),
                             column(4, numericInput("out_step27_cov_fig_width", label = "Figure width", value = 10, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step27_cov_fig_height", label = "Figure height", value = 7, min = 1, max = 20, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step27_cov_draw", label = "draw CoV plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step27_cov_scatter_plot", label = "Save CoV plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step27_cov_draw {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step27_cov_scatter_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "CoV",
                                  h4("Show the CoV table"),
                                  DT::dataTableOutput("out_step27_cov")),
                         
                         tabPanel(title = "flt_CoV",
                                  h4("Show the CoV table"),
                                  conditionalPanel(condition = ("input.act_step27_cov_filter > 0"),
                                                   DT::dataTableOutput("out_step27_cov_filtered") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "Fitted_CoV",
                                  h4("Show the fitted CoV results."),
                                  conditionalPanel(condition = ("input.act_step27_cov_fitted > 0"),
                                                   DT::dataTableOutput("out_step27_cov_fitted_results") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "CoV plot",
                                  h4("Show the CoV plot."),
                                  conditionalPanel(condition = ("input.act_step27_cov_draw > 0"),
                                                   plotOutput("out_step27_cov_scatter_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),

      ### step27 draw the cumulative of coefficient of variation (CoV) #########################################################
      tabItem(tabName = "step27_cov_cdf",
              h2("Draw the cumulative of CoV"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),

                # set the input panel
                column(5,
                       ## click to import the CoV eCDF table
                       box(title = "Import the CoV eCDF", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,

                           fileInput(inputId = "in_step27_cov_cdf", label = "CoV eCDF", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "CoV eCDF file"),
                           
                           actionButton("act_step27_cov_cdf_import", label = "Import the CoV eCDF", icon = icon("file")),

                           tags$style(HTML("#act_step27_cov_cdf_import {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Set the groups", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, width = 14,
                           
                           selectizeInput("in_step27_cov_cdf_gene", label = "Gene name", selected = "", 
                                          choices = NULL, multiple = FALSE),

                           selectizeInput("in_step27_cov_cdf_column", label = "Column", choices = NULL, 
                                          selected = NULL, multiple = FALSE),
                           
                           selectizeInput("in_step27_cov_cdf_groups", label = "Groups", choices = NULL, 
                                          selected = NULL, multiple = TRUE),
                           
                           
                           actionButton("act_step27_cov_cdf_filter", label = "Filter the CoV eCDF", icon = icon("calculator")),
                           
                           tags$style(HTML("#act_step27_cov_cdf_filter {background-color: #75aadb; color: black;}"))
                           
                       ),

                       box(title = "Draw the CDF of CoV", status = "primary", collapsible = TRUE, collapsed = FALSE,
                           solidHeader = TRUE, width = NULL,

                           fluidRow(
                             column(4, numericInput("in_step27_cov_cdf_x_min", label = "Ribo X min", value = 1, min = 1, step = 1)),
                             column(4, numericInput("in_step27_cov_cdf_x_max", label = "Ribo X max", value = NA, min = 1, step = 1)),
                             column(4, numericInput("in_step27_cov_cdf_x_break", label = "Ribo X breaks", value = 5, min = 3, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step27_cov_cdf_y_min", label = "Ribo Y min", value = 0, min = 0, step = 1)),
                             column(4, numericInput("in_step27_cov_cdf_y_max", label = "Ribo Y max", value = NA, min = 0, step = 1)),
                             column(4, numericInput("in_step27_cov_cdf_y_break", label = "Ribo Y breaks", value = 2, min = 2, step = 1))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step27_cov_cdf_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(6, sliderInput("in_step27_cov_cdf_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step27_cov_cdf_grid_width", label = "Grid width", value = 0.1, min = 0, max = 1, step = 0.05)),
                             column(4, sliderInput("in_step27_cov_cdf_line_width", label = "Line width", value = 0.8, min = 0, max = 5, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step27_cov_cdf_legend", label = "Legend", selected = "none",
                                                      choices = c("none", "right", "bottom", "top"))),
                             column(6, selectizeInput("in_step27_cov_cdf_y_scale", label = "Y scale", selected = "none",
                                                      choices = c("none", "sqrt", "log")))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step27_cov_cdf_color", label = "Color", selected = "Paired", choices = color_list)),
                             column(6, checkboxInput("in_step27_cov_cdf_facet", label = "Facet", value = TRUE))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step27_cov_cdf_xlabel", label = "Xlabel", value = "Nucleotide position")),
                             column(6, textInput("in_step27_cov_cdf_ylabel", label = "Ylabel", value = "Cumulative CoV"))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step27_cov_cdf_out_name", label = "Output", value = "CoV-CDF", placeholder = "output prefix")),
                             column(4, numericInput("out_step27_cov_cdf_fig_width", label = "Figure width", value = 8, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step27_cov_cdf_fig_height", label = "Figure height", value = 8, min = 1, max = 20, step = 0.5))
                           ),

                           fluidRow(
                             column(6, actionButton("act_step27_cov_cdf_draw", label = "draw CoV CDF", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step27_cov_cdf_plot", label = "Save CoV CDF",
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step27_cov_cdf_draw {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step27_cov_cdf_plot {background-color: #grey; color: black; margin-top: 0px;}"))
                       )
                ),

                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "CoV_CDF",
                                  h4("Show the CoV CDF table"),
                                  DT::dataTableOutput("out_step27_cov_cdf")),

                         tabPanel(title = "Filtered_CoV_CDF",
                                  h4("Show the filtered CoV CDF table"),
                                  DT::dataTableOutput("out_step27_cov_cdf_filtered")),
                         
                         tabPanel(title = "CoV CDF plot",
                                  h4("Show the CoV CDF plot."),
                                  conditionalPanel(condition = ("input.act_step27_cov_cdf_draw > 0"),
                                                   plotOutput("out_step27_cov_cdf_plot") %>%
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step28 Meta Codon ###############################
      ### step28 draw the metaplot of codon pausing ##############################################################
      tabItem(tabName = "step28_meta_codon_plot",
              h2("Draw the meta codon plot"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the meta codon table
                       box(title = "Import the meta codon density", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step28_meta_codon", label = "meta codon", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "meta codon table"),
                           
                           actionButton("act_step28_import_meta_codon", label = "import the density", icon = icon("file")),
                           tags$style(HTML("#act_step28_import_meta_codon {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Set the groups", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           selectizeInput("in_step28_meta_codon_column", label = "Design groups", choices = NULL, selected = NULL, multiple = FALSE),
                           selectizeInput("in_step28_meta_codon_groups", label = "Groups", choices = NULL, selected = NULL, multiple = TRUE)
                       ),
                       
                       box(title = "Draw meta codon", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step28_meta_codon_site", label = "Position", selected = "Codon", 
                                                      choices = c("Codon", "Nucleotide"))),
                             column(4, selectizeInput("in_step28_meta_codon_scaled", label = "Scaled", selected = "Density", 
                                                      choices = c("Scaled", "Density")))
                           ),
                           
                           fluidRow(
                             column(4, checkboxInput("in_step28_meta_codon_wrap", label = "Wrap plot", value = TRUE)),
                             column(4, numericInput("in_step28_meta_codon_wrap_num", label = "Wrap num", value = 3, step = 1)),
                             column(4, numericInput("in_step28_meta_codon_label", label = "Label", value = 0))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step28_meta_codon_grid_width", label = "Grid width", value = 0.2, min = 0, step = 0.05)),
                             column(4, textInput("in_step28_meta_codon_grid_color", label = "Grid color", value = "grey85")),
                             column(4, numericInput("in_step28_meta_codon_grid_width", label = "Grid width", value = 0.2, min = 0, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step28_meta_codon_x_min", label = "X min", value = NA, step = 1)),
                             column(4, numericInput("in_step28_meta_codon_x_max", label = "X max", value = NA, step = 1)),
                             column(4, numericInput("in_step28_meta_codon_x_breaks", label = "X breaks", value = 5, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step28_meta_codon_y_min", label = "Y min", value = NA, step = 0.1)),
                             column(4, numericInput("in_step28_meta_codon_y_max", label = "Y max", value = NA, step = 0.1)),
                             column(4, numericInput("in_step28_meta_codon_y_breaks", label = "Y breaks", value = 3, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step28_meta_codon_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(4, sliderInput("in_step28_meta_codon_dot_size", label = "Dot size", value = 1, min = 0, max = 10, step = 0.5)),
                             column(4, sliderInput("in_step28_meta_codon_line_size", label = "Line size", value = 0.8, min = 0, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step28_meta_codon_line_color", label = "Line color", selected = "Paired",
                                                      choices = color_list)),
                             column(6, selectizeInput("in_step28_meta_codon_heat_color", label = "Heat color", selected = "Blues",
                                                      choices = color_list)),
                             column(6, sliderInput("in_step28_meta_codon_alpha", label = "Alpha", value = 0.9, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step28_meta_codon_out_name", label = "Output", value = "codon", placeholder = "output prefix")),
                             column(4, numericInput("out_step28_meta_codon_fig_width", label = "Figure width", value = 8, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step28_meta_codon_fig_height", label = "Figure height", value = 4, min = 1, max = 20, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(4, actionButton("act_step28_meta_codon_draw_plot", label = "draw meta plot", icon = icon("chart-line"))),
                             column(4, downloadButton("save_step28_meta_codon_line", label = "Save line plot", 
                                                      class = "download-btn", icon = icon("download"))),
                             column(4, downloadButton("save_step28_meta_codon_heat", label = "Save heat map", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step28_meta_codon_draw_plot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step28_meta_codon_line {background-color: #grey; color: black; margin-top: 0px;}")),
                           tags$style(HTML("#save_step28_meta_codon_heat {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "MetaCodon",
                                  h4("Show the meta codon table"),
                                  DT::dataTableOutput("out_step28_meta_codon_dst")),
                         
                         tabPanel(title = "LinePlot",
                                  h4("Show the line plot."),
                                  conditionalPanel(condition = ("input.act_step28_meta_codon_draw_plot > 0"),
                                                   plotOutput("out_step28_meta_codon_line_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1))),
                         
                         tabPanel(title = "HeatMap",
                                  h4("Show the heat map."),
                                  conditionalPanel(condition = ("input.act_step28_meta_codon_draw_plot > 0"),
                                                   plotOutput("out_step28_meta_codon_heat_map") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      tabItem(tabName = "step28_meta_codon_seqlogo",
              h2("Draw the meta codon seq logo"),
              fluidRow(
                tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the meta codon table
                       box(title = "Import the meta codon sequence", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step28_meta_seq", label = "meta codon sequence", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "meta codon table"),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step28_meta_seq_seqtype", label = "Seq type", 
                                                      choices = c('dna', 'rna', 'aa'), selected = "dna")),
                             column(6, selectizeInput("in_step28_meta_seq_code", label = "Genetic code", selected = "SGC0 - Standard", 
                                                      choices = genetic_code_list))
                           ),
                           
                           actionButton("act_step28_import_meta_seq", label = "import the sequence", icon = icon("file")),
                           
                           tags$style(HTML("#act_step28_import_meta_seq {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Draw meta codon seqlogo", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(6, selectizeInput("in_step28_meta_seq_method", label = "Method", 
                                                      choices = c('bits', 'probability', 'custom'), selected = "probability")),
                             column(6, selectizeInput("in_step28_meta_seq_font_family", label = "Font family", 
                                                      choices = seqlogo_font_family, selected = "helvetica_regular")),
                             conditionalPanel(condition = "input.in_step28_meta_seq_method == 'custom'",
                                              column(4, selectizeInput("in_step28_meta_seq_custom", label = "Custom", selected = "foldchange", 
                                                                       choices = c("delta", "foldchange", "enrichment"))))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step28_meta_seq_font_stack", label = "Stack width", value = 0.9, min = 0, max = 1, step = 0.1)),
                             column(6, sliderInput("in_step28_meta_seq_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step28_meta_seq_fill_color", label = "Fill color", selected = "auto", choices = seqlogo_colors)),
                             column(6, sliderInput("in_step28_meta_seq_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step28_meta_seq_x_min", label = "X min", value = NA, step = 1)),
                             column(6, numericInput("in_step28_meta_seq_x_max", label = "X max", value = NA, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step28_meta_seq_out_name", label = "Output", value = "codon", placeholder = "output prefix")),
                             column(4, numericInput("out_step28_meta_seq_fig_width", label = "Figure width", value = 10, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step28_meta_seq_fig_height", label = "Figure height", value = 6, min = 1, max = 20, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step28_meta_seq_draw_plot", label = "draw meta plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step28_meta_seq_logo", label = "Save seq logo", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step28_meta_seq_draw_plot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step28_meta_seq_logo {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "MetaCodon",
                                  h4("Show the meta codon sequence"),
                                  DT::dataTableOutput("out_step28_meta_seq")),
                         
                         tabPanel(title = "SeqLogo",
                                  h4("Show the seq logo."),
                                  conditionalPanel(condition = ("input.act_step28_meta_seq_draw_plot > 0"),
                                                   plotOutput("out_step28_meta_seq_logo") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      ############################################################################################################
      ## step29 SeRP enrichment ###############################
      ### step29 draw the selective Ribo-seq enrichment ##########################################################
      tabItem(tabName = "step29_serp_enrich",
              h2("Draw the SeRP enrichment"),
              fluidRow(
                # tags$style(HTML(".checkbox label {display: flex; align-items: center; margin-top: 25px;}")),
                
                # set the input panel
                column(5,
                       ## click to import the metalot table
                       box(title = "Import the meta gene", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step29_serp_enrich_meta", label = "meta gene", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "meta gene table"),
                           
                           actionButton("act_step29_import_serp_enrich_meta", label = "import the metagene", icon = icon("file")),
                           
                           tags$style(HTML("#act_step29_import_serp_enrich_meta {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Calculate the enrichment", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(4, textInput("in_step29_serp_enrich_group", label = "Group", value = "Sample")),
                             column(4, selectizeInput("in_step29_serp_enrich_method", label = "Method",
                                                      choices = c("mix", 'individual'), selected = 'individual'))
                           ),
                           p("The 'mix' method only suitable for the mixed samples and will calculate the enrichment for the whole group,
                           The 'individual' method will calculate the enrichment for paired samples."),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step29_serp_enrich_input", label = "Input", 
                                                      choices = NULL, selected = NULL, multiple = TRUE)),
                             column(4, selectizeInput("in_step29_serp_enrich_ip", label = "Flag-IP", 
                                                      choices = NULL, selected = NULL, multiple = TRUE)),
                             column(4, selectizeInput("in_step29_serp_enrich_mock", label = "Mock-IP", 
                                                      choices = NULL, selected = NULL, multiple = TRUE))
                           ),
                           
                           actionButton("act_step29_calc_serp_enrich", label = "calculate the enrichment", icon = icon("calculator")),
                           
                           tags$style(HTML("#act_step29_calc_serp_enrich {background-color: #e78c45; color: black;}"))
                       ),
                       
                       box(title = "Draw the enrichment", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fluidRow(
                             column(6, selectizeInput("in_step29_serp_enrich_unit", label = "Unit", choices = c("aa", "nt"), selected = "aa")),
                             column(6, checkboxInput("in_step29_serp_enrich_facet", label = "Facet", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step29_serp_enrich_tis_min", label = "TIS min", value = -10, step = 1)),
                             column(6, numericInput("in_step29_serp_enrich_tis_max", label = "TIS max", value = 100, step = 1))
                           ),
                           
                           fluidRow(
                             column(6, numericInput("in_step29_serp_enrich_tts_min", label = "TTS min", value = -100, step = 1)),
                             column(6, numericInput("in_step29_serp_enrich_tts_max", label = "TTS max", value = 10, step = 1))
                           ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step29_serp_enrich_color", label = "Color", selected = "Paired", choices = color_list)),
                             column(6, sliderInput("in_step29_serp_enrich_alpha", label = "Alpha", value = 0.2, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step29_serp_enrich_font_size", label = "Font size", value = 12, min = 0, max = 20, step = 1)),
                             column(6, sliderInput("in_step29_serp_enrich_width", label = "Line width", value = 0.8, min = 0, max = 3, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step29_serp_enrich_out_name", label = "Output", value = "codon", placeholder = "output prefix")),
                             column(4, numericInput("out_step29_serp_enrich_fig_width", label = "Figure width", value = 10, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step29_serp_enrich_fig_height", label = "Figure height", value = 6, min = 1, max = 20, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step29_serp_enrich_draw_plot", label = "draw enrich plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step29_serp_enrich_meta", label = "Save meta plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step29_serp_enrich_draw_plot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step29_serp_enrich_meta {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Meta Gene",
                                  h4("Show the meta gene"),
                                  DT::dataTableOutput("out_step29_serp_enrich_meta")),
                         
                         tabPanel(title = "Enrichment",
                                  h4("Show the meta gene enrichment"),
                                  DT::dataTableOutput("out_step29_serp_enrich_meta_enrich")),
                         
                         tabPanel(title = "Meta plot",
                                  h4("Show the meta plot."),
                                  conditionalPanel(condition = ("input.act_step29_serp_enrich_draw_plot > 0"),
                                                   plotOutput("out_step29_serp_enrich_metaplot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      
      ############################################################################################################
      ## step30 SeRP peaks ###############################
      ### step30 draw the selective Ribo-seq binding peaks #######################################################
      tabItem(tabName = "step30_serp_peaks",
              h2("Draw the SeRP binding peaks"),
              fluidRow(
                # set the input panel
                column(5,
                       ## click to import the metaplot table
                       box(title = "Import the meta gene", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step30_serp_enrich_peaks", label = "enrich peaks", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "enrich peaks table"),
                           
                           actionButton("act_step30_import_serp_enrich_peaks", label = "import the enrich peaks", icon = icon("file")),
                           
                           tags$style(HTML("#act_step30_import_serp_enrich_peaks {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Draw the SeRP peaks", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           # fluidRow(
                           #   column(6, selectizeInput("in_step30_serp_enrich_unit", label = "Unit", choices = c("aa", "nt"), selected = "aa")),
                           #   column(6, checkboxInput("in_step30_serp_enrich_facet", label = "Facet", value = FALSE))
                           # ),
                           
                           fluidRow(
                             column(4, selectizeInput("in_step30_serp_enrich_gene", label = "Gene name", selected = "", 
                                                      choices = NULL, multiple = FALSE))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step30_serp_enrich_x_min", label = "RNA X min", value = 0, min = 1, step = 1)),
                             column(4, numericInput("in_step30_serp_enrich_x_max", label = "RNA X max", value = NA, min = 1, step = 1)),
                             column(4, numericInput("in_step30_serp_enrich_x_break", label = "RNA X breaks", value = 5, min = 3, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step30_serp_enrich_y_min", label = "RNA Y min", value = 0, min = 0, step = 1)),
                             column(4, numericInput("in_step30_serp_enrich_y_max", label = "RNA Y max", value = NA, min = 0, step = 1)),
                             column(4, numericInput("in_step30_serp_enrich_y_break", label = "RNA Y breaks", value = 2, min = 2, step = 1))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step30_serp_enrich_grid_width", label = "Grid width", value = 0.2, min = 0, max = 1, step = 0.1)),
                             column(4, textInput("in_step30_serp_enrich_grid_color", label = "Grid color", value = "grey85")),
                             column(4, checkboxInput("in_step30_serp_enrich_sqrt", label = "Sqrt", value = FALSE))
                           ),
                           
                           fluidRow(
                             column(4, textInput("in_step30_serp_enrich_enrich_color1", label = "Line color1", value = "#0072bd")),
                             column(4, textInput("in_step30_serp_enrich_enrich_color2", label = "Line color2", value = "#edb120")),
                             column(4, textInput("in_step30_serp_enrich_enrich_color3", label = "Line color3", value = "#d95319"))
                           ),
                           
                           fluidRow(
                             column(4, sliderInput("in_step30_serp_enrich_font_size", label = "Font size", value = 12, min = 0, max = 20, step = 1)),
                             column(4, sliderInput("in_step30_serp_enrich_line_width", label = "Line width", value = 0.8, min = 0, max = 3, step = 0.1)),
                             column(4, sliderInput("in_step30_serp_enrich_alpha", label = "Alpha", value = 0.8, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(6, textInput("in_step30_serp_enrich_xlabel", label = "Xlabel", value = "Gene position [AA]")),
                             column(6, textInput("in_step30_serp_enrich_ylabel", label = "Ylabel", value = "Enrichment [a.u.]"))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step30_serp_enrich_out_name", label = "Output", value = "codon", placeholder = "output prefix")),
                             column(4, numericInput("out_step30_serp_enrich_fig_width", label = "Figure width", value = 10, min = 1, max = 20, step = 0.5)),
                             column(4, numericInput("out_step30_serp_enrich_fig_height", label = "Figure height", value = 6, min = 1, max = 20, step = 0.5))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step30_serp_enrich_draw_plot", label = "draw peaks plot", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step30_serp_enrich_peaks", label = "Save peaks plot", 
                                                      class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step30_serp_enrich_draw_plot {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step30_serp_enrich_peaks {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Enrichment",
                                  h4("Show the gene enrichment peaks."),
                                  DT::dataTableOutput("out_step30_serp_enrich_peaks")),
                         
                         tabPanel(title = "Enrich peaks plot",
                                  h4("Show the enrich peaks plot."),
                                  conditionalPanel(condition = ("input.act_step30_serp_enrich_draw_plot > 0"),
                                                   plotOutput("out_step30_serp_enrich_peaks_plot") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      ),
      
      
      ############################################################################################################
      ## step31 SeRP motif ###############################
      ### draw the SeRP motif ####################################################################################
      tabItem(tabName = "step31_serp_motif",
              h2("Draw the SeRP seqlogo"),
              fluidRow(
                # set the input panel
                column(5,
                       ## click to import the metalot table
                       box(title = "Import the sequence data", status = "primary", collapsible = TRUE, collapsed = FALSE, 
                           solidHeader = TRUE, width = NULL,
                           
                           fileInput(inputId = "in_step31_serp_sequence", label = "peaks sequence", accept = c(".txt", ".TXT"),
                                     multiple = TRUE, placeholder = "enrich peaks sequence"),
                           
                           fluidRow(
                             column(4, numericInput("in_step31_serp_upstream_length", label = "Upstream length", value = 10, min = 1, max = 30, step = 1)),
                             column(4, numericInput("in_step31_serp_peak_length", label = "Peak length", value = 10, min = 1, max = 30, step = 1)),
                             column(4, numericInput("in_step31_serp_downstream_length", label = "Downstream length", value = 10, min = 1, max = 30, step = 1))
                           ),
                           
                           actionButton("act_step31_import_serp_sequence", label = "import the peaks sequence", icon = icon("file")),
                           
                           tags$style(HTML("#act_step31_import_serp_sequence {background-color: #75aadb; color: black;}"))
                       ),
                       
                       box(title = "Calc ulatethe consensus matrix.", collapsed = FALSE, collapsible = TRUE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           
                           fluidRow(
                             column(6, selectizeInput("in_step31_serp_method", label = "Method", 
                                                      choices = c('bits', 'probability', 'custom'), selected = "probability")),
                             conditionalPanel(condition = "input.in_step31_serp_method == 'custom'",
                                              column(4, selectizeInput("in_step31_serp_custom", label = "Custom", selected = "foldchange", 
                                                                       choices = c("delta", "foldchange", "enrichment"))))
                           ),
                           
                           actionButton("act_step31_calc_serp_consensus_matrix", label = "calc the consensus matrix", icon = icon("calculator")),
                           
                           tags$style(HTML("#act_step31_calc_serp_consensus_matrix {background-color: #75aadb; color: black;}"))
                           
                       ),
                       
                       box(title = "Draw the SeRP seqlogo.", collapsed = FALSE, collapsible = TRUE, 
                           width = NULL, status = "primary", solidHeader = TRUE,
                           
                           fluidRow(
                             column(4, selectizeInput("in_step31_serp_region", label = "Region", selected = "upstream",
                                                      choices = c('upstream', 'downstream'))),
                             column(4, selectizeInput("in_step31_serp_seqtype", label = "SeqType", selected = "aa",
                                                      choices = c('rna', 'aa'))),
                             column(4, selectizeInput("in_step31_serp_font_family", label = "Font family", selected = "helvetica_regular", 
                                                      choices = seqlogo_font_family))
                           ),
                           
                           fluidRow(
                             column(4, numericInput("in_step31_serp_seqlogo_x_min", label = "X min", value = NA, min = -30, max = 30, step = 1)),
                             column(4, numericInput("in_step31_serp_seqlogo_x_max", label = "X max", value = NA, min = -30, max = 30, step = 1)),
                             column(4, numericInput("in_step31_serp_seqlogo_x_breaks", label = "Breaks", value = 1, min = 1, max = 5, step = 1))
                           ),
                           
                           # fluidRow(
                           #   column(6, numericInput("in_step31_serp_seqlogo_y_min", label = "Y min", value = -10, min = -30, max = 30, step = 1)),
                           #   column(6, numericInput("in_step31_serp_seqlogo_y_max", label = "Y max", value = 10, min = -30, max = 30, step = 1))
                           # ),
                           
                           fluidRow(
                             column(6, selectizeInput("in_step31_serp_fill_color", label = "Fill color", selected = "auto", choices = seqlogo_colors)),
                             column(6, sliderInput("in_step31_serp_fill_alpha", label = "Fill alpha", value = 0.9, min = 0, max = 1, step = 0.05))
                           ),
                           
                           fluidRow(
                             column(6, sliderInput("in_step31_serp_font_size", label = "Font size", value = 15, min = 5, max = 25, step = 0.5)),
                             column(6, sliderInput("in_step31_serp_font_stack", label = "Stack width", value = 0.9, min = 0, max = 1, step = 0.1))
                           ),
                           
                           fluidRow(
                             column(4, textInput("out_step31_serp_fig_name", label = "Output", value = "SeRP-peaks", placeholder = "file name prefix")),
                             column(4, numericInput("out_step31_serp_fig_width", label = "Figure width", value = 10, min = 1)),
                             column(4, numericInput("out_step31_serp_fig_height", label = "Figure height", value = 5, min = 1))
                           ),
                           
                           fluidRow(
                             column(6, actionButton("act_step31_draw_serp_seqlogo", label = "draw the seqlogo", icon = icon("chart-line"))),
                             column(6, downloadButton("save_step31_serp_seqlogo", label = "Save seqlogo", class = "download-btn", icon = icon("download")))
                           ),
                           
                           tags$style(HTML("#act_step31_draw_serp_seqlogo {background-color: #e78c45; color: black;}")),
                           tags$style(HTML("#save_step31_serp_seqlogo {background-color: #grey; color: black; margin-top: 0px;}"))
                           
                       )
                ),
                
                column(7,
                       tabsetPanel(
                         ## set the output panel
                         tabPanel(title = "Sequence",
                                  h4("Show the sequence of gene enrichment peaks."),
                                  DT::dataTableOutput("out_step31_serp_peaks_sequence")),
                         
                         tabPanel(title = "Upstream_nt",
                                  h4("Show the consensus matrix of Upstream_nt."),
                                  DT::dataTableOutput("out_step31_serp_peaks_upstream_nt")),
                         tabPanel(title = "Downstream_nt",
                                  h4("Show the consensus matrix of Downstream_nt."),
                                  DT::dataTableOutput("out_step31_serp_peaks_downstream_nt")),
                         tabPanel(title = "Upstream_aa",
                                  h4("Show the consensus matrix of Upstream_aa."),
                                  DT::dataTableOutput("out_step31_serp_peaks_upstream_aa")),
                         tabPanel(title = "Downstream_aa",
                                  h4("Show the consensus matrix of Downstream_aa."),
                                  DT::dataTableOutput("out_step31_serp_peaks_downstream_aa")),
                         
                         tabPanel(title = "Seqlogo",
                                  h4("Show the seqlogo of enrich peaks."),
                                  conditionalPanel(condition = ("input.act_step31_draw_serp_seqlogo > 0"),
                                                   plotOutput("out_step31_serp_peaks_seqlogo") %>% 
                                                     withSpinner(color="#3c8dbc", type = 5, size = 1)))
                       )
                )
              )
      )
      
      
      
      ############################################################################################################
      ## step32 smORF ###############################
      ############################################################################################################
      
      
      
      ############################################################################################################
      ## step33 smORF expression ###############################
      ############################################################################################################
      
      
      
      ############################################################################################################
      ## step34 smORF mORF association ###############################
      ############################################################################################################
      
      
      
      
      
    )
  )
)
