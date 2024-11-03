#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#


function(input, output, session) {
  ## set the input file size
  options(shiny.maxRequestSize = 500*1024^2)
  
  ############################################################
  # 0 figure of overview #################
  ############################################################
  
  ## 0.1 figure of frame work #################
  fig_framework <- renderImage({
    filename <- normalizePath(file.path('images/', 'fig1_flowchart_v4_compress.png'))
    list(src = filename,
         contentType = 'png',
         width = NULL,
         height = 600,
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  ## 0.2 figure of offset detection #################
  fig_offset_detection <- renderImage({
    filename <- normalizePath(file.path('images/', 'fig2_offset_detection_v2_compress.png'))
    list(src = filename,
         contentType = 'png',
         width = NULL,
         height = 600,
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  ## 0.3 figure of quality control #################
  fig_quality_control <- renderImage({
    filename <- normalizePath(file.path('images/', 'fig3_quality_control_v1_compress.png'))
    list(src = filename,
         contentType = 'png',
         width = NULL,
         height = 600,
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  ## 0.4 figure of codon level #################
  fig_codon_level <- renderImage({
    filename <- normalizePath(file.path('images/', 'fig4_codon_level_v1_compress.png'))
    list(src = filename,
         contentType = 'png',
         width = NULL,
         height = 600,
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  ## 0.5 figure of gene level #################
  fig_gene_level <- renderImage({
    filename <- normalizePath(file.path('images/', 'fig5_gene_level_v1_compress.png'))
    list(src = filename,
         contentType = 'png',
         width = NULL,
         height = 600,
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  ## 0.4 figure of serp #################
  fig_serp <- renderImage({
    filename <- normalizePath(file.path('images/', 'fig6_serp_v1_compress.png'))
    list(src = filename,
         contentType = 'png',
         width = NULL,
         height = 600,
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  ## 0.5 figure of smorf-filter-model.png #################
  fig_smorf <- renderImage({
    filename <- normalizePath(file.path('images/', 'fig7_smorf-filter-model_compress.png'))
    list(src = filename,
         contentType = 'png',
         width = NULL,
         height = 600,
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  output$Plot_fig_framework <- fig_framework
  output$Plot_fig_offset_detection <- fig_offset_detection
  output$Plot_fig_quality_control <- fig_quality_control
  output$Plot_fig_codon_level <- fig_codon_level
  output$Plot_fig_gene_level <- fig_gene_level
  output$Plot_fig_serp <- fig_serp
  output$Plot_fig_smorf <- fig_smorf
  
  ############################################################
  # # get the client IP address #################
  # 
  # client_ip <- reactive({
  #   # browser()
  #   
  #   ip <- tryCatch({
  #     session$request$HTTP_HOST
  #     session$clientData$url_hostname
  #     
  #     # paste0(session$clientData$url_protocol, 
  #     #        session$clientData$url_pathname, session$clientData$url_pathname, 
  #     #        session$clientData$url_hostname, ':', session$clientData$url_port)
  #     
  #   }, error = function(e) {
  #     "Unknown IP"
  #   })
  #   ip
  # })
  # 
  # # show the IP address
  # output$ip_address <- renderText({
  #   paste("Your IP address is:", client_ip())
  # })
  # 
  

  ############################################################
  # 1. figure of detabase #################
  
  ## 1.1 show the table of example design #################
  ## import the example design
  step1_example <- reactive({
    req(input$act_step1_show_example)
    design_df <- read.xlsx("data/example_design.xlsx", sheet = 1, colNames = TRUE, rowNames = F) %>% 
      dplyr::mutate_all(as.character)
    
    return(design_df)
  })
  
  observeEvent(input$act_step1_show_example, {
    ## output the table
    output$out_step1_example <- DT::renderDataTable(server = F, {
      DT::datatable(step1_example(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10, 
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "example_table"), 
                                     list(extend = 'excel', filename = "example_table")))) %>% 
        formatStyle(names(step1_example()), textAlign = "right") 
    })
  })
  
  ## 1.2 show the table of design #################
  ## import the design
  import_design_clicked <- reactiveVal(FALSE)
  
  step1_design <- reactive({
    req(input$act_step1_import_design)
    import_design_clicked(TRUE)
    
    in_step1_design <- input$in_step1_design
    if (is.null(in_step1_design)) {return(NULL)}
    
    # import the sort the data with rank
    design_df <- read.xlsx(input$in_step1_design$datapath, 
                           sheet = input$in_step1_design_sheet, colNames = TRUE, rowNames = F)
    
    # check the column names
    # browser()
    
    # contains these columns: "Sample", "SeqType", "Rank", "Group"
    if (!all(c("Sample", "SeqType", "Rank", "Group") %in% colnames(design_df))) {
      # set the warning
      warning("The design file should contains four columns: 'Sample', 'SeqType', 'Rank', 'Group'")
      # command  to clear the input
      updatefileInput(session, "in_step1_design", label = "Upload the design file")
      
      return(NULL)
    }
    
    design_df <- design_df %>% 
      dplyr::arrange(SeqType, Rank) %>%
      dplyr::mutate(Sample = factor(Sample, levels = unique(Sample)),
                    Group = factor(Group, levels = unique(Group))) %>% 
      dplyr::mutate_all(as.character)
    
    return(design_df)
  })
  
  ## output the table
  observeEvent(input$act_step1_import_design, {
    output$out_step1_design <- DT::renderDataTable(server = F, { 
      DT::datatable(step1_design(), rownames = FALSE, filter = "top", extensions = 'Buttons', 
                    options = list(lengthMenu = c(10, 50, 100, 500),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtip',
                                   buttons = list(
                                     list(extend = 'csv', filename = "design_table"), 
                                     list(extend = 'excel', filename = "design_table")))) %>% 
        formatStyle(names(step1_design()), textAlign = "right") 
    })
  })
  
  ## 1.3 check and create the output directory #################
  # function for list files
  # list_files_tree <- function(path = ".", depth = 0) {
  #   files <- list.files(path, full.names = TRUE)
  #   result <- c(paste(rep("  ", depth), basename(path), "\n", sep = ""))
  #   for (file in files) {
  #     if (file.info(file)$isdir) {
  #       result <- c(result, paste0(rep("  ", depth + 1), "|--", basename(file), "/\n"))
  #     } else {
  #       result <- c(result, paste0(rep("  ", depth + 1), "|--", basename(file), "\n"))
  #     }
  #   }
  #   return(result)
  # }
  
  ## function for list files
  list_files_tree <- function(path = ".", depth = 0) {
    files <- list.files(path, full.names = TRUE)
    result <- data.frame(folder_structure = basename(path))
    for (file in files) {
      if (file.info(file)$isdir) {
        result <- rbind(result, paste0(rep("  ", depth + 1), "|--", basename(file), "/\n"))
      } else {
        result <- rbind(result, paste0(rep("  ", depth + 1), "|--", basename(file), "\n"))
      }
    }
    return(result)
  }
  
  ## check the output directory
  step1_output_dir <- reactive({
    req(input$act_create_folder)
    
    output_folder <- input$output_folder
    if (nchar(output_folder) == 0) {return(NULL)}
    
    if (!dir.exists(output_folder)) {
      dir.create(output_folder)
    }
    
    return(output_folder)
  })
  
  ## output the file list
  observeEvent(input$act_create_folder, {
    output$out_step1_folder <- renderTable({
      # browser()
      list_files_tree(step1_output_dir())
      # list.files(step1_output_dir())
    })
  })
  
  ## 1.4 show the annotation file #################
  ## import the annotation
  import_anno_clicked <- reactiveVal(FALSE)
  step1_annotation <- reactive({
    req(input$act_step1_import_anno)
    import_anno_clicked(TRUE)
    
    in_step1_anno <- input$in_step1_anno
    if (is.null(in_step1_anno)) {return(NULL)}
    
    annotation_file <- input$in_step1_anno$datapath
    annotation_df <- read.xlsx(annotation_file, sheet = input$in_step1_anno_sheet, colNames = TRUE, rowNames = F)
    
    if (input$in_step1_anno_longest){
      annotation_df <- annotation_df %>% 
        filter(rep_transcript == TRUE)
    }
    
    return(annotation_df)
  })
  
  ## output the table
  observeEvent(input$act_step1_import_anno, {
    # browser()
    ## output the table
    output$out_step1_annotation <- DT::renderDataTable(server = T, {
      DT::datatable(step1_annotation(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6, 
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "gene_table"), 
                                     list(extend = 'excel', filename = "gene_table"))))
    })
    
  })
  
  ## draw the length distribution
  step1_hist_plot <- reactive({
    req(input$act_step1_draw_histplot)
    
    # options(warn = -1)
    source("R/draw_cds_len_distr.R")
    
    if (input$in_step1_anno_type == 'cds_length'){
      gene_labels <- 'cds_length'
      plot_title = "CDS length distribution"
    } else if (input$in_step1_anno_type == 'utr5_length'){
      gene_labels <- 'utr5_length'
      plot_title = "5'-UTR length distribution"
    } else if (input$in_step1_anno_type == 'utr3_length'){
      gene_labels <- 'utr3_length'
      plot_title = "3'-UTR length distribution"
    } 
    
    hist_plot <- draw_cds_len_distr(distr = step1_annotation(), 
                                    label = gene_labels,
                                    title = plot_title,
                                    edge_color = input$in_step1_edge_color,
                                    fill_color= input$in_step1_fill_color,
                                    fill_alpha = input$in_step1_fill_alpha,
                                    font_size = input$in_step1_font_size,
                                    binwidth = input$in_step1_bin_range, 
                                    xmin = input$in_step1_length_xlim[1],
                                    xmax = input$in_step1_length_xlim[2],
                                    sqrty = input$in_step1_sqrty)
    
    return(list(hist_plot = hist_plot, gene_labels = gene_labels))
  })
  
  observeEvent(input$act_step1_draw_histplot, {
    # browser()
    
    ## output the length distribution
    output$out_step1_hist_plot <- renderPlot(
      width = input$out_step1_hist_width * 100,
      height = input$out_step1_hist_height * 100,
      {step1_hist_plot()$hist_plot})
    
    ## download the figure
    output$save_step1_hist <- downloadHandler(
      filename = function() {
        paste0(step1_hist_plot()$gene_labels, "-Histogram-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step1_hist_plot(), filename = file, width = input$out_step1_hist_width, height = input$out_step1_hist_height)
      }
    )
    
    # options(warn = 0)
  })
  
  ## 1.5 create the orgdb file #################
  ### 1.5.0 create the orgdb file #################
  ### retrieve the annotationhub file
  step1_hub <- reactive({
    req(input$act_step1_retrieve_orgdb)
    
    hub <- AnnotationHub::AnnotationHub(localHub = input$in_step1_orgdb_local)
    
    hub_table <- data.frame(species = hub$species,
                            taxonomy = hub$taxonomyid,
                            genome = hub$genome,
                            ah_id = hub$ah_id,
                            title = hub$title,
                            description = hub$description,
                            class = hub$rdataclass,
                            date = hub$rdatadateadded)
    
    return(list(hub = hub, hub_table = hub_table))
  })
  
  observeEvent(input$act_step1_retrieve_orgdb, {
    ## output the table
    output$out_step1_annotationhub <- DT::renderDataTable(server = T, {
      # browser()
      
      DT::datatable(step1_hub()$hub_table, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6, 
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "annotationhub_table"), 
                                     list(extend = 'excel', filename = "annotationhub_table"))))
      
    })
    
  })
  
  ### 1.5.1 install the orgdb #############################
  ## retrieve the orgdb file
  step1_retrieve_orgdb <- reactive({
    req(input$act_step1_filter_orgdb)
    
    if (is.null(step1_hub())) {return(NULL)}
    
    # browser()
    
    if (nchar(input$in_step1_orgdb_genus) != 0 & nchar(input$in_step1_orgdb_species) != 0) {
      my_hub <- AnnotationHub::query(step1_hub()$hub, 
                                     paste(input$in_step1_orgdb_genus, input$in_step1_orgdb_species))
      
      hub_table <- data.frame(species = my_hub$species,
                              taxonomy = my_hub$taxonomyid,
                              genome = my_hub$genome,
                              ah_id = my_hub$ah_id,
                              title = my_hub$title,
                              description = my_hub$description,
                              class = my_hub$rdataclass,
                              date = my_hub$rdatadateadded)
      
    } else if(ncahr(input$in_step1_orgdb_taxonomy) != 0) {
      my_hub <- AnnotationHub::query(step1_hub()$hub, input$in_step1_orgdb_taxonomy)
      
      hub_table <- data.frame(species = my_hub$species,
                              taxonomy = my_hub$taxonomyid,
                              genome = my_hub$genome,
                              ah_id = my_hub$ah_id,
                              title = my_hub$title,
                              description = my_hub$description,
                              class = my_hub$rdataclass,
                              date = my_hub$rdatadateadded)
      
    } else {
      return(NULL)
    }
    
    return(list(my_hub = my_hub, hub_table = hub_table))
  })
  
  ## retrieve the orgdb file
  observeEvent(input$act_step1_filter_orgdb, {
    # browser()
    ## output the table
    output$out_step1_orgdb_preview <- DT::renderDataTable(server = T, {
      DT::datatable(step1_retrieve_orgdb()$hub_table, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6, 
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "orgdb_table"), 
                                     list(extend = 'excel', filename = "orgdb_table"))))
    })
    
  })
  
  ### 1.5.2 install the orgdb #############################
  ## create the orgdb file
  step1_create_orgdb <- reactive({
    req(input$act_step1_create_orgdb)
    
    if (is.null(step1_retrieve_orgdb())) {return(NULL)}
    
    # browser()
    
    # retrieve the orgdb file
    my_orgdb <- step1_retrieve_orgdb()$my_hub[[input$in_step1_orgdb_records]]
    
    # make the gson go
    gson_go <- gson_GO(OrgDb = my_orgdb, keytype = 'ENTREZID', ont = 'ALL')
    
    # retrieve the gene information
    orgdb_column <- gsub(' ', '', input$in_step1_orgdb_column)
    orgdb_column <- gsub('\n', ',', orgdb_column)
    orgdb_column <- unlist(str_split(orgdb_column, ','))
    
    gene_mess <- AnnotationDbi::select(my_orgdb,
                                       keys = keys(my_orgdb, keytype = 'ENTREZID'),
                                       columns = orgdb_column)
    
    # browser()
    gene_column <- columns(my_orgdb)
    
    gene_table <- gene_mess %>% 
      dplyr::select(any_of(c('GID', 'SYMBOL', 'GENENAME'))) %>% 
      # tidyr::drop_na(GID) %>% 
      na.omit() %>% 
      unique()
    
    go_table <- gene_mess %>% 
      dplyr::select(any_of(c('GID', 'GO', 'EVIDENCE'))) %>% 
      # tidyr::drop_na(GID) %>% 
      na.omit() %>% 
      unique()
    
    
    # browser()
    
    org_genus <- gsub(' ', '', input$in_step1_orgdb_genus)
    org_species <- gsub(' ', '', input$in_step1_orgdb_species)
    
    # check the directory exists and create it if not
    if (!dir.exists(input$out_step1_orgdb_dir)) {
      dir.create(input$out_step1_orgdb_dir, recursive = TRUE)
    }
    
    orgdb_file = makeOrgPackage(version = "0.1",
                                gene_info = gene_table,
                                go = go_table,
                                maintainer = "rensc <rensc@163.com>",
                                author = "rensc <rensc@163.com>",
                                outputDir = input$out_step1_orgdb_dir,
                                tax_id = input$in_step1_orgdb_taxonomy,
                                genus = org_genus,
                                species = org_species,
                                goTable = 'go', 
                                verbose = TRUE)
    
    devtools::build(pkg = paste0(input$out_step1_orgdb_dir,
                                 "/",
                                 "org.", 
                                 substr(org_genus[1], 1, 1), 
                                 org_species,
                                 ".eg.db"))
    
    return(list(orgdb = my_orgdb, gene_mess = gene_mess, gson_go = gson_go))
  })
  
  observeEvent(input$act_step1_create_orgdb, {
    # browser()
    ## output the table
    output$out_step1_orgdb_mess <- DT::renderDataTable(server = T, {
      DT::datatable(step1_create_orgdb()$gene_mess, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6, 
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "orgdb_gene_table"), 
                                     list(extend = 'excel', filename = "orgdb_gene_table"))))
    })
    
    ## output the gson go
    output$save_step1_gson_go <- downloadHandler(
      filename = function() {
        paste0(input$out_step1_orgdb_name, "_go_", Sys.Date(), ".gson")
      },
      content = function(file) {
        gson::write.gson(step1_create_orgdb()$gson_go, file = file)
      }
    )
    
    ## output the go message
    output$save_step1_gene_mess <- downloadHandler(
      filename = function() {
        paste0(input$out_step1_orgdb_name, "_gene_mess_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(step1_create_orgdb()$gene_mess, file = file, RowName = FALSE)
      }
    )
    
    ## output the orgdb
    output$save_step1_orgdb <- downloadHandler(
      filename = function() {
        paste0(input$out_step1_orgdb_name, "_orgdb_", Sys.Date(), ".RData")
      },
      content = function(file) {
        save.image(step1_create_orgdb()$gson_go, file = file)
      }
    )
    
  })
  
  ### 1.5.3 install the orgdb #############################
  ## install the orgdb
  step1_install_orgdb <- reactive({
    req(input$act_step1_install_orgdb)
    
    # browser()
    if (is.null(input$in_step1_orgdb_file)) {return(NULL)}
    
    orgdb_file <- input$in_step1_orgdb_file$datapath
    
    # install the orgdb
    install.packages(orgdb_file, repos = NULL, type = "source")
    
    orgdb_log <- "Install the orgdb successfully!"
    
    return(orgdb_log)
  })
  
  observeEvent(input$act_step1_install_orgdb, {
    # browser()
    
    ## output the table
    output$out_step1_orgdb_install <- renderText({
      step1_install_orgdb()
    })
    
  })
  
  ### 1.6 retrieve the kegg file #############################
  ### retrieve the annotationhub file
  step1_gson_kegg <- reactive({
    req(input$act_step1_retrieve_kegg)
    
    if (is.null(input$in_step1_kegg_species)) {return(NULL)}
    
    # browser()
    
    my_kegg <- gson_KEGG(input$in_step1_kegg_species)
    
    kegg_table <- my_kegg@gsid2gene %>% 
      dplyr::left_join(my_kegg@gsid2name, by = 'gsid') %>% 
      dplyr::mutate(version = my_kegg@version,
                    date = my_kegg@accessed_date,
                    database = my_kegg@gsname,
                    species = my_kegg@species) %>% 
      dplyr::select(species, gsid, gene, name, version, date, database)
    
    return(list(kegg = my_kegg, kegg_table = kegg_table))
  })
  
  observeEvent(input$act_step1_retrieve_kegg, {
    ## output the table
    output$out_step1_kegg_mess <- DT::renderDataTable(server = T, {
      # browser()
      DT::datatable(step1_gson_kegg()$kegg_table, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6, 
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "kegg_table"), 
                                     list(extend = 'excel', filename = "kegg_table"))))
      
    })
    
    # output the gson kegg
    output$save_step1_gson_kegg <- downloadHandler(
      filename = function() {
        paste0(input$out_step1_kegg_name, "_KEGG_", Sys.Date(), ".gson")
      },
      content = function(file) {
        gson::write.gson(step1_gson_kegg()$kegg, file = file)
      }
    )
    
  })
  
  
  ## 1.7 show the sequence file #################
  ## import the gene sequence
  
  step1_sequence <- reactive({
    req(input$act_step1_import_seq)
    
    if (is.null(input$in_step1_seq)) {return(NULL)}
    
    # import the sequence
    cds_file <- input$in_step1_seq$datapath
    cds <- readSet(file = cds_file)
    
    # readDNAStringSet(input$in_step1_seq$datapath, format="fasta",
    #                  nrec=-1L, skip=0L, seek.first.rec=FALSE,
    #                  use.names=TRUE, with.qualities=FALSE)
    
    # browser()
    cds_preview = c()
    for (i in c(1:6)) {
      cds_preview = c(cds_preview, paste0(">", names(cds[i]), "\n", toString(cds[i]), "\n"))
    }
    
    return(list(cds = cds, cds_preview = cds_preview))
    
  })
  
  ## output the sequence
  observeEvent(input$act_step1_import_seq, {
    
    ## output the sequence
    output$out_step1_sequence <- renderText({
      step1_sequence()$cds_preview
    })
    
  })
  
  
  ## 1.8 calculate the codon usage #################
  ## import the gene sequence
  
  step1_codon_usage <- reactive({
    req(input$act_step1_calc_cu)
    
    if (is.null(step1_sequence())) {return(NULL)}
    
    # browser()
    
    source("R/create_codon_usage_table.R")
    codon_table <- create_codon_usage_table(seq = "DNA")
    
    # calculate the codon count
    cds_table <- codonTable(step1_sequence()$cds)
    
    cds_data_frame <- data.frame(Gene = cds_table@ID, Length = cds_table@len, cds_table@counts)
    
    cds_data_frame_long <- cds_data_frame %>% 
      tidyr::gather(key = 'Codon', value = 'Count', AAA:TTT) %>% 
      left_join(codon_table %>% dplyr::select(Codon, AA), by = "Codon") %>% 
      group_by(Gene, Codon) %>%
      dplyr::mutate(Frequency = Count / Length * 1000) %>% 
      group_by(Gene, AA) %>%
      dplyr::mutate(RSCU = Frequency / mean(Frequency),
                    CAI = Frequency / max(Frequency)) %>% 
      ungroup() %>% 
      tidyr::replace_na(list(Frequency = 0, RSCU = 0, CAI = 0)) %>%
      dplyr::mutate(Frequency = round(Frequency, input$in_step1_float),
                    RSCU = round(RSCU, input$in_step1_float),
                    CAI = round(CAI, input$in_step1_float))
    
    # calculate the codon usage of each gene
    cds_freq_rscu_cai <- list(
      Count = cds_data_frame,
      Frequency = cds_data_frame_long %>% dplyr::select(Gene, Length, Codon, Frequency) %>% tidyr::pivot_wider(names_from = Codon, values_from = Frequency),
      RSCU = cds_data_frame_long %>% dplyr::select(Gene, Length, Codon, RSCU) %>% tidyr::pivot_wider(names_from = Codon, values_from = RSCU),
      CAI = cds_data_frame_long %>% dplyr::select(Gene, Length, Codon, CAI) %>% tidyr::pivot_wider(names_from = Codon, values_from = CAI)
    )
    
    # calculate the codon usage of whole genes
    cds_codon_usage <- cds_data_frame %>% 
      gather(key = 'Codon', value = 'Count', AAA:TTT) %>% 
      group_by(Codon) %>%
      reframe(Count = sum(Count)) %>%
      mutate(Frequency = Count / sum(Count) * 1000)
    
    cds_codon_usage <- codon_table %>%
      left_join(cds_codon_usage, by = "Codon") %>% 
      group_by(AA) %>%
      arrange(AA, desc(Frequency)) %>%
      mutate(RSCU = Frequency / mean(Frequency),
             CAI = Frequency / max(Frequency)) %>% 
      mutate(Frequency = round(Frequency, input$in_step1_float),
             RSCU = round(RSCU, input$in_step1_float),
             CAI = round(CAI, input$in_step1_float)) %>% 
      ungroup()
    
    rm(cds_data_frame)
    
    return(list(cds_codon_usage = cds_codon_usage,
                cds_freq_rscu_cai = cds_freq_rscu_cai))
  })
  
  ## output the codon usage table and plot
  observeEvent(input$act_step1_calc_cu, {
    # browser()
    
    ## output the table
    output$out_step1_codon_table <- DT::renderDataTable(server = T, {
      DT::datatable(step1_codon_usage()$cds_freq_rscu_cai$Count, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10, 
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "gene_codon_usage"), 
                                     list(extend = 'excel', filename = "gene_codon_usage"))))
    })
    
    ## save the codon count table
    output$save_step1_codon_table <- downloadHandler(
      filename = function() {
        paste0(input$out_step1_codon_count_name, "-gene-codon-usage-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step1_codon_usage()$cds_freq_rscu_cai, file = file, rowNames = FALSE)
      }
    )
    
    ## output the codon usage table
    output$out_step1_cu_table <- DT::renderDataTable(server = F, {
      DT::datatable(step1_codon_usage()$cds_codon_usage, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10, 
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "whole_codon_usage"), 
                                     list(extend = 'excel', filename = "whole_codon_usage"))))
    })
    
  })
  
  
  ## 1.9 draw the codon usage #################
  ## import the gene sequence
  
  step1_cu_plot <- reactive({
    req(input$act_step1_draw_cu)
    
    if (is.null(step1_codon_usage())) {return(NULL)}
    
    # browser()
    source("R/draw_codon_usage.R")
    
    ## draw the figure
    freq_plot <- draw_codon_usage(codon_usage = step1_codon_usage()$cds_codon_usage,
                                  cu_class = "Frequency", stop_codon = input$in_step1_stop,
                                  type = input$in_step1_cu_plot_type,
                                  grid_width = input$in_step1_cu_grid_width, grid_color = input$in_step1_cu_grid_color,
                                  font_size = input$in_step1_cu_font_size, dot_size = input$in_step1_cu_dot_size,
                                  color_alpha = input$in_step1_cu_color_alpha, line_width = input$in_step1_cu_line_width,
                                  color_map = input$in_step1_cu_color_map, legend_row = input$in_step1_cu_legend_row,
                                  wrap = input$in_step1_cu_wrap, wrap_row = input$in_step1_cu_wrap_row)
    
    rscu_plot <- draw_codon_usage(codon_usage = step1_codon_usage()$cds_codon_usage,
                                  cu_class = "RSCU", stop_codon = input$in_step1_stop,
                                  type = input$in_step1_cu_plot_type,
                                  grid_width = input$in_step1_cu_grid_width, grid_color = input$in_step1_cu_grid_color,
                                  font_size = input$in_step1_cu_font_size, dot_size = input$in_step1_cu_dot_size,
                                  color_alpha = input$in_step1_cu_color_alpha, line_width = input$in_step1_cu_line_width,
                                  color_map = input$in_step1_cu_color_map, legend_row = input$in_step1_cu_legend_row,
                                  wrap = input$in_step1_cu_wrap, wrap_row = input$in_step1_cu_wrap_row)
    
    cai_plot <- draw_codon_usage(codon_usage = step1_codon_usage()$cds_codon_usage,
                                 cu_class = "CAI", stop_codon = input$in_step1_stop,
                                 type = input$in_step1_cu_plot_type,
                                 grid_width = input$in_step1_cu_grid_width, grid_color = input$in_step1_cu_grid_color,
                                 font_size = input$in_step1_cu_font_size, dot_size = input$in_step1_cu_dot_size,
                                 color_alpha = input$in_step1_cu_color_alpha, line_width = input$in_step1_cu_line_width,
                                 color_map = input$in_step1_cu_color_map, legend_row = input$in_step1_cu_legend_row,
                                 wrap = input$in_step1_cu_wrap, wrap_row = input$in_step1_cu_wrap_row)
    
    return(list(freq_plot = freq_plot, rscu_plot = rscu_plot, cai_plot = cai_plot))
  })
  
  
  ## output the codon usage table and plot
  observeEvent(input$act_step1_draw_cu, {
    # browser()
    ## output the table
    output$out_step1_freq_plot <- renderPlot(
      width = input$out_step1_cu_fig_width * 100,
      height = input$out_step1_cu_fig_height * 100,
      {step1_cu_plot()$freq_plot})
    
    output$out_step1_rscu_plot <- renderPlot(
      width = input$out_step1_cu_fig_width * 100,
      height = input$out_step1_cu_fig_height * 100,
      {step1_cu_plot()$rscu_plot})
    
    output$out_step1_cai_plot <- renderPlot(
      width = input$out_step1_cu_fig_width * 100,
      height = input$out_step1_cu_fig_height * 100,
      {step1_cu_plot()$cai_plot})
    
    ## download the figure
    output$save_step1_cu_freq <- downloadHandler(
      filename = function() {
        paste0(input$out_step1_cu_fig_name, "-Frequency-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step1_cu_plot()$freq_plot, filename = file, 
               width = input$out_step1_cu_fig_width, height = input$out_step1_cu_fig_height)
      }
    )
    
    output$save_step1_cu_rscu <- downloadHandler(
      filename = function() {
        paste0(input$out_step1_cu_fig_name, "-RSCU-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step1_cu_plot()$rscu_plot, filename = file, 
               width = input$out_step1_cu_fig_width, height = input$out_step1_cu_fig_height)
      }
    )
    
    output$save_step1_cu_cai <- downloadHandler(
      filename = function() {
        paste0(input$out_step1_cu_fig_name, "-CAI-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step1_cu_plot()$cai_plot, filename = file, 
               width = input$out_step1_cu_fig_width, height = input$out_step1_cu_fig_height)
      }
    )
    
  })
  
  ############################################################
  # 2. figure of mapping #################
  
  ## 2.1 import the mapping data #################
  ## import the mapping data
  step2_align_table <- reactive({
    req(input$act_step2_import_align)
    
    # import the alignment table
    if (is.null(input$in_step2_align_file)) {return(NULL)}
    
    align_table <- read.table(input$in_step2_align_file$datapath, sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>% 
      dplyr::mutate(Database = factor(Database, unique(Database)))
    
    # browser()
    # calculate the average of database
    if (input$in_step2_x == 'Sample') {
      align_table <- align_table
      
    } else if (import_design_clicked()) {
      align_table <- align_table %>%
        dplyr::left_join(step1_design(), by = 'Sample') %>%
        na.omit() %>% 
        dplyr::group_by(!!sym(input$in_step2_x), Database) %>% 
        dplyr::reframe(Count = mean(Count), Ratio = mean(Ratio)) %>% 
        # dplyr::rename(Sample = Group) %>% 
        dplyr::mutate(Database = factor(Database, unique(Database))) %>% 
        dplyr::ungroup() %>% 
        as.data.frame()
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    # format the data
    align_table$Label = paste0(round(align_table$Count / 1000000, 1), 'M')
    
    return(align_table)
  })
  
  ## output the table
  observeEvent(input$act_step2_import_align, {
    # browser()
    ## output the table
    output$out_step2_align_table <- DT::renderDataTable(server = F, {
      DT::datatable(step2_align_table(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "alignment_table"), 
                                     list(extend = 'excel', filename = "alignment_table"))))
    })
    
  })
  
  ## 2.2 draw the mapping data #################
  ## import the mapping data
  
  step2_align_plot <- reactive({
    req(input$act_step2_draw_align)
    
    if (is.null(step2_align_table())) {return(NULL)}
    
    # browser()
    source("R/draw_alignment.R")
    
    align_count <- draw_alignment(align = step2_align_table(),
                                  x = input$in_step2_x, y = "Count", 
                                  fill_color = input$in_step2_fill_color, edge_color = input$in_step2_edge_color,
                                  fill_alpha = input$in_step2_fill_alpha, bar_width = input$in_step2_bar_width,
                                  font_size = input$in_step2_font_size, label_size = input$in_step2_label_size,
                                  trans = input$in_step2_coord, text_label = input$in_step2_label)
    
    align_ratio <- draw_alignment(align = step2_align_table(),
                                  x = input$in_step2_x, y = "Ratio", 
                                  fill_color = input$in_step2_fill_color, edge_color = input$in_step2_edge_color,
                                  fill_alpha = input$in_step2_fill_alpha, bar_width = input$in_step2_bar_width,
                                  font_size = input$in_step2_font_size, label_size = input$in_step2_label_size,
                                  trans = input$in_step2_coord, text_label = input$in_step2_label)
    
    return(list(align_count = align_count, align_ratio = align_ratio))
    
  })
  
  ## output the alignment plot
  observeEvent(input$act_step2_draw_align, {
    # browser()
    ## show the alignment plot
    output$out_step2_count_plot <- renderPlot(
      width = input$out_step2_fig_width * 100,
      height = input$out_step2_fig_height * 100,
      {step2_align_plot()$align_count})
    
    output$out_step2_ratio_plot <- renderPlot(
      width = input$out_step2_fig_width * 100,
      height = input$out_step2_fig_height * 100,
      {step2_align_plot()$align_ratio})
    
    # save the alignment plot
    output$save_step2_count_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step2_fig_name, "-count-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step2_align_plot()$align_count, filename = file, 
               width = input$out_step2_fig_width, height = input$out_step2_fig_height)
      }
    )
    
    output$save_step2_ratio_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step2_fig_name, "-ratio-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step2_align_plot()$align_ratio, filename = file, 
               width = input$out_step2_fig_width, height = input$out_step2_fig_height)
      }
    )
    
  })
  
  
  ############################################################
  # 3. figure of length distribution #################
  
  ## 3.1 import the RNA-seq distribution data #################
  ## import the distribution data
  step3_rna_len_table <- reactive({
    req(input$act_step3_import_rna_length)
    
    # import the data
    if (is.null(input$in_step3_rnaseq)) {return(NULL)}
    
    rna_length_table <- read.table(input$in_step3_rnaseq$datapath, 
                                   sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    # browser()
    # calculate the average of group
    if (input$in_step3_rna_x == 'Sample') {
      rna_length_table <- rna_length_table
      
    } else if (import_design_clicked()) {
      rna_length_table <- rna_length_table %>%
        dplyr::left_join(step1_design(), by = 'Sample') %>%
        na.omit() %>% 
        dplyr::group_by(!!sym(input$in_step3_rna_x), Length) %>% 
        dplyr::reframe(Plus_Count = mean(Plus_Count), Minus_Count = mean(Minus_Count), 
                       Plus_Ratio = mean(Plus_Ratio), Minus_Ratio = mean(Minus_Ratio)) %>% 
        # dplyr::rename(Sample = Group) %>% 
        dplyr::ungroup() %>% 
        as.data.frame()
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    return(rna_length_table)
  })
  
  ## output the table
  observeEvent(input$act_step3_import_rna_length, {
    # browser()
    ## output the table
    output$out_step3_rna_distr <- DT::renderDataTable(server = F, {
      DT::datatable(step3_rna_len_table(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "rna_length_table"), 
                                     list(extend = 'excel', filename = "rnaa_length_table"))))
      
    })
  })
  
  ## 3.2 draw the RNA-seq distribution data #################
  ## import the distribution data
  
  step3_rna_distr_plot <- reactive({
    req(input$act_step3_draw_rna_length)
    
    if (is.null(step3_rna_len_table())) {return(NULL)}
    
    source("R/draw_dot_line.R")
    
    # format the data
    item_name <- paste(input$in_step3_rna_strand, input$in_step3_rna_type, sep = '_')
    distr <- step3_rna_len_table() %>% 
      dplyr::filter(Length >= input$in_step3_rna_xlim[1] & Length <= input$in_step3_rna_xlim[2])
    
    # browser()
    line_plot <- draw_dot_line(distr = distr, x = "Length", y = item_name,
                               group = input$in_step3_rna_x, facet = input$in_step3_rna_wrap,
                               dot_size = input$in_step3_rna_dot_size, font_size = input$in_step3_rna_font_size,
                               line_color = input$in_step3_rna_line_color, line_width = input$in_step3_rna_line_width,
                               xstart = input$in_step3_rna_xlim[1], xend = input$in_step3_rna_xlim[2],
                               type = "line")
    
    heat_plot <- draw_dot_line(distr = distr, x = "Length", y = item_name, 
                               group = input$in_step3_rna_x, facet = input$in_step3_rna_wrap,
                               font_size = input$in_step3_rna_font_size,
                               fill_color = input$in_step3_rna_fill_color, fill_alpha = input$in_step3_rna_fill_alpha, 
                               xstart = input$in_step3_rna_xlim[1], xend = input$in_step3_rna_xlim[2],
                               type = "heat")
    
    return(list(line_plot = line_plot, heat_plot = heat_plot))
    
  })
  
  ## output the alignment plot
  observeEvent(input$act_step3_draw_rna_length, {
    # browser()
    ## show the alignment plot
    output$out_step3_rna_line <- renderPlot(
      width = input$out_step3_rna_width * 100,
      height = input$out_step3_rna_height * 100,
      {step3_rna_distr_plot()$line_plot})
    
    output$out_step3_rna_heat <- renderPlot(
      width = input$out_step3_rna_width * 100,
      height = input$out_step3_rna_height * 100,
      {step3_rna_distr_plot()$heat_plot})
    
    # save the alignment plot
    output$save_step3_rna_line <- downloadHandler(
      filename = function() {
        paste0(input$out_step3_rna_fig_name, "-line-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step3_rna_distr_plot()$line_plot, filename = file, 
               width = input$out_step3_rna_width, height = input$out_step3_rna_height)
      }
    )
    
    output$save_step3_rna_heat <- downloadHandler(
      filename = function() {
        paste0(input$out_step3_rna_fig_name, "-heat-plot-", Sys.Date(), ".pdf")
      },
      
      content = function(file) {
        ggsave(plot = step3_rna_distr_plot()$heat_plot, filename = file, 
               width = input$out_step3_rna_width, height = input$out_step3_rna_height)
      }
    )
    
  })
  
  
  ## 3.3 import the Ribo-seq distribution data #################
  ## import the distribution data
  step3_ribo_len_table <- reactive({
    req(input$act_step3_import_ribo_length)
    
    # import the data
    if (is.null(input$in_step3_riboseq)) {return(NULL)}
    
    ribo_length_table <- read.table(input$in_step3_riboseq$datapath, 
                                    sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    # browser()
    # calculate the average of group
    if (input$in_step3_ribo_x == 'Sample') {
      ribo_length_table <- ribo_length_table
      
    } else if (import_design_clicked()) {
      ribo_length_table <- ribo_length_table %>%
        dplyr::left_join(step1_design(), by = 'Sample') %>%
        na.omit() %>% 
        dplyr::group_by(!!sym(input$in_step3_ribo_x), Length) %>% 
        dplyr::reframe(Plus_Count = mean(Plus_Count), Minus_Count = mean(Minus_Count), 
                       Plus_Ratio = mean(Plus_Ratio), Minus_Ratio = mean(Minus_Ratio)) %>% 
        # dplyr::rename(Sample = Group) %>% 
        dplyr::ungroup() %>% 
        as.data.frame()
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    return(ribo_length_table)
  })
  
  ## output the table
  observeEvent(input$act_step3_import_ribo_length, {
    # browser()
    ## output the table
    output$out_step3_ribo_distr <- DT::renderDataTable(server = F, {
      DT::datatable(step3_ribo_len_table(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "ribo_length_table"), 
                                     list(extend = 'excel', filename = "ribo_length_table"))))
      
    })
  })
  
  ## 3.4 draw the Ribo-seq distribution data #################
  ## import the distribution data
  
  step3_ribo_distr_plot <- reactive({
    req(input$act_step3_draw_ribo_length)
    
    if (is.null(step3_ribo_len_table())) {return(NULL)}
    
    source("R/draw_dot_line.R")
    
    # format the data
    item_name <- paste(input$in_step3_ribo_strand, input$in_step3_ribo_type, sep = '_')
    distr <- step3_ribo_len_table() %>% 
      dplyr::filter(Length >= input$in_step3_ribo_xlim[1] & Length <= input$in_step3_ribo_xlim[2])
    
    # browser()
    line_plot <- draw_dot_line(distr = distr, x = "Length", y = item_name, 
                               group = input$in_step3_ribo_x, facet = input$in_step3_ribo_wrap,
                               dot_size = input$in_step3_ribo_dot_size, font_size = input$in_step3_ribo_font_size,
                               line_color = input$in_step3_ribo_line_color, line_width = input$in_step3_ribo_line_width,
                               xstart = input$in_step3_ribo_xlim[1], xend = input$in_step3_ribo_xlim[2],
                               type = "line")
    
    heat_plot <- draw_dot_line(distr = distr, x = "Length", y = item_name, 
                               group = input$in_step3_ribo_x, facet = input$in_step3_ribo_wrap,
                               font_size = input$in_step3_ribo_font_size,
                               fill_color = input$in_step3_ribo_fill_color, fill_alpha = input$in_step3_ribo_fill_alpha, 
                               xstart = input$in_step3_ribo_xlim[1], xend = input$in_step3_ribo_xlim[2],
                               type = "heat")
    
    return(list(line_plot = line_plot, heat_plot = heat_plot))
    
  })
  
  ## output the alignment plot
  observeEvent(input$act_step3_draw_ribo_length, {
    # browser()
    ## show the alignment plot
    output$out_step3_ribo_line <- renderPlot(
      width = input$out_step3_ribo_width * 100,
      height = input$out_step3_ribo_height * 100,
      {step3_ribo_distr_plot()$line_plot})
    
    output$out_step3_ribo_heat <- renderPlot(
      width = input$out_step3_ribo_width * 100,
      height = input$out_step3_ribo_height * 100,
      {step3_ribo_distr_plot()$heat_plot})
    
    # save the alignment plot
    output$save_step3_ribo_line <- downloadHandler(
      filename = function() {
        paste0(input$out_step3_ribo_fig_name, "-line-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step3_ribo_distr_plot()$line_plot, filename = file, 
               width = input$out_step3_ribo_width, height = input$out_step3_ribo_height)
      }
    )
    
    output$save_step3_ribo_heat <- downloadHandler(
      filename = function() {
        paste0(input$out_step3_ribo_fig_name, "-heat-plot-", Sys.Date(), ".pdf")
      },
      
      content = function(file) {
        ggsave(plot = step3_ribo_distr_plot()$heat_plot, filename = file, 
               width = input$out_step3_ribo_width, height = input$out_step3_ribo_height)
      }
    )
    
  })
  
  
  ############################################################
  # 4. figure of gene saturation #################
  
  ## 4.1 import the gene saturation data #################
  ## import the gene saturation
  step4_rna_saturation <- reactive({
    req(input$act_step4_import_rna_saturation)
    
    # import the data
    if (is.null(input$in_step4_rnaseq)) {return(NULL)}
    
    rna_gene_table <- read.table(input$in_step4_rnaseq$datapath, 
                                 sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    # browser()
    # calculate the average of group
    if (input$in_step4_rna_x == 'Sample') {
      rna_gene_table <- rna_gene_table
      
    } else if (import_design_clicked()) {
      rna_gene_table <- rna_gene_table %>%
        dplyr::left_join(step1_design(), by = 'Sample') %>%
        na.omit() %>% 
        dplyr::group_by(!!sym(input$in_step4_rna_x), Tercile) %>% 
        dplyr::reframe(Covered = mean(Covered), Uncovered = mean(Uncovered)) %>% 
        dplyr::ungroup() %>% 
        as.data.frame()
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    return(rna_gene_table)
  })
  
  ## output the table
  observeEvent(input$act_step4_import_rna_saturation, {
    # browser()
    ## output the table
    output$out_step4_rna_saturation <- DT::renderDataTable({
      DT::datatable(step4_rna_saturation(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "rna_saturation_table"), 
                                     list(extend = 'excel', filename = "rna_saturation_table"))))
      
    })
  })
  
  ## 4.2 draw the RNA-seq gene saturation #################
  ## import the gene saturation data
  
  step4_rna_saturation_plot <- reactive({
    req(input$act_step4_draw_rna_saturation)
    
    if (is.null(step4_rna_saturation())) {return(NULL)}
    
    # browser()
    source("R/draw_dot_line.R")
    
    item_name <- input$in_step4_rna_x
    
    # format the data
    if (input$in_step4_rna_label) {
      ypeak = max(step4_rna_saturation()$Covered)
    } else {
      ypeak = NULL
    }
    
    rna_saturation <- step4_rna_saturation() %>% 
      dplyr::filter(Tercile != 0) %>% 
      tidyr::pivot_longer(cols = c("Covered", "Uncovered"),
                          names_to = "Class", 
                          values_to = "Count") %>% 
      dplyr::mutate(Class = factor(Class, levels = c("Covered", "Uncovered")))
    
    # browser()
    line_plot <- draw_dot_line(distr = rna_saturation,
                               x = "Tercile", y = "Count", group = input$in_step4_rna_x, 
                               dot_size = input$in_step4_rna_dot_size, font_size = input$in_step4_rna_font_size, 
                               line_color = input$in_step4_rna_line_color, line_width = input$in_step4_rna_line_width,
                               xstart = 0, xend = 90, title = "Gene Saturation",
                               ypeak = ypeak, type = "line") +
      ggplot2::facet_wrap(~Class, ncol = 2, scales = "free_y") +
      ggplot2::theme(strip.background = ggplot2::element_blank())
    
    heat_plot <- draw_dot_line(distr = rna_saturation,
                               x = "Tercile", y = "Count", group = input$in_step4_rna_x,
                               fill_color = input$in_step4_rna_fill_color, fill_alpha = input$in_step4_rna_fill_alpha,
                               xstart = 0, xend = 90, title = "Gene Saturation",
                               type = "heat") +
      ggplot2::facet_wrap(~Class, ncol = 2, scales = "free_y") +
      ggplot2::theme(strip.background = ggplot2::element_blank())
    
    return(list(line_plot = line_plot, heat_plot = heat_plot))
    
  })
  
  ## output the alignment plot
  observeEvent(input$act_step4_draw_rna_saturation, {
    # browser()
    ## show the alignment plot
    output$out_step4_rna_line <- renderPlot(
      width = input$out_step4_rna_width * 100,
      height = input$out_step4_rna_height * 100,
      {step4_rna_saturation_plot()$line_plot})
    
    output$out_step4_rna_heat <- renderPlot(
      width = input$out_step4_rna_width * 100,
      height = input$out_step4_rna_height * 100,
      {step4_rna_saturation_plot()$heat_plot})
    
    # save the alignment plot
    output$save_step4_rna_line <- downloadHandler(
      filename = function() {
        paste0(input$out_step4_rna_fig_name, "-saturation-line-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step4_rna_saturation_plot()$line_plot, filename = file,
               width = input$out_step4_rna_width, height = input$out_step4_rna_height)
      }
    )
    
    output$save_step4_rna_heat <- downloadHandler(
      filename = function() {
        paste0(input$out_step4_rna_fig_name, "-saturation-heat-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step4_rna_saturation_plot()$heat_plot, filename = file,
               width = input$out_step4_rna_width, height = input$out_step4_rna_height)
      }
    )
    
  })
  
  
  ## 4.3 import the Ribo-seq gene saturation #################
  ## import the gene saturation data
  step4_ribo_saturation <- reactive({
    req(input$act_step4_import_ribo_saturation)
    
    # import the data
    if (is.null(input$in_step4_riboseq)) {return(NULL)}
    
    ribo_gene_table <- read.table(input$in_step4_riboseq$datapath, 
                                  sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    # browser()
    # calculate the average of group
    if (input$in_step4_ribo_x == 'Sample') {
      ribo_gene_table <- ribo_gene_table
      
    } else if (import_design_clicked()) {
      ribo_gene_table <- ribo_gene_table %>%
        dplyr::left_join(step1_design(), by = 'Sample') %>%
        na.omit() %>% 
        dplyr::group_by(!!sym(input$in_step4_ribo_x), Tercile) %>% 
        dplyr::reframe(Covered = mean(Covered), Uncovered = mean(Uncovered)) %>% 
        dplyr::ungroup() %>% 
        as.data.frame()
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    return(ribo_gene_table)
  })
  
  ## output the table
  observeEvent(input$act_step4_import_ribo_saturation, {
    # browser()
    ## output the table
    output$out_step4_ribo_saturation <- DT::renderDataTable({
      DT::datatable(step4_ribo_saturation(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "ribo_saturation_table"), 
                                     list(extend = 'excel', filename = "ribo_saturation_table"))))
      
    })
  })
  
  ## 4.4 draw the ribo-seq gene saturation #################
  ## import the gene saturation data
  
  step4_ribo_saturation_plot <- reactive({
    req(input$act_step4_draw_ribo_saturation)
    
    if (is.null(step4_ribo_saturation())) {return(NULL)}
    
    # browser()
    source("R/draw_dot_line.R")
    
    item_name <- input$in_step4_ribo_x
    
    # format the data
    if (input$in_step4_ribo_label) {
      ypeak = max(step4_ribo_saturation()$Covered)
    } else {
      ypeak = NULL
    }
    
    ribo_saturation <- step4_ribo_saturation() %>% 
      dplyr::filter(Tercile != 0) %>% 
      tidyr::pivot_longer(cols = c("Covered", "Uncovered"),
                          names_to = "Class", 
                          values_to = "Count") %>% 
      dplyr::mutate(Class = factor(Class, levels = c("Covered", "Uncovered")))
    
    # browser()
    line_plot <- draw_dot_line(distr = ribo_saturation,
                               x = "Tercile", y = "Count", group = input$in_step4_ribo_x, 
                               dot_size = input$in_step4_ribo_dot, font_size = input$in_step4_ribo_font_size, 
                               line_color = input$in_step4_ribo_line_color, line_width = input$in_step4_ribo_line_width,
                               xstart = 0, xend = 90, title = "Gene Saturation",
                               ypeak = ypeak, type = "line") +
      ggplot2::facet_wrap(~Class, ncol = 2, scales = "free_y") +
      ggplot2::theme(strip.background = ggplot2::element_blank())
    
    heat_plot <- draw_dot_line(distr = ribo_saturation,
                               x = "Tercile", y = "Count", group = input$in_step4_ribo_x,
                               fill_color = input$in_step4_ribo_fill_color, fill_alpha = input$in_step4_ribo_fill_alpha,
                               xstart = 0, xend = 90, title = "Gene Saturation",
                               type = "heat") +
      ggplot2::facet_wrap(~Class, ncol = 2, scales = "free_y") +
      ggplot2::theme(strip.background = ggplot2::element_blank())
    
    return(list(line_plot = line_plot, heat_plot = heat_plot))
    
  })
  
  ## output the alignment plot
  observeEvent(input$act_step4_draw_ribo_saturation, {
    # browser()
    ## show the alignment plot
    output$out_step4_ribo_line <- renderPlot(
      width = input$out_step4_ribo_width * 100,
      height = input$out_step4_ribo_height * 100,
      {step4_ribo_saturation_plot()$line_plot})
    
    output$out_step4_ribo_heat <- renderPlot(
      width = input$out_step4_ribo_width * 100,
      height = input$out_step4_ribo_height * 100,
      {step4_ribo_saturation_plot()$heat_plot})
    
    # save the alignment plot
    output$save_step4_ribo_line <- downloadHandler(
      filename = function() {
        paste0(input$out_step4_ribo_fig_name, "-saturation-line-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step4_ribo_saturation_plot()$line_plot, filename = file,
               width = input$out_step4_ribo_width, height = input$out_step4_ribo_height)
      }
    )
    
    output$save_step4_ribo_heat <- downloadHandler(
      filename = function() {
        paste0(input$out_step4_ribo_fig_name, "-saturation-heat-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step4_ribo_saturation_plot()$heat_plot, filename = file,
               width = input$out_step4_ribo_width, height = input$out_step4_ribo_height)
      }
    )
    
  })
  
  ############################################################
  # 5. figure of reads digestion #################
  
  ## 5.1 import the reads digestion data #################
  ## import the reads digestion
  step5_ribo_digestion <- reactive({
    req(input$act_step5_import_ribo_digestion)
    
    # import the data
    if (is.null(input$in_step5_riboseq)) {return(NULL)}
    
    digest_table <- read.table(input$in_step5_riboseq$datapath, 
                               sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    # browser()
    # calculate the average of group
    if (import_design_clicked()) {
      flt_design <- step1_design() %>% 
        dplyr::filter(SeqType == "RIBO") %>% 
        dplyr::select(Sample, !!sym(input$in_step5_ribo_x))
      
      digest_table <- digest_table %>%
        dplyr::left_join(flt_design, by = 'Sample') %>%
        na.omit() %>% 
        tidyr::pivot_longer(cols = -c(input$in_step5_ribo_x, Sample, Digest, Site), names_to = "Base", values_to = "PWM") %>%
        dplyr::group_by(!!sym(input$in_step5_ribo_x), Digest, Site, Base) %>% 
        dplyr::reframe(PWM = mean(PWM)) %>% 
        tidyr::pivot_wider(names_from = Base, values_from = PWM)
      
    }
    
    return(digest_table)
  })
  
  ## output the table
  observeEvent(input$act_step5_import_ribo_digestion, {
    # browser()
    ## output the table
    output$out_step5_ribo_digestion <- DT::renderDataTable(server = T, {
      DT::datatable(step5_ribo_digestion(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "digestion_table"), 
                                     list(extend = 'excel', filename = "digestion_table"))))
      
    })
  })
  
  ## 5.2 draw the Ribo-seq reads digestion #################
  ## import the reads digestion data
  
  step5_ribo_digestion_plot <- reactive({
    req(input$act_step5_draw_ribo_digestion)
    
    if (is.null(step5_ribo_digestion())) {return(NULL)}
    
    # browser()
    source("R/draw_seqlogo.R")
    
    # format the data
    extract_columns <- function(format_df) {
      format_df <- format_df %>% 
        dplyr::select(-c(input$in_step5_ribo_x, Digest)) %>% 
        tibble::remove_rownames() %>% 
        tibble::column_to_rownames("Site") %>%
        as.data.frame() %>% 
        t() %>% 
        as.matrix()
      
      return(format_df)
    }
    
    End_3p_df <- step5_ribo_digestion() %>% dplyr::filter(Digest == "End_3p")
    End_3p_list <- lapply(split(End_3p_df, End_3p_df[input$in_step5_ribo_x]), extract_columns)
    
    End_5p_df <- step5_ribo_digestion() %>% dplyr::filter(Digest == "End_5p")
    End_5p_list <- lapply(split(End_5p_df, End_5p_df[input$in_step5_ribo_x]), extract_columns)
    
    # browser()
    end5_plot <- draw_seqlogo(pwm = End_5p_list,
                              method = input$in_step5_ribo_method, seq_type = 'dna',
                              col_scheme = input$in_step5_ribo_fill_color, col_alpha = input$in_step5_ribo_fill_alpha,
                              stack_width = input$in_step5_ribo_font_stack, 
                              font_family = input$in_step5_ribo_font_family, font_size = input$in_step5_ribo_font_size,
                              xstart = -5, xend = 10)
    
    end3_plot <- draw_seqlogo(pwm = End_3p_list,
                              method = input$in_step5_ribo_method, seq_type = 'dna',
                              col_scheme = input$in_step5_ribo_fill_color, col_alpha = input$in_step5_ribo_fill_alpha,
                              stack_width = input$in_step5_ribo_font_stack, 
                              font_family = input$in_step5_ribo_font_family, font_size = input$in_step5_ribo_font_size,
                              xstart = -10, xend = 5)
    
    return(list(end5_plot = end5_plot, end3_plot = end3_plot))
    
  })
  
  ## output the alignment plot
  observeEvent(input$act_step5_draw_ribo_digestion, {
    # browser()
    ## show the alignment plot
    output$out_step5_ribo_seqlogo_end_5p <- renderPlot(
      width = input$out_step5_ribo_width * 100,
      height = input$out_step5_ribo_height * 100,
      {step5_ribo_digestion_plot()$end5_plot})
    
    output$out_step5_ribo_seqlogo_end_3p <- renderPlot(
      width = input$out_step5_ribo_width * 100,
      height = input$out_step5_ribo_height * 100,
      {step5_ribo_digestion_plot()$end3_plot})
    
    # save the alignment plot
    output$save_step5_ribo_seqlogo_end_5p <- downloadHandler(
      filename = function() {
        paste0(input$out_step5_ribo_fig_name, "-5end-seqlogo-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step5_ribo_digestion_plot()$end5_plot, filename = file,
               width = input$out_step5_ribo_width, height = input$out_step5_ribo_height)
      }
    )
    
    output$save_step5_ribo_seqlogo_end_3p <- downloadHandler(
      filename = function() {
        paste0(input$out_step5_ribo_fig_name, "-3end-seqlogo-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step5_ribo_digestion_plot()$end3_plot, filename = file,
               width = input$out_step5_ribo_width, height = input$out_step5_ribo_height)
      }
    )
    
  })
  
  
  ############################################################
  # 6. figure of reads offset #################
  
  ## 6.1 import the tis offset data #################
  ## import the tis offset
  step6_frame_offset <- reactive({
    req(input$act_step6_import_frame_offset)
    
    # import the data
    if (is.null(input$in_step6_frame)) {return(NULL)}
    
    frame_offset_table <- read.table(input$in_step6_frame$datapath, 
                                     sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>% 
      dplyr::select(-Ribo)
    
    # browser()
    # calculate the average of group
    if (input$in_step6_frame_x == 'Sample') {
      frame_offset_table <- frame_offset_table
      
    } else if (import_design_clicked()) {
      frame_offset_table <- frame_offset_table %>%
        dplyr::left_join(step1_design(), by = 'Sample') %>%
        na.omit() %>% 
        dplyr::group_by(!!sym(input$in_step6_frame_x), Model, Length) %>% 
        dplyr::reframe(Frame0 = mean(Frame0),
                       Rpfs0 = mean(Rpfs0),
                       Frame1 = mean(Frame1),
                       Rpfs1 = mean(Rpfs1),
                       Frame2 = mean(Frame2),
                       Rpfs2 = mean(Rpfs2),
                       P_site = max(P_site[which.max(Rpfs0)]),
                       Rpfs = mean(Rpfs),
                       Periodicity = mean(Periodicity)) %>% 
        dplyr::ungroup() %>% 
        as.data.frame()
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    return(frame_offset_table)
  })
  
  ## output the table
  observeEvent(input$act_step6_import_frame_offset, {
    # browser()
    ## output the table
    output$out_step6_frame_offset <- DT::renderDataTable(server = F, {
      DT::datatable(step6_frame_offset(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "frame_offset_table"), 
                                     list(extend = 'excel', filename = "frame_offset_table"))))
      
    })
  })
  
  ## 6.2 draw the frame offset plot #################
  
  step6_frame_offset_plot <- reactive({
    req(input$act_step6_draw_frame_offset)
    
    if (is.null(step6_frame_offset())) {return(NULL)}
    
    # browser()
    source("R/draw_offset.R")
    
    # format the data
    frame_offset <- step6_frame_offset() %>%
      dplyr::select(-Model, -Frame0, -Frame1, -Frame2, -Rpfs, -Periodicity) %>%
      tidyr::pivot_longer(cols = c(Rpfs0, Rpfs1, Rpfs2), 
                          names_to = 'Frame', 
                          values_to = 'Count') %>%
      dplyr::group_by(!!sym(input$in_step6_frame_x), Length) %>%
      dplyr::mutate(Ratio = Count / sum(Count) * 100) %>%
      dplyr::filter(Length >= input$in_step6_frame_xlim[1] & Length <= input$in_step6_frame_xlim[2])
    
    # browser()
    frame_bar <- draw_offset(offset = frame_offset, plot_type = 'bar',
                             x = 'Length', y = input$in_step6_frame_value,
                             wrap_group = input$in_step6_frame_x, fill_group = 'Frame',
                             fill_color = input$in_step6_frame_fill_color,
                             fill_alpha = input$in_step6_frame_fill_alpha,
                             font_size = input$in_step6_frame_font_size,
                             bar_width = input$in_step6_frame_bar_width,
                             facet = input$in_step6_frame_warp,
                             xstart = input$in_step6_frame_xlim[1], xend = input$in_step6_frame_xlim[2])
    
    frame_heat <- draw_offset(offset = frame_offset, plot_type = 'heat',
                              x = 'Length', y = input$in_step6_frame_value,
                              wrap_group = input$in_step6_frame_x, fill_group = 'Frame',
                              fill_color = input$in_step6_frame_fill_color,
                              fill_alpha = input$in_step6_frame_fill_alpha,
                              font_size = input$in_step6_frame_font_size,
                              bar_width = input$in_step6_frame_bar_width,
                              facet = input$in_step6_frame_warp,
                              xstart = input$in_step6_frame_xlim[1], xend = input$in_step6_frame_xlim[2])
    
    return(list(frame_bar = frame_bar, frame_heat = frame_heat))
    
  })
  
  ## output the frame offset plot
  observeEvent(input$act_step6_draw_frame_offset, {
    # browser()
    ## show the alignment plot
    output$out_step6_frame_bar <- renderPlot(
      width = input$out_step6_frame_width * 100,
      height = input$out_step6_frame_height * 100,
      {step6_frame_offset_plot()$frame_bar})
    
    output$out_step6_frame_heat <- renderPlot(
      width = input$out_step6_frame_width * 100,
      height = input$out_step6_frame_height * 100,
      {step6_frame_offset_plot()$frame_heat})
    
    # save the alignment plot
    output$save_step6_frame_bar <- downloadHandler(
      filename = function() {
        paste0(input$out_step6_frame_fig_name, "-frame-offset-bar-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step6_frame_offset_plot()$frame_bar, filename = file,
               width = input$out_step6_frame_width, height = input$out_step6_frame_height)
      }
    )
    
    output$save_step6_frame_heat <- downloadHandler(
      filename = function() {
        paste0(input$out_step6_frame_fig_name, "-frame-offset-heat-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step6_frame_offset_plot()$frame_heat, filename = file,
               width = input$out_step6_frame_width, height = input$out_step6_frame_height)
      }
    )
    
  })
  
  ## 6.3 import the tis offset data #################
  step6_tis_offset <- reactive({
    req(input$act_step6_import_tis_offset)
    
    # import the data
    if (is.null(input$in_step6_tis)) {return(NULL)}
    
    tis_offset_table <- read.table(input$in_step6_tis$datapath, 
                                   sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    # browser()
    # calculate the average of group
    if (input$in_step6_tis_x == 'Sample') {
      tis_offset_table <- tis_offset_table
      
    } else if (import_design_clicked()) {
      tis_offset_table <- tis_offset_table %>%
        dplyr::left_join(step1_design(), by = 'Sample') %>%
        na.omit() %>% 
        dplyr::group_by(!!sym(input$in_step6_tis_x), Model, Length, Ribo) %>% 
        dplyr::reframe(Frame0 = mean(Frame0),
                       Rpfs0 = mean(Rpfs0),
                       Frame1 = mean(Frame1),
                       Rpfs1 = mean(Rpfs1),
                       Frame2 = mean(Frame2),
                       Rpfs2 = mean(Rpfs2),
                       P_site = max(P_site[which.max(Rpfs0)]),
                       Rpfs = mean(Rpfs),
                       Periodicity = mean(Periodicity)) %>% 
        dplyr::ungroup() %>% 
        as.data.frame()
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    return(tis_offset_table)
  })
  
  ## output the table
  observeEvent(input$act_step6_import_tis_offset, {
    # browser()
    ## output the table
    output$out_step6_tis_offset <- DT::renderDataTable(server = F, {
      DT::datatable(step6_tis_offset(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "tis_offset_table"), 
                                     list(extend = 'excel', filename = "tis_offset_table"))))
      
    })
  })
  
  ## 6.4 draw the tis offset plot #################
  
  step6_tis_offset_plot <- reactive({
    req(input$act_step6_draw_tis_offset)
    
    if (is.null(step6_tis_offset())) {return(NULL)}
    
    # browser()
    source("R/draw_offset.R")
    
    # format the data
    tis_offset <- step6_tis_offset() %>% 
      dplyr::select(-Model, -Frame0, -Frame1, -Frame2, -Rpfs, -Periodicity) %>% 
      tidyr::pivot_longer(cols = c(Rpfs0, Rpfs1, Rpfs2), 
                          names_to = 'TIS', 
                          values_to = 'Count') %>% 
      dplyr::group_by(!!sym(input$in_step6_tis_x), Length) %>%
      dplyr::mutate(Ratio = Count / sum(Count)) %>%
      dplyr::filter(Length >= input$in_step6_tis_xlim[1] & Length <= input$in_step6_tis_xlim[2])
    
    # browser()
    tis_bar <- draw_offset(offset = tis_offset, plot_type = 'bar',
                           x = 'Length', y = input$in_step6_tis_value,
                           wrap_group = input$in_step6_tis_x, fill_group = 'TIS',
                           fill_color = input$in_step6_tis_fill_color,
                           fill_alpha = input$in_step6_tis_fill_alpha,
                           font_size = input$in_step6_tis_font_size,
                           bar_width = input$in_step6_tis_bar_width,
                           facet = input$in_step6_tis_warp,
                           xstart = input$in_step6_tis_xlim[1], xend = input$in_step6_tis_xlim[2])
    
    tis_heat <- draw_offset(offset = tis_offset, plot_type = 'heat',
                            x = 'Length', y = input$in_step6_tis_value,
                            wrap_group = input$in_step6_tis_x, fill_group = 'TIS',
                            fill_color = input$in_step6_tis_fill_color,
                            fill_alpha = input$in_step6_tis_fill_alpha,
                            font_size = input$in_step6_tis_font_size,
                            bar_width = input$in_step6_tis_bar_width,
                            facet = input$in_step6_tis_warp,
                            xstart = input$in_step6_tis_xlim[1], xend = input$in_step6_tis_xlim[2])
    
    return(list(tis_bar = tis_bar, tis_heat = tis_heat))
    
  })
  
  ## output the tis offset plot
  observeEvent(input$act_step6_draw_tis_offset, {
    # browser()
    ## show the alignment plot
    output$out_step6_tis_bar <- renderPlot(
      width = input$out_step6_tis_width * 100,
      height = input$out_step6_tis_height * 100,
      {step6_tis_offset_plot()$tis_bar})
    
    output$out_step6_tis_heat <- renderPlot(
      width = input$out_step6_tis_width * 100,
      height = input$out_step6_tis_height * 100,
      {step6_tis_offset_plot()$tis_heat})
    
    # save the alignment plot
    output$save_step6_tis_bar <- downloadHandler(
      filename = function() {
        paste0(input$out_step6_tis_fig_name, "-tis-offset-bar-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step6_tis_offset_plot()$tis_bar, filename = file,
               width = input$out_step6_tis_width, height = input$out_step6_tis_height)
      }
    )
    
    output$save_step6_tis_heat <- downloadHandler(
      filename = function() {
        paste0(input$out_step6_tis_fig_name, "-tts-offset-heat-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step6_tis_offset_plot()$tis_heat, filename = file,
               width = input$out_step6_tis_width, height = input$out_step6_tis_height)
      }
    )
    
  })
  
  ## 6.5 import the end of tis offset data #################
  ## import the end of tis offset
  step6_end_offset <- reactive({
    req(input$act_step6_import_end_offset)
    
    # import the data
    if (is.null(input$in_step6_end)) {return(NULL)}
    
    end_offset_table <- read.delim(input$in_step6_end$datapath, sep = '\t', header = TRUE,
                                   check.names = F, stringsAsFactors = FALSE)
    
    # browser()
    
    # calculate the average of group
    if (input$in_step6_end_group == 'Sample') {
      end_offset_table <- end_offset_table %>% 
        dplyr::group_by(Sample, Length, Site, End) %>% 
        dplyr::mutate(Log2 = log2(Count + 1),
                      Scaled = scale(Count)[, 1]) %>% 
        dplyr::ungroup() %>% 
        as.data.frame()
      
    } else if (import_design_clicked()) {
      end_offset_table <- end_offset_table %>%
        dplyr::left_join(step1_design(), by = 'Sample') %>%
        na.omit() %>% 
        dplyr::group_by(!!sym(input$in_step6_end_group), Length, Site, End) %>% 
        # dplyr::select_if(~ !is.character(.)) %>%
        dplyr::mutate(Count = mean(Count), 
                      Log2 = log2(Count + 1),
                      Scaled = scale(Count)[, 1]) %>% 
        dplyr::ungroup() %>% 
        as.data.frame()
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    return(end_offset_table)
  })
  
  ## output the table
  observeEvent(input$act_step6_import_end_offset, {
    # browser()
    ## output the table
    output$out_step6_end_offset <- DT::renderDataTable({
      DT::datatable(step6_end_offset(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "rpfs_end_offset_table"), 
                                     list(extend = 'excel', filename = "rpfs_end_offset_table"))))
      
    })
  })
  
  ## 6.6 draw the end of tis offset plot #################
  
  step6_end_offset_plot <- reactive({
    req(input$act_step6_draw_end_offset)
    
    if (is.null(step6_end_offset())) {return(NULL)}
    
    # browser()
    
    source("R/draw_offset_end.R")
    
    # format the data
    if (isTRUE(input$in_step6_end_wrap)) {
      end_offset <- step6_end_offset()
      
    } else {
      
      end_offset <- step6_end_offset() %>% 
        dplyr::select(-!!sym(input$in_step6_end_group)) %>% 
        dplyr::mutate(Offset = as.numeric(Offset)) %>%
        dplyr::group_by(Length, Site, End, Offset) %>%
        dplyr::reframe(Count = mean(Count)) %>% 
        dplyr::mutate(Log2 = log2(Count + 1)) %>% 
        dplyr::mutate(Scaled = scale(Count)[, 1]) %>%
        dplyr::ungroup()
      
    }
    
    end_offset <- end_offset %>% 
      dplyr::filter(Offset >= input$in_step6_end_x_min & Offset <= input$in_step6_end_x_max) %>% 
      dplyr::filter(Length >= input$in_step6_end_y_min & Length <= input$in_step6_end_x_max) %>% 
      dplyr::mutate(Length = factor(Length, levels = sort(unique(Length), decreasing = TRUE))) %>%
      dplyr::mutate(End = factor(End, c("5p-end", "3p-end"))) %>% 
      dplyr::mutate(Site = factor(Site, levels = c('TIS', 'TTS')))
    
    # browser()
    
    tis_heat_plot <- draw_offset_end(offset = end_offset %>% dplyr::filter(Site == 'TIS'),
                                     x = "Offset", 
                                     y = "Length", 
                                     
                                     xstart = input$in_step6_end_x_min,
                                     xend = input$in_step6_end_x_max,
                                     ystart = input$in_step6_end_y_min,
                                     yend = input$in_step6_end_y_max,
                                     
                                     xlabel = "End of RPFs",
                                     ylabel = "RPFs Length",
                                     title = "TIS Offset",
                                     
                                     font_size = input$in_step6_end_font_size,
                                     
                                     wrap_split = input$in_step6_end_wrap,
                                     wrap_group = input$in_step6_end_group,
                                     
                                     edge_color = input$in_step6_end_edge_color,
                                     edge_width = input$in_step6_end_edge_width,
                                     
                                     # xbreaks = input$in_step6_end_x_breaks,
                                     # ybreaks = input$in_step6_end_y_breaks,
                                     
                                     fill_group = input$in_step6_end_value,
                                     fill_color = input$in_step6_end_fill_color,
                                     fill_alpha = input$in_step6_end_fill_alpha)
    
    tts_heat_plot <- draw_offset_end(offset = end_offset %>% dplyr::filter(Site == 'TTS'),
                                     x = "Offset", 
                                     y = "Length", 
                                     
                                     xstart = input$in_step6_end_x_min,
                                     xend = input$in_step6_end_x_max,
                                     ystart = input$in_step6_end_y_min,
                                     yend = input$in_step6_end_y_max,
                                     
                                     xlabel = "End of RPFs",
                                     ylabel = "RPFs Length",
                                     title = "TTS Offset",
                                     
                                     font_size = input$in_step6_end_font_size,
                                     
                                     wrap_split = input$in_step6_end_wrap,
                                     wrap_group = input$in_step6_end_group,
                                     
                                     edge_color = input$in_step6_end_edge_color,
                                     edge_width = input$in_step6_end_edge_width,
                                     
                                     # xbreaks = input$in_step6_end_x_breaks,
                                     # ybreaks = input$in_step6_end_y_breaks,
                                     
                                     fill_group = input$in_step6_end_value,
                                     fill_color = input$in_step6_end_fill_color,
                                     fill_alpha = input$in_step6_end_fill_alpha)
    
    return(list(tis_heat_plot = tis_heat_plot, tts_heat_plot = tts_heat_plot))
    
  })
  
  ## output the tis offset plot
  observeEvent(input$act_step6_draw_end_offset, {
    # browser()
    ## show the alignment plot
    output$out_step6_tis_end_heat <- renderPlot(
      width = input$out_step6_end_width * 100,
      height = input$out_step6_end_height * 100,
      {step6_end_offset_plot()$tis_heat_plot})
    
    output$out_step6_tts_end_heat <- renderPlot(
      width = input$out_step6_end_width * 100,
      height = input$out_step6_end_height * 100,
      {step6_end_offset_plot()$tts_heat_plot})
    
    # save the alignment plot
    output$save_step6_tis_end_heat <- downloadHandler(
      filename = function() {
        paste0(input$out_step6_end_fig_name, "-tis-end-offset-heat-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step6_end_offset_plot()$tis_heat_plot, filename = file,
               width = input$out_step6_end_width, height = input$out_step6_end_height)
      }
    )
    
    output$save_step6_tts_end_heat <- downloadHandler(
      filename = function() {
        paste0(input$out_step6_end_fig_name, "-tts-end-offset-heat-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step6_end_offset_plot()$tts_heat_plot, filename = file,
               width = input$out_step6_end_width, height = input$out_step6_end_height)
      }
    )
    
  })
  
  
  ############################################################
  # 7. figure of 3-nt periodicity #################
  
  ## 7.1 import the 3-nt periodicity table  #################
  step7_rna_period <- reactive({
    req(input$act_step7_import_rna_period)
    
    # import the data
    if (is.null(input$in_step7_rnaseq)) {return(NULL)}
    
    rna_period_table <- read.table(input$in_step7_rnaseq$datapath, 
                                   sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    # browser()
    # calculate the average of group
    if (input$in_step7_rna_x == 'Sample') {
      rna_period_table <- rna_period_table
      
    } else if (import_design_clicked()) {
      rna_period_table <- rna_period_table %>%
        dplyr::left_join(step1_design(), by = 'Sample') %>%
        na.omit() %>% 
        dplyr::group_by(!!sym(input$in_step7_rna_x), Frame) %>% 
        dplyr::reframe(Count = mean(Count),
                       Ratio = mean(Ratio)) %>% 
        dplyr::ungroup() %>%
        as.data.frame()
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    return(rna_period_table)
  })
  
  ## output the table
  observeEvent(input$act_step7_import_rna_period, {
    # browser()
    ## output the table
    output$out_step7_rna_period <- DT::renderDataTable(server = F, {
      DT::datatable(step7_rna_period(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "periodicity_table"), 
                                     list(extend = 'excel', filename = "periodicity_table"))
                    ))
      
    })
  })
  
  ## 7.2 draw the RNA-seq 3-nt periodicity plot #################
  
  step7_rna_period_plot <- reactive({
    req(input$act_step7_draw_rna_period)
    
    if (is.null(step7_rna_period())) {return(NULL)}
    
    period_df = step7_rna_period() %>% 
      dplyr::mutate(Ratio = round(Ratio, 2)) %>% 
      dplyr::mutate(Label = paste0(round(Count / 1000000, 1), 'M')) %>% 
      dplyr::mutate(Frame = factor(Frame, levels = c('2', '1', '0')))
    
    # browser()
    source("R/draw_period.R")
    perid_count_plot <- draw_period(period = period_df,
                                    x = input$in_step7_rna_x, y = "Count",
                                    fill_group = 'Frame',
                                    fill_color = input$in_step7_rna_fill_color,
                                    fill_alpha = input$in_step7_rna_fill_alpha,
                                    edge_color = input$in_step7_rna_edge_color,
                                    edge_width = input$in_step7_rna_edge_width,
                                    bar_width = input$in_step7_rna_bar_width,
                                    font_size = input$in_step7_rna_font_size,
                                    label_size = input$in_step7_rna_label_size,
                                    trans = input$in_step7_rna_trans,
                                    text_label = input$in_step7_rna_text)
    
    perid_ratio_plot <- draw_period(period = period_df,
                                    x = input$in_step7_rna_x, y = "Ratio",
                                    fill_group = 'Frame',
                                    fill_color = input$in_step7_rna_fill_color,
                                    fill_alpha = input$in_step7_rna_fill_alpha,
                                    edge_color = input$in_step7_rna_edge_color,
                                    edge_width = input$in_step7_rna_edge_width,
                                    bar_width = input$in_step7_rna_bar_width,
                                    font_size = input$in_step7_rna_font_size,
                                    label_size = input$in_step7_rna_label_size,
                                    trans = input$in_step7_rna_trans,
                                    text_label = input$in_step7_rna_text)
    
    return(list(perid_count_plot = perid_count_plot, perid_ratio_plot = perid_ratio_plot))
    
  })
  
  ## output the RNA-seq 3-nt periodicity plot
  observeEvent(input$act_step7_draw_rna_period, {
    # browser()
    ## show the 3-nt periodicity plot
    output$out_step7_rna_count <- renderPlot(
      width = input$out_step7_rna_width * 100,
      height = input$out_step7_rna_height * 100,
      {step7_rna_period_plot()$perid_count_plot})
    
    output$out_step7_rna_ratio <- renderPlot(
      width = input$out_step7_rna_width * 100,
      height = input$out_step7_rna_height * 100,
      {step7_rna_period_plot()$perid_ratio_plot})
    
    # save the 3-nt periodicity plot
    output$save_step7_rna_count <- downloadHandler(
      filename = function() {
        paste0(input$out_step7_rna_fig_name, "-period-count-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step7_rna_period_plot()$perid_count_plot, filename = file,
               width = input$out_step7_rna_width, height = input$out_step7_rna_height)
      }
    )
    
    output$save_step7_rna_ratio <- downloadHandler(
      filename = function() {
        paste0(input$out_step7_rna_fig_name, "-period-ratio-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step7_rna_period_plot()$perid_ratio_plot, filename = file,
               width = input$out_step7_rna_width, height = input$out_step7_rna_height)
      }
    )
    
  })
  
  
  ## 7.3 import the Ribo-seq 3-nt periodicity plot #################
  step7_ribo_period <- reactive({
    req(input$act_step7_import_ribo_period)
    
    # import the data
    if (is.null(input$in_step7_riboseq)) {return(NULL)}
    
    ribo_period_table <- read.table(input$in_step7_riboseq$datapath, 
                                    sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    # browser()
    # calculate the average of group
    if (input$in_step7_ribo_x == 'Sample') {
      ribo_period_table <- ribo_period_table
      
    } else if (import_design_clicked()) {
      ribo_period_table <- ribo_period_table %>%
        dplyr::left_join(step1_design(), by = 'Sample') %>%
        na.omit() %>% 
        dplyr::group_by(!!sym(input$in_step7_ribo_x), Frame) %>% 
        dplyr::reframe(Count = mean(Count),
                       Ratio = mean(Ratio)) %>% 
        dplyr::ungroup() %>%
        as.data.frame()
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    return(ribo_period_table)
  })
  
  ## output the table
  observeEvent(input$act_step7_import_ribo_period, {
    # browser()
    ## output the table
    output$out_step7_ribo_period <- DT::renderDataTable(server = F, {
      DT::datatable(step7_ribo_period(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "periodicity_table"), 
                                     list(extend = 'excel', filename = "periodicity_table"))
                    ))
      
    })
  })
  
  ## 7.4 draw the Ribo-seq 3-nt periodicity plot #################
  
  step7_ribo_period_plot <- reactive({
    req(input$act_step7_draw_ribo_period)
    
    if (is.null(step7_ribo_period())) {return(NULL)}
    
    period_df = step7_ribo_period() %>% 
      dplyr::mutate(Ratio = round(Ratio, 2)) %>% 
      dplyr::mutate(Label = paste0(round(Count / 1000000, 1), 'M')) %>% 
      dplyr::mutate(Frame = factor(Frame, levels = c('2', '1', '0')))
    
    # browser()
    source("R/draw_period.R")
    perid_count_plot <- draw_period(period = period_df,
                                    x = input$in_step7_ribo_x, y = "Count",
                                    fill_group = 'Frame',
                                    fill_color = input$in_step7_ribo_fill_color,
                                    fill_alpha = input$in_step7_ribo_fill_alpha,
                                    edge_color = input$in_step7_ribo_edge_color,
                                    edge_width = input$in_step7_ribo_edge_width,
                                    bar_width = input$in_step7_ribo_bar_width,
                                    font_size = input$in_step7_ribo_font_size,
                                    label_size = input$in_step7_ribo_label_size,
                                    trans = input$in_step7_ribo_trans,
                                    text_label = input$in_step7_ribo_text)
    
    perid_ratio_plot <- draw_period(period = period_df,
                                    x = input$in_step7_ribo_x, y = "Ratio",
                                    fill_group = 'Frame',
                                    fill_color = input$in_step7_ribo_fill_color,
                                    fill_alpha = input$in_step7_ribo_fill_alpha,
                                    edge_color = input$in_step7_ribo_edge_color,
                                    edge_width = input$in_step7_ribo_edge_width,
                                    bar_width = input$in_step7_ribo_bar_width,
                                    font_size = input$in_step7_ribo_font_size,
                                    label_size = input$in_step7_ribo_label_size,
                                    trans = input$in_step7_ribo_trans,
                                    text_label = input$in_step7_ribo_text)
    
    return(list(perid_count_plot = perid_count_plot, perid_ratio_plot = perid_ratio_plot))
    
  })
  
  ## output the 3-nt periodicity plot
  observeEvent(input$act_step7_draw_ribo_period, {
    # browser()
    ## show the 3-nt periodicity plot
    output$out_step7_ribo_count <- renderPlot(
      width = input$out_step7_ribo_width * 100,
      height = input$out_step7_ribo_height * 100,
      {step7_ribo_period_plot()$perid_count_plot})
    
    output$out_step7_ribo_ratio <- renderPlot(
      width = input$out_step7_ribo_width * 100,
      height = input$out_step7_ribo_height * 100,
      {step7_ribo_period_plot()$perid_ratio_plot})
    
    # save the 3-nt periodicity plot
    output$save_step7_ribo_count <- downloadHandler(
      filename = function() {
        paste0(input$out_step7_ribo_fig_name, "-period-count-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step7_ribo_period_plot()$perid_count_plot, filename = file,
               width = input$out_step7_ribo_width, height = input$out_step7_ribo_height)
      }
    )
    
    output$save_step7_ribo_ratio <- downloadHandler(
      filename = function() {
        paste0(input$out_step7_ribo_fig_name, "-period-ratio-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step7_ribo_period_plot()$perid_ratio_plot, filename = file,
               width = input$out_step7_ribo_width, height = input$out_step7_ribo_height)
      }
    )
    
  })
  
  
  ############################################################
  # 8. figure of metaplot #################
  
  ## 8.1 import the RNA-seq metaplot table  #################
  step8_rna_meta <- reactive({
    
    req(input$act_step8_import_rna_meta)
    
    if (is.null(input$in_step8_rnaseq)) {return(NULL)}
    
    rna_meta_table <- read.table(input$in_step8_rnaseq$datapath, 
                                 sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    # browser()
    # scale the data
    if (input$in_step8_rna_scale) {
      rna_meta_table <- rna_meta_table %>% 
        dplyr::group_by(!!sym(input$in_step8_rna_x)) %>%
        dplyr::mutate(Density = Density / mean(Density)) %>% 
        dplyr::ungroup()
    }
    
    # calculate the average of group
    if (input$in_step8_rna_x == 'Sample') {
      rna_meta_table <- rna_meta_table
      
    } else if (import_design_clicked()) {
      rna_meta_table <- rna_meta_table %>%
        dplyr::left_join(step1_design(), by = 'Sample') %>%
        na.omit() %>% 
        dplyr::group_by(!!sym(input$in_step8_rna_x), Meta, Nucleotide, Codon, Frame) %>% 
        dplyr::reframe(Density = mean(Density)) %>% 
        dplyr::ungroup() %>%
        as.data.frame()
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    return(rna_meta_table)
  })
  
  ## output the table
  observeEvent(input$act_step8_import_rna_meta, {
    # browser()
    
    ## output the table
    output$out_step8_rna_meta <- DT::renderDataTable({
      DT::datatable(step8_rna_meta(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "metagene_table", 
                                          exportOptions = list(modifier = list(page = 'all'))), 
                                     list(extend = 'excel', filename = "metagene_table"))
                    ))
      
    })
  })
  
  ## 8.2 draw the RNA-seq metaplot #################
  
  step8_rna_meta_plot <- reactive({
    req(input$act_step8_draw_rna_meta)
    
    if (is.null(step8_rna_meta())) {return(NULL)}
    
    # browser()
    source("R/draw_metaplot.R")
    meta_df = step8_rna_meta() %>% 
      dplyr::filter(Meta == "TIS" & Codon >= input$in_step8_rna_tis_x_min & Codon <= input$in_step8_rna_tis_x_max |
                      Meta == "TTS" & Codon >= input$in_step8_rna_tts_x_min & Codon <= input$in_step8_rna_tts_x_max) %>% 
      dplyr::mutate(Frame = factor(Frame, levels = c("0", "1", "2")))
    
    meta_bar_plot <- draw_metaplot(meta = meta_df, plot_type = "Bar",
                                   x = "Nucleotide", y = "Density",
                                   group = "Frame", font_size = input$in_step8_rna_font_size,
                                   fill_color = input$in_step8_rna_fill_color,
                                   fill_alpha = input$in_step8_rna_fill_alpha,
                                   facet = input$in_step8_rna_wrap,
                                   wrap_group = input$in_step8_rna_x)
    
    meta_line_plot <- draw_metaplot(meta = meta_df, plot_type = "Line",
                                    x = "Nucleotide", y = "Density",
                                    group = input$in_step8_rna_x, font_size = input$in_step8_rna_font_size,
                                    line_color = input$in_step8_rna_line_color,
                                    line_width = input$in_step8_rna_line_width,
                                    facet = input$in_step8_rna_wrap,
                                    wrap_group = input$in_step8_rna_x)
    
    meta_heat_plot <- draw_metaplot(meta = meta_df, plot_type = "Heat",
                                    x = "Nucleotide", y = input$in_step8_rna_x,
                                    group = 'Density', font_size = input$in_step8_rna_font_size,
                                    fill_color = input$in_step8_rna_fill_color,
                                    fill_alpha = input$in_step8_rna_fill_alpha,
                                    facet = input$in_step8_rna_wrap,
                                    wrap_group = input$in_step8_rna_x)
    
    return(list(meta_bar_plot = meta_bar_plot, meta_line_plot = meta_line_plot, meta_heat_plot = meta_heat_plot))
    
  })
  
  ## output the 3-nt periodicity plot
  observeEvent(input$act_step8_draw_rna_meta, {
    # browser()
    ## show the 3-nt periodicity plot
    output$out_step8_rna_bar <- renderPlot(
      width = input$out_step8_rna_width * 100,
      height = input$out_step8_rna_height * 100,
      {step8_rna_meta_plot()$meta_bar_plot})
    
    output$out_step8_rna_line <- renderPlot(
      width = input$out_step8_rna_width * 100,
      height = input$out_step8_rna_height * 100,
      {step8_rna_meta_plot()$meta_line_plot})
    
    output$out_step8_rna_heat <- renderPlot(
      width = input$out_step8_rna_width * 100,
      height = input$out_step8_rna_height * 100,
      {step8_rna_meta_plot()$meta_heat_plot})
    
    # save the 3-nt periodicity plot
    output$save_step8_rna_bar <- downloadHandler(
      filename = function() {
        paste0(input$out_step8_rna_fig_name, "-meta-bar-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step8_rna_meta_plot()$meta_bar_plot, filename = file,
               width = input$out_step8_rna_width, height = input$out_step8_rna_height)
      }
    )
    
    output$save_step8_rna_line <- downloadHandler(
      filename = function() {
        paste0(input$out_step8_rna_fig_name, "-meta-line-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step8_rna_meta_plot()$meta_line_plot, filename = file,
               width = input$out_step8_rna_width, height = input$out_step8_rna_height)
      }
    )
    
    output$save_step8_rna_heat <- downloadHandler(
      filename = function() {
        paste0(input$out_step8_rna_fig_name, "-meta-heat-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step8_rna_meta_plot()$meta_heat_plot, filename = file,
               width = input$out_step8_rna_width, height = input$out_step8_rna_height)
      }
    )
    
  })

  
  ## 8.3 import the Ribo-seq metaplot table  #################
  step8_ribo_meta <- reactive({
    
    req(input$act_step8_import_ribo_meta)
    
    if (is.null(input$in_step8_riboseq)) {return(NULL)}
    
    ribo_meta_table <- read.table(input$in_step8_riboseq$datapath, 
                                  sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    # browser()
    # scale the data
    if (input$in_step8_ribo_scale) {
      ribo_meta_table <- ribo_meta_table %>% 
        dplyr::group_by(!!sym(input$in_step8_ribo_x)) %>%
        dplyr::mutate(Density = Density / mean(Density)) %>% 
        dplyr::ungroup()
    }
    
    # calculate the average of group
    if (input$in_step8_ribo_x == 'Sample') {
      ribo_meta_table <- ribo_meta_table

    } else if (import_design_clicked()) {
      ribo_meta_table <- ribo_meta_table %>%
        dplyr::left_join(step1_design(), by = 'Sample') %>%
        na.omit() %>% 
        dplyr::group_by(!!sym(input$in_step8_ribo_x), Meta, Nucleotide, Codon, Frame) %>% 
        dplyr::reframe(Density = mean(Density)) %>% 
        dplyr::ungroup() %>%
        as.data.frame()
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    return(ribo_meta_table)
  })
  
  ## output the table
  observeEvent(input$act_step8_import_ribo_meta, {
    # browser()
    ## output the table
    output$out_step8_ribo_meta <- DT::renderDataTable({
      DT::datatable(step8_ribo_meta(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "metagene_table", 
                                          exportOptions = list(modifier = list(page = 'all'))), 
                                     list(extend = 'excel', filename = "metagene_table"))
                    ))
      
    })
  })

  
  ## 8.4 draw the Ribo-seq metaplot #################
  
  step8_ribo_meta_plot <- reactive({
    req(input$act_step8_draw_ribo_meta)
    
    if (is.null(step8_ribo_meta())) {return(NULL)}
    
    # browser()
    if (input$in_step8_ribo_x_type == "Nucleotide") {
      meta_df = step8_ribo_meta() %>% 
        dplyr::filter(Meta == "TIS" & Codon >= input$in_step8_ribo_tis_x_min & Codon <= input$in_step8_ribo_tis_x_max |
                        Meta == "TTS" & Codon >= input$in_step8_ribo_tts_x_min & Codon <= input$in_step8_ribo_tts_x_max) %>% 
        dplyr::mutate(Frame = factor(Frame, levels = c("0", "1", "2")))
      bar_group <- "Frame"

    } else if (input$in_step8_ribo_x_type == "Codon") {
      meta_df = step8_ribo_meta() %>% 
        dplyr::filter(Meta == "TIS" & Codon >= input$in_step8_ribo_tis_x_min & Codon <= input$in_step8_ribo_tis_x_max |
                        Meta == "TTS" & Codon >= input$in_step8_ribo_tts_x_min & Codon <= input$in_step8_ribo_tts_x_max) %>% 
        dplyr::group_by(!!sym(input$in_step8_ribo_x), Meta, Codon) %>%
        dplyr::reframe(Density = mean(Density))
      
      bar_group <- "Meta"
      
      if (input$in_step8_ribo_smooth > 0) {
        meta_df <- meta_df %>% 
          dplyr::group_by(!!sym(input$in_step8_ribo_x), Meta) %>%
          # dplyr::mutate(Density = loess(Density ~ seq_along(Density), span = input$in_step8_ribo_loess)$fitted) %>%
          dplyr::mutate(Density = rollmean(Density, k = input$in_step8_ribo_smooth, fill = NA)) %>%
          dplyr::ungroup()
      }
      
    }

    source("R/draw_metaplot.R")
    
    meta_bar_plot <- draw_metaplot(meta = meta_df, plot_type = "Bar",
                                   x = input$in_step8_ribo_x_type, y = "Density",
                                   ystart = input$in_step8_ribo_y_min, yend = input$in_step8_ribo_y_max,
                                   group = bar_group, font_size = input$in_step8_ribo_font_size,
                                   fill_color = input$in_step8_ribo_fill_color,
                                   fill_alpha = input$in_step8_ribo_fill_alpha,
                                   legend_position = input$in_step8_ribo_legend,
                                   facet = input$in_step8_ribo_wrap,
                                   wrap_group = input$in_step8_ribo_x)
    
    meta_line_plot <- draw_metaplot(meta = meta_df, plot_type = "Line",
                                    x = input$in_step8_ribo_x_type, y = "Density",
                                    ystart = input$in_step8_ribo_y_min, yend = input$in_step8_ribo_y_max,
                                    group = input$in_step8_ribo_x, font_size = input$in_step8_ribo_font_size,
                                    line_color = input$in_step8_ribo_line_color,
                                    line_width = input$in_step8_ribo_line_width,
                                    line_alpha = input$in_step8_ribo_fill_alpha,
                                    legend_position = input$in_step8_ribo_legend,
                                    facet = input$in_step8_ribo_wrap,
                                    wrap_group = input$in_step8_ribo_x)
    
    meta_heat_plot <- draw_metaplot(meta = meta_df, plot_type = "Heat",
                                    x = input$in_step8_ribo_x_type, y = input$in_step8_ribo_x,
                                    group = 'Density', font_size = input$in_step8_ribo_font_size,
                                    fill_color = input$in_step8_ribo_fill_color,
                                    fill_alpha = input$in_step8_ribo_fill_alpha,
                                    legend_position = input$in_step8_ribo_legend,
                                    facet = input$in_step8_ribo_wrap,
                                    wrap_group = input$in_step8_ribo_x)
    
    return(list(meta_bar_plot = meta_bar_plot, meta_line_plot = meta_line_plot, meta_heat_plot = meta_heat_plot))
    
  })
  
  ## output the 3-nt periodicity plot
  observeEvent(input$act_step8_draw_ribo_meta, {
    # browser()
    ## show the 3-nt periodicity plot
    output$out_step8_ribo_bar <- renderPlot(
      width = input$out_step8_ribo_width * 100,
      height = input$out_step8_ribo_height * 100,
      {step8_ribo_meta_plot()$meta_bar_plot})
    
    output$out_step8_ribo_line <- renderPlot(
      width = input$out_step8_ribo_width * 100,
      height = input$out_step8_ribo_height * 100,
      {step8_ribo_meta_plot()$meta_line_plot})
    
    output$out_step8_ribo_heat <- renderPlot(
      width = input$out_step8_ribo_width * 100,
      height = input$out_step8_ribo_height * 100,
      {step8_ribo_meta_plot()$meta_heat_plot})
    
    # save the 3-nt periodicity plot
    output$save_step8_ribo_bar <- downloadHandler(
      filename = function() {
        paste0(input$out_step8_ribo_fig_name, "-meta-bar-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step8_ribo_meta_plot()$meta_bar_plot, filename = file,
               width = input$out_step8_ribo_width, height = input$out_step8_ribo_height)
      }
    )
    
    output$save_step8_ribo_line <- downloadHandler(
      filename = function() {
        paste0(input$out_step8_ribo_fig_name, "-meta-line-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step8_ribo_meta_plot()$meta_line_plot, filename = file,
               width = input$out_step8_ribo_width, height = input$out_step8_ribo_height)
      }
    )
    
    output$save_step8_ribo_heat <- downloadHandler(
      filename = function() {
        paste0(input$out_step8_ribo_fig_name, "-meta-heat-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step8_ribo_meta_plot()$meta_heat_plot, filename = file,
               width = input$out_step8_ribo_width, height = input$out_step8_ribo_height)
      }
    )
    
  })
  
  
  ############################################################
  # 9. figure of coverage #################
  
  ## 9.1 import the RNA-seq coverage data #################
  step9_rna_coverage <- reactive({
    
    req(input$act_step9_import_rna_coverage)
    
    if (is.null(input$in_step9_rnaseq)) {return(NULL)}
    
    rna_coverage_table <- read.table(input$in_step9_rnaseq$datapath, 
                                     sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    # browser()
    # calculate the average of group
    if (input$in_step9_rna_x == 'Sample') {
      rna_coverage_table <- rna_coverage_table
      
      asymmetry_score <- rna_coverage_table %>%
        dplyr::filter(Region == "CDS") %>%
        dplyr::group_by(Sample) %>%
        dplyr::arrange(Bins) %>%
        summarise(
          first_half = sum(Density[1:(n()/2)]),
          second_half = sum(Density[(n()/2 + 1):n()]),
          asymmetry_score = log2(second_half / first_half)
        )
      
    } else if (import_design_clicked()) {
      rna_coverage_table <- rna_coverage_table %>%
        dplyr::left_join(step1_design(), by = 'Sample') %>%
        na.omit() %>% 
        dplyr::group_by(!!sym(input$in_step9_rna_x), Region, Bins) %>% 
        dplyr::reframe(Density = mean(Density)) %>% 
        dplyr::ungroup() %>%
        as.data.frame()
      
      asymmetry_score <- rna_coverage_table %>%
        dplyr::filter(Region == "CDS") %>%
        dplyr::group_by(Sample) %>%
        dplyr::arrange(Bins) %>%
        summarise(
          first_half = sum(Density[1:(n()/2)]),
          second_half = sum(Density[(n()/2 + 1):n()]),
          asymmetry_score = log2(second_half / first_half)
        )
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    # browser()
    # scale the data
    if (input$in_step9_rna_scale) {
      rna_coverage_table <- rna_coverage_table %>% 
        dplyr::group_by(!!sym(input$in_step9_rna_x)) %>%
        dplyr::mutate(Density = Density / mean(Density)) %>% 
        dplyr::ungroup()
    }
    
    return(list(rna_coverage_table = rna_coverage_table,
                asymmetry_score = asymmetry_score))
  })
  
  ## output the table
  observeEvent(input$act_step9_import_rna_coverage, {
    # browser()
    ## output the table
    output$out_step9_rna_coverage <- DT::renderDataTable(server = T, {
      DT::datatable(step9_rna_coverage()$rna_coverage_table, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "coverage_table"), 
                                     list(extend = 'excel', filename = "coverage_table"))
                    ))
      
    })
    
    ## output the table
    output$out_step9_rna_as_score <- DT::renderDataTable(server = T, {
      DT::datatable(step9_rna_coverage()$asymmetry_score, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "asymmetry_score_table"), 
                                     list(extend = 'excel', filename = "asymmetry_score_table"))
                    ))
      
    })
    
  })
  
  ## 9.2 draw the RNA-seq coverage #################
  
  step9_rna_coverage_plot <- reactive({
    req(input$act_step9_draw_rna_coverage)
    
    if (is.null(step9_rna_coverage()$rna_coverage_table)) {return(NULL)}
    
    # browser()
    source("R/draw_coverage.R")
    coverage_df = step9_rna_coverage()$rna_coverage_table %>% 
      dplyr::mutate(Region = factor(Region, levels = c("5-UTR", "CDS", "3-UTR")))
    
    coverage_line_plot <- draw_coverage(coverage = coverage_df, plot_type = "Line",
                                        x = "Bins", y = "Density",
                                        group = input$in_step9_rna_x, font_size = input$in_step9_rna_font_size,
                                        line_color = input$in_step9_rna_line_color,
                                        fill_alpha = input$in_step9_rna_fill_alpha,
                                        line_width = input$in_step9_rna_line_width,
                                        facet = input$in_step9_rna_wrap,
                                        wrap_group = input$in_step9_rna_x)
    
    coverage_heat_plot <- draw_coverage(coverage = coverage_df, plot_type = "Heat",
                                        x = "Bins", y = input$in_step9_rna_x,
                                        group = 'Density', font_size = input$in_step9_rna_font_size,
                                        fill_color = input$in_step9_rna_fill_color,
                                        fill_alpha = input$in_step9_rna_fill_alpha,
                                        facet = input$in_step9_rna_wrap,
                                        wrap_group = input$in_step9_rna_x)
    
    # browser()
    
    # draw the schematic of gene body
    utr5_length <- sum(str_count(step9_rna_coverage()$rna_coverage_table$Region, '5-UTR'))
    cds_length <- sum(str_count(step9_rna_coverage()$rna_coverage_table$Region, 'CDS'))
    utr3_length <- sum(str_count(step9_rna_coverage()$rna_coverage_table$Region, '3-UTR'))
    gene_length <- max(step9_rna_coverage()$rna_coverage_table$Bins)
    
    source("R/draw_Isoforms_Schematic.R")
    isoforms_rna_schematic <- draw_isoforms_schematic(gene_name = "Gene body",
                                                      utr5_length = utr5_length,
                                                      cds_length = cds_length,
                                                      utr3_length = utr3_length,
                                                      gene_color = c("#75aadb", "#264abd", "#75aadb"),
                                                      gene_width = c(2, 5, 2),
                                                      fill_alpha = input$in_step9_rna_fill_alpha,
                                                      xmin = NA,
                                                      xmax = NA,
                                                      xbreaks = 4,
                                                      arrow_length = 0.1,
                                                      arrow_width = 0.8,
                                                      arrow_color = "white",
                                                      font_size = input$in_step9_rna_font_size) +
      theme_void() +
      theme(legend.position = "bottom",
            legend.title = element_blank())
    
    # merge the figures
    coverage_line_plot <- cowplot::plot_grid(coverage_line_plot, isoforms_rna_schematic, 
                                             ncol = 1, align = 'v', rel_heights = c(4,1))
    
    coverage_heat_plot <- cowplot::plot_grid(coverage_heat_plot, isoforms_rna_schematic, 
                                             ncol = 1, align = 'v', rel_heights = c(4,1))
    
    return(list(coverage_line_plot = coverage_line_plot, coverage_heat_plot = coverage_heat_plot))
    
  })
  
  ## output the 3-nt periodicity plot
  observeEvent(input$act_step9_draw_rna_coverage, {
    # browser()
    ## show the 3-nt periodicity plot
    output$out_step9_rna_line <- renderPlot(
      width = input$out_step9_rna_width * 100,
      height = input$out_step9_rna_height * 100,
      {step9_rna_coverage_plot()$coverage_line_plot})
    
    output$out_step9_rna_heat <- renderPlot(
      width = input$out_step9_rna_width * 100,
      height = input$out_step9_rna_height * 100,
      {step9_rna_coverage_plot()$coverage_heat_plot})
    
    # save the 3-nt periodicity plot
    output$save_step9_rna_line <- downloadHandler(
      filename = function() {
        paste0(input$out_step9_rna_fig_name, "-coverage-line-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step9_rna_coverage_plot()$coverage_line_plot, filename = file,
               width = input$out_step9_rna_width, height = input$out_step9_rna_height)
      }
    )
    
    output$save_step9_rna_heat <- downloadHandler(
      filename = function() {
        paste0(input$out_step9_rna_fig_name, "-coverage-heat-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step9_rna_coverage_plot()$coverage_heat_plot, filename = file,
               width = input$out_step9_rna_width, height = input$out_step9_rna_height)
      }
    )
    
  })
  
  
  ## 9.3 import the Ribo-seq coverage data #################
  step9_ribo_coverage <- reactive({
    
    req(input$act_step9_import_ribo_coverage)
    
    if (is.null(input$in_step9_riboseq)) {return(NULL)}
    
    ribo_coverage_table <- read.table(input$in_step9_riboseq$datapath, 
                                      sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    # browser()
    
    # calculate the average of group
    if (input$in_step9_ribo_x == 'Sample') {
      ribo_coverage_table <- ribo_coverage_table
      
      asymmetry_score <- ribo_coverage_table %>%
        dplyr::filter(Region == "CDS") %>%
        dplyr::group_by(Sample) %>%
        dplyr::arrange(Bins) %>%
        summarise(
          first_half = sum(Density[1:(n()/2)]),
          second_half = sum(Density[(n()/2 + 1):n()]),
          asymmetry_score = log2(second_half / first_half)
        )
      
    } else if (import_design_clicked()) {
      ribo_coverage_table <- ribo_coverage_table %>%
        dplyr::left_join(step1_design(), by = 'Sample') %>%
        na.omit() %>% 
        dplyr::group_by(!!sym(input$in_step9_ribo_x), Region, Bins) %>% 
        dplyr::reframe(Density = mean(Density)) %>% 
        dplyr::ungroup() %>%
        as.data.frame()
      
      asymmetry_score <- ribo_coverage_table %>%
        dplyr::filter(Region == "CDS") %>%
        dplyr::group_by(!!sym(input$in_step9_ribo_x)) %>%
        dplyr::arrange(Bins) %>%
        summarise(
          first_half = sum(Density[1:(n()/2)]),
          second_half = sum(Density[(n()/2 + 1):n()]),
          asymmetry_score = log2(second_half / first_half)
        )
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    # browser()
    # scale the data
    if (input$in_step9_ribo_scale) {
      ribo_coverage_table <- ribo_coverage_table %>% 
        dplyr::group_by(!!sym(input$in_step9_ribo_x)) %>%
        dplyr::mutate(Density = Density / mean(Density)) %>% 
        dplyr::ungroup()
    }
    
    return(list(ribo_coverage_table = ribo_coverage_table,
                asymmetry_score = asymmetry_score))
  })
  
  ## output the table
  observeEvent(input$act_step9_import_ribo_coverage, {
    # browser()
    ## output the table
    output$out_step9_ribo_coverage <- DT::renderDataTable(server = T, {
      DT::datatable(step9_ribo_coverage()$ribo_coverage_table, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "coverage_table"), 
                                     list(extend = 'excel', filename = "coverage_table"))
                    ))
      
    })
    
    ## output the table
    output$out_step9_ribo_as_score <- DT::renderDataTable(server = F, {
      DT::datatable(step9_ribo_coverage()$asymmetry_score, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "asymmetry_score_table"), 
                                     list(extend = 'excel', filename = "asymmetry_score_table"))
                    ))
      
    })
  })
  
  ## 9.4 draw the Ribo-seq coverage #################
  
  step9_ribo_coverage_plot <- reactive({
    
    req(input$act_step9_draw_ribo_coverage)
    
    if (is.null(step9_ribo_coverage()$ribo_coverage_table)) {return(NULL)}
    
    # browser()
    source("R/draw_coverage.R")
    coverage_df = step9_ribo_coverage()$ribo_coverage_table
    
    coverage_line_plot <- draw_coverage(coverage = coverage_df, plot_type = "Line",
                                        x = "Bins", y = "Density",
                                        group = input$in_step9_ribo_x, font_size = input$in_step9_ribo_font_size,
                                        line_color = input$in_step9_ribo_line_color,
                                        fill_alpha = input$in_step9_ribo_fill_alpha,
                                        line_width = input$in_step9_ribo_line_width,
                                        facet = input$in_step9_ribo_wrap,
                                        wrap_group = input$in_step9_ribo_x)
    
    coverage_heat_plot <- draw_coverage(coverage = coverage_df, plot_type = "Heat",
                                        x = "Bins", y = input$in_step9_ribo_x,
                                        group = 'Density', font_size = input$in_step9_ribo_font_size,
                                        fill_color = input$in_step9_ribo_fill_color,
                                        fill_alpha = input$in_step9_ribo_fill_alpha,
                                        facet = input$in_step9_ribo_wrap,
                                        wrap_group = input$in_step9_ribo_x)
    
    # browser()
    
    # draw the schematic of gene body
    utr5_length <- sum(str_count(step9_ribo_coverage()$ribo_coverage_table$Region, '5-UTR'))
    cds_length <- sum(str_count(step9_ribo_coverage()$ribo_coverage_table$Region, 'CDS'))
    utr3_length <- sum(str_count(step9_ribo_coverage()$ribo_coverage_table$Region, '3-UTR'))
    gene_length <- max(step9_ribo_coverage()$ribo_coverage_table$Bins)
    
    source("R/draw_Isoforms_Schematic.R")
    isoforms_ribo_schematic <- draw_isoforms_schematic(gene_name = "Gene body",
                                                       utr5_length = utr5_length,
                                                       cds_length = cds_length,
                                                       utr3_length = utr3_length,
                                                       gene_color = c("#75aadb", "#264abd", "#75aadb"),
                                                       gene_width = c(2, 5, 2),
                                                       fill_alpha = input$in_step9_ribo_fill_alpha,
                                                       xmin = NA,
                                                       xmax = NA,
                                                       xbreaks = 4,
                                                       arrow_length = 0.1,
                                                       arrow_width = 0.8,
                                                       arrow_color = "white",
                                                       font_size = input$in_step9_ribo_font_size) +
      theme_void() +
      theme(legend.position = "bottom",
            legend.title = element_blank())
    
    # merge the figures
    coverage_line_plot <- cowplot::plot_grid(coverage_line_plot, isoforms_ribo_schematic, 
                                             ncol = 1, align = 'v', rel_heights = c(4,1))
    
    coverage_heat_plot <- cowplot::plot_grid(coverage_heat_plot, isoforms_ribo_schematic, 
                                             ncol = 1, align = 'v', rel_heights = c(4,1))
    
    return(list(coverage_line_plot = coverage_line_plot, coverage_heat_plot = coverage_heat_plot))
    
  })
  
  ## output the 3-nt periodicity plot
  observeEvent(input$act_step9_draw_ribo_coverage, {
    # browser()
    ## show the 3-nt periodicity plot
    output$out_step9_ribo_line <- renderPlot(
      width = input$out_step9_ribo_width * 100,
      height = input$out_step9_ribo_height * 100,
      {step9_ribo_coverage_plot()$coverage_line_plot})
    
    output$out_step9_ribo_heat <- renderPlot(
      width = input$out_step9_ribo_width * 100,
      height = input$out_step9_ribo_height * 100,
      {step9_ribo_coverage_plot()$coverage_heat_plot})
    
    # save the 3-nt periodicity plot
    output$save_step9_ribo_line <- downloadHandler(
      filename = function() {
        paste0(input$out_step9_ribo_fig_name, "-coverage-line-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step9_ribo_coverage_plot()$coverage_line_plot, filename = file,
               width = input$out_step9_ribo_width, height = input$out_step9_ribo_height)
      }
    )
    
    output$save_step9_ribo_heat <- downloadHandler(
      filename = function() {
        paste0(input$out_step9_ribo_fig_name, "-coverage-heat-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step9_ribo_coverage_plot()$coverage_heat_plot, filename = file,
               width = input$out_step9_ribo_width, height = input$out_step9_ribo_height)
      }
    )
    
  })
  
  
  ############################################################
  # 10. figure of abundance #################
  
  ## 10.1 import the RNA-seq abundance #################
  step10_rna_count <- reactive({
    
    req(input$act_step10_import_rna_count)
    
    if (is.null(input$in_step10_rnaseq)) {return(NULL)}
    
    # import the raw data table
    suffix <- c(".genes.results", ".isoforms.results", "_cds_rpf", "_cds_rpm", "_cds_rpkm", "_cds_tpm", "Aligned.sortedByCoord.out.bam")
    ftcount_column <- c("Chr", "Start", "End", "Strand", "Length")
    
    rna_count_table <- fread(input$in_step10_rnaseq$datapath, sep = '\t', check.names = FALSE,
                             header="auto", stringsAsFactors = FALSE) %>% 
      dplyr::rename_with(~ gsub(paste(suffix, collapse = '|'), '', .), everything()) %>%
      dplyr::rename(Gene = colnames(.)[1]) %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::select(-contains(ftcount_column)) %>% 
      dplyr::mutate(across(where(is.numeric), as.double)) %>% 
      dplyr::filter(rowSums(.) > input$in_step10_rna_filter)
    
    rna_count_table <- rna_count_table %>%
      dplyr::rename_all(~basename(names(rna_count_table)))
    
    # browser()
    
    # select the samples
    if (import_design_clicked()) {
      rna_count_table <- rna_count_table %>% 
        dplyr::select(any_of(step1_design()$Sample)) %>% 
        dplyr::filter(rowSums(.) > input$in_step10_rna_filter)
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    rna_count_table <- rna_count_table %>% 
      tibble::rownames_to_column("Gene")
    
    return(rna_count_table)
  })
  
  ## output the table
  observeEvent(input$act_step10_import_rna_count, {
    # browser()
    
    ## output the table
    output$out_step10_rna_count <- DT::renderDataTable(server = T, {
      DT::datatable(step10_rna_count(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "RNA_seq_count"), 
                                     list(extend = 'excel', filename = "RNA_seq_count"))
                    ))
    })
  })
  
  ## 10.2 normalize the RNA-seq abundance #################
  step10_rna_norm <- reactive({
    
    req(input$act_step10_norm_rna_count)
    
    if (is.null(step10_rna_count())) {return(NULL)}
    
    # browser()
    source("R/run_ruvseq_norm.R")
    
    if (import_design_clicked()) {
      # get the sample information
      phenodata <- step1_design() %>% 
        dplyr::filter(Sample %in% colnames(step10_rna_count())) %>%
        dplyr::filter(SeqType == input$in_step10_rna_seq) %>%
        dplyr::mutate(Group = factor(!!sym(input$in_step10_rna_group), 
                                     levels = unique(!!sym(input$in_step10_rna_group)))) %>% 
        tibble::column_to_rownames(var = "Sample")
      
      # normalize the RNA-seq abundance
      rna_norm <- run_ruvseq_norm(count = step10_rna_count() %>% tibble::column_to_rownames("Gene"),
                                  phenodata = phenodata,
                                  group = input$in_step10_rna_group,
                                  method = input$in_step10_rna_method,
                                  fill_color = input$in_step10_rna_fill_color,
                                  fill_alpha = input$in_step10_rna_fill_alpha,
                                  profile = input$in_step10_rna_seq)
      
      return(list(rna_norm = rna_norm$count_set2_norm, 
                  rna_rpm = rna_norm$count_set2_rpm,
                  phenodata = phenodata,
                  outmess = NULL))
      
    } else {
      return(list(rna_norm = NULL, 
                  rna_rpm = NULL,
                  phenodata = NULL,
                  outmess = "Step 1 Design has not been imported."))
    }
    
  })
  
  ## output the table
  observeEvent(input$act_step10_norm_rna_count, {
    # browser()
    output$out_step10_rna_norm <- DT::renderDataTable(server = T, {
      DT::datatable(step10_rna_norm()$rna_norm, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc', 
                                   buttons = list(
                                     list(extend = 'csv', filename = "RNA_seq_norm", page = "all"), 
                                     list(extend = 'excel', filename = "RNA_seq_norm", page = "all"))
                    ))
    })
    
    output$out_step10_rna_rpm <- DT::renderDataTable(server = T, {
      DT::datatable(step10_rna_norm()$rna_rpm, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "RNA_seq_rpm"),
                                     list(extend = 'excel', filename = "RNA_seq_rpm"))
                    ))
    })
    
    # save the expression level table
    output$out_step10_rna_norm_info <- renderText({
      step10_rna_norm()$outmess
    })
    
    output$save_step10_rna_norm <- downloadHandler(
      filename = function() {
        paste0(input$out_step10_rna_name, "-normalized-", Sys.Date(), input$out_step10_rna_format)
      },
      
      content = function(file) {
        if (input$out_step10_rna_format == ".txt") {
          write.table(x = step10_rna_norm()$rna_norm, file = file, row.names = FALSE, sep = '\t', quote = F)
        } else if (input$out_step10_rna_format == ".xlsx") {
          write.xlsx(x = step10_rna_norm()$rna_norm, file = file, rowNames = FALSE)
        }
      }
    )
    
    output$save_step10_rna_rpm <- downloadHandler(
      filename = function() {
        paste0(input$out_step10_rna_name, "-RPM-", Sys.Date(), input$out_step10_rna_format)
      },
      content = function(file) {
        if (input$out_step10_rna_format == ".txt") {
          write.table(x = step10_rna_norm()$rna_rpm, file = file, row.names = FALSE, sep = '\t', quote = F)
        } else if (input$out_step10_rna_format == ".xlsx") {
          write.xlsx(x = step10_rna_norm()$rna_rpm, file = file, rowNames = FALSE)
        }
      }
    )
    
  })
  
  ## 10.3 draw the RNA-seq abundance #################
  step10_rna_plot <- reactive({
    
    req(input$act_step10_draw_rna_expr)
    
    if (is.null(step10_rna_count())) {return(NULL)}
    if (is.null(step10_rna_norm())) {return(NULL)}
    
    source("R/draw_RLE.R")
    source("R/calc_RLE.R")
    source("R/draw_PCA.R")
    source("R/draw_eCDF.R")
    
    group <- input$in_step10_rna_group
    fill_color <- input$in_step10_rna_fill_color
    fill_alpha <- input$in_step10_rna_fill_alpha
    profile <- input$in_step10_rna_seq
    line_color <- input$in_step10_rna_line_color
    line_width <- input$in_step10_rna_line_width
    facet <- input$in_step10_rna_wrap
    log2_transform <- input$in_step10_rna_log2
    # xrange <- input$in_step10_rna_xlim
    # yrange <- input$in_step10_rna_ylim
    font_size <- input$in_step10_rna_font_size
    
    x_min <- input$in_step10_ribo_xlim_min
    x_max <- input$in_step10_ribo_xlim_max
    
    y_min <- input$in_step10_ribo_ylim_min
    y_max <- input$in_step10_ribo_ylim_max
    
    # browser()
    phenodata <- step10_rna_norm()$phenodata %>% tibble::rownames_to_column(var = "Sample")
    
    # draw the RLE plot
    raw_rle_mat <- calculate_rle(step10_rna_count() %>% tibble::column_to_rownames("Gene"))
    raw_rle_longer <- make_rle_mat(raw_rle_mat, phenodata)
    raw_rle_longer <- raw_rle_longer %>% dplyr::mutate(Sample = factor(Sample, levels = unique(phenodata$Sample)))
    raw_rle_plot <- draw_rle(rle_mat = raw_rle_longer, x = "Sample", y = "RLE", title = paste(profile, "(Raw count)"),
                             ymin = y_min, ymax = y_max, group = group, font_size = font_size,
                             fill_color = fill_color, fill_alpha = fill_alpha)
    
    norm_rle_mat <- calculate_rle(step10_rna_norm()$rna_norm %>% tibble::column_to_rownames("Gene"))
    norm_rle_longer <- make_rle_mat(norm_rle_mat, phenodata)
    norm_rle_longer <- norm_rle_longer %>% dplyr::mutate(Sample = factor(Sample, levels = unique(phenodata$Sample)))
    norm_rle_plot <- draw_rle(rle_mat = norm_rle_longer, x = "Sample", y = "RLE", title = paste(profile, "(Normalized)"),
                              ymin = y_min, ymax = y_max, group = group, font_size = font_size,
                              fill_color = fill_color, fill_alpha = fill_alpha)
    
    # draw the PCA plot
    raw_pca_res <- draw_pca(in_data = step10_rna_count() %>% tibble::column_to_rownames("Gene"), 
                            phenodata = phenodata, font_size = font_size,
                            group = group, title = paste(profile, "(Raw count)"),
                            fill_color = fill_color, fill_alpha = fill_alpha)
    
    norm_pca_res <- draw_pca(in_data = step10_rna_norm()$rna_norm %>% tibble::column_to_rownames("Gene"), 
                             phenodata = phenodata, font_size = font_size,
                             group = group, title = paste(profile, "(Normalized)"),
                             fill_color = fill_color, fill_alpha = fill_alpha)
    
    # browser()
    # draw the eCDF
    count_longer <-  step10_rna_count() %>% 
      tidyr::pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Count") %>%
      dplyr::left_join(phenodata, by = "Sample")
    
    norm_longer <- step10_rna_norm()$rna_norm %>% 
      tidyr::pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Normal") %>%
      dplyr::left_join(phenodata, by = "Sample")
    
    # log2 transform
    if (log2_transform) {
      count_longer <- count_longer %>% 
        dplyr::mutate(Count = log2(Count + 1)) %>% 
        dplyr::filter(Count > x_min, Count < x_max)
      
      norm_longer <- norm_longer %>% 
        dplyr::mutate(Normal = log2(Normal + 1)) %>% 
        dplyr::filter(Normal > x_min, Normal < x_max)
    }
    
    raw_ecdf_plot <- draw_ecdf(in_data = count_longer, 
                               x = "Count", 
                               group = group, 
                               title = paste(profile, "(Raw count)"),
                               xlabel = "Log2(Count + 1)",
                               facet = facet, 
                               font_size = font_size,
                               line_color = line_color, 
                               line_width = line_width)
    
    norm_ecdf_plot <- draw_ecdf(in_data = norm_longer, 
                                x = "Normal", 
                                group = group, 
                                title = paste(profile, "(Normalized)"),
                                xlabel = "Log2(Normalized + 1)",
                                facet = facet, 
                                font_size = font_size,
                                line_color = line_color, 
                                line_width = line_width)
    
    # format the plots
    rle_plot <- cowplot::plot_grid(raw_rle_plot, norm_rle_plot, ncol = 1)
    pca_plot <- cowplot::plot_grid(raw_pca_res$pca_scatter, norm_pca_res$pca_scatter, ncol = 2)
    ecdf_plot <- cowplot::plot_grid(raw_ecdf_plot, norm_ecdf_plot, ncol = 2)
    
    # browser()
    return(list(rle_plot = rle_plot,
                pca_plot = pca_plot,
                ecdf_plot = ecdf_plot))
  })
  
  ## output the table
  observeEvent(input$act_step10_draw_rna_expr, {
    # browser()
    ## show the expression level plot
    output$out_step10_rna_box <- renderPlot(
      width = input$out_step10_rna_width * 75,
      height = input$out_step10_rna_height * 150,
      {step10_rna_plot()$rle_plot})
    
    output$out_step10_rna_pca <- renderPlot(
      width = input$out_step10_rna_width * 100,
      height = input$out_step10_rna_height * 100,
      {step10_rna_plot()$pca_plot})
    
    output$out_step10_rna_ecdf <- renderPlot(
      width = input$out_step10_rna_width * 100,
      height = input$out_step10_rna_height * 100,
      {step10_rna_plot()$ecdf_plot})
    
    # save the expression level plot
    output$save_step10_rna_box <- downloadHandler(
      filename = function() {
        paste0(input$out_step10_rna_name, "-RLE-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step10_rna_plot()$rle_plot, filename = file,
               width = input$out_step10_rna_width, height = input$out_step10_rna_height)
      }
    )
    
    output$save_step10_rna_pca <- downloadHandler(
      filename = function() {
        paste0(input$out_step10_rna_name, "-temp-PCA-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step10_rna_plot()$pca_plot, filename = file,
               width = input$out_step10_rna_width, height = input$out_step10_rna_height)
      }
    )
    
    output$save_step10_rna_ecdf <- downloadHandler(
      filename = function() {
        paste0(input$out_step10_rna_name, "-eCDF-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step10_rna_plot()$ecdf_plot, filename = file,
               width = input$out_step10_rna_width, height = input$out_step10_rna_height)
      }
    )
    
  })
  
  
  ## 10.4 import the Ribo-seq abundance #################
  step10_ribo_count <- reactive({
    
    req(input$act_step10_import_ribo_count)
    
    if (is.null(input$in_step10_riboseq)) {return(NULL)}
    
    # import the raw data table
    suffix <- c(".genes.results", ".isoforms.results", "_cds_rpf", "_cds_rpm", "_cds_rpkm", "_cds_tpm", "Aligned.sortedByCoord.out.bam")
    ftcount_column <- c("Chr", "Start", "End", "Strand", "Length")
    
    ribo_count_table <- fread(input$in_step10_riboseq$datapath, sep = '\t', check.names = FALSE,
                              header="auto", stringsAsFactors = FALSE) %>% 
      dplyr::rename_with(~ gsub(paste(suffix, collapse = '|'), '', .), everything()) %>%
      dplyr::rename(Gene = colnames(.)[1]) %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::select(-contains(ftcount_column)) %>% 
      dplyr::mutate(across(where(is.numeric), as.double)) %>% 
      dplyr::filter(rowSums(.) > input$in_step10_ribo_filter)
    
    ribo_count_table <- ribo_count_table %>%
      dplyr::rename_all(~basename(names(ribo_count_table)))
    
    # select the samples
    if (import_design_clicked()) {
      ribo_count_table <- ribo_count_table %>% 
        dplyr::select(any_of(step1_design()$Sample)) %>% 
        dplyr::filter(rowSums(.) > input$in_step10_ribo_filter)
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    ribo_count_table <- ribo_count_table %>% 
      tibble::rownames_to_column("Gene")
    
    return(ribo_count_table)
  })
  
  ## output the table
  observeEvent(input$act_step10_import_ribo_count, {
    # browser()
    
    ## output the table
    output$out_step10_ribo_count <- DT::renderDataTable(server = T, {
      DT::datatable(step10_ribo_count(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "RNA_seq_count"), 
                                     list(extend = 'excel', filename = "RNA_seq_count"))
                    ))
    })
  })
  
  ## 10.5 normalize the Ribo-seq abundance #################
  
  step10_ribo_norm <- reactive({
    
    req(input$act_step10_norm_ribo_count)
    
    if (is.null(step10_ribo_count())) {return(NULL)}
    
    # browser()
    source("R/run_ruvseq_norm.R")
    
    if (import_design_clicked()) {
      # get the sample information
      phenodata <- step1_design() %>%
        dplyr::filter(Sample %in% colnames(step10_ribo_count())) %>% 
        dplyr::filter(SeqType == input$in_step10_ribo_seq) %>%
        dplyr::mutate(Group = factor(!!sym(input$in_step10_ribo_group), 
                                     levels = unique(!!sym(input$in_step10_ribo_group)))) %>% 
        tibble::column_to_rownames(var = "Sample")
      
      # normalize the RNA-seq abundance
      ribo_norm <- run_ruvseq_norm(count = step10_ribo_count() %>% tibble::column_to_rownames("Gene"),
                                   phenodata = phenodata,
                                   group = input$in_step10_ribo_group,
                                   method = input$in_step10_ribo_method,
                                   fill_color = input$in_step10_ribo_fill_color,
                                   fill_alpha = input$in_step10_ribo_fill_alpha,
                                   profile = input$in_step10_ribo_seq)
      
      return(list(ribo_norm = ribo_norm$count_set2_norm, 
                  ribo_rpm = ribo_norm$count_set2_rpm,
                  phenodata = phenodata,
                  outmess = NULL))
      
    } else {
      return(list(ribo_norm = ribo_norm$count_set2_norm, 
                  ribo_rpm = ribo_norm$count_set2_rpm,
                  phenodata = phenodata,
                  outmess = "Step 1 Design has not been imported."))
    }
    
  })
  
  ## output the table
  output$out_step10_rna_norm_info <- renderText({
    step10_ribo_norm()$outmess
  })
  
  observeEvent(input$act_step10_norm_ribo_count, {
    # browser()
    output$out_step10_ribo_norm <- DT::renderDataTable(server = T, {
      DT::datatable(step10_ribo_norm()$ribo_norm, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "Ribo_seq_rpm", page = "all"), 
                                     list(extend = 'excel', filename = "Ribo_seq_rpm", page = "all"))
                    ))
    })
    
    output$out_step10_ribo_rpm <- DT::renderDataTable(server = T, {
      DT::datatable(step10_ribo_norm()$ribo_rpm, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "Ribo_seq_rpm"), 
                                     list(extend = 'excel', filename = "Ribo_seq_rpm"))
                    ))
    })
    
    # save the expression level table
    output$save_step10_ribo_norm <- downloadHandler(
      filename = function() {
        paste0(input$out_step10_ribo_name, "-normalized-", Sys.Date(), input$out_step10_ribo_format)
      },
      content = function(file) {
        if (input$out_step10_ribo_format == ".txt") {
          write.table(x = step10_ribo_norm()$ribo_norm, file = file, row.names = FALSE, sep = '\t', quote = F)
        } else if (input$out_step10_ribo_format == ".csv") {
          write.xlsx(x = step10_ribo_norm()$ribo_norm, file = file, row.names = FALSE)
        } else if (input$out_step10_ribo_format == ".xlsx") {
          write.xlsx(x = step10_ribo_norm()$ribo_norm, file = file, rowNames = FALSE)
        }
      }
    )
    
    output$save_step10_ribo_rpm <- downloadHandler(
      filename = function() {
        paste0(input$out_step10_ribo_name, "-RPM-", Sys.Date(), input$out_step10_ribo_format)
      },
      content = function(file) {
        if (input$out_step10_ribo_format == ".txt") {
          write.table(x = step10_ribo_norm()$ribo_rpm, file = file, row.names = FALSE, sep = '\t', quote = F)
        } else if (input$out_step10_ribo_format == ".csv") {
          write.xlsx(x = step10_ribo_norm()$ribo_rpm, file = file, row.names = FALSE)
        } else if (input$out_step10_ribo_format == ".xlsx") {
          write.xlsx(x = step10_ribo_norm()$ribo_rpm, file = file, rowNames = FALSE)
        }
      }
    )
    
  })
  
  ## 10.6 draw the Ribo-seq abundance #################
  step10_ribo_plot <- reactive({
    
    req(input$act_step10_draw_ribo_expr)
    
    if (is.null(step10_ribo_count())) {return(NULL)}
    if (is.null(step10_ribo_norm())) {return(NULL)}
    
    source("R/draw_RLE.R")
    source("R/calc_RLE.R")
    source("R/draw_PCA.R")
    source("R/draw_eCDF.R")
    
    group <- input$in_step10_ribo_group
    fill_color <- input$in_step10_ribo_fill_color
    fill_alpha <- input$in_step10_ribo_fill_alpha
    profile <- input$in_step10_ribo_seq
    line_color <- input$in_step10_ribo_line_color
    line_width <- input$in_step10_ribo_line_width
    facet <- input$in_step10_ribo_wrap
    log2_transform <- input$in_step10_ribo_log2
    
    x_min <- input$in_step10_ribo_xlim_min
    x_max <- input$in_step10_ribo_xlim_max
    
    y_min <- input$in_step10_ribo_ylim_min
    y_max <- input$in_step10_ribo_ylim_max
    
    font_size = input$in_step10_ribo_font_size
    
    # browser()
    phenodata <- step10_ribo_norm()$phenodata %>% tibble::rownames_to_column(var = "Sample")
    
    # draw the RLE plot
    raw_rle_mat <- calculate_rle(step10_ribo_count() %>% tibble::column_to_rownames("Gene"))
    raw_rle_longer <- make_rle_mat(raw_rle_mat, phenodata)
    raw_rle_longer <- raw_rle_longer %>% dplyr::mutate(Sample = factor(Sample, levels = unique(phenodata$Sample)))
    raw_rle_plot <- draw_rle(rle_mat = raw_rle_longer, 
                             x = "Sample", y = "RLE", title = paste(profile, "(Raw count)"),
                             ymin = y_min, ymax = y_max, group = group, font_size = font_size,
                             fill_color = fill_color, fill_alpha = fill_alpha)
    
    norm_rle_mat <- calculate_rle(step10_ribo_norm()$ribo_norm %>% tibble::column_to_rownames("Gene"))
    norm_rle_longer <- make_rle_mat(norm_rle_mat, phenodata)
    norm_rle_longer <- norm_rle_longer %>% dplyr::mutate(Sample = factor(Sample, levels = unique(phenodata$Sample)))
    norm_rle_plot <- draw_rle(rle_mat = norm_rle_longer, x = "Sample", y = "RLE", title = paste(profile, "(Normalized)"),
                              ymin = y_min, ymax = y_max, group = group, font_size = font_size,
                              fill_color = fill_color, fill_alpha = fill_alpha)
    
    # draw the PCA plot
    raw_pca_res <- draw_pca(in_data = step10_ribo_count() %>% tibble::column_to_rownames("Gene"), 
                            phenodata = phenodata, font_size = font_size,
                            group = group, title = paste(profile, "(Raw count)"),
                            fill_color = fill_color, fill_alpha = fill_alpha)
    
    norm_pca_res <- draw_pca(in_data = step10_ribo_norm()$ribo_norm %>% tibble::column_to_rownames("Gene"), 
                             phenodata = phenodata, font_size = font_size,
                             group = group, title = paste(profile, "(Normalized)"),
                             fill_color = fill_color, fill_alpha = fill_alpha)
    
    # browser()
    # draw the eCDF
    count_longer <- step10_ribo_count() %>% 
      tidyr::pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Count") %>%
      dplyr::left_join(phenodata, by = "Sample")
    
    norm_longer <- step10_ribo_norm()$ribo_norm %>% 
      tidyr::pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Normal") %>%
      dplyr::left_join(phenodata, by = "Sample")
    
    # log2 transform
    if (log2_transform) {
      count_longer <- count_longer %>% 
        dplyr::mutate(Count = log2(Count + 1)) %>% 
        dplyr::filter(Count >= x_min, Count <= x_max)
      
      norm_longer <- norm_longer %>% 
        dplyr::mutate(Normal = log2(Normal + 1)) %>% 
        dplyr::filter(Normal >= x_min, Normal <= x_max)
    }
    
    raw_ecdf_plot <- draw_ecdf(in_data = count_longer, 
                               x = "Count", 
                               group = group, 
                               title = paste(profile, "(Raw count)"),
                               xlabel = "log2(Count + 1)",
                               facet = facet, font_size = font_size,
                               line_color = line_color, line_width = line_width)
    
    norm_ecdf_plot <- draw_ecdf(in_data = norm_longer,
                                x = "Normal", 
                                group = group, 
                                title = paste(profile, "(Normalized)"),
                                xlabel = "log2(Normalized + 1)",
                                facet = facet, font_size = font_size,
                                line_color = line_color, line_width = line_width)
    
    # format the plots
    rle_plot <- cowplot::plot_grid(raw_rle_plot, norm_rle_plot, ncol = 1)
    pca_plot <- cowplot::plot_grid(raw_pca_res$pca_scatter, norm_pca_res$pca_scatter, ncol = 2)
    ecdf_plot <- cowplot::plot_grid(raw_ecdf_plot, norm_ecdf_plot, ncol = 2)
    
    # browser()
    return(list(rle_plot = rle_plot,
                pca_plot = pca_plot,
                ecdf_plot = ecdf_plot))
  })
  
  ## output the table
  observeEvent(input$act_step10_draw_ribo_expr, {
    
    ## show the expression level plot
    output$out_step10_ribo_box <- renderPlot(
      width = input$out_step10_ribo_width * 70,
      height = input$out_step10_ribo_height * 150,
      {step10_ribo_plot()$rle_plot})
    
    output$out_step10_ribo_pca <- renderPlot(
      width = input$out_step10_ribo_width * 100,
      height = input$out_step10_ribo_height * 100,
      {step10_ribo_plot()$pca_plot})
    
    output$out_step10_ribo_ecdf <- renderPlot(
      width = input$out_step10_ribo_width * 100,
      height = input$out_step10_ribo_height * 100,
      {step10_ribo_plot()$ecdf_plot})
    
    # save the expression level plot
    output$save_step10_ribo_box <- downloadHandler(
      filename = function() {
        paste0(input$out_step10_ribo_name, "-RLE-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step10_ribo_plot()$rle_plot, filename = file,
               width = input$out_step10_ribo_width, height = input$out_step10_ribo_height)
      }
    )
    
    output$save_step10_ribo_pca <- downloadHandler(
      filename = function() {
        paste0(input$out_step10_ribo_name, "-temp-PCA-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step10_ribo_plot()$pca_plot, filename = file,
               width = input$out_step10_ribo_width, height = input$out_step10_ribo_height)
      }
    )
    
    output$save_step10_ribo_ecdf <- downloadHandler(
      filename = function() {
        paste0(input$out_step10_ribo_name, "-eCDF-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step10_ribo_plot()$ecdf_plot, filename = file,
               width = input$out_step10_ribo_width, height = input$out_step10_ribo_height)
      }
    )
    
  })
  
  
  ## 10.7 calculate the TE abundance #################
  
  step10_te_norm <- reactive({
    
    req(input$act_step10_import_te_calc)
    
    if (is.null(input$in_step10_te_rna)) {return(NULL)}
    if (is.null(input$in_step10_te_ribo)) {return(NULL)}
    
    # browser()
    source("R/fill_empty.R")
    source("R/calc_TE.R")
    
    # browser()
    
    # import the RNA-seq and Ribo-seq abundance
    suffix <- c(".genes.results", ".isoforms.results", "_cds_rpf", "_cds_rpm", "_cds_rpkm", "_cds_tpm", "Aligned.sortedByCoord.out.bam")
    ftcount_column <- c("Chr", "Start", "End", "Strand", "Length")
    
    rna_norm <- fread(input$in_step10_te_rna$datapath, sep = '\t', check.names = FALSE,
                      header="auto", stringsAsFactors = FALSE) %>% 
      dplyr::rename_with(~ gsub(paste(suffix, collapse = '|'), '', .), everything()) %>%
      dplyr::rename(Gene = colnames(.)[1]) %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::mutate(across(where(is.numeric), as.double)) %>% 
      dplyr::filter(rowSums(.) > input$in_step10_te_filter)
    
    ribo_norm <- fread(input$in_step10_te_ribo$datapath, sep = '\t', check.names = FALSE,
                       header="auto", stringsAsFactors = FALSE) %>% 
      dplyr::rename_with(~ gsub(paste(suffix, collapse = '|'), '', .), everything()) %>%
      dplyr::rename(Gene = colnames(.)[1]) %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::mutate(across(where(is.numeric), as.double)) %>%
      dplyr::filter(rowSums(.) > input$in_step10_te_filter)
    
    rna_norm <- rna_norm %>%
      dplyr::rename_all(~basename(names(rna_norm)))
    ribo_norm <- ribo_norm %>%
      dplyr::rename_all(~basename(names(ribo_norm)))
    
    rna_norm_fill0 <- fill_empty(expr_df = rna_norm, fill = input$in_step10_te_fill)
    ribo_norm_fill0 <- fill_empty(expr_df = ribo_norm, fill = input$in_step10_te_fill)
    
    # browser()
    # calculate the TE abundance
    if (import_design_clicked()) {
      
      # calculate TE
      te_res <- calc_te(rna = rna_norm_fill0, 
                        ribo = ribo_norm_fill0, 
                        phenodata = step1_design(), 
                        linked = input$in_step10_te_linked,
                        clean = input$in_step10_te_clean)
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    # browser()
    return(list(te_norm = te_res$te, phenodata = te_res$phenodata))
  })
  
  ## output the table
  observeEvent(input$act_step10_import_te_calc, {
    # browser()
    output$out_step10_te_norm <- DT::renderDataTable({
      DT::datatable(step10_te_norm()$te_norm, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "TE", page = "all"), 
                                     list(extend = 'excel', filename = "TE", page = "all"))
                    ))
    })
    
    # save the TE table
    output$save_step10_te_norm <- downloadHandler(
      filename = function() {
        paste0(input$out_step10_te_name, "-normalized-", Sys.Date(), ".txt")
      },
      content = function(file) {
        write.table(x = step10_te_norm()$te_norm, file = file, row.names = FALSE, sep = '\t', quote = F)
      }
    )
    
    # output$save_step10_te_norm <- downloadHandler(
    #   filename = function() {
    #     paste0(input$out_step10_te_name, "-normalized-", Sys.Date(), input$out_step10_te_format)
    #   },
    #   content = function(file) {
    #     if (input$out_step10_te_format == ".txt") {
    #       write.table(x = step10_te_norm()$te_norm, file = file, row.names = FALSE, sep = '\t', quote = F)
    #     } else if (input$out_step10_te_format == ".csv") {
    #       write.xlsx(x = step10_te_norm()$te_norm, file = file, row.names = FALSE)
    #     } else if (input$out_step10_te_format == ".xlsx") {
    #       write.xlsx(x = step10_te_norm()$te_norm, file = file, rowNames = FALSE)
    #     }
    #   }
    # )
    
  })
  
  
  ## 10.8 draw the TE abundance #################
  
  step10_te_plot <- reactive({
    
    req(input$act_step10_draw_te_expr)
    
    if (is.null(step10_te_norm())) {return(NULL)}
    
    source("R/draw_PCA.R")
    source("R/draw_eCDF.R")
    
    group <- input$in_step10_te_group
    fill_color <- input$in_step10_te_fill_color
    fill_alpha <- input$in_step10_te_fill_alpha
    profile <- "TE"
    line_color <- input$in_step10_te_line_color
    line_width <- input$in_step10_te_line_width
    facet <- input$in_step10_te_wrap
    log2_transform <- input$in_step10_te_log2
    xrange <- input$in_step10_te_xlim
    font_size <- input$in_step10_te_font_size
    
    # browser()
    phenodata <- step10_te_norm()$phenodata %>% 
      dplyr::filter(SeqType == "TE")
    
    te_norm <- step10_te_norm()$te_norm %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::select(any_of(phenodata$Sample)) %>% 
      na.omit()
    
    # draw the PCA plot
    te_pca_res <- draw_pca(in_data = te_norm, 
                           phenodata = phenodata, 
                           font_size = font_size,
                           group = group, 
                           title = paste(profile, "(Normalized)"),
                           fill_color = fill_color, 
                           fill_alpha = fill_alpha)
    
    # browser()
    # draw the eCDF
    te_longer <- te_norm %>% 
      tidyr::pivot_longer(cols = phenodata$Sample, names_to = "Sample", values_to = "TE") %>%
      dplyr::left_join(phenodata, by = "Sample")
    
    # log2 transform
    if (log2_transform) {
      te_longer <- te_longer %>% 
        dplyr::mutate(TE = log2(TE)) %>% 
        dplyr::filter(TE >= xrange[1], TE <= xrange[2])
    }
    
    te_ecdf_plot <- draw_ecdf(in_data = te_longer, 
                              x = "TE",
                              group = group, 
                              title = paste(profile, "(Normalized)"),
                              xlabel = "log2(TE)",
                              facet = facet, 
                              font_size = font_size,
                              line_color = line_color, 
                              line_width = line_width)
    
    return(list(pca_plot = te_pca_res$pca_scatter,
                ecdf_plot = te_ecdf_plot))
  })
  
  ## output the table
  observeEvent(input$act_step10_draw_te_expr, {
    # browser()
    ## show the expression level plot
    output$out_step10_te_pca <- renderPlot(
      width = input$out_step10_te_width * 100,
      height = input$out_step10_te_height * 100,
      {step10_te_plot()$pca_plot})
    
    output$out_step10_te_ecdf <- renderPlot(
      width = input$out_step10_te_width * 100,
      height = input$out_step10_te_height * 100,
      {step10_te_plot()$ecdf_plot})
    
    # save the expression level plot
    output$save_step10_te_pca <- downloadHandler(
      filename = function() {
        paste0(input$out_step10_te_name, "-temp-PCA-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step10_te_plot()$pca_plot, filename = file,
               width = input$out_step10_te_width, height = input$out_step10_te_height)
      }
    )
    
    output$save_step10_te_ecdf <- downloadHandler(
      filename = function() {
        paste0(input$out_step10_te_name, "-eCDF-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step10_te_plot()$ecdf_plot, filename = file,
               width = input$out_step10_te_width, height = input$out_step10_te_height)
      }
    )
    
  })
  
  ############################################################
  ## 11 draw the heatmap #################
  
  ### 11.1 import the RNA-seq table #################
  step11_rna_norm <- reactive({
    
    req(input$act_step11_import_rna_norm)
    
    if (is.null(input$in_step11_rnaseq)) {return(NULL)}
    
    # import the raw data table
    suffix <- c(".genes.results", ".isoforms.results", "_cds_rpf", "_cds_rpm", "_cds_rpkm", "_cds_tpm", "Aligned.sortedByCoord.out.bam")
    ftcount_column <- c("Chr", "Start", "End", "Strand", "Length")
    
    rna_norm_table <- fread(input$in_step11_rnaseq$datapath, sep = '\t', check.names = FALSE,
                            header="auto", stringsAsFactors = FALSE) %>% 
      dplyr::rename_with(~ gsub(paste(suffix, collapse = '|'), '', .), everything()) %>%
      dplyr::rename(Gene = colnames(.)[1]) %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::mutate(across(where(is.numeric), as.double))
    
    rna_norm_table <- rna_norm_table %>%
      dplyr::rename_all(~basename(names(rna_norm_table)))
    
    # browser()
    
    # select the samples
    if (import_design_clicked()) {
      # group the samples
      rna_norm_table <- rna_norm_table %>% 
        dplyr::select(any_of(step1_design()$Sample)) %>% 
        tibble::rownames_to_column("Gene") %>% 
        tidyr::pivot_longer(cols = -Gene, 
                            names_to = "Sample",
                            values_to = "Norm") %>%
        dplyr::left_join(step1_design(), by = "Sample") %>% 
        dplyr::group_by(Gene, !!sym(input$in_step11_rna_group)) %>%
        dplyr::reframe(Norm = mean(Norm)) %>%
        tidyr::pivot_wider(names_from = !!sym(input$in_step11_rna_group),
                           values_from = Norm) %>% 
        tibble::column_to_rownames("Gene")
      
      groups_order <- step1_design()[, input$in_step11_rna_group]
      
      rna_norm_table <- rna_norm_table %>% 
        dplyr::select(any_of(c("Gene", groups_order)))
      
      outmess <- NULL
    } else {
      outmess <- "Step 1 Design has not been imported."
    }
    
    # browser()
    
    # filter the gene average
    rna_norm_table <- rna_norm_table %>% 
      dplyr::filter(rowSums(.) > input$in_step11_rna_rowsum)
    
    # filter the samples
    rna_norm_table <- rna_norm_table %>% 
      dplyr::mutate(Effective_Sample = rowSums(. > 0)) %>% 
      dplyr::mutate(Effective_Ratio = Effective_Sample / (ncol(.) - 1) * 100) %>% 
      dplyr::filter(Effective_Ratio >= input$in_step11_rna_colsum) %>% 
      dplyr::select(-c(Effective_Sample, Effective_Ratio))
    
    # filter the specified genes
    if (input$in_step11_rna_gene_name != "") {
      heat_gene_name <- gsub(' ', ',', input$in_step11_rna_gene_name)
      heat_gene_name <- unlist(strsplit(heat_gene_name, ","))
      
      rna_norm_table <- rna_norm_table %>% 
        dplyr::filter(rownames(.) %in% heat_gene_name)
    }
    
    # filter the top expression genes
    if (input$in_step11_rna_top > 0) {
      rna_norm_table <- rna_norm_table %>% 
        dplyr::mutate(rowsum = rowSums(.)) %>% 
        dplyr::arrange(desc(rowsum)) %>% 
        head(input$in_step11_rna_top) %>%
        dplyr::select(-rowsum)
    }
    
    # calculate the correlation
    rna_norm_corr <- cor(rna_norm_table)
    
    # convert the rowname to column
    rna_norm_table <- rna_norm_table %>% 
      tibble::rownames_to_column("Gene") %>% 
      dplyr::arrange(Gene)
    
    return(list(rna_norm_table = rna_norm_table, 
                rna_norm_corr = rna_norm_corr,
                outmess = outmess))
  })
  
  ## output the table
  observeEvent(input$act_step11_import_rna_norm, {
    # browser()
    
    ## output the table
    output$out_step11_rna_norm <- DT::renderDataTable(server = T, {
      DT::datatable(step11_rna_norm()$rna_norm_table, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "RNA_seq_filter"), 
                                     list(extend = 'excel', filename = "RNA_seq_filter"))
                    ))
    })
    
    output$out_step11_rna_corr_tab <- DT::renderDataTable(server = T, {
      DT::datatable(step11_rna_norm()$rna_norm_corr, rownames = TRUE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "RNA_seq_corr"), 
                                     list(extend = 'excel', filename = "RNA_seq_corr"))
                    ))
    })
    
    output$save_step11_rna_corr <- downloadHandler(
      filename = function() {
        paste0(input$out_step11_rna_name, "-correlation-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step11_rna_norm()$rna_norm_corr, file = file, rowNames = FALSE)
      }
    )
    
  })
  
  ### 11.2 draw the RNA-seq correlation heatmap #################
  step11_rna_corr_plot <- reactive({
    
    req(input$act_step11_draw_rna_corr)
    
    if (is.null(step11_rna_norm())) {return(NULL)}
    
    my_color <- rev(colorRampPalette(RColorBrewer::brewer.pal(8, input$in_step11_rna_corr_fill_color))(100))
    my_color <- alpha(my_color, input$in_step11_rna_corr_fill_alpha)
    
    my_breaks <- seq(input$in_step11_rna_corr_min, input$in_step11_rna_corr_max,
                     (input$in_step11_rna_corr_max - input$in_step11_rna_corr_min) / 100)
    
    # browser()
    # draw the correlation heatmap
    corr_heat_plot <- pheatmap(mat = step11_rna_norm()$rna_norm_corr %>% as.matrix(),
                               scale = "none",
                               cluster_cols = F,
                               cluster_rows = F,
                               show_rownames = T,
                               show_colnames = T,
                               angle_col = input$in_step11_rna_corr_angle,
                               fontsize = input$in_step11_rna_corr_font_size,
                               
                               gaps_row = unlist(strsplit(input$in_step11_rna_corr_row_gap, ",")) %>% as.numeric(),
                               gaps_col = unlist(strsplit(input$in_step11_rna_corr_col_gap, ",")) %>% as.numeric(),
                               
                               breaks = my_breaks,
                               color = my_color,
                               border_color = input$in_step11_rna_corr_border_color,
                               na_col = "grey",
                               
                               display_numbers = input$in_step11_rna_corr_number,
                               number_format = paste0("%.", input$in_step11_rna_corr_num_format, "f"),
                               fontsize_number = input$in_step11_rna_corr_number_size,
                               number_color = input$in_step11_rna_corr_num_color,
                               
                               xlab = "",
                               ylab = "",
                               main = "RNA-seq (correlation)",
                               silent = T)
    
    
    # browser()
    return(corr_heat_plot)
  })
  
  ## output the table
  observeEvent(input$act_step11_draw_rna_corr, {
    
    ## show the expression level plot
    output$out_step11_rna_corr_plot <- renderPlot(
      width = input$out_step11_rna_corr_width * 80,
      height = input$out_step11_rna_corr_height * 80,
      {step11_rna_corr_plot()})
    
    # save the expression level plot
    output$save_step11_rna_corr_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step11_rna_corr_name, "-correlation-heat-map-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step11_rna_corr_plot() %>% as.ggplot(), filename = file,
               width = input$out_step11_rna_corr_width, height = input$out_step11_rna_corr_height)
      }
    )
    
  })
  
  
  ### 11.3 draw the RNA-seq cluster heatmap #################
  step11_rna_plot <- reactive({
    source("R/draw_heatmap.R")
    
    req(input$act_step11_draw_rna_expr)
    
    if (is.null(step11_rna_norm())) {return(NULL)}
    
    fill_color <- input$in_step11_rna_fill_color
    fill_alpha <- input$in_step11_rna_fill_alpha
    border_color <- input$in_step11_rna_border_color
    
    cluster_method <- input$in_step11_rna_method
    cluster_distance <- input$in_step11_rna_distance
    
    heat_gap_row <- unlist(strsplit(input$in_step11_rna_row_gap, ",")) %>% as.numeric()
    heat_gap_col <- unlist(strsplit(input$in_step11_rna_col_gap, ",")) %>% as.numeric()
    
    row_clust <- input$in_step11_rna_row_clust
    col_clust <- input$in_step11_rna_col_clust
    
    log2_transform <- input$in_step11_rna_log2
    
    show_rownames <- input$in_step11_rna_rowname
    show_colnames <- input$in_step11_rna_colname
    
    font_size <- input$in_step11_rna_font_size
    
    # browser()
    # log2 transform
    if (log2_transform) {
      rna_norm <- step11_rna_norm()$rna_norm_table %>% 
        tibble::column_to_rownames("Gene") %>%
        # convert the norm to log2(norm + 1)
        dplyr::mutate(across(where(is.numeric), ~log2(. + 1)))
      
    } else {
      rna_norm <- step11_rna_norm()$rna_norm_table %>% 
        tibble::column_to_rownames("Gene")
    }
    
    scale_heat_plot <- draw_heat(in_dat = rna_norm, 
                                 color = fill_color, alpha = fill_alpha, border_color = border_color,
                                 scale = "row",
                                 gaps_row = heat_gap_row, gaps_col = heat_gap_col,
                                 cluster_cols = col_clust, cluster_rows = row_clust,
                                 show_rownames = show_rownames, show_colnames = show_colnames,
                                 cluster_distance = cluster_distance, cluster_method = cluster_method,
                                 angle_col = 90, silent = T, 
                                 fontsize = font_size, title = "RNA-seq (scaled)")
    
    unscale_heat_plot <- draw_heat(in_dat = rna_norm,
                                   color = fill_color, alpha = fill_alpha, border_color = border_color,
                                   scale = "none",
                                   gaps_row = heat_gap_row, gaps_col = heat_gap_col,
                                   cluster_cols = col_clust, cluster_rows = row_clust,
                                   show_rownames = show_rownames, show_colnames = show_colnames,
                                   cluster_distance = cluster_distance, cluster_method = cluster_method,
                                   angle_col = 90, silent = T, 
                                   fontsize = font_size, title = "RNA-seq (un-scaled)")
    
    # browser()
    return(list(scale_heat_plot = scale_heat_plot,
                unscale_heat_plot = unscale_heat_plot))
  })
  
  ## output the table
  observeEvent(input$act_step11_draw_rna_expr, {
    # browser()
    ## show the expression level plot
    output$out_step11_rna_scale <- renderPlot(
      width = input$out_step11_rna_width * 100,
      height = input$out_step11_rna_height * 100,
      {step11_rna_plot()$scale_heat_plot})
    
    output$out_step11_rna_unscale <- renderPlot(
      width = input$out_step11_rna_width * 100,
      height = input$out_step11_rna_height * 100,
      {step11_rna_plot()$unscale_heat_plot})
    
    # save the expression level plot
    output$save_step11_rna_scale <- downloadHandler(
      filename = function() {
        paste0(input$out_step11_rna_name, "-scaled-heat-map-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step11_rna_plot()$scale_heat_plot, filename = file,
               width = input$out_step11_rna_width, height = input$out_step11_rna_height)
      }
    )
    
    output$save_step11_rna_unscale <- downloadHandler(
      filename = function() {
        paste0(input$out_step11_rna_name, "-un-scaled-heat-map-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step11_rna_plot()$unscale_heat_plot, filename = file,
               width = input$out_step11_rna_width, height = input$out_step11_rna_height)
      }
    )
    
  })
  
  
  ### 11.4 import the Ribo-seq table #################
  step11_ribo_norm <- reactive({
    
    req(input$act_step11_import_ribo_norm)
    
    if (is.null(input$in_step11_riboseq)) {return(NULL)}
    
    # import the raw data table
    suffix <- c(".genes.results", ".isoforms.results", "_cds_rpf", "_cds_rpm", "_cds_rpkm", "_cds_tpm", "Aligned.sortedByCoord.out.bam")
    ftcount_column <- c("Chr", "Start", "End", "Strand", "Length")
    
    ribo_norm_table <- fread(input$in_step11_riboseq$datapath, sep = '\t', check.names = FALSE,
                             header="auto", stringsAsFactors = FALSE) %>% 
      dplyr::rename_with(~ gsub(paste(suffix, collapse = '|'), '', .), everything()) %>%
      dplyr::rename(Gene = colnames(.)[1]) %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::mutate(across(where(is.numeric), as.double))
    
    ribo_norm_table <- ribo_norm_table %>%
      dplyr::rename_all(~basename(names(ribo_norm_table)))
    
    # browser()
    
    # select the samples
    if (import_design_clicked()) {
      # group the samples
      ribo_norm_table <- ribo_norm_table %>% 
        dplyr::select(any_of(step1_design()$Sample)) %>% 
        tibble::rownames_to_column("Gene") %>% 
        tidyr::pivot_longer(cols = -Gene, 
                            names_to = "Sample",
                            values_to = "Norm") %>%
        dplyr::left_join(step1_design(), by = "Sample") %>% 
        dplyr::group_by(Gene, !!sym(input$in_step11_ribo_group)) %>%
        dplyr::reframe(Norm = mean(Norm)) %>%
        tidyr::pivot_wider(names_from = !!sym(input$in_step11_ribo_group),
                           values_from = Norm) %>% 
        tibble::column_to_rownames("Gene")
      
      groups_order <- step1_design()[, input$in_step11_ribo_group]
      
      ribo_norm_table <- ribo_norm_table %>% 
        dplyr::select(any_of(c("Gene", groups_order)))
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    # browser()
    
    # filter the gene average
    ribo_norm_table <- ribo_norm_table %>% 
      dplyr::filter(rowSums(.) > input$in_step11_ribo_rowsum)
    
    # filter the samples
    ribo_norm_table <- ribo_norm_table %>% 
      dplyr::mutate(Effective_Sample = rowSums(. > 0)) %>% 
      dplyr::mutate(Effective_Ratio = Effective_Sample / (ncol(.) - 1) * 100) %>% 
      dplyr::filter(Effective_Ratio >= input$in_step11_ribo_colsum) %>% 
      dplyr::select(-c(Effective_Sample, Effective_Ratio))
    
    # filter the specified genes
    if (input$in_step11_ribo_gene_name != "") {
      heat_gene_name <- gsub(' ', ',', input$in_step11_ribo_gene_name)
      heat_gene_name <- unlist(strsplit(heat_gene_name, ","))
      
      ribo_norm_table <- ribo_norm_table %>% 
        dplyr::filter(rownames(.) %in% heat_gene_name)
    }
    
    # filter the top expression genes
    if (input$in_step11_ribo_top > 0) {
      ribo_norm_table <- ribo_norm_table %>% 
        dplyr::mutate(rowsum = rowSums(.)) %>% 
        dplyr::arrange(desc(rowsum)) %>% 
        head(input$in_step11_ribo_top) %>%
        dplyr::select(-rowsum)
    }
    
    # calculate the correlation
    ribo_norm_corr <- cor(ribo_norm_table)
    
    # convert the rowname to column
    ribo_norm_table <- ribo_norm_table %>% 
      tibble::rownames_to_column("Gene") %>% 
      dplyr::arrange(Gene)
    
    
    return(list(ribo_norm_table = ribo_norm_table, 
                ribo_norm_corr = ribo_norm_corr))
  })
  
  ## output the table
  observeEvent(input$act_step11_import_ribo_norm, {
    # browser()
    
    ## output the table
    output$out_step11_ribo_norm <- DT::renderDataTable(server = T, {
      DT::datatable(step11_ribo_norm()$ribo_norm_table, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "Ribo_seq_filter"), 
                                     list(extend = 'excel', filename = "Ribo_seq_filter"))
                    ))
    })
    
    output$out_step11_ribo_corr_tab <- DT::renderDataTable(server = T, {
      DT::datatable(step11_ribo_norm()$ribo_norm_corr, rownames = TRUE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "Ribo_seq_corr"), 
                                     list(extend = 'excel', filename = "Ribo_seq_corr"))
                    ))
    })
    
    output$save_step11_ribo_corr <- downloadHandler(
      filename = function() {
        paste0(input$out_step11_ribo_name, "-correlation-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step11_ribo_norm()$ribo_norm_corr, file = file, rowNames = FALSE)
      }
    )
    
  })
  
  
  ### 11.5 draw the Ribo-seq correlation heatmap #################
  step11_ribo_corr_plot <- reactive({
    
    req(input$act_step11_draw_ribo_corr)
    
    if (is.null(step11_ribo_norm())) {return(NULL)}
    
    my_color <- rev(colorRampPalette(RColorBrewer::brewer.pal(8, input$in_step11_ribo_corr_fill_color))(100))
    my_color <- alpha(my_color, input$in_step11_ribo_corr_fill_alpha)
    
    my_breaks <- seq(input$in_step11_ribo_corr_min, input$in_step11_ribo_corr_max,
                     (input$in_step11_ribo_corr_max - input$in_step11_ribo_corr_min) / 100)
    
    # browser()
    # draw the correlation heatmap
    corr_heat_plot <- pheatmap(mat = step11_ribo_norm()$ribo_norm_corr %>% as.matrix(),
                               scale = "none",
                               cluster_cols = F,
                               cluster_rows = F,
                               show_rownames = T,
                               show_colnames = T,
                               angle_col = input$in_step11_ribo_corr_angle,
                               fontsize = input$in_step11_ribo_corr_font_size,
                               
                               gaps_row = unlist(strsplit(input$in_step11_ribo_corr_row_gap, ",")) %>% as.numeric(),
                               gaps_col = unlist(strsplit(input$in_step11_ribo_corr_col_gap, ",")) %>% as.numeric(),
                               
                               breaks = my_breaks,
                               color = my_color,
                               border_color = input$in_step11_ribo_corr_border_color,
                               na_col = "grey",
                               
                               display_numbers = input$in_step11_ribo_corr_number,
                               number_format = paste0("%.", input$in_step11_ribo_corr_num_format, "f"),
                               fontsize_number = input$in_step11_ribo_corr_number_size,
                               number_color = input$in_step11_ribo_corr_num_color,
                               
                               xlab = "",
                               ylab = "",
                               main = "Ribo-seq (correlation)",
                               silent = T)
    
    # browser()
    return(corr_heat_plot)
  })
  
  ## output the table
  observeEvent(input$act_step11_draw_ribo_corr, {
    
    ## show the expression level plot
    output$out_step11_ribo_corr_plot <- renderPlot(
      width = input$out_step11_ribo_corr_width * 80,
      height = input$out_step11_ribo_corr_height * 80,
      {step11_ribo_corr_plot()})
    
    # save the expression level plot
    output$save_step11_ribo_corr_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step11_ribo_corr_name, "-correlation-heat-map-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step11_ribo_corr_plot() %>% as.ggplot(), filename = file,
               width = input$out_step11_ribo_corr_width, height = input$out_step11_ribo_corr_height)
      }
    )
    
  })
  
  
  ### 11.6 draw the Ribo-seq heatmap #################
  step11_ribo_plot <- reactive({
    source("R/draw_heatmap.R")
    
    req(input$act_step11_draw_ribo_expr)
    
    if (is.null(step11_ribo_norm())) {return(NULL)}
    
    fill_color <- input$in_step11_ribo_fill_color
    fill_alpha <- input$in_step11_ribo_fill_alpha
    border_color <- input$in_step11_ribo_border_color
    
    cluster_method <- input$in_step11_ribo_method
    cluster_distance <- input$in_step11_ribo_distance
    
    heat_gap_row <- unlist(strsplit(input$in_step11_ribo_row_gap, ",")) %>% as.numeric()
    heat_gap_col <- unlist(strsplit(input$in_step11_ribo_col_gap, ",")) %>% as.numeric()
    
    row_clust <- input$in_step11_ribo_row_clust
    col_clust <- input$in_step11_ribo_col_clust
    
    log2_transform <- input$in_step11_ribo_log2
    
    show_rownames <- input$in_step11_ribo_rowname
    show_colnames <- input$in_step11_ribo_colname
    
    font_size <- input$in_step11_ribo_font_size
    
    # browser()
    # log2 transform
    if (log2_transform) {
      ribo_norm <- step11_ribo_norm()$ribo_norm_table %>% 
        tibble::column_to_rownames("Gene") %>%
        # convert the norm to log2(norm + 1)
        dplyr::mutate(across(where(is.numeric), ~log2(. + 1)))
      
    } else {
      ribo_norm <- step11_ribo_norm()$ribo_norm_table %>% 
        tibble::column_to_rownames("Gene")
    }
    
    scale_heat_plot <- draw_heat(in_dat = ribo_norm, 
                                 color = fill_color, alpha = fill_alpha, border_color = border_color,
                                 scale = "row",
                                 gaps_row = heat_gap_row, gaps_col = heat_gap_col,
                                 cluster_cols = col_clust, cluster_rows = row_clust,
                                 show_rownames = show_rownames, show_colnames = show_colnames,
                                 cluster_distance = cluster_distance, cluster_method = cluster_method,
                                 angle_col = 90, silent = T, 
                                 fontsize = font_size, title = "Ribo-seq (scaled)")
    
    unscale_heat_plot <- draw_heat(in_dat = ribo_norm,
                                   color = fill_color, alpha = fill_alpha, border_color = border_color,
                                   scale = "none",
                                   gaps_row = heat_gap_row, gaps_col = heat_gap_col,
                                   cluster_cols = col_clust, cluster_rows = row_clust,
                                   show_rownames = show_rownames, show_colnames = show_colnames,
                                   cluster_distance = cluster_distance, cluster_method = cluster_method,
                                   angle_col = 90, silent = T, 
                                   fontsize = font_size, title = "Ribo-seq (un-scaled)")
    
    return(list(scale_heat_plot = scale_heat_plot,
                unscale_heat_plot = unscale_heat_plot))
  })
  
  ## output the table
  observeEvent(input$act_step11_draw_ribo_expr, {
    # browser()
    ## show the expression level plot
    output$out_step11_ribo_scale <- renderPlot(
      width = input$out_step11_ribo_width * 100,
      height = input$out_step11_ribo_height * 100,
      {step11_ribo_plot()$scale_heat_plot})
    
    output$out_step11_ribo_unscale <- renderPlot(
      width = input$out_step11_ribo_width * 100,
      height = input$out_step11_ribo_height * 100,
      {step11_ribo_plot()$unscale_heat_plot})
    
    # save the expression level plot
    output$save_step11_ribo_scale <- downloadHandler(
      filename = function() {
        paste0(input$out_step11_ribo_name, "-scaled-heat-map-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step11_ribo_plot()$scale_heat_plot, filename = file,
               width = input$out_step11_ribo_width, height = input$out_step11_ribo_height)
      }
    )
    
    output$save_step11_ribo_unscale <- downloadHandler(
      filename = function() {
        paste0(input$out_step11_ribo_name, "-un-scaled-heat-map-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step11_ribo_plot()$unscale_heat_plot, filename = file,
               width = input$out_step11_ribo_width, height = input$out_step11_ribo_height)
      }
    )
    
  })
  
  
  ### 11.7 import the TE table #################
  step11_te_norm <- reactive({
    
    req(input$act_step11_import_te_norm)
    
    if (is.null(input$in_step11_te)) {return(NULL)}
    
    # import the raw data table
    suffix <- c(".genes.results", ".isoforms.results", "_cds_rpf", "_cds_rpm", "_cds_rpkm", "_cds_tpm", "Aligned.sortedByCoord.out.bam")
    ftcount_column <- c("Chr", "Start", "End", "Strand", "Length")
    
    te_norm_table <- fread(input$in_step11_te$datapath, sep = '\t', check.names = FALSE,
                           header="auto", stringsAsFactors = FALSE) %>% 
      dplyr::rename_with(~ gsub(paste(suffix, collapse = '|'), '', .), everything()) %>%
      dplyr::rename(Gene = colnames(.)[1]) %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::mutate(across(where(is.numeric), as.double))
    
    te_norm_table <- te_norm_table %>%
      dplyr::rename_all(~basename(names(te_norm_table)))
    
    # browser()
    
    # select the samples
    if (import_design_clicked()) {
      # group the samples
      te_norm_table <- te_norm_table %>% 
        dplyr::select(any_of(step1_design()$Sample)) %>% 
        tibble::rownames_to_column("Gene") %>% 
        tidyr::pivot_longer(cols = -Gene, 
                            names_to = "Sample",
                            values_to = "Norm") %>%
        dplyr::left_join(step1_design(), by = "Sample") %>% 
        dplyr::group_by(Gene, !!sym(input$in_step11_te_group)) %>%
        dplyr::reframe(Norm = mean(Norm)) %>%
        tidyr::pivot_wider(names_from = !!sym(input$in_step11_te_group),
                           values_from = Norm) %>% 
        tibble::column_to_rownames("Gene")
      
      groups_order <- step1_design()[, input$in_step11_te_group]
      
      te_norm_table <- te_norm_table %>% 
        dplyr::select(any_of(c("Gene", groups_order)))
      
    } else {
      print("Step 1 Design has not been imported.")
    }
    
    # browser()
    
    # filter the gene average
    te_norm_table <- te_norm_table %>% 
      dplyr::filter(rowSums(.) > input$in_step11_te_rowsum)
    
    # filter the samples
    te_norm_table <- te_norm_table %>% 
      dplyr::mutate(Effective_Sample = rowSums(. > 0)) %>% 
      dplyr::mutate(Effective_Ratio = Effective_Sample / (ncol(.) - 1) * 100) %>% 
      dplyr::filter(Effective_Ratio >= input$in_step11_te_colsum) %>% 
      dplyr::select(-c(Effective_Sample, Effective_Ratio))
    
    # filter the specified genes
    if (input$in_step11_te_gene_name != "") {
      heat_gene_name <- gsub(' ', ',', input$in_step11_te_gene_name)
      heat_gene_name <- unlist(strsplit(heat_gene_name, ","))
      
      te_norm_table <- te_norm_table %>% 
        dplyr::filter(rownames(.) %in% heat_gene_name)
    }
    
    # filter the top expression genes
    if (input$in_step11_te_top > 0) {
      te_norm_table <- te_norm_table %>% 
        dplyr::mutate(rowsum = rowSums(.)) %>% 
        dplyr::arrange(desc(rowsum)) %>% 
        head(input$in_step11_te_top) %>%
        dplyr::select(-rowsum)
    }
    
    # calculate the correlation
    te_norm_corr <- cor(te_norm_table)
    
    # convert the rowname to column
    te_norm_table <- te_norm_table %>% 
      tibble::rownames_to_column("Gene") %>% 
      dplyr::arrange(Gene)
    
    
    return(list(te_norm_table = te_norm_table, 
                te_norm_corr = te_norm_corr))
  })
  
  ## output the table
  observeEvent(input$act_step11_import_te_norm, {
    # browser()
    
    ## output the table
    output$out_step11_te_norm <- DT::renderDataTable(server = T, {
      DT::datatable(step11_te_norm()$te_norm_table, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "TE_filter"), 
                                     list(extend = 'excel', filename = "TE_filter"))
                    ))
    })
    
    output$out_step11_te_corr_tab <- DT::renderDataTable(server = T, {
      DT::datatable(step11_te_norm()$te_norm_corr, rownames = TRUE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "TE_corr"), 
                                     list(extend = 'excel', filename = "TE_corr"))
                    ))
    })
    
    # download the correlation
    output$save_step11_te_corr <- downloadHandler(
      filename = function() {
        paste0(input$out_step11_te_name, "-correlation-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step11_te_norm()$te_norm_corr, file = file, rowNames = FALSE)
      }
    )
    
  })
  
  
  ### 11.8 draw the TE correlation heatmap #################
  step11_te_corr_plot <- reactive({
    
    req(input$act_step11_draw_te_corr)
    
    if (is.null(step11_te_norm())) {return(NULL)}
    
    my_color <- rev(colorRampPalette(RColorBrewer::brewer.pal(8, input$in_step11_te_corr_fill_color))(100))
    my_color <- alpha(my_color, input$in_step11_te_corr_fill_alpha)
    
    my_breaks <- seq(input$in_step11_te_corr_min, input$in_step11_te_corr_max,
                     (input$in_step11_te_corr_max - input$in_step11_te_corr_min) / 100)
    
    # browser()
    # draw the correlation heatmap
    corr_heat_plot <- pheatmap(mat = step11_te_norm()$te_norm_corr %>% as.matrix(),
                               scale = "none",
                               cluster_cols = F,
                               cluster_rows = F,
                               show_rownames = T,
                               show_colnames = T,
                               angle_col = input$in_step11_te_corr_angle,
                               fontsize = input$in_step11_te_corr_font_size,
                               
                               gaps_row = unlist(strsplit(input$in_step11_te_corr_row_gap, ",")) %>% as.numeric(),
                               gaps_col = unlist(strsplit(input$in_step11_te_corr_col_gap, ",")) %>% as.numeric(),
                               
                               breaks = my_breaks,
                               color = my_color,
                               border_color = input$in_step11_te_corr_border_color,
                               na_col = "grey",
                               
                               display_numbers = input$in_step11_te_corr_number,
                               number_format = paste0("%.", input$in_step11_te_corr_num_format, "f"),
                               fontsize_number = input$in_step11_te_corr_number_size,
                               number_color = input$in_step11_te_corr_num_color,
                               
                               xlab = "",
                               ylab = "",
                               main = "TE (correlation)",
                               silent = T)
    
    return(corr_heat_plot)
  })
  
  ## output the table
  observeEvent(input$act_step11_draw_te_corr, {
    
    ## show the expression level plot
    output$out_step11_te_corr_plot <- renderPlot(
      width = input$out_step11_te_corr_width * 80,
      height = input$out_step11_te_corr_height * 80,
      {step11_te_corr_plot()})
    
    # save the expression level plot
    output$save_step11_te_corr_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step11_te_corr_name, "-correlation-heat-map-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step11_te_corr_plot() %>% as.ggplot(), filename = file,
               width = input$out_step11_te_corr_width, height = input$out_step11_te_corr_height)
      }
    )
    
  })
  
  
  ### 11.9 draw the TE heatmap #################
  step11_te_plot <- reactive({
    source("R/draw_heatmap.R")
    
    req(input$act_step11_draw_te_expr)
    
    if (is.null(step11_te_norm())) {return(NULL)}
    
    fill_color <- input$in_step11_te_fill_color
    fill_alpha <- input$in_step11_te_fill_alpha
    border_color <- input$in_step11_te_border_color
    
    cluster_method <- input$in_step11_te_method
    cluster_distance <- input$in_step11_te_distance
    
    heat_gap_row <- unlist(strsplit(input$in_step11_te_row_gap, ",")) %>% as.numeric()
    heat_gap_col <- unlist(strsplit(input$in_step11_te_col_gap, ",")) %>% as.numeric()
    
    row_clust <- input$in_step11_te_row_clust
    col_clust <- input$in_step11_te_col_clust
    
    log2_transform <- input$in_step11_te_log2
    
    show_rownames <- input$in_step11_te_rowname
    show_colnames <- input$in_step11_te_colname
    
    font_size <- input$in_step11_te_font_size
    
    # browser()
    # log2 transform
    if (log2_transform) {
      te_norm <- step11_te_norm()$te_norm_table %>% 
        tibble::column_to_rownames("Gene") %>%
        # convert the norm to log2(norm + 1)
        dplyr::mutate(across(where(is.numeric), ~log2(. + 1)))
      
    } else {
      te_norm <- step11_te_norm()$te_norm_table %>% 
        tibble::column_to_rownames("Gene")
    }
    
    scale_heat_plot <- draw_heat(in_dat = te_norm, 
                                 color = fill_color, alpha = fill_alpha, border_color = border_color,
                                 scale = "row",
                                 gaps_row = heat_gap_row, gaps_col = heat_gap_col,
                                 cluster_cols = col_clust, cluster_rows = row_clust,
                                 show_rownames = show_rownames, show_colnames = show_colnames,
                                 cluster_distance = cluster_distance, cluster_method = cluster_method,
                                 angle_col = 90, silent = T, 
                                 fontsize = font_size, title = "TE (scaled)")
    
    unscale_heat_plot <- draw_heat(in_dat = te_norm,
                                   color = fill_color, alpha = fill_alpha, border_color = border_color,
                                   scale = "none",
                                   gaps_row = heat_gap_row, gaps_col = heat_gap_col,
                                   cluster_cols = col_clust, cluster_rows = row_clust,
                                   show_rownames = show_rownames, show_colnames = show_colnames,
                                   cluster_distance = cluster_distance, cluster_method = cluster_method,
                                   angle_col = 90, silent = T, 
                                   fontsize = font_size, title = "TE (un-scaled)")
    
    return(list(scale_heat_plot = scale_heat_plot,
                unscale_heat_plot = unscale_heat_plot))
  })
  
  ## output the table
  observeEvent(input$act_step11_draw_te_expr, {
    # browser()
    ## show the expression level plot
    output$out_step11_te_scale <- renderPlot(
      width = input$out_step11_te_width * 100,
      height = input$out_step11_te_height * 100,
      {step11_te_plot()$scale_heat_plot})
    
    output$out_step11_te_unscale <- renderPlot(
      width = input$out_step11_te_width * 100,
      height = input$out_step11_te_height * 100,
      {step11_te_plot()$unscale_heat_plot})
    
    # save the expression level plot
    output$save_step11_te_scale <- downloadHandler(
      filename = function() {
        paste0(input$out_step11_te_name, "-scaled-heat-map-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step11_te_plot()$scale_heat_plot, filename = file,
               width = input$out_step11_te_width, height = input$out_step11_te_height)
      }
    )
    
    output$save_step11_te_unscale <- downloadHandler(
      filename = function() {
        paste0(input$out_step11_te_name, "-un-scaled-heat-map-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step11_te_plot()$unscale_heat_plot, filename = file,
               width = input$out_step11_te_width, height = input$out_step11_te_height)
      }
    )
    
  })
  
  
  ############################################################
  ## 12 draw the PCA plot #################
  
  ### 12.1 import the RNA-seq exprs #######################
  step12_rna_norm <- reactive({
    
    req(input$act_step12_import_rna_norm)
    
    if (is.null(input$in_step12_rnaseq)) {return(NULL)}
    
    # import the raw data table
    suffix <- c(".genes.results", ".isoforms.results", "_cds_rpf", "_cds_rpm", "_cds_rpkm", "_cds_tpm", "Aligned.sortedByCoord.out.bam")
    ftcount_column <- c("Chr", "Start", "End", "Strand", "Length")
    
    rna_norm <- fread(input$in_step12_rnaseq$datapath, sep = '\t', check.names = FALSE,
                      header="auto", stringsAsFactors = FALSE) %>% 
      dplyr::rename_with(~ gsub(paste(suffix, collapse = '|'), '', .), everything()) %>%
      dplyr::rename(Gene = colnames(.)[1]) %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::mutate(across(where(is.numeric), as.double))
    
    rna_norm <- rna_norm %>%
      dplyr::rename_all(~basename(names(rna_norm)))
    
    # browser()
    
    # select the samples
    # if (import_design_clicked()) {
    #   # group the samples
    #   rna_norm <- rna_norm %>% 
    #     dplyr::select(any_of(step1_design()$Sample)) %>% 
    #     tibble::rownames_to_column("Gene") %>% 
    #     tidyr::pivot_longer(cols = -Gene, 
    #                         names_to = "Sample",
    #                         values_to = "Norm") %>%
    #     dplyr::left_join(step1_design(), by = "Sample") %>% 
    #     dplyr::group_by(Gene, !!sym(input$in_step12_rna_group)) %>%
    #     dplyr::reframe(Norm = mean(Norm)) %>%
    #     tidyr::pivot_wider(names_from = !!sym(input$in_step12_rna_group),
    #                        values_from = Norm) %>% 
    #     tibble::column_to_rownames("Gene")
    #   
    # } else {
    #   print("Step 1 Design has not been imported.")
    # }
    
    # browser()
    
    # filter the gene average
    rna_norm <- rna_norm %>% 
      dplyr::filter(rowSums(.) > input$in_step12_rna_rowsum)
    
    # filter the samples
    rna_norm <- rna_norm %>% 
      dplyr::mutate(Effective_Sample = rowSums(. > 0)) %>% 
      dplyr::mutate(Effective_Ratio = Effective_Sample / (ncol(.) - 1) * 100) %>% 
      dplyr::filter(Effective_Ratio >= input$in_step12_rna_colsum) %>% 
      dplyr::select(-c(Effective_Sample, Effective_Ratio))
    
    # filter the specified genes
    if (input$in_step12_rna_gene_name != "") {
      heat_gene_name <- gsub(' ', ',', input$in_step12_rna_gene_name)
      heat_gene_name <- unlist(strsplit(heat_gene_name, ","))
      
      rna_norm <- rna_norm %>% 
        dplyr::filter(rownames(.) %in% heat_gene_name)
    }
    
    # filter the top expression genes
    if (input$in_step12_rna_top > 0) {
      rna_norm <- rna_norm %>% 
        dplyr::mutate(rowsum = rowSums(.)) %>% 
        dplyr::arrange(desc(rowsum)) %>% 
        head(input$in_step12_rna_top) %>%
        dplyr::select(-rowsum)
    }
    
    # browser()
    # convert the norm to log2(norm  + 1)
    if (input$in_step12_rna_log2) {
      rna_norm <- log2(rna_norm + 1)
    }
    
    # convert the rowname to column
    rna_norm <- rna_norm %>% 
      tibble::rownames_to_column("Gene") %>% 
      dplyr::arrange(Gene)
    
    
    return(rna_norm)
  })
  
  ## output the table
  observeEvent(input$act_step12_import_rna_norm, {
    # browser()
    
    ## output the table
    output$out_step12_rna_norm <- DT::renderDataTable({
      DT::datatable(step12_rna_norm(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "RNA_seq_filter"), 
                                     list(extend = 'excel', filename = "RNA_seq_filter"))
                    ))
    })
    
  })
  
  ### 12.2 draw the RNA-seq PCA plot #################
  step12_rna_plot <- reactive({
    
    req(input$act_step12_draw_rna_pca)
    
    if (is.null(step12_rna_norm())) {return(NULL)}
    
    # browser()
    source("R/draw_PCA.R")
    
    if (import_design_clicked()) {
      # make the pheno table
      phenodata <- step1_design() %>% 
        dplyr::filter(Sample %in% colnames(step12_rna_norm())) %>%
        dplyr::filter(SeqType == "RNA") %>%
        dplyr::mutate(Group = factor(!!sym(input$in_step12_rna_group), 
                                     levels = unique(!!sym(input$in_step12_rna_group)))) %>% 
        dplyr::mutate(row_name = Sample) %>%  
        tibble::column_to_rownames(var = "row_name")
      
    } else {
      print("Step 1 Design has not been imported.")
      {return(NULL)}
    }
    
    group <- input$in_step12_rna_group
    
    dot_color <- input$in_step12_rna_dot_color
    dot_alpha <- input$in_step12_rna_dot_alpha
    
    line_color <- input$in_step12_rna_line_color
    line_width <- input$in_step12_rna_line_width
    
    xmax <- input$in_step12_rna_xmax
    xmin <- input$in_step12_rna_xmin
    
    ymax <- input$in_step12_rna_ymax
    ymin <- input$in_step12_rna_ymin
    
    gene_num <- 50
    
    dot_size <- input$in_step12_rna_dot_size
    font_size <- input$in_step12_rna_font_size
    
    label_sample <- input$in_step12_rna_label_sample
    label_size <- input$in_step12_rna_label_size
    # browser()
    
    # draw the PCA plot
    norm_pca_res <- draw_pca(in_data = step12_rna_norm() %>% tibble::column_to_rownames("Gene"), 
                             phenodata = phenodata, 
                             x_min = xmin, x_max = xmax, 
                             y_min = ymin, y_max = ymax, 
                             gene_num = gene_num,
                             label_sample = label_sample,
                             label_size = label_size,
                             font_size = font_size, 
                             dot_size = dot_size,
                             fill_color = dot_color, 
                             fill_alpha = dot_alpha, 
                             group = group, 
                             title = "RNA-seq")
    
    # browser()
    return(list(eig_plot = norm_pca_res$eig_plot,
                pca_plot = norm_pca_res$pca_scatter,
                var_bar_plot = norm_pca_res$pca_contrib_var,
                ind_bar_plot = norm_pca_res$pca_contrib_ind,
                pca_res = norm_pca_res$pca_res))
    
  })
  
  ## output the table
  observeEvent(input$act_step12_draw_rna_pca, {
    
    # browser()
    output$out_step12_rna_eigen <- DT::renderDataTable({
      DT::datatable(step12_rna_plot()$pca_res$var.coord, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "RNA_seq_filter"), 
                                     list(extend = 'excel', filename = "RNA_seq_filter"))
                    ))
    })
    
    output$out_step12_rna_var_contrib <- DT::renderDataTable({
      DT::datatable(step12_rna_plot()$pca_res$var.contrib, rownames = TRUE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "RNA_seq_filter"), 
                                     list(extend = 'excel', filename = "RNA_seq_filter"))
                    ))
    })
    
    output$out_step12_rna_ind_contrib <- DT::renderDataTable({
      DT::datatable(step12_rna_plot()$pca_res$ind.contrib, rownames = TRUE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "RNA_seq_filter"), 
                                     list(extend = 'excel', filename = "RNA_seq_filter"))
                    ))
    })
    
    output$save_step12_rna_eigen <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_rna_name, "-PCA-results-", Sys.Date(), '.xlsx')
      },
      content = function(file) {
        write.xlsx(x = step12_rna_plot()$pca_res, file = file, rowNames = TRUE)
      }
    )
    
    ## show the PCA plot
    output$out_step12_rna_eigen_plot <- renderPlot(
      width = input$out_step12_rna_width * 80,
      height = input$out_step12_rna_height * 60,
      {step12_rna_plot()$eig_plot})
    
    output$out_step12_rna_scatter <- renderPlot(
      width = input$out_step12_rna_width * 100,
      height = input$out_step12_rna_height * 80,
      {step12_rna_plot()$pca_plot})
    
    output$out_step12_rna_variances <- renderPlot(
      width = input$out_step12_rna_width * 80,
      height = input$out_step12_rna_height * 100,
      {step12_rna_plot()$var_bar_plot})
    
    output$out_step12_rna_individual <- renderPlot(
      width = input$out_step12_rna_width * 80,
      height = input$out_step12_rna_height * 100,
      {step12_rna_plot()$ind_bar_plot})
    
    # save the expression level plot
    output$save_step12_rna_eigen_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_rna_name, "-eigen-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step12_rna_plot()$eig_var_plot, filename = file,
               width = input$out_step12_rna_width, height = input$out_step12_rna_height)
      }
    )
    
    output$save_step12_rna_scatter <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_rna_name, "-PCA-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step12_rna_plot()$pca_plot, filename = file,
               width = input$out_step12_rna_width, height = input$out_step12_rna_height)
      }
    )
    
    output$save_step12_rna_variances <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_rna_name, "-variances-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step12_rna_plot()$var_bar_plot, filename = file,
               width = input$out_step12_rna_width, height = input$out_step12_rna_height)
      }
    )
    
    output$save_step12_rna_individual <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_rna_name, "-individual-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step12_rna_plot()$ind_bar_plot, filename = file,
               width = input$out_step12_rna_width, height = input$out_step12_rna_height)
      }
    )
    
  })
  
  ### 12.3 import the Ribo-seq exprs #################
  step12_ribo_norm <- reactive({
    
    req(input$act_step12_import_ribo_norm)
    
    if (is.null(input$in_step12_riboseq)) {return(NULL)}
    
    # import the raw data table
    suffix <- c(".genes.results", ".isoforms.results", "_cds_rpf", "_cds_rpm", "_cds_rpkm", "_cds_tpm", "Aligned.sortedByCoord.out.bam")
    ftcount_column <- c("Chr", "Start", "End", "Strand", "Length")
    
    ribo_norm <- fread(input$in_step12_riboseq$datapath, sep = '\t', check.names = FALSE,
                       header="auto", stringsAsFactors = FALSE) %>% 
      dplyr::rename_with(~ gsub(paste(suffix, collapse = '|'), '', .), everything()) %>%
      dplyr::rename(Gene = colnames(.)[1]) %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::mutate(across(where(is.numeric), as.double))
    
    ribo_norm <- ribo_norm %>%
      dplyr::rename_all(~basename(names(ribo_norm)))
    
    # browser()
    
    # filter the gene average
    ribo_norm <- ribo_norm %>% 
      dplyr::filter(rowSums(.) > input$in_step12_ribo_rowsum)
    
    # filter the samples
    ribo_norm <- ribo_norm %>% 
      dplyr::mutate(Effective_Sample = rowSums(. > 0)) %>% 
      dplyr::mutate(Effective_Ratio = Effective_Sample / (ncol(.) - 1) * 100) %>% 
      dplyr::filter(Effective_Ratio >= input$in_step12_ribo_colsum) %>% 
      dplyr::select(-c(Effective_Sample, Effective_Ratio))
    
    # filter the specified genes
    if (input$in_step12_ribo_gene_name != "") {
      heat_gene_name <- gsub(' ', ',', input$in_step12_ribo_gene_name)
      heat_gene_name <- unlist(strsplit(heat_gene_name, ","))
      
      ribo_norm <- ribo_norm %>% 
        dplyr::filter(rownames(.) %in% heat_gene_name)
    }
    
    # filter the top expression genes
    if (input$in_step12_ribo_top > 0) {
      ribo_norm <- ribo_norm %>% 
        dplyr::mutate(rowsum = rowSums(.)) %>% 
        dplyr::arrange(desc(rowsum)) %>% 
        head(input$in_step12_ribo_top) %>%
        dplyr::select(-rowsum)
    }
    
    # browser()
    # convert the norm to log2(norm  + 1)
    if (input$in_step12_ribo_log2) {
      ribo_norm <- log2(ribo_norm + 1)
    }
    
    # convert the rowname to column
    ribo_norm <- ribo_norm %>% 
      tibble::rownames_to_column("Gene") %>% 
      dplyr::arrange(Gene)
    
    
    return(ribo_norm)
  })
  
  ## output the table
  observeEvent(input$act_step12_import_ribo_norm, {
    # browser()
    
    ## output the table
    output$out_step12_ribo_norm <- DT::renderDataTable({
      DT::datatable(step12_ribo_norm(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "ribo_seq_filter"), 
                                     list(extend = 'excel', filename = "ribo_seq_filter"))
                    ))
    })
    
  })
  
  ### 12.4 draw the Ribo-seq PCA plot #################
  step12_ribo_plot <- reactive({
    
    req(input$act_step12_draw_ribo_pca)
    
    if (is.null(step12_ribo_norm())) {return(NULL)}
    
    # browser()
    source("R/draw_PCA.R")
    
    if (import_design_clicked()) {
      # make the pheno table
      phenodata <- step1_design() %>% 
        dplyr::filter(Sample %in% colnames(step12_ribo_norm())) %>%
        dplyr::filter(SeqType == "RIBO") %>%
        dplyr::mutate(Group = factor(!!sym(input$in_step12_ribo_group), 
                                     levels = unique(!!sym(input$in_step12_ribo_group)))) %>% 
        dplyr::mutate(row_name = Sample) %>% 
        tibble::column_to_rownames(var = "row_name")
      
    } else {
      print("Step 1 Design has not been imported.")
      {return(NULL)}
    }
    
    group <- input$in_step12_ribo_group
    
    dot_color <- input$in_step12_ribo_dot_color
    dot_alpha <- input$in_step12_ribo_dot_alpha
    
    line_color <- input$in_step12_ribo_line_color
    line_width <- input$in_step12_ribo_line_width
    
    xmax <- input$in_step12_ribo_xmax
    xmin <- input$in_step12_ribo_xmin
    
    ymax <- input$in_step12_ribo_ymax
    ymin <- input$in_step12_ribo_ymin
    
    gene_num <- 50
    
    dot_size <- input$in_step12_ribo_dot_size
    font_size <- input$in_step12_ribo_font_size
    
    label_sample <- input$in_step12_ribo_label_sample
    label_size <- input$in_step12_ribo_label_size
    # browser()
    
    # draw the PCA plot
    norm_pca_res <- draw_pca(in_data = step12_ribo_norm() %>% tibble::column_to_rownames("Gene"), 
                             phenodata = phenodata, 
                             x_min = xmin, x_max = xmax, 
                             y_min = ymin, y_max = ymax, 
                             gene_num = gene_num,
                             label_sample = label_sample,
                             label_size = label_size,
                             font_size = font_size, 
                             dot_size = dot_size,
                             fill_color = dot_color, 
                             fill_alpha = dot_alpha, 
                             group = group, 
                             title = "Ribo-seq")
    
    # browser()
    return(list(eig_plot = norm_pca_res$eig_plot,
                pca_plot = norm_pca_res$pca_scatter,
                var_bar_plot = norm_pca_res$pca_contrib_var,
                ind_bar_plot = norm_pca_res$pca_contrib_ind,
                pca_res = norm_pca_res$pca_res))
    
  })
  
  ## output the table
  observeEvent(input$act_step12_draw_ribo_pca, {
    
    # browser()
    output$out_step12_ribo_eigen <- DT::renderDataTable({
      DT::datatable(step12_ribo_plot()$pca_res$var.coord, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "Ribo_seq_filter"), 
                                     list(extend = 'excel', filename = "Ribo_seq_filter"))
                    ))
    })
    
    output$out_step12_ribo_var_contrib <- DT::renderDataTable({
      DT::datatable(step12_ribo_plot()$pca_res$var.contrib, rownames = TRUE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "Ribo_seq_filter"), 
                                     list(extend = 'excel', filename = "Ribo_seq_filter"))
                    ))
    })
    
    output$out_step12_ribo_ind_contrib <- DT::renderDataTable({
      DT::datatable(step12_ribo_plot()$pca_res$ind.contrib, rownames = TRUE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "Ribo_seq_filter"), 
                                     list(extend = 'excel', filename = "Ribo_seq_filter"))
                    ))
    })
    
    output$save_step12_ribo_eigen <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_ribo_name, "-PCA-results-", Sys.Date(), '.xlsx')
      },
      content = function(file) {
        write.xlsx(x = step12_ribo_plot()$pca_res, file = file, rowNames = TRUE)
      }
    )
    
    ## show the PCA plot
    output$out_step12_ribo_eigen_plot <- renderPlot(
      width = input$out_step12_ribo_width * 80,
      height = input$out_step12_ribo_height * 60,
      {step12_ribo_plot()$eig_plot})
    
    output$out_step12_ribo_scatter <- renderPlot(
      width = input$out_step12_ribo_width * 100,
      height = input$out_step12_ribo_height * 80,
      {step12_ribo_plot()$pca_plot})
    
    output$out_step12_ribo_variances <- renderPlot(
      width = input$out_step12_ribo_width * 80,
      height = input$out_step12_ribo_height * 100,
      {step12_ribo_plot()$var_bar_plot})
    
    output$out_step12_ribo_individual <- renderPlot(
      width = input$out_step12_ribo_width * 80,
      height = input$out_step12_ribo_height * 100,
      {step12_ribo_plot()$ind_bar_plot})
    
    # save the expression level plot
    output$save_step12_ribo_eigen_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_ribo_name, "-eigen-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step12_ribo_plot()$eig_var_plot, filename = file,
               width = input$out_step12_ribo_width, height = input$out_step12_ribo_height)
      }
    )
    
    output$save_step12_ribo_scatter <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_ribo_name, "-PCA-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step12_ribo_plot()$pca_plot, filename = file,
               width = input$out_step12_ribo_width, height = input$out_step12_ribo_height)
      }
    )
    
    output$save_step12_ribo_variances <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_ribo_name, "-variances-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step12_ribo_plot()$var_bar_plot, filename = file,
               width = input$out_step12_ribo_width, height = input$out_step12_ribo_height)
      }
    )
    
    output$save_step12_ribo_individual <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_ribo_name, "-individual-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step12_ribo_plot()$ind_bar_plot, filename = file,
               width = input$out_step12_ribo_width, height = input$out_step12_ribo_height)
      }
    )
    
  })
  
  ### 12.5 import the TE exprs #################
  step12_te_norm <- reactive({
    
    req(input$act_step12_import_te_norm)
    
    if (is.null(input$in_step12_te)) {return(NULL)}
    
    # import the raw data table
    suffix <- c(".genes.results", ".isoforms.results", "_cds_rpf", "_cds_rpm", "_cds_rpkm", "_cds_tpm", "Aligned.sortedByCoord.out.bam")
    ftcount_column <- c("Chr", "Start", "End", "Strand", "Length")
    
    te_norm <- fread(input$in_step12_te$datapath, sep = '\t', check.names = FALSE,
                     header="auto", stringsAsFactors = FALSE) %>% 
      dplyr::rename_with(~ gsub(paste(suffix, collapse = '|'), '', .), everything()) %>%
      dplyr::rename(Gene = colnames(.)[1]) %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::mutate(across(where(is.numeric), as.double))
    
    te_norm <- te_norm %>%
      dplyr::rename_all(~basename(names(te_norm)))
    
    # browser()
    
    # filter the gene average
    te_norm <- te_norm %>% 
      dplyr::filter(rowSums(.) > input$in_step12_te_rowsum)
    
    # filter the samples
    te_norm <- te_norm %>% 
      dplyr::mutate(Effective_Sample = rowSums(. > 0)) %>% 
      dplyr::mutate(Effective_Ratio = Effective_Sample / (ncol(.) - 1) * 100) %>% 
      dplyr::filter(Effective_Ratio >= input$in_step12_te_colsum) %>% 
      dplyr::select(-c(Effective_Sample, Effective_Ratio))
    
    # filter the specified genes
    if (input$in_step12_te_gene_name != "") {
      heat_gene_name <- gsub(' ', ',', input$in_step12_te_gene_name)
      heat_gene_name <- unlist(strsplit(heat_gene_name, ","))
      
      te_norm <- te_norm %>% 
        dplyr::filter(rownames(.) %in% heat_gene_name)
    }
    
    # filter the top expression genes
    if (input$in_step12_te_top > 0) {
      te_norm <- te_norm %>% 
        dplyr::mutate(rowsum = rowSums(.)) %>% 
        dplyr::arrange(desc(rowsum)) %>% 
        head(input$in_step12_te_top) %>%
        dplyr::select(-rowsum)
    }
    
    # browser()
    # convert the norm to log2(norm  + 1)
    if (input$in_step12_te_log2) {
      te_norm <- log2(te_norm + 1)
    }
    
    # convert the rowname to column
    te_norm <- te_norm %>% 
      tibble::rownames_to_column("Gene") %>% 
      dplyr::arrange(Gene)
    
    
    return(te_norm)
  })
  
  ## output the table
  observeEvent(input$act_step12_import_te_norm, {
    # browser()
    
    ## output the table
    output$out_step12_te_norm <- DT::renderDataTable({
      DT::datatable(step12_te_norm(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "TE_filter"), 
                                     list(extend = 'excel', filename = "TE_filter"))
                    ))
    })
    
  })
  
  ### 12.6 draw the TE PCA plot #################
  step12_te_plot <- reactive({
    
    req(input$act_step12_draw_te_pca)
    
    if (is.null(step12_te_norm())) {return(NULL)}
    
    # browser()
    source("R/draw_PCA.R")
    
    if (import_design_clicked()) {
      # make the pheno table
      phenodata <- step1_design() %>% 
        dplyr::filter(Sample %in% colnames(step12_te_norm())) %>%
        dplyr::filter(SeqType == "TE") %>%
        dplyr::mutate(Group = factor(!!sym(input$in_step12_te_group), 
                                     levels = unique(!!sym(input$in_step12_te_group)))) %>% 
        dplyr::mutate(row_name = Sample) %>% 
        tibble::column_to_rownames(var = "row_name")
      
    } else {
      print("Step 1 Design has not been imported.")
      {return(NULL)}
    }
    
    group <- input$in_step12_te_group
    
    dot_color <- input$in_step12_te_dot_color
    dot_alpha <- input$in_step12_te_dot_alpha
    
    line_color <- input$in_step12_te_line_color
    line_width <- input$in_step12_te_line_width
    
    xmax <- input$in_step12_te_xmax
    xmin <- input$in_step12_te_xmin
    
    ymax <- input$in_step12_te_ymax
    ymin <- input$in_step12_te_ymin
    
    gene_num <- 50
    
    dot_size <- input$in_step12_te_dot_size
    font_size <- input$in_step12_te_font_size
    
    label_sample <- input$in_step12_te_label_sample
    label_size <- input$in_step12_te_label_size
    # browser()
    
    # draw the PCA plot
    norm_pca_res <- draw_pca(in_data = step12_te_norm() %>% tibble::column_to_rownames("Gene"), 
                             phenodata = phenodata, 
                             x_min = xmin, x_max = xmax, 
                             y_min = ymin, y_max = ymax, 
                             gene_num = gene_num,
                             label_sample = label_sample,
                             label_size = label_size,
                             font_size = font_size, 
                             dot_size = dot_size,
                             fill_color = dot_color, 
                             fill_alpha = dot_alpha, 
                             group = group, 
                             title = "TE")
    
    # browser()
    return(list(eig_plot = norm_pca_res$eig_plot,
                pca_plot = norm_pca_res$pca_scatter,
                var_bar_plot = norm_pca_res$pca_contrib_var,
                ind_bar_plot = norm_pca_res$pca_contrib_ind,
                pca_res = norm_pca_res$pca_res))
    
  })
  
  ## output the table
  observeEvent(input$act_step12_draw_te_pca, {
    
    # browser()
    output$out_step12_te_eigen <- DT::renderDataTable({
      DT::datatable(step12_te_plot()$pca_res$var.coord, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "TE_filter"), 
                                     list(extend = 'excel', filename = "TE_filter"))
                    ))
    })
    
    output$out_step12_te_var_contrib <- DT::renderDataTable({
      DT::datatable(step12_te_plot()$pca_res$var.contrib, rownames = TRUE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "TE_filter"), 
                                     list(extend = 'excel', filename = "TE_filter"))
                    ))
    })
    
    output$out_step12_te_ind_contrib <- DT::renderDataTable({
      DT::datatable(step12_te_plot()$pca_res$ind.contrib, rownames = TRUE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "TE_filter"), 
                                     list(extend = 'excel', filename = "TE_filter"))
                    ))
    })
    
    output$save_step12_te_eigen <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_te_name, "-PCA-results-", Sys.Date(), '.xlsx')
      },
      content = function(file) {
        write.xlsx(x = step12_te_plot()$pca_res, file = file, rowNames = TRUE)
      }
    )
    
    ## show the PCA plot
    output$out_step12_te_eigen_plot <- renderPlot(
      width = input$out_step12_te_width * 80,
      height = input$out_step12_te_height * 60,
      {step12_te_plot()$eig_plot})
    
    output$out_step12_te_scatter <- renderPlot(
      width = input$out_step12_te_width * 100,
      height = input$out_step12_te_height * 80,
      {step12_te_plot()$pca_plot})
    
    output$out_step12_te_variances <- renderPlot(
      width = input$out_step12_te_width * 80,
      height = input$out_step12_te_height * 100,
      {step12_te_plot()$var_bar_plot})
    
    output$out_step12_te_individual <- renderPlot(
      width = input$out_step12_te_width * 80,
      height = input$out_step12_te_height * 100,
      {step12_te_plot()$ind_bar_plot})
    
    # save the expression level plot
    output$save_step12_te_eigen_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_te_name, "-eigen-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step12_te_plot()$eig_var_plot, filename = file,
               width = input$out_step12_te_width, height = input$out_step12_te_height)
      }
    )
    
    output$save_step12_te_scatter <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_te_name, "-PCA-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step12_te_plot()$pca_plot, filename = file,
               width = input$out_step12_te_width, height = input$out_step12_te_height)
      }
    )
    
    output$save_step12_te_variances <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_te_name, "-variances-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step12_te_plot()$var_bar_plot, filename = file,
               width = input$out_step12_te_width, height = input$out_step12_te_height)
      }
    )
    
    output$save_step12_te_individual <- downloadHandler(
      filename = function() {
        paste0(input$out_step12_te_name, "-individual-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step12_te_plot()$ind_bar_plot, filename = file,
               width = input$out_step12_te_width, height = input$out_step12_te_height)
      }
    )
    
  })
  
  
  ############################################################
  ## 13 run the deltaTE #################
  
  ### 13.1 import the design for edgeR analysis ##############
  
  observeEvent(input$act_step1_import_design, {
    col_name <- names(step1_design())
    updateSelectizeInput(session, "in_step13_edger_group", choices = col_name, selected = col_name[3], server = TRUE)
    # browser()
    observe({
      req(input$in_step13_edger_group)
      select_choices <- step1_design() %>% dplyr::filter(SeqType == input$in_step13_edger_type)
      select_choices <- select_choices[[input$in_step13_edger_group]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step13_edger_group1", choices = select_choices, selected = select_choices[1], server = TRUE)
      updateSelectizeInput(session, "in_step13_edger_group2", choices = select_choices, selected = select_choices[2], server = TRUE)
    })
  })
  
  ### 13.2 import the annotation for edgeR analysis ##########
  observeEvent(input$act_step1_import_anno, {
    col_name <- names(step1_annotation())
    updateSelectizeInput(session, "in_step13_edger_anno_lab", choices = col_name, selected = col_name[2], server = TRUE)
    updateSelectizeInput(session, "in_step13_edger_anno_col", choices = col_name, selected = col_name[3], server = TRUE)
  })
  
  ### 13.3 import the exprs data for edgeR analysis ##########
  step13_edger_exprs <- reactive({
    
    req(input$act_step13_import_edger)
    
    source("R/retrieve_group.R")
    
    if (is.null(input$in_step13_edger)) {return(NULL)}
    
    # import the raw data table
    suffix <- c(".genes.results", ".isoforms.results", "_cds_rpf", "_cds_rpm", "_cds_rpkm", "_cds_tpm", "Aligned.sortedByCoord.out.bam")
    ftcount_column <- c("Chr", "Start", "End", "Strand", "Length")
    
    edger_exprs <- fread(input$in_step13_edger$datapath, sep = '\t', check.names = FALSE,
                         header = "auto", stringsAsFactors = FALSE) %>% 
      dplyr::rename_with(~ gsub(paste(suffix, collapse = '|'), '', .), everything()) %>%
      dplyr::rename(Gene = colnames(.)[1]) %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::select(-contains(ftcount_column)) %>% 
      dplyr::mutate(across(where(is.numeric), as.double))
    
    edger_exprs <- edger_exprs %>%
      dplyr::rename_all(~basename(names(edger_exprs)))
    
    # select the samples
    if (import_design_clicked()) {
      
      # browser()
      if (is.null(input$in_step13_edger_type)) {return(NULL)}
      if (nchar(input$in_step13_edger_group) == 0) {return(NULL)}
      if (nchar(input$in_step13_edger_group1) == 0) {return(NULL)}
      if (nchar(input$in_step13_edger_group2) == 0) {return(NULL)}
      
      edger_exprs <- retrieve_group(exprs = edger_exprs,
                                    design_table = step1_design() %>% 
                                      dplyr::filter(SeqType == input$in_step13_edger_type),
                                    sample_flt = "Sample",
                                    group_flt = input$in_step13_edger_group,
                                    group1 = input$in_step13_edger_group1,
                                    group2 = input$in_step13_edger_group2,
                                    rowsum = input$in_step13_edger_rowsum,
                                    stdev = input$in_step13_edger_stdev)
      outmess <- NULL
      return(list(edger_exprs = edger_exprs$gene_exprs,
                  flt_design = edger_exprs$flt_design,
                  outmess = outmess))
      
    } else {
      # convert the rowname to column
      edger_exprs <- edger_exprs %>% 
        tibble::rownames_to_column("Gene") %>% 
        dplyr::arrange(Gene)
      outmess <- c("Step 1 Design has not been imported!")
      
      return(list(edger_exprs = edger_exprs, 
                  flt_design = NULL,
                  outmess = outmess))
    }
    
  })
  
  ## output the table
  observeEvent(input$act_step13_import_edger, {
    
    ## output the table
    output$out_step13_edger_exprs <- DT::renderDataTable(server = F, {
      DT::datatable(step13_edger_exprs()$edger_exprs, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "edgeR_exprs"), 
                                     list(extend = 'excel', filename = "edgeR_exprs"))
                    ))
    })
    
    output$out_step13_edger_design <- DT::renderDataTable(server = T, {
      DT::datatable(step13_edger_exprs()$flt_design, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "edgeR_design"), 
                                     list(extend = 'excel', filename = "edgeR_design"))
                    ))
    })
    
    output$out_step13_edger_exprs_info <- renderText({
      step13_edger_exprs()$outmess
    })
    
  })
  
  ### 13.4 run edgeR analysis ################################
  step13_run_edger <- reactive({
    
    req(input$act_step13_run_edger)
    
    source("R/run_edgeR.R")
    
    if (is.null(step1_design())) {return(NULL)}
    if (is.null(step13_edger_exprs())) {return(NULL)}
    
    # browser()
    
    # filter the gene annotation
    if (!is.null(input$in_step13_edger_anno_lab) & !is.null(input$in_step13_edger_anno_col)) 
    {
      flt_columns <- unique(c(input$in_step13_edger_anno_lab,
                              input$in_step13_edger_anno_col))
      
      anno <- step1_annotation() %>% 
        dplyr::select(any_of(flt_columns))
    } else {
      anno <- NULL
    }
    
    # browser()
    
    # run the edger analysis
    edger_res <- run_edger(
      exprs = step13_edger_exprs()$edger_exprs,
      design = step13_edger_exprs()$flt_design, # filtered design table
      method = input$in_step13_edger_test,
      bcv = input$in_step13_edger_bcv, # bcv for filter the DEGs
      padj_cutoff = input$in_step13_edger_padj, # padj for filter the DEGs
      log2fc_cutoff = input$in_step13_edger_log2fc, # log2fc for filter the DEGs
      flt_sig = F, # filter the sig diff genes
      anno = anno, # annotation table
      join_flag = input$in_step13_edger_anno_lab # join flag of the annotation table
    )
    
    # summary_table <- summary(edger_res$DTEGs)
    
    return(list(
      DEGs = edger_res$DEGs,
      summary_table = edger_res$summary_table
    ))
  })
  
  
  ## output the table
  observeEvent(input$act_step13_run_edger, {
    
    ## output the table
    output$out_step13_edger_res <- DT::renderDataTable(server = T, {
      DT::datatable(step13_run_edger()$DEGs, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "DEGs"), 
                                     list(extend = 'excel', filename = "DEGs"))
                    ))
    })
    
    output$out_step13_edger_summary <- DT::renderDataTable(server = T, {
      DT::datatable(step13_run_edger()$summary_table, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "edger_summary"),
                                     list(extend = 'excel', filename = "edger_summary"))
                    ))
    })
    
    output$save_step13_edger_res <- downloadHandler(
      filename = function() {
        paste0(input$out_step13_edger_name, "-edger-results-", Sys.Date(), '.xlsx')
      },
      content = function(file) {
        write.xlsx(x = step13_run_edger()$DEGs, sheetName = input$out_step13_edger_sheet,
                   file = file, rowNames = FALSE)
      }
    )
    
  })
  
  
  ### 13.5 import the design for DESeq2 analysis ##############
  observeEvent(input$act_step1_import_design, {
    col_name <- names(step1_design())
    updateSelectizeInput(session, "in_step13_deseq2_group", choices = col_name, selected = col_name[3], server = TRUE)
    # browser()
    observe({
      req(input$in_step13_deseq2_group)
      select_choices <- step1_design() %>% dplyr::filter(SeqType == input$in_step13_deseq2_type)
      select_choices <- select_choices[[input$in_step13_deseq2_group]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step13_deseq2_group1", choices = select_choices, selected = select_choices[1], server = TRUE)
      updateSelectizeInput(session, "in_step13_deseq2_group2", choices = select_choices, selected = select_choices[2], server = TRUE)
    })
  })
  
  ### 13.6 import the annotation for DESeq2 analysis ##########
  observeEvent(input$act_step1_import_anno, {
    col_name <- names(step1_annotation())
    updateSelectizeInput(session, "in_step13_deseq2_anno_lab", choices = col_name, selected = col_name[2], server = TRUE)
    updateSelectizeInput(session, "in_step13_deseq2_anno_col", choices = col_name, selected = col_name[3], server = TRUE)
  })
  
  ### 13.7 import the exprs data for DESeq2 analysis #########
  step13_deseq2_exprs <- reactive({
    
    req(input$act_step13_import_deseq2)
    
    source("R/retrieve_group.R")
    
    if (is.null(input$in_step13_deseq2)) {return(NULL)}
    
    # import the raw data table
    # browser()
    
    suffix <- c(".genes.results", ".isoforms.results", "_cds_rpf", "_cds_rpm", "_cds_rpkm", "_cds_tpm", "Aligned.sortedByCoord.out.bam")
    ftcount_column <- c("Chr", "Start", "End", "Strand", "Length")
    
    gene_exprs <- fread(input$in_step13_deseq2$datapath, sep = '\t', check.names = FALSE,
                        header = "auto", stringsAsFactors = FALSE) %>% 
      dplyr::rename_with(~ gsub(paste(suffix, collapse = '|'), '', .), everything()) %>%
      dplyr::rename(Gene = colnames(.)[1]) %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::select(-contains(ftcount_column)) %>% 
      dplyr::mutate(across(where(is.numeric), as.double))
    
    gene_exprs <- gene_exprs %>%
      dplyr::rename_all(~basename(names(gene_exprs)))
    
    # select the samples
    if (import_design_clicked()) {
      
      # browser()
      if (is.null(input$in_step13_deseq2_type)) {return(NULL)}
      if (nchar(input$in_step13_deseq2_group) == 0) {return(NULL)}
      if (nchar(input$in_step13_deseq2_group1) == 0) {return(NULL)}
      if (nchar(input$in_step13_deseq2_group2) == 0) {return(NULL)}
      
      gene_exprs <- retrieve_group(exprs = gene_exprs,
                                   design_table = step1_design() %>% 
                                     dplyr::filter(SeqType == input$in_step13_deseq2_type),
                                   sample_flt = "Sample",
                                   group_flt = input$in_step13_deseq2_group,
                                   group1 = input$in_step13_deseq2_group1,
                                   group2 = input$in_step13_deseq2_group2,
                                   rowsum = input$in_step13_deseq2_rowsum,
                                   stdev = input$in_step13_deseq2_stdev)
      outmess <- NULL
      return(list(gene_exprs = gene_exprs$gene_exprs, 
                  flt_design = gene_exprs$flt_design, 
                  outmess = outmess))
    } else {
      # convert the rowname to column
      gene_exprs <- gene_exprs %>% 
        tibble::rownames_to_column("Gene") %>% 
        dplyr::arrange(Gene)
      outmess <- c("Step 1 Design has not been imported!")
      
      return(list(gene_exprs = gene_exprs, 
                  flt_design = NULL, 
                  outmess = outmess))
    }
    
  })
  
  ## output the table
  observeEvent(input$act_step13_import_deseq2, {
    
    ## output the table
    output$out_step13_deseq2_exprs <- DT::renderDataTable(server = T, {
      DT::datatable(step13_deseq2_exprs()$gene_exprs, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "deseq2_filter"), 
                                     list(extend = 'excel', filename = "deseq2_filter"))
                    ))
    })
    
    output$out_step13_deseq2_design <- DT::renderDataTable(server = T, {
      DT::datatable(step13_deseq2_exprs()$flt_design, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "deseq2_filter"), 
                                     list(extend = 'excel', filename = "deseq2_filter"))
                    ))
    })
    
    output$out_step13_deseq2_exprs_info <- renderText({
      step13_deseq2_exprs()$outmess
    })
    
  })
  
  ### 13.8 run DESeq2 analysis ###############################
  step13_run_deseq2 <- reactive({
    
    req(input$act_step13_run_deseq2)
    
    source("R/run_DESeq2.R")
    
    if (is.null(step1_design())) {return(NULL)}
    if (is.null(step13_deseq2_exprs())) {return(NULL)}
    
    # browser()
    
    # filter the gene annotation
    if (!is.null(input$in_step13_deseq2_anno_lab) & !is.null(input$in_step13_deseq2_anno_col)) 
    {
      flt_columns <- unique(c(input$in_step13_deseq2_anno_lab,
                              input$in_step13_deseq2_anno_col))
      
      anno <- step1_annotation() %>% 
        dplyr::select(any_of(flt_columns))
    } else {
      anno <- NULL
    }
    
    # browser()
    
    # run the deseq2 analysis
    deseq2_res <- run_deseq2(
      exprs = step13_deseq2_exprs()$gene_exprs,
      design = step13_deseq2_exprs()$flt_design, # filtered design table
      method = input$in_step13_deseq2_test,
      padj_cutoff = input$in_step13_deseq2_padj, # padj for filter the DEGs
      log2fc_cutoff = input$in_step13_deseq2_log2fc, # log2fc for filter the DEGs
      flt_sig = F, # filter the sig diff genes
      anno = anno, # annotation table
      join_flag = input$in_step13_deseq2_anno_lab # join flag of the annotation table
    )
    
    # summary_table <- summary(deseq2_res$DTEGs)
    
    return(list(
      DEGs = deseq2_res$DEGs,
      summary_table = deseq2_res$summary_table
    ))
  })
  
  
  ## output the table
  observeEvent(input$act_step13_run_deseq2, {
    
    ## output the table
    output$out_step13_deseq2_res <- DT::renderDataTable(server = T, {
      DT::datatable(step13_run_deseq2()$DEGs, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "DEGs"), 
                                     list(extend = 'excel', filename = "DEGs"))
                    ))
    })
    
    output$out_step13_deseq2_summary <- DT::renderDataTable(server = T, {
      DT::datatable(step13_run_deseq2()$summary_table, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "deseq2_summary"),
                                     list(extend = 'excel', filename = "deseq2_summary"))
                    ))
    })
    
    output$save_step13_deseq2_res <- downloadHandler(
      filename = function() {
        paste0(input$out_step13_deseq2_name, "-deseq2-results-", Sys.Date(), '.xlsx')
      },
      content = function(file) {
        write.xlsx(x = step13_run_deseq2()$DEGs, sheetName = input$out_step13_deseq2_sheet,
                   file = file, rowNames = FALSE)
      }
    )
    
  })
  
  
  ### 13.9 import the design for deltaTE analysis ############
  observeEvent(input$act_step1_import_design, {
    col_name <- names(step1_design())
    updateSelectizeInput(session, "in_step13_deltate_rna_group", choices = col_name, selected = col_name[3], server = TRUE)
    updateSelectizeInput(session, "in_step13_deltate_ribo_group", choices = col_name, selected = col_name[3], server = TRUE)
    
    # browser()
    observe({
      req(input$in_step13_deltate_rna_group)
      select_choices <- step1_design() %>% dplyr::filter(SeqType == "RNA")
      select_choices <- select_choices[[input$in_step13_deltate_rna_group]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step13_deltate_group1", choices = select_choices, selected = select_choices[1], server = TRUE)
      updateSelectizeInput(session, "in_step13_deltate_group2", choices = select_choices, selected = select_choices[2], server = TRUE)
    })
    
    observe({
      req(input$in_step13_deltate_ribo_group)
      select_choices <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      select_choices <- select_choices[[input$in_step13_deltate_ribo_group]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step13_deltate_group3", choices = select_choices, selected = select_choices[1], server = TRUE)
      updateSelectizeInput(session, "in_step13_deltate_group4", choices = select_choices, selected = select_choices[2], server = TRUE)
    })
  })
  
  
  ### 13.10 import the annotation for deltaTE analysis #######
  observeEvent(input$act_step1_import_anno, {
    observe({
      req(input$act_step1_import_anno)
      col_name <- names(step1_annotation())
      updateSelectizeInput(session, "in_step13_deltate_anno_lab", choices = col_name, selected = col_name[2], server = TRUE)
      updateSelectizeInput(session, "in_step13_deltate_anno_col", choices = col_name, selected = col_name[3], server = TRUE)
    })
  })
  
  ### 13.11 import the exprs data for deltaTE analysis ######## 
  step13_deltate_exprs <- reactive({
    
    req(input$act_step13_import_deltate)
    
    source("R/retrieve_group.R")
    
    if (is.null(input$in_step13_deltate_rna)) {return(NULL)}
    if (is.null(input$in_step13_deltate_ribo)) {return(NULL)}
    
    # import the raw data table
    suffix <- c(".genes.results", ".isoforms.results", "_cds_rpf", "_cds_rpm", "_cds_rpkm", "_cds_tpm", "Aligned.sortedByCoord.out.bam")
    ftcount_column <- c("Chr", "Start", "End", "Strand", "Length")
    
    rna_exprs <- fread(input$in_step13_deltate_rna$datapath, sep = '\t', check.names = FALSE,
                       header = "auto", stringsAsFactors = FALSE) %>% 
      dplyr::rename_with(~ gsub(paste(suffix, collapse = '|'), '', .), everything()) %>%
      dplyr::rename(Gene = colnames(.)[1]) %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::select(-contains(ftcount_column)) %>% 
      dplyr::mutate(across(where(is.numeric), as.double)) %>% 
      round()
    
    ribo_exprs <- fread(input$in_step13_deltate_ribo$datapath, sep = '\t', check.names = FALSE,
                        header = "auto", stringsAsFactors = FALSE) %>%  
      dplyr::rename_with(~ gsub(paste(suffix, collapse = '|'), '', .), everything()) %>%
      dplyr::rename(Gene = colnames(.)[1]) %>% 
      tibble::column_to_rownames("Gene") %>%
      dplyr::mutate(across(where(is.numeric), as.double)) %>% 
      round()
    
    rna_exprs <- rna_exprs %>%
      dplyr::rename_all(~basename(names(rna_exprs)))
    
    ribo_exprs <- ribo_exprs %>%
      dplyr::rename_all(~basename(names(ribo_exprs)))
    
    # browser()
    
    # select the samples
    if (import_design_clicked()) {
      
      # browser()
      if (is.null(rna_exprs)) {return(NULL)}
      if (is.null(ribo_exprs)) {return(NULL)}
      
      if (nchar(input$in_step13_deltate_rna_group) == 0) {return(NULL)}
      if (nchar(input$in_step13_deltate_ribo_group) == 0) {return(NULL)}
      
      if (nchar(input$in_step13_deltate_group1) == 0) {return(NULL)}
      if (nchar(input$in_step13_deltate_group2) == 0) {return(NULL)}
      if (nchar(input$in_step13_deltate_group3) == 0) {return(NULL)}
      if (nchar(input$in_step13_deltate_group4) == 0) {return(NULL)}
      
      # select the design
      rna_exprs <- retrieve_group(exprs = rna_exprs,
                                  design_table = step1_design() %>% dplyr::filter(SeqType == "RNA"),
                                  sample_flt = "Sample",
                                  group_flt = input$in_step13_deltate_rna_group,
                                  group1 = input$in_step13_deltate_group1,
                                  group2 = input$in_step13_deltate_group2,
                                  rowsum = input$in_step13_deltate_rowsum,
                                  stdev = input$in_step13_deltate_stdev)
      
      ribo_exprs <- retrieve_group(exprs = ribo_exprs,
                                   design_table = step1_design() %>% dplyr::filter(SeqType == "RIBO"),
                                   sample_flt = "Sample",
                                   group_flt = input$in_step13_deltate_ribo_group,
                                   group1 = input$in_step13_deltate_group3,
                                   group2 = input$in_step13_deltate_group4,
                                   rowsum = input$in_step13_deltate_rowsum,
                                   stdev = input$in_step13_deltate_stdev)
      
      flt_design <- rbind(rna_exprs$flt_design, ribo_exprs$flt_design) %>% 
        dplyr::distinct_all()
      
      outmess <- NULL
      
      return(list(rna_exprs = rna_exprs$gene_exprs,
                  ribo_exprs = ribo_exprs$gene_exprs,
                  flt_design = flt_design,
                  outmess = outmess))
      
    } else {
      
      # convert the rowname to column
      rna_exprs <- rna_exprs %>%
        tibble::rownames_to_column("Gene") %>%
        dplyr::arrange(Gene)
      
      ribo_exprs <- ribo_exprs %>%
        tibble::rownames_to_column("Gene") %>%
        dplyr::arrange(Gene)
      
      outmess <- c("Step 1 Design has not been imported!")
    }
    
    return(list(rna_exprs = rna_exprs,
                ribo_exprs = ribo_exprs,
                flt_design = NULL,
                outmess = outmess))
  })
  
  ## output the table
  observeEvent(input$act_step13_import_deltate, {
    
    ## output the table
    output$out_step13_deltate_rna_exprs <- DT::renderDataTable(server = T, {
      DT::datatable(step13_deltate_exprs()$rna_exprs, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "deltate_rna_filter"), 
                                     list(extend = 'excel', filename = "deltate_rna_filter"))
                    ))
    })
    
    output$out_step13_deltate_ribo_exprs <- DT::renderDataTable(server = T, {
      DT::datatable(step13_deltate_exprs()$ribo_exprs, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "deltate_ribo_filter"), 
                                     list(extend = 'excel', filename = "deltate_ribo_filter"))
                    ))
    })
    
    output$out_step13_deltate_design <- DT::renderDataTable(server = T, {
      DT::datatable(step13_deltate_exprs()$flt_design, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "deltate_design"), 
                                     list(extend = 'excel', filename = "deltate_design"))
                    ))
    })
    
    output$out_step13_deltate_rna_info <- renderText({
      step13_deltate_exprs()$outmess
    })
    
    output$out_step13_deltate_ribo_info <- renderText({
      step13_deltate_exprs()$outmess
    })
    
  })
  
  ### 13.12 run deltaTE analysis ##############################
  step13_run_deltate <- reactive({
    
    req(input$act_step13_run_deltate)
    
    source("R/run_deltaTE.R")
    
    if (is.null(step1_design())) {return(NULL)}
    if (is.null(step13_deltate_exprs()$rna_exprs)) {return(NULL)}
    if (is.null(step13_deltate_exprs()$ribo_exprs)) {return(NULL)}
    
    # browser()
    
    # filter the gene annotation
    if (!is.null(input$in_step13_deltate_anno_lab) & !is.null(input$in_step13_deltate_anno_col)) 
    {
      flt_columns <- unique(c(input$in_step13_deltate_anno_lab,
                              input$in_step13_deltate_anno_col))
      
      anno <- step1_annotation() %>% 
        dplyr::select(any_of(flt_columns))
    } else {
      anno <- NULL
    }
    
    # browser()
    
    # run the deltaTE analysis
    deltate_res <- run_delta_te(
      dt_ribo = step13_deltate_exprs()$ribo_exprs,
      dt_mrna = step13_deltate_exprs()$rna_exprs,
      dt_design = step13_deltate_exprs()$flt_design, # filtered design table
      padj_cutoff = input$in_step13_deltate_padj, # padj for filter the DEGs
      log2fc_cutoff = input$in_step13_deltate_log2fc, # log2fc for filter the DEGs
      flt_sig = F, # filter the sig diff genes
      anno = anno, # annotation table
      join_flag = input$in_step13_deltate_anno_lab # join flag of the annotation table
    )
    
    # summary_table <- summary(deltate_res$DTEGs)
    
    return(list(
      DTEGs = deltate_res$DTEGs,
      rna_DEGs = deltate_res$rna_DEGs,
      ribo_DEGs = deltate_res$ribo_DEGs,
      forwarded = deltate_res$forwarded,
      exclusive = deltate_res$exclusive,
      intensified = deltate_res$intensified,
      buffered = deltate_res$buffered,
      summary_table = deltate_res$summary_table
    ))
  })
  
  
  ## output the table
  observeEvent(input$act_step13_run_deltate, {
    
    ## output the table
    output$out_step13_deltate_res <- DT::renderDataTable(server = T, {
      DT::datatable(step13_run_deltate()$DTEGs, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "DTEGs"), 
                                     list(extend = 'excel', filename = "DTEGs"))
                    ))
    })
    
    output$out_step13_deltate_rna_degs <- DT::renderDataTable(server = T, {
      DT::datatable(step13_run_deltate()$rna_DEGs, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "rna_DEGs"), 
                                     list(extend = 'excel', filename = "rna_DEGs"))
                    ))
    })
    
    output$out_step13_deltate_ribo_degs <- DT::renderDataTable(server = T, {
      DT::datatable(step13_run_deltate()$ribo_DEGs, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "ribo_DEGs"), 
                                     list(extend = 'excel', filename = "ribo_DEGs"))
                    ))
    })
    
    output$out_step13_deltate_summary <- DT::renderDataTable({
      DT::datatable(step13_run_deltate()$summary_table, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "deltate_summary"),
                                     list(extend = 'excel', filename = "deltate_summary"))
                    ))
    })
    
    output$save_step13_deltate_res <- downloadHandler(
      filename = function() {
        paste0(input$out_step13_deltate_name, "-deltaTE-results-", Sys.Date(), '.xlsx')
      },
      content = function(file) {
        write.xlsx(x = step13_run_deltate(), file = file, rowNames = FALSE)
      }
    )
    
  })
  
  ############################################################
  ## 14 draw the volcano plot #################
  
  ### 14.1 import the degs data ##############################
  step14_vol_degs <- reactive({
    
    req(input$act_step14_import_vol_degs)
    
    if (is.null(input$in_step14_vol)) {return(NULL)}
    
    # import the raw data table
    vol_degs <- read.xlsx(input$in_step14_vol$datapath,
                          sheet = input$in_step14_vol_sheet,
                          colNames = T, rowNames = F)
    
    # label the down/up/ns genes again with specified padj and log2fc
    flt_padj <- input$in_step14_vol_padj
    flt_logfc <- input$in_step14_vol_logfc
    
    logfc_name <- input$in_step14_vol_x
    padj_name <- input$in_step14_vol_y
    
    # browser()
    
    if (!is.null(flt_padj) & !is.null(flt_logfc)) {
      flt_logfc <- abs(flt_logfc)
      vol_degs <- vol_degs %>%
        dplyr::mutate(DEGs = if_else(!!sym(padj_name) < flt_padj & !!sym(logfc_name) >= flt_logfc, "UP",
                                     if_else(!!sym(padj_name) < flt_padj & !!sym(logfc_name) <= -flt_logfc, "DOWN",
                                             "NS")))
    }
    
    
    return(vol_degs)
  })
  
  ## output the table
  observeEvent(input$act_step14_import_vol_degs, {
    
    ## output the table
    output$out_step14_vol_degs <- DT::renderDataTable(server = T, {
      DT::datatable(step14_vol_degs(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "vol_degs_filter"), 
                                     list(extend = 'excel', filename = "vol_degs_filter"))
                    ))
    })
    
  })
  
  ### 14.2 draw the volcano plot ############################
  step14_vol_plot <- reactive({
    
    req(input$act_step14_draw_volcano)
    
    # browser()
    source("R/draw_volcano.R")
    
    if (is.null(step14_vol_degs())) {return(NULL)}
    
    # filter the gene annotation
    vol_label_list <- unlist(strsplit(gsub(' ', '', input$in_step14_vol_label_degs), ","))
    
    # browser()
    
    # draw the volcano plot
    volcano_plot <- draw_volcano(degs = step14_vol_degs(),
                                 log2fc_column = input$in_step14_vol_x,
                                 pvalue_column = input$in_step14_vol_y,
                                 log2fc = input$in_step14_vol_logfc, 
                                 pvalue = input$in_step14_vol_padj,
                                 class = "DEGs",
                                 up_color = input$in_step14_vol_up_color,
                                 down_color = input$in_step14_vol_down_color,
                                 ns_color = input$in_step14_vol_ns_color,
                                 dot_size = input$in_step14_vol_dot_size,
                                 font_size = input$in_step14_vol_font_size,
                                 x_max = input$in_step14_vol_xlim_max,
                                 x_min = input$in_step14_vol_xlim_min,
                                 
                                 fill_0 = input$in_step14_vol_pvalue0,
                                 
                                 vol_sqrt = input$in_step14_vol_sqrt,
                                 vol_log2 = input$in_step14_vol_log2,
                                 
                                 vol_title = input$in_step14_vol_title,
                                 vol_xlab = input$in_step14_vol_xlab,
                                 vol_ylab = input$in_step14_vol_ylab,
                                 
                                 label_type = input$in_step14_vol_label_type,
                                 label_col = input$in_step14_vol_label_col,
                                 label_color = input$in_step14_vol_label_color,
                                 label_size = input$in_step14_vol_label_size,
                                 label_list = vol_label_list,
                                 label_num = input$in_step14_vol_label_num,
                                 overlaps_num = input$in_step14_vol_overlaps_num
    )
    # draw the MA plot
    source("R/draw_MAplot.R")
    ma_plot <- draw_MAplot(degs = step14_vol_degs(),
                           log2fc_column = input$in_step14_vol_x,
                           basemean_column = input$in_step14_vol_basemean,
                           pvalue_column = input$in_step14_vol_y,
                           log2fc = input$in_step14_vol_logfc, 
                           pvalue = input$in_step14_vol_padj,
                           class = "DEGs",
                           
                           up_color = input$in_step14_vol_up_color,
                           down_color = input$in_step14_vol_down_color,
                           ns_color = input$in_step14_vol_ns_color,
                           dot_size = input$in_step14_vol_dot_size,
                           font_size = input$in_step14_vol_font_size,
                           
                           x_max = input$in_step14_vol_xlim_max,
                           x_min = input$in_step14_vol_xlim_min,
                           
                           fill_0 = input$in_step14_vol_pvalue0,
                           
                           vol_sqrt = input$in_step14_vol_sqrt,
                           vol_log2 = input$in_step14_vol_log2,
                           
                           vol_title = input$in_step14_vol_title,
                           vol_xlab = input$in_step14_vol_ma_xlab,
                           vol_ylab = input$in_step14_vol_ma_ylab,
                           
                           label_type = input$in_step14_vol_label_type,
                           
                           label_col = input$in_step14_vol_label_col,
                           label_color = input$in_step14_vol_label_color,
                           label_size = input$in_step14_vol_label_size,
                           label_list = vol_label_list,
                           label_num = input$in_step14_vol_label_num,
                           overlaps_num = input$in_step14_vol_overlaps_num
    )
    
    return(list(volcano_plot = volcano_plot, ma_plot = ma_plot))
    
  })
  
  
  ## output the figure
  observeEvent(input$act_step14_draw_volcano, {
    
    ## output the figure
    output$out_step14_vol_plot <- renderPlot(
      width = input$out_step14_vol_width * 100,
      height = input$out_step14_vol_height * 100,
      {step14_vol_plot()$volcano_plot})
    
    output$out_step14_vol_maplot <- renderPlot(
      width = input$out_step14_vol_width * 100,
      height = input$out_step14_vol_height * 80,
      {step14_vol_plot()$ma_plot})
    
    # save the volcano plot
    output$save_step14_vol_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step14_vol_out_name, "-volcano-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step14_vol_plot()$volcano_plot, filename = file,
               width = input$out_step14_vol_width, height = input$out_step14_vol_height)
      }
    )
    
    
    output$save_step14_vol_maplot <- downloadHandler(
      filename = function() {
        paste0(input$out_step14_vol_out_name, "-MA-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step14_vol_plot()$ma_plot, filename = file,
               width = input$out_step14_vol_width, height = input$out_step14_vol_height)
      }
    )
    
  })
  
  
  ### 14.3 import the degs data ##############################
  step14_quadra_degs <- reactive({
    
    req(input$act_step14_import_quadra_degs)
    
    if (is.null(input$in_step14_quadra_degs_x)) {return(NULL)}
    if (is.null(input$in_step14_quadra_degs_y)) {return(NULL)}
    
    source("R/read_quadra.R")
    # browser()
    
    degs_x <- input$in_step14_quadra_degs_x_sheet
    degs_y <- input$in_step14_quadra_degs_y_sheet
    
    if (degs_x == "rna_DEGs" & degs_y == "ribo_DEGs") {
      category <- "RNA-Ribo"
    } else if (degs_x == "rna_DEGs" & degs_y == "DTEGs") {
      category <- "RNA-TE"
    } else if (degs_x == "ribo_DEGs" & degs_y == "DTEGs") {
      category <- "Ribo-TE"
    }
    
    # import the raw data table
    quadra_degs <- read_quadra(
      data1_file = input$in_step14_quadra_degs_x$datapath,
      data1_sheet = input$in_step14_quadra_degs_x_sheet,
      data2_file = input$in_step14_quadra_degs_y$datapath,
      data2_sheet = input$in_step14_quadra_degs_y_sheet,
      category = category,
      gene_column = "Gene",
      label_column = input$in_step14_quadra_label_col,
      logfc_column = input$in_step14_quadra_x,
      padj_column = input$in_step14_quadra_y,
      logfc_cutoff = input$in_step14_quadra_logfc,
      pvalue_cutoff = input$in_step14_quadra_padj
    )
    
    return(quadra_degs)
  })
  
  ## output the table
  observeEvent(input$act_step14_import_quadra_degs, {
    
    ## output the table
    output$out_step14_quadra_degs <- DT::renderDataTable(server = T, {
      DT::datatable(step14_quadra_degs(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "quadra_degs_filter"), 
                                     list(extend = 'excel', filename = "quadra_degs_filter"))
                    ))
    })
    
    ## save the table
    output$save_step14_quadra_degs <- downloadHandler(
      filename = function() {
        paste0(input$out_step14_quadra_out_name, "-class-DEGs-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(step14_quadra_degs(), file = file, rowNames = FALSE, sheetName = "quadra_DEGs")
      }
    )
    
  })
  
  ### 14.4 draw the quadra plot ##############################
  step14_quadra_plot <- reactive({
    
    req(input$act_step14_draw_quadra_plot)
    
    # browser()
    source("R/draw_quadra.R")
    
    if (is.null(step14_quadra_degs())) {return(NULL)}
    
    # filter the gene annotation
    quadra_label_list <- unlist(strsplit(gsub(' ', '', input$in_step14_quadra_label_degs), ","))
    
    # browser()
    degs_x <- input$in_step14_quadra_degs_x_sheet
    degs_y <- input$in_step14_quadra_degs_y_sheet
    
    if (degs_x == "rna_DEGs" & degs_y == "ribo_DEGs") {
      xais <- "RNA_log2FC"
      yais <- "Ribo_log2FC"
    } else if (degs_x == "rna_DEGs" & degs_y == "DTEGs") {
      xais <- "RNA_log2FC"
      yais <- "TE_log2FC"
    } else if (degs_x == "ribo_DEGs" & degs_y == "DTEGs") {
      xais <- "Ribo_log2FC"
      yais <- "TE_log2FC"
    }
    
    # annotate the gene class
    if (input$in_step14_quadra_ns){
      quadra_degs <- step14_quadra_degs() %>% 
        dplyr::filter(DEGs != "NS")
      
      sections_count <- quadra_degs %>% 
        dplyr::group_by(Sections) %>% 
        dplyr::summarise(Count = n())
      
      quadra_degs <- quadra_degs %>% 
        dplyr::left_join(sections_count, by = "Sections") %>% 
        dplyr::mutate(Sections = paste0(Sections, " (", Count, ")")) %>% 
        dplyr::select(-Count)
      
    } else {

      sections_count <- step14_quadra_degs() %>% 
        dplyr::group_by(Sections) %>% 
        dplyr::summarise(Count = n())
      
      quadra_degs <- step14_quadra_degs() %>% 
        dplyr::left_join(sections_count, by = "Sections") %>% 
        dplyr::mutate(Sections = paste0(Sections, " (", Count, ")")) %>% 
        dplyr::select(-Count)
    }
    
    # draw the volcano plot
    quadra_plot <- draw_quadra(degs = quadra_degs,
                               x = xais,
                               y = yais,
                               x_logfc = input$in_step14_quadra_logfc,
                               y_logfc = input$in_step14_quadra_logfc,
                               
                               title = input$in_step14_quadra_title,
                               xlabel = input$in_step14_quadra_xlabel,
                               ylabel = input$in_step14_quadra_ylabel,
                               
                               remove_ns = input$in_step14_quadra_ns, 
                               
                               dot_size = input$in_step14_quadra_dot_size,
                               gene_class = "Sections",
                               color = input$in_step14_quadra_color,
                               font_size = input$in_step14_quadra_font_size,
                               
                               x_max = input$in_step14_quadra_xlim_max,
                               x_min = input$in_step14_quadra_xlim_min,
                               y_max = input$in_step14_quadra_ylim_max,
                               y_min = input$in_step14_quadra_ylim_min,
                               
                               label_column = input$in_step14_quadra_label_col,
                               label_color = input$in_step14_quadra_label_color,
                               label_size = input$in_step14_quadra_label_size,
                               label_list = quadra_label_list,
                               overlaps_num = input$in_step14_quadra_overlaps_num
    )
    
    return(quadra_plot)
    
  })
  
  
  ## output the figure
  observeEvent(input$act_step14_draw_quadra_plot, {
    
    ## output the figure
    output$out_step14_quadra_plot <- renderPlot(
      width = input$out_step14_quadra_width * 100,
      height = input$out_step14_quadra_height * 100,
      {step14_quadra_plot()})
    
    # save the quadrant plot
    output$save_step14_quadra_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step14_quadra_out_name, "-quadrant-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step14_quadra_plot(), filename = file,
               width = input$out_step14_quadra_width, height = input$out_step14_quadra_height)
      }
    )
    
  })
  
  ### 14.5 import the degs data ##############################
  step14_delta_degs <- reactive({
    
    req(input$act_step14_import_delta_degs)
    
    if (is.null(input$in_step14_delta_rna)) {return(NULL)}
    if (is.null(input$in_step14_delta_ribo)) {return(NULL)}
    
    source("R/read_quadra.R")
    # browser()
    
    # import the raw data table
    delta_degs <- read_quadra(
      data1_file = input$in_step14_delta_rna$datapath,
      data1_sheet = input$in_step14_delta_rna_sheet,
      data2_file = input$in_step14_delta_ribo$datapath,
      data2_sheet = input$in_step14_delta_ribo_sheet,
      category = "RNA-Ribo",
      delta_column = input$in_step14_delta_column,
      gene_column = input$ in_step14_delta_gene_column,
      label_column = input$in_step14_delta_label_col,
      logfc_column = 'log2FoldChange',
      padj_column = 'padj',
      logfc_cutoff = input$in_step14_delta_logfc,
      pvalue_cutoff = input$in_step14_delta_padj
    )
    
    return(delta_degs)
  })
  
  ## output the table
  observeEvent(input$act_step14_import_delta_degs, {
    
    ## output the table
    output$out_step14_delta_degs <- DT::renderDataTable(server = T, {
      DT::datatable(step14_delta_degs(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "delta_degs_filter"), 
                                     list(extend = 'excel', filename = "delta_degs_filter"))
                    ))
    })
    
    ## save the table
    output$save_step14_delta_degs <- downloadHandler(
      filename = function() {
        paste0(input$out_step14_delta_out_name, "-class-DEGs-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(step14_delta_degs(), file = file, rowNames = FALSE, sheetName = "delta_DEGs")
      }
    )
    
  })
  
  ### 14.6 draw the deltaTE plot ##############################
  step14_delta_plot <- reactive({
    
    req(input$act_step14_draw_delta_plot)
    
    source("R/draw_delta.R")
    
    if (is.null(step14_delta_degs())) {return(NULL)}
    
    # filter the gene annotation
    # browser()
    delta_label_list <- unlist(strsplit(gsub(' ', '', input$in_step14_delta_label_degs), ","))
    
    # browser()

    # filter out the 'others' and 'NS' genes
    delta_degs <- step14_delta_degs()
    
    if (input$in_step14_delta_others) {
      delta_degs <- delta_degs %>% 
        dplyr::filter(Delta != "others")
    }
    
    if (input$in_step14_delta_ns) {
      delta_degs <- delta_degs %>% 
        dplyr::filter(DEGs != "NS")
    }
    
    # annotate the gene class
    delta_count <- delta_degs %>% 
      dplyr::group_by(Delta) %>% 
      dplyr::summarise(Count = n())
    
    delta_degs <- delta_degs %>% 
      dplyr::left_join(delta_count, by = "Delta") %>% 
      dplyr::mutate(Delta = paste0(Delta, " (", Count, ")")) %>% 
      dplyr::select(-Count)

    # draw the volcano plot
    delta_plot <- draw_delta(degs = delta_degs,
                             x = "RNA_log2FC",
                             y = "Ribo_log2FC",
                             x_logfc = input$in_step14_delta_logfc,
                             y_logfc = input$in_step14_delta_logfc,
                             Delta = "Delta", # input$in_step14_delta_column,
                             
                             title = input$in_step14_delta_title,
                             xlabel = input$in_step14_delta_xlabel,
                             ylabel = input$in_step14_delta_ylabel,
                             
                             remove_others = input$in_step14_delta_others, 
                             remove_ns = input$in_step14_delta_ns, 
                             
                             dot_size = input$in_step14_delta_dot_size,
                             color = input$in_step14_delta_color,
                             font_size = input$in_step14_delta_font_size,
                             
                             x_max = input$in_step14_delta_xlim_max,
                             x_min = input$in_step14_delta_xlim_min,
                             y_max = input$in_step14_delta_ylim_max,
                             y_min = input$in_step14_delta_ylim_min,

                             label_column = input$in_step14_delta_label_col,
                             label_color = input$in_step14_delta_label_color,
                             label_size = input$in_step14_delta_label_size,
                             label_list = delta_label_list,
                             overlaps_num = input$in_step14_delta_overlaps_num
    )
    
    return(delta_plot)
    
  })
  
  
  ## output the figure
  observeEvent(input$act_step14_draw_delta_plot, {
    
    ## output the figure
    output$out_step14_delta_plot <- renderPlot(
      width = input$out_step14_delta_width * 100,
      height = input$out_step14_delta_height * 100,
      {step14_delta_plot()})
    
    # save the delta plot
    output$save_step14_delta_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step14_delta_out_name, "-delta-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step14_delta_plot(), filename = file,
               width = input$out_step14_delta_width, height = input$out_step14_delta_height)
      }
    )
    
  })
  
  
  ############################################################
  ## 15 draw the DEGs plot #################
  
  ### 15.1 import the datasets ##############################
  step15_venn_data_sets <- reactive({
    
    req(input$act_step15_venn_calc)
    
    # browser()
    
    # check the set 1 input data
    if(nchar(input$in_step15_venn_set_1) == 0){
      set1 <- NULL
    } else {
      set1 <- scan(text = input$in_step15_venn_set_1, what = character())
      set1 <- unique(set1)
    }
    
    # check the set 2 input data
    if(nchar(input$in_step15_venn_set_2) == 0){
      set2 <- NULL
    } else {
      set2 <- scan(text = input$in_step15_venn_set_2, what = character())
      set2 <- unique(set2)
    }
    
    # check the set 3 input data
    if(nchar(input$in_step15_venn_set_3) == 0){
      set3 <- NULL
    } else {
      set3 <- scan(text = input$in_step15_venn_set_3, what = character())
      set3 <- unique(set3)
    }
    
    # check the set 4 input data
    if(nchar(input$in_step15_venn_set_4) == 0){
      set4 <- NULL
    } else {
      set4 <- scan(text = input$in_step15_venn_set_4, what = character())
      set4 <- unique(set4)
    }
 
    # browser()
    
    # remove the duplicated data
    merge_data_set <- unique(c(set1, set2, set3, set4))
    venn_data_set <- UpSetR::fromList(list(set1, set2, set3, set4)) %>% as.data.frame()
    
    # rename the columns
    colnames(venn_data_set) <- c(input$in_step15_venn_set_1_title,
                                 input$in_step15_venn_set_2_title,
                                 input$in_step15_venn_set_3_title,
                                 input$in_step15_venn_set_4_title)

    # remove the columns with all zeros
    col_sums <- colSums(venn_data_set)
    venn_data_set <- venn_data_set[, col_sums != 0]
    # rownames(venn_data_set)[venn_data_set$Set3 == 1]
    
    # browser()
    
    # label the venn data set
    venn_data_set <- venn_data_set %>% 
      dplyr::mutate(across(everything(), ~replace(., . == 1, which(colnames(venn_data_set) == cur_column())))) %>% 
      dplyr::mutate(Flag = apply(., 1, paste0, collapse = "")) %>% 
      dplyr::mutate(Flag = str_replace_all(Flag, "0", "")) %>% 
      dplyr::mutate(IDs = merge_data_set)

    
    venn_set_list <- venn_data_set %>% 
      dplyr::select(-Flag) %>% 
      pivot_longer(cols = -IDs, names_to = "Set", values_to = "Count") %>% 
      filter(Count > 0) %>% 
      group_by(Set) %>%
      summarize(Genes = list(IDs)) %>% 
      deframe()
    
    
    return(list(venn_data_set = venn_data_set, set_list = venn_set_list))
  })
  
  ## output the table
  observeEvent(input$act_step15_venn_calc, {
    
    ## output the tables
    output$out_step15_venn_data_sets <- DT::renderDataTable(server = T, {
      DT::datatable(step15_venn_data_sets()$venn_data_set, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "venn_data_set"), 
                                     list(extend = 'excel', filename = "venn_data_set"))
                    ))
    })
    
    # save the venn data set
    output$save_step15_venn_set <- downloadHandler(
      filename = function() {
        paste0(input$out_step15_venn_set_title, "-merged-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step15_venn_data_sets()$venn_data_set, file = file, rowNames = FALSE, colNames = TRUE)
      }
    )
    
  })
 
  
  ### 15.2 draw the venn diagram ##############################
  step15_venn_diagram <- reactive({
    
    req(input$act_step15_venn_diagram)
    
    if (is.null(step15_venn_data_sets())) {return(NULL)}
    
    # browser()

    color_list <- c(input$in_step15_venn_set_1_color, 
                    input$in_step15_venn_set_2_color,
                    input$in_step15_venn_set_3_color,
                    input$in_step15_venn_set_4_color)

    # draw the venn plot
    venn_diagram <- venn.diagram(
      step15_venn_data_sets()$set_list,
      filename = NULL,
      output = TRUE,
      
      scaled = input$in_step15_venn_scale,
      alpha = input$in_step15_venn_alpha,
      
      # circle
      lwd = input$in_step15_venn_line_width,
      lty = 1,
      col = "black",
      fill = color_list[1:length(step15_venn_data_sets()$set_list)],
      
      # number
      cex = input$in_step15_venn_number_size,
      fontface = input$in_step15_venn_fontface,
      fontfamily = input$in_step15_venn_fontfamily,
      
      # category
      cat.cex = input$in_step15_venn_font_size,
      cat.col = "black",
      cat.fontface = input$in_step15_venn_fontface,
      cat.fontfamily = input$in_step15_venn_fontfamily,
      cat.default.pos = "outer"
      # cat.pos = c(-1, 1, 0)
      # cat.dist = c(0.2, 0.1, 0.1),
      # rotation = 1

    )
    
    return(venn_diagram)
    
  })
  
  ## output the figure
  observeEvent(input$act_step15_venn_diagram, {
    
    ## output the figure
    output$out_step15_venn_diagram <- renderPlot(
      width = input$in_step15_venn_width * 100,
      height = input$in_step15_venn_height * 100,
      {grid.draw(step15_venn_diagram())})

    # save the venn plot
    output$save_step15_venn_diagram <- downloadHandler(
      filename = function() {
        paste0(input$out_step15_venn_set_title, "-venn-diagram-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        pdf(file, width = input$in_step15_venn_width, height = input$in_step15_venn_height)
        grid.draw(step15_venn_diagram())
        dev.off()
        # ggsave(plot = step15_venn_diagram(), filename = file,
        #        width = input$in_step15_venn_width, height = input$in_step15_venn_height)
      }
    )
    
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ### 15.1 import the degs data ##############################
  step15_degs_venn_degs <- reactive({
    
    req(input$act_step15_degs_venn_import)
    
    if (is.null(input$in_step15_degs_venn_rna)) {return(NULL)}
    if (is.null(input$in_step15_degs_venn_ribo)) {return(NULL)}
    
    # browser()
    rna_degs <- NULL
    ribo_degs <- NULL
    te_degs <- NULL
    
    gene_column <- input$in_step15_degs_venn_gene
    degs_column <- input$in_step15_degs_venn_degs
    
    # import the raw data table
    if (!is.null(input$in_step15_degs_venn_rna)) {
      rna_degs <- read.xlsx(input$in_step15_degs_venn_rna$datapath,
                            sheet = input$in_step15_degs_venn_rna_sheet,
                            colNames = T, rowNames = F) %>%
        dplyr::select(!!sym(gene_column), 
                      !!sym(degs_column)) %>%
        na.omit() %>%
        dplyr::mutate(Type = paste(input$in_step15_degs_venn_set_1,
                                   !!sym(degs_column), sep = "_")) %>%
        dplyr::select(-!!sym(degs_column))
    }
    
    if (!is.null(input$in_step15_degs_venn_ribo)) {
      ribo_degs <- read.xlsx(input$in_step15_degs_venn_ribo$datapath,
                             sheet = input$in_step15_degs_venn_ribo_sheet,
                             colNames = T, rowNames = F) %>%
        dplyr::select(!!sym(gene_column), 
                      !!sym(degs_column)) %>%
        na.omit() %>%
        dplyr::mutate(Type = paste(input$in_step15_degs_venn_set_2,
                                   !!sym(degs_column), sep = "_")) %>%
        dplyr::select(-!!sym(degs_column))
    }
    
    if (!is.null(input$in_step15_degs_venn_te)) {
      te_degs <- read.xlsx(input$in_step15_degs_venn_te$datapath,
                           sheet = input$in_step15_degs_venn_te_sheet,
                           colNames = T, rowNames = F) %>%
        dplyr::select(!!sym(gene_column),
                      !!sym(degs_column)) %>%
        na.omit() %>%
        dplyr::mutate(Type = paste(input$in_step15_degs_venn_set_3,
                                   !!sym(degs_column), sep = "_")) %>%
        dplyr::select(-!!sym(degs_column))
    }
    
    # browser()
    venn_degs <- rbind(rna_degs, ribo_degs, te_degs) %>%
      dplyr::filter(!str_detect(Type, "NS"))
    
    return(venn_degs)
  })
  
  ## output the table
  observeEvent(input$act_step15_degs_venn_import, {
    
    ## output the table
    output$out_step15_degs_venn_degs <- DT::renderDataTable(server = T, {
      DT::datatable(step15_degs_venn_degs(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "venn_degs_filter"), 
                                     list(extend = 'excel', filename = "venn_degs_filter"))
                    ))
    })
    
  })
  
  ### 15.2 draw the degs venn ##############################
  step15_degs_venn_plot <- reactive({
    
    req(input$act_step15_draw_venn_plot)
    
    if (is.null(step15_degs_venn_degs())) {return(NULL)}
    
    source("R/draw_venn.R")
    
    # make the venn list
    down_df <- step15_degs_venn_degs() %>%
      dplyr::filter(str_detect(Type, "DOWN"))
    
    up_df <- step15_degs_venn_degs() %>%
      dplyr::filter(str_detect(Type, "UP"))
    
    down_list <- split(down_df[[1]], down_df[[2]])
    up_list <- split(up_df[[1]], up_df[[2]])
    venn_list <- split(step15_degs_venn_degs()[[1]], step15_degs_venn_degs()[[2]])
    
    # browser()
    # draw the down venn plot
    down_venn_plot <- draw_venn(
      venn_data = down_list, 
      set_color = input$in_step15_degs_venn_set_color,
      set_size = input$in_step15_degs_venn_set_size,
      label = input$in_step15_degs_venn_label,
      label_alpha = input$in_step15_degs_venn_label_alpha,
      label_geom = input$in_step15_degs_venn_label_geom,
      label_color = input$in_step15_degs_venn_label_color,
      label_size = input$in_step15_degs_venn_label_size,
      label_txtWidth = input$in_step15_degs_venn_label_width,
      edge_lty = input$in_step15_degs_venn_edge_lty,
      edge_size = input$in_step15_degs_venn_edge_size,
      low_color = input$in_step15_degs_venn_low_color,
      high_color = input$in_step15_degs_venn_high_color,
      font_size = input$in_step15_degs_venn_font_size
    )
    
    # draw the vup enn plot
    up_venn_plot <- draw_venn(
      venn_data = up_list, 
      set_color = input$in_step15_degs_venn_set_color,
      set_size = input$in_step15_degs_venn_set_size,
      label = input$in_step15_degs_venn_label,
      label_alpha = input$in_step15_degs_venn_label_alpha,
      label_geom = input$in_step15_degs_venn_label_geom,
      label_color = input$in_step15_degs_venn_label_color,
      label_size = input$in_step15_degs_venn_label_size,
      label_txtWidth = input$in_step15_degs_venn_label_width,
      edge_lty = input$in_step15_degs_venn_edge_lty,
      edge_size = input$in_step15_degs_venn_edge_size,
      low_color = input$in_step15_degs_venn_low_color,
      high_color = input$in_step15_degs_venn_high_color,
      font_size = input$in_step15_degs_venn_font_size
    )
    
    # draw the upset venn plot # venn_list
    fsize <- input$in_step15_degs_venn_text_scale
    
    venn_plot <- upset(fromList(venn_list), 
                       keep.order = TRUE,
                       sets = names(venn_list),
                       # order.by = c("freq", "degree"),
                       sets.x.label = "Set Size",
                       mainbar.y.label = "Intersection Size",
                       nsets = nsets, 
                       nintersects = input$in_step15_degs_venn_set_num,
                       mb.ratio = c(input$in_step15_degs_venn_mb_ratio, 
                                    1 - input$in_step15_degs_venn_mb_ratio), 
                       matrix.color = input$in_step15_degs_venn_dot_color, 
                       line.size = input$in_step15_degs_venn_line_size, 
                       point.size = input$in_step15_degs_venn_dot_size,
                       matrix.dot.alpha = 0, shade.color = 'grey90',
                       main.bar.color = input$in_step15_degs_venn_up_color,
                       sets.bar.color = brewer.pal(length(venn_list), 
                                                   input$in_step15_degs_venn_left_color), 
                       text.scale = c(fsize, fsize, fsize, fsize, fsize, fsize))
    
    return(list(down_venn = down_venn_plot,
                up_venn = up_venn_plot,
                venn_plot = venn_plot))
    
  })
  
  ## output the figure
  observeEvent(input$act_step15_draw_venn_plot, {
    
    ## output the figure
    output$out_step15_up_venn_plot <- renderPlot(
      width = input$out_step15_degs_venn_width * 100,
      height = input$out_step15_degs_venn_height * 100,
      {step15_degs_venn_plot()$down_venn})
    
    output$out_step15_down_venn_plot <- renderPlot(
      width = input$out_step15_degs_venn_width * 100,
      height = input$out_step15_degs_venn_height * 100,
      {step15_degs_venn_plot()$up_venn})
    
    output$out_step15_degs_venn_plot <- renderPlot(
      width = input$out_step15_degs_venn_width * 120,
      height = input$out_step15_degs_venn_height * 80,
      {step15_degs_venn_plot()$venn_plot})
    
    # save the venn plot
    output$save_step15_down_venn_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step15_degs_venn_out_name, "-down-venn-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step15_degs_venn_plot()$down_venn, filename = file,
               width = input$out_step15_degs_venn_width, height = input$out_step15_degs_venn_height)
      }
    )
    
    output$save_step15_up_venn_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step15_degs_venn_out_name, "-up-venn-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step15_degs_venn_plot()$up_venn, filename = file,
               width = input$out_step15_degs_venn_width, height = input$out_step15_degs_venn_height)
      }
    )
    
    output$save_step15_degs_venn_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step15_degs_venn_out_name, "-upset-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        export::graph2pdf(step15_degs_venn_plot()$venn_plot, file = file,
                          width = input$out_step15_degs_venn_width, height = input$out_step15_degs_venn_height)
        # ggsave(plot = step15_degs_venn_plot()$venn_plot, filename = file,
        #        width = input$out_step15_degs_venn_width, height = input$out_step15_degs_venn_height)
      }
    )
    
  })
  
  ### 15.3 update the degs list ##############################
  observeEvent(input$in_step15_merge_rna, {
    observe({
      req(input$in_step15_merge_rna)
      file_name <- paste(input$in_step15_merge_rna$name, collapse = "\n")
      re_name <- gsub(".xlsx|.XLSX", "", file_name)
      
      updateTextAreaInput(session, "in_step15_merge_rna_file", value = file_name)
      updateTextAreaInput(session, "in_step15_merge_rna_name", value = re_name)
    })
  })
  
  observeEvent(input$in_step15_merge_ribo, {
    observe({
      req(input$in_step15_merge_ribo)
      file_name <- paste(input$in_step15_merge_ribo$name, collapse = "\n")
      re_name <- gsub(".xlsx|.XLSX", "", file_name)
      
      updateTextAreaInput(session, "in_step15_merge_ribo_file", value = file_name)
      updateTextAreaInput(session, "in_step15_merge_ribo_name", value = re_name)
    })
  })
  
  observeEvent(input$in_step15_merge_te, {
    observe({
      req(input$in_step15_merge_te)
      file_name <- paste(input$in_step15_merge_te$name, collapse = "\n")
      re_name <- gsub(".xlsx|.XLSX", "", file_name)
      
      updateTextAreaInput(session, "in_step15_merge_te_file", value = file_name)
      updateTextAreaInput(session, "in_step15_merge_te_name", value = re_name)
    })
  })
  
  
  ### 15.3 import the degs list ##############################
  step15_merge_degs <- reactive({
    
    req(input$act_step15_merge_degs)
    
    if (is.null(input$in_step15_merge_rna) & is.null(input$in_step15_merge_ribo) & is.null(input$in_step15_merge_te)) {return(NULL)}
    
    rna_degs <- NULL
    ribo_degs <- NULL
    te_degs <- NULL
    
    degs_column <- gsub(',', '\n', input$in_step15_merge_degs_column)
    anno_column <- gsub(',', '\n', input$in_step15_merge_anno_column)
    
    degs_column <- unlist(str_split(degs_column, '\n'))
    anno_column <- unlist(str_split(anno_column, '\n'))
    
    selected_column <- c(degs_column, anno_column)
    
    # browser()
    # import the rna data table
    if (!is.null(input$in_step15_merge_rna)) {
      rna_name <- unlist(str_split(input$in_step15_merge_rna_name, '\n'))
      file_num <- nrow(input$in_step15_merge_rna)
      
      for (i in 1:file_num) {
        tryCatch({
          degs <- read.xlsx(input$in_step15_merge_rna$datapath[i], 
                            sheet = input$in_step15_merge_rna_sheet,
                            colNames = T, rowNames = F) %>%
            dplyr::select(selected_column) %>%
            dplyr::mutate(Groups = "rna_DEGs",
                          Files = rna_name[i])
          
        }, error = function(e) {
          message("Error: ", e)
          degs <- NULL
        })
        
        if (is.null(rna_degs)) {
          rna_degs <- degs
        } else {
          rna_degs <- rbind(rna_degs, degs)
        }
      }
    }
    
    # import the ribo data table
    if (!is.null(input$in_step15_merge_ribo)) {
      ribo_name <- unlist(str_split(input$in_step15_merge_ribo_name, '\n'))
      file_num <- nrow(input$in_step15_merge_ribo)
      
      for (i in 1:file_num) {
        tryCatch({
          degs <- read.xlsx(input$in_step15_merge_ribo$datapath[i],
                            sheet = input$in_step15_merge_ribo_sheet,
                            colNames = T, rowNames = F) %>%
            dplyr::select(selected_column) %>%
            dplyr::mutate(DEGs = DEGs,
                          Groups = "ribo_DEGs",
                          Files = ribo_name[i])
          
        }, error = function(e) {
          message("Error: ", e)
          degs <- NULL
        })
        
        if (is.null(ribo_degs)) {
          ribo_degs <- degs
        } else {
          ribo_degs <- rbind(ribo_degs, degs)
        }
      }
    }
    
    # import the te data table
    if (!is.null(input$in_step15_merge_te)) {
      te_name <- unlist(str_split(input$in_step15_merge_te_name, '\n'))
      file_num <- nrow(input$in_step15_merge_te)
      
      for (i in 1:file_num) {
        tryCatch({
          degs <- read.xlsx(input$in_step15_merge_te$datapath[i],
                            sheet = input$in_step15_merge_te_sheet,
                            colNames = T, rowNames = F) %>%
            dplyr::select(selected_column) %>%
            dplyr::mutate(DEGs = DEGs,
                          Groups = "DTEGs",
                          Files = te_name[i])
          
        }, error = function(e) {
          message("Error: ", e)
          degs <- NULL
        })
        
        if (is.null(te_degs)) {
          te_degs <- degs
        } else {
          te_degs <- rbind(te_degs, degs)
        }
      }
    }
    
    # browser()
    if (isTRUE(input$in_step15_merge_degs_sig)){
      merge_degs <- rbind(rna_degs, ribo_degs, te_degs) %>%
        dplyr::filter(DEGs != "NS") %>%
        dplyr::filter(DEGs != "ns")
      
    } else {
      merge_degs <- rbind(rna_degs, ribo_degs, te_degs)
    }
    
    return(merge_degs)
  })
  
  ## output the table
  observeEvent(input$act_step15_merge_degs, {
    
    ## output the table
    output$out_step15_merge_degs <- DT::renderDataTable(server = T, {
      DT::datatable(step15_merge_degs(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "merge_degs_filter"), 
                                     list(extend = 'excel', filename = "merge_degs_filter"))
                    ))
    })
    
    output$save_step15_merge_degs <- downloadHandler(
      filename = function() {
        paste0(input$out_step15_merge_out_name, "-DEGs-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step15_merge_degs(), file = file, rowNames = F, sheetName = "DEGs")
      }
    )
    
  })
  
  ### 15.4 import the degs ##############################
  step15_import_merge_degs <- reactive({
    
    req(input$act_step15_import_merge_degs)
    
    if (is.null(input$in_step15_merge_degs)) {return(NULL)}
    
    # browser()
    # import the merge degs table
    merge_degs <- read.xlsx(input$in_step15_merge_degs$datapath,
                            sheet = input$in_step15_merge_degs_sheet,
                            colNames = T, rowNames = F)
    
    return(merge_degs)
  })
  
  ## output the table
  observeEvent(input$act_step15_import_merge_degs, {
    
    ## output the table
    output$out_step15_import_degs <- DT::renderDataTable(server = T, {
      DT::datatable(step15_import_merge_degs(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "merge_degs_filter"), 
                                     list(extend = 'excel', filename = "merge_degs_filter"))
                    ))
    })
    
  })
  
  
  ### 15.5 draw the degs barplot ##############################
  step15_degs_barplot <- reactive({
    
    req(input$act_step15_draw_degs)
    
    if (is.null(step15_import_merge_degs())) {return(NULL)}
    
    # browser()
    source("R/draw_DEGs.R")
    
    # make the degs count
    degs_count <- step15_import_merge_degs() %>%
      dplyr::group_by(Files, Groups, DEGs) %>%
      dplyr::reframe(Count = n()) %>%
      ungroup() %>% 
      as.data.frame()
    
    # draw the degs bar plot
    degs_barplot <- draw_degs_bar(degs = degs_count, 
                                  x = input$in_step15_bar_x,
                                  y = input$in_step15_bar_y,
                                  fill_degs = input$in_step15_bar_degs,
                                  
                                  up_color = input$in_step15_bar_up_color,
                                  down_color = input$in_step15_bar_down_color,
                                  alpha = input$in_step15_bar_alpha,
                                  edge_color = 'white',
                                  
                                  bar_width = input$in_step15_bar_width,
                                  font_size = input$in_step15_bar_font_size,
                                  
                                  label = input$in_step15_bar_label,
                                  label_color = 'black',
                                  label_size = input$in_step15_bar_label_size,
                                  
                                  bar_degs = input$in_step15_bar_degs,
                                  bar_group = input$in_step15_bar_group,
                                  
                                  title = "",
                                  xlab = "Groups",
                                  ylab = "Count")
    
    return(degs_barplot)
  })
  
  ## output the figure
  observeEvent(input$act_step15_draw_degs, {
    
    ## output the figure
    output$out_step15_degs_bar <- renderPlot(
      width = input$out_step15_bar_width * 100,
      height = input$out_step15_bar_height * 100,
      {step15_degs_barplot()})
    
    output$save_step15_degs_bar <- downloadHandler(
      filename = function() {
        paste0(input$out_step15_bar_out_name, "-DEGs-bar-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step15_degs_barplot(), filename = file,
               width = input$out_step15_bar_width, height = input$out_step15_bar_height)
      }
    )
    
  })
  
  
  ############################################################
  ## 16 run the GO enrichment #################
  
  ### 16.1 import the degs data ##############################
  step16_go_degs <- reactive({
    
    req(input$act_step16_import_go_degs)
    
    if (is.null(input$in_step16_go_degs)) {return(NULL)}
    
    # browser()
    # import the raw data table
    degs <- tryCatch(
      {
        read.xlsx(xlsxFile = input$in_step16_go_degs$datapath, 
                  sheet = input$in_step16_go_degs_sheet, rowNames = FALSE) %>% 
          dplyr::filter(!str_detect(DEGs, "NS")) %>% 
          dplyr::filter(!str_detect(DEGs, "ns"))
      },
      error = function(e) {
        message("An error occurred while reading the Excel file: \n", e)
        return(NULL)
      }
    )
    
    return(degs)
  })
  
  ## output the table
  observeEvent(input$act_step16_import_go_degs, {
    
    ## output the table
    output$out_step16_go_degs <- DT::renderDataTable(server = T, {
      DT::datatable(step16_go_degs(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "degs_filter"), 
                                     list(extend = 'excel', filename = "degs_filter"))
                    ))
    })
    
    # show the installed packages
    select_choices <- installed.packages() %>%
      as.data.frame() %>% 
      dplyr::filter(str_detect(Package, 'org.')) 
    select_choices <- select_choices$Package
    
    updateSelectizeInput(session, "in_step16_go_orgdb", choices = select_choices, selected = select_choices[1], server = TRUE)
    
  })
  
  
  ### 16.2 run the go enrichment ##############################
  step16_go_enrichment <- reactive({
    
    req(input$act_step16_go_enrich)
    
    if (is.null(step16_go_degs())) {return(NULL)}
    
    # browser()
    source("R/run_GO_Terms.R")
    
    # set the degs list
    degs <- step16_go_degs() %>% 
      dplyr::select(input$in_step16_go_gene_column,
                    input$in_step16_go_degs_column) %>% 
      na.omit()
    
    if (input$in_step16_go_degs_group) {
      degs_list <- split(degs[[input$in_step16_go_gene_column]], 
                         degs[[input$in_step16_go_degs_column]])
    } else {
      degs_list <- degs[[input$in_step16_go_gene_column]]
    }
    
    # run the go enrichment
    go_terms <- run_go_terms(degs = degs_list,
                             degs_group = input$in_step16_go_degs_group,
                             database = input$in_step16_go_database,
                             orgdb = input$in_step16_go_orgdb,
                             gson_file = input$in_step16_go_gson,
                             keytype = input$in_step16_go_key_type,
                             ontology = "ALL",
                             pAdjustMethod = input$in_step16_go_pAdjustMethod,
                             pvalueCutoff = input$in_step16_go_pvalue,
                             qvalueCutoff = input$in_step16_go_qvalue,
                             readable = FALSE)
    
    return(go_terms)
  })
  
  ## output the table
  observeEvent(input$act_step16_go_enrich, {
    
    ## output the table
    output$out_step16_go_results <- DT::renderDataTable(server = T, {
      DT::datatable(step16_go_enrichment() %>% as.data.frame(), 
                    rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "go_enrichment"), 
                                     list(extend = 'excel', filename = "go_enrichment"))
                    ))
    })
    
    # download the go enrichment results
    output$save_step16_go_results <- downloadHandler(
      filename = function() {
        paste0(input$out_step16_go_out_name, "-GO-enrichment-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step16_go_enrichment(), file = file, rowNames = FALSE)
      }
    )
    
  })
  
  ### 16.3 simplify the go enrichment ##############################
  step16_go_enrichment_simplify <- reactive({
    
    req(input$act_step16_simplify_go_enrich)
    
    if (is.null(step16_go_enrichment())) {return(NULL)}
    
    # simplify the go enrichment
    # browser()
    
    go_terms_simplify <- clusterProfiler::simplify(
      x = step16_go_enrichment(),
      cutoff = input$in_step16_go_simplify_cutoff,
      by = input$in_step16_go_simplify_by,
      select_fun = min
    )
    
    return(go_terms_simplify)
  })
  
  ## output the table
  observeEvent(input$act_step16_simplify_go_enrich, {
    
    ## output the table
    output$out_step16_simplify_go_results <- DT::renderDataTable(server = T, {
      DT::datatable(step16_go_enrichment_simplify() %>% as.data.frame(), 
                    rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "simplify_go_enrichment"), 
                                     list(extend = 'excel', filename = "simplify_go_enrichment"))
                    ))
    })
    
    # download the go enrichment results
    output$save_step16_simplify_go_enrich_results <- downloadHandler(
      filename = function() {
        paste0(input$out_step16_go_out_name, "-simplify-GO-enrichment-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step16_go_enrichment_simplify(), file = file, rowNames = FALSE)
      }
    )
    
  })
  
  ### 16.4 draw the go enrichment ##############################
  step16_go_enrichment_plot <- reactive({
    
    req(input$act_step16_draw_go_enrich_plot)
    
    # browser()
    
    if (input$act_step16_simplify_go_enrich != 0 && !is.null(step16_go_enrichment_simplify())) {
      go_results <- step16_go_enrichment_simplify() %>% 
        dplyr::filter(pvalue < input$in_step16_go_pvalue_flt,
                      qvalue < input$in_step16_go_qvalue_flt,
                      Count >= input$in_step16_go_count_flt)
      
    } else if (input$act_step16_go_enrich != 0 && !is.null(step16_go_enrichment())) {
      go_results <- step16_go_enrichment() %>% 
        dplyr::filter(pvalue < input$in_step16_go_pvalue_flt,
                      qvalue < input$in_step16_go_qvalue_flt,
                      Count >= input$in_step16_go_count_flt)
      
    } else {
      return(NULL)
    }
    
    source("R/draw_GO_Terms.R")
    
    # browser()
    
    # draw the go enrichment
    go_enrich_plot <- draw_go_terms(go_results = go_results,
                                    x = "GeneRatio",
                                    color = "p.adjust",
                                    
                                    showCategory = input$in_step16_go_item_num,
                                    
                                    dotsize = input$in_step16_go_dot_size,
                                    chr_width = input$in_step16_go_item_width,
                                    
                                    dot_alpha = input$in_step16_go_dot_alpha,
                                    low_color = input$in_step16_go_low_color,
                                    high_color = input$in_step16_go_high_color,
                                    
                                    ontology = input$in_step16_go_ontology,
                                    
                                    facet_group1 = input$in_step16_go_facet_group1,
                                    facet_group2 = input$in_step16_go_facet_group2,
                                    
                                    gridwidth = input$in_step16_go_grid_width,
                                    font_size = input$in_step16_go_font_size,
                                    
                                    xmin = input$in_step16_go_xmin,
                                    xmax = input$in_step16_go_xmax,
                                    x_breaks = input$in_step16_go_break,
                                    
                                    clustered = input$in_step16_go_degs_group,
                                    
                                    title = input$in_step16_go_title,
                                    xlab = input$in_step16_go_xlab,
                                    ylab = input$in_step16_go_ylab)
    
    return(go_enrich_plot)
  })
  
  ## output the plot
  observeEvent(input$act_step16_draw_go_enrich_plot, {
    
    ## output the figure
    output$out_step16_go_dot_plot <- renderPlot(
      width = input$out_step16_go_width * 100,
      height = input$out_step16_go_height * 100,
      {step16_go_enrichment_plot()})
    
    output$save_step16_go_dot_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step16_go_out_name, "-GO-enrich-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step16_go_enrichment_plot(), filename = file,
               width = input$out_step16_go_width, height = input$out_step16_go_height)
      }
    )
  })
  
  ### 16.5 import the merged DEGs ##############################
  step16_group_go_degs <- reactive({
    
    req(input$act_step16_import_group_go_degs)
    
    if (is.null(input$in_step16_group_go_degs)) {return(NULL)}
    
    # browser()
    # import the merge degs table
    merge_degs <- tryCatch(
      {
        read.xlsx(input$in_step16_group_go_degs$datapath,
                  sheet = input$in_step16_group_go_degs_sheet,
                  colNames = T, rowNames = F)
      },
      error = function(e) {
        message("An error occurred while reading the Excel file: \n", e)
        return(NULL)
      }
    )
    
    return(merge_degs)
  })
  
  ## output the table
  observeEvent(input$act_step16_import_group_go_degs, {
    
    ## output the table
    output$out_step16_group_go_degs <- DT::renderDataTable(server = T, {
      DT::datatable(step16_group_go_degs(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "merge_degs_filter"), 
                                     list(extend = 'excel', filename = "merge_degs_filter"))
                    ))
    })
    
    # show the installed packages
    select_choices <- installed.packages() %>%
      as.data.frame() %>% 
      dplyr::filter(str_detect(Package, 'org.')) 
    select_choices <- select_choices$Package
    
    updateSelectizeInput(session, "in_step16_group_go_orgdb", choices = select_choices, selected = select_choices[1], server = TRUE)
    
    # show the DEGs groups
    degs_select_choices <- colnames(step16_group_go_degs())
    degs_select_choices <- grep("Gene|baseMean|log2FoldChange|pvalue|padj|GeneID|Symbol|Ensembl|Name|logFC|logCPM|LR|PValue|FDR", 
                                degs_select_choices, invert = T, value = T)
    updateSelectizeInput(session, "in_step16_group_go_group_column", choices = degs_select_choices, 
                         selected = degs_select_choices, server = TRUE)
    
  })
  
  ### 16.6 run the group go enrichment ##############################
  step16_group_go_enrichment <- reactive({
    
    req(input$act_step16_group_go_enrich)
    
    if (is.null(step16_group_go_degs())) {return(NULL)}
    
    source("R/run_Group_GO_Terms.R")
    
    # set the degs list
    # browser()
    group_column <- c(input$in_step16_group_go_gene_column,
                      input$in_step16_group_go_group_column)
    
    degs <- step16_group_go_degs() %>% 
      tidyr::drop_na(all_of(group_column))
    
    
    # run the go enrichment
    go_terms <- run_group_go_terms(degs = degs %>% remove_rownames(),
                                   degs_id = input$in_step16_group_go_gene_column,
                                   degs_group = input$in_step16_group_go_group_column,
                                   database = input$in_step16_group_go_database,
                                   gson_file = input$in_step16_group_go_gson,
                                   orgdb = input$in_step16_group_go_orgdb,
                                   keytype = input$in_step16_group_go_key_type,
                                   ontology = "ALL",
                                   pAdjustMethod = input$in_step16_group_go_pAdjustMethod,
                                   pvalueCutoff = input$in_step16_group_go_pvalue,
                                   qvalueCutoff = input$in_step16_group_go_qvalue)
    
    return(go_terms)
  })
  
  ## output the table
  observeEvent(input$act_step16_group_go_enrich, {
    
    ## output the table
    output$out_step16_group_go_results <- DT::renderDataTable(server = T, {
      DT::datatable(step16_group_go_enrichment() %>% as.data.frame(), 
                    rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "group_go_enrichment"), 
                                     list(extend = 'excel', filename = "group_go_enrichment"))
                    ))
    })
    
    # download the go enrichment results
    output$save_step16_group_go_enrich_results <- downloadHandler(
      filename = function() {
        paste0(input$out_step16_group_go_out_name, "-group-GO-enrichment-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step16_group_go_enrichment(), file = file, rowNames = FALSE)
      }
    )
    
    # update the select input
    go_select_choices <- c("ONTOLOGY", "Cluster", input$in_step16_group_go_group_column)
    updateSelectizeInput(session, "in_step16_group_go_facet_group1", choices = go_select_choices, 
                         selected = NULL, server = TRUE)
    updateSelectizeInput(session, "in_step16_group_go_facet_group2", choices = go_select_choices, 
                         selected = NULL, server = TRUE)
  })
  
  
  ### 16.7 simplify the go enrichment ##############################
  step16_group_go_enrichment_simplify <- reactive({
    
    req(input$act_step16_simplify_group_go_enrich)
    
    if (is.null(step16_group_go_enrichment())) {return(NULL)}
    
    # simplify the go enrichment
    # browser()
    
    go_terms_simplify <- clusterProfiler::simplify(
      x = step16_group_go_enrichment(),
      cutoff = input$in_step16_group_go_simplify_cutoff,
      by = input$in_step16_group_go_simplify_by,
      select_fun = min
    )
    
    return(go_terms_simplify)
  })
  
  ## output the table
  observeEvent(input$act_step16_simplify_group_go_enrich, {
    
    ## output the table
    output$out_step16_simplify_group_go_results <- DT::renderDataTable(server = T, {
      DT::datatable(step16_group_go_enrichment_simplify() %>% as.data.frame(), 
                    rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "simplify_go_enrichment"), 
                                     list(extend = 'excel', filename = "simplify_go_enrichment"))
                    ))
    })
    
    # download the go enrichment results
    output$save_step16_simplify_group_go_enrich_results <- downloadHandler(
      filename = function() {
        paste0(input$out_step16_group_go_out_name, "-simplify-group-GO-enrichment-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step16_group_go_enrichment_simplify(), file = file, rowNames = FALSE)
      }
    )
    
  })
  
  ### 16.8 draw the group go enrichment ##############################
  step16_group_go_enrichment_plot <- reactive({
    
    req(input$act_step16_draw_group_go_enrich_plot)
    
    if (input$act_step16_simplify_group_go_enrich != 0 && !is.null(step16_group_go_enrichment_simplify())) {
      go_results <- step16_group_go_enrichment_simplify() %>% 
        dplyr::filter(pvalue < input$in_step16_group_go_pvalue_flt,
                      qvalue < input$in_step16_group_go_qvalue_flt,
                      Count >= input$in_step16_group_go_count_flt)
      
    } else if (input$act_step16_group_go_enrich != 0 && !is.null(step16_group_go_enrichment())) {
      go_results <- step16_group_go_enrichment() %>% 
        dplyr::filter(pvalue < input$in_step16_group_go_pvalue_flt,
                      qvalue < input$in_step16_group_go_qvalue_flt,
                      Count >= input$in_step16_group_go_count_flt)
      
    } else {
      return(NULL)
    }
    
    source("R/draw_GO_Terms.R")
    
    # browser()
    
    # draw the go enrichment
    go_enrich_plot <- draw_go_terms(go_results = go_results,
                                    x = "GeneRatio",
                                    color = "p.adjust",
                                    
                                    showCategory = input$in_step16_group_go_item_num,
                                    
                                    dotsize = input$in_step16_group_go_dot_size,
                                    chr_width = input$in_step16_group_go_item_width,
                                    
                                    dot_alpha = input$in_step16_group_go_dot_alpha,
                                    low_color = input$in_step16_group_go_low_color,
                                    high_color = input$in_step16_group_go_high_color,
                                    
                                    ontology = input$in_step16_group_go_ontology,
                                    
                                    facet_group1 = input$in_step16_group_go_facet_group1,
                                    facet_group2 = input$in_step16_group_go_facet_group2,
                                    
                                    gridwidth = input$in_step16_group_go_grid_width,
                                    font_size = input$in_step16_group_go_font_size,
                                    
                                    clustered = FALSE,
                                    
                                    title = input$in_step16_group_go_title,
                                    xlab = input$in_step16_group_go_xlab,
                                    ylab = input$in_step16_group_go_ylab)
    
    return(go_enrich_plot)
  })
  
  ## output the plot
  observeEvent(input$act_step16_draw_group_go_enrich_plot, {
    
    ## output the figure
    output$out_step16_group_go_plot <- renderPlot(
      width = input$out_step16_group_go_width * 100,
      height = input$out_step16_group_go_height * 100,
      {step16_group_go_enrichment_plot()})
    
    output$save_step16_group_go_enrich_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step16_group_go_out_name, "-group-GO-enrich-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step16_group_go_enrichment_plot(), filename = file,
               width = input$out_step16_group_go_width, height = input$out_step16_group_go_height)
      }
    )
  })
  
  ############################################################
  ## 17 run the GO GSEA enrichment #################
  
  ### 17.1 import the degs data ##############################
  step17_go_gsea_degs <- reactive({
    
    req(input$act_step17_import_go_gsea_degs)
    
    if (is.null(input$in_step17_go_gsea_degs)) {return(NULL)}
    
    # browser()
    # import the raw data table
    
    degs <- tryCatch({
      read.xlsx(xlsxFile = input$in_step17_go_gsea_degs$datapath, 
                sheet = input$in_step17_go_gsea_degs_sheet, rowNames = FALSE)
    },
    error = function(e) {
      message("An error occurred while reading the Excel file: \n", e)
      return(NULL)
    }
    )
    
    return(degs)
  })
  
  ## output the table
  observeEvent(input$act_step17_import_go_gsea_degs, {
    
    ## output the table
    output$out_step17_go_gsea_degs <- DT::renderDataTable(server = T, {
      DT::datatable(step17_go_gsea_degs(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "degs_filter"), 
                                     list(extend = 'excel', filename = "degs_filter"))
                    ))
    })
    
    # show the installed packages
    select_choices <- installed.packages() %>%
      as.data.frame() %>% 
      dplyr::filter(str_detect(Package, 'org.')) 
    select_choices <- select_choices$Package
    
    updateSelectizeInput(session, "in_step17_go_gsea_orgdb", choices = select_choices, selected = select_choices[1], server = TRUE)
  })
  
  
  ### 17.2 create the genes list ##############################
  step17_go_gsea_degs_vector <- reactive({
    
    req(input$act_step17_create_go_gsea_list)
    
    if (is.null(step17_go_gsea_degs())) {return(NULL)}
    
    # browser()
    
    # make vector for gsea enrichment
    
    degs_log2fc <- step17_go_gsea_degs() %>% 
      tidyr::drop_na(!!sym(input$in_step17_go_gsea_gene_column)) %>%
      dplyr::arrange(desc(!!sym(input$in_step17_go_gsea_log2fc_column))) %>% 
      dplyr::distinct(!!sym(input$in_step17_go_gsea_gene_column), .keep_all = TRUE) %>%
      dplyr::pull(input$in_step17_go_gsea_log2fc_column, input$in_step17_go_gsea_gene_column)
    
    return(degs_log2fc)
  })
  
  ## output the table
  observeEvent(input$act_step17_create_go_gsea_list, {
    
    ## output the table
    output$out_step17_go_gsea_gene_list <- renderPrint({
      print("The top 10 genes for GO GSEA enrichment are:")
      step17_go_gsea_degs_vector() %>% head(10)
    })
    
  })
  
  ### 17.3 run the gsea enrichment ############################
  step17_go_gsea_enrichment <- reactive({
    
    req(input$act_step17_go_gsea_enrich)
    
    if (is.null(step17_go_gsea_degs_vector())) {return(NULL)}
    
    # browser()
    source("R/run_GO_GSEA_Terms.R")
    
    # run the go enrichment
    go_gsea_terms <- run_go_gsea_terms(degs = step17_go_gsea_degs_vector(),
                                       database = input$in_step17_go_gsea_database,
                                       orgdb = input$in_step17_go_gsea_orgdb,
                                       gson_file = input$in_step17_go_gsea_gson,
                                       keytype = input$in_step17_go_gsea_key_type,
                                       ontology = "ALL",
                                       eps = input$in_step17_go_gsea_eps,
                                       minGSSize = input$in_step17_go_gsea_minsize,
                                       maxGSSize = input$in_step17_go_gsea_maxsize,
                                       pAdjustMethod = input$in_step17_go_gsea_pAdjustMethod,
                                       pvalueCutoff = input$in_step17_go_gsea_pvalue)
    # browser()
    
    return(go_gsea_terms)
  })
  
  ## output the table
  observeEvent(input$act_step17_go_gsea_enrich, {
    
    ## output the table
    output$out_step17_go_gsea_results <- DT::renderDataTable(server = T, {
      DT::datatable(step17_go_gsea_enrichment() %>% as.data.frame(), 
                    rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "go_gsea_enrichment"), 
                                     list(extend = 'excel', filename = "go_gsea_enrichment"))
                    ))
    })
    
    # download the go enrichment results
    output$save_step17_go_gsea_results <- downloadHandler(
      filename = function() {
        paste0(input$out_step17_go_gsea_out_name, "-GO-GSEA-enrichment-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step17_go_gsea_enrichment(), file = file, rowNames = FALSE)
      }
    )
    
  })
  
  ### 17.4 draw the go gsea dotplot ##############################
  step17_go_gsea_enrichment_plot <- reactive({
    
    req(input$act_step17_draw_go_gsea_dot_plot)
    
    if (is.null(step17_go_gsea_enrichment())) {
      return(NULL)
    } else {
      go_gsea_results <- step17_go_gsea_enrichment() %>% 
        dplyr::filter(pvalue < input$in_step17_go_gsea_pvalue_flt,
                      qvalue < input$in_step17_go_gsea_qvalue_flt)
    }
    
    source("R/draw_GO_Terms.R")
    
    # browser()
    
    # draw the go gsea enrichment plot
    go_gsea_dot_plot <- draw_go_terms(go_results = go_gsea_results,
                                      x = "GeneRatio",
                                      color = "p.adjust",
                                      
                                      ontology = input$in_step17_go_gsea_ontology,
                                      
                                      showCategory = input$in_step17_go_gsea_item_num,
                                      
                                      dotsize = input$in_step17_go_gsea_dot_size,
                                      chr_width = input$in_step17_go_gsea_item_width,
                                      
                                      dot_alpha = input$in_step17_go_gsea_dot_alpha,
                                      low_color = input$in_step17_go_gsea_low_color,
                                      high_color = input$in_step17_go_gsea_high_color,
                                      
                                      facet_group1 = input$in_step17_go_gsea_facet_group1,
                                      facet_group2 = input$in_step17_go_gsea_facet_group2,
                                      
                                      gridwidth = input$in_step17_go_gsea_grid_width,
                                      font_size = input$in_step17_go_gsea_font_size,
                                      
                                      xmin = input$in_step17_go_gsea_xmin,
                                      xmax = input$in_step17_go_gsea_xmax,
                                      x_breaks = input$in_step17_go_gsea_break,
                                      
                                      clustered = FALSE,
                                      
                                      title = input$in_step17_go_gsea_title,
                                      xlab = input$in_step17_go_gsea_xlab,
                                      ylab = input$in_step17_go_gsea_ylab)
    
    return(go_gsea_dot_plot)
  })
  
  ## output the plot
  observeEvent(input$act_step17_draw_go_gsea_dot_plot, {
    
    ## output the figure
    output$out_step17_go_gsea_dot_plot <- renderPlot(
      width = input$out_step17_go_gsea_width * 100,
      height = input$out_step17_go_gsea_height * 100,
      {step17_go_gsea_enrichment_plot()})
    
    output$save_step17_go_gsea_dot_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step17_go_gsea_out_name, "-GO-GSEA-dotplot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step17_go_gsea_enrichment_plot(), filename = file,
               width = input$out_step17_go_gsea_width, height = input$out_step17_go_gsea_height)
      }
    )
  })
  
  
  ### 17.5 draw the go gseaplot ##############################
  step17_go_gsea_line_plot <- reactive({
    
    req(input$act_step17_draw_go_gsea_lineplot)
    
    if (is.null(step17_go_gsea_enrichment())) {
      return(NULL)
    } else {
      go_gsea_results <- step17_go_gsea_enrichment() %>% 
        dplyr::filter(pvalue < input$in_step17_go_gsea_pvalue_flt,
                      qvalue < input$in_step17_go_gsea_qvalue_flt)
    }
    
    source("R/draw_GSEA_Terms.R")
    
    # browser()
    
    geneSetID <- gsub(" ", "", input$in_step17_go_gsea_line_item_name)
    
    geneSetID <- unlist(str_split(geneSetID, ",|\n"))
    
    rel_heights <- c(input$in_step17_go_gsea_line_height1,
                     input$in_step17_go_gsea_line_height2,
                     input$in_step17_go_gsea_line_height3)
    
    subplots <- as.numeric(input$in_step17_go_gsea_line_subplot)
    
    # draw the go gsea enrichment plot
    go_gsea_line_plot <- draw_gsea_terms(gsea_results = go_gsea_results,
                                         geneSetID = geneSetID,
                                         title = input$in_step17_go_gsea_line_title,
                                         base_size = input$in_step17_go_gsea_line_fontsize,
                                         rel_heights = rel_heights,
                                         subplots = subplots,
                                         color = input$in_step17_go_gsea_line_color,
                                         pvalue_table = input$in_step17_go_gsea_line_pvalue,
                                         ES_geom = input$in_step17_go_gsea_line_es_geom)
    
    return(go_gsea_line_plot)
  })
  
  ## output the plot
  observeEvent(input$act_step17_draw_go_gsea_lineplot, {
    
    ## output the figure
    output$out_step17_go_gsea_line_plot <- renderPlot(
      width = input$out_step17_go_gsea_line_width * 100,
      height = input$out_step17_go_gsea_line_height * 100,
      {step17_go_gsea_line_plot()})
    
    output$save_step17_go_gsea_lineplot <- downloadHandler(
      filename = function() {
        paste0(input$out_step17_go_gsea_line_out_name, "-GO-gseaplot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step17_go_gsea_line_plot(), filename = file,
               width = input$out_step17_go_gsea_line_width, height = input$out_step17_go_gsea_line_height)
      }
    )
  })
  
  ############################################################
  ## 18 run the KEGG enrichment #################
  
  ### 18.1 import the degs data ##############################
  step18_kegg_degs <- reactive({
    
    req(input$act_step18_import_kegg_degs)
    
    if (is.null(input$in_step18_kegg_degs)) {return(NULL)}
    
    # browser()
    # import the raw data table
    
    degs <- tryCatch(
      {
        read.xlsx(xlsxFile = input$in_step18_kegg_degs$datapath, 
                  sheet = input$in_step18_kegg_degs_sheet, rowNames = FALSE) %>% 
          dplyr::filter(!str_detect(DEGs, "NS")) %>% 
          dplyr::filter(!str_detect(DEGs, "ns"))
      },
      error = function(e) {
        message("An error occurred while reading the Excel file: \n", e)
        return(NULL)
      }
    )
    
    return(degs)
  })
  
  ## output the table
  observeEvent(input$act_step18_import_kegg_degs, {
    
    ## output the table
    output$out_step18_kegg_degs <- DT::renderDataTable(server = T, {
      DT::datatable(step18_kegg_degs(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "degs_filter"), 
                                     list(extend = 'excel', filename = "degs_filter"))
                    ))
    })
    
  })
  
  
  ### 18.2 run the kegg enrichment ##############################
  step18_kegg_enrichment <- reactive({
    
    req(input$act_step18_kegg_enrich)
    
    if (is.null(step18_kegg_degs())) {return(NULL)}
    
    source("R/run_KEGG_Terms.R")
    
    # set the degs list
    degs <- step18_kegg_degs() %>% 
      dplyr::select(input$in_step18_kegg_gene_column, DEGs) %>% 
      na.omit()
    
    if (input$in_step18_kegg_degs_group) {
      degs_list <- split(degs[[input$in_step18_kegg_gene_column]], degs$DEGs)
    } else {
      degs_list <- degs[[input$in_step18_kegg_gene_column]]
    }
    
    # run the go enrichment
    # browser()
    kegg_terms <- run_kegg_terms(degs = degs_list,
                                 degs_group = input$in_step18_kegg_degs_group,
                                 database = input$in_step18_kegg_database,
                                 species = input$in_step18_kegg_species,
                                 gson_file = input$in_step18_kegg_gson$datapath,
                                 keytype = input$in_step18_kegg_key_type,
                                 pAdjustMethod = input$in_step18_kegg_pAdjustMethod,
                                 pvalueCutoff = input$in_step18_kegg_pvalue,
                                 qvalueCutoff = input$in_step18_kegg_qvalue,
                                 use_internal_data = input$in_step18_kegg_internal)
    
    return(kegg_terms)
  })
  
  ## output the table
  observeEvent(input$act_step18_kegg_enrich, {
    
    ## output the table
    output$out_step18_kegg_results <- DT::renderDataTable(server = T, {
      DT::datatable(step18_kegg_enrichment() %>% as.data.frame(), 
                    rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "kegg_enrichment"), 
                                     list(extend = 'excel', filename = "kegg_enrichment"))
                    ))
    })
    
    # download the go enrichment results
    output$save_step18_kegg_results <- downloadHandler(
      filename = function() {
        paste0(input$out_step18_kegg_out_name, "-KEGG-enrichment-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step18_kegg_enrichment(), file = file, rowNames = FALSE)
      }
    )
    
  })
  
  
  ### 18.3 draw the kegg enrichment ##############################
  step18_kegg_enrichment_plot <- reactive({
    
    req(input$act_step18_draw_enrich_plot)
    
    if (is.null(step18_kegg_enrichment())) {
      return(NULL)
    } else {
      kegg_results <- step18_kegg_enrichment() %>% 
        dplyr::filter(pvalue < input$in_step18_kegg_pvalue_flt,
                      qvalue < input$in_step18_kegg_qvalue_flt,
                      Count >= input$in_step18_kegg_count_flt)
    }
    
    source("R/draw_KEGG_Terms.R")
    
    # browser()
    
    # draw the kegg enrichment plot
    kegg_enrich_plot <- draw_kegg_terms(kegg_results = kegg_results,
                                        x = "GeneRatio",
                                        color = "p.adjust",
                                        
                                        showCategory = input$in_step18_kegg_item_num,
                                        
                                        dotsize = input$in_step18_kegg_dot_size,
                                        chr_width = input$in_step18_kegg_item_width,
                                        
                                        dot_alpha = input$in_step18_kegg_dot_alpha,
                                        low_color = input$in_step18_kegg_low_color,
                                        high_color = input$in_step18_kegg_high_color,
                                        
                                        facet_group1 = input$in_step18_kegg_facet_group1,
                                        facet_group2 = NULL,
                                        
                                        gridwidth = input$in_step18_kegg_grid_width,
                                        font_size = input$in_step18_kegg_font_size,
                                        
                                        xmin = input$in_step18_kegg_xmin,
                                        xmax = input$in_step18_kegg_xmax,
                                        x_breaks = input$in_step18_kegg_break,
                                        
                                        clustered = input$in_step18_kegg_degs_group,
                                        
                                        title = input$in_step18_kegg_title,
                                        xlab = input$in_step18_kegg_xlab,
                                        ylab = input$in_step18_kegg_ylab)
    
    return(kegg_enrich_plot)
  })
  
  ## output the plot
  observeEvent(input$act_step18_draw_enrich_plot, {
    
    ## output the figure
    output$out_step18_kegg_dot_plot <- renderPlot(
      width = input$out_step18_kegg_width * 100,
      height = input$out_step18_kegg_height * 100,
      {step18_kegg_enrichment_plot()})
    
    output$save_step18_kegg_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step18_kegg_out_name, "-KEGG-enrich-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step18_kegg_enrichment_plot(), filename = file,
               width = input$out_step18_kegg_width, height = input$out_step18_kegg_height)
      }
    )
  })
  
  ### 18.4 import the merged DEGs ##############################
  step18_group_kegg_degs <- reactive({
    
    req(input$act_step18_import_group_kegg_degs)
    
    if (is.null(input$in_step18_group_kegg_degs)) {return(NULL)}
    
    # browser()
    # import the merge degs table
    # merge_degs <- read.xlsx(input$in_step18_group_kegg_degs$datapath,
    #                         sheet = input$in_step18_group_kegg_degs_sheet,
    #                         colNames = T, rowNames = F)
    # 
    
    merge_degs <- tryCatch(
      {
        read.xlsx(input$in_step18_group_kegg_degs$datapath,
                  sheet = input$in_step18_group_kegg_degs_sheet,
                  colNames = TRUE, rowNames = FALSE)
      },
      error = function(e) {
        message("An error occurred while reading the Excel file: \n", e)
        return(NULL)
      }
    )
    
    return(merge_degs)
  })
  
  ## output the table
  observeEvent(input$act_step18_import_group_kegg_degs, {
    
    ## output the table
    output$out_step18_group_kegg_degs <- DT::renderDataTable(server = T, {
      DT::datatable(step18_group_kegg_degs(), rownames = FALSE, filter = "top", extensions = 'Buttons', selection = 'single',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "merge_degs_filter"), 
                                     list(extend = 'excel', filename = "merge_degs_filter"))
                    ))
    })
    
    # show the DEGs groups
    degs_select_choices <- colnames(step18_group_kegg_degs())
    degs_select_choices <- grep("Gene|baseMean|log2FoldChange|pvalue|padj|GeneID|Symbol|Ensembl|Name|logFC|logCPM|LR|PValue|FDR", 
                                degs_select_choices, invert = T, value = T)
    updateSelectizeInput(session, "in_step18_group_kegg_group_column", choices = degs_select_choices, selected = degs_select_choices, server = TRUE)
    
  })
  
  ### 18.5 run the group go enrichment ##############################
  step18_group_kegg_enrichment <- reactive({
    
    req(input$act_step18_group_kegg_enrich)
    
    if (is.null(step18_group_kegg_degs())) {return(NULL)}
    
    source("R/run_Group_KEGG_Terms.R")
    
    # set the degs list
    group_column <- c(input$in_step18_group_kegg_gene_column,
                      input$in_step18_group_kegg_group_column)
    
    degs <- step18_group_kegg_degs() %>% 
      tidyr::drop_na(all_of(group_column))
    
    # browser()
    # run the group kegg enrichment
    kegg_terms <- run_group_kegg_terms(degs = degs,
                                       degs_id = input$in_step18_group_kegg_gene_column,
                                       degs_group = input$in_step18_group_kegg_group_column,
                                       database = input$in_step18_group_kegg_database,
                                       species = input$in_step18_group_kegg_species,
                                       gson_file = input$in_step18_group_kegg_gson$datapath,
                                       keytype = input$in_step18_group_kegg_key_type,
                                       pAdjustMethod = input$in_step18_group_kegg_pAdjustMethod,
                                       pvalueCutoff = input$in_step18_group_kegg_pvalue,
                                       qvalueCutoff = input$in_step18_group_kegg_qvalue,
                                       use_internal_data = input$in_step18_group_kegg_internal)
    
    # update the select input
    kegg_select_choices <- c('Cluster', input$in_step18_group_kegg_group_column)
    updateSelectizeInput(session, "in_step18_group_kegg_facet_group1", choices = kegg_select_choices, 
                         selected = "", server = TRUE)
    updateSelectizeInput(session, "in_step18_group_kegg_facet_group2", choices = kegg_select_choices, 
                         selected = "", server = TRUE)
    
    return(kegg_terms)
  })
  
  ## output the table
  observeEvent(input$act_step18_group_kegg_enrich, {
    
    ## output the table
    output$out_step18_group_kegg_results <- DT::renderDataTable(server = T, {
      DT::datatable(step18_group_kegg_enrichment() %>% as.data.frame(), 
                    rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "kegg_enrichment"), 
                                     list(extend = 'excel', filename = "kegg_enrichment"))
                    ))
    })
    
    # download the kegg enrichment results
    output$save_step18_group_kegg_enrich_results <- downloadHandler(
      filename = function() {
        paste0(input$out_step18_group_kegg_out_name, "-group-KEGG-enrichment-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step18_group_kegg_enrichment(), file = file, rowNames = FALSE)
      }
    )
    
  })
  
  ### 18.6 draw the group go enrichment ##############################
  step18_group_kegg_enrichment_plot <- reactive({
    
    req(input$act_step18_draw_group_kegg_enrich_plot)
    
    if (is.null(step18_group_kegg_enrichment())) {
      return(NULL)
    } else {
      group_kegg_results <- step18_group_kegg_enrichment() %>% 
        dplyr::filter(pvalue < input$in_step18_group_kegg_pvalue_flt,
                      qvalue < input$in_step18_group_kegg_qvalue_flt,
                      Count >= input$in_step18_group_kegg_count_flt)
    }
    
    source("R/draw_KEGG_Terms.R")
    
    # browser()
    
    # draw the group kegg enrichment
    kegg_enrich_plot <- draw_kegg_terms(kegg_results = group_kegg_results,
                                        x = "GeneRatio",
                                        color = "p.adjust",
                                        
                                        showCategory = input$in_step18_group_kegg_item_num,
                                        
                                        dotsize = input$in_step18_group_kegg_dot_size,
                                        chr_width = input$in_step18_group_kegg_item_width,
                                        
                                        dot_alpha = input$in_step18_group_kegg_dot_alpha,
                                        low_color = input$in_step18_group_kegg_low_color,
                                        high_color = input$in_step18_group_kegg_high_color,
                                        
                                        facet_group1 = input$in_step18_group_kegg_facet_group1,
                                        facet_group2 = input$in_step18_group_kegg_facet_group2,
                                        
                                        gridwidth = input$in_step18_group_kegg_grid_width,
                                        font_size = input$in_step18_group_kegg_font_size,
                                        
                                        xmin = NA,
                                        xmax = NA,
                                        
                                        clustered = FALSE,
                                        
                                        title = input$in_step18_group_kegg_title,
                                        xlab = input$in_step18_group_kegg_xlab,
                                        ylab = input$in_step18_group_kegg_ylab)
    
    return(kegg_enrich_plot)
  })
  
  ## output the plot
  observeEvent(input$act_step18_draw_group_kegg_enrich_plot, {
    
    ## output the figure
    output$out_step18_group_kegg_plot <- renderPlot(
      width = input$out_step18_group_kegg_width * 100,
      height = input$out_step18_group_kegg_height * 100,
      {step18_group_kegg_enrichment_plot()})
    
    output$save_step18_group_kegg_enrich_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step18_group_kegg_out_name, "-group-KEGG-enrich-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step18_group_kegg_enrichment_plot(), filename = file,
               width = input$out_step18_group_kegg_width, height = input$out_step18_group_kegg_height)
      }
    )
  })
  
  
  ############################################################
  ## 19 run the KEGG GSEA enrichment #################
  
  ### 19.1 import the degs data ##############################
  step19_kegg_gsea_degs <- reactive({
    
    req(input$act_step19_import_kegg_gsea_degs)
    
    if (is.null(input$in_step19_kegg_gsea_degs)) {return(NULL)}
    
    # browser()
    # import the raw data table
    
    degs <- tryCatch({
      read.xlsx(xlsxFile = input$in_step19_kegg_gsea_degs$datapath, 
                sheet = input$in_step19_kegg_gsea_degs_sheet, rowNames = FALSE)
    },
    error = function(e) {
      message("An error occurred while reading the Excel file: \n", e)
      return(NULL)
    }
    )
    
    return(degs)
  })
  
  ## output the table
  observeEvent(input$act_step19_import_kegg_gsea_degs, {
    
    ## output the table
    output$out_step19_kegg_gsea_degs <- DT::renderDataTable(server = T, {
      DT::datatable(step19_kegg_gsea_degs(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "degs_filter"), 
                                     list(extend = 'excel', filename = "degs_filter"))
                    ))
    })
    
  })
  
  
  ### 19.2 create the genes list ##############################
  step19_kegg_gsea_degs_vector <- reactive({
    
    req(input$act_step19_create_kegg_gsea_list)
    
    if (is.null(step19_kegg_gsea_degs())) {return(NULL)}
    
    # browser()
    
    # make vector for gsea enrichment
    degs_log2fc <- step19_kegg_gsea_degs() %>% 
      tidyr::drop_na(!!sym(input$in_step19_kegg_gsea_gene_column)) %>%
      dplyr::arrange(desc(!!sym(input$in_step19_kegg_gsea_log2fc_column))) %>% 
      dplyr::distinct(!!sym(input$in_step19_kegg_gsea_gene_column), .keep_all = TRUE) %>%
      dplyr::pull(input$in_step19_kegg_gsea_log2fc_column, input$in_step19_kegg_gsea_gene_column)
    
    return(degs_log2fc)
  })
  
  ## output the table
  observeEvent(input$act_step19_create_kegg_gsea_list, {
    
    ## output the table
    output$out_step19_kegg_gsea_gene_list <- renderPrint({
      print("The top 10 genes for KEGG GSEA enrichment are:")
      step19_kegg_gsea_degs_vector() %>% head(10)
    })
    
  })
  
  ### 19.3 run the gsea enrichment ############################
  step19_kegg_gsea_enrichment <- reactive({
    
    req(input$act_step19_kegg_gsea_enrich)
    
    if (is.null(step19_kegg_gsea_degs_vector())) {return(NULL)}
    
    # browser()
    source("R/run_KEGG_GSEA_Terms.R")
    
    # run the go enrichment
    kegg_gsea_terms <- run_kegg_gsea_terms(degs = step19_kegg_gsea_degs_vector(),
                                           database = input$in_step19_kegg_gsea_database,
                                           gson_file = input$in_step19_kegg_gsea_gson,
                                           species = input$in_step19_kegg_gsea_species,
                                           kegg_type = input$in_step19_kegg_gsea_type,
                                           keytype = input$in_step19_kegg_gsea_key_type,
                                           eps = input$in_step19_kegg_gsea_eps,
                                           minGSSize = input$in_step19_kegg_gsea_minsize,
                                           maxGSSize = input$in_step19_kegg_gsea_maxsize,
                                           pAdjustMethod = input$in_step19_kegg_gsea_pAdjustMethod,
                                           pvalueCutoff = input$in_step19_kegg_gsea_pvalue)
    
    return(kegg_gsea_terms)
  })
  
  ## output the table
  observeEvent(input$act_step19_kegg_gsea_enrich, {
    
    ## output the table
    output$out_step19_kegg_gsea_results <- DT::renderDataTable(server = T, {
      DT::datatable(step19_kegg_gsea_enrichment() %>% as.data.frame(), 
                    rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "kegg_gsea_enrichment"), 
                                     list(extend = 'excel', filename = "kegg_gsea_enrichment"))
                    ))
    })
    
    # download the go enrichment results
    output$save_step19_kegg_gsea_results <- downloadHandler(
      filename = function() {
        paste0(input$out_step19_kegg_gsea_out_name, "-KEGG-GSEA-enrichment-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step19_kegg_gsea_enrichment(), file = file, rowNames = FALSE)
      }
    )
    
  })
  
  
  ### 19.4 draw the kegg gsea dot plot ##############################
  step19_kegg_gsea_enrichment_plot <- reactive({
    
    req(input$act_step19_draw_kegg_gsea_dot_plot)
    
    if (is.null(step19_kegg_gsea_enrichment())) {
      return(NULL)
    } else {
      kegg_gsea_results <- step19_kegg_gsea_enrichment() %>% 
        dplyr::filter(pvalue < input$in_step19_kegg_gsea_pvalue_flt,
                      qvalue < input$in_step19_kegg_gsea_qvalue_flt)
    }
    
    source("R/draw_KEGG_Terms.R")
    
    # browser()
    
    # draw the go gsea enrichment plot
    kegg_gsea_dot_plot <- draw_kegg_terms(kegg_results = kegg_gsea_results,
                                          x = "GeneRatio",
                                          color = "p.adjust",
                                          
                                          showCategory = input$in_step19_kegg_gsea_item_num,
                                          
                                          dotsize = input$in_step19_kegg_gsea_dot_size,
                                          chr_width = input$in_step19_kegg_gsea_item_width,
                                          
                                          dot_alpha = input$in_step19_kegg_gsea_dot_alpha,
                                          low_color = input$in_step19_kegg_gsea_low_color,
                                          high_color = input$in_step19_kegg_gsea_high_color,
                                          
                                          facet_group1 = input$in_step19_kegg_gsea_facet_group1,
                                          facet_group2 = NULL,
                                          
                                          gridwidth = input$in_step19_kegg_gsea_grid_width,
                                          font_size = input$in_step19_kegg_gsea_font_size,
                                          
                                          xmin = input$in_step19_kegg_gsea_xmin,
                                          xmax = input$in_step19_kegg_gsea_xmax,
                                          x_breaks = input$in_step19_kegg_gsea_break,
                                          
                                          clustered = TRUE,
                                          
                                          title = input$in_step19_kegg_gsea_title,
                                          xlab = input$in_step19_kegg_gsea_xlab,
                                          ylab = input$in_step19_kegg_gsea_ylab)
    
    return(kegg_gsea_dot_plot)
  })
  
  ## output the plot
  observeEvent(input$act_step19_draw_kegg_gsea_dot_plot, {
    
    ## output the figure
    output$out_step19_kegg_gsea_dot_plot <- renderPlot(
      width = input$out_step19_kegg_gsea_width * 100,
      height = input$out_step19_kegg_gsea_height * 100,
      {step19_kegg_gsea_enrichment_plot()})
    
    output$save_step19_kegg_gsea_dot_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step19_kegg_gsea_out_name, "-KEGG-GSEA-dotplot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step19_kegg_gsea_enrichment_plot(), filename = file,
               width = input$out_step19_kegg_gsea_width, height = input$out_step19_kegg_gsea_height)
      }
    )
  })
  
  ### 19.4 draw the kegg gsea line plot ##############################
  step19_kegg_gsea_line_plot <- reactive({
    
    req(input$act_step19_draw_kegg_gsea_lineplot)
    
    if (is.null(step19_kegg_gsea_enrichment())) {
      return(NULL)
    } else {
      kegg_gsea_results <- step19_kegg_gsea_enrichment() %>% 
        dplyr::filter(pvalue < input$in_step19_kegg_gsea_pvalue_flt,
                      qvalue < input$in_step19_kegg_gsea_qvalue_flt)
    }
    
    source("R/draw_GSEA_Terms.R")
    
    # browser()
    
    geneSetID <- gsub(" ", "", input$in_step19_kegg_gsea_line_item_name)
    
    geneSetID <- unlist(str_split(geneSetID, ",|\n"))
    
    rel_heights <- c(input$in_step19_kegg_gsea_line_height1,
                     input$in_step19_kegg_gsea_line_height2,
                     input$in_step19_kegg_gsea_line_height3)
    
    subplots <- as.numeric(input$in_step19_kegg_gsea_line_subplot)
    
    # draw the go gsea enrichment plot
    kegg_gsea_line_plot <- draw_gsea_terms(gsea_results = kegg_gsea_results,
                                           geneSetID = geneSetID,
                                           title = input$in_step19_kegg_gsea_line_title,
                                           base_size = input$in_step19_kegg_gsea_line_fontsize,
                                           rel_heights = rel_heights,
                                           subplots = subplots,
                                           color = input$in_step19_kegg_gsea_line_color,
                                           pvalue_table = input$in_step19_kegg_gsea_line_pvalue,
                                           ES_geom = input$in_step19_kegg_gsea_line_es_geom)
    
    return(kegg_gsea_line_plot)
  })
  
  ## output the plot
  observeEvent(input$act_step19_draw_kegg_gsea_lineplot, {
    
    ## output the figure
    output$out_step19_kegg_gsea_line_plot <- renderPlot(
      width = input$out_step19_kegg_gsea_line_width * 100,
      height = input$out_step19_kegg_gsea_line_height * 100,
      {step19_kegg_gsea_line_plot()})
    
    output$save_step19_kegg_gsea_lineplot <- downloadHandler(
      filename = function() {
        paste0(input$out_step19_kegg_gsea_line_out_name, "-KEGG-gseaplot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step19_kegg_gsea_line_plot(), filename = file,
               width = input$out_step19_kegg_gsea_line_width, height = input$out_step19_kegg_gsea_line_height)
      }
    )
  })
  
  
  ############################################################
  ## step20 Codon Enrich ###############################
  
  ### 20.1 import the codon count table ##############################
  step20_codon_usage_table <- reactive({
    
    req(input$act_step20_codon_import_table)
    
    if (is.null(input$in_step20_codon_table)) {return(NULL)}
    
    # browser()
    
    # import the codon count table
    codon_table <- read.xlsx(xlsxFile = input$in_step20_codon_table$datapath,
                             sheet = input$in_step20_codon_table_sheet)
    
    return(codon_table)
    
  })
  
  
  ## output the table
  observeEvent(input$act_step20_codon_import_table, {
    
    ## output the table
    output$out_step20_codon_table <- DT::renderDataTable(server = T, {
      DT::datatable(step20_codon_usage_table(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "codon_usage_table"),
                                     list(extend = 'excel', filename = "codon_usage_table"))
                    ))
    })
    
  })
  
  
  
  ### 20.2 import the gene table ##############################
  import_codon_gene_clicked <- reactiveVal(FALSE)
  
  step20_codon_gene_table <- reactive({
    
    req(input$act_step20_codon_import_gene_table)
    import_codon_gene_clicked(TRUE)
    
    if (is.null(input$in_step20_codon_gene_table)) {return(NULL)}
    
    # browser()
    
    # import the codon count table
    gene_table <- read.xlsx(xlsxFile = input$in_step20_codon_gene_table$datapath,
                             sheet = input$in_step20_codon_gene_table_sheet)
    
    return(gene_table)
    
  })
  
  
  ## output the table
  observeEvent(input$act_step20_codon_import_gene_table, {
    
    ## output the table
    output$out_step20_codon_gene_table <- DT::renderDataTable(server = T, {
      DT::datatable(step20_codon_gene_table(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "gene_class_table"),
                                     list(extend = 'excel', filename = "gene_class_table"))
                    ))
    })
    
  })
  
  
  
  ### 20.3 filter the codon table ##############################
  step20_codon_filtered_table <- reactive({
    
    req(input$act_step20_codon_filter_table)
    
    if (is.null(step20_codon_usage_table())) {return(NULL)}

    # browser()
    
    if (import_codon_gene_clicked()) {
      gene_table <- step20_codon_gene_table() %>% 
        dplyr::select(input$in_step20_codon_gene_table_column, input$in_step20_codon_gene_table_class) %>% 
        dplyr::rename("Gene" = input$in_step20_codon_gene_table_column, "Class" = input$in_step20_codon_gene_table_class)

      codon_table <- step20_codon_usage_table() %>%
        dplyr::left_join(gene_table, by = setNames("Gene", input$in_step20_codon_table_column)) %>%
        dplyr::mutate(Value = rowMeans(across(all_of(input$in_step20_codon_list)))) %>%
        dplyr::filter(Value >= input$in_step20_codon_value) %>% 
        dplyr::select(input$in_step20_codon_table_column, Class,
                      all_of(input$in_step20_codon_list), Value) %>% 
        dplyr::arrange(Class, Value) %>% 
        dplyr::mutate(Class = factor(Class, levels = unique(Class))) %>% 
        na.omit()

    } else {
      codon_table <- step20_codon_usage_table() %>%
        dplyr::mutate(Class = "Total") %>% 
        dplyr::mutate(Value = rowMeans(across(all_of(input$in_step20_codon_list)))) %>%
        dplyr::filter(Value >= input$in_step20_codon_value) %>% 
        dplyr::select(input$in_step20_codon_table_column, Class, 
                      all_of(input$in_step20_codon_list), Value) %>% 
        dplyr::arrange(Class, Value) %>% 
        na.omit()
      
    }
    
    class_count <- codon_table %>% 
      dplyr::group_by(Class) %>%
      summarise(Count = n())
    
    codon_table <- codon_table %>% 
      dplyr::left_join(class_count, by = "Class") %>% 
      dplyr::mutate(Class = paste0(Class, " (", Count, ")")) %>% 
      dplyr::select(-Count)

    return(codon_table)
  })
  
  ## output the table
  observeEvent(input$act_step20_codon_filter_table, {

    ## output the table
    output$out_step20_codon_filter_table <- DT::renderDataTable(server = T, {
      DT::datatable(step20_codon_filtered_table(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "codon_filter_table"),
                                     list(extend = 'excel', filename = "codon_filter_table"))
                    ))
    })
    
    ## save the table
    output$act_step20_codon_save_table <- downloadHandler(
      filename = function() {
        paste0(input$out_step20_codon_table_name, "-", 
               paste0(input$in_step20_codon_list, collapse = '-'), 
               "-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step20_codon_filtered_table(), file = file, rowNames = FALSE)
      }
    )
    
  })
  
  
  ### 20.4 draw the codon boxplot ##############################
  step20_codon_boxplot <- reactive({
    
    req(input$act_step20_codon_draw_box)
    
    if (is.null(step20_codon_filtered_table())) {return(NULL)}

    # browser()
    
    # draw the codon usage figure
    source("R/draw_codon_usage_box.R")
    
    cu_boxplot <- draw_cu_box(
      codon_usage = step20_codon_filtered_table(),
      
      x = "Class",
      y = "Value",
      
      box_color = input$in_step20_codon_box_color,
      box_alpha = input$in_step20_codon_box_alpha,
      box_width = input$in_step20_codon_box_width,
      
      line_width = input$in_step20_codon_box_line,
      line_color = input$in_step20_codon_box_line_color,
      
      outlier_show = input$in_step20_codon_box_outlier,
      outlier_color = input$in_step20_codon_box_outlier_color,
      outlier_size = input$in_step20_codon_box_outlier_size,
      
      font_size = input$in_step20_codon_box_fontsize,
      
      legend_position = input$in_step20_codon_box_legend,
      
      y_min = input$in_step20_codon_box_y_min,
      y_max = input$in_step20_codon_box_y_max,
      y_breaks = input$in_step20_codon_box_y_breaks,
      
      x_angle = input$in_step20_codon_box_x_angle,
      y_sqrt = input$in_step20_codon_box_y_sqrt,
      
      cu_title = input$out_step20_codon_box_title,
      cu_xlab = input$out_step20_codon_box_xlabel,
      cu_ylab = input$out_step20_codon_box_ylabel)
      
    return(cu_boxplot)
    
  })
  
  
  ## output the plot
  observeEvent(input$act_step20_codon_draw_box, {
    
    ## output the figure
    output$out_step20_codon_save_box <- renderPlot(
      width = input$out_step20_codon_box_width * 100,
      height = input$out_step20_codon_box_height * 100,
      {step20_codon_boxplot()})
    
    output$save_step20_codon_save_box <- downloadHandler(
      filename = function() {
        paste0(input$out_step20_codon_box_name,
               "-",
               paste0(input$in_step20_codon_list, collapse = '-'), 
               "-boxplot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step20_codon_boxplot(), filename = file,
               width = input$out_step20_codon_box_width, height = input$out_step20_codon_box_height)
      }
    )
  })
  
  
  
  ### 20.5 draw the codon CDF plot ##############################
  step20_codon_cdfplot <- reactive({
    
    req(input$act_step20_codon_draw_cdf)
    
    if (is.null(step20_codon_filtered_table())) {return(NULL)}
    
    # browser()
    
    # draw the codon usage figure
    source("R/draw_eCDF.R")
    
    cu_boxplot <- draw_ecdf(in_data = step20_codon_filtered_table(), 
                            x = "Value",
                            group = "Class", 
                            
                            x_min = input$in_step20_codon_cdf_x_min,
                            x_max = input$in_step20_codon_cdf_x_max,
                            x_breaks = input$in_step20_codon_cdf_x_breaks,
                            
                            title = input$in_step20_codon_cdf_title,
                            xlabel = input$in_step20_codon_cdf_xlabel,
                            ylabel = input$in_step20_codon_cdf_ylabel,
                            
                            legend_position = input$in_step20_codon_cdf_legend,
                            facet = input$in_step20_codon_cdf_facet, 
                            font_size = input$in_step20_codon_cdf_fontsize,
                            
                            log2_transform = TRUE,
                            
                            line_color = input$in_step20_codon_cdf_color,
                            line_alpha = input$in_step20_codon_cdf_alpha,
                            line_width = input$in_step20_codon_cdf_width)
    
    return(cu_boxplot)
    
  })
  
  
  ## output the plot
  observeEvent(input$act_step20_codon_draw_cdf, {
    
    ## output the figure
    output$out_step20_codon_save_cdf <- renderPlot(
      width = input$out_step20_codon_cdf_width * 100,
      height = input$out_step20_codon_cdf_height * 100,
      {step20_codon_cdfplot()})
    
    output$save_step20_codon_save_cdf <- downloadHandler(
      filename = function() {
        paste0(input$out_step20_codon_cdf_name,
               "-",
               paste0(input$in_step20_codon_list, collapse = '-'), 
               "-cdfplot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step20_codon_cdfplot(), filename = file,
               width = input$out_step20_codon_cdf_width, height = input$out_step20_codon_cdf_height)
      }
    )
  })
  
  
  ### 20.6 import the codon count table ##############################
  step20_go_codon_table <- reactive({
    
    req(input$act_step20_go_import_codon_table)
    
    if (is.null(input$in_step20_go_codon_table)) {return(NULL)}
    
    # browser()
    
    # import the codon count table
    codon_table <- read.xlsx(xlsxFile = input$in_step20_go_codon_table$datapath,
                             sheet = input$in_step20_go_codon_table_sheet)
    
    return(codon_table)
    
  })
  
  ## output the table
  observeEvent(input$act_step20_go_import_codon_table, {
    
    ## output the table
    output$out_step20_go_codon_table <- DT::renderDataTable(server = T, {
      DT::datatable(step20_go_codon_table(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "codon_count_table"),
                                     list(extend = 'excel', filename = "codon_count_table"))
                    ))
    })
    
  })
  
  
  ### 20.7 make the codon vector ##############################
  step20_go_codon_vector <- reactive({
    
    req(input$act_step20_go_create_codon_table_vector)
    
    # browser()
    
    if (is.null(step20_go_codon_table())) {return(NULL)}
    
    # make vector for gsea enrichment
    codon_vector <- step20_go_codon_table() %>% 
      dplyr::select(-Length) %>% 
      dplyr::select(input$in_step20_go_codon_gene_column, all_of(input$in_step20_go_codon_table_codon)) %>% 
      dplyr::mutate(Value = rowMeans(across(all_of(input$in_step20_go_codon_table_codon)))) %>%
      dplyr::filter(Value >= input$in_step20_go_codon_table_number) %>% 
      dplyr::mutate(Scaled = scale(Value)[, 1]) %>%
      dplyr::select(input$in_step20_go_codon_gene_column, Scaled) %>%
      na.omit() %>%
      dplyr::arrange(desc(Scaled)) %>% 
      dplyr::distinct(!!sym(input$in_step20_go_codon_gene_column), .keep_all = TRUE) %>%
      dplyr::pull(Scaled, input$in_step20_go_codon_gene_column)
    
    # show the installed packages
    select_choices <- installed.packages() %>%
      as.data.frame() %>% 
      dplyr::filter(str_detect(Package, 'org.')) 
    select_choices <- select_choices$Package
    
    updateSelectizeInput(session, "in_step20_go_codon_gsea_orgdb", choices = select_choices, selected = select_choices[1], server = TRUE)
    
    return(codon_vector)
  })
  
  ## output the table
  observeEvent(input$act_step20_go_create_codon_table_vector, {
    
    ## output the table
    output$out_step20_go_codon_table_vector <- renderPrint({
      print("The top 10 genes for GO Codon GSEA enrichment are:")
      step20_go_codon_vector() %>% head(10)
    })
    
  })
  
  
  ### 20.8 run the gsea enrichment ############################
  step20_go_codon_gsea_enrichment <- reactive({
    
    req(input$act_step20_go_codon_gsea_enrich)
    
    if (is.null(step20_go_codon_vector())) {return(NULL)}
    
    # browser()
    source("R/run_GO_GSEA_Terms.R")
    
    # run the go enrichment
    go_gsea_terms <- run_go_gsea_terms(degs = step20_go_codon_vector(),
                                       database = input$in_step20_go_codon_gsea_database,
                                       orgdb = input$in_step20_go_codon_gsea_orgdb,
                                       gson_file = input$in_step20_go_codon_gsea_gson,
                                       keytype = input$in_step20_go_codon_gsea_key_type,
                                       ontology = "ALL",
                                       eps = input$in_step20_go_codon_gsea_eps,
                                       minGSSize = input$in_step20_go_codon_gsea_minsize,
                                       maxGSSize = input$in_step20_go_codon_gsea_maxsize,
                                       pAdjustMethod = input$in_step20_go_codon_gsea_pAdjustMethod,
                                       pvalueCutoff = input$in_step20_go_codon_gsea_pvalue)
    
    # browser()
    
    return(go_gsea_terms)
  })
  
  ## output the table
  observeEvent(input$act_step20_go_codon_gsea_enrich, {
    
    ## output the table
    output$out_step20_go_codon_gsea_results <- DT::renderDataTable(server = T, {
      DT::datatable(step20_go_codon_gsea_enrichment() %>% as.data.frame(), 
                    rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "go_codon_gsea_enrichment"), 
                                     list(extend = 'excel', filename = "go_codon_gsea_enrichment"))
                    ))
    })
    
    # download the go enrichment results
    output$save_step20_go_codon_gsea_results <- downloadHandler(
      filename = function() {
        paste0(input$out_step20_go_codon_gsea_out_name, "-GO-Codon-GSEA-enrichment-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step20_go_codon_gsea_enrichment(), file = file, rowNames = FALSE)
      }
    )
    
  })
  
  
  ### 20.9 draw the go gsea dotplot ##############################
  step20_go_codon_gsea_enrichment_plot <- reactive({
    
    req(input$act_step20_draw_go_codon_gsea_dot_plot)
    
    if (is.null(step20_go_codon_gsea_enrichment())) {
      return(NULL)
    } else {
      go_codon_gsea_results <- step20_go_codon_gsea_enrichment() %>% 
        dplyr::filter(pvalue < input$in_step20_go_codon_gsea_pvalue_flt,
                      qvalue < input$in_step20_go_codon_gsea_qvalue_flt)
    }
    
    source("R/draw_GO_Terms.R")
    
    # browser()
    
    # draw the go gsea enrichment plot
    go_codon_gsea_dot_plot <- draw_go_terms(go_results = go_codon_gsea_results,
                                            x = "GeneRatio",
                                            color = "p.adjust",
                                            
                                            ontology = input$in_step20_go_codon_gsea_ontology,
                                            
                                            showCategory = input$in_step20_go_codon_gsea_item_num,
                                            
                                            dotsize = input$in_step20_go_codon_gsea_dot_size,
                                            chr_width = input$in_step20_go_codon_gsea_item_width,
                                            
                                            dot_alpha = input$in_step20_go_codon_gsea_dot_alpha,
                                            low_color = input$in_step20_go_codon_gsea_low_color,
                                            high_color = input$in_step20_go_codon_gsea_high_color,
                                            
                                            facet_group1 = input$in_step20_go_codon_gsea_facet_group1,
                                            facet_group2 = input$in_step20_go_codon_gsea_facet_group2,
                                            
                                            gridwidth = input$in_step20_go_codon_gsea_grid_width,
                                            font_size = input$in_step20_go_codon_gsea_font_size,
                                            
                                            xmin = input$in_step20_go_codon_gsea_xmin,
                                            xmax = input$in_step20_go_codon_gsea_xmax,
                                            x_breaks = input$in_step20_go_codon_gsea_break,
                                            
                                            clustered = FALSE,
                                            
                                            title = input$in_step20_go_codon_gsea_title,
                                            xlab = input$in_step20_go_codon_gsea_xlab,
                                            ylab = input$in_step20_go_codon_gsea_ylab)
    
    return(go_codon_gsea_dot_plot)
  })
  
  ## output the plot
  observeEvent(input$act_step20_draw_go_codon_gsea_dot_plot, {
    
    ## output the figure
    output$out_step20_go_codon_gsea_dot_plot <- renderPlot(
      width = input$out_step20_go_codon_gsea_width * 100,
      height = input$out_step20_go_codon_gsea_height * 100,
      {step20_go_codon_gsea_enrichment_plot()})
    
    output$save_step20_go_codon_gsea_dot_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step20_go_codon_gsea_out_name, "-GO-Codon-GSEA-dotplot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step20_go_codon_gsea_enrichment_plot(), filename = file,
               width = input$out_step20_go_codon_gsea_width, height = input$out_step20_go_codon_gsea_height)
      }
    )
  })
  
  
  ### 20.10 draw the go gseaplot ##############################
  step20_go_codon_gsea_line_plot <- reactive({
    
    req(input$act_step20_draw_go_codon_gsea_lineplot)
    
    if (is.null(step20_go_codon_gsea_enrichment())) {
      return(NULL)
    } else {
      go_codon_gsea_results <- step20_go_codon_gsea_enrichment() %>% 
        dplyr::filter(pvalue < input$in_step20_go_codon_gsea_pvalue_flt,
                      qvalue < input$in_step20_go_codon_gsea_qvalue_flt)
    }
    
    source("R/draw_GSEA_Terms.R")
    
    # browser()
    
    geneSetID <- gsub(" ", "", input$in_step20_go_codon_gsea_line_item_name)
    
    geneSetID <- unlist(str_split(geneSetID, ",|\n"))
    
    rel_heights <- c(input$in_step20_go_codon_gsea_line_height1,
                     input$in_step20_go_codon_gsea_line_height2,
                     input$in_step20_go_codon_gsea_line_height3)
    
    subplots <- as.numeric(input$in_step20_go_codon_gsea_line_subplot)
    
    # draw the go gsea enrichment plot
    go_codon_gsea_line_plot <- draw_gsea_terms(gsea_results = go_codon_gsea_results,
                                               geneSetID = geneSetID,
                                               title = input$in_step20_go_codon_gsea_line_title,
                                               base_size = input$in_step20_go_codon_gsea_line_fontsize,
                                               rel_heights = rel_heights,
                                               subplots = subplots,
                                               color = input$in_step20_go_codon_gsea_line_color,
                                               pvalue_table = input$in_step20_go_codon_gsea_line_pvalue,
                                               ES_geom = input$in_step20_go_codon_gsea_line_es_geom)
    
    return(go_codon_gsea_line_plot)
  })
  
  ## output the plot
  observeEvent(input$act_step20_draw_go_codon_gsea_lineplot, {
    
    ## output the figure
    output$out_step20_go_codon_gsea_line_plot <- renderPlot(
      width = input$out_step20_go_codon_gsea_line_width * 100,
      height = input$out_step20_go_codon_gsea_line_height * 100,
      {step20_go_codon_gsea_line_plot()})
    
    output$save_step20_go_codon_gsea_lineplot <- downloadHandler(
      filename = function() {
        paste0(input$out_step20_go_codon_gsea_line_out_name, "-GO-Codon-gseaplot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step20_go_codon_gsea_line_plot(), filename = file,
               width = input$out_step20_go_codon_gsea_line_width, height = input$out_step20_go_codon_gsea_line_height)
      }
    )
  })
  
  
  
  
  
  ### 20.11 import the codon count table ##############################
  step20_kegg_codon_table <- reactive({
    
    req(input$act_step20_kegg_import_codon_table)
    
    if (is.null(input$in_step20_kegg_codon_table)) {return(NULL)}
    
    # browser()
    
    # import the codon count table
    codon_table <- read.xlsx(xlsxFile = input$in_step20_kegg_codon_table$datapath,
                             sheet = input$in_step20_kegg_codon_table_sheet)
    
    return(codon_table)
    
  })
  
  ## output the table
  observeEvent(input$act_step20_kegg_import_codon_table, {
    
    ## output the table
    output$out_step20_kegg_codon_table <- DT::renderDataTable(server = T, {
      DT::datatable(step20_kegg_codon_table(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "codon_count_table"),
                                     list(extend = 'excel', filename = "codon_count_table"))
                    ))
    })
    
  })
  
  
  ### 20.12 make the codon vector ##############################
  step20_kegg_codon_vector <- reactive({
    
    req(input$act_step20_kegg_create_codon_table_vector)
    
    # browser()
    
    if (is.null(step20_kegg_codon_table())) {return(NULL)}
    
    # make vector for gsea enrichment
    codon_vector <- step20_kegg_codon_table() %>% 
      dplyr::select(-Length) %>% 
      dplyr::select(input$in_step20_kegg_codon_gene_column, all_of(input$in_step20_kegg_codon_table_codon)) %>% 
      dplyr::mutate(Value = rowMeans(across(all_of(input$in_step20_kegg_codon_table_codon)))) %>%
      dplyr::filter(Value >= input$in_step20_kegg_codon_table_number) %>% 
      dplyr::mutate(Scaled = scale(Value)[, 1]) %>%
      dplyr::select(input$in_step20_kegg_codon_gene_column, Scaled) %>%
      dplyr::arrange(desc(Scaled)) %>% 
      dplyr::distinct(!!sym(input$in_step20_kegg_codon_gene_column), .keep_all = TRUE) %>%
      dplyr::pull(Scaled, input$in_step20_kegg_codon_gene_column)
    
    return(codon_vector)
  })
  
  ## output the table
  observeEvent(input$act_step20_kegg_create_codon_table_vector, {
    
    ## output the table
    output$out_step20_kegg_codon_table_vector <- renderPrint({
      print("The top 10 genes for KEGG GSEA enrichment are:")
      step20_kegg_codon_vector() %>% head(10)
    })
    
  })
  
  
  ### 20.13 run the gsea enrichment ############################
  step20_kegg_codon_gsea_enrichment <- reactive({
    
    req(input$act_step20_kegg_codon_gsea_enrich)
    
    # browser()
    if (is.null(step20_kegg_codon_vector())) {return(NULL)}
    
    source("R/run_kEGG_GSEA_Terms.R")
    
    # run the kegg enrichment
    kegg_codon_gsea_terms <- run_kegg_gsea_terms(degs = step20_kegg_codon_vector(),
                                                 database = input$in_step20_kegg_codon_gsea_database,
                                                 gson_file = input$in_step20_kegg_codon_gsea_gson,
                                                 species = input$in_step20_kegg_codon_gsea_species,
                                                 kegg_type = input$in_step20_kegg_codon_gsea_type,
                                                 keytype = input$in_step20_kegg_codon_gsea_key_type,
                                                 eps = input$in_step20_kegg_codon_gsea_eps,
                                                 minGSSize = input$in_step20_kegg_codon_gsea_minsize,
                                                 maxGSSize = input$in_step20_kegg_codon_gsea_maxsize,
                                                 pAdjustMethod = input$in_step20_kegg_codon_gsea_pAdjustMethod,
                                                 pvalueCutoff = input$in_step20_kegg_codon_gsea_pvalue)
    
    return(kegg_codon_gsea_terms)
  })
  
  ## output the table
  observeEvent(input$act_step20_kegg_codon_gsea_enrich, {
    
    ## output the table
    output$out_step20_kegg_codon_gsea_results <- DT::renderDataTable(server = T, {
      DT::datatable(step20_kegg_codon_gsea_enrichment() %>% as.data.frame(), 
                    rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(6, 10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "kegg_codon_gsea_enrichment"), 
                                     list(extend = 'excel', filename = "kegg_codon_gsea_enrichment"))
                    ))
    })
    
    # download the go enrichment results
    output$save_step20_kegg_codon_gsea_results <- downloadHandler(
      filename = function() {
        paste0(input$out_step20_kegg_codon_gsea_out_name, "-KEGG-Codon-GSEA-enrichment-", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        write.xlsx(x = step20_kegg_codon_gsea_enrichment(), file = file, rowNames = FALSE)
      }
    )
    
  })
  
  ### 20.14 draw the kegg gsea enrichment ##############################
  step20_kegg_codon_gsea_enrichment_plot <- reactive({
    
    req(input$act_step20_draw_kegg_codon_gsea_dot_plot)
    
    if (is.null(step20_kegg_codon_gsea_enrichment())) {
      return(NULL)
    } else {
      kegg_codon_gsea_results <- step20_kegg_codon_gsea_enrichment() %>% 
        dplyr::filter(pvalue < input$in_step20_kegg_codon_gsea_pvalue_flt,
                      qvalue < input$in_step20_kegg_codon_gsea_qvalue_flt)
    }
    
    source("R/draw_KEGG_Terms.R")
    
    # browser()
    
    # draw the go gsea enrichment plot
    kegg_codon_gsea_dot_plot <- draw_kegg_terms(kegg_results = kegg_codon_gsea_results,
                                                x = "GeneRatio",
                                                color = "p.adjust",
                                                
                                                showCategory = input$in_step20_kegg_codon_gsea_item_num,
                                                
                                                dotsize = input$in_step20_kegg_codon_gsea_dot_size,
                                                chr_width = input$in_step20_kegg_codon_gsea_item_width,
                                                
                                                dot_alpha = input$in_step20_kegg_codon_gsea_dot_alpha,
                                                low_color = input$in_step20_kegg_codon_gsea_low_color,
                                                high_color = input$in_step20_kegg_codon_gsea_high_color,
                                                
                                                facet_group1 = input$in_step20_kegg_codon_gsea_facet_group1,
                                                facet_group2 = NULL,
                                                
                                                gridwidth = input$in_step20_kegg_codon_gsea_grid_width,
                                                font_size = input$in_step20_kegg_codon_gsea_font_size,
                                                
                                                xmin = input$in_step20_kegg_codon_gsea_xmin,
                                                xmax = input$in_step20_kegg_codon_gsea_xmax,
                                                x_breaks = input$in_step20_kegg_codon_gsea_break,
                                                
                                                clustered = TRUE,
                                                
                                                title = input$in_step20_kegg_codon_gsea_title,
                                                xlab = input$in_step20_kegg_codon_gsea_xlab,
                                                ylab = input$in_step20_kegg_codon_gsea_ylab)
    
    return(kegg_codon_gsea_dot_plot)
  })
  
  ## output the plot
  observeEvent(input$act_step20_draw_kegg_codon_gsea_dot_plot, {
    
    ## output the figure
    output$out_step20_kegg_codon_gsea_dot_plot <- renderPlot(
      width = input$out_step20_kegg_codon_gsea_width * 100,
      height = input$out_step20_kegg_codon_gsea_height * 100,
      {step20_kegg_codon_gsea_enrichment_plot()})
    
    output$save_step20_kegg_codon_gsea_dot_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step20_kegg_codon_gsea_out_name, "-KEGG-Codon-GSEA-dotplot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step20_kegg_codon_gsea_enrichment_plot(), filename = file,
               width = input$out_step20_kegg_codon_gsea_width, height = input$out_step20_kegg_codon_gsea_height)
      }
    )
  })
  
  ### 20.15 draw the kegg gsea enrichment ##############################
  step20_kegg_codon_gsea_line_plot <- reactive({
    
    req(input$act_step20_draw_kegg_codon_gsea_lineplot)
    
    if (is.null(step20_kegg_codon_gsea_enrichment())) {
      return(NULL)
    } else {
      kegg_codon_gsea_results <- step20_kegg_codon_gsea_enrichment() %>% 
        dplyr::filter(pvalue < input$in_step20_kegg_codon_gsea_pvalue_flt,
                      qvalue < input$in_step20_kegg_codon_gsea_qvalue_flt)
    }
    
    source("R/draw_GSEA_Terms.R")
    
    # browser()
    
    geneSetID <- gsub(" ", "", input$in_step20_kegg_codon_gsea_line_item_name)
    
    geneSetID <- unlist(str_split(geneSetID, ",|\n"))
    
    rel_heights <- c(input$in_step20_kegg_codon_gsea_line_height1,
                     input$in_step20_kegg_codon_gsea_line_height2,
                     input$in_step20_kegg_codon_gsea_line_height3)
    
    subplots <- as.numeric(input$in_step20_kegg_codon_gsea_line_subplot)
    
    # draw the go gsea enrichment plot
    kegg_codon_gsea_line_plot <- draw_gsea_terms(gsea_results = kegg_codon_gsea_results,
                                                 geneSetID = geneSetID,
                                                 title = input$in_step20_kegg_codon_gsea_line_title,
                                                 base_size = input$in_step20_kegg_codon_gsea_line_fontsize,
                                                 rel_heights = rel_heights,
                                                 subplots = subplots,
                                                 color = input$in_step20_kegg_codon_gsea_line_color,
                                                 pvalue_table = input$in_step20_kegg_codon_gsea_line_pvalue,
                                                 ES_geom = input$in_step20_kegg_codon_gsea_line_es_geom)
    
    return(kegg_codon_gsea_line_plot)
  })
  
  ## output the plot
  observeEvent(input$act_step20_draw_kegg_codon_gsea_lineplot, {
    
    ## output the figure
    output$out_step20_kegg_codon_gsea_line_plot <- renderPlot(
      width = input$out_step20_kegg_codon_gsea_line_width * 100,
      height = input$out_step20_kegg_codon_gsea_line_height * 100,
      {step20_kegg_codon_gsea_line_plot()})
    
    output$save_step20_kegg_codon_gsea_lineplot <- downloadHandler(
      filename = function() {
        paste0(input$out_step20_kegg_codon_gsea_line_out_name, "-KEGG-gseaplot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step20_kegg_codon_gsea_line_plot(), filename = file,
               width = input$out_step20_kegg_codon_gsea_line_width, height = input$out_step20_kegg_codon_gsea_line_height)
      }
    )
  })
  
  
  ############################################################
  ## step21 Gene plot ###############################
  
  ### 21.0 set the design ##############################
  observeEvent(input$act_step1_import_design, {
    col_name <- names(step1_design())
    updateSelectizeInput(session, "in_step21_isoforms_rna_column", choices = col_name, selected = col_name[3], server = TRUE)
    updateSelectizeInput(session, "in_step21_isoforms_ribo_column", choices = col_name, selected = col_name[3], server = TRUE)
    
    # browser()
    observe({
      req(input$in_step21_isoforms_rna_column)
      
      rna_design <- step1_design() %>% dplyr::filter(SeqType == "RNA")
      rna_select_choices <- unique(rna_design[[input$in_step21_isoforms_rna_column]])
      updateSelectizeInput(session, "in_step21_isoforms_rna_groups", choices = rna_select_choices, selected = rna_select_choices, server = TRUE)
    })
    
    observe({
      req(input$in_step21_isoforms_ribo_column)
      
      ribo_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      ribo_select_choices <- unique(ribo_design[[input$in_step21_isoforms_ribo_column]])
      updateSelectizeInput(session, "in_step21_isoforms_ribo_groups", choices = ribo_select_choices, selected = ribo_select_choices, server = TRUE)
    })
    
  })
  
  ### 21.1 import the rna-seq isoforms data ##############################
  step21_isoforms_rna_reads <- reactive({
    
    req(input$act_step21_isoforms_import_rna_seq)
    
    if (is.null(input$in_step21_isoforms_rna_seq)) {return(NULL)}
    
    # import the raw data table
    isoforms_rna_reads <- fread(file = input$in_step21_isoforms_rna_seq$datapath,
                                index = NULL, header="auto", check.names = FALSE, sep = '\t')
    
    # update the input$in_step21_isoforms_gene with unique gene names
    # if (nchar(input$in_step21_isoforms_gene) == 0 & (is.null(step21_isoforms_ribo_reads) | is.null(step21_isoforms_rna_reads))) {
    #   uniq_gene_names <- unique(isoforms_rna_reads$name)
    #   updateSelectizeInput(session, "in_step21_isoforms_gene", choices = uniq_gene_names, selected = uniq_gene_names[1], server = TRUE)
    # }
    
    uniq_gene_names <- unique(isoforms_rna_reads$name)
    updateSelectizeInput(session, "in_step21_isoforms_gene", choices = uniq_gene_names, selected = uniq_gene_names[1], server = TRUE)
    
    # browser()
    # calculate the average of the rna-seq isoforms data
    if (!is.null(input$in_step21_isoforms_rna_groups)) {
      
      if (is.null(step1_design())) {
        return(isoforms_rna_reads)
      }
      
      rna_design <- step1_design() %>% dplyr::filter(SeqType == "RNA")
      
      for (i in 1:length(input$in_step21_isoforms_rna_groups)) {
        now_group <- input$in_step21_isoforms_rna_groups[i]
        now_name <- rna_design %>%
          dplyr::filter(!!sym(input$in_step21_isoforms_rna_column) == now_group) %>% 
          dplyr::pull("Sample")
        
        isoforms_rna_reads <- isoforms_rna_reads %>%
          dplyr::mutate(!!now_group := rowMeans(dplyr::select(., all_of(now_name))))
      }
      
      isoforms_rna_reads <- isoforms_rna_reads %>%
        dplyr::select(any_of(c("name", "now_nt", "from_tis", "from_tts", "region", "codon", "frame", 
                               input$in_step21_isoforms_rna_groups))) %>%
        tidyr::pivot_longer(cols = input$in_step21_isoforms_rna_groups, 
                            names_to = "Sample", values_to = "Abundance")
      
      return(isoforms_rna_reads)
      
    } else {
      
      isoforms_rna_reads <- isoforms_rna_reads %>%
        tidyr::pivot_longer(cols = -c("name", "now_nt", "from_tis", "from_tts", "region", "codon", "frame"), 
                            names_to = "Sample", values_to = "Abundance")
      
      return(isoforms_rna_reads)
    }
  })
  
  ## output the table
  observeEvent(input$act_step21_isoforms_import_rna_seq, {
    
    ## output the table
    output$out_step21_isoforms_rna_seq <- DT::renderDataTable(server = T, {
      DT::datatable(step21_isoforms_rna_reads(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "isoforms_rna_seq"),
                                     list(extend = 'excel', filename = "isoforms_rna_seq"))
                    ))
    })
    
  })
  
  
  ### 21.2 import the ribo-seq isoforms data ##############################
  step21_isoforms_ribo_reads <- reactive({
    
    req(input$act_step21_isoforms_import_ribo_seq)
    
    if (is.null(input$in_step21_isoforms_ribo_seq)) {return(NULL)}
    
    # import the raw data table
    isoforms_ribo_reads <- fread(file = input$in_step21_isoforms_ribo_seq$datapath,
                                 index = NULL, header="auto", check.names = FALSE, sep = '\t')
    
    # update the input$in_step21_isoforms_gene with unique gene names
    # if (nchar(input$in_step21_isoforms_gene) == 0 & (is.null(step21_isoforms_ribo_reads) | is.null(step21_isoforms_rna_reads))) {
    #   uniq_gene_names <- unique(isoforms_ribo_reads$name)
    #   updateSelectizeInput(session, "in_step21_isoforms_gene", choices = uniq_gene_names, selected = uniq_gene_names[1], server = TRUE)
    # }
    
    uniq_gene_names <- unique(isoforms_ribo_reads$name)
    updateSelectizeInput(session, "in_step21_isoforms_gene", choices = uniq_gene_names, selected = uniq_gene_names[1], server = TRUE)
    
    # browser()
    # calculate the average of the ribo-seq isoforms data
    if (!is.null(input$in_step21_isoforms_ribo_groups)) {
      
      if (is.null(step1_design())) {
        return(isoforms_ribo_reads)
      }
      
      ribo_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      
      for (i in 1:length(input$in_step21_isoforms_ribo_groups)) {
        now_group <- input$in_step21_isoforms_ribo_groups[i]
        now_name <- ribo_design %>%
          dplyr::filter(!!sym(input$in_step21_isoforms_ribo_column) == now_group) %>% 
          dplyr::pull("Sample")
        
        isoforms_ribo_reads <- isoforms_ribo_reads %>%
          dplyr::mutate(!!now_group := rowMeans(dplyr::select(., all_of(now_name))))
      }
      
      isoforms_ribo_reads <- isoforms_ribo_reads %>%
        dplyr::select(any_of(c("name", "now_nt", "from_tis", "from_tts", "region", "codon", "frame", 
                               input$in_step21_isoforms_ribo_groups))) %>%
        tidyr::pivot_longer(cols = input$in_step21_isoforms_ribo_groups, 
                            names_to = "Sample", values_to = "Abundance")
      
      return(isoforms_ribo_reads)
      
    } else {
      
      isoforms_ribo_reads <- isoforms_ribo_reads %>%
        tidyr::pivot_longer(cols = -c("name", "now_nt", "from_tis", "from_tts", "region", "codon", "frame"), 
                            names_to = "Sample", values_to = "Abundance")
      
      return(isoforms_ribo_reads)
    }
  })
  
  ## output the table
  observeEvent(input$act_step21_isoforms_import_ribo_seq, {
    
    ## output the table
    output$out_step21_isoforms_ribo_seq <- DT::renderDataTable(server = T, {
      DT::datatable(step21_isoforms_ribo_reads(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE,
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "isoforms_ribo_seq"),
                                     list(extend = 'excel', filename = "isoforms_ribo_seq"))
                    ))
    })
    
  })
  
  
  ### 21.3 draw the rna-seq isoforms data ##############################
  step21_isoforms_rna_plot <- reactive({
    
    req(input$act_step21_isoforms_draw_rna_seq)
    
    if (is.null(step21_isoforms_rna_reads())) {return(NULL)}
    if (is.null(input$in_step21_isoforms_gene)) {return(NULL)}
    
    # filter the raw data table
    isoforms_rna_reads <- step21_isoforms_rna_reads() %>%
      dplyr::filter(name == input$in_step21_isoforms_gene)
    
    # draw the gene schematic
    isoforms_rna_uniq <- isoforms_rna_reads %>% 
      dplyr::filter(Sample == Sample[1])
    
    utr5_length <- sum(str_count(isoforms_rna_uniq$region, '5utr'))
    cds_length <- sum(str_count(isoforms_rna_uniq$region, 'cds'))
    utr3_length <- sum(str_count(isoforms_rna_uniq$region, '3utr'))
    gene_length <- max(isoforms_rna_uniq$now_nt)
    
    xmin <- ifelse(!is.na(input$in_step21_isoforms_rna_x_min), input$in_step21_isoforms_rna_x_min, NA)
    xmax <- ifelse(!is.na(input$in_step21_isoforms_rna_x_max), input$in_step21_isoforms_rna_x_max, NA)
    
    source("R/draw_Isoforms_Schematic.R")
    isoforms_rna_schematic <- draw_isoforms_schematic(gene_name = input$in_step21_isoforms_gene,
                                                      
                                                      utr5_length = utr5_length,
                                                      cds_length = cds_length,
                                                      utr3_length = utr3_length,
                                                      
                                                      gene_color = c("#75aadb", "#264abd", "#75aadb"),
                                                      gene_width = c(2, 5, 2),
                                                      fill_alpha = input$in_step21_isoforms_alpha,
                                                      
                                                      xmin = xmin,
                                                      xmax = xmax,
                                                      xbreaks = input$in_step21_isoforms_rna_x_break,
                                                      
                                                      arrow_length = 0.1,
                                                      arrow_width = 0.8,
                                                      arrow_color = "white",
                                                      
                                                      font_size = input$in_step21_isoforms_font_size)
    
    # browser()
    # draw the rna-seq isoforms data
    # c("name", "now_nt", "from_tis", "from_tts", "region", "codon", "frame", "Sample", "Abundance"), 
    source("R/draw_Isoforms_Mapped.R")
    
    ymin <- ifelse(!is.na(input$in_step21_isoforms_rna_y_min), input$in_step21_isoforms_rna_y_min, NA)
    ymax <- ifelse(!is.na(input$in_step21_isoforms_rna_y_max), input$in_step21_isoforms_rna_y_max, NA)
    
    if (input$in_step21_isoforms_y_scale == 'log') {
      isoforms_rna_reads <- isoforms_rna_reads %>%
        dplyr::mutate(Abundance = log2(Abundance + 1))
    }
    
    isoforms_rna_reads_plot <- draw_isoforms_mapped(gene_name = input$in_step21_isoforms_gene,
                                                    gene_reads = isoforms_rna_reads,
                                                    
                                                    x = 'now_nt',
                                                    y = 'Abundance',
                                                    color = 'Sample',
                                                    
                                                    xlabel = input$in_step21_isoforms_xlabel,
                                                    ylabel = input$in_step21_isoforms_ylabel,
                                                    title = paste0(input$in_step21_isoforms_gene, " (RNA-seq)"),
                                                    
                                                    gene_color = input$in_step21_isoforms_color,
                                                    
                                                    plot_type = input$in_step21_isoforms_type,
                                                    
                                                    line_width = input$in_step21_isoforms_line_width,
                                                    grid_width = input$in_step21_isoforms_grid_width,
                                                    fill_alpha = input$in_step21_isoforms_alpha,
                                                    
                                                    xmin = xmin,
                                                    xmax = xmax,
                                                    xbreaks = input$in_step21_isoforms_rna_x_break,
                                                    
                                                    ymin = ymin,
                                                    ymax = ymax,
                                                    ybreaks = input$in_step21_isoforms_rna_y_break,
                                                    
                                                    sqrt = input$in_step21_isoforms_y_scale,
                                                    
                                                    facet = input$in_step21_isoforms_facet,
                                                    show_legend = input$in_step21_isoforms_legend,
                                                    
                                                    font_size = input$in_step21_isoforms_font_size)
    
    # combine the gene schematic and the rna-seq data
    isoforms_rna_plot <- cowplot::plot_grid(isoforms_rna_reads_plot, isoforms_rna_schematic,
                                            align = 'hv', nrow = 2, axis = 'tblr',
                                            rel_heights = c(input$in_step21_isoforms_top_height, 
                                                            input$in_step21_isoforms_bottom_height))
    
    return(isoforms_rna_plot)
    
  })
  
  ## output the figure
  observeEvent(input$act_step21_isoforms_draw_rna_seq, {
    
    ## output the figure
    output$out_step21_isoforms_rna_seq_plot <- renderPlot(
      width = input$out_step21_isoforms_fig_width * 100,
      height = input$out_step21_isoforms_fig_height * 100,
      {step21_isoforms_rna_plot()})
    
    output$save_step21_isoforms_draw_rna_seq <- downloadHandler(
      filename = function() {
        paste0(input$out_step21_isoforms_out_name, "-", input$in_step21_isoforms_gene, "-RNA-seq-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step21_isoforms_rna_plot(), filename = file,
               width = input$out_step21_isoforms_fig_width, height = input$out_step21_isoforms_fig_height)
      }
    )
    
  })
  
  
  ### 21.4 draw the ribo-seq isoforms data ##############################
  step21_isoforms_ribo_plot <- reactive({
    
    req(input$act_step21_isoforms_draw_ribo_seq)
    
    if (is.null(step21_isoforms_ribo_reads())) {return(NULL)}
    if (is.null(input$in_step21_isoforms_gene)) {return(NULL)}
    
    # filter the raw data table
    isoforms_ribo_reads <- step21_isoforms_ribo_reads() %>%
      dplyr::filter(name == input$in_step21_isoforms_gene)
    
    # draw the gene schematic
    isoforms_ribo_uniq <- isoforms_ribo_reads %>% 
      dplyr::filter(Sample == Sample[1])
    
    utr5_length <- sum(str_count(isoforms_ribo_uniq$region, '5utr'))
    cds_length <- sum(str_count(isoforms_ribo_uniq$region, 'cds'))
    utr3_length <- sum(str_count(isoforms_ribo_uniq$region, '3utr'))
    gene_length <- max(isoforms_ribo_uniq$now_nt)
    
    xmin <- ifelse(!is.na(input$in_step21_isoforms_ribo_x_min), input$in_step21_isoforms_ribo_x_min, NA)
    xmax <- ifelse(!is.na(input$in_step21_isoforms_ribo_x_max), input$in_step21_isoforms_ribo_x_max, NA)
    
    source("R/draw_Isoforms_Schematic.R")
    isoforms_ribo_schematic <- draw_isoforms_schematic(gene_name = input$in_step21_isoforms_gene,
                                                       
                                                       utr5_length = utr5_length,
                                                       cds_length = cds_length,
                                                       utr3_length = utr3_length,
                                                       
                                                       gene_color = c("#75aadb", "#264abd", "#75aadb"),
                                                       gene_width = c(2, 5, 2),
                                                       fill_alpha = input$in_step21_isoforms_alpha,
                                                       
                                                       xmin = xmin,
                                                       xmax = xmax,
                                                       xbreaks = input$in_step21_isoforms_ribo_x_break,
                                                       
                                                       arrow_length = 0.1,
                                                       arrow_width = 0.8,
                                                       arrow_color = "white",
                                                       
                                                       font_size = input$in_step21_isoforms_font_size)
    
    # browser()
    # draw the ribo-seq isoforms data
    # c("name", "now_nt", "from_tis", "from_tts", "region", "codon", "frame", "Sample", "Abundance"), 
    source("R/draw_Isoforms_Mapped.R")
    
    ymin <- ifelse(!is.na(input$in_step21_isoforms_ribo_y_min), input$in_step21_isoforms_ribo_y_min, NA)
    ymax <- ifelse(!is.na(input$in_step21_isoforms_ribo_y_max), input$in_step21_isoforms_ribo_y_max, NA)
    
    if (input$in_step21_isoforms_y_scale == 'log') {
      isoforms_ribo_reads <- isoforms_ribo_reads %>%
        dplyr::mutate(Abundance = log2(Abundance + 1))
    }
    
    isoforms_ribo_reads_plot <- draw_isoforms_mapped(gene_name = input$in_step21_isoforms_gene,
                                                     gene_reads = isoforms_ribo_reads,
                                                     
                                                     x = 'now_nt',
                                                     y = 'Abundance',
                                                     color = 'Sample',
                                                     
                                                     xlabel = input$in_step21_isoforms_xlabel,
                                                     ylabel = input$in_step21_isoforms_ylabel,
                                                     title = paste0(input$in_step21_isoforms_gene, " (Ribo-seq)"),
                                                     
                                                     gene_color = input$in_step21_isoforms_color,
                                                     
                                                     plot_type = input$in_step21_isoforms_type,
                                                     
                                                     line_width = input$in_step21_isoforms_line_width,
                                                     grid_width = input$in_step21_isoforms_grid_width,
                                                     fill_alpha = input$in_step21_isoforms_alpha,
                                                     
                                                     xmin = xmin,
                                                     xmax = xmax,
                                                     xbreaks = input$in_step21_isoforms_ribo_x_break,
                                                     
                                                     ymin = ymin,
                                                     ymax = ymax,
                                                     ybreaks = input$in_step21_isoforms_ribo_y_break,
                                                     
                                                     sqrt = input$in_step21_isoforms_y_scale,
                                                     
                                                     facet = input$in_step21_isoforms_facet,
                                                     show_legend = input$in_step21_isoforms_legend,
                                                     
                                                     font_size = input$in_step21_isoforms_font_size)
    
    # combine the gene schematic and the rna-seq data
    isoforms_ribo_plot <- cowplot::plot_grid(isoforms_ribo_reads_plot, isoforms_ribo_schematic,
                                             align = 'hv', nrow = 2, axis = 'tblr',
                                             rel_heights = c(input$in_step21_isoforms_top_height, 
                                                             input$in_step21_isoforms_bottom_height))
    
    return(isoforms_ribo_plot)
    
  })
  
  ## output the figure
  observeEvent(input$act_step21_isoforms_draw_ribo_seq, {
    
    ## output the figure
    output$out_step21_isoforms_ribo_seq_plot <- renderPlot(
      width = input$out_step21_isoforms_fig_width * 100,
      height = input$out_step21_isoforms_fig_height * 100,
      {step21_isoforms_ribo_plot()})
    
    output$save_step21_isoforms_draw_ribo_seq <- downloadHandler(
      filename = function() {
        paste0(input$out_step21_isoforms_out_name, "-", input$in_step21_isoforms_gene, "-Ribo-seq-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step21_isoforms_ribo_plot(), filename = file,
               width = input$out_step21_isoforms_fig_width, height = input$out_step21_isoforms_fig_height)
      }
    )
    
  })
  
  
  ############################################################
  ## step22 Pausing ###############################
  
  ### 22.0 set the design ##############################
  observeEvent(input$act_step1_import_design, {
    col_name <- names(step1_design())
    updateSelectizeInput(session, "in_step22_column", choices = col_name, selected = col_name[3], server = TRUE)
    updateSelectizeInput(session, "in_step22_diff_column", choices = col_name, selected = col_name[3], server = TRUE)
    
    # browser()
    observe({
      req(input$in_step22_column)
      step1_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      select_choices <- step1_design[[input$in_step22_column]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step22_groups", choices = select_choices, selected = select_choices, server = TRUE)
    })
    
    observe({
      req(input$in_step22_diff_column)
      step1_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      select_choices <- step1_design[[input$in_step22_diff_column]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step22_diff_groups1", choices = select_choices, selected = select_choices[1], server = TRUE)
      updateSelectizeInput(session, "in_step22_diff_groups2", choices = select_choices, selected = select_choices[2], server = TRUE)
    })
    
  })
  
  ### 22.1 import the pausing score data ##############################
  step22_pausing_score <- reactive({
    
    req(input$act_step22_import_pausing)
    
    if (is.null(input$in_step22_pausing)) {return(NULL)}
    
    # import the raw data table
    pausing_score <- read.table(file = input$in_step22_pausing$datapath, 
                                row.names = NULL, header = T, check.names = F, sep = '\t')
    
    # browser()
    absolute_total_ps <- pausing_score %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_absolute_total_ps")) %>% 
      dplyr::rename_with(~str_remove(., "_absolute_total_ps"), contains("_absolute_total_ps")) %>% 
      dplyr::mutate(Class = "Total")
    absolute_valid_ps <- pausing_score %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_absolute_valid_ps")) %>% 
      dplyr::rename_with(~str_remove(., "_absolute_valid_ps"), contains("_absolute_valid_ps")) %>% 
      dplyr::mutate(Class = "Valid")
    absolute_ps <- rbind(absolute_valid_ps, absolute_total_ps)
    
    relative_total_pausing <- pausing_score %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_relative_total_ps")) %>% 
      dplyr::rename_with(~str_remove(., "_relative_total_ps"), contains("_relative_total_ps")) %>%
      dplyr::mutate(Class = "Total")
    relative_valid_pausing <- pausing_score %>%
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_relative_valid_ps")) %>% 
      dplyr::rename_with(~str_remove(., "_relative_valid_ps"), contains("_relative_valid_ps")) %>%
      dplyr::mutate(Class = "Valid")
    relative_ps <- rbind(relative_valid_pausing, relative_total_pausing)
    
    # calculate the average of the pausing score
    if (!is.null(input$in_step22_groups)) {
      
      if (is.null(step1_design())) {
        return(list(absolute_pausing = absolute_ps,
                    relative_pausing = relative_ps))
      }
      
      # browser()
      
      absolute_pausing <- absolute_ps %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr', 'Class'), names_to = "Sample", values_to = "PauseScore") %>%
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step22_column) %in% input$in_step22_groups) %>%
        dplyr::group_by(Codon, AA, Abbr, Class, !!sym(input$in_step22_column)) %>% 
        dplyr::reframe(PauseScore = mean(PauseScore)) %>% 
        dplyr::ungroup()
      
      relative_pausing <- relative_ps %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr', 'Class'), names_to = "Sample", values_to = "PauseScore") %>%
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step22_column) %in% input$in_step22_groups) %>%
        dplyr::group_by(Codon, AA, Abbr, Class, !!sym(input$in_step22_column)) %>% 
        dplyr::reframe(PauseScore = mean(PauseScore)) %>% 
        dplyr::ungroup()
      
      return(list(absolute_pausing = absolute_pausing,
                  relative_pausing = relative_pausing))
    } else {
      
      absolute_pausing <- absolute_ps %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr', 'Class'), names_to = "Sample", values_to = "PauseScore") 
      
      relative_pausing <- relative_ps %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr', 'Class'), names_to = "Sample", values_to = "PauseScore")
      
      return(list(absolute_pausing = absolute_pausing,
                  relative_pausing = relative_pausing))
    }
    
  })
  
  ## output the table
  observeEvent(input$act_step22_import_pausing, {
    
    ## output the table
    output$out_step22_absolute_pausing <- DT::renderDataTable(server = T, {
      DT::datatable(step22_pausing_score()$absolute_pausing, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "absolute_ps_filter"), 
                                     list(extend = 'excel', filename = "absolute_ps_filter"))
                    ))
    })
    
    output$out_step22_relative_pausing <- DT::renderDataTable(server = T, {
      DT::datatable(step22_pausing_score()$relative_pausing, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "relative_ps_filter"), 
                                     list(extend = 'excel', filename = "relative_ps_filter"))
                    ))
    })
    
  })
  
  
  ### 22.2 draw the pausing score plot ##############################
  step22_pausing_plot <- reactive({
    
    req(input$act_step22_draw_pausing_score)
    
    if (is.null(step22_pausing_score())) {return(NULL)}
    
    # browser()
    source("R/draw_codon.R")
    
    absolute_codon_table <- step22_pausing_score()$absolute_pausing %>% 
      dplyr::filter(Class == input$in_step22_valid)
    
    relative_codon_table <- step22_pausing_score()$relative_pausing %>% 
      dplyr::filter(Class == input$in_step22_valid)
    
    # set the dot color aes
    if (nchar(input$in_step22_column) == 0) {
      fill_color <- "Sample"
    } else {
      fill_color <- input$in_step22_column
    }
    
    # absolute pausing score
    absolute_pausing_plot <- draw_codon(codon = absolute_codon_table,
                                        x = "Codon",
                                        y = "PauseScore",
                                        xlab = "Codon",
                                        ylab = "Pause Score",
                                        title = "Codon pausing score",
                                        fill_color = fill_color,
                                        color_palette = input$in_step22_dot_color,
                                        color_alpha = input$in_step22_alpha,
                                        dot_size = input$in_step22_dot_size,
                                        line_size = input$in_step22_line_size,
                                        
                                        grid_width = input$in_step22_grid_width,
                                        grid_color = input$in_step22_grid_color,
                                        
                                        font_size = input$in_step22_font_size,
                                        font_color = 'black',
                                        legend_row = input$in_step22_legend_row,
                                        
                                        wrap = input$in_step22_wrap,
                                        wrap_group = "AA",
                                        wrap_row = input$in_step22_wrap_row,
                                        wrap_scales = input$in_step22_wrap_scales)
    
    # relative pausing score
    relative_pausing_plot <- draw_codon(codon = relative_codon_table,
                                        x = "Codon",
                                        y = "PauseScore",
                                        xlab = "Codon",
                                        ylab = "Relative Pause Score",
                                        title = "Codon pausing score",
                                        fill_color = fill_color,
                                        color_palette = input$in_step22_dot_color,
                                        color_alpha = input$in_step22_alpha,
                                        dot_size = input$in_step22_dot_size,
                                        line_size = input$in_step22_line_size,
                                        
                                        grid_width = input$in_step22_grid_width,
                                        grid_color = input$in_step22_grid_color,
                                        
                                        font_size = input$in_step22_font_size,
                                        font_color = 'black',
                                        legend_row = input$in_step22_legend_row,
                                        
                                        wrap = input$in_step22_wrap,
                                        wrap_group = "AA",
                                        wrap_row = input$in_step22_wrap_row,
                                        wrap_scales = input$in_step22_wrap_scales)
    
    return(list(absolute_pausing_plot = absolute_pausing_plot,
                relative_pausing_plot = relative_pausing_plot))
  })
  
  ## output the figure
  observeEvent(input$act_step22_draw_pausing_score, {
    
    ## output the figure
    output$out_step22_abs_pausing_plot <- renderPlot(
      width = input$out_step22_fig_width * 100,
      height = input$out_step22_fig_height * 100,
      {step22_pausing_plot()$absolute_pausing_plot})
    
    output$out_step22_rel_pausing_plot <- renderPlot(
      width = input$out_step22_fig_width * 100,
      height = input$out_step22_fig_height * 100,
      {step22_pausing_plot()$relative_pausing_plot})
    
    output$save_step22_abs_pausing_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step22_out_name, "-absolute-pausing-score-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step22_pausing_plot()$absolute_pausing_plot, filename = file,
               width = input$out_step22_fig_width, height = input$out_step22_fig_height)
      }
    )
    
    output$save_step22_rel_pausing_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step22_out_name, "-relative-pausing-score-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step22_pausing_plot()$relative_pausing_plot, filename = file,
               width = input$out_step22_fig_width, height = input$out_step22_fig_height)
      }
    )
    
  })
  
  
  ### 22.3 import the pausing score table ##############################
  step22_raw_pausing_score <- reactive({
    req(input$act_step22_diff_import_pausing)
    
    if (is.null(input$in_step22_diff_pausing)) {return(NULL)}
    
    # import the raw data table
    pausing_score <- read.table(file = input$in_step22_diff_pausing$datapath, 
                                row.names = NULL, header = T, check.names = F, sep = '\t')
    
    return(pausing_score)
  })
  
  ## output the table
  observeEvent(input$act_step22_diff_import_pausing, {
    
    ## output the table
    output$out_step22_raw_pausing <- DT::renderDataTable(server = T, {
      DT::datatable(step22_raw_pausing_score(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "raw_ps_filter"), 
                                     list(extend = 'excel', filename = "raw_ps_filter"))
                    ))
    })
    
  })
  
  ### 22.4 calculate the differential pausing score ##############################
  step22_diff_pausing_score <- reactive({
    
    req(input$act_step22_diff_calculate)
    
    if (is.null(step22_raw_pausing_score)) {return(NULL)}
    
    # browser()
    pausing_score <- step22_raw_pausing_score()
    
    absolute_total_ps <- pausing_score %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_absolute_total_ps")) %>% 
      dplyr::rename_with(~str_remove(., "_absolute_total_ps"), contains("_absolute_total_ps")) %>% 
      dplyr::mutate(Class = "Total")
    absolute_valid_ps <- pausing_score %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_absolute_valid_ps")) %>% 
      dplyr::rename_with(~str_remove(., "_absolute_valid_ps"), contains("_absolute_valid_ps")) %>% 
      dplyr::mutate(Class = "Valid")
    absolute_ps <- rbind(absolute_valid_ps, absolute_total_ps)
    
    relative_total_pausing <- pausing_score %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_relative_total_ps")) %>% 
      dplyr::rename_with(~str_remove(., "_relative_total_ps"), contains("_relative_total_ps")) %>%
      dplyr::mutate(Class = "Total")
    relative_valid_pausing <- pausing_score %>%
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_relative_valid_ps")) %>% 
      dplyr::rename_with(~str_remove(., "_relative_valid_ps"), contains("_relative_valid_ps")) %>%
      dplyr::mutate(Class = "Valid")
    relative_ps <- rbind(relative_valid_pausing, relative_total_pausing)
    
    # browser()
    
    # calculate the differential pausing score
    if (nchar(input$in_step22_diff_column) > 0) {
      
      if (is.null(step1_design())) {
        return(list(absolute_pausing = absolute_ps,
                    relative_pausing = relative_ps))
      }
      
      # browser()
      
      absolute_pausing <- absolute_ps %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr', 'Class'), names_to = "Sample", values_to = "PauseScore") %>% 
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step22_diff_column) %in% c(input$in_step22_diff_groups1, input$in_step22_diff_groups2)) %>%
        dplyr::group_by(Codon, AA, Abbr, Class, !!sym(input$in_step22_diff_column)) %>% 
        dplyr::reframe(PauseScore = mean(PauseScore)) %>% 
        dplyr::ungroup() %>% 
        tidyr::pivot_wider(names_from = !!sym(input$in_step22_diff_column), values_from = PauseScore) %>% 
        mutate(Delta = !!sym(input$in_step22_diff_groups1) - !!sym(input$in_step22_diff_groups2),
               FoldChange = !!sym(input$in_step22_diff_groups1) / !!sym(input$in_step22_diff_groups2))
      
      relative_pausing <- relative_ps %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr', 'Class'), names_to = "Sample", values_to = "PauseScore") %>% 
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step22_diff_column) %in% c(input$in_step22_diff_groups1, input$in_step22_diff_groups2)) %>%
        dplyr::group_by(Codon, AA, Abbr, Class, !!sym(input$in_step22_diff_column)) %>% 
        dplyr::reframe(PauseScore = mean(PauseScore)) %>% 
        dplyr::ungroup() %>% 
        tidyr::pivot_wider(names_from = !!sym(input$in_step22_diff_column), values_from = PauseScore) %>% 
        mutate(Delta = !!sym(input$in_step22_diff_groups1) - !!sym(input$in_step22_diff_groups2),
               FoldChange = !!sym(input$in_step22_diff_groups1) / !!sym(input$in_step22_diff_groups2))
      
      return(list(absolute_pausing = absolute_pausing,
                  relative_pausing = relative_pausing))
      
    } else {
      
      return(list(absolute_pausing = absolute_ps,
                  relative_pausing = relative_ps))
      
    }
    
  })
  
  ## output the table
  observeEvent(input$act_step22_diff_import_pausing, {
    
    ## output the table
    output$out_step22_diff_absolute_pausing <- DT::renderDataTable(server = T, {
      DT::datatable(step22_diff_pausing_score()$absolute_pausing, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "absolute_diff_ps_filter"), 
                                     list(extend = 'excel', filename = "absolute_diff_ps_filter"))
                    ))
    })
    
    output$out_step22_diff_relative_pausing <- DT::renderDataTable(server = T, {
      DT::datatable(step22_diff_pausing_score()$relative_pausing, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "relative_diff_ps_filter"), 
                                     list(extend = 'excel', filename = "relative_diff_ps_filter"))
                    ))
    })
    
  })
  
  
  ### 22.4 draw the differential pausing score plot ##############################
  step22_diff_pausing_plot <- reactive({
    
    req(input$act_step22_diff_draw_pausing_score)
    
    if (is.null(step22_diff_pausing_score())) {return(NULL)}
    
    # browser()
    source("R/draw_codon.R")
    
    absolute_codon_table <- step22_diff_pausing_score()$absolute_pausing %>% 
      dplyr::filter(Class == input$in_step22_diff_valid)
    
    relative_codon_table <- step22_diff_pausing_score()$relative_pausing %>% 
      dplyr::filter(Class == input$in_step22_diff_valid)
    
    # set the ylabel
    if (input$in_step22_diff_method == "Delta") {
      ylab_group <- paste0(input$in_step22_diff_groups1, " - ", input$in_step22_diff_groups2)
    } else if (input$in_step22_diff_method == "FoldChange") {
      ylab_group <- paste0(input$in_step22_diff_groups1, " / ", input$in_step22_diff_groups2)
    } else {
      ylab_group <- input$in_step22_diff_method
    }
    
    # absolute pausing score
    absolute_pausing_plot <- draw_codon(codon = absolute_codon_table,
                                        x = "Codon",
                                        y = input$in_step22_diff_method,
                                        xlab = "Codon",
                                        ylab = ylab_group,
                                        title = "Codon pausing score",
                                        fill_color = "AA",
                                        color_palette = input$in_step22_diff_dot_color,
                                        color_alpha = input$in_step24_diff_alpha,
                                        
                                        dot_size = input$in_step22_diff_dot_size,
                                        line_size = input$in_step22_diff_line_size,
                                        grid_width = input$in_step22_diff_grid_width,
                                        grid_color = input$in_step22_diff_grid_color,
                                        
                                        font_size = input$in_step22_diff_font_size,
                                        font_color = 'black',
                                        legend_row = input$in_step22_diff_legend_row,
                                        
                                        wrap = input$in_step22_diff_wrap,
                                        wrap_group = "AA",
                                        wrap_row = input$in_step22_diff_wrap_row,
                                        wrap_scales = input$in_step22_diff_wrap_scales)
    
    # relative pausing score
    relative_pausing_plot <- draw_codon(codon = relative_codon_table,
                                        x = "Codon",
                                        y = input$in_step22_diff_method,
                                        xlab = "Codon",
                                        ylab = ylab_group,
                                        title = "Relative codon pausing score",
                                        fill_color = "AA",
                                        color_palette = input$in_step22_diff_dot_color,
                                        color_alpha = input$in_step24_diff_alpha,
                                        
                                        dot_size = input$in_step22_diff_dot_size,
                                        line_size = input$in_step22_diff_line_size,
                                        grid_width = input$in_step22_diff_grid_width,
                                        grid_color = input$in_step22_diff_grid_color,
                                        
                                        font_size = input$in_step22_diff_font_size,
                                        font_color = 'black',
                                        legend_row = input$in_step22_diff_legend_row,
                                        
                                        wrap = input$in_step22_diff_wrap,
                                        wrap_group = "AA",
                                        wrap_row = input$in_step22_diff_wrap_row,
                                        wrap_scales = input$in_step22_diff_wrap_scales)
    
    return(list(absolute_pausing_plot = absolute_pausing_plot,
                relative_pausing_plot = relative_pausing_plot))
  })
  
  ## output the figure
  observeEvent(input$act_step22_diff_draw_pausing_score, {
    
    ## output the figure
    output$out_step22_diff_abs_pausing_plot <- renderPlot(
      width = input$out_step22_diff_fig_width * 100,
      height = input$out_step22_diff_fig_height * 100,
      {step22_diff_pausing_plot()$absolute_pausing_plot})
    
    output$out_step22_diff_rel_pausing_plot <- renderPlot(
      width = input$out_step22_diff_fig_width * 100,
      height = input$out_step22_diff_fig_height * 100,
      {step22_diff_pausing_plot()$relative_pausing_plot})
    
    output$save_step22_diff_abs_pausing_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step22_diff_out_name, "-diff-absolute-pausing-score-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step22_diff_pausing_plot()$absolute_pausing_plot, filename = file,
               width = input$out_step22_diff_fig_width, height = input$out_step22_diff_fig_height)
      }
    )
    
    output$save_step22_diff_rel_pausing_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step22_diff_out_name, "-diff-relative-pausing-score-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step22_diff_pausing_plot()$relative_pausing_plot, filename = file,
               width = input$out_step22_diff_fig_width, height = input$out_step22_diff_fig_height)
      }
    )
    
  })
  
  ############################################################
  ## step23 Occupancy ###############################
  
  ### 23.0 set the design ##############################
  observeEvent(input$act_step1_import_design, {
    col_name <- names(step1_design())
    updateSelectizeInput(session, "in_step23_column", choices = col_name, selected = col_name[3], server = TRUE)
    updateSelectizeInput(session, "in_step23_diff_column", choices = col_name, selected = col_name[3], server = TRUE)
    
    # browser()
    observe({
      req(input$in_step23_column)
      step1_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      select_choices <- step1_design[[input$in_step23_column]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step23_groups", choices = select_choices, selected = select_choices, server = TRUE)
    })
    
    observe({
      req(input$in_step23_diff_column)
      step1_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      select_choices <- step1_design[[input$in_step23_diff_column]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step23_diff_groups1", choices = select_choices, selected = select_choices[1], server = TRUE)
      updateSelectizeInput(session, "in_step23_diff_groups2", choices = select_choices, selected = select_choices[2], server = TRUE)
    })
    
  })
  
  ### 23.1 import the codon occupancy data ##############################
  step23_occupancy <- reactive({
    
    req(input$act_step23_import_occupancy)
    
    if (is.null(input$in_step23_occupancy)) {return(NULL)}
    
    # import the raw data table
    occupancy <- read.table(file = input$in_step23_occupancy$datapath, 
                            row.names = NULL, header = T, check.names = F, sep = '\t')
    
    # browser()
    absolute_occup <- occupancy %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_absolute_occupancy")) %>% 
      dplyr::rename_with(~str_remove(., "_absolute_occupancy"), contains("_absolute_occupancy"))
    
    relative_occup <- occupancy %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_relative_occupancy")) %>% 
      dplyr::rename_with(~str_remove(., "_relative_occupancy"), contains("_relative_occupancy"))
    
    # calculate the average of the occupancy
    if (!is.null(input$in_step23_groups) & !is.null(input$in_step23_column)) {
      
      # browser()
      absolute_occupancy <- absolute_occup %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "Occupancy") %>%
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step23_column) %in% input$in_step23_groups) %>%
        dplyr::group_by(Codon, AA, Abbr, !!sym(input$in_step23_column)) %>% 
        dplyr::reframe(Occupancy = mean(Occupancy)) %>% 
        dplyr::ungroup()
      
      relative_occupancy <- relative_occup %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "Occupancy") %>%
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step23_column) %in% input$in_step23_groups) %>%
        dplyr::group_by(Codon, AA, Abbr, !!sym(input$in_step23_column)) %>% 
        dplyr::reframe(Occupancy = mean(Occupancy)) %>% 
        dplyr::ungroup()
      
      return(list(absolute_occupancy = absolute_occupancy,
                  relative_occupancy = relative_occupancy))
    } else {
      
      absolute_occupancy <- absolute_occup %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "Occupancy") 
      
      relative_occupancy <- relative_occup %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "Occupancy")
      
      return(list(absolute_occupancy = absolute_occupancy,
                  relative_occupancy = relative_occupancy))
    }
    
  })
  
  ## output the table
  observeEvent(input$act_step23_import_occupancy, {
    
    ## output the table
    output$out_step23_absolute_occupancy <- DT::renderDataTable(server = T, {
      DT::datatable(step23_occupancy()$absolute_occupancy, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "occupancy_filter"), 
                                     list(extend = 'excel', filename = "occupancy_filter"))
                    ))
    })
    
    output$out_step23_relative_occupancy <- DT::renderDataTable(server = T, {
      DT::datatable(step23_occupancy()$relative_occupancy, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "relative_occupancy_filter"), 
                                     list(extend = 'excel', filename = "relative_occupancy_filter"))
                    ))
    })
    
  })
  
  
  ### 23.2 draw the occupancy plot ##############################
  step23_occupancy_plot <- reactive({
    
    req(input$act_step23_draw_occupancy)
    
    if (is.null(step23_occupancy())) {return(NULL)}
    
    # browser()
    source("R/draw_codon.R")
    
    # set the dot color aes
    if (nchar(input$in_step23_column) == 0) {
      fill_color <- "Sample"
    } else {
      fill_color <- input$in_step23_column
    }
    
    # absolute occupancy score
    absolute_occupancy_plot <- draw_codon(codon = step23_occupancy()$absolute_occupancy,
                                          x = "Codon",
                                          y = "Occupancy",
                                          xlab = "Codon",
                                          ylab = "Occupancy",
                                          title = "Codon occupancy",
                                          fill_color = fill_color,
                                          color_palette = input$in_step23_dot_color,
                                          color_alpha = input$in_step23_alpha,
                                          
                                          dot_size = input$in_step23_dot_size,
                                          line_size = input$in_step23_line_size,
                                          grid_width = input$in_step23_grid_width,
                                          grid_color = input$in_step23_grid_color,
                                          
                                          font_size = input$in_step23_font_size,
                                          font_color = 'black',
                                          legend_row = input$in_step23_legend_row,
                                          
                                          wrap = input$in_step23_wrap,
                                          wrap_group = "AA",
                                          wrap_row = input$in_step23_wrap_row,
                                          wrap_scales = input$in_step23_wrap_scales)
    
    # relative occupancy score
    relative_occupancy_plot <- draw_codon(codon = step23_occupancy()$relative_occupancy,
                                          x = "Codon",
                                          y = "Occupancy",
                                          xlab = "Codon",
                                          ylab = "Relative Occupancy",
                                          title = "Codon occupancy",
                                          fill_color = fill_color,
                                          color_palette = input$in_step23_dot_color,
                                          color_alpha = input$in_step23_alpha,
                                          
                                          dot_size = input$in_step23_dot_size,
                                          line_size = input$in_step23_line_size,
                                          grid_width = input$in_step23_grid_width,
                                          grid_color = input$in_step23_grid_color,
                                          
                                          font_size = input$in_step23_font_size,
                                          font_color = 'black',
                                          legend_row = input$in_step23_legend_row,
                                          
                                          wrap = input$in_step23_wrap,
                                          wrap_group = "AA",
                                          wrap_row = input$in_step23_wrap_row,
                                          wrap_scales = input$in_step23_wrap_scales)
    
    return(list(absolute_occupancy_plot = absolute_occupancy_plot,
                relative_occupancy_plot = relative_occupancy_plot))
  })
  
  ## output the figure
  observeEvent(input$act_step23_draw_occupancy, {
    
    ## output the figure
    output$out_step23_absolute_occupancy_plot <- renderPlot(
      width = input$out_step23_fig_width * 100,
      height = input$out_step23_fig_height * 100,
      {step23_occupancy_plot()$absolute_occupancy_plot})
    
    output$out_step23_relative_occupancy_plot <- renderPlot(
      width = input$out_step23_fig_width * 100,
      height = input$out_step23_fig_height * 100,
      {step23_occupancy_plot()$relative_occupancy_plot})
    
    output$save_step23_absolute_occupancy_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step23_out_name, "-occupancy-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step23_occupancy_plot()$absolute_occupancy_plot, filename = file,
               width = input$out_step23_fig_width, height = input$out_step23_fig_height)
      }
    )
    
    output$save_step23_relative_occupancy_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step23_out_name, "-relative-occupancy-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step23_occupancy_plot()$relative_occupancy_plot, filename = file,
               width = input$out_step23_fig_width, height = input$out_step23_fig_height)
      }
    )
    
  })
  
  
  ### 23.3 import the occupancy table ##############################
  step23_raw_occupancy <- reactive({
    req(input$act_step23_diff_import_occupancy)
    
    if (is.null(input$in_step23_diff_occupancy)) {return(NULL)}
    
    # import the raw data table
    occupancy <- read.table(file = input$in_step23_diff_occupancy$datapath, 
                            row.names = NULL, header = T, check.names = F, sep = '\t')
    
    return(occupancy)
  })
  
  ## output the table
  observeEvent(input$act_step23_diff_import_occupancy, {
    
    ## output the table
    output$out_step23_raw_occupancy <- DT::renderDataTable(server = T, {
      DT::datatable(step23_raw_occupancy(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "occupancy_filter"), 
                                     list(extend = 'excel', filename = "occupancy_filter"))
                    ))
    })
    
  })
  
  ### 23.4 calculate the differential occupancy ##############################
  step23_diff_occupancy <- reactive({
    
    req(input$act_step23_diff_calculate)
    
    if (is.null(step23_raw_occupancy)) {return(NULL)}
    
    # browser()
    occupancy <- step23_raw_occupancy()
    
    absolute_occup <- occupancy %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_absolute_occupancy")) %>% 
      dplyr::rename_with(~str_remove(., "_absolute_occupancy"), contains("_absolute_occupancy"))
    
    relative_occup <- occupancy %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_relative_occupancy")) %>% 
      dplyr::rename_with(~str_remove(., "_relative_occupancy"), contains("_relative_occupancy"))
    
    # browser()
    
    # calculate the differential occupancy score
    if (nchar(input$in_step23_diff_column) > 0) {
      
      if (is.null(step1_design())) {
        return(list(absolute_occupancy = absolute_occup,
                    relative_occupancy = relative_occup))
      }
      
      # browser()
      
      absolute_occupancy <- absolute_occup %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "Occupancy") %>% 
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step23_diff_column) %in% c(input$in_step23_diff_groups1, input$in_step23_diff_groups2)) %>%
        dplyr::group_by(Codon, AA, Abbr, !!sym(input$in_step23_diff_column)) %>% 
        dplyr::reframe(Occupancy = mean(Occupancy)) %>% 
        dplyr::ungroup() %>% 
        tidyr::pivot_wider(names_from = !!sym(input$in_step23_diff_column), values_from = Occupancy) %>% 
        mutate(Delta = !!sym(input$in_step23_diff_groups1) - !!sym(input$in_step23_diff_groups2),
               FoldChange = !!sym(input$in_step23_diff_groups1) / !!sym(input$in_step23_diff_groups2))
      
      relative_occupancy <- relative_occup %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "Occupancy") %>% 
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step23_diff_column) %in% c(input$in_step23_diff_groups1, input$in_step23_diff_groups2)) %>%
        dplyr::group_by(Codon, AA, Abbr, !!sym(input$in_step23_diff_column)) %>% 
        dplyr::reframe(Occupancy = mean(Occupancy)) %>% 
        dplyr::ungroup() %>% 
        tidyr::pivot_wider(names_from = !!sym(input$in_step23_diff_column), values_from = Occupancy) %>% 
        mutate(Delta = !!sym(input$in_step23_diff_groups1) - !!sym(input$in_step23_diff_groups2),
               FoldChange = !!sym(input$in_step23_diff_groups1) / !!sym(input$in_step23_diff_groups2))
      
      return(list(absolute_occupancy = absolute_occupancy,
                  relative_occupancy = relative_occupancy))
      
    } else {
      
      return(list(absolute_occupancy = absolute_occup,
                  relative_occupancy = relative_occup))
      
    }
    
  })
  
  ## output the table
  observeEvent(input$act_step23_diff_import_occupancy, {
    
    ## output the table
    output$out_step23_diff_absolute_occupancy <- DT::renderDataTable(server = T, {
      DT::datatable(step23_diff_occupancy()$absolute_occupancy, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "absolute_diff_occupancy_filter"), 
                                     list(extend = 'excel', filename = "absolute_diff_occupancy_filter"))
                    ))
    })
    
    output$out_step23_diff_relative_occupancy <- DT::renderDataTable(server = T, {
      DT::datatable(step23_diff_occupancy()$relative_occupancy, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "relative_diff_occupancy_filter"), 
                                     list(extend = 'excel', filename = "relative_diff_occupancy_filter"))
                    ))
    })
    
  })
  
  
  ### 23.4 draw the differential occupancy plot ##############################
  step23_diff_occupancy_plot <- reactive({
    
    req(input$act_step23_diff_draw_occupancy)
    
    if (is.null(step23_diff_occupancy())) {return(NULL)}
    
    source("R/draw_codon.R")
    
    # browser()
    
    # set the ylabel
    if (input$in_step23_diff_method == "Delta") {
      ylab_group <- paste0(input$in_step23_diff_groups1, " - ", input$in_step23_diff_groups2)
    } else if (input$in_step23_diff_method == "FoldChange") {
      ylab_group <- paste0(input$in_step23_diff_groups1, " / ", input$in_step23_diff_groups2)
    } else {
      ylab_group <- input$in_step23_diff_method
    }
    
    # absolute occupancy
    absolute_occupancy_plot <- draw_codon(codon = step23_diff_occupancy()$absolute_occupancy,
                                          x = "Codon",
                                          y = input$in_step23_diff_method,
                                          xlab = "Codon",
                                          ylab = ylab_group,
                                          title = "Codon occupancy",
                                          fill_color = "AA",
                                          color_palette = input$in_step23_diff_dot_color,
                                          color_alpha = input$in_step23_diff_alpha,
                                          
                                          dot_size = input$in_step23_diff_dot_size,
                                          line_size = input$in_step23_diff_line_size,
                                          grid_width = input$in_step23_diff_grid_width,
                                          grid_color = input$in_step23_diff_grid_color,
                                          
                                          font_size = input$in_step23_diff_font_size,
                                          font_color = 'black',
                                          legend_row = input$in_step23_diff_legend_row,
                                          
                                          wrap = input$in_step23_diff_wrap,
                                          wrap_group = "AA",
                                          wrap_row = input$in_step23_diff_wrap_row,
                                          wrap_scales = input$in_step23_diff_wrap_scales)
    
    # relative occupancy
    relative_occupancy_plot <- draw_codon(codon = step23_diff_occupancy()$relative_occupancy,
                                          x = "Codon",
                                          y = input$in_step23_diff_method,
                                          xlab = "Codon",
                                          ylab = ylab_group,
                                          title = "Relative codon occupancy",
                                          fill_color = "AA",
                                          color_palette = input$in_step23_diff_dot_color,
                                          color_alpha = input$in_step23_diff_alpha,
                                          
                                          dot_size = input$in_step23_diff_dot_size,
                                          line_size = input$in_step23_diff_line_size,
                                          grid_width = input$in_step23_diff_grid_width,
                                          grid_color = input$in_step23_diff_grid_color,
                                          
                                          font_size = input$in_step23_diff_font_size,
                                          font_color = 'black',
                                          legend_row = input$in_step23_diff_legend_row,
                                          
                                          wrap = input$in_step23_diff_wrap,
                                          wrap_group = "AA",
                                          wrap_row = input$in_step23_diff_wrap_row,
                                          wrap_scales = input$in_step23_diff_wrap_scales)
    
    return(list(absolute_occupancy_plot = absolute_occupancy_plot,
                relative_occupancy_plot = relative_occupancy_plot))
  })
  
  ## output the figure
  observeEvent(input$act_step23_diff_draw_occupancy, {
    
    ## output the figure
    output$out_step23_diff_abs_occupancy_plot <- renderPlot(
      width = input$out_step23_diff_fig_width * 100,
      height = input$out_step23_diff_fig_height * 100,
      {step23_diff_occupancy_plot()$absolute_occupancy_plot})
    
    output$out_step23_diff_rel_occupancy_plot <- renderPlot(
      width = input$out_step23_diff_fig_width * 100,
      height = input$out_step23_diff_fig_height * 100,
      {step23_diff_occupancy_plot()$relative_occupancy_plot})
    
    output$save_step23_diff_abs_occupancy_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step23_diff_out_name, "-diff-absolute-occupancy-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step23_diff_occupancy_plot()$absolute_occupancy_plot, filename = file,
               width = input$out_step23_diff_fig_width, height = input$out_step23_diff_fig_height)
      }
    )
    
    output$save_step23_diff_rel_occupancy_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step23_diff_out_name, "-diff-relative-occupancy-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step23_diff_occupancy_plot()$relative_occupancy_plot, filename = file,
               width = input$out_step23_diff_fig_width, height = input$out_step23_diff_fig_height)
      }
    )
    
  })
  
  ############################################################
  ## step24 CDT ###############################
  
  ### 24.0 set the design ##############################
  observeEvent(input$act_step1_import_design, {
    col_name <- names(step1_design())
    updateSelectizeInput(session, "in_step24_column", choices = col_name, selected = col_name[3], server = TRUE)
    updateSelectizeInput(session, "in_step24_diff_column", choices = col_name, selected = col_name[3], server = TRUE)
    
    # browser()
    observe({
      req(input$in_step24_column)
      step1_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      select_choices <- step1_design[[input$in_step24_column]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step24_groups", choices = select_choices, selected = select_choices, server = TRUE)
    })
    
    observe({
      req(input$in_step24_diff_column)
      step1_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      select_choices <- step1_design[[input$in_step24_diff_column]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step24_diff_groups1", choices = select_choices, selected = select_choices[1], server = TRUE)
      updateSelectizeInput(session, "in_step24_diff_groups2", choices = select_choices, selected = select_choices[2], server = TRUE)
    })
    
  })
  
  ### 24.1 import the codon cdt data ##############################
  step24_cdt <- reactive({
    
    req(input$act_step24_import_cdt)
    
    if (is.null(input$in_step24_cdt)) {return(NULL)}
    
    # import the raw data table
    cdt <- read.table(file = input$in_step24_cdt$datapath, 
                      row.names = NULL, header = T, check.names = F, sep = '\t')
    
    # browser()
    
    # filter the raw or norm data
    if (input$in_step24_norm == "Raw") {
      cdt <- cdt %>% dplyr::select(!contains("_norm"))
    } else if (input$in_step24_norm == "Normlized") {
      cdt <- cdt %>% dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_norm"))
    }
    
    # retrieve the absolute and relative cdt
    absolute_cdt <- cdt %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), !contains("relative")) %>% 
      dplyr::rename_with(~str_remove(., "_absolute_cdt"), contains("_absolute_cdt")) %>% 
      dplyr::rename_with(~str_remove(., "_norm_cdt"), contains("_norm_cdt"))
    
    relative_cdt <- cdt %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_relative_cdt")) %>% 
      dplyr::rename_with(~str_remove(., "_norm_relative_cdt"), contains("_norm_relative_cdt")) %>% 
      dplyr::rename_with(~str_remove(., "_relative_cdt"), contains("_relative_cdt"))
    
    # calculate the average of the cdt
    if (!is.null(input$in_step24_groups) & !is.null(input$in_step24_column)) {
      
      # browser()
      absolute_cdt <- absolute_cdt %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "CDT") %>%
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step24_column) %in% input$in_step24_groups) %>%
        dplyr::group_by(Codon, AA, Abbr, !!sym(input$in_step24_column)) %>% 
        dplyr::reframe(CDT = mean(CDT)) %>% 
        dplyr::ungroup()
      
      relative_cdt <- relative_cdt %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "CDT") %>%
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step24_column) %in% input$in_step24_groups) %>%
        dplyr::group_by(Codon, AA, Abbr, !!sym(input$in_step24_column)) %>% 
        dplyr::reframe(CDT = mean(CDT)) %>% 
        dplyr::ungroup()
      
      return(list(absolute_cdt = absolute_cdt,
                  relative_cdt = relative_cdt))
    } else {
      
      absolute_cdt <- absolute_cdt %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "CDT") 
      
      relative_cdt <- relative_cdt %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "CDT")
      
      return(list(absolute_cdt = absolute_cdt,
                  relative_cdt = relative_cdt))
    }
    
  })
  
  ## output the table
  observeEvent(input$act_step24_import_cdt, {
    
    ## output the table
    output$out_step24_absolute_cdt <- DT::renderDataTable(server = T, {
      DT::datatable(step24_cdt()$absolute_cdt, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "cdt_filter"), 
                                     list(extend = 'excel', filename = "cdt_filter"))
                    ))
    })
    
    output$out_step24_relative_cdt <- DT::renderDataTable(server = T, {
      DT::datatable(step24_cdt()$relative_cdt, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "relative_cdt_filter"), 
                                     list(extend = 'excel', filename = "relative_cdt_filter"))
                    ))
    })
    
  })
  
  
  ### 24.2 draw the cdt plot ##############################
  step24_cdt_plot <- reactive({
    
    req(input$act_step24_draw_cdt)
    
    if (is.null(step24_cdt())) {return(NULL)}
    
    # browser()
    source("R/draw_codon.R")
    
    # set the dot color aes
    if (nchar(input$in_step24_column) == 0) {
      fill_color <- "Sample"
    } else {
      fill_color <- input$in_step24_column
    }
    
    # absolute cdt score
    absolute_cdt_plot <- draw_codon(codon = step24_cdt()$absolute_cdt,
                                    x = "Codon",
                                    y = "CDT",
                                    xlab = "Codon",
                                    ylab = "CDT",
                                    title = "Codon decoding time",
                                    fill_color = fill_color,
                                    color_palette = input$in_step24_dot_color,
                                    color_alpha = input$in_step24_alpha,
                                    
                                    dot_size = input$in_step24_dot_size,
                                    line_size = input$in_step24_line_size,
                                    grid_width = input$in_step24_grid_width,
                                    grid_color = input$in_step24_grid_color,
                                    
                                    font_size = input$in_step24_font_size,
                                    font_color = 'black',
                                    legend_row = input$in_step24_legend_row,
                                    
                                    wrap = input$in_step24_wrap,
                                    wrap_group = "AA",
                                    wrap_row = input$in_step24_wrap_row,
                                    wrap_scales = input$in_step24_wrap_scales)
    
    # relative cdt score
    relative_cdt_plot <- draw_codon(codon = step24_cdt()$relative_cdt,
                                    x = "Codon",
                                    y = "CDT",
                                    xlab = "Codon",
                                    ylab = "Relative CDT",
                                    title = "Codon decoding time",
                                    fill_color = fill_color,
                                    color_palette = input$in_step24_dot_color,
                                    color_alpha = input$in_step24_diff_alpha,
                                    
                                    dot_size = input$in_step24_dot_size,
                                    line_size = input$in_step24_line_size,
                                    grid_width = input$in_step24_grid_width,
                                    grid_color = input$in_step24_grid_color,
                                    
                                    font_size = input$in_step24_font_size,
                                    font_color = 'black',
                                    legend_row = input$in_step24_legend_row,
                                    
                                    wrap = input$in_step24_wrap,
                                    wrap_group = "AA",
                                    wrap_row = input$in_step24_wrap_row,
                                    wrap_scales = input$in_step24_wrap_scales)
    
    return(list(absolute_cdt_plot = absolute_cdt_plot,
                relative_cdt_plot = relative_cdt_plot))
  })
  
  ## output the figure
  observeEvent(input$act_step24_draw_cdt, {
    
    ## output the figure
    output$out_step24_absolute_cdt_plot <- renderPlot(
      width = input$out_step24_fig_width * 100,
      height = input$out_step24_fig_height * 100,
      {step24_cdt_plot()$absolute_cdt_plot})
    
    output$out_step24_relative_cdt_plot <- renderPlot(
      width = input$out_step24_fig_width * 100,
      height = input$out_step24_fig_height * 100,
      {step24_cdt_plot()$relative_cdt_plot})
    
    output$save_step24_absolute_cdt_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step24_out_name, "-cdt-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step24_cdt_plot()$absolute_cdt_plot, filename = file,
               width = input$out_step24_fig_width, height = input$out_step24_fig_height)
      }
    )
    
    output$save_step24_relative_cdt_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step24_out_name, "-relative-cdt-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step24_cdt_plot()$relative_cdt_plot, filename = file,
               width = input$out_step24_fig_width, height = input$out_step24_fig_height)
      }
    )
    
  })
  
  
  ### 24.3 import the cdt table ##############################
  step24_raw_cdt <- reactive({
    req(input$act_step24_diff_import_cdt)
    
    if (is.null(input$in_step24_diff_cdt)) {return(NULL)}
    
    # import the raw data table
    cdt <- read.table(file = input$in_step24_diff_cdt$datapath, 
                      row.names = NULL, header = T, check.names = F, sep = '\t')
    
    return(cdt)
  })
  
  ## output the table
  observeEvent(input$act_step24_diff_import_cdt, {
    
    ## output the table
    output$out_step24_raw_cdt <- DT::renderDataTable(server = T, {
      DT::datatable(step24_raw_cdt(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "cdt_filter"), 
                                     list(extend = 'excel', filename = "cdt_filter"))
                    ))
    })
    
  })
  
  ### 24.4 calculate the differential cdt ##############################
  step24_diff_cdt <- reactive({
    
    req(input$act_step24_diff_calculate)
    
    if (is.null(step24_raw_cdt)) {return(NULL)}
    
    # browser()
    cdt <- step24_raw_cdt()
    
    absolute_cdt <- cdt %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_absolute_cdt")) %>% 
      dplyr::rename_with(~str_remove(., "_absolute_cdt"), contains("_absolute_cdt"))
    
    relative_cdt <- cdt %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_relative_cdt")) %>% 
      dplyr::rename_with(~str_remove(., "_relative_cdt"), contains("_relative_cdt"))
    
    # browser()
    
    # calculate the differential cdt score
    if (nchar(input$in_step24_diff_column) > 0) {
      
      if (is.null(step1_design())) {
        return(list(absolute_cdt = absolute_cdt,
                    relative_cdt = relative_cdt))
      }
      
      # browser()
      
      absolute_cdt <- absolute_cdt %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "cdt") %>% 
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step24_diff_column) %in% c(input$in_step24_diff_groups1, input$in_step24_diff_groups2)) %>%
        dplyr::group_by(Codon, AA, Abbr, !!sym(input$in_step24_diff_column)) %>% 
        dplyr::reframe(cdt = mean(cdt)) %>% 
        dplyr::ungroup() %>% 
        tidyr::pivot_wider(names_from = !!sym(input$in_step24_diff_column), values_from = cdt) %>% 
        mutate(Delta = !!sym(input$in_step24_diff_groups1) - !!sym(input$in_step24_diff_groups2),
               FoldChange = !!sym(input$in_step24_diff_groups1) / !!sym(input$in_step24_diff_groups2))
      
      relative_cdt <- relative_cdt %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "cdt") %>% 
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step24_diff_column) %in% c(input$in_step24_diff_groups1, input$in_step24_diff_groups2)) %>%
        dplyr::group_by(Codon, AA, Abbr, !!sym(input$in_step24_diff_column)) %>% 
        dplyr::reframe(cdt = mean(cdt)) %>% 
        dplyr::ungroup() %>% 
        tidyr::pivot_wider(names_from = !!sym(input$in_step24_diff_column), values_from = cdt) %>% 
        mutate(Delta = !!sym(input$in_step24_diff_groups1) - !!sym(input$in_step24_diff_groups2),
               FoldChange = !!sym(input$in_step24_diff_groups1) / !!sym(input$in_step24_diff_groups2))
      
      return(list(absolute_cdt = absolute_cdt,
                  relative_cdt = relative_cdt))
      
    } else {
      
      return(list(absolute_cdt = absolute_cdt,
                  relative_cdt = relative_cdt))
      
    }
    
  })
  
  ## output the table
  observeEvent(input$act_step24_diff_import_cdt, {
    
    ## output the table
    output$out_step24_diff_absolute_cdt <- DT::renderDataTable(server = T, {
      DT::datatable(step24_diff_cdt()$absolute_cdt, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "absolute_diff_cdt_filter"), 
                                     list(extend = 'excel', filename = "absolute_diff_cdt_filter"))
                    ))
    })
    
    output$out_step24_diff_relative_cdt <- DT::renderDataTable(server = T, {
      DT::datatable(step24_diff_cdt()$relative_cdt, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "relative_diff_cdt_filter"), 
                                     list(extend = 'excel', filename = "relative_diff_cdt_filter"))
                    ))
    })
    
  })
  
  
  ### 24.4 draw the differential cdt plot ##############################
  step24_diff_cdt_plot <- reactive({
    
    req(input$act_step24_diff_draw_cdt)
    
    if (is.null(step24_diff_cdt())) {return(NULL)}
    
    source("R/draw_codon.R")
    
    # browser()
    
    # absolute cdt
    absolute_cdt_plot <- draw_codon(codon = step24_diff_cdt()$absolute_cdt,
                                    x = "Codon",
                                    y = input$in_step24_diff_method,
                                    xlab = "Codon",
                                    ylab = "cdt",
                                    title = "Codon cdt",
                                    fill_color = "AA",
                                    color_palette = input$in_step24_diff_dot_color,
                                    color_alpha = input$in_step24_diff_alpha,
                                    
                                    dot_size = input$in_step24_diff_dot_size,
                                    line_size = input$in_step24_diff_line_size,
                                    grid_width = input$in_step24_diff_grid_width,
                                    grid_color = input$in_step24_diff_grid_color,
                                    
                                    font_size = input$in_step24_diff_font_size,
                                    font_color = 'black',
                                    legend_row = input$in_step24_diff_legend_row,
                                    
                                    wrap = input$in_step24_diff_wrap,
                                    wrap_group = "AA",
                                    wrap_row = input$in_step24_diff_wrap_row,
                                    wrap_scales = input$in_step24_diff_wrap_scales)
    
    # relative cdt
    relative_cdt_plot <- draw_codon(codon = step24_diff_cdt()$relative_cdt,
                                    x = "Codon",
                                    y = input$in_step24_diff_method,
                                    xlab = "Codon",
                                    ylab = "Relative cdt",
                                    title = "Codon cdt",
                                    fill_color = "AA",,
                                    color_palette = input$in_step24_diff_dot_color,
                                    color_alpha = input$in_step24_diff_alpha,
                                    
                                    dot_size = input$in_step24_diff_dot_size,
                                    line_size = input$in_step24_diff_line_size,
                                    grid_width = input$in_step24_diff_grid_width,
                                    grid_color = input$in_step24_diff_grid_color,
                                    
                                    font_size = input$in_step24_diff_font_size,
                                    font_color = 'black',
                                    legend_row = input$in_step23_diff_legend_row,
                                    
                                    wrap = input$in_step24_diff_wrap,
                                    wrap_group = "AA",
                                    wrap_row = input$in_step24_diff_wrap_row,
                                    wrap_scales = input$in_step24_diff_wrap_scales)
    
    return(list(absolute_cdt_plot = absolute_cdt_plot,
                relative_cdt_plot = relative_cdt_plot))
  })
  
  ## output the figure
  observeEvent(input$act_step24_diff_draw_cdt, {
    
    ## output the figure
    output$out_step24_diff_abs_cdt_plot <- renderPlot(
      width = input$out_step24_diff_fig_width * 100,
      height = input$out_step24_diff_fig_height * 100,
      {step24_diff_cdt_plot()$absolute_cdt_plot})
    
    output$out_step24_diff_rel_cdt_plot <- renderPlot(
      width = input$out_step24_diff_fig_width * 100,
      height = input$out_step24_diff_fig_height * 100,
      {step24_diff_cdt_plot()$relative_cdt_plot})
    
    output$save_step24_diff_abs_cdt_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step24_diff_out_name, "-diff-absolute-cdt-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step24_diff_cdt_plot()$absolute_cdt_plot, filename = file,
               width = input$out_step24_diff_fig_width, height = input$out_step24_diff_fig_height)
      }
    )
    
    output$save_step24_diff_rel_cdt_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step24_diff_out_name, "-diff-relative-cdt-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step24_diff_cdt_plot()$relative_cdt_plot, filename = file,
               width = input$out_step24_diff_fig_width, height = input$out_step24_diff_fig_height)
      }
    )
    
  })
  
  
  ############################################################
  ## step25 CST ###############################
  
  ### 25.0 set the design ##############################
  observeEvent(input$act_step1_import_design, {
    col_name <- names(step1_design())
    updateSelectizeInput(session, "in_step25_column", choices = col_name, selected = col_name[3], server = TRUE)
    updateSelectizeInput(session, "in_step25_diff_column", choices = col_name, selected = col_name[3], server = TRUE)
    
    # browser()
    observe({
      req(input$in_step25_column)
      step1_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      select_choices <- step1_design[[input$in_step25_column]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step25_groups", choices = select_choices, selected = select_choices, server = TRUE)
    })
    
    observe({
      req(input$in_step25_diff_column)
      step1_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      select_choices <- step1_design[[input$in_step25_diff_column]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step25_diff_groups1", choices = select_choices, selected = select_choices[1], server = TRUE)
      updateSelectizeInput(session, "in_step25_diff_groups2", choices = select_choices, selected = select_choices[2], server = TRUE)
    })
    
  })
  
  ### 25.1 import the codon cst data ##############################
  step25_cst <- reactive({
    
    req(input$act_step25_import_cst)
    
    if (is.null(input$in_step25_cst)) {return(NULL)}
    
    # import the raw data table
    cst <- read.table(file = input$in_step25_cst$datapath, 
                      row.names = NULL, header = T, check.names = F, sep = '\t')
    
    # browser()
    absolute_cst <- cst %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_absolute_cst")) %>% 
      dplyr::rename_with(~str_remove(., "_absolute_cst"), contains("_absolute_cst"))
    
    relative_cst <- cst %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_relative_cst")) %>% 
      dplyr::rename_with(~str_remove(., "_relative_cst"), contains("_relative_cst"))
    
    # calculate the average of the cst
    if (!is.null(input$in_step25_groups) & !is.null(input$in_step25_column)) {
      
      # browser()
      absolute_cst <- absolute_cst %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "CST") %>%
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step25_column) %in% input$in_step25_groups) %>%
        dplyr::group_by(Codon, AA, Abbr, !!sym(input$in_step25_column)) %>% 
        dplyr::reframe(CST = mean(CST)) %>% 
        dplyr::ungroup()
      
      relative_cst <- relative_cst %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "CST") %>%
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step25_column) %in% input$in_step25_groups) %>%
        dplyr::group_by(Codon, AA, Abbr, !!sym(input$in_step25_column)) %>% 
        dplyr::reframe(CST = mean(CST)) %>% 
        dplyr::ungroup()
      
      return(list(absolute_cst = absolute_cst,
                  relative_cst = relative_cst))
    } else {
      
      absolute_cst <- absolute_cst %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "CST") 
      
      relative_cst <- relative_cst %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "CST")
      
      return(list(absolute_cst = absolute_cst,
                  relative_cst = relative_cst))
    }
    
  })
  
  ## output the table
  observeEvent(input$act_step25_import_cst, {
    
    ## output the table
    output$out_step25_absolute_cst <- DT::renderDataTable(server = T, {
      DT::datatable(step25_cst()$absolute_cst, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "cst_filter"), 
                                     list(extend = 'excel', filename = "cst_filter"))
                    ))
    })
    
    output$out_step25_relative_cst <- DT::renderDataTable(server = T, {
      DT::datatable(step25_cst()$relative_cst, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "relative_cst_filter"), 
                                     list(extend = 'excel', filename = "relative_cst_filter"))
                    ))
    })
    
  })
  
  
  ### 25.2 draw the cst plot ##############################
  step25_cst_plot <- reactive({
    
    req(input$act_step25_draw_cst)
    
    if (is.null(step25_cst())) {return(NULL)}
    
    # browser()
    source("R/draw_codon.R")
    
    # set the dot color aes
    if (nchar(input$in_step25_column) == 0) {
      fill_color <- "Sample"
    } else {
      fill_color <- input$in_step25_column
    }
    
    # absolute cst score
    absolute_cst_plot <- draw_codon(codon = step25_cst()$absolute_cst,
                                    x = "Codon",
                                    y = "CST",
                                    xlab = "Codon",
                                    ylab = "CST",
                                    title = "Codon selection time",
                                    fill_color = fill_color,
                                    color_palette = input$in_step25_dot_color,
                                    color_alpha = input$in_step25_alpha,
                                    
                                    dot_size = input$in_step25_dot_size,
                                    line_size = input$in_step25_line_size,
                                    grid_width = input$in_step25_grid_width,
                                    grid_color = input$in_step25_grid_color,
                                    
                                    font_size = input$in_step25_font_size,
                                    font_color = 'black',
                                    legend_row = input$in_step25_legend_row,
                                    
                                    wrap = input$in_step25_wrap,
                                    wrap_group = "AA",
                                    wrap_row = input$in_step25_wrap_row,
                                    wrap_scales = input$in_step25_wrap_scales)
    
    # relative cst score
    relative_cst_plot <- draw_codon(codon = step25_cst()$relative_cst,
                                    x = "Codon",
                                    y = "CST",
                                    xlab = "Codon",
                                    ylab = "CST",
                                    title = "Codon selection time",
                                    fill_color = fill_color,
                                    color_palette = input$in_step25_dot_color,
                                    color_alpha = input$in_step25_alpha,
                                    
                                    dot_size = input$in_step25_dot_size,
                                    line_size = input$in_step25_line_size,
                                    grid_width = input$in_step25_grid_width,
                                    grid_color = input$in_step25_grid_color,
                                    
                                    font_size = input$in_step25_font_size,
                                    font_color = 'black',
                                    legend_row = input$in_step25_legend_row,
                                    
                                    wrap = input$in_step25_wrap,
                                    wrap_group = "AA",
                                    wrap_row = input$in_step25_wrap_row,
                                    wrap_scales = input$in_step25_wrap_scales)
    
    return(list(absolute_cst_plot = absolute_cst_plot,
                relative_cst_plot = relative_cst_plot))
  })
  
  ## output the figure
  observeEvent(input$act_step25_draw_cst, {
    
    ## output the figure
    output$out_step25_absolute_cst_plot <- renderPlot(
      width = input$out_step25_fig_width * 100,
      height = input$out_step25_fig_height * 100,
      {step25_cst_plot()$absolute_cst_plot})
    
    output$out_step25_relative_cst_plot <- renderPlot(
      width = input$out_step25_fig_width * 100,
      height = input$out_step25_fig_height * 100,
      {step25_cst_plot()$relative_cst_plot})
    
    output$save_step25_absolute_cst_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step25_out_name, "-cst-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step25_cst_plot()$absolute_cst_plot, filename = file,
               width = input$out_step25_fig_width, height = input$out_step25_fig_height)
      }
    )
    
    output$save_step25_relative_cst_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step25_out_name, "-relative-cst-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step25_cst_plot()$relative_cst_plot, filename = file,
               width = input$out_step25_fig_width, height = input$out_step25_fig_height)
      }
    )
    
  })
  
  
  ### 25.3 import the cst table ##############################
  step25_raw_cst <- reactive({
    req(input$act_step25_diff_import_cst)
    
    if (is.null(input$in_step25_diff_cst)) {return(NULL)}
    
    # import the raw data table
    cst <- read.table(file = input$in_step25_diff_cst$datapath, 
                      row.names = NULL, header = T, check.names = F, sep = '\t')
    
    return(cst)
  })
  
  ## output the table
  observeEvent(input$act_step25_diff_import_cst, {
    
    ## output the table
    output$out_step25_raw_cst <- DT::renderDataTable(server = T, {
      DT::datatable(step25_raw_cst(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "cst_filter"), 
                                     list(extend = 'excel', filename = "cst_filter"))
                    ))
    })
    
  })
  
  ### 25.4 calculate the differential cst ##############################
  step25_diff_cst <- reactive({
    
    req(input$act_step25_diff_calculate)
    
    if (is.null(step25_raw_cst)) {return(NULL)}
    
    # browser()
    cst <- step25_raw_cst()
    
    absolute_cst <- cst %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_absolute_cst")) %>% 
      dplyr::rename_with(~str_remove(., "_absolute_cst"), contains("_absolute_cst"))
    
    relative_cst <- cst %>% 
      dplyr::select(c('Codon', 'AA', 'Abbr'), contains("_relative_cst")) %>% 
      dplyr::rename_with(~str_remove(., "_relative_cst"), contains("_relative_cst"))
    
    # browser()
    
    # calculate the differential cst score
    if (nchar(input$in_step25_diff_column) > 0) {
      
      if (is.null(step1_design())) {
        return(list(absolute_cst = absolute_cst,
                    relative_cst = relative_cst))
      }
      
      # browser()
      
      absolute_cst <- absolute_cst %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "cst") %>% 
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step25_diff_column) %in% c(input$in_step25_diff_groups1, input$in_step25_diff_groups2)) %>%
        dplyr::group_by(Codon, AA, Abbr, !!sym(input$in_step25_diff_column)) %>% 
        dplyr::reframe(cst = mean(cst)) %>% 
        dplyr::ungroup() %>% 
        tidyr::pivot_wider(names_from = !!sym(input$in_step25_diff_column), values_from = cst) %>% 
        mutate(Delta = !!sym(input$in_step25_diff_groups1) - !!sym(input$in_step25_diff_groups2),
               FoldChange = !!sym(input$in_step25_diff_groups1) / !!sym(input$in_step25_diff_groups2))
      
      relative_cst <- relative_cst %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), names_to = "Sample", values_to = "cst") %>% 
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step25_diff_column) %in% c(input$in_step25_diff_groups1, input$in_step25_diff_groups2)) %>%
        dplyr::group_by(Codon, AA, Abbr, !!sym(input$in_step25_diff_column)) %>% 
        dplyr::reframe(cst = mean(cst)) %>% 
        dplyr::ungroup() %>% 
        tidyr::pivot_wider(names_from = !!sym(input$in_step25_diff_column), values_from = cst) %>% 
        mutate(Delta = !!sym(input$in_step25_diff_groups1) - !!sym(input$in_step25_diff_groups2),
               FoldChange = !!sym(input$in_step25_diff_groups1) / !!sym(input$in_step25_diff_groups2))
      
      return(list(absolute_cst = absolute_cst,
                  relative_cst = relative_cst))
      
    } else {
      
      return(list(absolute_cst = absolute_cst,
                  relative_cst = relative_cst))
      
    }
    
  })
  
  ## output the table
  observeEvent(input$act_step25_diff_import_cst, {
    
    ## output the table
    output$out_step25_diff_absolute_cst <- DT::renderDataTable(server = T, {
      DT::datatable(step25_diff_cst()$absolute_cst, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "absolute_diff_cst_filter"), 
                                     list(extend = 'excel', filename = "absolute_diff_cst_filter"))
                    ))
    })
    
    output$out_step25_diff_relative_cst <- DT::renderDataTable(server = T, {
      DT::datatable(step25_diff_cst()$relative_cst, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "relative_diff_cst_filter"), 
                                     list(extend = 'excel', filename = "relative_diff_cst_filter"))
                    ))
    })
    
  })
  
  
  ### 25.4 draw the differential cst plot ##############################
  step25_diff_cst_plot <- reactive({
    
    req(input$act_step25_diff_draw_cst)
    
    if (is.null(step25_diff_cst())) {return(NULL)}
    
    source("R/draw_codon.R")
    
    # browser()
    
    # absolute cst
    absolute_cst_plot <- draw_codon(codon = step25_diff_cst()$absolute_cst,
                                    x = "Codon",
                                    y = input$in_step25_diff_method,
                                    xlab = "Codon",
                                    ylab = "CST",
                                    title = "Codon selection time",
                                    fill_color = "AA",
                                    color_palette = input$in_step25_diff_dot_color,
                                    color_alpha = input$in_step25_diff_alpha,
                                    
                                    dot_size = input$in_step25_diff_dot_size,
                                    line_size = input$in_step25_diff_line_size,
                                    grid_width = input$in_step25_diff_grid_width,
                                    grid_color = input$in_step25_diff_grid_color,
                                    
                                    font_size = input$in_step25_diff_font_size,
                                    font_color = 'black',
                                    legend_row = input$in_step25_diff_legend_row,
                                    
                                    wrap = input$in_step25_diff_wrap,
                                    wrap_group = "AA",
                                    wrap_row = input$in_step25_diff_wrap_row,
                                    wrap_scales = input$in_step25_diff_wrap_scales)
    
    # relative cst
    relative_cst_plot <- draw_codon(codon = step25_diff_cst()$relative_cst,
                                    x = "Codon",
                                    y = input$in_step25_diff_method,
                                    xlab = "Codon",
                                    ylab = "Relative CST",
                                    title = "Codon selection time",
                                    fill_color = "AA",,
                                    color_palette = input$in_step25_diff_dot_color,
                                    color_alpha = input$in_step25_diff_alpha,
                                    
                                    dot_size = input$in_step25_diff_dot_size,
                                    line_size = input$in_step25_diff_line_size,
                                    grid_width = input$in_step25_diff_grid_width,
                                    grid_color = input$in_step25_diff_grid_color,
                                    
                                    font_size = input$in_step25_diff_font_size,
                                    font_color = 'black',
                                    legend_row = input$in_step25_diff_legend_row,
                                    
                                    wrap = input$in_step25_diff_wrap,
                                    wrap_group = "AA",
                                    wrap_row = input$in_step25_diff_wrap_row,
                                    wrap_scales = input$in_step25_diff_wrap_scales)
    
    return(list(absolute_cst_plot = absolute_cst_plot,
                relative_cst_plot = relative_cst_plot))
  })
  
  ## output the figure
  observeEvent(input$act_step25_diff_draw_cst, {
    
    ## output the figure
    output$out_step25_diff_abs_cst_plot <- renderPlot(
      width = input$out_step25_diff_fig_width * 100,
      height = input$out_step25_diff_fig_height * 100,
      {step25_diff_cst_plot()$absolute_cst_plot})
    
    output$out_step25_diff_rel_cst_plot <- renderPlot(
      width = input$out_step25_diff_fig_width * 100,
      height = input$out_step25_diff_fig_height * 100,
      {step25_diff_cst_plot()$relative_cst_plot})
    
    output$save_step25_diff_abs_cst_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step25_diff_out_name, "-diff-absolute-cst-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step25_diff_cst_plot()$absolute_cst_plot, filename = file,
               width = input$out_step25_diff_fig_width, height = input$out_step25_diff_fig_height)
      }
    )
    
    output$save_step25_diff_rel_cst_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step25_diff_out_name, "-diff-relative-cst-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step25_diff_cst_plot()$relative_cst_plot, filename = file,
               width = input$out_step25_diff_fig_width, height = input$out_step25_diff_fig_height)
      }
    )
    
  })
  
  
  ############################################################
  ## step26 Odd Ratio ###############################
  
  ### 26.0 import the codon odd_ratio data ##############################
  step26_odd_ratio <- reactive({
    
    req(input$act_step26_import_odd_ratio)
    
    if (is.null(input$in_step26_odd_ratio)) {return(NULL)}
    
    # import the raw data table
    odd_ratio <- read.table(file = input$in_step26_odd_ratio$datapath, 
                            row.names = NULL, header = T, check.names = F, sep = '\t')
    
    # check the file format
    
    # browser()
    
    if (input$in_step26_class == "number") {
      # retrieve the odd ratio number
      odd_ratio <- odd_ratio %>% 
        dplyr::select(c('Codon', 'AA', 'Abbr', 'codon_sum', 'control_codon_sum', 'treat_codon_sum', 
                        'control_mean_odd_ratio', 'treat_mean_odd_ratio',
                        'control_number', 'treat_number', 'control_number_relative', 
                        'treat_number_relative', 'control_group', 'treat_group')) %>% 
        dplyr::arrange(control_group, treat_group, Abbr)
      
      ck_tr_group <- odd_ratio %>% 
        dplyr::select(c('Codon', 'AA', 'Abbr', 'codon_sum', 'control_group', 'treat_group')) %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr', 'codon_sum'), 
                            names_to = "Group1", values_to = "Sample")
      
      ck_tr_rel <- odd_ratio %>% 
        dplyr::select(c('Codon', 'AA', 'Abbr'), contains("number") & contains("relative")) %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), 
                            names_to = "Group2", values_to = "Odd_Ratio")
      
      ck_tr_abs <- odd_ratio %>% 
        dplyr::select(c('Codon', 'AA', 'Abbr'), contains("number") & !contains("relative")) %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), 
                            names_to = "Group2", values_to = "Odd_Ratio")
      
    } else if (input$in_step26_class == "proportion") {
      # retrieve the odd ratio proportion
      odd_ratio <- odd_ratio %>% 
        dplyr::select(c('Codon', 'AA', 'Abbr', 'codon_sum', 'control_codon_sum', 'treat_codon_sum', 
                        'control_mean_odd_ratio', 'treat_mean_odd_ratio',
                        'control_proportion', 'treat_proportion', 'control_proportion_relative', 
                        'treat_proportion_relative', 'control_group', 'treat_group')) %>% 
        dplyr::arrange(control_group, treat_group, Abbr)
      
      ck_tr_group <- odd_ratio %>% 
        dplyr::select(c('Codon', 'AA', 'Abbr', 'codon_sum', 'control_group', 'treat_group')) %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr', 'codon_sum'), 
                            names_to = "Group1", values_to = "Sample")
      
      ck_tr_rel <- odd_ratio %>% 
        dplyr::select(c('Codon', 'AA', 'Abbr'), contains("proportion") & contains("relative")) %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), 
                            names_to = "Group2", values_to = "Odd_Ratio")
      
      ck_tr_abs <- odd_ratio %>% 
        dplyr::select(c('Codon', 'AA', 'Abbr'), contains("proportion") & !contains("relative")) %>% 
        tidyr::pivot_longer(cols = -c('Codon', 'AA', 'Abbr'), 
                            names_to = "Group2", values_to = "Odd_Ratio")
      
    }
    
    absolute_odd_ratio <- cbind(ck_tr_group, ck_tr_abs[, c('Group2', 'Odd_Ratio')])
    relative_odd_ratio <- cbind(ck_tr_group, ck_tr_rel[, c('Group2', 'Odd_Ratio')])
    
    return(list(absolute_odd_ratio = absolute_odd_ratio,
                relative_odd_ratio = relative_odd_ratio))
  })
  
  ## output the table
  observeEvent(input$act_step26_import_odd_ratio, {
    
    ## output the table
    output$out_step26_absolute_odd_ratio <- DT::renderDataTable(server = T, {
      DT::datatable(step26_odd_ratio()$absolute_odd_ratio, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "odd_ratio_filter"), 
                                     list(extend = 'excel', filename = "odd_ratio_filter"))
                    ))
    })
    
    output$out_step26_relative_odd_ratio <- DT::renderDataTable(server = T, {
      DT::datatable(step26_odd_ratio()$relative_odd_ratio, rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "relative_odd_ratio_filter"), 
                                     list(extend = 'excel', filename = "relative_odd_ratio_filter"))
                    ))
    })
    
  })
  
  
  ### 26.1 set the design ##############################
  observeEvent(input$act_step26_import_odd_ratio, {
    
    observe({
      req(input$act_step26_import_odd_ratio)
      
      if (is.null(step26_odd_ratio())) {return(NULL)}
      
      uniq_group <- step26_odd_ratio()$absolute_odd_ratio %>% 
        dplyr::select(c("Group1", "Sample")) %>% 
        dplyr::distinct()
      
      ck_group <- uniq_group %>% 
        dplyr::filter(Group1 == 'control_group')
      tr_group <- uniq_group %>% 
        dplyr::filter(Group1 == 'treat_group')
      
      updateSelectizeInput(session, "in_step26_control", choices = ck_group$Sample, selected = ck_group$Sample[1], server = TRUE)
      updateSelectizeInput(session, "in_step26_treat", choices = tr_group$Sample, selected = tr_group$Sample[1], server = TRUE)
      
    })
    
  })
  
  
  ### 26.2 draw the odd_ratio plot ##############################
  step26_odd_ratio_plot <- reactive({
    
    req(input$act_step26_draw_odd_ratio)
    
    if (is.null(step26_odd_ratio())) {return(NULL)}
    
    # browser()
    source("R/draw_codon.R")
    
    absolute_odd_ratio <- step26_odd_ratio()$absolute_odd_ratio %>% 
      dplyr::filter(Sample %in% c(input$in_step26_control, input$in_step26_treat))
    
    relative_odd_ratio <- step26_odd_ratio()$relative_odd_ratio %>% 
      dplyr::filter(Sample %in% c(input$in_step26_control, input$in_step26_treat))
    
    # absolute odd_ratio score
    absolute_odd_ratio_plot <- draw_codon(codon = absolute_odd_ratio,
                                          x = "Codon",
                                          y = "Odd_Ratio",
                                          xlab = "Codon",
                                          ylab = "Odd Ratio",
                                          title = "Codon Odd Ratio",
                                          fill_color = "Sample",
                                          dot_size = input$in_step26_dot_size,
                                          color_palette = input$in_step26_dot_color,
                                          color_alpha = input$in_step26_alpha,
                                          
                                          line_size = input$in_step26_line_size,
                                          grid_width = input$in_step26_grid_width,
                                          grid_color = input$in_step26_grid_color,
                                          
                                          font_size = input$in_step26_font_size,
                                          font_color = 'black',
                                          legend_row = input$in_step26_legend_row,
                                          
                                          wrap = input$in_step26_wrap,
                                          wrap_group = "AA",
                                          wrap_row = input$in_step26_wrap_row,
                                          wrap_scales = input$in_step26_wrap_scales)
    
    # relative odd_ratio score
    relative_odd_ratio_plot <- draw_codon(codon = relative_odd_ratio,
                                          x = "Codon",
                                          y = "Odd_Ratio",
                                          xlab = "Codon",
                                          ylab = "Odd Ratio",
                                          title = "Codon Odd Ratio",
                                          fill_color = "Sample",
                                          dot_size = input$in_step26_dot_size,
                                          color_palette = input$in_step26_dot_color,
                                          color_alpha = input$in_step26_alpha,
                                          
                                          line_size = input$in_step26_line_size,
                                          grid_width = input$in_step26_grid_width,
                                          grid_color = input$in_step26_grid_color,
                                          
                                          font_size = input$in_step26_font_size,
                                          font_color = 'black',
                                          legend_row = input$in_step26_legend_row,
                                          
                                          wrap = input$in_step26_wrap,
                                          wrap_group = "AA",
                                          wrap_row = input$in_step26_wrap_row,
                                          wrap_scales = input$in_step26_wrap_scales)
    
    return(list(absolute_odd_ratio_plot = absolute_odd_ratio_plot,
                relative_odd_ratio_plot = relative_odd_ratio_plot))
  })
  
  ## output the figure
  observeEvent(input$act_step26_draw_odd_ratio, {
    
    ## output the figure
    output$out_step26_absolute_odd_ratio_plot <- renderPlot(
      width = input$out_step26_fig_width * 100,
      height = input$out_step26_fig_height * 100,
      {step26_odd_ratio_plot()$absolute_odd_ratio_plot})
    
    output$out_step26_relative_odd_ratio_plot <- renderPlot(
      width = input$out_step26_fig_width * 100,
      height = input$out_step26_fig_height * 100,
      {step26_odd_ratio_plot()$relative_odd_ratio_plot})
    
    output$save_step26_absolute_odd_ratio_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step26_out_name, "-odd_ratio-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step26_odd_ratio_plot()$absolute_odd_ratio_plot, filename = file,
               width = input$out_step26_fig_width, height = input$out_step26_fig_height)
      }
    )
    
    output$save_step26_relative_odd_ratio_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step26_out_name, "-relative-odd_ratio-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step26_odd_ratio_plot()$relative_odd_ratio_plot, filename = file,
               width = input$out_step26_fig_width, height = input$out_step26_fig_height)
      }
    )
    
  })
  
  
  ############################################################
  ## step27 CoV ###############################
  
  ### 27.0 set the design ##############################
  observeEvent(input$act_step1_import_design, {
    col_name <- names(step1_design())
    updateSelectizeInput(session, "in_step27_cov_column", choices = col_name, selected = col_name[3], server = TRUE)
    
    # browser()
    observe({
      req(input$in_step27_cov_column)
      step1_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      select_choices <- step1_design[[input$in_step27_cov_column]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step27_cov_groups", choices = select_choices, selected = select_choices, server = TRUE)
      
    })
    
    
    updateSelectizeInput(session, "in_step27_cov_cdf_column", choices = col_name, selected = col_name[3], server = TRUE)
    
    observe({
      req(input$in_step27_cov_cdf_column)
      step1_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      select_choices <- step1_design[[input$in_step27_cov_cdf_column]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step27_cov_cdf_groups", choices = select_choices, selected = select_choices, server = TRUE)
      
    })
    
  })
  
  ### 27.1 import the cov data ##############################
  step27_cov_table <- reactive({
    
    req(input$act_step27_import_cov)
    
    if (is.null(input$in_step27_cov)) {return(NULL)}
    
    # import the raw data table
    cov_table <- fread(file = input$in_step27_cov$datapath, header = T, check.names = F, sep = '\t') %>% 
      na.omit()
    
    return(cov_table)
    
  })
  
  ## output the table
  observeEvent(input$act_step27_import_cov, {
    
    ## output the table
    output$out_step27_cov <- DT::renderDataTable(server = T, {
      DT::datatable(step27_cov_table(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "cov_table"), 
                                     list(extend = 'excel', filename = "cov_table"))
                    ))
    })
    
  })
  
  
  ### 27.2 filter the cov data ##############################
  step27_cov_flt_table <- reactive({
    
    req(input$act_step27_cov_filter)
    
    if (is.null(step27_cov_table())) {return(NULL)}
    
    # browser()
    
    # filter the speicific samples
    if (import_design_clicked()) {
      
      flt_step1_design <- step1_design() %>% 
        dplyr::select(Sample, !!sym(input$in_step27_cov_column))
      
      selected_samples <- flt_step1_design %>% 
        dplyr::filter(!!sym(input$in_step27_cov_column) %in% input$in_step27_cov_groups) %>% 
        dplyr::pull(Sample)
      
      cov_table <- step27_cov_table() %>% 
        dplyr::select(name, contains(selected_samples)) %>% 
        filter(if_any(contains("_sum"), ~is.numeric(.) & . >= input$in_step27_cov_sum)) %>% 
        filter(if_any(contains("_mean"), ~is.numeric(.) & . >= input$in_step27_cov_mean))
      
      if (!is.null(input$in_step27_cov_column) & !is.null(input$in_step27_cov_groups)) {
        
        flt_sum_table <- cov_table %>%
          dplyr::select(name, contains("_sum")) %>%
          dplyr::rename_with(~str_remove(., "_sum"), contains("_sum")) %>%
          tidyr::pivot_longer(cols = -name, names_to = "Sample", values_to = "Sum") %>%
          dplyr::left_join(flt_step1_design, by = c('Sample' = 'Sample')) %>%
          dplyr::group_by(name, !!sym(input$in_step27_cov_column)) %>%
          dplyr::reframe(Sum = mean(Sum)) %>% 
          dplyr::mutate(log2Sum = log2(Sum))
        
        flt_mean_table <- cov_table %>% 
          dplyr::select(name, contains("_mean")) %>%
          dplyr::rename_with(~str_remove(., "_mean"), contains("_mean")) %>% 
          tidyr::pivot_longer(cols = -name, names_to = "Sample", values_to = "Mean") %>%
          dplyr::left_join(flt_step1_design, by = c('Sample' = 'Sample')) %>% 
          dplyr::group_by(name, !!sym(input$in_step27_cov_column)) %>% 
          dplyr::reframe(Mean = mean(Mean)) %>% 
          dplyr::mutate(log2Mean = log2(Mean))
        
        flt_cov_table <- cov_table %>% 
          dplyr::select(name, contains("_CoV")) %>%
          dplyr::rename_with(~str_remove(., "_CoV"), contains("_CoV")) %>% 
          tidyr::pivot_longer(cols = -name, names_to = "Sample", values_to = "CoV") %>%
          dplyr::left_join(flt_step1_design, by = c('Sample' = 'Sample')) %>% 
          dplyr::group_by(name, !!sym(input$in_step27_cov_column)) %>% 
          dplyr::reframe(CoV = mean(CoV)) %>% 
          dplyr::mutate(log2CoV = log2(CoV))
        
        melt_cov_table <- flt_sum_table %>%
          dplyr::left_join(flt_mean_table, by = c('name', input$in_step27_cov_column)) %>% 
          dplyr::left_join(flt_cov_table, by = c('name', input$in_step27_cov_column))
        
        return(melt_cov_table)
        
      } else {
        return(NULL)
        
      }
      
    } else {
      
      cov_table <- step27_cov_table() %>% 
        filter(if_any(contains("_sum"), ~is.numeric(.) & . >= input$in_step27_cov_sum)) %>% 
        filter(if_any(contains("_mean"), ~is.numeric(.) & . >= input$in_step27_cov_mean))
      
      flt_sum_table <- cov_table %>% 
        dplyr::select(name, contains("_sum")) %>%
        dplyr::rename_with(~str_remove(., "_sum"), contains("_sum")) %>% 
        tidyr::pivot_longer(cols = -name, names_to = "Sample", values_to = "Sum") %>%
        dplyr::group_by(name, Sample) %>% 
        dplyr::reframe(Sum = mean(Sum)) %>% 
        dplyr::mutate(log2Sum = log2(Sum))
      
      flt_mean_table <- cov_table %>% 
        dplyr::select(name, contains("_mean")) %>%
        dplyr::rename_with(~str_remove(., "_mean"), contains("_mean")) %>% 
        tidyr::pivot_longer(cols = -name, names_to = "Sample", values_to = "Mean") %>%
        dplyr::group_by(name, Sample) %>% 
        dplyr::reframe(Mean = mean(Mean)) %>% 
        dplyr::mutate(log2Mean = log2(Mean))
      
      flt_cov_table <- cov_table %>% 
        dplyr::select(name, contains("_CoV")) %>%
        dplyr::rename_with(~str_remove(., "_CoV"), contains("_CoV")) %>% 
        tidyr::pivot_longer(cols = -name, names_to = "Sample", values_to = "CoV") %>%
        dplyr::group_by(name, Sample) %>% 
        dplyr::reframe(CoV = mean(CoV)) %>% 
        dplyr::mutate(log2CoV = log2(CoV))
      
      melt_cov_table <- flt_sum_table %>%
        dplyr::left_join(flt_mean_table, by = c('name', 'Sample')) %>% 
        dplyr::left_join(flt_cov_table, by = c('name', 'Sample')) %>% 
        na.omit()
      
      return(melt_cov_table)
    }
    
  })
  
  ## output the table
  observeEvent(input$act_step27_cov_filter, {
    
    ## output the table
    output$out_step27_cov_filtered <- DT::renderDataTable(server = T, {
      DT::datatable(step27_cov_flt_table(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "cov_table"), 
                                     list(extend = 'excel', filename = "cov_table"))
                    ))
    })
    
  })
  
  ### 27.3 fit the CoV ##############################
  step27_cov_fit_table <- reactive({
    
    req(input$act_step27_cov_fitted)
    
    if (is.null(step27_cov_flt_table())) {return(NULL)}
    
    # browser()
    
    # use the nchar to replace the is.null, is.null and is.na can't be used here
    if (nchar(input$in_step27_cov_column) == 0){
      group_column <- "Sample"
    } else {
      group_column <- input$in_step27_cov_column
    }
    
    # browser()
    now_groups <- step27_cov_flt_table() %>% 
      dplyr::select(!!sym(group_column)) %>% 
      dplyr::distinct() %>% 
      dplyr::pull(!!sym(group_column))
    
    # fit all the groups and return the results to dataframe
    # set the CoV model
    nonlinear_model <- function(mu, alpha, beta) {
      return(0.5 * log2((beta / mu) + alpha))
    }
    
    now_fitted_results <- list()
    
    for (i in 1:length(now_groups)) {
      Mean <- step27_cov_flt_table() %>% 
        dplyr::filter(!!sym(group_column) == now_groups[i]) %>% 
        dplyr::pull(Mean)
      CV <- step27_cov_flt_table() %>% 
        dplyr::filter(!!sym(group_column) == now_groups[i]) %>% 
        dplyr::pull(CoV)
      
      cov_fit <- nlsLM(log2(CV) ~ nonlinear_model(Mean, alpha, beta), start = list(alpha = 1, beta = 1))
      
      cov_fit_summary <- summary(cov_fit)$coefficients
      fitted_alpha <- cov_fit_summary[1, 1]
      fitted_beta <- cov_fit_summary[2, 1]
      
      # browser()
      cov_fit_table <- cov_fit_summary %>% 
        as.data.frame() %>% 
        dplyr::mutate(Group = now_groups[i]) %>% 
        dplyr::mutate(Formula = paste0("log2(CoV) = 0.5 * log2[(", round(fitted_beta, 4), "/mean) + ", round(fitted_alpha, 4), "]")) %>% 
        tibble::rownames_to_column(var = "Parameters")
      
      now_fitted_results[[i]] <- cov_fit_table
    }
    
    fitted_results <- do.call(rbind, now_fitted_results)
    
    return(fitted_results)
  })
  
  ## output the table
  observeEvent(input$act_step27_cov_fitted, {
    
    ## output the table
    output$out_step27_cov_fitted_results <- DT::renderDataTable(server = T, {
      DT::datatable(step27_cov_fit_table(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "cov_fitted_table"), 
                                     list(extend = 'excel', filename = "cov_fitted_table"))
                    ))
    })
    
  })
  
  ### 27.4 draw the scatter cov data ##############################
  step27_cov_scatter_plot <- reactive({
    
    req(input$act_step27_cov_draw)
    
    if (is.null(step27_cov_flt_table())) {return(NULL)}
    
    # browser()
    source("R/draw_CoV.R")
    
    if (nchar(input$in_step27_cov_column) == 0){
      group_color <- "Sample"
    } else {
      group_color <- input$in_step27_cov_column
    }
    
    # browser()
    # cov scatter plot
    cov_fit_plot <- draw_cov_line(cov_table = step27_cov_flt_table(),
                                  x = 'log2Mean', y = 'log2CoV', group = group_color,
                                  
                                  facet = input$in_step27_cov_facet, 
                                  wrap_num = input$in_step27_cov_wrap_num,
                                  
                                  dot_color = input$in_step27_cov_dot_color,
                                  dot_alpha = input$in_step27_cov_dot_alpha,
                                  dot_size = input$in_step27_cov_dot_size,
                                  
                                  line_color = input$in_step27_cov_line_color,
                                  line_alpha = input$in_step27_cov_line_alpha,
                                  line_width = input$in_step27_cov_line_width,
                                  
                                  font_size = input$in_step27_cov_font_size,
                                  xstart = input$in_step27_cov_x_min, xend = input$in_step27_cov_x_max,
                                  ystart = input$in_step27_cov_y_min, yend = input$in_step27_cov_y_max,
                                  xbreak = input$in_step27_cov_x_breaks, ybreak = input$in_step27_cov_y_breaks,
                                  
                                  type = input$in_step27_cov_figure_type,
                                  
                                  xlabel = input$in_step27_cov_xlabel,
                                  ylabel = input$in_step27_cov_ylabel,
                                  title = input$in_step27_cov_title
    )
    
    return(cov_fit_plot)
  })
  
  ## output the figure
  observeEvent(input$act_step27_cov_draw, {
    
    ## output the figure
    output$out_step27_cov_scatter_plot <- renderPlot(
      width = input$out_step27_cov_fig_width * 100,
      height = input$out_step27_cov_fig_height * 100,
      {step27_cov_scatter_plot()})
    
    output$save_step27_cov_scatter_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step27_cov_out_name, "-CoV-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step27_cov_scatter_plot(), filename = file,
               width = input$out_step27_cov_fig_width, height = input$out_step27_cov_fig_height)
      }
    )
    
  })
  
  
  ### 27.5 import the cov cdf data ##############################
  step27_cov_cdf_table <- reactive({
    
    req(input$act_step27_cov_cdf_import)
    
    if (is.null(input$in_step27_cov_cdf)) {return(NULL)}
    
    # import the raw data table
    cov_cdf_table <- fread(file = input$in_step27_cov_cdf$datapath, header = T, check.names = F, sep = '\t') %>% 
      na.omit()
    
    
    # browser()
    
    # filter the speicific samples
    uniq_gene_names <- unique(cov_cdf_table$name)
    updateSelectizeInput(session, "in_step27_cov_cdf_gene", choices = uniq_gene_names, selected = uniq_gene_names[1], server = TRUE)
    
    return(cov_cdf_table)
    
  })
  
  ## output the table
  observeEvent(input$act_step27_cov_cdf_import, {
    
    ## output the table
    output$out_step27_cov_cdf <- DT::renderDataTable(server = T, {
      DT::datatable(step27_cov_cdf_table(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "cov_cdf_table"), 
                                     list(extend = 'excel', filename = "cov_cdf_table"))
                    ))
    })
    
  })
  
  
  ### 27.6 filter the CoV CDF data ##############################
  step27_cov_cdf_flt_table <- reactive({
    
    req(input$act_step27_cov_cdf_filter)
    
    if (is.null(step27_cov_cdf_table())) {return(NULL)}

    # browser()
    
    # filter the raw data table
    cov_cdf_table <- step27_cov_cdf_table() %>%
      dplyr::filter(name == input$in_step27_cov_cdf_gene)
    
    # calculate the average of the CoV CDF data
    if (nchar(input$in_step27_cov_cdf_column) != 0) {

      ribo_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      
      for (i in 1:length(input$in_step27_cov_cdf_groups)) {
        now_group <- input$in_step27_cov_cdf_groups[i]
        now_name <- ribo_design %>%
          dplyr::filter(!!sym(input$in_step27_cov_cdf_column) == now_group) %>% 
          dplyr::pull("Sample")
        
        cov_cdf_table <- cov_cdf_table %>%
          dplyr::mutate(!!now_group := rowMeans(dplyr::select(., all_of(now_name))))
      }
      
      cov_cdf_table <- cov_cdf_table %>%
        dplyr::select(any_of(c("name", "now_nt", "from_tis", "from_tts", "region", "codon", "frame", 
                               input$in_step27_cov_cdf_groups))) %>%
        tidyr::pivot_longer(cols = input$in_step27_cov_cdf_groups, 
                            names_to = "Sample", values_to = "CoV")
      
      return(cov_cdf_table)
      
    } else {
      
      cov_cdf_table <- cov_cdf_table %>%
        tidyr::pivot_longer(cols = -c("name", "now_nt", "from_tis", "from_tts", "region", "codon", "frame"), 
                            names_to = "Sample", values_to = "CoV")
      
      return(cov_cdf_table)
    }
    
  })
  
  ## output the table
  observeEvent(input$act_step27_cov_cdf_filter, {
    
    ## output the table
    output$out_step27_cov_cdf_filtered <- DT::renderDataTable(server = T, {
      DT::datatable(step27_cov_cdf_flt_table(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "cov_cdf_flt_table"), 
                                     list(extend = 'excel', filename = "cov_cdf_flt_table"))
                    ))
    })
    
  })
  
  
  ### 27.7 draw the CoV CDF ##############################
  step27_cov_cdf_plot <- reactive({
    
    req(input$act_step27_cov_cdf_draw)
    
    if (is.null(step27_cov_cdf_flt_table())) {return(NULL)}

    # browser()
    
    # draw the CoV CDF
    # c("name", "now_nt", "from_tis", "from_tts", "region", "codon", "frame", "Sample", "CoV"), 
    source("R/draw_Isoforms_Mapped.R")
    
    xmin <- ifelse(!is.na(input$in_step27_cov_cdf_x_min), input$in_step27_cov_cdf_x_min, NA)
    xmax <- ifelse(!is.na(input$in_step27_cov_cdf_x_max), input$in_step27_cov_cdf_x_max, NA)
    ymin <- ifelse(!is.na(input$in_step27_cov_cdf_y_min), input$in_step27_cov_cdf_y_min, NA)
    ymax <- ifelse(!is.na(input$in_step27_cov_cdf_y_max), input$in_step27_cov_cdf_y_max, NA)
    
    if (input$in_step27_cov_cdf_y_scale == 'log') {
      cov_cdf_flt_table <- step27_cov_cdf_flt_table() %>%
        dplyr::mutate(CoV = log2(CoV + 1))
    } else {
      cov_cdf_flt_table <- step27_cov_cdf_flt_table()
    }
    
    # browser()
    
    cov_cdf_plot <- draw_isoforms_mapped(gene_name = input$in_step27_cov_cdf_gene,
                                         gene_reads = cov_cdf_flt_table,
                                         
                                         x = 'now_nt',
                                         y = 'CoV',
                                         color = 'Sample',
                                         
                                         xlabel = input$in_step27_cov_cdf_xlabel,
                                         ylabel = input$in_step27_cov_cdf_ylabel,
                                         title = paste0(input$in_step27_cov_cdf_gene, " (CoV)"),
                                         
                                         gene_color = input$in_step27_cov_cdf_color,
                                         
                                         plot_type = "line",
                                         
                                         line_width = input$in_step27_cov_cdf_line_width,
                                         grid_width = input$in_step27_cov_cdf_grid_width,
                                         fill_alpha = input$in_step27_cov_cdf_alpha,
                                         
                                         xmin = xmin,
                                         xmax = xmax,
                                         xbreaks = input$in_step27_cov_cdf_x_break,
                                         
                                         ymin = ymin,
                                         ymax = ymax,
                                         ybreaks = input$in_step27_cov_cdf_y_break,
                                         
                                         sqrt = input$in_step27_cov_cdf_y_scale,
                                         
                                         facet = input$in_step27_cov_cdf_facet,
                                         show_legend = input$in_step27_cov_cdf_legend,
                                         
                                         font_size = input$in_step27_cov_cdf_font_size)
    
    return(cov_cdf_plot)
    
  })
  
  ## output the figure
  observeEvent(input$act_step27_cov_cdf_draw, {
    
    ## output the figure
    output$out_step27_cov_cdf_plot <- renderPlot(
      width = input$out_step27_cov_cdf_fig_width * 100,
      height = input$out_step27_cov_cdf_fig_height * 100,
      {step27_cov_cdf_plot()})
    
    output$save_step27_cov_cdf_plot <- downloadHandler(
      filename = function() {
        paste0(input$out_step27_cov_cdf_out_name, "-", input$in_step27_cov_cdf_gene, "-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step27_cov_cdf_plot(), filename = file,
               width = input$out_step27_cov_cdf_fig_width, height = input$out_step27_cov_cdf_fig_height)
      }
    )
    
  })
  
  
  ############################################################
  ## step28 Meta Codon ###############################
  
  ### 28.0 set the design ##############################
  observeEvent(input$act_step1_import_design, {
    col_name <- names(step1_design())
    updateSelectizeInput(session, "in_step28_meta_codon_column", choices = col_name, selected = col_name[3], server = TRUE)
    
    # browser()
    observe({
      req(input$in_step28_meta_codon_column)
      step1_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      select_choices <- step1_design[[input$in_step28_meta_codon_column]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step28_meta_codon_groups", choices = select_choices, selected = select_choices, server = TRUE)
      
    })
    
  })
  
  ### 28.1 import the meta codon data ##############################
  step28_meta_codon_dst <- reactive({
    
    req(input$act_step28_import_meta_codon)
    
    if (is.null(input$in_step28_meta_codon)) {return(NULL)}
    
    # import the raw data table
    meta_codon_dst <- read.table(file = input$in_step28_meta_codon$datapath, 
                                 row.names = NULL, header = T, check.names = F, sep = '\t')
    
    # browser()
    
    # calculate the average of the meta codon density
    if (nchar(input$in_step28_meta_codon_column) != 0 && nchar(input$in_step28_meta_codon_groups[1]) != 0) {
      
      # browser()
      
      meta_codon_dst <- meta_codon_dst %>% 
        tidyr::pivot_longer(cols = -Codon:-Frame, names_to = "Sample", values_to = "Density") %>%
        dplyr::left_join(step1_design(), by = c('Sample' = 'Sample')) %>%
        na.omit() %>% 
        dplyr::filter(!!sym(input$in_step28_meta_codon_column) %in% c(input$in_step28_meta_codon_groups)) %>%
        dplyr::group_by(Codon, Nucleotide, Frame, !!sym(input$in_step28_meta_codon_column)) %>% 
        dplyr::reframe(Density = mean(Density)) %>% 
        dplyr::ungroup()
      
    } else {
      
      meta_codon_dst <- meta_codon_dst %>% 
        tidyr::pivot_longer(cols = -Codon:-Frame, names_to = "Sample", values_to = "Density")
      
    }
    
    # scale the density
    # browser()
    
    meta_codon_dst <- meta_codon_dst %>% 
      dplyr::group_by(!!sym(input$in_step28_meta_codon_column)) %>% 
      dplyr::mutate(Scaled = Density / mean(Density)) %>% 
      dplyr::ungroup()
    
    return(meta_codon_dst)
    
  })
  
  ## output the table
  observeEvent(input$act_step28_import_meta_codon, {
    
    ## output the table
    output$out_step28_meta_codon_dst <- DT::renderDataTable(server = T, {
      DT::datatable(step28_meta_codon_dst(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "meta_codon_dst"), 
                                     list(extend = 'excel', filename = "meta_codon_dst"))
                    ))
    })
    
  })
  
  
  ### 28.2 draw the meta codon plot ##############################
  step28_meta_codon_plot <- reactive({
    
    req(input$act_step28_meta_codon_draw_plot)
    
    if (is.null(step28_meta_codon_dst())) {return(NULL)}
    
    # browser()
    source("R/draw_dot_line.R")
    if (nchar(input$in_step28_meta_codon_column) == 0){
      group_value = "Sample"
    } else{
      group_value = input$in_step28_meta_codon_column
    }
    
    meta_codon_dst <- step28_meta_codon_dst() %>% 
      dplyr::filter(!!sym(input$in_step28_meta_codon_site) >= input$in_step28_meta_codon_x_min) %>% 
      dplyr::filter(!!sym(input$in_step28_meta_codon_site) <= input$in_step28_meta_codon_x_max)
    
    # browser()
    # meta codon line plot
    line_plot <- draw_dot_line(distr = step28_meta_codon_dst(), 
                               x = input$in_step28_meta_codon_site, y = input$in_step28_meta_codon_scaled, 
                               titles = c("Meta Codon Density"), group = group_value, 
                               facet = input$in_step28_meta_codon_wrap, wrap_num = input$in_step28_meta_codon_wrap_num,
                               dot_size = input$in_step28_meta_codon_dot_size, font_size = input$in_step28_meta_codon_font_size,
                               line_color = input$in_step28_meta_codon_line_color, line_width = input$in_step28_meta_codon_line_size,
                               fill_alpha = input$in_step28_meta_codon_alpha, 
                               xstart = input$in_step28_meta_codon_x_min, xend = input$in_step28_meta_codon_x_max,
                               ystart = input$in_step28_meta_codon_y_min, yend = input$in_step28_meta_codon_y_max,
                               xbreak = input$in_step28_meta_codon_x_breaks, ybreak = input$in_step28_meta_codon_y_breaks,
                               grid_width = input$in_step28_meta_codon_grid_width,
                               grid_color = input$in_step28_meta_codon_grid_color,
                               xpeak = input$in_step28_meta_codon_label,
                               type = "line")
    
    # meta codon heatmap
    heat_plot <- draw_dot_line(distr = step28_meta_codon_dst(), 
                               x = input$in_step28_meta_codon_site, y = input$in_step28_meta_codon_scaled, 
                               titles = c("Meta Codon Density"), group = group_value, facet = FALSE,
                               font_size = input$in_step28_meta_codon_font_size,
                               fill_color = input$in_step28_meta_codon_heat_color, fill_alpha = input$in_step28_meta_codon_alpha, 
                               xstart = input$in_step28_meta_codon_x_min, xend = input$in_step28_meta_codon_x_max,
                               xbreak = input$in_step28_meta_codon_x_breaks,
                               type = "heat")
    
    return(list(line_plot = line_plot,
                heat_plot = heat_plot))
  })
  
  ## output the figure
  observeEvent(input$act_step28_meta_codon_draw_plot, {
    
    ## output the figure
    output$out_step28_meta_codon_line_plot <- renderPlot(
      width = input$out_step28_meta_codon_fig_width * 100,
      height = input$out_step28_meta_codon_fig_height * 100,
      {step28_meta_codon_plot()$line_plot})
    
    output$out_step28_meta_codon_heat_map <- renderPlot(
      width = input$out_step28_meta_codon_fig_width * 100,
      height = input$out_step28_meta_codon_fig_height * 100,
      {step28_meta_codon_plot()$heat_plot})
    
    output$save_step28_meta_codon_line <- downloadHandler(
      filename = function() {
        paste0(input$out_step28_meta_codon_out_name, "-line-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step28_meta_codon_plot()$line_plot, filename = file,
               width = input$out_step28_meta_codon_fig_width, height = input$out_step28_meta_codon_fig_height)
      }
    )
    
    output$save_step28_meta_codon_heat <- downloadHandler(
      filename = function() {
        paste0(input$out_step28_meta_codon_out_name, "-heat-map-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step28_meta_codon_plot()$heat_plot, filename = file,
               width = input$out_step28_meta_codon_fig_width, height = input$out_step28_meta_codon_fig_height)
      }
    )
    
  })
  
  
  ### 28.3 import the meta codon sequence ##############################
  step28_meta_codon_sequence <- reactive({
    
    req(input$act_step28_import_meta_seq)
    
    if (is.null(input$in_step28_meta_seq)) {return(NULL)}
    
    # browser()
    # import the raw data table
    meta_codon_seq_table <- fread(file = input$in_step28_meta_seq$datapath, 
                                  header = T, check.names = F, sep = '\t') %>% 
      as_tibble()
    
    # translate the codon to AA with biostrings
    genetic_code <- unlist(str_split(input$in_step28_meta_seq_code, " - "))[1]
    
    translate_codon <- function(codon_seq, genetic_code) {
      protein_seq <- translate(DNAStringSet(codon_seq), 
                               no.init.codon = T,
                               genetic.code = getGeneticCode(genetic_code)) %>% as.character()
      return(protein_seq)
    }
    
    # browser()
    
    # convert the sequence to the desired type
    if (input$in_step28_meta_seq_seqtype == "dna") {
      meta_codon_seq_table <- meta_codon_seq_table
      
    } else if (input$in_step28_meta_seq_seqtype == "rna") {
      meta_codon_seq_table <- meta_codon_seq_table %>% 
        dplyr::mutate_at(vars(!contains(c('Position', 'name'))), ~str_replace_all(.x, "T", "U"))
      
    } else if (input$in_step28_meta_seq_seqtype == "aa") {
      meta_codon_seq_table <- meta_codon_seq_table %>% 
        dplyr::mutate_at(vars(!contains(c('Position', 'name'))), ~translate_codon(.x, genetic_code))
    }
    
    return(meta_codon_seq_table)
    
  })
  
  ## output the table
  observeEvent(input$act_step28_import_meta_seq, {
    
    ## output the table
    output$out_step28_meta_seq <- DT::renderDataTable(server = T, {
      DT::datatable(step28_meta_codon_sequence(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "meta_codon_dst"), 
                                     list(extend = 'excel', filename = "meta_codon_dst"))
                    ))
    })
    
  })
  
  
  ### 28.4 draw the seqlogo ##############################
  step28_meta_codon_seqlogo <- reactive({
    req(input$act_step28_meta_seq_draw_plot)
    
    if (is.null(step28_meta_codon_sequence())) {return(NULL)}
    
    # browser()
    source("R/draw_seqlogo.R")
    
    # filter the range of the sequence
    meta_codon_sequence <- step28_meta_codon_sequence() %>% 
      dplyr::select(-Codon, -name)
    
    column_names <- colnames(meta_codon_sequence) %>% 
      as.numeric()
    
    min_seq_len <- min(column_names)
    max_seq_len <- max(column_names)
    
    xstart <- ifelse(is.na(input$in_step28_meta_seq_x_min), min_seq_len, as.numeric(input$in_step28_meta_seq_x_min))
    xend <- ifelse(is.na(input$in_step28_meta_seq_x_max), max_seq_len, as.numeric(input$in_step28_meta_seq_x_max))
    xrange <- paste(c(xstart:xend))
    
    merged_sequence <- meta_codon_sequence %>% 
      dplyr::select(any_of(xrange)) %>%
      dplyr::mutate(sequence = do.call(paste0, c(., sep = ""))) %>%
      dplyr::pull(sequence)
    
    # convert the sequence to the biostrings object
    if (input$in_step28_meta_seq_seqtype == "dna") {
      seq_freq <- 0.25
      merged_sequence_set <- Biostrings::DNAStringSet(merged_sequence)
      xstart <- xstart * 3
      xend <- xend * 3 + 2
      meta_codon_sequence_pfm <- consensusMatrix(merged_sequence_set, baseOnly = TRUE) %>% 
        subset(., rowSums(.) > 0) %>% 
        as.data.frame() %>% 
        dplyr::rename_all(~paste(seq(xstart, xend, 1))) %>% 
        as.matrix()
      
    } else if (input$in_step28_meta_seq_seqtype == "rna") {
      seq_freq <- 0.25
      merged_sequence_set <- Biostrings::RNAStringSet(merged_sequence)
      xstart <- xstart * 3
      xend <- xend * 3 + 2
      meta_codon_sequence_pfm <- consensusMatrix(merged_sequence_set, baseOnly = TRUE) %>% 
        subset(., rowSums(.) > 0) %>% 
        as.data.frame() %>% 
        dplyr::rename_all(~paste(seq(xstart, xend, 1))) %>% 
        as.matrix()
      
    } else if (input$in_step28_meta_seq_seqtype == "aa") {
      seq_freq <- 0.05
      merged_sequence_set <- AAStringSet(merged_sequence)
      meta_codon_sequence_pfm <- consensusMatrix(merged_sequence_set, baseOnly = TRUE) %>% 
        subset(., rowSums(.) > 0) %>% 
        as.data.frame() %>% 
        dplyr::rename_all(~paste(seq(xstart, xend, 1))) %>% 
        as.matrix()
    }
    
    # format the pfm/pwm/ppm
    if (input$in_step28_meta_seq_method == "custom") {
      if (input$in_step28_meta_seq_custom == "raw") {
        meta_codon_sequence_mat <- meta_codon_sequence_pfm
        meta_codon_sequence_ppm <- meta_codon_sequence_pfm / colSums(meta_codon_sequence_pfm)
        meta_codon_sequence_pwm <- log2(meta_codon_sequence_ppm / seq_freq)
        
      } else if (input$in_step28_meta_seq_custom == "foldchange"){
        meta_codon_sequence_mat <- log2(meta_codon_sequence_pfm / rowMeans(meta_codon_sequence_pfm))
        
      } else if (input$in_step28_meta_seq_custom == "delta") {
        meta_codon_sequence_mat <- meta_codon_sequence_pfm - rowMeans(meta_codon_sequence_pfm)
        
      } else if (input$in_step28_meta_seq_custom == "enrichment") {
        meta_codon_sequence_mat <- (meta_codon_sequence_pfm - rowMeans(meta_codon_sequence_pfm)) / rowMeans(meta_codon_sequence_pfm)
      } 
    } else {
      meta_codon_sequence_mat <- meta_codon_sequence_pfm
    }
    
    # browser()
    meta_codon_sequence_plot <- draw_seqlogo(pwm = meta_codon_sequence_mat,
                                             method = input$in_step28_meta_seq_method, seq_type = input$in_step28_meta_seq_seqtype,
                                             col_scheme = input$in_step28_meta_seq_fill_color, col_alpha = input$in_step28_meta_seq_fill_alpha,
                                             stack_width = input$in_step28_meta_seq_font_stack, 
                                             font_family = input$in_step28_meta_seq_font_family, font_size = input$in_step28_meta_seq_font_size,
                                             xstart = xstart, xend = xend,
                                             xlabel = "site", ylabel = input$in_step28_meta_seq_custom, title = "meta-codon seqlogo")
    
    return(meta_codon_sequence_plot)
    
  })
  
  ## output the alignment plot
  observeEvent(input$act_step28_meta_seq_draw_plot, {
    # browser()
    ## show the alignment plot
    output$out_step28_meta_seq_logo <- renderPlot(
      width = input$out_step28_meta_seq_fig_width * 100,
      height = input$out_step28_meta_seq_fig_height * 100,
      {step28_meta_codon_seqlogo()})
    
    # save the alignment plot
    output$save_step28_meta_seq_logo <- downloadHandler(
      filename = function() {
        paste0(input$out_step28_meta_seq_out_name, "-meta-codon-seqlogo-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step28_meta_codon_seqlogo(), filename = file,
               width = input$out_step28_meta_seq_fig_width, height = input$out_step28_meta_seq_fig_height)
      }
    )
    
  })
  
  
  ############################################################
  ## step29 SeRP enrichment ###############################
  
  ### 29.0 set the design ##############################
  observeEvent(input$act_step1_import_design, {
    col_name <- names(step1_design())
    updateSelectizeInput(session, "in_step28_meta_codon_column", choices = col_name, selected = col_name[3], server = TRUE)
    
    # browser()
    observe({
      req(input$in_step28_meta_codon_column)
      step1_design <- step1_design() %>% dplyr::filter(SeqType == "RIBO")
      select_choices <- step1_design[[input$in_step28_meta_codon_column]]
      select_choices <- unique(select_choices)
      updateSelectizeInput(session, "in_step28_meta_codon_groups", choices = select_choices, selected = select_choices, server = TRUE)
      
    })
    
  })
  
  ### 29.1 import the Ribo-seq metaplot table  #################
  step29_ribo_meta <- reactive({
    
    req(input$act_step29_import_serp_enrich_meta)
    
    if (is.null(input$in_step29_serp_enrich_meta)) {return(NULL)}
    
    ribo_meta_table <- read.table(input$in_step29_serp_enrich_meta$datapath, 
                                  sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    
    select_choices <- unique(ribo_meta_table$Sample)
    
    updateSelectizeInput(session, "in_step29_serp_enrich_input", choices = select_choices, selected = select_choices[1], server = TRUE)
    updateSelectizeInput(session, "in_step29_serp_enrich_ip", choices = select_choices, selected = select_choices[1], server = TRUE)
    updateSelectizeInput(session, "in_step29_serp_enrich_mock", choices = select_choices, selected = select_choices[1], server = TRUE)
    
    return(ribo_meta_table)
  })
  
  ## output the table
  observeEvent(input$act_step29_import_serp_enrich_meta, {
    # browser()
    ## output the table
    output$out_step29_serp_enrich_meta <- DT::renderDataTable({
      DT::datatable(step29_ribo_meta(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "metagene_table"), 
                                     list(extend = 'excel', filename = "metagene_table"))
                    ))
      
    })
  })
  
  
  ### 29.2 calculate the enrichment #################
  step29_enrich_meta <- reactive({
    
    req(input$act_step29_calc_serp_enrich)
    
    if (is.null(step29_ribo_meta())) {return(NULL)}
    
    # browser()
    ribo_meta_table <- step29_ribo_meta()
    
    source("R/calc_Enrichment.R")
    # calculate the enrichment and save to ip_wt_enrich/ip_mock_enrich table
    ip_wt_enrich = calc_enrichment(meta_table = ribo_meta_table,
                                   group = "Sample",
                                   group1 = input$in_step29_serp_enrich_ip,
                                   group2 = input$in_step29_serp_enrich_input,
                                   category = "IP/INPUT",
                                   method = input$in_step29_serp_enrich_method)
    
    mock_wt_enrich = calc_enrichment(meta_table = ribo_meta_table,
                                     group = "Sample",
                                     group1 = input$in_step29_serp_enrich_mock,
                                     group2 = input$in_step29_serp_enrich_input,
                                     category = "MOCK/INPUT",
                                     method = input$in_step29_serp_enrich_method)
    
    meta_enrich <- rbind(ip_wt_enrich, mock_wt_enrich) %>% 
      dplyr::arrange(Category, Paired, Meta, Nucleotide, Codon) %>% 
      dplyr::select(Category, Paired, Meta, Nucleotide, Codon, Enrichment, Samples1, Samples2, Density1, Density2)
    
    # calculate the average of codon
    # if (input$in_step29_serp_enrich_unit == 'aa') {
    #   meta_enrich <- meta_enrich %>%
    #     dplyr::group_by(Meta, Codon, Samples1, Samples2, Density1, Density2, Category, Paired) %>% 
    #     dplyr::reframe(Nucleotide = mean(Nucleotide),
    #                    Enrichment = mean(Enrichment)) %>% 
    #     dplyr::ungroup() %>% 
    #     dplyr::arrange(Category, Paired, Meta, Nucleotide, Codon) %>% 
    #     dplyr::select(Category, Paired, Meta, Nucleotide, Codon, Enrichment, Samples1, Samples2, Density1, Density2)
    #   
    # } else if (input$in_step29_serp_enrich_unit == 'nt') {
    #   meta_enrich <- meta_enrich %>% 
    #     dplyr::select(Category, Paired, Meta, Nucleotide, Codon, Enrichment, Samples1, Samples2, Density1, Density2)
    # }
    
    return(meta_enrich)
  })
  
  ## output the table
  observeEvent(input$act_step29_calc_serp_enrich, {
    # browser()
    ## output the table
    output$out_step29_serp_enrich_meta_enrich <- DT::renderDataTable({
      DT::datatable(step29_enrich_meta(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 6,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "enrich_metagene_table"), 
                                     list(extend = 'excel', filename = "enrich_metagene_table"))
                    ))
    })
  })
  
  ### 29.3 draw the Ribo-seq metaplot #################
  
  step29_meta_enrich_plot <- reactive({
    req(input$act_step29_serp_enrich_draw_plot)
    
    if (is.null(step29_enrich_meta())) {return(NULL)}
    
    # browser()
    
    # filter the range of meta enrich data
    meta_enrich_table <- step29_enrich_meta() %>% 
      dplyr::filter(Meta == "TIS" & Codon >= input$in_step29_serp_enrich_tis_min & Codon <= input$in_step29_serp_enrich_tis_max |
                      Meta == "TTS" & Codon >= input$in_step29_serp_enrich_tts_min & Codon <= input$in_step29_serp_enrich_tts_max)
    
    # browser()
    source("R/draw_meta_enrich.R")
    
    if (input$in_step29_serp_enrich_unit == 'aa') {
      x_axis <- "Codon"
    } else if (input$in_step29_serp_enrich_unit == 'nt') {
      x_axis <- "Nucleotide"
    }
    
    meta_enrich_plot <- draw_meta_enrich(meta = meta_enrich_table,
                                         x = x_axis, y = "Enrichment",
                                         font_size = input$in_step29_serp_enrich_font_size,
                                         line_color = input$in_step29_serp_enrich_color,
                                         line_width = input$in_step29_serp_enrich_width,
                                         cl_alpha = input$in_step29_serp_enrich_alpha,
                                         group = "Category", wrap_group = "Category",
                                         facet = input$in_step29_serp_enrich_facet,
                                         xlabel = x_axis, ylabel = "Enrichment [A.U.]")
    
    
    return(meta_enrich_plot)
    
  })
  
  ## output the enrichment plot
  observeEvent(input$act_step29_serp_enrich_draw_plot, {
    # browser()
    ## show the enrichment plot
    output$out_step29_serp_enrich_metaplot <- renderPlot(
      width = input$out_step29_serp_enrich_fig_width * 100,
      height = input$out_step29_serp_enrich_fig_height * 100,
      {step29_meta_enrich_plot()})
    
    # save the enrichment plot
    output$save_step29_serp_enrich_meta <- downloadHandler(
      filename = function() {
        paste0(input$out_step29_serp_enrich_out_name, "-meta-enrich-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step29_meta_enrich_plot(), filename = file,
               width = input$out_step29_serp_enrich_fig_width, height = input$out_step29_serp_enrich_fig_height)
      }
    )
    
  })
  
  
  
  ############################################################
  ## step30 SeRP peaks ###############################
  
  ### 30.1 import the enrich peaks table  #################
  step30_enrich_peaks <- reactive({
    
    req(input$act_step30_import_serp_enrich_peaks)
    
    if (is.null(input$in_step30_serp_enrich_peaks)) {return(NULL)}
    
    enrich_peaks_table <- read.table(input$in_step30_serp_enrich_peaks$datapath, 
                                     sep = '\t', header = TRUE, stringsAsFactors = FALSE) %>% 
      dplyr::mutate(gene_start = from_tis - min(from_tis))
    
    uniq_gene_names <- unique(enrich_peaks_table$name)
    updateSelectizeInput(session, "in_step30_serp_enrich_gene", choices = uniq_gene_names, selected = uniq_gene_names[1], server = TRUE)
    
    return(enrich_peaks_table)
  })
  
  ## output the table
  observeEvent(input$act_step30_import_serp_enrich_peaks, {
    # browser()
    ## output the table
    output$out_step30_serp_enrich_peaks <- DT::renderDataTable({
      DT::datatable(step30_enrich_peaks(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "enrich_peaks"), 
                                     list(extend = 'excel', filename = "enrich_peaks"))
                    ))
      
    })
  })
  
  
  ### 30.2 draw the enrich peaks plot #################
  step30_enrich_peaks_plot <- reactive({
    
    req(input$act_step30_serp_enrich_draw_plot)
    
    if (is.null(step30_enrich_peaks())) {return(NULL)}
    
    # filter the gene data
    enrich_gene_name <- input$in_step30_serp_enrich_gene
    enrich_gene_data <- step30_enrich_peaks() %>% dplyr::filter(name == enrich_gene_name)
    
    # browser()
    # set the min and max
    if (is.na(input$in_step30_serp_enrich_x_min)) {
      xmin <- min(enrich_gene_data$gene_start)
    } else {
      xmin <- input$in_step30_serp_enrich_x_min
    }
    
    if (is.na(input$in_step30_serp_enrich_x_max)){
      xmax <- max(enrich_gene_data$gene_start)
    } else {
      xmax <- input$in_step30_serp_enrich_x_max
    }
    
    # set the gene length
    utr5_length <- enrich_gene_data %>% dplyr::filter(region == "5utr") %>% NROW()
    cds_length <- enrich_gene_data %>% dplyr::filter(region == "cds") %>% NROW()
    utr3_length <- enrich_gene_data %>% dplyr::filter(region == "3utr") %>% NROW()
    
    # draw the isoforms schematic
    source("R/draw_Isoforms_Schematic.R")
    isoforms_ribo_schematic <- draw_isoforms_schematic(gene_name = enrich_gene_name,
                                                       utr5_length = utr5_length,
                                                       cds_length = cds_length,
                                                       utr3_length = utr3_length,
                                                       gene_color = c("#75aadb", "#264abd", "#75aadb"),
                                                       gene_width = c(2, 5, 2),
                                                       fill_alpha = input$in_step30_serp_enrich_alpha,
                                                       xmin = xmin,
                                                       xmax = xmax,
                                                       xbreaks = input$in_step30_serp_enrich_x_break,
                                                       arrow_length = 0.1,
                                                       arrow_width = 0.8,
                                                       arrow_color = "white",
                                                       font_size = input$in_step30_serp_enrich_font_size) +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    
    # browser()
    source("R/draw_SeRP_Peaks.R")
    
    # filter the range of enrich peaks data
    enrich_peaks_table <- enrich_gene_data %>%
      dplyr::filter(gene_start >= xmin) %>%
      dplyr::filter(gene_start <= xmax)
    
    line_color = c(input$in_step30_serp_enrich_enrich_color1,
                   input$in_step30_serp_enrich_enrich_color2,
                   input$in_step30_serp_enrich_enrich_color3)
    
    enrich_peaks_plot <- draw_serp_peaks(gene_name = enrich_gene_name,
                                         gene_peaks = enrich_peaks_table,
                                         
                                         x = "gene_start",
                                         enrich = "enrich",
                                         edge = "edge",
                                         bound = "bound",
                                         
                                         xlabel = input$in_step30_serp_enrich_xlabel,
                                         ylabel = input$in_step30_serp_enrich_ylabel,
                                         title = input$in_step30_serp_enrich_title,
                                         
                                         line_color = line_color,
                                         
                                         line_width = input$in_step30_serp_enrich_line_width,
                                         fill_alpha = input$in_step30_serp_enrich_alpha,
                                         
                                         xmin = xmin,
                                         xmax = xmax,
                                         xbreaks = input$in_step30_serp_enrich_x_break,
                                         
                                         ymin = input$in_step30_serp_enrich_y_min,
                                         ymax = input$in_step30_serp_enrich_y_max,
                                         ybreaks = input$in_step30_serp_enrich_y_break,
                                         
                                         grid_width = input$in_step30_serp_enrich_grid_width,
                                         grid_color = input$in_step30_serp_enrich_grid_color,
                                         sqrt = input$in_step30_serp_enrich_sqrt,
                                         
                                         font_size = input$in_step30_serp_enrich_font_size)
    
    # combine the gene schematic and the enrichment peaks
    isoforms_enrich_peaks_plot <- cowplot::plot_grid(enrich_peaks_plot, isoforms_ribo_schematic,
                                                     align = 'hv', nrow = 2, axis = 'tblr',
                                                     rel_heights = c(3, 1))
    
    return(isoforms_enrich_peaks_plot)
    
  })
  
  ## output the enrichment plot
  observeEvent(input$act_step30_serp_enrich_draw_plot, {
    # browser()
    ## show the enrichment plot
    output$out_step30_serp_enrich_peaks_plot <- renderPlot(
      width = input$out_step30_serp_enrich_fig_width * 100,
      height = input$out_step30_serp_enrich_fig_height * 100,
      {step30_enrich_peaks_plot()})
    
    # save the enrichment plot
    output$save_step30_serp_enrich_peaks <- downloadHandler(
      filename = function() {
        paste0(input$out_step30_serp_enrich_out_name, "-enrich-peaks-plot-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step30_enrich_peaks_plot(), filename = file,
               width = input$out_step30_serp_enrich_fig_width, height = input$out_step30_serp_enrich_fig_height)
      }
    )
    
  })
  
  
  ############################################################
  ## step31 SeRP motif ###############################
  
  ### 31.1 import the enrich peaks table  #################
  step31_serp_sequence <- reactive({
    
    req(input$act_step31_import_serp_sequence)
    
    # import the serp sequence
    if (is.null(input$in_step31_serp_sequence)) {return(NULL)}
    
    enrich_peaks_seq <- read.table(input$in_step31_serp_sequence$datapath,
                                   sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    # 
    # enrich_peaks_seq <- fread(input$in_step31_serp_sequence$datapath, header = TRUE, stringsAsFactors = FALSE) %>% 
    #   as.data.frame()
    
    # browser()
    
    # filter the minimum length of the sequence
    enrich_peaks_seq <- enrich_peaks_seq %>%
      dplyr::filter(nchar(upstream_aa) >= input$in_step31_serp_upstream_length) %>%
      dplyr::filter(nchar(peak_aa) >= input$in_step31_serp_peak_length) %>% 
      dplyr::filter(nchar(downstream_aa) >= input$in_step31_serp_downstream_length)
    
    
    # # trim the length of the upstream sequence
    # if (input$in_step31_serp_upstream_length > 0) {
    #   enrich_peaks_seq <- enrich_peaks_seq %>%
    #     dplyr::mutate(upstream_nt = substr(upstream_nt, nchar(upstream_nt) - input$in_step31_serp_upstream_length * 3 + 1, nchar(upstream_nt))) %>%
    #     dplyr::mutate(upstream_aa = substr(upstream_aa, nchar(upstream_aa) - input$in_step31_serp_upstream_length + 1, nchar(upstream_aa))) %>% 
    #     dplyr::mutate(upstream_merge_nt_site = paste0(-input$in_step31_serp_upstream_length * 3 + 1, input$in_step31_serp_peak_length * 3))
    #   
    # } else if (input$in_step31_serp_upstream_length <= 0) {
    #   enrich_peaks_seq <- enrich_peaks_seq %>%
    #     dplyr::mutate(upstream_nt = "") %>%
    #     dplyr::mutate(upstream_aa = "")
    # }
    # 
    # # trim the length of the upstream sequence
    # if (input$in_step31_serp_downstream_length > 0) {
    #   enrich_peaks_seq <- enrich_peaks_seq %>%
    #     dplyr::mutate(downstream_nt = substr(downstream_nt, 1, input$in_step31_serp_downstream_length * 3))%>%
    #     dplyr::mutate(downstream_aa = substr(downstream_aa, 1, input$in_step31_serp_downstream_length))
    #   
    # } else if (input$in_step31_serp_upstream_length <= 0) {
    #   enrich_peaks_seq <- enrich_peaks_seq %>%
    #     dplyr::mutate(downstream_nt = "") %>%
    #     dplyr::mutate(downstream_aa = "")
    # }
    # 
    # 
    # # trim the length of the peak sequence
    # if (input$in_step31_serp_peak_length > 0) {
    #   enrich_peaks_seq <- enrich_peaks_seq %>%
    #     dplyr::mutate(right_peak_nt = substr(peak_nt, nchar(peak_nt) - input$in_step31_serp_peak_length * 3 + 1, nchar(peak_nt))) %>%
    #     dplyr::mutate(right_peak_aa = substr(peak_aa, nchar(peak_aa) - input$in_step31_serp_peak_length + 1, nchar(peak_aa))) %>%
    #     
    #     dplyr::mutate(left_peak_nt = substr(peak_nt, 1, input$in_step31_serp_peak_length * 3)) %>%
    #     dplyr::mutate(left_peak_aa = substr(peak_aa, 1, input$in_step31_serp_peak_length))
    #   
    # } else if (input$in_step31_serp_peak_length <= 0) {
    #   enrich_peaks_seq <- enrich_peaks_seq %>%
    #     dplyr::mutate(right_peak_nt = "") %>%
    #     dplyr::mutate(right_peak_aa = "") %>%
    #   
    #     dplyr::mutate(left_peak_nt = "") %>%
    #     dplyr::mutate(left_peak_aa = "")
    # }
    # 
    # 
    # # merge the length of the down-stream and up-stream peak sequence
    # enrich_peaks_seq <- enrich_peaks_seq %>%
    #   dplyr::mutate(upstream_merge_nt_site = paste0(-input$in_step31_serp_upstream_length * 3 + 1, input$in_step31_serp_peak_length * 3)) %>% 
    #   dplyr::mutate(upstream_merge_nt = paste0(upstream_nt, substr(peak_nt, 1, input$in_step31_serp_peak_length * 3))) %>% 
    #   
    #   dplyr::mutate(downstream_merge_nt_site = paste0(-input$in_step31_serp_peak_length * 3, input$in_step31_serp_downstream_length * 3 - 1)) %>% 
    #   dplyr::mutate(downstream_merge_nt = paste0(substr(peak_nt, nchar(peak_nt) - input$in_step31_serp_peak_length * 3 + 1, nchar(peak_nt)), downstream_nt)) %>%
    #   
    #   dplyr::mutate(upstream_merge_aa = paste0(upstream_aa, substr(peak_aa, 1, input$in_step31_serp_peak_length))) %>%
    #   dplyr::mutate(downstream_merge_aa = paste0(substr(peak_aa, nchar(peak_aa) - input$in_step31_serp_peak_length + 1, nchar(peak_aa)), downstream_aa))
    # 
    

    enrich_peaks_seq <- enrich_peaks_seq %>%
      dplyr::mutate(upstream_nt = substr(upstream_nt, nchar(upstream_nt) - input$in_step31_serp_upstream_length * 3 + 1, nchar(upstream_nt))) %>%
      dplyr::mutate(downstream_nt = substr(downstream_nt, 1, input$in_step31_serp_downstream_length * 3)) %>% 
      
      dplyr::mutate(upstream_aa = substr(upstream_aa, nchar(upstream_aa) - input$in_step31_serp_upstream_length + 1, nchar(upstream_aa))) %>%
      dplyr::mutate(downstream_aa = substr(downstream_aa, 1, input$in_step31_serp_downstream_length)) %>% 
      
      # dplyr::mutate(upstream_merge_nt_site = paste0(-input$in_step31_serp_upstream_length * 3 + 1, input$in_step31_serp_peak_length * 3)) %>% 
      dplyr::mutate(upstream_merge_nt = paste0(upstream_nt, substr(peak_nt, 1, input$in_step31_serp_peak_length * 3))) %>% 
      
      # dplyr::mutate(downstream_merge_nt_site = paste0(-input$in_step31_serp_peak_length * 3, input$in_step31_serp_downstream_length * 3 - 1)) %>% 
      dplyr::mutate(downstream_merge_nt = paste0(substr(peak_nt, nchar(peak_nt) - input$in_step31_serp_peak_length * 3 + 1, nchar(peak_nt)), downstream_nt)) %>%
      
      dplyr::mutate(upstream_merge_aa = paste0(upstream_aa, substr(peak_aa, 1, input$in_step31_serp_peak_length))) %>%
      dplyr::mutate(downstream_merge_aa = paste0(substr(peak_aa, nchar(peak_aa) - input$in_step31_serp_peak_length + 1, nchar(peak_aa)), downstream_aa))
    
    # drop the stop codon 
    enrich_peaks_seq <- enrich_peaks_seq %>%
      dplyr::filter(!str_detect(upstream_aa, "[*]")) %>% 
      dplyr::filter(!str_detect(downstream_aa, "[*]")) %>% 
      dplyr::filter(!str_detect(peak_aa, "[*]"))
    
    return(enrich_peaks_seq)
  })
  
  ## output the table
  observeEvent(input$act_step31_import_serp_sequence, {
    # browser()
    ## output the table
    output$out_step31_serp_peaks_sequence <- DT::renderDataTable({
      DT::datatable(step31_serp_sequence(), rownames = FALSE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "enrich_peak_sequence"), 
                                     list(extend = 'excel', filename = "enrich_peak_sequence"))
                    ))
      
    })
  })
  
  ### 31.2 calculate the consensus matrix #################
  step31_serp_consensus_matrix <- reactive({
    req(input$act_step31_calc_serp_consensus_matrix)
    
    if (is.null(step31_serp_sequence())) {return(NULL)}
    
    # browser()
    
    # convert the nucleotide sequence to the biostrings object
    upstream_merge_nt <- step31_serp_sequence() %>% dplyr::select(upstream_merge_nt) %>% pull()
    downstream_merge_nt <- step31_serp_sequence() %>% dplyr::select(downstream_merge_nt) %>% pull()
    
    # set the range of the sequence
    upstream_merge_nt_site = seq(-input$in_step31_serp_upstream_length * 3 + 1, input$in_step31_serp_peak_length * 3, 1)
    downstream_merge_nt_site = seq(-input$in_step31_serp_peak_length * 3, input$in_step31_serp_downstream_length * 3 - 1, 1)

    upstream_merge_nt_set_pfm <- consensusMatrix(DNAStringSet(upstream_merge_nt), baseOnly = TRUE) %>% 
      subset(., rowSums(.) > 0) %>% 
      as.data.frame() %>% 
      dplyr::rename_all(~as.character(upstream_merge_nt_site)) %>% 
      as.matrix()
    
    downstream_merge_nt_set_pfm <- consensusMatrix(DNAStringSet(downstream_merge_nt), baseOnly = TRUE) %>% 
      subset(., rowSums(.) > 0) %>% 
      as.data.frame() %>% 
      dplyr::rename_all(~as.character(downstream_merge_nt_site)) %>% 
      as.matrix()
    
    # convert the amino acid sequence to the biostrings object
    upstream_merge_aa <- step31_serp_sequence() %>% dplyr::select(upstream_merge_aa) %>% pull()
    downstream_merge_aa <- step31_serp_sequence() %>% dplyr::select(downstream_merge_aa) %>% pull()
    
    # set the range of the sequence
    upstream_merge_aa_site = seq(-input$in_step31_serp_upstream_length + 1, input$in_step31_serp_peak_length, 1)
    downstream_merge_aa_site = seq(-input$in_step31_serp_peak_length, input$in_step31_serp_downstream_length - 1, 1)
    
    upstream_merge_aa_set_pfm <- consensusMatrix(AAStringSet(upstream_merge_aa), baseOnly = TRUE) %>% 
      subset(., rowSums(.) > 0) %>% 
      as.data.frame() %>% 
      dplyr::rename_all(~as.character(upstream_merge_aa_site)) %>% 
      as.matrix()
    
    downstream_merge_aa_set_pfm <- consensusMatrix(AAStringSet(downstream_merge_aa), baseOnly = TRUE) %>% 
      subset(., rowSums(.) > 0) %>% 
      as.data.frame() %>% 
      dplyr::rename_all(~as.character(downstream_merge_aa_site)) %>% 
      as.matrix()
    
    # browser()
    source("R/calc_PFM_PPM_PWM.R")
    
    # format the pfm/pwm/ppm
    if (input$in_step31_serp_method == "custom") {
      upstream_merge_nt_set_list <- calc_pfm_ppm_pwm(seq_pfm = upstream_merge_nt_set_pfm, seq_freq = 0.25, 
                                                     method = input$in_step31_serp_custom)
      downstream_merge_nt_set_list <- calc_pfm_ppm_pwm(seq_pfm = downstream_merge_nt_set_pfm, seq_freq = 0.25, 
                                                       method = input$in_step31_serp_custom)
      upstream_merge_aa_set_list <- calc_pfm_ppm_pwm(seq_pfm = upstream_merge_aa_set_pfm, seq_freq = 0.05, 
                                                     method = input$in_step31_serp_custom)
      downstream_merge_aa_set_list <- calc_pfm_ppm_pwm(seq_pfm = downstream_merge_aa_set_pfm, seq_freq = 0.05, 
                                                       method = input$in_step31_serp_custom)
      
      upstream_merge_nt_set_mat <- upstream_merge_nt_set_list$seq_mat
      downstream_merge_nt_set_mat <- downstream_merge_nt_set_list$seq_mat
      upstream_merge_aa_set_mat <- upstream_merge_aa_set_list$seq_mat
      downstream_merge_aa_set_mat <- downstream_merge_aa_set_list$seq_mat
      
    } else {
      upstream_merge_nt_set_mat <- upstream_merge_nt_set_pfm
      downstream_merge_nt_set_mat <- downstream_merge_nt_set_pfm
      
      upstream_merge_aa_set_mat <- upstream_merge_aa_set_pfm
      downstream_merge_aa_set_mat <- downstream_merge_aa_set_pfm
    }
    
    return(list(upstream_merge_nt_set_mat = upstream_merge_nt_set_mat, 
                downstream_merge_nt_set_mat = downstream_merge_nt_set_mat, 
                upstream_merge_aa_set_mat = upstream_merge_aa_set_mat,
                downstream_merge_aa_set_mat = downstream_merge_aa_set_mat))
    
  })
  
  
  ## output the table
  observeEvent(input$act_step31_calc_serp_consensus_matrix, {
    # browser()
    ## output the table
    output$out_step31_serp_peaks_upstream_nt <- DT::renderDataTable({
      DT::datatable(step31_serp_consensus_matrix()$upstream_merge_nt_set_mat, rownames = TRUE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "upstream_merge_nt_set_mat"), 
                                     list(extend = 'excel', filename = "upstream_merge_nt_set_mat"))
                    ))
    })
    
    output$out_step31_serp_peaks_downstream_nt <- DT::renderDataTable({
      DT::datatable(step31_serp_consensus_matrix()$downstream_merge_nt_set_mat, rownames = TRUE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "downstream_merge_nt_set_mat"), 
                                     list(extend = 'excel', filename = "downstream_merge_nt_set_mat"))
                    ))
    })
    
    output$out_step31_serp_peaks_upstream_aa <- DT::renderDataTable({
      DT::datatable(step31_serp_consensus_matrix()$upstream_merge_aa_set_mat, rownames = TRUE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "upstream_merge_aa_set_mat"), 
                                     list(extend = 'excel', filename = "upstream_merge_aa_set_mat"))
                    ))
    })
    
    output$out_step31_serp_peaks_downstream_aa <- DT::renderDataTable({
      DT::datatable(step31_serp_consensus_matrix()$downstream_merge_aa_set_mat, rownames = TRUE, filter = "top", extensions = 'Buttons',
                    options = list(pageLength = 10,
                                   lengthMenu = c(10, 50, 100, 500, 1000),
                                   scrollX = TRUE, 
                                   dom = 'lBfrtipc',
                                   buttons = list(
                                     list(extend = 'csv', filename = "downstream_merge_aa_set_mat"), 
                                     list(extend = 'excel', filename = "downstream_merge_aa_set_mat"))
                    ))
    })
  })
  
  
  
  ### 31.3 draw the seqlogo ##############################
  step31_serp_seqlogo <- reactive({
    req(input$act_step31_draw_serp_seqlogo)
    
    if (is.null(step31_serp_consensus_matrix())) {return(NULL)}
    
    # browser()
    
    source("R/draw_seqlogo.R")
    
    # retrieve the matrix
    if (input$in_step31_serp_region == "upstream" && input$in_step31_serp_seqtype == "rna") {
      consensus_matrix <- step31_serp_consensus_matrix()$upstream_merge_nt_set_mat
      
    } else if (input$in_step31_serp_region == "upstream" && input$in_step31_serp_seqtype == "aa") {
      consensus_matrix <- step31_serp_consensus_matrix()$upstream_merge_aa_set_mat
      
    } else if (input$in_step31_serp_region == "downstream" && input$in_step31_serp_seqtype == "rna") {
      consensus_matrix <- step31_serp_consensus_matrix()$downstream_merge_nt_set_mat
      
    } else if (input$in_step31_serp_region == "downstream" && input$in_step31_serp_seqtype == "aa") {
      consensus_matrix <- step31_serp_consensus_matrix()$downstream_merge_aa_set_mat
      
    }
    
    xrange <- as.numeric(colnames(consensus_matrix))
    
    xstart <- ifelse(is.na(input$in_step31_serp_seqlogo_x_min), min(xrange), input$in_step31_serp_seqlogo_x_min)
    xend <- ifelse(is.na(input$in_step31_serp_seqlogo_x_max), max(xrange), input$in_step31_serp_seqlogo_x_max)
    
    consensus_matrix <- consensus_matrix %>% 
      as.data.frame() %>% 
      dplyr::select(as.character(xstart:xend)) %>% 
      as.matrix()
    
    # browser()
    serp_sequence_plot <- draw_seqlogo(pwm = consensus_matrix,
                                       method = input$in_step31_serp_method, seq_type = input$in_step31_serp_seqtype,
                                       col_scheme = input$in_step31_serp_fill_color, col_alpha = input$in_step31_serp_fill_alpha,
                                       stack_width = input$in_step31_serp_font_stack, 
                                       font_family = input$in_step31_serp_font_family, font_size = input$in_step31_serp_font_size,
                                       xstart = xstart, xend = xend, x_breaks = input$in_step31_serp_seqlogo_x_breaks,
                                       xlabel = "site", ylabel = input$in_step31_serp_custom, title = "SeRP seqlogo")
    return(serp_sequence_plot)
    
  })
  
  
  ## output the alignment plot
  observeEvent(input$act_step31_draw_serp_seqlogo, {
    # browser()
    ## show the alignment plot
    output$out_step31_serp_peaks_seqlogo <- renderPlot(
      width = input$out_step31_serp_fig_width * 100,
      height = input$out_step31_serp_fig_height * 100,
      {step31_serp_seqlogo()})
    
    # save the alignment plot
    output$save_step31_serp_seqlogo <- downloadHandler(
      filename = function() {
        paste0(input$out_step31_serp_fig_name, "-peaks-seqlogo-", Sys.Date(), ".pdf")
      },
      content = function(file) {
        ggsave(plot = step31_serp_seqlogo(), filename = file,
               width = input$out_step31_serp_fig_width, height = input$out_step31_serp_fig_height)
      }
    )
    
  })
  
  
  ############################################################
  ## step32 smORF ###############################
  ############################################################
  
  
  
  
  ############################################################
  ## step33 smORF expression ###############################
  ############################################################
  
  
  
  
  ############################################################
  ## step34 smORF and mORF association ###############################
  ############################################################
  
  
  
  
  
  
  
  
  
  
  
}
