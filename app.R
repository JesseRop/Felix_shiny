source("leish_lib_files_loading.R")


ui = navbarPage(title = "Leishmania scRNA-seq visualization", theme = shinytheme("cerulean"), 
               #tags$head(tags$style("div.dataTables_scrollHead span {color: black}")),
               #h2("Leishmania scRNA-seq visualization"),
               tabPanel("Gene expression",wellPanel(
                 fluidRow(
                   #h3("Gene expression"),
                   column(8,
                          wellPanel(
                            fluidRow(
                              column(
                                4,
                                selectInput(
                                  inputId = "marker_genes",
                                  label = strong("Select genes:"),
                                  choices = all_genes,
                                  multiple = T,
                                  selected = "LmxM.33.3645"
                                ),
                              )
                            ),
                            fluidRow(
                              column(
                                6,
                                # box(
                                h4("Violin plot"),
                                #   width = NULL,
                                #   solidHeader = TRUE,
                                withSpinner(plotOutput("vln_plot"))
                                # )
                              ),
                              column(
                                6,
                                # box(
                                h4("UMAP"),
                                # width = NULL,
                                # solidHeader = TRUE,
                                withSpinner(plotOutput("marker_umap"))
                                # )
                              )
                              
                            )
                          )
                   ),
                   
                   column(
                     4,
                     wellPanel(
                       # box(
                       h4("Dotplot"),
                       #   width = NULL,
                       #   solidHeader = TRUE,
                       
                       selectInput(
                         inputId = "marker_genes_dotplot",
                         label = strong("Select genes:"),
                         choices = all_genes,
                         multiple = T,
                         selected = c("LmxM.33.3645","LmxM.15.1240", "LmxM.08.0810", "LmxM.10.0870")
                       ),
                       withSpinner(plotOutput("dotplot"))
                       # )
                     )
                   )
                 )
               )),  
               tabPanel("Clustering",wellPanel(fluidRow(
                 
                 #h3("Clustering"),
                 column(
                   5,
                   wellPanel(
                     
                     # box(
                     h4("UMAP Clusters"),
                     # width = NULL,
                     # solidHeader = TRUE,
                     sliderInput(
                       inputId = "UMAP",
                       label = strong("Clustering resolution"),
                       value = res1,
                       min = res1,
                       max = res2,
                       step = diff_res
                     ),
                     withSpinner(plotOutput("labelled_umap", height = 250))
                     
                     # )
                     ,
                     # box(title = "Cluster labelling based on marker genes",
                     #     width = NULL,
                     #     solidHeader = TRUE,
                     uiOutput("cluster_annot")
                     # )
                   )),
                 
                 column(
                   7, 
                   wellPanel(
                     # box(
                     h4("Cluster markers"),
                     # width = NULL,
                     # solidHeader = TRUE,
                     withSpinner(DTOutput("marker_table"))
                     
                     # )
                   )
                 )
               )
               ))
               
               
)

server = function(input, output) {
  
  leishmania_dta = reactive({
    leishmania_umap_list_res[[(input$UMAP * 10)-1]]
    
  })
  
  
  cluster_marker_table = reactive({
    leishmania_marker_tables_res[[(input$UMAP * 10)-1]] %>% group_by(cluster) %>% select("gene",
                                                                                          "cluster",
                                                                                          "avg_logFC",
                                                                                          "pct.1",
                                                                                          "pct.2",
                                                                                          "p_val_adj") %>% mutate_at(vars(matches("p_val|pval")), ~ formatC(., format = "e", digits = 2))
    
  })
  
  
  output$marker_table = DT::renderDataTable({
    numeric_cols =
      colnames(data.frame(cluster_marker_table()))[which_numeric_cols(data.frame(cluster_marker_table()))]
    
    #   # Javascript-enabled table.
    #   datatable(
    DT::datatable(
      cluster_marker_table(),
      colnames = c(
        "Gene",
        "scRNA-seq cluster",
        "Avg_logFC",
        "% in cluster 1",
        "% in cluster 2",
        "Adjusted P"
      ),
      selection = "single",
      rownames = FALSE,
      filter = list(position = "top", plain = TRUE),
      style = "default",
      extensions = "Scroller",
      options = list(
        deferRender = TRUE,
        scrollY = 350,
        scroller = TRUE,
        lengthMenu = FALSE,
        autoWidth = FALSE
      )
    ) %>%
      DT::formatSignif(columns = numeric_cols, digits = 3)
  }, server = TRUE)
  
  
  ######################
  ### Labelling clusters under subtitle 2.3
  ######################
  #Generating dynamic fields for labelling UMAP clusters and initializing the fields with placeholders
  output$cluster_annot = renderUI({
    #umap_cluster_modified1 = leishmania_dta()
    
    if(length(unique(leishmania_dta()[["seurat_clusters"]][["seurat_clusters"]])) == length(cluster_names)){
      do.call(flowLayout, 
              lapply(0:(length(unique(leishmania_dta()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
                textInput(inputId = paste("labeller", x), label = strong(paste("Input cluster", x, "label")), value = cluster_names[x+1])
                
              })
      )
      
    } else {
      do.call(flowLayout,
              lapply(0:(length(unique(leishmania_dta()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
                textInput(inputId = paste("labeller", x), label = strong(paste("Input cluster", x, "label")), value = x)
                
              })
      )
    }
  })
  
  ##Storing names decided on by the researchers based on optimal clustering to start off the differential expression visualization
  annotation = reactiveValues(annot = cluster_names)
  
  ##Observer to allow updating of cluster names dynamically as typed in
  observe({
    
    req(unlist(lapply(0:(length(unique(leishmania_dta()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
      new_cluster_name = input[[paste("labeller",x)]]
    })))
    annotation$annot = unlist(lapply(0:(length(unique(leishmania_dta()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
      new_cluster_name = input[[paste("labeller",x)]]
    }))
  })
  
  
  ##Renaming clusters
  umap_cluster_labelled = reactive({
    umap_names = annotation$annot
    
    leishmania_dta = leishmania_dta()
    if(length(unique(leishmania_dta()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
      names(umap_names) = levels(leishmania_dta)
      leishmania_dta = RenameIdents(leishmania_dta, umap_names)
      
    } else {
      leishmania_dta
    }
    
  })
  
  output$labelled_umap = renderPlot({
    #HoverLocator(DimPlot(umap_cluster_labelled(), label = TRUE, label.size = 6), information = FetchData(umap_cluster_labelled(), vars = c("ident", "nFeature_RNA")))
    DimPlot(umap_cluster_labelled(), label = TRUE, label.size = 6)
    
  })
  ######################
  ###END Labelling clusters under subtitle 2.3
  ######################
  
  
  ######################
  ### Visualizition of top cluster markers under subtitle 2.3
  ######################
  
  
  output$vln_plot = renderPlot({
    VlnPlot(umap_cluster_labelled(), features = input$marker_genes)
  })
  output$marker_umap = renderPlot({
    FeaturePlot(umap_cluster_labelled(), features = input$marker_genes)
  })
  
  output$dotplot = renderPlot({
    DotPlot(umap_cluster_labelled(), features = input$marker_genes_dotplot, dot.scale = 6) + RotatedAxis()
  })
  ######################
  ### END of visualization of top cluster markers under subtitle 2.3
  ######################
  
  
}



shinyApp(ui = ui, server = server)
