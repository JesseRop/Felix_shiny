setwd("D:/GCRF_UoG/Felix_shiny/")
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
                              
                            ),
                            tags$hr(),
                          wellPanel(style = "background:LightBlue",
                                    tags$hr(),
                                    tags$p(style = "font-family:Comic Sans MS;color:black",
                                           "Visualization of gene expression of one or more genes in the different cell populations using Violin plots and UMAP feature plots."
                                    ),
                                    tags$hr()
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
                       withSpinner(plotOutput("dotplot")),
                       tags$hr(),
                       wellPanel(style = "background:LightBlue",
                                 tags$hr(),
                                 tags$p(style = "font-family:Comic Sans MS;color:black",
                                        #paste("Comparison of", input$select_markers_dotplot[1], input$select_markers_dotplot[2], "expression between", conditions[1], "and", conditions[2], "across clusters using a dotplot. Red for Wildtype and Blue for CD18 KO cells with the increasing in intensity of the respective colour correlating with the level of gene expression in the cluster"),
                                        # paste(unlist(for (i in 1:length(input$select_markers_dotplot)) {
                                        #   print(paste(input$select_markers_dotplot[i]))
                                        # }
                                        paste("Comparison of gene expression across cell populations using a dotplot. The genes are on the x-axis and the populations on the y-axis. The increase in intensity from grey to blue correlates with the increase in the average level of gene expression across all cells in the cluster. The size of the dot corresponds to the percentage of cells in the cluster expressing the gene.")
                                        
                                 ),
                                 tags$hr()
                       )
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
                       value = res_default,
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
                     uiOutput("cluster_annot"),
                     # )
                     tags$hr(),
                   wellPanel(style = "background:LightBlue",
                             tags$hr(),
                             tags$p(style = "font-family:Comic Sans MS;color:black",
                                    #paste("Comparison of", input$select_markers_dotplot[1], input$select_markers_dotplot[2], "expression between", conditions[1], "and", conditions[2], "across clusters using a dotplot. Red for Wildtype and Blue for CD18 KO cells with the increasing in intensity of the respective colour correlating with the level of gene expression in the cluster"),
                                    # paste(unlist(for (i in 1:length(input$select_markers_dotplot)) {
                                    #   print(paste(input$select_markers_dotplot[i]))
                                    # }
                                    paste("Labelling of cell populations based on top cluster markers in the 'Cluster marker' table on the left.")
                                    
                             ),
                             tags$hr()
                   )
                   )
                   ),
                 
                 column(
                   7, 
                   wellPanel(
                     # box(
                     h4("Cluster markers"),
                     # width = NULL,
                     # solidHeader = TRUE,
                     withSpinner(DTOutput("marker_table")),
                     tags$hr(),
                   wellPanel(style = "background:LightBlue",
                             tags$hr(),
                             tags$p(style = "font-family:Comic Sans MS;color:black",
                                    #paste("Comparison of", input$select_markers_dotplot[1], input$select_markers_dotplot[2], "expression between", conditions[1], "and", conditions[2], "across clusters using a dotplot. Red for Wildtype and Blue for CD18 KO cells with the increasing in intensity of the respective colour correlating with the level of gene expression in the cluster"),
                                    # paste(unlist(for (i in 1:length(input$select_markers_dotplot)) {
                                    #   print(paste(input$select_markers_dotplot[i]))
                                    # }
                                    paste("Table of markers of cluster/s selected from the 'ScRNA-seq cluster' column. Here, markers are genes highly expressed in a cluster as compared to all other clusters. A set of top marker genes in cluster can be used to infer its identity.")
                                    
                             ),
                             tags$hr()
                   )
                     
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
        "% in interest cluster",
        "% in other clusters",
        "Adjusted P"
      ),
      selection = "single",
      rownames = FALSE,
      filter = list(position = "top", plain = TRUE),
      style = "default",
      extensions = c("Buttons","Scroller"),
      options = list(
        deferRender = TRUE,
        scrollY = 350,
        #scroller = TRUE,
        #lengthMenu = FALSE,
        autoWidth = FALSE,
        dom = "Blfrtip",
        buttons = 
          list(list(
            extend = "collection",
            buttons = c("csv", "pdf"),
            text = "Download"
          )),  # end of buttons customization
        
        # customize the length menu
        lengthMenu = list( c(10, 20, -1), c(10, 20, "All")), # end of lengthMenu customization
        pageLength = 10
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
