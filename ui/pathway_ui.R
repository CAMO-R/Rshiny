
pathway_ui <- function(id, label = "Aggregated Pathway") {
  ns <- NS(id)
  tabPanel("Pathway-based analysis", value=id,
           sidebarLayout(
             sidebarPanel(
               h2("Pathway c-scores & d-scores calculation"),
               selectInput(ns("measure"), label = "Score type", 
                           choices = list("F-measure (recommended)" = "Fmeasure",
                                          "Youden index" = "youden", 
                                          "Geometric mean" = "geo.mean")),
               numericInput(ns("permNumPathway"), 
                            label = "Number of permutations", 
                            value = 100),
               checkboxInput(ns("parallel"), 'Use parallel computation', FALSE),
               uiOutput(ns("para_cores")),
               br(),
               radioButtons(ns("select_pathwayfile"), label = "Pathway database",
                       choices = list("Select from existing pathway database"="exist",
                                      "Upload a list of pathways"="upload"),
                       selected = "upload"),
               uiOutput(ns("exist_pathwayfile")),
               uiOutput(ns("upload_pathwayfile")),

               
               ## pathway selection advanced setting
               bsCollapsePanel("Advanced settings for pathway selection",
                               NULL,
                               numericInput(ns("pathwaysizeLowerCut"), 
                                            label = "Minimum pathway size", value = 5),
                               numericInput(ns("pathwaysizeUpperCut"), 
                                            label = "Maximum pathway size", value = 200),
                               numericInput(ns("overlapsizeCut"), 
                                            label = "Lower bound of the minimum number of overlapping genes across studies", value = 5),
                               numericInput(ns("medDECut"), 
                                            label = "Lower bound of the median number of overlapping DE genes across studies", value = 3),
                               numericInput(ns("minDECut"), 
                                            label = "Lower bound of the minimum number of overlapping DE genes across studies", value = 0),
                               numericInput(ns("qfisherCut"), 
                                            label = "Upper bound of the Fisher combination q-value", value = 0.05),
                               textInput(ns("topPathNnum"), 
                                            label = "If only top pathways in at least one study are considered (use this if Fisher combination q-value is too stringent), please input the number of top pathways.", value = ""),
                               style="primary"
               ),
               
               actionButton(ns('ACS_ADS'), 'Pathway c-scores & d-scores', 
                            class="btn-success", icon = icon("play")),
               
               ##tunning K
               tags$hr(),
               h2("Pathway clustering"),
               actionButton(ns('tuneK_ACS'), 
                            'Scree plot to determine the number of pathway clusters', 
                            class="btn-success",
                            icon = icon("play")),
               hr(),
               numericInput(ns("K_ACS"), 
                            label = "Optimal number of pathway clusters K", value = 1),
               selectInput(ns('select_hashtbfile'), 'Select a noun-pathway matrix for text mining', 
                           c("Noun pharases from KEGG and Reactome (homo sapiens)"="hsa_text",
                             "Noun pharases from KEGG and Reactome (mus musculus)"="mmu_text",
                             "Noun pharases from KEGG and Reactome (rattus norvegicus)"="rno_text",
                             "Noun pharases from KEGG and Reactome (caenorhabditis elegans)"="cel_text",
                             "Noun pharases from KEGG and Reactome (drosophila melanogaster)"="dme_text",
                             "Upload a noun-pathway matrix file"="upload",
                             "Skip text mining"="no_textmining"),
                           selected = NULL),
               uiOutput(ns("upload_hashtbfile")),

               ## pathway clustering advanced setting
               bsCollapsePanel("Advanced settings for pathway clustering", 
                               NULL, 
                               numericInput(ns("silCut"), "Silhouette index cutoff to control scatterness", 0.1),
                               #numericInput(ns("tmThres"), "FDR control for key words from text mining", 0.2),
                               numericInput(ns("comProbCut"), "Lower bound of co-membership proportion shown in heatmap", 0.7),
                               style="primary"
               ),
               actionButton(ns('pathclust_ACS'), 
                            'Pathway clustering', 
                            class="btn-success",
                            icon = icon("play")),
               hr(),
               h2("Visualization of pathway DE evidence and c/d-scores"),
               uiOutput(ns("ACS_DEgroup")),
               actionButton(ns('plotACS_DE'), 
                            'DE evidence and c/d-scores plot', 
                            class="btn-success",
                            icon = icon("play"))
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel("Table of pathway-level c-scores and d-scores",
                          h3("Pathway-level c-scores"),
                          DT::dataTableOutput(ns("pathwayACS_Table")),
                          h3("Pathway-level d-scores"),
                          DT::dataTableOutput(ns("pathwayADS_Table"))),
                 tabPanel("Scree plot",
                          h3("Scree plot to select the optimal K"),
                          imageOutput(ns("tuneKFig"), height = 500)),
                 tabPanel("Visualization of pathway clusters",
                          h3("MDS map of all pathways based on pathway-level c-scores colored by cluster labels"),
                          imageOutput(ns("mdsPathway"), height = 500),
                          h3("Heatmap of pathway-level c-score pvalues ordered by cluster labels"),
                          imageOutput(ns("heatmapPathway"),height = 500)),
                 tabPanel("Text mining on each pathway cluster",
                          h3("Table of enriched key words & Co-membership heatmap of each cluster"),
                          #h5("Key words table: noun phrases in each pathway cluster"),
                          h5("Co-membership probability: the proportion of pathways in which the two studies were assigned together (based on study similarity for each pathway in individual pathway analysis) in each pathway cluster."),
                          uiOutput(ns("comemPlots"))),
                 tabPanel("DE evidence and c/d-scores plot",
                          h3("c-scores/d-scores - DE evidence plot"),
                          numericInput(ns("numNearPath"),
                                       "Select the number of nearest pathways to display after clicking",
                                       3,min = 1, step = 1),
                          textOutput(ns("nearPathwayName")),
                          uiOutput(ns("clickedPathway")),
                          uiOutput(ns("ACS_DEplots")))
               )
             )
           )
  )
}
