indPathway_ui <- function(id, label = "Individual Pathway") {
  ns <- NS(id)
  tabPanel("Individual pathway analysis", value=id,
           sidebarLayout(
             sidebarPanel(
               h2("Individual pathway results browser"),
               h3("(a) Pathway-level c-scores & d-scores"),
               uiOutput(ns("selectPathway")),
               hr(),
               h3("(b) Visualization of individual pathway"),
               uiOutput(ns("choosePathway")),
               #radioButtons(ns("browser_scoreType"), "Select a score type", choices = list("c-scores" = "c_scores", "d-scores" = "d_scores"), selected = "c_scores"),
               #checkboxInput(ns("browser_useACS"), 'Use c-scores to generate results', TRUE),
               checkboxGroupInput(ns('browserVisualizations'), 'Select types of graph to display',
                                  c("mdsModel - MDS map of studies for the selected pathway"="mdsModel",
                                    "clustModel - Clustering of studies for the selected pathway"="clustModel",
                                    "genePM - Heatmap of posterior probability of differential expression for all genes in the selected pathway"="genePM", 
                                    "keggView - Topological visualization of a KEGG pathway"="keggView",
                                    "reactomeView - Topological visualization of a Reactome pathway"="reactomeView")),
               #checkboxInput(ns("showKegg"), 'Generate and show Kegg topological plots of selected pathway(s)', TRUE),
               uiOutput(ns("KEGG_panel")),
               #checkboxInput(ns("showReactome"), 'Generate and show Reactome topological plots of selected pathway(s)', TRUE),
               uiOutput(ns("Reactome_panel")),
               actionButton(ns("showVisual"), 'Generate visualizations', class = "btn-success", icon = icon("play")),
               hr(),
               h3("(c) Topological gene module detection for KEGG pathways"),
               uiOutput(ns('MD_chooseKEGGname')),
               h5("Select a study pair:"),
               fluidRow(
                 column(width = 3,
                        uiOutput(ns('MD_dataPair1'))
                 ),
                 column(width = 3, 
                        uiOutput(ns('MD_dataPair2'))
                 )
               ),
               textInput(ns("MD_KEGGspecies"), "KEGG organism code", "hsa"),
               radioButtons(ns("MD_searchType"), "Select a searching algorithm", choices = list("Simulated Annealing" = "SA", "Exhaustive Searching" = "Exhaustive"), selected = "Exhaustive"),
               radioButtons(ns("MD_geneType"), "Select a module type", choices = list("Concordant" = "concordant", "Discordant" = "discordant"), selected = "discordant"),
               numericInput(ns("DE_PM_cut"), "Lower bound of the posterior DE probability in both studies", 0.2, min = 0),
               radioButtons(ns("MD_KEGG.dataGisTopologyG"), "Gene nomenclature in the study matches with the KEGG pathway database?", choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "FALSE"),
               uiOutput(ns('MD_choosekeggGenes_dat_map')),
               
               bsCollapsePanel("Advanced settings for topological gene module detection",
                               NULL,
                               textInput(ns("MD_pathwayID"), "KEGG pathway ID without the organism prefix. If not provided, ID will be retrieved from KEGGREST."),
                               numericInput(ns("minM"), "Minimum module size", 4, min = 1),
                               textInput(ns("maxM"), "Maximum module size. If not provided, the algorithm will search up the number of concordant/discordant genes.", ""),
                               numericInput(ns("MD_P"), "Number of permutations", 1000, min = 1),
                               numericInput(ns("MD_cores"), "Number of parallel computing cores", 1, min = 1),
                               numericInput(ns("MD_seed"), "Random seed", 12345, min = 1),
                               textInput(ns("node_sep"), "Separation string to concatenate multiple genes in one node", "-"),
                               hr(),
                               p("Hyperparameters for Simulated Annealing:"),
                               numericInput(ns("Tm0"), "Initial temperature", 10, min = 1),
                               numericInput(ns("mu"), "Temperature multiplier", 0.95),
                               numericInput(ns("epsilon"), "Final temperature", 1e-5),
                               numericInput(ns("N"), "Number of maximum annealing times", 3000, min = 1),
                               numericInput(ns("reps_eachM"), "Number of searching repetitions at each module size", 1, min = 1),
                               numericInput(ns("topG_from_previous"), "Number of top module results stored as initials for next module size", 1, min = 1),
                               style="primary"
               ),
               actionButton(ns('tuneMD'), 'Elbow plot to determine the optimal module size', class="btn-success", icon = icon("play")),
               # p("Indicate module size(s) to plot"),
               textInput(ns("MD_whichToDraw"), "Select module size(s) to show in topological plot", "all"),
               helpIcon(ns('which_modele'), "Please type the module sizes of interests (check elbow plot) and use \",\" as separator (e.g. 1,2,3). Leave it as default if you want all of them."),
               actionButton(ns('plotMD'), 'Topological plot with the detected module(s) highlighted', class="btn-success", icon = icon("play")),
               hr(),
               h2("Save results for all pathways"),
               helpIcon(ns('individual_analysis'), "Results will be saved to the working directory selected."),
               #radioButtons(ns("save_scoreType"), "Select a score type", choices = list("c-scores" = "c_scores", "d-scores" = "d_scores"), selected = "c_scores"),
               
               #checkboxInput(ns("useADS"), 'Use c-scores to generate results', TRUE),
               #uiOutput(ns('selectVisualizations')),
               checkboxGroupInput(ns('selectVisualizations'), 'Select types of graph to save',
                                  c("mdsModel - MDS map of studies for the selected pathway"="mdsModel",
                                    "clustModel - Clustering of studies for the selected pathway"="clustModel",
                                    "genePM - Heatmap of posterior probability of differential expression for all genes in the selected pathway"="genePM", 
                                    "keggView - Topological visualization of a KEGG pathway"="keggView",
                                    "reactomeView - Topological visualization of a Reactome pathway"="reactomeView")),
               uiOutput(ns('keggViewSelect')),
               uiOutput(ns('reactomeViewSelect')),
               actionButton(ns('plotAll'), 'Run and Save', class="btn-success", icon = icon("play"))
             ),
             mainPanel(
               tabsetPanel(
               tabPanel("Table of pathway-level c-scores and d-scores",
                        h3("Table of pathway-level c-scores and d-scores of the selected pathway"),
                        textOutput(ns("Show_selected_pathway_name")),
                        textOutput(ns("Pathway_ACS_ADS_note")),
                        DT::dataTableOutput(ns("selectedPathwayACS_ADS_Table"))
               ),
               tabPanel("Visualization for individual pathway",
                        h3("MDS map, Clustering heatmap, Heatmap of posterior DE probability, Topological plot (KEGG and Reactome pathways)"),
                        uiOutput(ns("inputPathway"))
               ),
               tabPanel("Topological gene module detection for KEGG pathways - elbow plot",
                        h3("Elbow plot to select the optimal module size"),
                        plotOutput(ns("MDelbowPlot"),width="auto",height="auto"),
                        textOutput(ns("MDbestSize"))
               ),
               tabPanel("Topological gene module detection for KEGG pathways - topological plot",
                        h3("KEGG topological plot(s) with highlighted module"),
                        uiOutput(ns("MDplots"))
               )
               )
             )
           )
  )
}
