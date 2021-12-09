global_ui <- function(id, label = "Global") {
  ns <- NS(id)
  tabPanel("Genome-wide analysis", value=id,
           sidebarLayout(
             sidebarPanel(
               h2("Genome-wide c-scores & d-scores calculation"),
               #uiOutput(ns("compType")),
               ## Global ARS/ADS
               selectInput(ns("measure"), label = "Score type", 
                           choices = list("F-measure (recommended)" = "Fmeasure",
                                          "Youden index" = "youden", 
                                          "Geometric mean" = "geo.mean")),
               
               numericInput(ns("permNumGlobal"), 
                            label = "Number of permutations", 
                            value = 100),

               actionButton(ns('ACS_ADS'), 'Genome-wide c-scores & d-scores', 
                            class="btn-success",icon = icon("play")), 
               
               hr(),
               actionButton(ns("plotGlobalMDS"),
                            'Genome-wide MDS map',
                            class="btn-success",
                            icon = icon("play")),   
               tags$hr()
             ),
             mainPanel(
               tabsetPanel(tabPanel("Table of genome-wide c-scores and d-scores",
                                    h3("Genome-wide c-scores & d-scores"),
                                    textOutput(ns("Global_ACS_ADS_note")),
                                    DT::dataTableOutput(ns("globalACS_ADSTable"))),
                           tabPanel("Multidimensional scaling (MDS) map",
                                    h3("MDS map of all studies based on genome-wide c-scores"),
                                    plotOutput(ns("globalMdsFig"), height = 500))
               )
               
               #plotOutput(ns("ARS_DE"), height = 500, click = clickOpts("plot_click"),
               #hover = hoverOpts("plot_hover")),
               #verbatimTextOutput(ns("click_info")),
               #verbatimTextOutput(ns("hover_info")),
               #h3("KEGG pathway viewer"),
               #imageOutput(ns("KEGGtopo"), height = 500)))
               
             )
           )
  )
}