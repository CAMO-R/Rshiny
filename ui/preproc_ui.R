preproc_ui <- function(id, label= "preprocessing data") {
  ns <- NS(id)

  tabPanel("Data Uploading and Preprocessing", value=id,
    sidebarLayout(
      sidebarPanel(
        useShinyjs(),
        ##########################
        # Upload Data 
        ##########################
        h2("Upload data"),
        #### chose species
        selectInput(ns('select_species'), 'Select a species', 
                    c("homo sapiens (hs)"="hs","mus musculus (mm)"="mm",
                      "rattus norvegicus (rn)"="rn","caenorhabditis elegans (ce)"="ce",
                      "drosophila melanogaster (dm)"="dm","other"="other")),
        uiOutput(ns("other_species")),
        # textInput(ns("species"), label="Species name:", value = "", width = NULL,
        #           placeholder = NULL),
        # radioButtons(ns("species"), label = "Chose species",
        #         choices = list("human" = "human", "mouse" = "mouse"),
        #         selected = "human"),
        #### chose input data type
        radioButtons(ns("inputType"), label = "Input data type",
                choices = list("Pre-calculated p-value and log fold change (logFC) from gene-level differential expression analysis" = 1, 
                               "Preprocessed gene expression data" = 2),
                selected = 2),

        #### if input data are p-values
        conditionalPanel(
          condition = paste0("input['", ns("inputType"), "'] == '1' "),
          fileInput(ns("pvalfile"), 'Upload p-value and logFC file (.csv)',
            accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
          ),

        #### if input data are raw data perform single DE first
        conditionalPanel(
          condition = paste0("input['", ns("inputType"), "'] == '2' "),
          fileInput(ns("exprfile"), 'Upload gene expression data file (.csv)',
            accept=c('text/csv', 'text/comma-separated-values,text/plain', 
              '.csv')
          ),

          fileInput(ns("clinical"), 'Upload clinical file (.csv)',
            accept=c('text/csv', 'text/comma-separated-values,text/plain', 
              '.csv')
          ),
          uiOutput(ns('caseName')),
          radioButtons(ns("dtype"), label = "Differential analysis method",
                       choices = list("LIMMA (for microarray intensities  / log2-transformed RNA-Seq normalized counts)" = "microarray", 
                                      "DEseq2 (for RNA-Seq counts)" = "RNAseq"),
                       selected = "microarray")
        ),
        # #### calculate signed delta and q-values
        actionButton(ns('SingleDE'), 'Bayesian differential analysis', class="btn-success",
                     icon = icon("play")),
        
        # h2("Preprocess data"),
        # br(),
        # actionButton(ns('preprocSingleStudy'), 'Preprocess single study', class="btn-success"),
        # tags$hr(),

        ##########################
        # Save and Metadata      #
        ##########################
        h2("Save single study"),
        textInput(ns("studyName"), "Study name (Please do not use '_' in study name):", value = ""),
        actionButton(ns('saveStudy'), 'Save', icon=icon("save"), class="btn-success")
      ),

      mainPanel(
        h3("Study summary"),
        DT::dataTableOutput(ns("studySummary")),
        hr(),
        h3("Bayesian differential analysis summary"),
        textOutput(ns('DEdescription')),
        br(),
        DT::dataTableOutput(ns("DESummary"))
      )
    )
  )
}
