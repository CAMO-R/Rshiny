saved_data_ui <- function(id, label = "saved data of single study or multiple study") {
  ns <- NS(id)
  tabPanel("Saved Data", value=id,
    sidebarLayout(
      sidebarPanel(
        p("Selected studies:"), helpIcon(ns('merge_select'), "Click on the rows of the list of saved studies to select studies to merge/delete"),
        textOutput(ns("selected")),
        br(),
        selectInput(ns('select_orthologous'), 'Cross-species ortholog matching file', 
                    c("homo sapiens (hs) vs mus musculus (mm)"="hs_mm_orth",
                      "homo sapiens (hs) vs rattus norvegicus (rn)"="hs_rs_orth",
                      "homo sapiens (hs) vs caenorhabditis elegans (ce)"="hs_ce_orth",
                      "homo sapiens (hs) vs drosophila melanogaster (dm)"="hs_dm_orth",
                      "caenorhabditis elegans (ce) vs drosophila melanogaster (dm)"="ce_dm_orth",
                      "Upload an ortholog file"="upload"),
                    selected = NULL),
        uiOutput(ns("upload_orthologous")),
        # fileInput(ns("orthologous"), 'Upload orthologous file',
        #   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
        uiOutput(ns('reference')),
        actionButton(ns("merge"), 'Match and merge', class="btn-success", icon=icon("rocket")),
        hr(),
        hr(),
        #p("Selected datasets"), helpIcon(ns('merge_select'), HELP.select.datasets),
        actionButton(ns('delete'), 'Delete selected studies', icon=icon("trash"), 
                     class="btn-danger")
      ),
      mainPanel(
        h3("List of saved studies"),
        DT::dataTableOutput(ns("table")),
        hr(),
        h3("Ortholog matching file selected"),
        DT::dataTableOutput(ns("table_orth")),
        hr(),
        h3("List of merged studies"),
        DT::dataTableOutput(ns("table_merge"))
        
      )
    )
  )
}
