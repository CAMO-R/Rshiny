
save_result_ui <- function(id, label = "Save results") {
  ns <- NS(id)
  tabPanel("Save results", value=id,
           sidebarLayout(
             sidebarPanel(
               h3("Generate and save Kegg plots for study pair"),
               uiOutput(ns('keggStudy1')),
               uiOutput(ns('keggStudy2')),
               actionButton(ns("KEGG"),"KEGG plots",class = "btn-success"),
               uiOutput(ns('selectPathwayName')),
               hr(),
               h3("Generate and save visulizations"),
               uiOutput(ns('selectVisualizations')),
               actionButton(ns('plotAll'), 'Generate visualizations', 
                            class="btn-success"),
               tags$hr()
             ),
             mainPanel()
           )
  )
}
