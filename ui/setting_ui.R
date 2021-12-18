setting_ui <- function(id, label = "global settings") {
  ns <- NS(id)
  
  tabPanel("Setting", value=id,
           h2("Welcome to CAMO",align = "middle",style="primary"),
           tags$hr(),
           #HTML('<center><img src="Rshiny_flowchart.png" width="500" </center>'),
           fluidRow(
             img(src='Shiny_flowchart.jpg',align="center",width="1500"),
             p("CAMO is an analytical software with R-Shiny based graphical user interface (GUI) for evaluating omics congruence of model organisms. It performs threshold-free Bayesian differential analysis and generates quantitative concordance and discordance scores (c-scores and d-scores) both genome-wide and at pathway level for all pair-wise studies. Based on the c-scores/d-scores, CAMO conducts a series of downstream machine learning and bioinformatics analysis with interactive visualization for pathway knowledge retrieval and topological gene module detection. Outputs from the tool will provide foundations for hypothesis generation and subsequent translational investigations."), 
             p("Our tool is available for download on github: ",a(strong("CAMO."), href="https://github.com/weiiizong/CAMO_Rshiny",target="_blank"), 
             "For detailed implementation of the tool, please refer to our ",a(strong("Tutorials."), href="https://github.com/metaOmicsCAMO/tutorial/
blob/master/CAMO_turtorial.pdf",target="_blank")), 
             p("CAMO is developed and maintained by ", a("Dr. George Tseng's group ",href="http://tsenglab.biostat.pitt.edu",target="_blank"),"(Department of Biostatistics, University of Pittsburgh) and ", a("Dr. Tianzhou Ma's group ",href="https://matianzhou.github.io",target="_blank"), "(Department of Epidemiology and Biostatistics, University of Maryland College Park)."),
             p("We recommend users to use R (>=3.5.0) to implement our tool. If you are using R 3.4, 
               you may encounter errors in installing dependencies of the modules. "),
             style="text-indent: 20px; font-size: 16px"),
           tags$hr(),
    mainPanel(
      h2("Session Information"),
      verbatimTextOutput(ns("urlText")),
      helpIcon("working_dir_help", "During the computation, some output files or images are automatically saved to this directory."),
      h2("Saving directory:", style="display:inline"),
      directoryInput(ns('directory'), label='select a directory')
    )
  )
}
