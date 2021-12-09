setting_server <- function(input, output, session) {
  
  ns <- NS("Setting")
  
  ##########################
  # Reactive Values        #
  ##########################
  PACKAGES <- reactiveValues(installed=installed.packages()[,"Package"])
  
  installed <- function(pkg) {
    pkg %in% PACKAGES$installed
  }
  
  ##########################
  # Validation             #
  ##########################
  
  ##########################
  # Observers              #
  ##########################
  observeEvent(input$tabChange, {
    PACKAGES$installed <- installed.packages()[,"Package"]
  })
  
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$directory
    },
    handlerExpr = {
      if (input$directory > 0) {
        #path = choose.dir(default = readDirectoryInput(session, 'directory'))
        path = rstudioapi::selectDirectory(path = readDirectoryInput(session, 'directory'))
        print(path)
        if(length(path) == 0) {
          sendErrorMessage(session, MSG.no.working.dir)
        } else {
          updateDirectoryInput(session, 'directory', value = path)
          DB.set.working.dir(db, path)
          output$working.dir <- renderText({DB.load.working.dir(db)})
        }
      }
    }
  )
  
  ##########################
  # Render output/UI       #
  ##########################
  # working directory
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 DB.load.working.dir(db)
               },
               handlerExpr = {
                 output$working.dir <- renderText({
                   try({ DB.load.working.dir(db) }, session)
                 })
               })
  # output$working.dir <- renderText({
  #   try({ DB.load.working.dir(db) }, session)
  # })
  # session information
  output$urlText <- renderText({
    server.type <- ""
    if (session$clientData$url_hostname == "127.0.0.1")
      server.type <- "local"
    else
      server.type <- "remote"
    paste(sep = "",
          "protocol: ", session$clientData$url_protocol, "\n",
          "hostname: ", session$clientData$url_hostname, "\n",
          "port: ",     session$clientData$url_port,     "\n",
          "server type: ", server.type,     "\n"
    )
  })
  
}
