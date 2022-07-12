saved_data_server <- function(input, output, session) {
  
  ns <- NS("saved_data")
  
  DB <- reactiveValues(meta=meta(db), 
                       all_studies=DB.load(db, list.files(path=db@dir)))
  
  observeEvent(input$tabChange, {
    DB$all_studies <- DB.load(db, list.files(path=db@dir))
    DB$meta <- meta(db)  
    output$table_merge <- DT::renderDataTable(
      if(file.exists(paste(DB.load.working.dir(db), 
                           "MergedDB.rds", sep="/"))){
        DT::datatable({
          MergedSpecies=MergedSpecies.load(db) 
          MergedStudyNames=MergedStudyNames.load(db) 
          data.frame(MergedSpecies,MergedStudyNames)
        })
      }
    )
  })
  output$table <- DT::renderDataTable({
    DT::datatable({DB$meta})
  })


  # observeEvent(input$orthologous, {
  #   if (!is.null(input$orthologous)) {
  #     inFile <- input$orthologous
  #     DB$full_ortholog <- read.csv(inFile$datapath, stringsAsFactors = F,
  #                                  header=T)
  #   }
  # })
  
  
  observeEvent(c(input$select_orthologous,input$orthologous),{
    path_old <- getwd()
    
    try({
      if(input$select_orthologous == "no_orth"){
        DB$full_ortholog = NULL
      }else if(input$select_orthologous == "upload"){
        if(!is.null(input$orthologous)){
          file <- input$orthologous
          ext <- tools::file_ext(file$datapath)
          req(file)
          validate(need(ext == "RData"|ext == "rda", "Please upload a .RData or .rda file"))
          DB$full_ortholog = get(load(file$datapath))
        }
        # 
        # if(!is.null(input$orthologous)){
        #   inFile <- input$orthologous
        #   ext <- tools::file_ext(inFile$datapath)
        #   req(inFile)
        #   validate(need(ext == "csv", "Please upload a .csv"))
        #   DB$full_ortholog = read.csv(inFile$datapath, stringsAsFactors = F,
        #                               header=T)
        # }
      }else if(input$select_orthologous == "hs_mm_orth"){
        data(hs_mm_orth, package = "CAMO")
        DB$full_ortholog = hs_mm_orth
      }else if(input$select_orthologous == "hs_rn_orth"){
        data(hs_rn_orth, package = "CAMO")
        DB$full_ortholog = hs_rn_orth
      }else if(input$select_orthologous == "hs_ce_orth"){
        data(hs_ce_orth, package = "CAMO")
        DB$full_ortholog = hs_ce_orth
      }else if(input$select_orthologous == "hs_dm_orth"){
        data(hs_dm_orth, package = "CAMO")
        DB$full_ortholog = hs_dm_orth
      }else if(input$select_orthologous == "ce_dm_orth"){
        data(ce_dm_orth, package = "CAMO")
        DB$full_ortholog = ce_dm_orth
      }
      output$table_orth <- DT::renderDataTable(DT::datatable({
        DB$full_ortholog
      }))
    },session)
    setwd(path_old)
    done(session)
  })

  
  
  output$selected <- renderText({
    selected <- input$table_rows_selected
    if(length(selected) == 0){
      "You haven't select any study yet"
    }else{
      paste(DB$meta[selected,"studyName"], sep=", ")
      
    }
  })
  
  observeEvent(input$delete, {
    selected <- input$table_rows_selected
    if(length(selected) != 0 & !is.null(DB$meta)){
      selected <- rownames(DB$meta[selected,])
      DB.delete(db, selected)
      DB$all_studies <- DB.load(db, list.files(path=db@dir)) #DB$all_studies order and meta(db) order does not match
      studies <- DB.load(db, list.files(path=db@dir))
      db.meta  <- lapply(DB$all_studies, function(study_use) meta(study_use))
      DB$meta <- do.call(rbind, db.meta)
      sendSuccessMessage(session, paste(selected, "deleted"))
    }else{
      sendErrorMessage(session, "You haven't select any study yet")
    }
  })
  
  
  observeEvent(input$merge, {
    wait(session, "Match and merge")
    try({
      selected <- input$table_rows_selected
      print(selected)
      print(DB$meta)
      if(length(selected) != 0 & !is.null(DB$meta)){
        selected <- as.numeric(rownames(DB$meta[selected,]))
        species <- sapply(1:length(DB$all_studies), function(x) 
          DB$all_studies[[x]]@species)[selected]
        print(species)
        studyNames <- sapply(1:length(DB$all_studies), function(x) 
          DB$all_studies[[x]]@studyName)[selected]
        print(studyNames)
        
        mcmc.list <- lapply(1:length(DB$all_studies), function(x) 
          DB$all_studies[[x]]@MCMC)[selected]
        # if(is.null(DB$full_ortholog)){
        #   data(hs_mm_orth, package = "CAMO")
        #   DB$full_ortholog <- hs_mm_orth
        # }
        #print(paste("ref number :", which(species == input$reference)[1], sep=""))
        mcmc.merge.list <- CAMO::merge(mcmc.list, species = species,
                                       ortholog.db = DB$full_ortholog, 
                                       reference=which(species == input$reference)[1])
        names(mcmc.merge.list) = studyNames
        #print("finish merge")
        
        saveRDS(mcmc.merge.list,
                file=paste(DB.load.working.dir(db), 
                           "MergedDB.rds", sep="/"))
        saveRDS(species, 
                file=paste(DB.load.working.dir(db), 
                           "MergedSpecies.rds", sep="/"))
        saveRDS(studyNames, 
                file=paste(DB.load.working.dir(db), 
                           "MergedStudyNames.rds", sep="/"))
        
        output$table_merge <- DT::renderDataTable(      
          DT::datatable({
          MergedSpecies=MergedSpecies.load(db) 
          MergedStudyNames=MergedStudyNames.load(db) 
          data.frame(MergedSpecies,MergedStudyNames)
        }))
        
        message = paste("Data are successfully merged")
        sendSuccessMessage(session, message)
      }else{
        sendErrorMessage(session, "You haven't select any study yet")
      }
    },session)
    done(session)
  }, label="save study")
  

  ##########################
  # Render output/UI       #
  ##########################
  output$reference = renderUI({
    if(length(DB$all_studies) != 0 & input$select_orthologous != "no_orth"){
      selectInput(ns('reference'), 'Reference species', 
                  as.character(unique(sapply(1:length(DB$all_studies), function(x) 
                    DB$all_studies[[x]]@species))),
                  selected=as.character(unique(sapply(1:length(DB$all_studies), 
                                                      function(x) DB$all_studies[[x]]@species))) [1])
      
    }
  })
  
  
  output$upload_orthologous = renderUI({
    if(input$select_orthologous == "upload"){
      fileInput(ns("orthologous"), "Upload an ortholog file (.RData/.rda)", accept = c(".RData",".rda"))
    }
  })
}
