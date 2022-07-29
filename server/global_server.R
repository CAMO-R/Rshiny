global_server <- function(input, output, session){
  ns <- NS("global")
  
  ##########################
  # Reactive Values        #
  ##########################
  DB <- reactiveValues(
    MergedDB=MergedDB.load(db), 
    MergedSpecies=MergedSpecies.load(db), 
    MergedStudyNames=MergedStudyNames.load(db), 
    pathway.list=NULL,
    ACS_ADS_global=NULL,
    ACS_ADS_pathway=NULL
  )
  ClustDEevid = reactiveValues(res=NULL,DEevid.pos=NULL,DEevid.neg=NULL)
  ##########################
  # Observers              #
  ##########################
  # tab change to load merged data
  observeEvent(input$tabChange, {
    try({
      DB$MergedDB <- MergedDB.load(db)
      DB$MergedSpecies <- MergedSpecies.load(db)
      DB$MergedStudyNames <- MergedStudyNames.load(db)
      print(DB$MergedStudyNames)
      
      if(length(DB$MergedDB)<2) {
        stop("At least two studies are needed")
      }
      # if(length(unique(DB$MergedSpecies))<2) {
      #   stop("At least two species are neeeded")
      # }
      output$globalACS_ADSTable <- DT::renderDataTable({
        if (!file.exists(paste(DB.load.working.dir(db),
                               "ACS_ADS_global.RData", sep="/"))){
          data.frame(NULL)
        }else{
          load(paste(DB.load.working.dir(db),"ACS_ADS_global.RData", sep="/"))
          DB$ACS_ADS_global <- ACS_ADS_global
          if(length(DB$ACS_ADS_global$ACS)==1){
            res <- matrix(1, length(DB$MergedDB), 
                          length(DB$MergedDB))
            colnames(res) <- rownames(res) <- DB$MergedStudyNames
            res[1,2] <- paste(round(DB$ACS_ADS_global$ACS,3), " (p-val= ", 
                              round(DB$ACS_ADS_global$ACSpvalue,3), ")", sep="")
            res[2,1] <- paste(round(DB$ACS_ADS_global$ADS,3), " (p-val= ", 
                              round(DB$ACS_ADS_global$ADSpvalue,3), ")", sep="")
            res
            print(res)
          }else{
            res <- DB$ACS_ADS_global$ACS
            for(i in 1:nrow(res)){
              for(j in 1:ncol(res)){
                if(j>=i){
                  res[i,j] <- paste(
                    round(DB$ACS_ADS_global$ACS[i,j],3), 
                    " (pval=", 
                    round(DB$ACS_ADS_global$ACSpvalue[i,j],3), ")", sep="")
                }
                else{
                  res[i,j] <- paste(
                    round(DB$ACS_ADS_global$ADS[i,j],3), 
                    " (pval=", 
                    round(DB$ACS_ADS_global$ADSpvalue[i,j],3), ")", sep="")
                }
              }
            }
            res
          }
        }
      })
      
    }, session)
    print(paste("saving directory is: ", DB.load.working.dir(db), sep=""))
  })
  
  # # select comparison type
  # observeEvent(input$compType, {
  #   try({
  #     if(input$compType != DB$compType){
  #       stop(paste("Please select ", DB$compType, sep=""))
  #     }
  #   }, session)
  # })
  
  # global ACS/ADS
  # output  global ACS/ADS
  output$Global_ACS_ADS_note <- renderText("Upper triangular matrix: genome-wide c-scores (p-value). 
                                           Lower triangular matrix: genome-wide d-scores (p-value).")
  observeEvent(input$ACS_ADS, {
    wait(session, "Calculating global c-scores & d-scores, may take a while")
    path_old <- getwd()
    try({
      path_old <- getwd()
      setwd(DB.load.working.dir(db))
      ACS_ADS_global <- multi_ACS_ADS_global(DB$MergedDB, DB$MergedStudyNames,
                                             measure=input$measure, B=input$permNumGlobal)
      save(ACS_ADS_global, file="ACS_ADS_global.RData")
      DB$ACS_ADS_global <- ACS_ADS_global
      print("Genome-wide c-scores and d-scores are saved as 'ACS_ADS_global'.")
      
      load("ACS_ADS_global.RData")
      DB$ACS_ADS_global <- ACS_ADS_global
      setwd(path_old)
      
      # ## global and pathway ACS/ADS
      # # single
      # if (DB$compType == "single"){
      #   path_old <- getwd()
      #   setwd(DB.load.working.dir(db))
      #   
      #   dat1 <- DB$MergedDB[[1]]
      #   dat2 <- DB$MergedDB[[2]]
      #   deIndex1 <- attr(DB$MergedDB[[1]], "DEindex")
      #   deIndex2 <- attr(DB$MergedDB[[2]], "DEindex")
      #   permOut <- perm_global(dat1,dat2,measure=input$measure,B=input$permNumGlobal)
      #   ACS_ADS_global$ACS <- ACS_global(dat1, dat2, deIndex1, deIndex2, measure=input$measure)
      #   print("global ACS finished")
      #   ACS_ADS_global$ACSpvalue <- pACS_global(dat1, dat2, deIndex1, deIndex2, 
      #                                           input$measure, DB$ACS_ADS_global$ACS, permOut)
      #   print("glocal pACS finished")
      #   ACS_ADS_global$ADS <- ADS_global(dat1, dat2, deIndex1, deIndex2, measure=input$measure)
      #   print("glocal ADS finished")
      #   ACS_ADS_global$ADSpvalue <- pADS_global(dat1, dat2, deIndex1, deIndex2,
      #                                           input$measure, DB$ACS_ADS_global$ADS, permOut)
      #   print("glocal pADS finished")
      #   DB$ACS_ADS_global <- ACS_ADS_global
      #   save(ACS_ADS_global, file="ACS_ADS_global.RData")
      #   print("global ACS/ADS saved")
      #   
      #   setwd(path_old)
      #   
      #   ## multiple comparison
      # }else if(DB$compType == "multiple"){
      #   
      #   ## multiple global 
      #   path_old <- getwd()
      #   setwd(DB.load.working.dir(db))
      #   ACS_ADS_global <- multi_ACS_ADS_global(DB$MergedDB, DB$MergedStudyNames,
      #                                          measure=input$measure, B=input$permNumGlobal)
      #   save(ACS_ADS_global, file="ACS_ADS_global.RData")
      #   DB$ACS_ADS_global <- ACS_ADS_global
      #   print("ACS_ADS_global saved")
      #   
      #   load("ACS_ADS_global.RData")
      #   DB$ACS_ADS_global <- ACS_ADS_global
      #   setwd(path_old)
      # }else{
      #   stop("At least one study is needed in each species")
      # }
      
      # output  global ACS/ADS
    output$globalACS_ADSTable <- DT::renderDataTable({
        if (is.null(DB$ACS_ADS_global$ACS)){
          data.frame(NULL)
        }else{
          if(length(DB$ACS_ADS_global$ACS)==1){
            res <- matrix(1, length(DB$MergedDB), 
                          length(DB$MergedDB))
            colnames(res) <- rownames(res) <- DB$MergedStudyNames
            res[1,2] <- paste(round(DB$ACS_ADS_global$ACS,3), " (p-val= ", 
                              round(DB$ACS_ADS_global$ACSpvalue,3), ")", sep="")
            res[2,1] <- paste(round(DB$ACS_ADS_global$ADS,3), " (p-val= ", 
                              round(DB$ACS_ADS_global$ADSpvalue,3), ")", sep="")
            res
            print(res)
          }else{
            res <- DB$ACS_ADS_global$ACS
            for(i in 1:nrow(res)){
              for(j in 1:ncol(res)){
                if(j>=i){
                  res[i,j] <- paste(
                    round(DB$ACS_ADS_global$ACS[i,j],3), 
                    " (pval=", 
                    round(DB$ACS_ADS_global$ACSpvalue[i,j],3), ")", sep="")
                }
                else{
                  res[i,j] <- paste(
                    round(DB$ACS_ADS_global$ADS[i,j],3), 
                    " (pval=", 
                    round(DB$ACS_ADS_global$ADSpvalue[i,j],3), ")", sep="")
                }
              }
            }
            res
          }
        }
      })
      sendSuccessMessage(session, "Genome-wide c-scores & d-scores are saved as 'ACS_ADS_global'.")
    }, session)
    setwd(path_old)
    done(session)
  })
  
  
  observeEvent(input$plotGlobalMDS, {
    wait(session, "Generating MDS using global c-scores")
    path_old <- getwd()
    try({
      library(utils)
      path_old <- getwd()
      setwd(DB.load.working.dir(db))
      load("ACS_ADS_global.RData")
      DB$ACS_ADS_global <- ACS_ADS_global
      res <- mdsGlobal(DB$ACS_ADS_global$ACS,DB$MergedStudyNames, sep="_",
                       file="globalMDS.pdf")
      print("Genome-wide MDS map is done")
      output$globalMdsFig <- renderPlot({
        res
      }, width = 600)
      setwd(path_old)
      sendSuccessMessage(session, "Genome-wide MDS map is done")
    }, session)
    setwd(path_old)
    done(session)
  })
  
  # # comparison type
  # output$compType = renderUI({
  #   selectInput(ns('compType'), 'Selelct comparison type:', 
  #               c("single", "multiple"),
  #               selected=as.character(DB$compType))
  # })
  
}

