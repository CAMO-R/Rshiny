preproc_server <- function(input, output, session) {

  ns <- NS("preproc")
  ##########################
  # Reactive Values        #
  ##########################
  DB   <- reactiveValues(names=DB.ls(db))
  STUDY <- reactiveValues(action="", 
    update=0, ori=NULL, preview=NULL,  
    expr=NULL, clinical=NULL, 
    studyName=NULL, species=NULL, genes=NULL,
    DE_p=NULL, DE_lfc=NULL, MCMC=NULL)
  SUMMARY <- reactiveValues(studySummary=data.frame(NULL), 
    DESummary=data.frame(NULL), DESummaryPM=data.frame(NULL))

  ##########################
  # Validation             #
  ##########################
  validate.study <- function(study) {
    if(is.null(STUDY$MCMC))  {
      stop(MSG.datasetInput.noinput)
    }
    studyName <- STUDY$studyName
    if(is.null(studyName) || studyName == "") {
      stop(MSG.study.noname)
    }
    if(studyName %in% DB$names) stop(MSG.study.duplicate(studyName))  
  }

  ##########################
  # Observers              #
  ##########################
  # watch for tab change, get the newest list of all data
  observeEvent(input$tabChange, {DB$names <- DB.ls(db)}, 
               label="tab change")
  
  
  species0 = eventReactive(c(input$select_species,input$species),{
    ifelse(input$select_species == "other",input$species,input$select_species)
  })
  # watch for pvalue files upload
  observeEvent(input$pvalfile, {
    if (!is.null(input$pvalfile)) {
      STUDY$action <- STUDY.pval.upload
      STUDY$update <- STUDY$update + 1
      STUDY$species <- species0()
    }
  }, label="p-value file upload")

  # watch for expression file upload
  observeEvent(input$exprfile, {
    if (!is.null(input$exprfile)) {
      STUDY$action <- STUDY.expression.upload
      STUDY$update <- STUDY$update + 1
      STUDY$species <- species0()
    }
  }, label="expression file upload")

  # watch for clinical file upload
  observeEvent(input$clinical, {
    if (!is.null(input$clinical)) {
      STUDY$action <- STUDY.clinical.upload
      STUDY$update <- STUDY$update + 1
    }
  }, label="clinical file upload")

  # Setting study from file type
  observeEvent(STUDY$update, {
    if (STUDY$action == STUDY.pval.upload) {
      wait(session, "Parsing pvalue File...")
      inFile <- input$pvalfile
      tryCatch({
        pvalfile <- read.pData(inFile$datapath)
        STUDY$genes <- rownames(pvalfile)
        STUDY$DE_p <- pvalfile[,"pvalue"]
        STUDY$DE_lfc <- pvalfile[,"logFC"]
        studySummary <- data.frame(species=STUDY$species, 
          NumGenes=length(STUDY$DE_p))
        SUMMARY$studySummary = studySummary
        SUMMARY$DESummary <- cbind(STUDY$DE_p, STUDY$DE_lfc)
        rownames(SUMMARY$DESummary) <- STUDY$genes
        colnames(SUMMARY$DESummary) <- c("pvalue", "logFC")
      }, error=function(error) {
        sendErrorMessage(session, MSG.file.corrupted)
      })
      done(session)
    } else if (STUDY$action == STUDY.expression.upload) {
      wait(session, "Parsing Expression File...")
      inFile <- input$exprfile
      tryCatch({
        STUDY$expr <- read.rawData(inFile$datapath)
        STUDY$genes <- rownames(STUDY$expr)
      }, error=function(error) {
        sendErrorMessage(session, MSG.file.corrupted)
      })
      done(session)
    }  else if (STUDY$action == STUDY.clinical.upload) {
      wait(session, "Parsing Clinical File...")
      inFile <- input$clinical
      tryCatch({
        STUDY$clinical <- read.groupData(inFile$datapath)
        group_labels <- sort(unique(STUDY$clinical))
        studySummary <- data.frame(species=STUDY$species, 
          NumSamples=ncol(STUDY$expr), NumGenes=nrow(STUDY$expr), 
          NumGroup1=sum(STUDY$clinical==group_labels[1]),
          NumGroup2=sum(STUDY$clinical==group_labels[2]))
        colnames(studySummary)[4] <- paste("NumGroup_", group_labels[1], 
          sep="")
        colnames(studySummary)[5] <- paste("NumGroup_", group_labels[2], 
          sep="")
        SUMMARY$studySummary = studySummary
      }, error=function(error) {
        sendErrorMessage(session, MSG.file.corrupted)
      })
      done(session)
    }
  }, label="setting STUDY$ori from file upload or selection")
  
  # DE analysis
  observeEvent(input$SingleDE, {
    if(is.null(SUMMARY$DESummary)){
      DEmethod = ifelse(input$dtype == "microarray", "Limma", "DEseq2")
      wait(session, paste0("Running DE analysis using ",DEmethod,"..."))
      try({
        SingleDERes <- indDE(data=STUDY$expr, 
                             group=as.factor(STUDY$clinical), 
                             data.type=input$dtype, case.label=input$caseName, 
                             ctrl.label=setdiff(unique(STUDY$clinical), input$caseName))
        DE_lfc <- SingleDERes[,"logFC"]
        DE_p <- SingleDERes[,"pvalue"]
        STUDY$DE_p <- DE_p
        STUDY$DE_lfc <- DE_lfc
        SUMMARY$DESummary <- cbind(STUDY$DE_p, STUDY$DE_lfc)
        rownames(SUMMARY$DESummary) <- STUDY$genes
        colnames(SUMMARY$DESummary) <- c("pvalue", "logFC")
      }, session)
    }
    
    wait(session, "Running threshold-free Bayesian differential analysis to generate posterior DE probabilities...")
    try({
      STUDY$MCMC <- bayes(SUMMARY$DESummary, seed=12345)
      signPM.vec <- apply(STUDY$MCMC,1,mean)
      SUMMARY$DESummaryPM <- cbind(SUMMARY$DESummary, signPM.vec)
      rownames(SUMMARY$DESummaryPM) <- STUDY$genes
      colnames(SUMMARY$DESummaryPM) <- c("pvalue", "logFC","Posterior DE probability")
    }, session)
  
    done(session)
  })

  # # preprocess single studies
  # observeEvent(input$preprocSingleStudy, {
  #   wait(session, "Preprocessing, may take a while")
  #   try({
  #     # cat("test")
  #     # bayes <- function(pData, seed=12345){
  #     #   
  #     #   set.seed(seed)
  #     #   
  #     #   ## deSelect part
  #     #   
  #     #   DEindex <- deSelect(pData)
  #     #   
  #     #   ## bayesP part   
  #     #   p <- pData[,1]
  #     #   lfc <- pData[,2]
  #     #   G <- nrow(pData)
  #     #   z <- PtoZ(p,lfc)
  #     #   iteration <- 5000
  #     #   burnin <- 3000
  #     #   thin <- 10
  #     #   names(z) <- rownames(pData)
  #     #   prop <- SelectGamma(p)
  #     #   if(prop <= 0.3) {
  #     #     gamma <- G*0.3
  #     #   } else{
  #     #     gamma <- G*prop
  #     #   }
  #     #   MCMCout <- MCMC(z, iteration, gamma) 
  #     #   signdelta <- MCMCout$Y[,-c(1:burnin)]
  #     #   
  #     #   signdelta.sub <- signdelta[,seq(1,ncol(signdelta),by=thin)]
  #     #   rownames(signdelta.sub) <- names(z)
  #     #   
  #     #   attr(signdelta.sub,"DEindex") <- DEindex
  #     #   
  #     #   return(signdelta.sub) #the full sign delta for a dataset (subsample 500)
  #     # }  
  #     
  #     STUDY$MCMC <- bayes(SUMMARY$DESummary, seed=12345)
  #     signPM.vec <- apply(STUDY$MCMC,1,mean)
  #     SUMMARY$DESummaryPM <- cbind(SUMMARY$DESummary, signPM.vec)
  #     rownames(SUMMARY$DESummaryPM) <- STUDY$genes
  #     colnames(SUMMARY$DESummaryPM) <- c("pvalue", "logFC","Bayes posterior DE probability")
  #   }, session)
  #   done(session)
  # })

  # Save and Metadata
  observeEvent(input$saveStudy, {
    wait(session, "Saving study")
    try({
      if(length(grep("_",input$studyName)) != 0){
        stop("Please do not use '_' in study name")
      }
      STUDY$studyName <- input$studyName

      validate.study(STUDY)
      study_use <- new("Study1",
        studyName=STUDY$studyName,
        species=STUDY$species,
        MCMC=STUDY$MCMC)
      
      DB.save(db, study_use)
      sendSuccessMessage(session, paste("Study", study_use@studyName, "saved."))
      DB$names <- DB.ls(db)
      reset("species")
      reset("inputType")
      reset("pvalfile")
      reset("exprfile")
      reset("clinical")
      reset("studyName")
      reset("PCut")
      reset("DENum")
      SUMMARY$studySummary <- data.frame(NULL)
      SUMMARY$DESummary <- data.frame(NULL)
      SUMMARY$DESummaryPM <- data.frame(NULL)
    }, session)
    done(session)
  }, label="save study")

  ##########################
  # Render output/UI       #
  ##########################
  output$caseName = renderUI({
      selectInput(ns('caseName'), 'Case name', 
        as.character(sort(unique(STUDY$clinical))),
       selected=as.character(sort(unique(STUDY$clinical)))[2])
  })
  
  output$other_species = renderUI({
    if(input$select_species == "other"){
      textInput(ns("species"), label="Specify the species name:", value = "", width = NULL,
                placeholder = NULL)
    }
  })

  # studySummary
  output$studySummary <- DT::renderDataTable({
    SUMMARY$studySummary
  })

  # DESummary
  output$DESummary <- DT::renderDataTable({
    if (nrow(SUMMARY$DESummaryPM) != 0){
      SUMMARY$DESummaryPM
    }else{
      SUMMARY$DESummary
    }
  })
  
  #DEdescription
  DEmethod <- eventReactive(input$dtype,{
    ifelse(input$dtype == "microarray", "Limma", "DEseq2")
  })
  output$DEdescription <- renderText({
    if (nrow(SUMMARY$DESummaryPM) != 0 | nrow(SUMMARY$DESummary) != 0 ){
      paste0("P-values and logFC are generated from the classical differential analysis (DE) using ", DEmethod(),". Posterior probabilities of DE assignment are generated from the Bayesian Gaussian mixture model.")
    }
  })
}
