indPathway_server <- function(input, output, session){
  ns <- NS("indPathway")
  
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
    #keggSpecies="", keggGenesfile=NULL, 
    #reactomeSpecies="", reactomeGenesfile=NULL,
    #KEGGpathwayID=NULL, gene_type=NULL, 
    #plotall=NULL,
    #MD_data=NULL, MD_size=NULL,
    #MD_whichToDraw=NULL
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
      if(length(DB$MergedDB)<2) {
        stop("At least two studies are needed")
      }
      # if(length(unique(DB$MergedSpecies))<2) {
      #   stop("At least two species are needed")
      # }
      # if(length(DB$MergedDB)==2 & length(unique(DB$MergedSpecies))==2){
      #   DB$compType <- "single"
      # }else{
      #   DB$compType <- "multiple"
      # }
      print(DB$MergedStudyNames)
      
      if(file.exists(paste0(DB.load.working.dir(db),"/ACS_ADS_pathway.RData"))){
        load(paste0(DB.load.working.dir(db),"/ACS_ADS_pathway.RData"))
        DB$ACS_ADS_pathway = ACS_ADS_pathway
      }
      if(file.exists(paste0(DB.load.working.dir(db),"/select_pathway.RData"))){
        load(paste0(DB.load.working.dir(db),"/select_pathway.RData"))
        DB$pathway.list = pathway.list
      }
      print(paste("saving directory is: ", DB.load.working.dir(db), sep=""))
    }, session)
  })
  
  # Generate all visualization plots
  observe({
    if(is.null(input$selectVisualizations)) {
      output$keggViewSelect = renderUI({NULL})
      output$reactomeViewSelect = renderUI({NULL})
    }
    if ("keggView" %in% input$selectVisualizations) {
      output$keggViewSelect = renderUI({
        bsCollapsePanel("Settings for topological visualization",
                        hr(),
                        textInput(ns("keggSpecies"), "KEGG organisms code", "hsa"),
                        radioButtons(ns("KEGG.dataGisEntrezID"), "Gene names in data matrix are Entrez IDs?", choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "FALSE"),
                        fileInput(ns("keggGenesfile"), "Upload a data frame mapping (merged) gene names in data [column1] to 
Entrez IDs [column2] (.RData/.rda). If not provided & choose FALSE in the question above & KEGG organisms code is one of 'hsa', 'mmu','rno, 'cel' or 'dme', gene symbols will be automatically mapped to Entrez IDs by Bioconductor packages
                                    'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db', 'org.Ce.eg.db' or 'org.Dm.eg.db'", accept = c(".RData","rda")),
                        fileInput(ns("keggIDfile"), "Upload a list mapping pathway names [content] to KEGG IDs [name] (.RData/.rda). If not provided, ID will be retrieved from KEGGREST.", accept = c(".RData","rda")),
                        checkboxGroupInput(ns('keggViewSelected'), 'Select a subset of studies to generalize KEGG topology plots',
                                           DB$MergedStudyNames, selected = DB$MergedStudyNames))
      })
    } else {
      output$keggViewSelect = renderUI({NULL})
    }
    if ("reactomeView" %in% input$selectVisualizations) {
      output$reactomeViewSelect = renderUI({
        bsCollapsePanel("Settings for topological visualization",
                        hr(),
                        textInput(ns("reactomeSpecies"), "Reactome organisms code", "HSA"),
                        fileInput(ns("reactomeGenesfile"), "Upload a data frame mapping (merged) gene names in data [column1] to to gene names in Reactome topology [column2] (.RData/.rda). 
                                  If not provided, it assumes the gene names in data matrix and topology plots are the same type.", accept = c(".RData","rda")),
                        fileInput(ns("reactomeIDfile"), "Upload a list mapping pathway names [content] to Reactome IDs [name] (.RData/.rda). If not provided, ID will be retrieved from reactome.db", accept = c(".RData","rda")),
                        checkboxGroupInput(ns('keggViewSelected'), 'Select a subset of studies to generalize Reactome topology plots',
                                           DB$MergedStudyNames, selected = DB$MergedStudyNames)
        )
      })
    } else {
      output$reactomeViewSelect = renderUI({NULL})
    }
    if(("keggView" %in% input$selectVisualizations) && ("reactomeView" %in% input$selectVisualizations)) {
      output$keggViewSelect = renderUI({
        bsCollapsePanel("Settings for topological visualization",
                        hr(),
                        textInput(ns("keggSpecies"), "KEGG organism code", "hsa"),
                        radioButtons(ns("KEGG.dataGisEntrezID"), "Gene names in data matrix are Entrez IDs?", choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "FALSE"),
                        fileInput(ns("keggGenesfile"), "Upload a data frame mapping (merged) gene names in data [column1] to 
Entrez IDs [column2] (.RData/.rda). If not provided & choose FALSE in the question above & KEGG organism code is one of 'hsa', 'mmu','rno, 'cel' or 'dme', gene symbols will be automatically mapped to Entrez IDs by Bioconductor packages
                                    'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db', 'org.Ce.eg.db' or 'org.Dm.eg.db'", accept = c(".RData","rda")),
                        fileInput(ns("keggIDfile"), "Upload a list mapping pathway names [content] to KEGG IDs [name] (.RData/.rda). If not provided, ID will be retrieved from KEGGREST.", accept = c(".RData","rda")),
                        hr(),
                        textInput(ns("reactomeSpecies"), "Reactome organism code", "HSA"),
                        fileInput(ns("reactomeGenesfile"), "Upload a data frame mapping (merged) gene names in data [column1] to to gene names in Reactome topology [column2] (.RData/.rda). 
                                  If not provided, it assumes the gene names in data matrix and topology plots are the same type.", accept = c(".RData","rda")),
                        fileInput(ns("reactomeIDfile"), "Upload a list mapping pathway names [content] to Reactome IDs [name] (.RData/.rda). If not provided, ID will be retrieved from reactome.db", accept = c(".RData","rda")),
                        checkboxGroupInput(ns('keggViewSelected'), 'Select a subset of studies to generalize KEGG/Reactome topology plots',
                                           DB$MergedStudyNames, selected = DB$MergedStudyNames))
      })
      output$reactomeViewSelect = renderUI({NULL})
    }
  })
  observeEvent(input$plotAll, {
    wait(session, "Generating individual pathway plots, may take a while...")
    path_old = getwd()
    try({ 
      path_old = getwd()
      if(file.exists(paste0(DB.load.working.dir(db),"/ACS_ADS_pathway.RData"))){
        load(paste0(DB.load.working.dir(db),"/ACS_ADS_pathway.RData"))
        DB$ACS_ADS_pathway = ACS_ADS_pathway
      }
      if(file.exists(paste0(DB.load.working.dir(db),"/select_pathway.RData"))){
        load(paste0(DB.load.working.dir(db),"/select_pathway.RData"))
        DB$pathway.list = pathway.list
      }
      
      ### keggView
      if(!is.null(input$keggGenesfile)){
        file <- input$keggGenesfile
        ext <- tools::file_ext(file$datapath)
        req(file)
        validate(need(ext == "RData"|ext == "rda", "Please upload a .RData or .rda file"))
        keggGenesfile = get(load(file$datapath))
      }
      else{
        keggGenesfile=NULL
      }
      if(!is.null(input$keggIDfile)){
        file <- input$keggIDfile
        ext <- tools::file_ext(file$datapath)
        req(file)
        validate(need(ext == "RData"|ext == "rda", "Please upload a .RData or .rda file"))
        keggIDfile = get(load(file$datapath))
      }else{
        keggIDfile=NULL
      }
      #DB$plotall = input$selectVisualization
      #DB$keggGenesfile=keggGenesfile
      ### reactomeView
      if(!is.null(input$reactomeGenesfile)){
        file <- input$reactomeGenesfile
        ext <- tools::file_ext(file$datapath)
        req(file)
        validate(need(ext == "RData"|ext == "rda", "Please upload a .RData or .rda file"))
        reactomeGenesfile = get(load(file$datapath))
      }
      else{
        reactomeGenesfile=NULL
      }
      if(!is.null(input$reactomeIDfile)){
        file <- input$reactomeIDfile
        ext <- tools::file_ext(file$datapath)
        req(file)
        validate(need(ext == "RData"|ext == "rda", "Please upload a .RData or .rda file"))
        reactomeIDfile = get(load(file$datapath))
      }else{
        reactomeIDfile=NULL
      }
      #DB$reactomeSpecies = input$reactomeSpecies
      #DB$reactomeGenesfile=reactomeGenesfile
      # load("res_clustPath.RData")
      # res$CluterLabelwithoutScatter
      keggViewSelect <- DB$MergedStudyNames[grep(paste(input$keggViewSelected,collapse="|"), 
                             DB$MergedStudyNames, value=FALSE)]
      print(keggViewSelect)
      
      setwd(DB.load.working.dir(db))

      multiOutput(mcmc.merge.list=DB$MergedDB, 
                  dataset.names=DB$MergedStudyNames, 
                  select.pathway.list=DB$pathway.list,
                  DB$ACS_ADS_pathway,
                  output=input$selectVisualizations,
                  use_ADS = FALSE,
                  ViewPairSelect = keggViewSelect, 
                  kegg.species = input$keggSpecies, KEGG.dataGisEntrezID=ifelse(input$KEGG.dataGisEntrezID=="TRUE",TRUE,FALSE),
                  KEGG.dataG2EntrezID = keggGenesfile,KEGG.pathID2name=keggIDfile,
                  reactome.species=input$reactomeSpecies, Reactome.dataG2TopologyGtype= reactomeGenesfile,Reactome.pathID2name=reactomeIDfile)
      setwd(path_old)
      message = paste("Individual pathway plots saved.")
      sendSuccessMessage(session, message)
    }, session) 
    setwd(path_old)
    done(session)
  })
  
  
  # visulize ACS/ADS for selected individual pathway
  observeEvent(input$selectedPathway, {
    if(input$selectedPathway!=""){
      wait(session, "Generating ACS/ADS table for slected pathway...")
      path_old = getwd()
      try({
        path_old = getwd()
        if(file.exists(paste0(DB.load.working.dir(db),"/ACS_ADS_pathway.RData"))){
          load(paste0(DB.load.working.dir(db),"/ACS_ADS_pathway.RData"))
          DB$ACS_ADS_pathway = ACS_ADS_pathway
        }
        if(file.exists(paste0(DB.load.working.dir(db),"/ACS_ADS_global.RData"))){
          load(paste0(DB.load.working.dir(db),"/ACS_ADS_global.RData"))
          DB$ACS_ADS_global = ACS_ADS_global
        }
        selectedPathwayACS <- DB$ACS_ADS_pathway$ACS.mat[input$selectedPathway,,drop=F]
        selectedPathwaypACS <- DB$ACS_ADS_pathway$ACSpvalue.mat[input$selectedPathway,,drop=F]
        selectedPathwayADS <- DB$ACS_ADS_pathway$ADS.mat[input$selectedPathway,,drop=F]
        selectedPathwaypADS <- DB$ACS_ADS_pathway$ADSpvalue.mat[input$selectedPathway,,drop=F]
        res <- DB$ACS_ADS_global$ACS
        c = 1
        d = 1
        for(i in 1:(nrow(DB$ACS_ADS_global$ACS))){
          for(j in 1:ncol(DB$ACS_ADS_global$ACS)){
            if(j>i){
              res[i,j] <- paste0(round(selectedPathwayACS[1,c],3), " (pval=",
                                 round(selectedPathwaypACS[1,c],3), ")")
              c <- c+1
            }
            else if(j<i){
              res[i,j] <- paste0(round(selectedPathwayADS[1,d],3), " (pval=",
                                 round(selectedPathwaypADS[1,d],3), ")")
              d <- d+1
            }
          }
        }
        setwd(path_old)
        output$Show_selected_pathway_name <- renderText(paste0("Pathway name: ", row.names(selectedPathwayACS)))
        output$Pathway_ACS_ADS_note <- renderText("Upper triangular matrix: pathway c-scores (p-value). 
                                                  Lower triangular matrix: pathway d-scores (p-value).")
        output$selectedPathwayACS_ADS_Table <- DT::renderDataTable({
          if (!(is.null(res))) {
            res
          }
        })
        
      }, session)
      setwd(path_old)
      }
    done(session)
  })
  
  #Individual visualization
  observeEvent(input$showVisual, {
    path_old = getwd()
    try({
      if(!is.null(input$chosenPathway)){
        selectPathwayName = input$chosenPathway
        
        path_old = getwd()
        if(file.exists(paste0(DB.load.working.dir(db),"/ACS_ADS_pathway.RData"))){
          load(paste0(DB.load.working.dir(db),"/ACS_ADS_pathway.RData"))
          DB$ACS_ADS_pathway = ACS_ADS_pathway
        }
        if(file.exists(paste0(DB.load.working.dir(db),"/select_pathway.RData"))){
          load(paste0(DB.load.working.dir(db),"/select_pathway.RData"))
          DB$pathway.list = pathway.list
        }

        selectedPathwayAS = DB$ACS_ADS_pathway$ACS.mat[selectPathwayName,]
        selectedPathwaypAS = DB$ACS_ADS_pathway$ACSpvalue.mat[selectPathwayName,]
        
        # if(input$browser_scoreType == "c_scores"){
        #   selectedPathwayAS = DB$ACS_ADS_pathway$ACS.mat[selectPathwayName,]
        #   selectedPathwaypAS = DB$ACS_ADS_pathway$ACSpvalue.mat[selectPathwayName,]
        # }else{
        #   selectedPathwayAS = DB$ACS_ADS_pathway$ADS.mat[selectPathwayName,]
        #   selectedPathwaypAS = DB$ACS_ADS_pathway$ADSpvalue.mat[selectPathwayName,]
        # }
        
        if(!is.null(input$keggGenesfile)){
          file <- input$keggGenesfile
          ext <- tools::file_ext(file$datapath)
          req(file)
          validate(need(ext == "RData"|ext == "rda", "Please upload a .RData or .rda file"))
          keggGenesfile = get(load(file$datapath))
        }
        else{
          keggGenesfile=NULL
        }
        #DB$keggGenesfile=keggGenesfile
        if(length(DB$MergedDB) > 1) {
          wait(session, "Generating analysis results for selected pathway(s)...")
          
          lapply(1:length(selectPathwayName), function(n){
            ### clear output
            # studyPairs.index = combn(length(DB$MergedStudyNames),2)
            # lapply(1:ncol(studyPairs.index), function(k){
            #   output[[paste0("keggViewDat",n,k)]] = renderText({})
            #   output[[paste0("keggView",n,k)]] = renderText({
            #     validate(
            #       FALSE
            #     )
            #   })
            #   output[[paste0("reactomeViewDat",n,k)]] = renderText({})
            #   output[[paste0("reactomeView",n,k)]] = renderText({
            #     validate(
            #       FALSE
            #     )
            #   })
            # })
            
            if(grepl("/", selectPathwayName[n])){
              selectPathwayName[n] <- sub("/","-",selectPathwayName[n])
            }
            
            output[[paste0("clickedPathwayName",n)]] = renderText({paste0(n,". ",selectPathwayName[n])})
            
            #mdsModel
            if("mdsModel" %in% input$browserVisualizations){
              if(file.exists(paste(DB.load.working.dir(db),"/mdsModel/", selectPathwayName[n], ".png", sep=""))){
                output[[paste0("mdsPathway",n)]] = renderImage({
                  img.src <- paste(DB.load.working.dir(db), 
                                   "/mdsModel/", selectPathwayName[n], ".png", sep="")
                  list(src=img.src, contentType='image/png', alt="module",width=400)
                }, deleteFile = FALSE)
                print("mdsModel done")
              }else{
                temp.dir = paste0(DB.load.working.dir(db),"/mdsModel")
                dir.create(temp.dir,recursive = T)
                setwd(temp.dir)
                res = mdsModel(unlist(c(selectedPathwayAS[n,])),DB$MergedStudyNames,selectPathwayName[n],sep="_")
                setwd(path_old)
                output[[paste0("mdsPathway",n)]] = renderImage({
                  img.src <- paste(DB.load.working.dir(db), 
                                   "/mdsModel/", selectPathwayName[n], ".png", sep="")
                  list(src=img.src, contentType='image/png', alt="module",width=400)
                }, deleteFile = FALSE)
                print("mdsModel done")
              }
            }
            
            #clustModel
            if("clustModel" %in% input$browserVisualizations){
              if(file.exists(paste(DB.load.working.dir(db),"/clustModel/", selectPathwayName[n], ".png", sep=""))){
                output[[paste0("clustPathway",n)]] = renderImage({
                  img.src = paste(DB.load.working.dir(db), "/clustModel/", selectPathwayName[n],'.png',sep="")
                  list(src=img.src, contentType='image/png', alt="module",width=400)
                }, deleteFile = FALSE)
                print("clustModel done")
                
              }else{
                temp.dir = paste0(DB.load.working.dir(db),"/clustModel")
                dir.create(temp.dir,recursive = T)
                setwd(temp.dir)
                cluster.assign.path = SA_algo(unlist(c(selectedPathwaypAS[n,])),DB$MergedStudyNames,sep="_")
                if(length(unique(cluster.assign.path))>1 && class(cluster.assign.path) != "try-error" ){
                  res = clustModel(unlist(c(selectedPathwaypAS[n,])),DB$MergedStudyNames,cluster.assign.path,
                                   selectPathwayName[n],sep="_")
                }
                if(length(unique(cluster.assign.path))==1 || class(cluster.assign.path) == "try-error" ){
                  warning("clustModel only identifies one cluster")
                  res = clustModelOne(unlist(c(selectedPathwaypAS[n,])),DB$MergedStudyNames,
                                      selectPathwayName[n],sep="_")
                }
                setwd(path_old)
                output[[paste0("clustPathway",n)]] = renderImage({
                  img.src <- paste(DB.load.working.dir(db), 
                                   "/clustModel/", selectPathwayName[n], ".png", sep="")
                  list(src=img.src, contentType='image/png', alt="module",width=400)
                }, deleteFile = FALSE)
                print("clustModel done")
              }
            }
            
            #genePM
            if("genePM" %in% input$browserVisualizations){
              if(file.exists(paste(DB.load.working.dir(db),"/genePM/", selectPathwayName[n], ".png", sep=""))){
                output[[paste0("genePM",n)]] = renderImage({
                  img.src <- paste(DB.load.working.dir(db), "/genePM/", selectPathwayName[n],'.png',sep="")
                  list(src=img.src, alt="module",height=400)
                }, deleteFile = FALSE)
                print("genePM done")
                
              }else{
                temp.dir = paste0(DB.load.working.dir(db),"/genePM")
                dir.create(temp.dir,recursive = T)
                setwd(temp.dir)
                
                signPM.list = lapply(DB$MergedDB,function(x) apply(x,1,mean))
                names(signPM.list) = DB$MergedStudyNames
                hm = genePM(signPM.list, pathway.genes=DB$pathway.list[[selectPathwayName[n]]],
                            pathway.name=selectPathwayName[n])
                setwd(path_old)
                output[[paste0("genePM",n)]] = renderImage({
                  img.src <- paste(DB.load.working.dir(db), "/genePM/", selectPathwayName[n],'.png',sep="")
                  list(src=img.src, alt="module",height=400)
                }, deleteFile = FALSE)
                print("genePM done")
              }
            }
            
            ViewSelect2 = grep(paste(input$keggViewSelected2,collapse="|"), 
                                    DB$MergedStudyNames, value=FALSE)
            if(is.null(ViewSelect2)){
              data.pair = combn(length(DB$MergedStudyNames),2)
            }else{
              data.pair = combn(ViewSelect2,2)
            }
            
            #### visulize Kegg topology plots-------------
            if(grepl("KEGG",selectPathwayName[n]) && ("keggView" %in% input$browserVisualizations)){
              kegg.species = input$keggSpecies2
              
              ##match KEGG pathway names to KEGG IDs
              if(!is.null(input$keggIDfile2)){
                file <- input$keggIDfile2
                ext <- tools::file_ext(file$datapath)
                req(file)
                validate(need(ext == "RData"|ext == "rda", "Please upload a .RData or .rda file"))
                KEGG.pathID2name = get(load(file$datapath))
              }
              else{
                message("Retrieve KEGG pathway IDs from KEGGREST...")
                KEGG.pathID2name = lapply(KEGGREST::keggList("pathway",input$keggSpecies2),function(x){
                  strsplit(x," - ")[[1]][-length(strsplit(x," - ")[[1]])]
                })
                names(KEGG.pathID2name) = gsub("path:","",names(KEGG.pathID2name))
              }
              keggk.name.clean = gsub("KEGG ","",selectPathwayName[n])
              pathwayID = gsub(kegg.species,"",names(KEGG.pathID2name)[which(KEGG.pathID2name == keggk.name.clean)])

              ##construct signPMmat for matched pathwayIDs
              if(length(pathwayID) == 0){
                print(paste0("No KEGG pathway ID found for ",selectPathwayName[n],". Please check R package 'KEGGREST' for correct pathway name."))
                }else{
                keggID = paste(kegg.species,pathwayID,sep="")
                kegg.name0 = gsub("/","_",selectPathwayName[n],fixed = T)
                
                dir.path = paste0(DB.load.working.dir(db),"/keggView/", gsub("/","_",selectPathwayName[n],fixed = T))
                
                lapply(1:ncol(data.pair), function(c){
                  ds1 = DB$MergedStudyNames[data.pair[1,c]]
                  ds2 = DB$MergedStudyNames[data.pair[2,c]]
                  
                  img.src = paste(dir.path,"/",kegg.name0,"_",ds1,"_",ds2,".png",sep="")
                  if(file.exists(img.src)){
                    output[[paste0("keggView",n,c)]] = renderImage({
                      list(src=img.src, contentType='image/png', alt="module",width = 800)
                    }, deleteFile = FALSE)
                    output[[paste0("keggViewDat",n,c)]] = renderText({
                      paste0(selectPathwayName[n], " topological plot for data pair ",ds1," and ",ds2)
                    })
                  }else{
                    dir.create(dir.path,recursive = T)
                    setwd(dir.path)
                    
                    ds1 = DB$MergedStudyNames[data.pair[1,c]]
                    ds2 = DB$MergedStudyNames[data.pair[2,c]]
                    
                    dat1 = DB$MergedDB[[data.pair[1,c]]]
                    dat2 = DB$MergedDB[[data.pair[2,c]]]
                    
                    overlap.genes <- intersect(rownames(dat1),DB$pathway.list[[selectPathwayName[n]]])
                    signPM.mat <- cbind(apply(dat1[overlap.genes,],1,mean),
                                        apply(dat2[overlap.genes,],1,mean))
                    colnames(signPM.mat) <- c(ds1,ds2)
                    
                    std.genes <- rownames(signPM.mat)
                    
                    if(input$KEGG.dataGisEntrezID2 == "TRUE"){
                      message("'KEGG.dataGisEntrezID == TRUE', gene names are used for KEGG pathview directly.")
                      
                    }else if(!is.null(input$keggGenesfile2)){
                      file <- input$keggGenesfile2
                      ext <- tools::file_ext(file$datapath)
                      req(file)
                      validate(need(ext == "RData"|ext == "rda", "Please upload a .RData or .rda file"))
                      KEGG.dataG2EntrezID = get(load(file$datapath))
                      entrezID = KEGG.dataG2EntrezID[match(rownames(signPM.mat),KEGG.dataG2EntrezID[,1]),2]
                      rownames(signPM.mat) = entrezID
                      
                      message("Gene names are converted to EntrezID by the provided data frame KEGG.dataG2EntrezID.")
                
                    }else if(kegg.species == "hsa"){
                      map.ls = as.list(as.list(org.Hs.eg.db::org.Hs.egALIAS2EG))
                      entrezID = sapply(rownames(signPM.mat),function(g) map.ls[[g]][[1]])
                      rownames(signPM.mat) = entrezID
                      
                      message("Gene names are converted to EntrezID by org.Hs.eg.db::org.Hs.egALIAS2EG.")
                      
                    }else if(kegg.species == "mmu"){
                      map.ls = as.list(as.list(org.Mm.eg.db::org.Mm.egALIAS2EG))
                      entrezID = sapply(rownames(signPM.mat),function(g) map.ls[[g]][[1]])
                      rownames(signPM.mat) = entrezID
                      
                      message("Gene names are converted to EntrezID by org.Hs.eg.db::org.Mm.egALIAS2EG.")
                      
                    }else if(kegg.species == "rno"){
                      map.ls = as.list(as.list(org.Rn.eg.db::org.Rn.egALIAS2EG))
                      entrezID = sapply(rownames(signPM.mat),function(g) map.ls[[g]][[1]])
                      rownames(signPM.mat) = entrezID
                      
                      message("Gene names are converted to EntrezID by org.Hs.eg.db::org.Rn.egALIAS2EG")
                      
                    }else if(kegg.species == "cel"){
                      map.ls = as.list(as.list(org.Ce.eg.db::org.Ce.egALIAS2EG))
                      entrezID = sapply(rownames(signPM.mat),function(g) map.ls[[g]][[1]])
                      rownames(signPM.mat) = entrezID
                      
                      message("Gene names are converted to EntrezID by org.Hs.eg.db::org.Ce.egALIAS2EG")
                      
                    }else if(kegg.species == "dme"){
                      map.ls = as.list(as.list(org.Dm.eg.db::org.Dm.egALIAS2EG))
                      entrezID = sapply(rownames(signPM.mat),function(g) map.ls[[g]][[1]])
                      rownames(signPM.mat) = entrezID
                      
                      message("Gene names are converted to EntrezID by org.Hs.eg.db::org.Dm.egALIAS2EG")
                      
                    }
                    
                    
                    
                    if(kegg.species == "hsa"|kegg.species == "mmu"){
                      #Average DE genes in each node
                      #Only applicable to human and mouse because both their XML genes and pathview package use EntrezID
                      download.kegg(pathway.id = pathwayID, kegg.species, kegg.dir = ".", file.type="xml")
                      parsePathway = KEGGgraph::parseKGML(paste(getwd(), "/",kegg.species, pathwayID,".xml",sep = ""))
                      parsePathway = KEGGgraph::splitKEGGgroup(parsePathway)
                      
                      entries = KEGGgraph::nodes(parsePathway)
                      types = sapply(entries, KEGGgraph::getType)
                      entryNames = as.list(sapply(entries, KEGGgraph::getName))
                      if(any(types == "group") || any(types=="map")){
                        entryNames = entryNames[!(types %in% c("group","map"))]
                      }
                      entryIds = names(entryNames)
                      entryNames = lapply(1:length(entryNames), function(i) paste(entryNames[[i]],collapse="_"))
                      names(entryNames) = entryIds
                      
                      entryNames.unique = unique(entryNames)
                      
                      xmlG = gsub(paste0(kegg.species,":"),"",unlist(entryNames.unique))
                      xmlG.ls = lapply(xmlG, function(x){
                        strsplit(x,"_")[[1]]
                      })
                      
                      mergePMls = lapply(1:length(xmlG.ls), function(x){
                        genes = xmlG.ls[[x]]
                        cmG = intersect(rownames(signPM.mat),genes)
                        if(length(cmG) !=0){
                          sub.signPM.mat = matrix(signPM.mat[cmG,],ncol = 2)
                          avgPM = apply(sub.signPM.mat, 2, mean)
                          return(avgPM)
                        }
                      })
                      names(mergePMls) = xmlG
                      mergePMmat = do.call(rbind,mergePMls)
                      
                      row.names(mergePMmat) = sapply(row.names(mergePMmat), function(x) strsplit(x,"_")[[1]][1])
                      signPM.mat = mergePMmat
                    }
                    res = pathview(gene.data = signPM.mat, pathway.id = pathwayID,
                                   species = kegg.species, out.suffix = "", kegg.native = T,
                                   key.pos = "bottomright", map.null=T,cex = 0.15)
                    file.rename(paste(keggID,"..multi.png",sep=""), 
                                paste(kegg.name0,"_",ds1,"_",ds2,".png",sep=""))
                    file.remove(paste(keggID,".xml",sep=""))
                    file.remove(paste(keggID,".png",sep=""))
                    
                    output[[paste0("keggView",n,c)]] = renderImage({
                      list(src=img.src, contentType='image/png', alt="module",width = 800)
                    }, deleteFile = FALSE)
                    output[[paste0("keggViewDat",n,c)]] = renderText({
                      paste0(selectPathwayName[n], "  topological plot for data pair ",ds1," and ",ds2)
                    })
                    
                    setwd(path_old)
                  }
                })
              }
            } 
            # else{
            #   lapply(1:ncol(data.pair), function(c){
            #     output[[paste0("keggViewDat",n,c)]] = renderText({})
            #     output[[paste0("keggView",n,c)]] = renderText({
            #       validate(
            #         need(grepl("KEGG",selectPathwayName[n]),'')
            #       )
            #     })
            #   })
            # }
            
            ### visualize Reactome topology plots--------------------
            if(grepl("Reactome",selectPathwayName[n]) && ("reactomeView" %in% input$browserVisualizations)){
              reactome.species = input$reactomeSpecies2
              mcmc.merge.list = DB$MergedDB
              
              ##match Reactome pathway names to Reactome IDs
              if(!is.null(input$reactomeIDfile2)){
                file <- input$reactomeIDfile2
                ext <- tools::file_ext(file$datapath)
                req(file)
                validate(need(ext == "RData"|ext == "rda", "Please upload a .RData or .rda file"))
                Reactome.pathID2name = get(load(file$datapath))
              }
              else{
                message("Retrieve Reactome pathway IDs from reactome.db...")
                pathid2name = as.list(reactomePATHID2NAME)
                pathid2name_species = pathid2name[grep(reactome.species,names(pathid2name))]
                Reactome.pathID2name = lapply(pathid2name_species, function(x) strsplit(x,": ")[[1]][2])
              }
              
              ##match gene names in pathway files to reactome gene names
              if(!is.null(input$reactomeGenesfile2)){
                file <- input$reactomeGenesfile2
                ext <- tools::file_ext(file$datapath)
                req(file)
                validate(need(ext == "RData"|ext == "rda", "Please upload a .RData or .rda file"))
                Reactome.dataG2TopologyGtype = get(load(file$datapath))
                
                #match mcmc.merge.list genes
                mcmcG = row.names(mcmc.merge.list[[1]])
                mcmcG.topology = Reactome.dataG2TopologyGtype[match(mcmcG,Reactome.dataG2TopologyGtype[,1]),2]
                na.index = which(is.na(mcmcG.topology))
                if(length(na.index) !=0){
                  mcmcG.topology.rmna = mcmcG.topology[-na.index]
                  mcmc.merge.list.topology = lapply(mcmc.merge.list, function(x) {
                    x_rm = x[-na.index,]
                    row.names(x_rm) = mcmcG.topology.rmna
                    return(x_rm)
                  })
                }else{
                  mcmc.merge.list.topology = lapply(mcmc.merge.list, function(x) {
                    row.names(x) = mcmcG.topology
                    return(x)
                  })
                }
                #match pathway genes
                select.pathway.list.topology = lapply(DB$pathway.list, function(x){
                  pathwayG.topology = Reactome.dataG2TopologyGtype[match(x,Reactome.dataG2TopologyGtype[,1]),2]
                  pathwayG.topology = pathwayG.topology[-which(is.na(pathwayG.topology)|pathwayG.topology == "")]
                  return(pathwayG.topology)
                })
              }else{
                message("Reactome.dataG2TopologyGtype is not provided, gene names are used to hightlight genes Reactome topology plots directly.")
                
                mcmc.merge.list.topology = mcmc.merge.list
                select.pathway.list.topology = DB$pathway.list
                }
              source_python(system.file("ImageProcess.py", package = "CAMO"))

              aname1 = selectPathwayName[n]
              aname0 = gsub("Reactome ", "", aname1)
              aID = names(Reactome.pathID2name[Reactome.pathID2name==aname0])
              
              if(length(aID) == 0){
                print(paste0("No Reactome pathway ID found for ",selectPathwayName[n]))
              }else{
                reactome.df = data.frame(ID = as.character(aID),
                                         Genes = as.character(sapply(select.pathway.list.topology[aname1], function(x) paste(x,collapse = " "))))
                
                dir.path = paste0(DB.load.working.dir(db),"/reactomeView/", gsub("/","_",aname1,fixed = T))
                
                lapply(1:ncol(data.pair), function(c){
                  ds1 = DB$MergedStudyNames[data.pair[1,c]]
                  ds2 = DB$MergedStudyNames[data.pair[2,c]]
                  
                  img.src <- paste(dir.path,"/",aID,"_",ds1,"_",ds2,".jpeg",sep="")
                  if(file.exists(img.src)){
                    output[[paste0("reactomeView",n,c)]] = renderImage({
                      list(src=img.src, contentType='image/png', alt="module",width = 800)
                    }, deleteFile = FALSE)
                    output[[paste0("reactomeViewDat",n,c)]] = renderText({
                      paste0(selectPathwayName[n], " reactomeView for data pair ",ds1," and ",ds2)
                    })
                  }else{
                    dir.create(dir.path,recursive = T)
                    setwd(dir.path)
                    file.copy(from = system.file("pallete.jpeg", package = "CAMO"), to = getwd())
                    
                    dat1 = mcmc.merge.list.topology[[data.pair[1,c]]]
                    dat2 = mcmc.merge.list.topology[[data.pair[2,c]]]
                    overlap.genes0 = intersect(rownames(dat1),select.pathway.list.topology[[aname1]])
                    if(length(overlap.genes0) != 0){
                      signPM.mat = cbind(apply(dat1[overlap.genes0,],1,mean),
                                         apply(dat2[overlap.genes0,],1,mean))
                      signPM = data.frame(Genes = row.names(signPM.mat),signPM.mat)
                      colnames(signPM) = c("Genes","dat1","dat2")
                      
                      CallFromR(signPM=signPM,ReactomePath=reactome.df,datadir=paste0(getwd(),"/"),pathwayID=unlist(aID, use.names=FALSE))
                      
                      file.rename(paste(aID,"New.jpeg",sep=""),
                                  paste(aID,"_",ds1,"_",ds2,".jpeg",sep=""))
                      file.remove(paste(aID,".jpeg",sep=""))
                      file.remove(paste(aID,"(1).jpeg",sep=""))
                      file.remove(paste(aID,".sbgn",sep=""))
                    }
                    
                    output[[paste0("reactomeView",n,c)]] = renderImage({
                      list(src=img.src, contentType='image/png', alt="module",width = 800)
                    }, deleteFile = FALSE)
                    output[[paste0("reactomeViewDat",n,c)]] = renderText({
                      paste0(selectPathwayName[n], " reactomeView for data pair ",ds1," and ",ds2)
                    })
                    setwd(path_old)
                    
                  }
                  
                })
              }
              setwd(path_old)
            } 
            # else{
            #   lapply(1:ncol(data.pair), function(c){
            #     output[[paste0("reactomeViewDat",n,c)]] = renderText({})
            #     output[[paste0("reactomeView",n,c)]] = renderText({
            #       validate(
            #         need(grepl("Reactome",selectPathwayName[n]),'')
            #       )
            #     })
            #   })
            # }
            
          })
        }
      }
    },session)
    setwd(path_old)
    done(session)
  })
  
  
  observeEvent(input$tuneMD, {
    wait(session, "Running module detection at each module size, may take a while...")
    path_old = getwd()
    try({    
    path_old <- getwd()
    print(paste0("wd:",getwd()))
    
    names(DB$MergedDB) = DB$MergedStudyNames
    kegg.species = input$MD_KEGGspecies
    keggk.name.clean = gsub("KEGG ","",input$MD_KEGGname)
    keggk.name0 = gsub(" / ","_",input$MD_KEGGname,fixed = T)
    
    if(input$MD_pathwayID != ""){
      pathwayID = input$MD_pathwayID
    }else{
      KEGG.pathID2name = lapply(KEGGREST::keggList("pathway",kegg.species),function(x) strsplit(x," - ")[[1]][-length(strsplit(x," - ")[[1]])])
      names(KEGG.pathID2name) = gsub("path:","",names(KEGG.pathID2name))
      pathwayID = gsub(kegg.species,"",names(KEGG.pathID2name)[which(KEGG.pathID2name == keggk.name.clean)])
    }
    
    
    # kegg_id_name <- unlist(as.list(KEGGPATHNAME2ID))
    # keggk.name <- gsub("KEGG ","", input$chosenKeggID)
    # KEGGpathwayID = kegg_id_name[keggk.name]
    KEGGpathwayID = paste0(kegg.species,pathwayID)
    gene_type = input$MD_geneType
    if(input$maxM == ""){
      maxM = NULL
    }else if(as.numeric(input$maxM) < input$minM){
      showNotification(
        ui = "maxM shold not smaller than minM.",
        type = "message",
        duration = 3
      )
    }else{
      maxM = as.numeric(input$maxM)
    }
    if (input$MD_data1 == input$MD_data2 || (is.null(input$MD_data1) & is.null(input$MD_data2))){
      showNotification(
        ui = "Please choose a pair (2) of studies.",
        type = "message",
        duration = 3
      )
    } else {
      MD_data <- c(input$MD_data1,input$MD_data2)
      dir.path <- paste0(DB.load.working.dir(db),"/keggModule/",KEGGpathwayID,"_",keggk.name0)
      
      if(!file.exists(paste(dir.path,"/KEGG_module_", KEGGpathwayID,"_", input$MD_data1, "_",
                           input$MD_data2, "_", gene_type,".RData", sep=""))){
        dir.create(dir.path,recursive = T)
        setwd(dir.path)
        if(input$MD_KEGG.dataGisTopologyG == "TRUE"){
          MD_res = KEGG_module(mcmc.merge.list = DB$MergedDB, 
                               dataset.names = DB$MergedStudyNames,
                               KEGGspecies=input$MD_KEGGspecies,
                               KEGGpathwayID = pathwayID, 
                               KEGG.dataGisTopologyG = TRUE,
                               KEGG.dataG2topologyG = NULL,
                               data.pair = c(input$MD_data1,input$MD_data2),
                               gene_type = gene_type, 
                               DE_PM_cut = input$DE_PM_cut, minM = input$minM, maxM = maxM,
                               B = input$MD_P, cores = input$MD_cores, search_method = input$MD_searchType,
                               reps_eachM = input$reps_eachM, topG_from_previous = input$topG_from_previous, Tm0 = input$Tm0,
                               mu = input$mu, epsilon = input$epsilon, N = input$N, Elbow_plot = T,
                               seed = input$MD_seed,
                               sep = input$node_sep)
          
        }else if(!is.null(input$keggGenes_dat_map)){
          inFile <- input$keggGenes_dat_map
          ext <- tools::file_ext(inFile$datapath)
          req(inFile)
          validate(need(ext == "RData"|ext == "rda", "Please upload a .RData or .rda file"))
          KEGG.dataG2topologyG = get(load(file$datapath))
          
          MD_res = KEGG_module(mcmc.merge.list = DB$MergedDB, 
                               dataset.names = DB$MergedStudyNames,
                               KEGGspecies=input$MD_KEGGspecies,
                               KEGGpathwayID = pathwayID, 
                               KEGG.dataGisTopologyG = FALSE,
                               KEGG.dataG2topologyG = KEGG.dataG2topologyG,
                               data.pair = c(input$MD_data1,input$MD_data2),
                               gene_type = gene_type, 
                               DE_PM_cut = input$DE_PM_cut, minM = input$minM, maxM = maxM,
                               B = input$MD_P, cores = input$MD_cores, search_method = input$MD_searchType,
                               reps_eachM = input$reps_eachM, topG_from_previous = input$topG_from_previous, Tm0 = input$Tm0,
                               mu = input$mu, epsilon = input$epsilon, N = input$N, Elbow_plot = T,
                               seed = input$MD_seed,
                               sep = input$node_sep)
          
        }else if(input$MD_KEGGspecies %in% c("hsa","mmu","rno","cel","dme")){
          MD_res = KEGG_module(mcmc.merge.list = DB$MergedDB, 
                               dataset.names = DB$MergedStudyNames,
                               KEGGspecies=input$MD_KEGGspecies,
                               KEGGpathwayID = pathwayID, 
                               KEGG.dataGisTopologyG = FALSE,
                               KEGG.dataG2topologyG = NULL,
                               data.pair = c(input$MD_data1,input$MD_data2),
                               gene_type = gene_type, 
                               DE_PM_cut = input$DE_PM_cut, minM = input$minM, maxM = maxM,
                               B = input$MD_P, cores = input$MD_cores, search_method = input$MD_searchType,
                               reps_eachM = input$reps_eachM, topG_from_previous = input$topG_from_previous, Tm0 = input$Tm0,
                               mu = input$mu, epsilon = input$epsilon, N = input$N, Elbow_plot = T,
                               seed = input$MD_seed,
                               sep = input$node_sep)
          
        }

        save(MD_res, file = paste("KEGG_module_", KEGGpathwayID,"_", input$MD_data1, "_",
                                  input$MD_data2, "_", gene_type,".RData", sep=""))
        setwd(path_old)
      }else{
        load(paste(dir.path,"/KEGG_module_", KEGGpathwayID,"_", input$MD_data1, "_",
                   input$MD_data2, "_", gene_type,".RData", sep=""))
      }
      output$MDelbowPlot <- renderImage({
        img.src <- paste(dir.path,"/",KEGGpathwayID,"_", gene_type,"_", input$MD_data1, "_", input$MD_data2, "_", input$MD_searchType, "_elbow_plot.png", sep="")
        list(src=img.src, contentType='image/png', alt="module")
      }, deleteFile = FALSE)
      output$MDbestSize <- renderText({
        paste0("Best module size: ", gsub("minG", "", MD_res$bestSize))
      })
      
      message = paste("Module detections are done.")
      sendSuccessMessage(session, message)
    }
    setwd(path_old)
    },session)
    setwd(path_old)
    done(session)
  })
  
  
  observeEvent(input$plotMD, {
    wait(session, "Generating KEGG topological plots with highlighted modules, may take a while...")
    path_old <- getwd()
    try({
      path_old <- getwd()
      names(DB$MergedDB) = DB$MergedStudyNames
      
      kegg.species = input$MD_KEGGspecies
      keggk.name.clean = gsub("KEGG ","",input$MD_KEGGname)
      keggk.name0 = gsub(" / ","_",input$MD_KEGGname,fixed = T)
      
      if(input$MD_pathwayID != ""){
        pathwayID = input$MD_pathwayID
      }else{
        KEGG.pathID2name = lapply(KEGGREST::keggList("pathway",kegg.species),function(x) strsplit(x," - ")[[1]][-length(strsplit(x," - ")[[1]])])
        names(KEGG.pathID2name) = gsub("path:","",names(KEGG.pathID2name))
        pathwayID = gsub(kegg.species,"",names(KEGG.pathID2name)[which(KEGG.pathID2name == keggk.name.clean)])
      }
      
      
      # kegg_id_name <- unlist(as.list(KEGGPATHNAME2ID))
      # keggk.name <- gsub("KEGG ","", input$chosenKeggID)
      # KEGGpathwayID = kegg_id_name[keggk.name]
      KEGGpathwayID = paste0(kegg.species,pathwayID)
      gene_type = input$MD_geneType
      
      if(input$MD_data1 == input$MD_data2){
        stop("Please choose two different studies.")
      }
      MD_data = c(input$MD_data1,input$MD_data2)
      
      MD_size = input$minM
      

      
      dir.path =  paste0(DB.load.working.dir(db),"/keggModule/",KEGGpathwayID,"_",keggk.name0)
      print(dir.path)
      if(file.exists(dir.path)){
        setwd(dir.path)
        MD_res_file = paste(dir.path,"/KEGG_module_", KEGGpathwayID,"_", MD_data[1], "_",
                            MD_data[2], "_", gene_type,".RData",sep = "")
        if(file.exists(MD_res_file)){
          MD_res = get(load(MD_res_file))
        }else{
          stop("Please run module detection algorithm first.")
        }
      }else{
        stop("Please run module detection algorithm first.")
      }

      if(input$MD_whichToDraw != "all"){
        MD_whichToDraw_index = as.numeric(unlist(strsplit(input$MD_whichToDraw, split=",")))-input$minM+1
        MD_whichToDraw_size = as.numeric(unlist(strsplit(input$MD_whichToDraw, split=",")))
      } else{
        MD_whichToDraw_index = c(1:length(MD_res[["minG.ls"]]))
        MD_whichToDraw_size = "all"
      }
      
      output$MDplots = renderUI({
        lapply(1:length(MD_whichToDraw_index), function(i){
          index = MD_whichToDraw_index[i]
          topologyG = MD_res$minG.ls[[index]]$minG
          if(is.null(dim(topologyG))){
            ncol = 1
          }else{
            ncol = ncol(topologyG)
          }
          lapply(1:ncol, function(j){
            list(span(textOutput(ns(paste0("keggMDDat",i,j))), style="color:red;font-size: 20px;"),
                 DT::dataTableOutput(ns(paste0("keggMDDat2",i,j))),
                 plotOutput(ns(paste0("keggMD",i,j)),width="auto",height="auto"))
          })
        })
      })
      
      
      # output$MDplots = renderUI({
      #   lapply(1:count, function(i){
      #     list(textOutput(ns(paste0("keggMDDat",i))),
      #          textOutput(ns(paste0("keggMDDat2",i))),
      #          plotOutput(ns(paste0("keggMD",i)),width="auto",height="auto"))
      #   })
      # })
      
      res = KEGG_module_topology_plot(MD_res, MD_whichToDraw_size)
      lapply(1:length(MD_whichToDraw_index), function(i){
        index = MD_whichToDraw_index[i]
        topologyG = MD_res$minG.ls[[index]]$minG
        if(is.null(dim(topologyG))){
          signPM.mat = MD_res$mergePMmat[topologyG,]
          node.ls = sapply(row.names(signPM.mat), function(x) strsplit(x,input$node_sep))
          id.ls = lapply(node.ls, function(n) paste0(kegg.species,":", n))
          # name.ls = lapply(1:length(id.ls), function(n){
          #   #node_name = paste(id.ls[[n]],collapse = "_")
          #   genes_ID_name = sapply(strsplit(keggList(id.ls[[n]]),";"),function(x) x[1])
          #   node_genes = paste(sapply(1:length(genes_ID_name), function(x){
          #     paste0(gsub(":","",names(genes_ID_name)[x])," (",genes_ID_name[x],")")
          #   }),collapse = "; ")
          #   
          #   paste0("node",n,": ",node_genes)
          # })
          # name.ls = lapply(1:length(id.ls), function(m) {
          #   tmp.ls = strsplit(keggList(id.ls[[m]]),";")
          #   lapply(tmp.ls, function(x) x[1])
          # })

          name.df = data.frame(Node = paste0("node",1:length(id.ls)), 
                               Genes = sapply(1:length(id.ls), function(n){
                                 #node_name = paste(id.ls[[n]],collapse = "_")
                                 genes_ID_name = sapply(strsplit(KEGGREST::keggList(id.ls[[n]]),";"),function(x) x[1])
                                 paste(sapply(1:length(genes_ID_name), function(x){
                                   paste0(gsub(":","",names(genes_ID_name)[x])," (",genes_ID_name[x],")")
                                 }),collapse = "; ")
                               }))
          
          img.src <- paste(getwd(),"/", MD_res$KEGGspecies, MD_res$KEGGpathwayID,"_", MD_res$data.pair[1], "_", 
                           MD_res$data.pair[2], "_", names(MD_res$minG.ls)[index],"_",gene_type, ".png", sep="")
          if(file.exists(img.src)){
            output[[paste0("keggMDDat",i,1)]] = renderText({
              paste0("Module size: ", index+input$minM-1)
            })
            # output[[paste0("keggMDDat2",i,1)]] = renderText({
            #   #paste(unlist(name.ls),sep=" ")
            #   name.ls
            # })
            output[[paste0("keggMDDat2",i,1)]] = DT::renderDataTable({
              #paste(unlist(name.ls),sep=" ")
              name.df
            })
            output[[paste0("keggMD",i,1)]] = renderImage({
              list(src=img.src, contentType='image/png', alt="module",width = 800)
            }, deleteFile = FALSE)
          }
        }else{
          lapply(1:ncol(topologyG), function(j){
            topologyG0 = topologyG[,j]
            signPM.mat = MD_res$mergePMmat[topologyG0,]
            node.ls = sapply(row.names(signPM.mat), function(x) strsplit(x,input$node_sep))
            id.ls = lapply(node.ls, function(n) paste0(MD_res$KEGGspecies,":", n))
            name.df = data.frame(Node = paste0("node",1:length(id.ls)), 
                                 Genes = sapply(1:length(id.ls), function(n){
                                   #node_name = paste(id.ls[[n]],collapse = "_")
                                   genes_ID_name = sapply(strsplit(KEGGREST::keggList(id.ls[[n]]),";"),function(x) x[1])
                                   paste(sapply(1:length(genes_ID_name), function(x){
                                     paste0(gsub(":","",names(genes_ID_name)[x])," (",genes_ID_name[x],")")
                                   }),collapse = "; ")
                                 }))

            img.src <- paste(getwd(),"/",MD_res$KEGGspecies, MD_res$KEGGpathwayID,"_", MD_res$data.pair[1], "_", 
                             MD_res$data.pair[2], "_", names(MD_res$minG.ls)[index], "_", j,"_",gene_type, ".png", sep="")
            if(file.exists(img.src)){
              output[[paste0("keggMDDat",i,j)]] = renderText({
                paste0("Module size: ", index+input$minM-1, " (alternative",j,")")
              })
              output[[paste0("keggMDDat2",i,j)]] = DT::renderDataTable({
                name.df
              })
              output[[paste0("keggMD",i,j)]] = renderImage({
                list(src=img.src, contentType='image/png', alt="module",width = 800)
              }, deleteFile = FALSE)
            }
          })
        }
      })
      setwd(path_old)
      
      message = paste("KEGG topology plots with highlighted local communities are done.")
      sendSuccessMessage(session, message)
    },session)
    setwd(path_old)
    done(session)
  })
  
  
  ##########################
  # Render output/UI       #
  ##########################
  # upload gene names mapping
  # output$MD_choosekeggGenes_dat_map = renderUI({
  #   if(input$MD_KEGG.dataGisTopologyG=='FALSE' & input$MD_KEGGspecies %in% c('hsa','mmu','rno')){
  #     list(selectInput(ns('select_keggGenes_dat_map'), 'Matching options', 
  #                 c("Match gene symbols to Entrez ID by Bioconductor package"="bioc",
  #                   "Upload a file"="upload"),
  #                 selected = NULL),
  #     uiOutput(ns("upload_keggGenes_dat_maps")))
  #   }else if(input$MD_KEGG.dataGisTopologyG=='FALSE' & input$MD_KEGGspecies %in% c('cel')){
  #     list(selectInput(ns('select_keggGenes_dat_map'), 'Matching options',
  #                      c("Match gene symbols to to WormBase sequence name by Bioconductor package 'biomaRt'"="bioc",
  #                        "Upload a file"="upload"),
  #                      selected = NULL),
  #          uiOutput(ns("upload_keggGenes_dat_maps")))
  #   }else if(input$MD_KEGG.dataGisTopologyG=='FALSE' & input$MD_KEGGspecies %in% c('dme')){
  #     list(selectInput(ns('select_keggGenes_dat_map'), 'Matching options',
  #                      c("Match gene symbols to FlyBase CG IDs by Bioconductor package 'org.Dm.eg.db'"="bioc",
  #                        "Upload a file"="upload"),
  #                      selected = NULL),
  #          uiOutput(ns("upload_keggGenes_dat_maps")))
  #   }else if(input$MD_KEGG.dataGisTopologyG=='FALSE' & !input$MD_KEGGspecies %in% c('hsa','mmu','rno','cel','dme')){
  #     uiOutput(ns("upload_keggGenes_dat_maps"))
  #   }
  # })

  output$MD_choosekeggGenes_dat_map = renderUI({
    if(input$MD_KEGG.dataGisTopologyG=='FALSE' & input$MD_KEGGspecies %in% c('hsa','mmu','rno','cel','dme')){
      list(selectInput(ns('select_keggGenes_dat_map'), 'Matching options', 
                       c("Mapping by Bioconductor"="bioc",
                         "Upload mapping file"="upload"),
                       selected = NULL),
           uiOutput(ns("upload_keggGenes_dat_maps")))
    }else if(input$MD_KEGG.dataGisTopologyG=='FALSE' & !input$MD_KEGGspecies %in% c('hsa','mmu','rno','cel','dme')){
      uiOutput(ns("upload_keggGenes_dat_maps"))
    }
  })
  
  output$upload_keggGenes_dat_maps = renderUI({
    if(input$select_keggGenes_dat_map == "upload"){
      fileInput(ns("keggGenes_dat_map"), "Upload a data frame to match gene names in studies to KEGG topological plots (.RData/.rda)", accept = c(".RData",".rda"))
    }
  })
  
  # select pathway to visual pathway ACS/ADS
  output$selectPathway = renderUI({
    selectInput(ns('selectedPathway'), 'Select a pathway:',
                rownames(DB$ACS_ADS_pathway$ACS.mat))
  })
  
  
  output$choosePathway = renderUI({
    selectInput(ns('chosenPathway'), 'Select a pathway:',
                rownames(DB$ACS_ADS_pathway$ACS.mat))
  })
  
  
  output$inputPathway = renderUI({
    numPath = length(input$chosenPathway)
    ViewSelect2 = grep(paste(input$keggViewSelected2,collapse="|"), 
                       DB$MergedStudyNames, value=FALSE)
    if(is.null(ViewSelect2)){
      data.pair = combn(length(DB$MergedStudyNames),2)
    }else{
      data.pair = combn(ViewSelect2,2)
    }
    unlist(lapply(1:numPath, function(n){
      list(#span(textOutput(ns(paste0("clickedPathwayName",n))),style="color:red"),
           plotOutput(ns(paste0("mdsPathway",n)),width="auto",height="auto"),
           plotOutput(ns(paste0("clustPathway",n)),width="auto",height="auto"),
           plotOutput(ns(paste0("genePM",n)),width="auto",height="auto"),
           lapply(1:ncol(data.pair), function(c){
             list(textOutput(ns(paste0("keggViewDat",n,c))),
                  plotOutput(ns(paste0("keggView",n,c)),width="auto",height="auto"),
                  textOutput(ns(paste0("reactomeViewDat",n,c))),
                  plotOutput(ns(paste0("reactomeView",n,c)),width="auto",height="auto"))
           })
      )
    }),recursive = F)
    
    
    #lapply(1:numPath, function(n){
    #DT::dataTableOutput(ns(paste0("clickedPathway",n)))
    #})
  })
  
  # select a kegg pathway to do module detection
  output$MD_chooseKEGGname = renderUI({
    allPathwayNames <- rownames(DB$ACS_ADS_pathway$ACS.mat)
    keggPathwayNames <- allPathwayNames[grep("KEGG", allPathwayNames)]
    selectInput(ns('MD_KEGGname'), 'Select a KEGG pathway:', keggPathwayNames)
  })
  
  # select a pair of datasets for module detection within selected pathway
  # output$MD_dataPair = renderUI({
  #   checkboxGroupInput(ns('MD_data'), "Select a pair of studies for module detection",DB$MergedStudyNames)
  # })
  output$MD_dataPair1 = renderUI({
    selectInput(ns('MD_data1'), "study1",DB$MergedStudyNames)
  })
  output$MD_dataPair2 = renderUI({
    selectInput(ns('MD_data2'), "study2",DB$MergedStudyNames)
  })
  # observeEvent(input$MD_searchType, {
  #   if(input$MD_searchType == "SA") {
  #     output$MD_condPanel <- renderUI({
  #       list(numericInput(ns("reps_eachM"), "Number of searching repetitions at each module sizle", 1, min = 1),
  #            numericInput(ns("topG_from_previous"), "Number of top module results stored as initials for next module sizle", 1, min = 1),
  #            numericInput(ns("Tm0"), "Number of searching repetitions at each module sizle", 10, min = 1),
  #            numericInput(ns("mu"), "Temperature multiplier", 0.95),
  #            numericInput(ns("epsilon"), "Final temperature", 1e-5, ),
  #            numericInput(ns("N"), "Number of maximum annealing times", 3000, min = 1)
  #       )
  #     })
  #   } else {
  #     output$MD_condPanel <- renderUI({NULL})
  #   }
  # })
  
  observe({
    if(is.null(input$browserVisualizations)) {
      output$KEGG_panel = renderUI({NULL})
      output$Reactome_panel = renderUI({NULL})
    }
    if ("keggView" %in% input$browserVisualizations) {
      output$KEGG_panel <- renderUI({
        bsCollapsePanel("Settings for topological visualization",
                        hr(),
                        textInput(ns("keggSpecies2"), "KEGG organism code", "hsa"),
                        radioButtons(ns("KEGG.dataGisEntrezID2"), "Gene names in data matrix are Entrez IDs?", choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "FALSE"),
                        fileInput(ns("keggGenesfile2"), "Upload a data frame mapping (merged) gene names in data [column1] to 
Entrez IDs [column2] (.RData/.rda). If not provided & choose FALSE in the question above & KEGG organism code is one of 'hsa', 'mmu','rno, 'cel' or 'dme', gene symbols will be automatically mapped to Entrez IDs by Bioconductor packages
                                    'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db', 'org.Ce.eg.db' or 'org.Dm.eg.db'", accept = c(".RData","rda")),
                        fileInput(ns("keggIDfile2"), "Upload a list mapping pathway names [content] to KEGG IDs [name] (.RData/.rda). If not provided, ID will be retrieved from KEGGREST.", accept = c(".RData","rda")),
                        checkboxGroupInput(ns('keggViewSelected2'), 'Select a subset of studies to generate topological plots',
                                           DB$MergedStudyNames, selected = DB$MergedStudyNames))
      })
    } else {
      output$KEGG_panel = renderUI({NULL})
    }
    if ("reactomeView" %in% input$browserVisualizations) {
      output$Reactome_panel <- renderUI({
        bsCollapsePanel("Settings for topological visualization",
                        hr(),
                        textInput(ns("reactomeSpecies2"), "Reactome organism code", "HSA"),
                        fileInput(ns("reactomeGenesfile2"), "Upload a data frame mapping (merged) gene names in data [column1] to to gene names in Reactome topology [column2] (.RData/.rda). 
                                  If not provided, it assumes the gene names in data matrix and topology plots are the same type.", accept = c(".RData","rda")),
                        fileInput(ns("reactomeIDfile2"), "Upload a list mapping pathway names [content] to Reactome IDs [name] (.RData/.rda). If not provided, ID will be retrieved from reactome.db", accept = c(".RData","rda")),
                        checkboxGroupInput(ns('keggViewSelected2'), 'Select a subset of studies to generate topological plots',
                                           DB$MergedStudyNames, selected = DB$MergedStudyNames))
        
      })
    } else {
      output$Reactome_panel = renderUI({NULL})
    }
    if(("keggView" %in% input$browserVisualizations) && ("reactomeView" %in% input$browserVisualizations)) {
      output$KEGG_panel = renderUI({
        bsCollapsePanel("Settings for topological visualization",
                        hr(),
                        textInput(ns("keggSpecies2"), "KEGG organism code", "hsa"),
                        radioButtons(ns("KEGG.dataGisEntrezID2"), "Gene names in data matrix are Entrez IDs?", choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "FALSE"),
                        fileInput(ns("keggGenesfile2"), "Upload a data frame mapping (merged) gene names in data [column1] to 
Entrez IDs [column2] (.RData/.rda). If not provided & choose FALSE in the question above & KEGG organism code is one of 'hsa', 'mmu','rno, 'cel' or 'dme', gene symbols will be automatically mapped to Entrez IDs by Bioconductor packages
                                  'org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db', 'org.Ce.eg.db' or 'org.Dm.eg.db'", accept = c(".RData","rda")),
                        fileInput(ns("keggIDfile2"), "Upload a list mapping pathway names [content] to KEGG IDs [name] (.RData/.rda). If not provided, ID will be retrieved from KEGGREST.", accept = c(".RData","rda")),
                        hr(),
                        textInput(ns("reactomeSpecies2"), "Reactome organism code", "HSA"),
                        fileInput(ns("reactomeGenesfile2"), "Upload a data frame mapping (merged) gene names in data [column1] to to gene names in Reactome topology [column2] (.RData/.rda). 
                                  If not provided, it assumes the gene names in data matrix and topology plots are the same type.", accept = c(".RData","rda")),
                        fileInput(ns("reactomeIDfile2"), "Upload a list mapping pathway names [content] to Reactome IDs [name] (.RData/.rda). If not provided, ID will be retrieved from reactome.db", accept = c(".RData","rda")),
                        checkboxGroupInput(ns('keggViewSelected2'), 'Select a subset of studies to generate topological plots',
                                           DB$MergedStudyNames, selected = DB$MergedStudyNames))
      })
      output$Reactome_panel = renderUI({NULL})
    }
  })
}