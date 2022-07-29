pathway_server <- function(input, output, session){
  ns <- NS("pathway")
  
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
      #   stop("At least two species are needed")
      # }
      
      output$pathwayACS_Table <- DT::renderDataTable({
        if (!file.exists(paste(DB.load.working.dir(db),
                               "ACS_ADS_pathway.RData", sep="/"))){
          data.frame(NULL)
        }else{
          load(paste(DB.load.working.dir(db),
                     "ACS_ADS_pathway.RData", sep="/"))
          DB$ACS_ADS_pathway <- ACS_ADS_pathway
          res_ACS <- ACS_ADS_pathway$ACS.mat
          for(i in 1:nrow(res_ACS)){
            for(j in 1:ncol(res_ACS)){
              res_ACS[i,j] <-
                paste(round(ACS_ADS_pathway$ACS.mat[i,j],3), " (pval=",
                      round(ACS_ADS_pathway$ACSpvalue.mat[i,j], 3), ")", sep="")
            }
          }
          res_ACS
        }
      })
      output$pathwayADS_Table <- DT::renderDataTable({
        if (!file.exists(paste(DB.load.working.dir(db),
                               "ACS_ADS_pathway.RData", sep="/"))){
          data.frame(NULL)
        }else{
          load(paste(DB.load.working.dir(db),
                     "ACS_ADS_pathway.RData", sep="/"))
          DB$ACS_ADS_pathway <- ACS_ADS_pathway
          res_ADS <- ACS_ADS_pathway$ADS.mat
          for(i in 1:nrow(res_ADS)){
            for(j in 1:ncol(res_ADS)){
              res_ADS[i,j] <-
                paste(round(ACS_ADS_pathway$ADS.mat[i,j],3), " (pval=",
                      round(ACS_ADS_pathway$ADSpvalue.mat[i,j], 3), ")", sep="")
            }
          }
          res_ADS
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
  
  # global and pathway ACS/ADS
  # output pathway ACS/ADS table

  
  
  pathway.list0 = eventReactive(c(input$pathwayfile,input$pathwayfile_ex),{
    if(input$select_pathwayfile == "upload"){
      if(!is.null(input$pathwayfile)){
        file <- input$pathwayfile
        ext <- tools::file_ext(file$datapath)
        req(file)
        validate(need(ext == "RData"|ext == "rda", "Please upload a .RData or .rda file"))
        pathway.list = get(load(file$datapath))
      }
    }else {
      print(paste0("Selected pathway database: ",input$pathwayfile_ex))
      data(input$pathwayfile_ex)
      pathway.list = unlist(lapply(input$pathwayfile_ex, function(x) eval(parse(text = x))),recursive = F)
      
      # full.pathway.list.short = list("KEGG"=kegg.pathway.list,
      #                                "Reactome"=reactome.pathway.list,
      #                                "GO"=go.pathway.list,
      #                                "Biocarta"=biocarta.pathway.list,
      #                                "PID"=pid.pathway.list,
      #                                "Wiki"=wiki.pathway.list,
      #                                "KEGG_cel"=worm.pathway.list[grep("KEGG",names(worm.pathway.list))],
      #                                "Reactome_CEL"=worm.pathway.list[grep("Reactome",names(worm.pathway.list))])
      # pathway.list = unlist(full.pathway.list.short[input$pathwayfile_ex],recursive = F)
      # names(pathway.list) = sapply(strsplit(names(pathway.list),".",fixed = T),function(x) paste(x[-1],collapse = "_"))
    }
    pathway.list
    
  })
  
  observeEvent(input$ACS_ADS, {
    wait(session, "Calculating pathway c-scores & d-scores, may take a while")
    path_old <- getwd()
    try({ 
      path_old <- getwd()
      setwd(DB.load.working.dir(db))
      print(paste0("Total number of pathways: ",length(pathway.list0())))
      
      if(input$topPathNnum == ""){
        topPath.indStudy.num = NULL
      }else{
        topPath.indStudy.num = as.numeric(input$topPathNnum)
      }
      select.pathway <- pathSelect(DB$MergedDB, pathway.list0(),
                                   pathwaysize.lower.cut = input$pathwaysizeLowerCut,
                                   pathwaysize.upper.cut = input$pathwaysizeUpperCut,
                                   overlapsize.cut = input$overlapsizeCut,
                                   med.de.cut = input$medDECut,
                                   min.de.cut = input$minDECut,
                                   qfisher.cut = input$qfisherCut,
                                   topPath.indStudy.num = topPath.indStudy.num)
      
      print(paste0("Number of selected pathways: ",length(select.pathway)))
      
      pathway.list <- pathway.list0()[select.pathway]
      DB$pathway.list <- pathway.list
      save(pathway.list, file="select_pathway.RData")
      
      load("select_pathway.RData")
      n.cores = reactive({
        ifelse(input$parallel == TRUE, input$cores, 1)
      })
      ACS_ADS_pathway <- multi_ACS_ADS_pathway(mcmc.merge.list=DB$MergedDB,
                                               dataset.names=DB$MergedStudyNames,
                                               select.pathway.list=DB$pathway.list,
                                               measure=input$measure,
                                               B=input$permNumPathway,
                                               parallel=input$parallel,
                                               n.cores=n.cores)
      save(ACS_ADS_pathway, file="ACS_ADS_pathway.RData")
      load("ACS_ADS_pathway.RData")
      setwd(path_old)
      
      DB$ACS_ADS_pathway <- ACS_ADS_pathway
      
      res_ACS <- ACS_ADS_pathway$ACS.mat
      for(i in 1:nrow(res_ACS)){
        for(j in 1:ncol(res_ACS)){
          res_ACS[i,j] <-
            paste(round(ACS_ADS_pathway$ACS.mat[i,j],3), " (pval=",
                  round(ACS_ADS_pathway$ACSpvalue.mat[i,j], 3), ")", sep="")
        }
      }
      
      res_ADS <- ACS_ADS_pathway$ADS.mat
      for(i in 1:nrow(res_ADS)){
        for(j in 1:ncol(res_ADS)){
          res_ADS[i,j] <-
            paste(round(ACS_ADS_pathway$ADS.mat[i,j],3), " (pval=",
                  round(ACS_ADS_pathway$ADSpvalue.mat[i,j], 3), ")", sep="")
        }
      }
      
      # output pathway ACS/ADS table
      DB$pathACS_summary <- res_ACS
      output$pathwayACS_Table <- DT::renderDataTable({
        if (!(is.null(DB$pathACS_summary))) {
          DB$pathACS_summary
        }
      })
      DB$pathADS_summary <- res_ADS
      output$pathwayADS_Table <- DT::renderDataTable({
        if (!(is.null(DB$pathADS_summary))) {
          DB$pathADS_summary
        }
      })
      
      sendSuccessMessage(session, "Pathway-level c-scores & d-scores are saved as 'ACS_ADS_pathway'.")
      
    }, session)
    setwd(path_old)
    done(session)
  })
  
  
  # Tune pathway number of clusters using ACS
  observeEvent(input$tuneK_ACS, {
    wait(session,"Consensus clustering to generate scree plot...")
    path_old <- getwd()
    try({ 
      if(length(DB$MergedDB) > 1) {
        #wait(session, "Pathway clustering")
        path_old = getwd()
        if(file.exists(paste0(DB.load.working.dir(db),"/ACS_ADS_pathway.RData"))){
          load(paste0(DB.load.working.dir(db),"/ACS_ADS_pathway.RData"))
          DB$ACS_ADS_pathway = ACS_ADS_pathway
        }
        
        temp.dir = paste0(DB.load.working.dir(db),"/ConsensusClusteringScreePlots")
        if (!dir.exists(temp.dir)) dir.create(temp.dir,recursive = T)
        setwd(temp.dir)
        
        # if(input$clustACS == "c_scores"){
        #   results = clustPathway(DB$ACS_ADS_pathway$ACSpvalue.mat)
        # }else{
        #   results = clustPathway(DB$ACS_ADS_pathway$ADSpvalue.mat)
        # }
        results = clustPathway(DB$ACS_ADS_pathway$ACSpvalue.mat)
        
        setwd(path_old)
        
        output$tuneKFig <- renderImage({
          img.src <- paste(DB.load.working.dir(db), 
                           "/ConsensusClusteringScreePlots/Pathway clustering/consensus012.png", sep="")
          list(src=img.src, contentType='image/png', alt="module")
        }, deleteFile = FALSE)
        sendSuccessMessage(session, "Scree plot for selecting K is done")
      }
    }, session)
    setwd(path_old)
    done(session)
  })
  
  
  #Pathway Clustering
  hashtb0 = eventReactive(c(input$select_hashtbfile,input$hashtbfile),{
    if(input$select_hashtbfile == "upload"){
      if(!is.null(input$hashtbfile)){
        file <- input$hashtbfile
        ext <- tools::file_ext(file$datapath)
        req(file)
        validate(need(ext == "RData"|ext == "rda", "Please upload a .RData or .rda file"))
        load(file$datapath)
        hashtb = get(load(file$datapath))
      }
    }else if(input$select_hashtbfile == "hsa_text"){
      data(hashtb_hsa)
    }else if(input$select_hashtbfile == "cel_text"){
      data(hashtb_cel)
    }else if(input$select_hashtbfile == "mmu_text"){
      data(hashtb_mmu)
    }else if(input$select_hashtbfile == "rno_text"){
      data(hashtb_rno)
    }else if(input$select_hashtbfile == "dme_text"){
      data(hashtb_dme)
    }else{
      hashtb = NULL
    }
    hashtb
  })
  
  observeEvent(input$pathclust_ACS,{
    wait(session, "Clustering pathways, text mining and generating visualizations...")
    path_old <- getwd()
    try({
      path_old <- getwd()
      setwd(DB.load.working.dir(db))
      load("ACS_ADS_pathway.RData")
      DB$ACS_ADS_pathway <- ACS_ADS_pathway
      load("select_pathway.RData")
      DB$pathway.list <- pathway.list
      
      #data(hashtbR0216_phrase0.05) ###allow select later
      DB$hashtb = hashtb0()
      names(DB$MergedDB) <- DB$MergedStudyNames
      # if(file.exists(paste0("clustPathway/Clustering_Summary_K", input$K_ACS, ".csv"))){
      #   print("found csv file")
      #   file.remove(paste0("clustPathway/Clustering_Summary_K", input$K_ACS,".csv"))
      # }

      res = multiOutput(mcmc.merge.list=DB$MergedDB, 
                        dataset.names=DB$MergedStudyNames, 
                        select.pathway.list=DB$pathway.list,
                        DB$ACS_ADS_pathway,
                        output=c("clustPathway"),
                        optK=input$K_ACS, sil_cut=input$silCut,
                        hashtb=DB$hashtb,
                        use_ADS = F,keywords_cut=1,
                        text.permutation = "all", comemberProb_cut=input$comProbCut)
      
      output$mdsPathway <- renderImage({
        img.src <- paste(DB.load.working.dir(db), 
                         "/clustPathway/mdsPathway_K_", input$K_ACS, ".jpeg", sep="")
        list(src=img.src, contentType='image/png', alt="module")
      }, deleteFile = FALSE)
      output$heatmapPathway <- renderImage({
        img.src <- paste(DB.load.working.dir(db), 
                         "/clustPathway/heatmapPathway_K_", input$K_ACS, ".jpeg", sep="")
        list(src=img.src, contentType='image/png', alt="module")
      }, deleteFile = FALSE)
      
      ClustDEevid$res = res
      
      setwd(path_old)
      
      setwd(DB.load.working.dir(db))
      save(res, file="res_clustPath.RData")
      message = paste0("Pathway clustering results are saved. ",res$msg)
      sendSuccessMessage(session, message)
      
      setwd(path_old)
      if(length(DB$MergedDB) > 1){
        CluterLabelwithScatter = res$CluterLabelwithoutScatter[-res$scatter.index]
        rmClust = setdiff(unique(res$CluterLabelwithoutScatte),unique(CluterLabelwithScatter))
        
        lapply(1:input$K_ACS, function(n){
          if(length(rmClust) == 0 |((length(rmClust) != 0) & (!n %in% rmClust))){
            output[[paste0("clustInfo", n)]] = renderText({
              paste0("Cluster ", n,":")
            })
            output[[paste0("keywords",n)]] = DT::renderDataTable({
              res$Tm_filtered[[as.numeric(n)]]
            })
          }else{
            output[[paste0("clustInfo", n)]] = renderText({
              paste0("Cluster ", n," is removed because of scatterness. Please consider a smaller cluster nunmber.")
            }) 
          }
        })
        
        if(length(DB$MergedDB)>2) {
          lapply(1:input$K_ACS, function(n){
            if(length(rmClust) == 0 |((length(rmClust) != 0) & (!n %in% rmClust))){
              img.src <- paste0(DB.load.working.dir(db),"/comemberPlot/ComemMat_cluster_",n,"_threshold_",input$comProbCut, ".jpeg")
              output[[paste0("comemPlot", n)]] = renderImage({
                list(src=img.src, contentType='image/png', alt="module")
              },deleteFile=FALSE)
            }
          })
        }
      }

    },session)
    setwd(path_old)
    done(session)
  })

  #Plot ACS/ADS-DE plot
  observeEvent(input$plotACS_DE, {
    wait(session, "Generating DE evidence plot...")
    path_old = getwd()
    try({ 
      if(length(DB$MergedDB) > 1  & length(input$ACS_DEgroup) >1) {
        path_old <- getwd()
        if(file.exists(paste0(DB.load.working.dir(db),"/ACS_ADS_pathway.RData"))){
          load(paste0(DB.load.working.dir(db),"/ACS_ADS_pathway.RData"))
          DB$ACS_ADS_pathway = ACS_ADS_pathway
        }
        if(file.exists(paste0(DB.load.working.dir(db),"/select_pathway.RData"))){
          load(paste0(DB.load.working.dir(db),"/select_pathway.RData"))
          DB$pathway.list = pathway.list
        }
        if(file.exists(paste0(DB.load.working.dir(db),"/res_clustPath.RData"))){
          load(paste0(DB.load.working.dir(db),"/res_clustPath.RData"))
          ClustDEevid$res <- res
          print("load res_clustPath successfully")
        }else{
          ClustDEevid$res = NULL
          print("Pathway clustering results not found")
        }
        
        #calculate DE mat
        ACS_pvalue = DB$ACS_ADS_pathway$ACSpvalue.mat
        ACSlog10p.mat = -log10(ACS_pvalue)
        ADS_pvalue = DB$ACS_ADS_pathway$ADSpvalue.mat
        ADSlog10p.mat = -log10(ADS_pvalue)
        
        all(row.names(ACSlog10p.mat) == row.names(ADSlog10p.mat))
        P <- nrow(ACSlog10p.mat)
        
        allgenes <- rownames(DB$MergedDB[[1]])
        pm.list <- lapply(1:length(DB$MergedDB), function(x)
          apply(DB$MergedDB[[x]],1,mean))
        names(pm.list) <- DB$MergedStudyNames
        
        DEevid <- matrix(NA,P,length(DB$MergedStudyNames))
        rownames(DEevid) = row.names(ACS_pvalue)
        colnames(DEevid) = DB$MergedStudyNames
        for(j in 1:P){
          pathj <- rownames(ACSlog10p.mat)[j]
          genej <- DB$pathway.list[[pathj]]
          intergenej <- intersect(genej,allgenes)
          
          for(ds in DB$MergedStudyNames){
            DEevid[j,ds] <- mean(abs(pm.list[[ds]])[intergenej],na.rm = T)
          }
        }
        ClustDEevid$DEevid = DEevid
        # 
        # save(DEevid, file = paste0(DB.load.working.dir(db),"/DEevidMat.RData"))
        # 
        # if(file.exists(paste0(DB.load.working.dir(db),"/DEevidMat.RData"))){
        #   load(paste0(DB.load.working.dir(db),"/DEevidMat.RData"))
        #   ClustDEevid$DEevid = DEevid
        #   print("Finish DE matrix")
        # }
        
        #Plot
        if(!is.null(ClustDEevid$res)){
          cluster.lb = ClustDEevid$res$CluterLabelwithoutScatter #withCluster indicator will be used to denote color
          cluster.lb[ClustDEevid$res$scatter.index] = "scatter"
        }else{
          cluster.lb = NULL
        }
        
        dpairs = combn(input$ACS_DEgroup,2)
        dpairs.index = combn(length(input$ACS_DEgroup),2)
        Plist = lapply(1:ncol(dpairs), function(x) {
          ds1 <- dpairs[1,x]
          ds2 <- dpairs[2,x]
          DEevid1 = DEevid[,ds1]
          DEevid2 = DEevid[,ds2]
          ACSp <- ACS_pvalue[,paste(ds1,ds2,sep="_")]
          ADSp <- ADS_pvalue[,paste(ds1,ds2,sep="_")]
          size.scale = 1/length(input$ACS_DEgroup)
          plist <- ACS_ADS_DE(ds1,ds2,DEevid1,DEevid2, ACSp, ADSp, cluster.lb, size.scale = size.scale)
          #save(plist, file = paste0(DB.load.working.dir(db),"/plist.RData"))
          return(plist)
        })
        
        Pcordi = cbind(rep(seq(1:length(input$ACS_DEgroup)),each = length(input$ACS_DEgroup)),
                       rep(seq(1:length(input$ACS_DEgroup)),times = length(input$ACS_DEgroup)))
        
        lapply(1:nrow(Pcordi), function(i){
          r = Pcordi[i,1]
          c = Pcordi[i,2]
          if (r < c){
            plist.index = which(sapply(1:ncol(dpairs.index), function(y){
              all(dpairs.index[,y] == Pcordi[i,])
            }))
            output[[paste0("ACS_DE",r,c)]] = renderPlot({
              #eval(parse(text = paste0("output$ACS_DE",r,c))) = renderPlot({
              Plist[[plist.index]][[1]]
            })
            output[[paste0("ACS_DEtext",r,c)]] = renderText({
              if(!is.null(input[[paste0("plot_hover",r,c)]])){
                paste("DE_",dpairs[1,plist.index],"=",round(input[[paste0("plot_hover",r,c)]]$x,3),
                      "\nDE_",dpairs[2,plist.index],"=",round(input[[paste0("plot_hover",r,c)]]$y,3),sep = "")
              }
            })
          } else if (r > c){
            plist.index = which(sapply(1:ncol(dpairs.index), function(y){
              (dpairs.index[1,y] == c)&(dpairs.index[2,y] == r)
            }))
            output[[paste0("ACS_DE",r,c)]] = renderPlot({
              #eval(parse(text = paste0("output$ACS_DE",r,c))) = renderPlot({
              Plist[[plist.index]][[2]]
            })
            output[[paste0("ACS_DEtext",r,c)]] = renderText({
              if(!is.null(input[[paste0("plot_hover",r,c)]])){
                paste("DE_",dpairs[2,plist.index],"=",round(input[[paste0("plot_hover",r,c)]]$x,3),
                      "\nDE_",dpairs[1,plist.index],"=",round(input[[paste0("plot_hover",r,c)]]$y,3),sep = "")
              }
            })
          } else {
            #Plist.org[[i]] = rectGrob(gp=gpar(fill="white"))
            output[[paste0("ACS_DE",r,c)]] = renderPlot({
              ggplot()+
                theme_void() +
                coord_fixed(ylim=c(0,-1),xlim=c(0,1)) +
                labs(x="",y="") +
                theme(legend.title = element_blank(),
                      axis.text.x = element_text(size = 32),
                      axis.text.y = element_text(size = 32),
                      panel.border = element_blank())
            })
          }
        })
        
        print("c/d scores - DE evidence plot is done")
        setwd(path_old)
        sendSuccessMessage(session, "c/d scores-DE evidence plot is done")
        
      }else{
        stop("At least two studies should be provided and selected.")
      }
    }, session)
    setwd(path_old)
    done(session)
  })
  
  #Reactive pathway table
  observeEvent(input$plotACS_DE,{
    lapply(1:input$numNearPath,function(n) {
      output[[paste0("clickedPathwayName",n)]] = renderText({})
      #output[[paste0("clickedPathway",n)]] = DT::renderDataTable({})
    })
  })#to reset table
  
  click_info_func = function(){
    clickVec = list()
    for (r in 1:length(input$ACS_DEgroup)) {
      for (c in 1:length(input$ACS_DEgroup)) {
        clickVec = c(clickVec,input[[paste0("plot_click",r,c)]])
      }
    }
    return(clickVec)
  }
  click_info = reactive(click_info_func())
  
  #Generate infoMat and pathways
  Info_func = function(){
    if (!is.null(ClustDEevid$DEevid)){
      #cluster = ClustDEevid$res$CluterLabelwithScatter #withCluster indicator will be used to denote color
      DEevid = ClustDEevid$DEevid
      #tmWd = sapply(1:length(ClustDEevid$res$Tm_filtered), function(x) paste(row.names(ClustDEevid$res$Tm_filtered[[x]])[1:10],collapse = ", "))
      #tmDf = data.frame(cluster = 1:length(unique(cluster)),KeyWds = tmWd)
      dpairs = combn(input$ACS_DEgroup,2)
      dpairs.index = combn(length(input$ACS_DEgroup),2)
      
      pathway = NULL
      Pcordi = cbind(rep(seq(1:length(input$ACS_DEgroup)),each = length(input$ACS_DEgroup)),
                     rep(seq(1:length(input$ACS_DEgroup)),times = length(input$ACS_DEgroup)))
      tbRes = lapply(1:nrow(Pcordi), function(i){
        r = Pcordi[i,1]
        c = Pcordi[i,2]
        
        if(!is.null(input[[paste0("plot_click",r,c)]])){
          print("start pathway searching")
          ds1 = input$ACS_DEgroup[r]
          ds2 = input$ACS_DEgroup[c]
          
          if (r < c){
            dataT = nearPoints(data.frame(DEevid), input[[paste0("plot_click",r,c)]],
                               xvar = ds2,
                               yvar = ds1,
                               threshold = Inf, maxpoints = input$numNearPath,
                               addDist = TRUE)
            
            
          } else if (r > c){
            dataT = nearPoints(data.frame(DEevid), input[[paste0("plot_click",r,c)]],
                               xvar = ds2,
                               yvar = ds1,
                               threshold = Inf, maxpoints = input$numNearPath,
                               addDist = TRUE)
          }
          pathway = row.names(dataT)
          
          infoMat_func  = function(p){
            infoMat = matrix(nrow = length(input$ACS_DEgroup),
                             ncol = length(input$ACS_DEgroup))
            row.names(infoMat) = input$ACS_DEgroup
            colnames(infoMat) = input$ACS_DEgroup
            for (i in 1:length(input$ACS_DEgroup)) {
              for (j in 1:length(input$ACS_DEgroup)) {
                ds1 = input$ACS_DEgroup[i]
                ds2 = input$ACS_DEgroup[j]
                if(i<j){
                  pathway.index = which(row.names(DB$ACS_ADS_pathway[["ACS.mat"]])==p)
                  pair.index = intersect(grep(ds1,colnames(DB$ACS_ADS_pathway[["ACS.mat"]])),
                                         grep(ds2,colnames(DB$ACS_ADS_pathway[["ACS.mat"]])))
                  
                  aACS = round(DB$ACS_ADS_pathway[["ACS.mat"]][pathway.index,pair.index],2)
                  aACSp = round(DB$ACS_ADS_pathway[["ACSpvalue.mat"]][pathway.index,pair.index],2)
                  
                  pathway.indexDE = which(row.names(DEevid)==p)
                  pair.indexDE1 = which(colnames(DEevid)==ds1)
                  #pair.indexDE2 = which(colnames(DEevid)==ds2)
                  # DE = round(DEevid[pathway.indexDE, pair.indexDE1],2)
                  
                  infoMat[i,j] = paste0("ACS=",aACS,"\n(pval=",aACSp,")")
                }else if(i>j){
                  pathway.index = which(row.names(DB$ACS_ADS_pathway[["ADS.mat"]])==p)
                  pair.index = intersect(grep(ds1,colnames(DB$ACS_ADS_pathway[["ADS.mat"]])),
                                         grep(ds2,colnames(DB$ACS_ADS_pathway[["ADS.mat"]])))
                  aADS = round(DB$ACS_ADS_pathway[["ADS.mat"]][pathway.index,pair.index],2)
                  aADSp = round(DB$ACS_ADS_pathway[["ADSpvalue.mat"]][pathway.index,pair.index],2)
                  
                  pathway.indexDE = which(row.names(DEevid)==p)
                  pair.indexDE1 = which(colnames(DEevid)==ds1)
                  #pair.indexDE2 = which(colnames(DEevid)==ds2)
                  #DE = round(DEevid[pathway.indexDE, pair.indexDE1],2)
                  #DEneg2 = round(DEevid[pathway.indexDE, pair.indexDE2],2)
                  
                  infoMat[i,j] = paste0("ADS=",aADS,"(pval=",aADSp,")")
                }else{
                  infoMat[i,j] = ""
                }
              }
              DE = round(DEevid[pathway.indexDE, pair.indexDE1],2)
              row.names(infoMat)[i] = paste0(input$ACS_DEgroup[i], " (DE=", DE, ")")
            }
            return(infoMat)
          }
          infoMat_ls = lapply(pathway, infoMat_func)
          return(list(infoMat_ls = infoMat_ls,pathway=pathway))
        }
      })
      return(unlist(tbRes,recursive = F))
    }else{
      return(NULL)
    }
    
  }
  Info = eventReactive(list(click_info(),input$numNearPath),Info_func())
  
  #Generate table info
  observeEvent(Info(),{
    lapply(1:input$numNearPath,function(n) {
      print(n)
      pathway = Info()$pathway[[n]]
      infoMat = Info()$infoMat_ls[[n]]
      
      print(pathway)

      output[[paste0("clickedPathwayName",n)]] = renderText({paste0(n,". ",pathway)})
      #output[[paste0("clickedPathway",n)]] = DT::renderDataTable({infoMat})
    })
  })
  
  #Replot after clicking
  observeEvent(Info(),{
    #res = ClustDEevid$res
    DEevid = ClustDEevid$DEevid
    #Plot
    # cluster.lb = ClustDEevid$res$CluterLabelwithoutScatter #withCluster indicator will be used to denote color
    # cluster.lb[ClustDEevid$res$scatter.index] = "scatter"
    
    dpairs = combn(input$ACS_DEgroup,2)
    dpairs.index = combn(length(input$ACS_DEgroup),2)
    
    wait(session,"Updating ACS-ADS-DE plots...")
    path_old <- getwd()
    #setwd(DB.load.working.dir(db))
    pathway = Info()$pathway
    clickedVec = ifelse(row.names(DB$ACS_ADS_pathway[["ACS.mat"]]) %in% pathway,"red","grey")
    
    Plist2 = lapply(1:ncol(dpairs), function(x) {
      ds1 <- dpairs[1,x]
      ds2 <- dpairs[2,x]
      DEevid1 = DEevid[,ds1]
      DEevid2 = DEevid[,ds2]
      ACSp <- DB$ACS_ADS_pathway[["ACSpvalue.mat"]][,paste(ds1,ds2,sep="_")]
      ADSp <- DB$ACS_ADS_pathway[["ADSpvalue.mat"]][,paste(ds1,ds2,sep="_")]
      size.scale = 1/length(input$ACS_DEgroup)
      plist <- DEevid_ACS_plot_clicked(ds1,ds2,DEevid1,DEevid2,
                                       ACSp, ADSp, clickedVec, size.scale=size.scale)
      save(plist, file = paste0(DB.load.working.dir(db),"/plist_clicked.RData"))
      return(plist)
    })
    print("Finish updated plot list")
    setwd(path_old)
    
    print("start replot")
    Pcordi = cbind(rep(seq(1:length(input$ACS_DEgroup)),each = length(input$ACS_DEgroup)),
                   rep(seq(1:length(input$ACS_DEgroup)),times = length(input$ACS_DEgroup)))
    Plist.org = list()
    lapply(1:nrow(Pcordi), function(i){
      #print(i)
      r = Pcordi[i,1]
      c = Pcordi[i,2]
      if (r < c){
        plist.index = which(sapply(1:ncol(dpairs.index), function(y){
          all(dpairs.index[,y] == Pcordi[i,])
        }))
        Plist.org[[i]] = Plist2[[plist.index]][[1]]
        output[[paste0("ACS_DE",r,c)]] = renderPlot({
          #eval(parse(text = paste0("output$ACS_DE",r,c))) = renderPlot({
          Plist2[[plist.index]][[1]]
        })  
        output[[paste0("ACS_DEtext",r,c)]] = renderText({
          if(!is.null(input[[paste0("plot_hover",r,c)]])){
            paste("DE_",dpairs[1,plist.index],"=",round(input[[paste0("plot_hover",r,c)]]$x,3),
                  "\nDE_",dpairs[2,plist.index],"=",round(input[[paste0("plot_hover",r,c)]]$y,3),sep = "")
          }
        })  
        
      } else if (r > c){
        plist.index = which(sapply(1:ncol(dpairs.index), function(y){
          (dpairs.index[1,y] == c)&(dpairs.index[2,y] == r)
        }))
        Plist.org[[i]] = Plist2[[plist.index]][[2]]
        output[[paste0("ACS_DE",r,c)]] = renderPlot({
          #eval(parse(text = paste0("output$ACS_DE",r,c))) = renderPlot({
          Plist2[[plist.index]][[2]]
        })
        output[[paste0("ACS_DEtext",r,c)]] = renderText({
          if(!is.null(input[[paste0("plot_hover",r,c)]])){
            paste("DE_",dpairs[2,plist.index],"=",round(input[[paste0("plot_hover",r,c)]]$x,3),
                  "\nDE_",dpairs[1,plist.index],"=",round(input[[paste0("plot_hover",r,c)]]$y,3),sep = "")
          }
        })  
      } else {
        Plist.org[[i]] = rectGrob(gp=gpar(fill="white"))
        output[[paste0("ACS_DE",r,c)]] = renderPlot({
          ggplot()+
            theme_void() +
            coord_fixed(ylim=c(0,-1),xlim=c(0,1)) + 
            labs(x="",y="") +
            theme(legend.title = element_blank(),
                  axis.text.x = element_text(size = 8),
                  axis.text.y = element_text(size = 8),
                  panel.border = element_blank())
        })
      }
      shinyjs::reset(paste0("plot_click",r,c))
    })
    
    pathway = NULL#reset clicked pathway
    done(session)
    
  })
  
  ##########################
  # Render output/UI       #
  ##########################
  output$para_cores = renderUI({
    if(input$parallel == TRUE){
      numericInput(ns("cores"), "Number of cores", 1)
    }
  })
  
  
  observeEvent(input$K_ACS,{
    output$comemPlots = renderUI({
      lapply(1:input$K_ACS, function(n){
        list(
          span(textOutput(ns(paste0("clustInfo",n))), style="color:red;font-size: 20px;"),
          DT::dataTableOutput(ns(paste0("keywords",n))),
          plotOutput(ns(paste0("comemPlot",n))),
          hr()
        )

      })
    })

  })
  
  
  # select datasets for ACS-DE plot
  output$ACS_DEgroup = renderUI({
    checkboxGroupInput(ns('ACS_DEgroup'), "Select studies",DB$MergedStudyNames,selected = DB$MergedStudyNames)
  })
  
  # ACS-DE plots layout
  output$ACS_DEplots = renderUI({
    nCol = nRow = length(input$ACS_DEgroup)
    if(nCol == 0 | nRow == 0){
      Colwidth = 2
      Rowwidth = "200px"
    }else{
      Colwidth = floor(12/nCol)
      Rowwidth = 1000/nRow
    }
    fluidRow(
      #tagList(
      lapply(1:nCol, function(c){
        column(Colwidth,lapply(1:nRow, function(r) {#change width later
          list(plotOutput(ns(paste0("ACS_DE",r,c)),height = Rowwidth,
                          click = ns(paste0("plot_click",r,c)),
                          hover = ns(paste0("plot_hover",r,c))),
               verbatimTextOutput(ns(paste0("ACS_DEtext",r,c)),placeholder = F)
          )
        }))
      })
      #)
    )
  })
  
  output$clickedPathway = renderUI({
    numPath = input$numNearPath
    studyPairs.index = combn(length(DB$MergedStudyNames),2)
    unlist(lapply(1:numPath, function(n){
      list(span(textOutput(ns(paste0("clickedPathwayName",n))),style="color:red")
      )
    }),recursive = F)
  })
  
  output$upload_pathwayfile = renderUI({
    if(input$select_pathwayfile == "upload"){
      fileInput(ns("pathwayfile"), "Upload pathway list file (.RData/.rda)", accept = c(".RData",".rda"))
    }
  })
  
  output$exist_pathwayfile = renderUI({
    if(input$select_pathwayfile == "exist"){
      checkboxGroupInput(ns("pathwayfile_ex"), "Choose from:",
                         c("KEGG homo sapiens (gene symbols)"="kegg.pathway.list_hsa",
                           "Reactome homo sapiens (gene symbols)"="reactome.pathway.list_hsa",
                           "Gene Ontology homo sapiens (gene symbols)"="go.pathway.list_hsa",
                           "Biocarta homo sapiens (gene symbols)"="biocarta.pathway.list_hsa",
                           "Pathway Interaction Database homo sapiens (gene symbols)"="pid.pathway.list_hsa",
                           "WikiPathways homo sapiens (gene symbols)"="wiki.pathway.list_hsa",
                           
                           "KEGG caenorhabditis elegans (sequence names)"="kegg.pathway.list_cel",
                           "Reactome caenorhabditis elegans (sequence names)"="kegg.pathway.list_cel",
                           "KEGG caenorhabditis elegans (gene symbols)"="kegg.pathway.list_cel_GeneNames",
                           "Reactome caenorhabditis elegans (gene symbols)"="reactome.pathway.list_cel_GeneNames",
                           
                           "KEGG drosophila melanogaster (gene symbols)"="kegg.pathway.list_dme",
                           "Reactome drosophila melanogaster (gene symbols)"="reactome.pathway.list_dme",
                           
                           "KEGG mus musculus (gene symbols)"="kegg.pathway.list_mmu",
                           "Reactome mus musculus (gene symbols)"="reactome.pathway.list_mmu",
                           
                           "KEGG rattus norvegicus (gene symbols)"="kegg.pathway.list_rno",
                           "Reactome rattus norvegicus (gene symbols)"="reactome.pathway.list_rno"))
    }
  })
  
  output$upload_hashtbfile = renderUI({
    if(input$select_hashtbfile == "upload"){
      fileInput(ns("hashtbfile"), "Upload a noun-pathway file for text mining (.RData/.rda)", accept = c(".RData",".rda"))
    }
  })
  
}

