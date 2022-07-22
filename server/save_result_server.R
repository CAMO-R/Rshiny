save_result_server <- function(input, output, session){
  ns <- NS("save_result")
  
  ##########################
  # Reactive Values        #
  ##########################
  DB <- reactiveValues(
    MergedDB=MergedDB.load(db), 
    MergedSpecies=MergedSpecies.load(db), 
    MergedStudyNames=MergedStudyNames.load(db), 
    globalACS=NULL, globalACSpvalue=NULL,
    pathwayACS=NULL, pathwayACSpvalue=NULL,
    globalADS=NULL, globalADSpvalue=NULL,
    pathwayADS=NULL, pathwayADSpvalue=NULL,
    compType=NULL, pathway.list=NULL,
    pathACS_summary=data.frame(NULL),
    pathADS_summary=data.frame(NULL),
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
      if(length(DB$MergedDB)==2 & length(unique(DB$MergedSpecies))==2){
        DB$compType <- "single"
      }else{
        DB$compType <- "multiple"
      }
    }, session)
    print(paste("saving directory is:, ", DB.load.working.dir(db), sep=""))
  })
  
  # select comparison type
  observeEvent(input$compType, {
    try({
      if(input$compType != DB$compType){
        stop(paste("Please select ", DB$compType, sep=""))
      }
    }, session)
  })
  
  # global and pathway ACS/ADS
  observeEvent(input$ACS_ADS, {
    wait(session, "Calculating pathway ACS/ADS, may take a while")
    try({ 
      ## pathway select
      data(pathway.list) ## include pathway.list
      select.pathway <- pathSelect(DB$MergedDB, pathway.list, 
                                   pathwaysize.lower.cut = input$pathwaysizeLowerCut, 
                                   pathwaysize.upper.cut = input$pathwaysizeUpperCut, 
                                   overlapsize.cut = input$overlapsizeCut, 
                                   med.de.cut = input$medDECut, 
                                   qfisher.cut = input$qfisherCut)
      save(select.pathway, file="select_pathway.RData")
      load("select_pathway.RData")
      DB$pathway.list <- pathway.list[select.pathway]
      
      # single
      if (DB$compType == "single"){
        path_old <- getwd()
        setwd(DB.load.working.dir(db))
        
        ## 1vs1 pathway
        marginOut <- margin_pathway(dat1, dat2, DB$pathway.list, measure=input$measure)
        permOut_pathway <- perm_pathway(dat1,dat2,DB$pathway.list,measure=input$measure,
                                        B=input$permNumPathway,parallel=F,n.cores=n.cores)
        pathwayACS <- ACS_pathway(dat1, dat2, deIndex1, deIndex2, DB$pathway.list, 
                                  measure=input$measure, marginOut)
        print("pathway ACS finished")
        pathwayACSpvalue <- pACS_pathway(dat1, dat2, deIndex1, deIndex2, DB$pathway.list, 
                                         input$measure, pathwayACS, permOut_pathway, marginOut)
        print("global pACS finished")
        pathwayADS <- ADS_pathway(dat1, dat2, deIndex1, deIndex2, DB$pathway.list, 
                                  measure=input$measure, marginOut)
        print("global ADS finished")
        pathwayADSpvalue <- pADS_pathway(dat1, dat2, deIndex1, deIndex2, DB$pathway.list, 
                                         input$measure, pathwayADS, permOut_pathway, marginOut)
        print("global pADS finished")
        ACS_ADS_pathway$ACS <- pathwayACS
        ACS_ADS_pathway$ADS <- pathwayADS
        ACS_ADS_pathway$ACSpvalue <- pathwayACSpvalue
        ACS_ADS_pathway$ADSpvalue <- pathwayADSpvalue
        DB$ACS_ADS_pathway <- ACS_ADS_pathway
        save(ACS_ADS_pathway, file="ACS_ADS_pathway.RData")
        print("pathway ACS/ADS saved")
        
        res_ACS <- cbind(pathwayACS, pathwayACSpvalue)
        rownames(res_ACS) <- names(DB$pathway.list)
        colnames(res_ACS) <- c("ACS","ACS p-val")
        res_ADS <- cbind(pathwayADS, pathwayADSpvalue)
        rownames(res_ADS) <- names(DB$pathway.list)
        colnames(res_ADS) <- c("ADS","ADS p-val")
        
        setwd(path_old)
        
        ## multiple comparison
      }else if(DB$compType == "multiple"){
        
        ## multiple pathway
        path_old <- getwd()
        setwd(DB.load.working.dir(db))
        ACS_ADS_pathway <- multi_ACS_ADS_pathway(mcmc.merge.list=DB$MergedDB, 
                                                 dataset.names=DB$MergedStudyNames, 
                                                 select.pathway.list=DB$pathway.list, 
                                                 measure=input$measure,
                                                 B=input$permNumPathway, 
                                                 parallel=F, n.cores=4)
        save(ACS_ADS_pathway, file="ACS_ADS_pathway.RData")
        
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
      }else{
        stop("At least one study is needed in each species")
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
    }, session)
    done(session)
  })
  
  
  # visulize ACS/ADS for selected pathway
  observeEvent(input$selectedPathway, {
    load(paste0(DB.load.working.dir(db),"/ACS_ADS_pathway.RData"))
    DB$ACS_ADS_pathway <- ACS_ADS_pathway
    load(paste0(DB.load.working.dir(db),"/ACS_ADS_global.RData"))
    DB$ACS_ADS_global <- ACS_ADS_global
    wait(session, "Generating ACS/ADS table for slected pathway...")
    try({
      if(input$selectedPathway!=""){
        path_old <- getwd()
        selectedPathwayACS <- DB$pathACS_summary[input$selectedPathway,]
        selectedPathwayADS <- DB$pathADS_summary[input$selectedPathway,]
        if(length(DB$ACS_ADS_global$ACS)==1){
          res <- matrix(1, length(DB$MergedDB), 
                        length(DB$MergedDB))
          colnames(res) <- rownames(res) <- DB$MergedStudyNames
          print(length(selectedPathwayACS))
          res[1,2] <- paste(round(selectedPathwayACS[1],3), " (p-val= ", 
                            round(selectedPathwayACS[2],3), ")", sep="")
          res[2,1] <- paste(round(selectedPathwayADS[1],3), " (p-val= ", 
                            round(selectedPathwayADS[2],3), ")", sep="")
        } else{
          res <- DB$ACS_ADS_global$ACS
          c = 1
          d = 1
          for(i in 1:(nrow(DB$ACS_ADS_global$ACS))){
            for(j in 1:ncol(DB$ACS_ADS_global$ACS)){
              if(j>i){
                res[i,j] <- selectedPathwayACS[1,c]
                c <- c+1
              }
              else if(j<i){
                res[i,j] <- selectedPathwayADS[1,d]
                d <- d+1
              }
            }
          }
        }
        setwd(path_old)
        output$Show_selected_pathway_name <- renderText(paste0("The pathway you selected is:", row.names(selectedPathwayACS)))
        output$Pathway_ACS_ADS_note <- renderText("The upper triangular matrix contains the ACS scores with their p-values, 
                                                  while the lower triangular matrix contains the ADS scores with their p-values.")
        output$selectedPathwayACS_ADS_Table <- DT::renderDataTable({
          if (!(is.null(res))) {
            res
          }
        })
      }
    }, session)
    done(session)
  })
  
  
  # Tune pathway number of clusters using ACS
  observeEvent(input$tuneK_ACS, {
    print("start_tuneK")
    wait(session, "Tuning ACS K")
    try({ 
      #load(paste0(DB.load.working.dir(db),"/ACS_ADS_pathway.RData"))
      #DB$ACS_ADS_pathway <- ACS_ADS_pathway
      if(length(DB$MergedDB) > 2 & 
         length(unique(DB$MergedSpecies)) == 2) {
        wait(session, "ACS Pathway clustering")
        path_old <- getwd()
        setwd(DB.load.working.dir(db))
        clustDiag(DB$ACS_ADS_pathway, TRUE)
        output$tuneKFig_ACS <- renderImage({
          img.src <- paste(DB.load.working.dir(db), 
                           "/clustDiag/Pathway clustering/consensus012.png", sep="")
          list(src=img.src, contentType='image/png', alt="module")
        }, deleteFile = FALSE)
        setwd(path_old)
      }
    }, session)
    done(session)
  })
  
  # Tune pathway number of clusters using ADS
  observeEvent(input$tuneK_ADS, {
    wait(session, "Tuning ADS K")
    try({ 
      #load("C:/Users/Fancygirl_zyj/Box Sync/CAMO/output/ACS_ADS_pathway.RData")
      #DB$ACS_ADS_pathway <- ACS_ADS_pathway
      if(length(DB$MergedDB) > 2 & 
         length(unique(DB$MergedSpecies)) == 2) {
        wait(session, "ADS Pathway clustering")
        path_old <- getwd()
        setwd(DB.load.working.dir(db))
        clustDiag(DB$ACS_ADS_pathway, FALSE)
        output$tuneKFig_ADS <- renderImage({
          img.src <- paste(DB.load.working.dir(db), 
                           "/clustDiag/Pathway clustering/consensus012.png", sep="")
          list(src=img.src, contentType='image/png', alt="module")
        }, deleteFile = FALSE)
        setwd(path_old)
      }
    }, session)
    done(session)
  })
  
  #Pathway Clustering
  observeEvent(input$pathclust_ACS,{
    wait(session, "ACS Pathway clustering and text mining...")
    path_old <- getwd()
    setwd(DB.load.working.dir(db))
    ## test, remember to comment
    load("ACS_ADS_pathway.RData")
    DB$ACS_ADS_pathway <- ACS_ADS_pathway
    data(pathway.list) ## include pathway.list
    select.pathway <- pathSelect(DB$MergedDB, pathway.list,
                                 pathwaysize.lower.cut = input$pathwaysizeLowerCut,
                                 pathwaysize.upper.cut = input$pathwaysizeUpperCut,
                                 overlapsize.cut = input$overlapsizeCut,
                                 med.de.cut = input$medDECut,
                                 qfisher.cut = input$qfisherCut)
    save(select.pathway, file="select_pathway.RData")
    load("select_pathway.RData")
    DB$pathway.list <- pathway.list[select.pathway]
    ## end test
    
    data(hashtbR0216_phrase0.05) ###allow select later
    res = multiOutput(mcmc.merge.list=DB$MergedDB, 
                      dataset.names=DB$MergedStudyNames, 
                      select.pathway.list=DB$pathway.list,
                      DB$ACS_ADS_pathway,
                      output=c("clustPathway","mdsModel","clustModel"),
                      hashtb=hashtb,pathways=pathways,keggViewSelect = c(1,2),optK=input$K_ACS,
                      kegg_pathname=NULL,hs_gene_id=NULL,thres=0.05,text.permutation = "all",
                      sil_cut=0.1, ADS=FALSE)
    # res = pathwayResRShiny(DB$MergedDB, DB$MergedStudyNames, DB$pathway.list,
    #                        DB$pathwayACS, hashtb=hashtb, pathways=pathways, optK=input$K_ACS)
    ClustDEevid$res = res
    save(res, file="res_ACS.RData")
    print("finish multiple ACS pathway analysis")
    setwd(path_old)
    message = paste("ACS Pathway clustering results saved.")
    sendSuccessMessage(session, message)
    done(session)
  })
  
  #Pathway Clustering
  observeEvent(input$pathclust_ADS,{
    wait(session, "ADS Pathway clustering and text mining...")
    ## test, remember to comment
    # load("C:/Users/Fancygirl_zyj/Box Sync/CAMO/output/ACS_ADS_pathway.RData")
    # DB$ACS_ADS_pathway <- ACS_ADS_pathway
    data(pathway.list) ## include pathway.list
    # select.pathway <- pathSelect(DB$MergedDB, pathway.list, 
    #                              pathwaysize.lower.cut = input$pathwaysizeLowerCut, 
    #                              pathwaysize.upper.cut = input$pathwaysizeUpperCut, 
    #                              overlapsize.cut = input$overlapsizeCut, 
    #                              med.de.cut = input$medDECut, 
    #                              qfisher.cut = input$qfisherCut)
    #save(select.pathway, file="select_pathway.RData")
    ## end test
    path_old <- getwd()
    setwd(DB.load.working.dir(db))
    load("select_pathway.RData")
    DB$pathway.list <- pathway.list[select.pathway]
    data(hashtb_human) ###allow select later
    res <- multiOutput(mcmc.merge.list=DB$MergedDB, 
                       dataset.names=DB$MergedStudyNames, 
                       select.pathway.list=DB$pathway.list,
                       DB$ACS_ADS_pathway,
                       output=c("clustPathway","mdsModel","clustModel", "genePM"),
                       hashtb=hashtb,pathways=pathways,keggViewSelect = c(1,2),optK=input$K_ADS,
                       kegg_pathname=NULL,hs_gene_id=NULL,thres=0.05,text.permutation = "all",
                       sil_cut=0.1, ADS=TRUE)
    # res = pathwayResRShiny(DB$MergedDB, DB$MergedStudyNames, DB$pathway.list,
    #                        DB$pathwayADS, hashtb=hashtb, pathways=pathways, optK=input$K_ADS)
    ClustDEevid$res = res
    save(res, file="res_ADS.RData") # also working
    print("finish multiple ADS pathway analysis")
    setwd(path_old)
    message = paste("ADS Pathway clustering results saved.")
    sendSuccessMessage(session, message)
    done(session)
  })
  
  
  #Plot ACS-DE plot
  observeEvent(input$plotACS_DE, {
    wait(session, "Plotting ACS-ADS-DE plot...")
    try({ 
      data(pathway.list)
      ## test, remember to comment
      # select.pathway <- pathSelect(DB$MergedDB, pathway.list,
      #                              pathwaysize.lower.cut = input$pathwaysizeLowerCut,
      #                              pathwaysize.upper.cut = input$pathwaysizeUpperCut,
      #                              overlapsize.cut = input$overlapsizeCut,
      #                              med.de.cut = input$medDECut,
      #                              qfisher.cut = input$qfisherCut)
      # save(select.pathway, file="select_pathway.RData")
      ## end test
      if(length(DB$MergedDB) > 2 & 
         length(unique(DB$MergedSpecies)) == 2) {
        path_old <- getwd()
        setwd(DB.load.working.dir(db))
        load("ACS_ADS_pathway.RData")
        DB$ACS_ADS_pathway <- ACS_ADS_pathway
        load("select_pathway.RData")
        DB$pathway.list <- pathway.list
        load("res_ACS.RData")
        ClustDEevid$res <- res
        print("load res_ACS successfully")
        
        #calculate DE mat
        ACS_pvalue = DB$ACS_ADS_pathway$ACSpvalue.mat
        ACSlog10p.mat = -log10(ACS_pvalue)
        ADS_pvalue = DB$ACS_ADS_pathway$ADSpvalue.mat
        ADSlog10p.mat = -log10(ADS_pvalue)
        
        all(row.names(ACSlog10p.mat) == row.names(ADSlog10p.mat))
        P <- nrow(ACSlog10p.mat)
        
        # allgenes <- rownames(DB$MergedDB[[1]])
        # pm.list <- lapply(1:length(DB$MergedDB), function(x) 
        #   apply(DB$MergedDB[[x]],1,mean))
        # names(pm.list) <- DB$MergedStudyNames
        # 
        # DEevid <- matrix(NA,P,length(DB$MergedStudyNames))
        # rownames(DEevid) = row.names(ACS_pvalue)
        # colnames(DEevid) = DB$MergedStudyNames
        # for(j in 1:P){
        #   print(j)
        #   pathj <- rownames(ACSlog10p.mat)[j]
        #   print(pathj)
        #   genej <- DB$pathway.list[[pathj]]
        #   intergenej <- intersect(genej,allgenes)
        #   
        #   for(ds in DB$MergedStudyNames){
        #     DEevid[j,ds] <- deevid(abs(pm.list[[ds]])[intergenej])
        #   }
        # }
        # save(DEevid, file = "DEevidMat.RData")
        load("DEevidMat.RData")
        ClustDEevid$DEevid = DEevid
        print("Finish DE matrix")
        
        #Plot
        cluster.lb = ClustDEevid$res$CluterLabelwithoutScatter #withCluster indicator will be used to denote color
        cluster.lb[ClustDEevid$res$scatter.index] = "scatter"
        #print(cluster.lb)
        
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
          
          save(plist, file = "plist.RData")
          return(plist)
        })
        print("Finish Plist")
        #output$ACS_DE11 = renderPlot(ggplot(cars))
        
        Pcordi = cbind(rep(seq(1:length(input$ACS_DEgroup)),each = length(input$ACS_DEgroup)),
                       rep(seq(1:length(input$ACS_DEgroup)),times = length(input$ACS_DEgroup)))
        print("Pcordi is okay")
        
        lapply(1:nrow(Pcordi), function(i){
          print(i)
          r = Pcordi[i,1]
          c = Pcordi[i,2]
          if (r < c){
            print("start r<c")
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
            print("r<c is okay")
          } else if (r > c){
            print("start r>c")
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
            print("r>c is okay")
          } else {
            start("else")
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
            print("else is okay")
          }
        })
        
        print("DE plot is ready")
        setwd(path_old)
      }
    }, session)
    done(session)
  })
  
  
  #Reactive coordinate
  #output$ACS_DEplots_info = renderPrint({
  #for (r in 1:length(input$ACS_DEgroup)) {
  #for (c in 1:length(input$ACS_DEgroup)) {
  #if(!is.null(input[[paste0("plot_hover",r,c)]])){
  #str(input[[paste0("plot_hover",r,c)]])
  #}
  #}
  #}
  #})
  
  
  #Reactive pathway table
  observeEvent(input$plotACS_DE,{
    lapply(1:input$numNearPath,function(n) {
      output[[paste0("clickedPathwayName",n)]] = renderText({})
      output[[paste0("clickedPathway",n)]] = DT::renderDataTable({})
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
    if ((!is.null(ClustDEevid$res)) & (!is.null(ClustDEevid$DEevid))){
      cluster = ClustDEevid$res$CluterLabelwithScatter #withCluster indicator will be used to denote color
      DEevid = ClustDEevid$DEevid
      tmWd = sapply(1:length(ClustDEevid$res$Tm_filtered), function(x) paste(row.names(ClustDEevid$res$Tm_filtered[[x]])[1:10],collapse = ", "))
      tmDf = data.frame(cluster = 1:length(unique(cluster)),KeyWds = tmWd)
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
                               xvar = ds1,
                               yvar = ds2,
                               threshold = Inf, maxpoints = input$numNearPath,
                               addDist = TRUE)
            
            
          } else if (r > c){
            dataT = nearPoints(data.frame(DEevid), input[[paste0("plot_click",r,c)]],
                               xvar = ds1,
                               yvar = ds2,
                               threshold = Inf, maxpoints = input$numNearPath,
                               addDist = TRUE)
          }
          pathway = row.names(dataT)
          print(pathway)
          
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
              # print(row.names(infoMat)[i])
              # print(paste0(input$ACS_DEgroup[i], " (DE=", DE, ")"))
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
  # ClustDEevid$res <- res
  # ClustDEevid$DEevid <- DEevid
  Info = eventReactive(list(click_info(),input$numNearPath),Info_func())
  
  #Generate table info
  observeEvent(Info(),{
    print(Info())
    print(Info()$pathway[[1]])
    print(Info()$infoMat_ls[[1]])
    
    lapply(1:input$numNearPath,function(n) {
      pathway = Info()$pathway[[n]]
      infoMat = Info()$infoMat_ls[[n]]
      
      output[[paste0("clickedPathwayName",n)]] = renderText({paste0(n,". ",pathway)})
      output[[paste0("clickedPathway",n)]] = DT::renderDataTable({infoMat})
      
    })
  })
  
  
  #Replot
  observeEvent(Info(),{
    res = ClustDEevid$res
    DEevid = ClustDEevid$DEevid
    #Plot
    # cluster.lb = ClustDEevid$res$CluterLabelwithoutScatter #withCluster indicator will be used to denote color
    # cluster.lb[ClustDEevid$res$scatter.index] = "scatter"
    
    dpairs = combn(input$ACS_DEgroup,2)
    dpairs.index = combn(length(input$ACS_DEgroup),2)
    
    wait(session,"Updating ACS-ADS-DE plots and table...")
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
      save(plist, file = "plist_clicked.RData")
      return(plist)
    })
    
    print("Finish Plist2")
    
    Pcordi = cbind(rep(seq(1:length(input$ACS_DEgroup)),each = length(input$ACS_DEgroup)),
                   rep(seq(1:length(input$ACS_DEgroup)),times = length(input$ACS_DEgroup)))
    Plist.org = list()
    lapply(1:nrow(Pcordi), function(i){
      print(i)
      r = Pcordi[i,1]
      c = Pcordi[i,2]
      print("start replot")
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
    
    #pathway = NULL#reset clicked pathway
    done(session)
    
  })
  
  
  # Generate all visualization plots
  observeEvent(input$plotAll, {
    wait(session, "Generating individual pathway plots, may take a while...")
    #try({ 
    path_old <- getwd()
    setwd(DB.load.working.dir(db))
    load("C:/Users/Fancygirl_zyj/Box Sync/CAMO/output/ACS_ADS_pathway.RData")
    DB$ACS_ADS_pathway <- ACS_ADS_pathway
    load("select_pathway.RData")
    DB$pathway.list <- pathway.list
    
    data(hashtbR0216_phrase0.05)
    
    ## exclude kegg pathway 
    if("KEGG Metabolic pathways" %in% names(DB$pathway.list)){
      index <- which(names(DB$pathway.list) == "KEGG Metabolic pathways")
      DB$pathway.list <- DB$pathway.list[-index]
      print(paste("number of pathways:",length(DB$pathway.list)))
    }
    if("KEGG D-Arginine and D-ornithine metabolism" %in% names(DB$pathway.list)){
      index <- which(names(DB$pathway.list) == "KEGG D-Arginine and D-ornithine metabolism")
      DB$pathway.list <- DB$pathway.list[-index]
      print(paste("number of pathways:",length(DB$pathway.list)))
    }
    
    library("org.Hs.eg.db")
    # 1 vs 1
    if(DB$compType == "single"){
      #hs_gene_id <- unlist(mget(x=rownames(DB$MergedDB[[1]]),envir=org.Hs.egALIAS2EG))
      
      CS::singleOutput(DB$MergedDB, DB$MergedStudyNames, DB$pathway.list, 
                       output=c("genePM"))
      
    }   # m vs m
    else if(length(DB$MergedDB) > 2 & 
            length(unique(DB$MergedSpecies)) == 2) {
      
      multiOutput(mcmc.merge.list=DB$MergedDB, 
                  dataset.names=DB$MergedStudyNames, 
                  select.pathway.list=DB$pathway.list,
                  DB$ACS_ADS_pathway,
                  output=input$selectVisualizations,
                  hashtb=hashtb,pathways=pathways,keggViewSelect = NULL,optK=input$K_ACS,
                  kegg_pathname=NULL,hs_gene_id=NULL,thres=0.05,text.permutation = "all",
                  sil_cut=0.1, ADS=FALSE)
    }
    
    studyPairs = combn(length(DB$MergedStudyNames),2)
    for (i in 1:ncol(studyPairs)) {
      if(studyPairs[1,i] != studyPairs[2,i]){
        kegg_shiny(DB$MergedDB, DB$MergedStudyNames, 
                   DB$pathway.list,
                   keggViewSelect = c(studyPairs[1,i],studyPairs[2,i]))
      }
    }
    setwd(path_old)
    
    message = paste("Individual pathway plots saved.")
    sendSuccessMessage(session, message)
    #}, session) 
    done(session)
  })
  
  observeEvent(input$KEGG,{
    wait(session, "Generating pathway KEGG plots, may take a while...")
    path_old <- getwd()
    setwd(DB.load.working.dir(db))
    load("ACS_ADS_pathway.RData")
    DB$ACS_ADS_pathway <- ACS_ADS_pathway
    load("select_pathway.RData")
    DB$pathway.list <- pathway.list
    kegg.ind <- grep("KEGG",names(DB$pathway.list))
    data(hashtbR0216_phrase0.05)
    
    
    
    if("KEGG Metabolic pathways" %in% names(DB$pathway.list)){
      index <- which(names(DB$pathway.list) == "KEGG Metabolic pathways")
      DB$pathway.list <- DB$pathway.list[-index]
      print(paste("number of pathways:",length(DB$pathway.list)))
    }
    if("KEGG D-Arginine and D-ornithine metabolism" %in% names(DB$pathway.list)){
      index <- which(names(DB$pathway.list) == "KEGG D-Arginine and D-ornithine metabolism")
      DB$pathway.list <- DB$pathway.list[-index]
      print(paste("number of pathways:",length(DB$pathway.list)))
    }
    library("org.Hs.eg.db")
    
    study1 = which(DB$MergedStudyNames==input$keggStudy1)
    study2 = which(DB$MergedStudyNames==input$keggStudy2)
    
    kegg_shiny(DB$MergedDB, DB$MergedStudyNames, 
               DB$pathway.list,
               keggViewSelect = c(study1, study2))
    # for (i in 1:ncol(studyPairs)) {
    #   if(studyPairs[1,i] != studyPairs[2,i]){
    #     kegg_shiny(DB$MergedDB, DB$MergedStudyNames, 
    #                DB$pathway.list,
    #                keggViewSelect = c(study1, study2))
    #   }
    # }
    setwd(path_old)
  })
  
  
  
  observeEvent(Info(), {
    selectPathwayName = Info()$pathway
    studyPairs.index = combn(length(DB$MergedStudyNames),2)
    
    if(length(DB$MergedDB) > 1) {
      #wait(session, "Model analysis")
      #try({ 
      lapply(1:length(selectPathwayName), function(n){
        if(grepl("/", selectPathwayName[n])){
          selectPathwayName[n] <- sub("/","-",selectPathwayName[n])
        }
        print(paste0("plotP",n))
        if(dir.exists(paste0(DB.load.working.dir(db), 
                             "/mdsModel/"))){
          output[[paste0("mdsPathway",n)]] = renderImage({
            img.src <- paste(DB.load.working.dir(db), 
                             "/mdsModel/", selectPathwayName[n], ".png", sep="")
            list(src=img.src, contentType='image/png', alt="module",width=300)
          }, deleteFile = FALSE)
        }
        if(dir.exists(paste0(DB.load.working.dir(db), "/clustModel/"))){
          output[[paste0("clustPathway",n)]] = renderImage({
            img.src <- paste(DB.load.working.dir(db), "/clustModel/", selectPathwayName[n],'_.png',sep="")
            list(src=img.src, contentType='image/png', alt="module",width=300)
          }, deleteFile = FALSE)
        }
        if(dir.exists(paste0(DB.load.working.dir(db), "/genePM/"))){
          output[[paste0("genePM",n)]] = renderImage({
            img.src <- paste(DB.load.working.dir(db), "/genePM/", selectPathwayName[n],'.jpeg',sep="")
            list(src=img.src, alt="module",height=300)
          }, deleteFile = FALSE)
        }
        
        # if(!grepl("KEGG",selectPathwayName[n])){
        #   lapply(1:ncol(studyPairs.index), function(c){
        #     output[[paste0("keggView",n,c)]] = renderImage(list())
        #     output[[paste0("keggViewDat",n,c)]] = renderText("")
        #   })
        # }
        
        
        if(grepl("KEGG",selectPathwayName[n])){
          kegg_pathname <- unlist(as.list(KEGGPATHID2NAME)) ## KEGG pathway name <-> ID
          path_old <- getwd()
          print(path_old)
          setwd(DB.load.working.dir(db))
          print(getwd())
          load("select_pathway.RData")
          DB$select.pathway.list <- DB$pathway.list[select.pathway]
          
          dir.path <- "keggView"
          dir.output <- getwd()
          if (!file.exists(dir.path)) dir.create(dir.path)
          setwd(paste(dir.output,"/",dir.path,sep=""))
          
          lapply(1:ncol(studyPairs.index), function(c){
            ds1 = DB$MergedStudyNames[studyPairs.index[1,c]]
            ds2 = DB$MergedStudyNames[studyPairs.index[2,c]]
            
            dat1 = DB$MergedDB[[studyPairs.index[1,c]]]
            dat2 = DB$MergedDB[[studyPairs.index[2,c]]]
            
            overlap.genes <- intersect(rownames(dat1),DB$select.pathway.list[[selectPathwayName[n]]])
            signPM.mat <- cbind(apply(dat1[overlap.genes,],1,mean),
                                apply(dat2[overlap.genes,],1,mean))
            colnames(signPM.mat) <- c(ds1,ds2)
            keggk.name1 <- gsub("KEGG ","",selectPathwayName[n]) 
            print(keggk.name1)
            pathwayID <- names(kegg_pathname)[which(kegg_pathname==keggk.name1)]
            #pathwayID <- gsub("hsa","",names(kegg_pathname)[which(kegg_pathname==keggk.name1)])
            print(pathwayID)
            res <- keggView(mat=signPM.mat,pathwayID)
            hsaName <- paste("hsa",pathwayID,sep="")
            print(hsaName)
            print(paste(selectPathwayName[n],"_",ds1,"_",ds2,".png",sep=""))
            file.rename(paste(hsaName,"..multi.png",sep=""), 
                        paste(selectPathwayName[n],"_",ds1,"_",ds2,".png",sep=""))
            file.remove(paste(hsaName,".xml",sep=""))
            file.remove(paste(hsaName,".png",sep=""))   
            
            img.src <- paste(dir.output,"/",dir.path,"/",selectPathwayName[n],"_",ds1,"_",ds2,".png",sep="")
            if(file.exists(img.src)){
              output[[paste0("keggView",n,c)]] = renderImage({
                list(src=img.src, contentType='image/png', alt="module",width = 300)
              }, deleteFile = FALSE)
              output[[paste0("keggViewDat",n,c)]] = renderText({
                paste0("KEGG pathview for data pair ",ds1," and ",ds2)
              })
            }
            
            # if(dir.exists(paste0("/keggView_",ds1,"_",ds2,"/"))){
            #   output[[paste0("keggView",n,c)]] = renderImage({
            #     img.src <- paste(DB.load.working.dir(db), "/keggView_",ds1,"_",ds2,"/", selectPathwayName[n],'.png',sep="")
            #     list(src=img.src, contentType='image/png', alt="module",width = 300)
            #   }, deleteFile = FALSE)
            # }
          })
          setwd(path_old)
        } 
        # else{
        #   lapply(1:ncol(studyPairs.index), function(c){
        #   shinyjs::reset(paste0("keggView",n,c))
        #   output[[paste0("keggViewDat",n,c)]] = renderText({})
        #   })
        # }
        
      })
      
      ### KEGG viewer
      #if(grepl("KEGG",input$selectPathwayName) == TRUE){
      #output$KEGGtopo <- renderImage({
      #img.src <- paste(DB.load.working.dir(db), "/keggView/", 
      #input$selectPathwayName, ".png", sep="")
      #list(src=img.src, contentType='image/png', alt="module")
      #}, deleteFile = FALSE)
      #}
      # }, session)
      done(session)
    }
  })
  
  ##########################
  # Render output/UI       #
  ##########################
  # select pathway to visual pathway ACS/ADS
  output$selectPathway = renderUI({
    selectInput(ns('selectedPathway'), 'Select a pathway to look at ACS/ADS table',
                rownames(DB$pathACS_summary), selected = "")
  })
  
  # select study to visual KEGG
  output$keggStudy1 = renderUI({
    selectInput(ns('keggStudy1'), 'Study 1 for KEGG pathway viewer', 
                DB$MergedStudyNames[DB$MergedSpecies == unique(DB$MergedSpecies)[1]],
                selected=as.character((DB$MergedStudyNames[DB$MergedSpecies == unique(DB$MergedSpecies)[1]])[1]))
  })
  
  # select study to visual KEGG
  output$keggStudy2 = renderUI({
    selectInput(ns('keggStudy2'), 'Study 2 for KEGG pathway viewer', 
                DB$MergedStudyNames[DB$MergedSpecies == unique(DB$MergedSpecies)[2]],
                selected=as.character((DB$MergedStudyNames[DB$MergedSpecies == unique(DB$MergedSpecies)[2]])[1]))
  })
  
  # select type of cisulizations
  output$selectVisualizations = renderUI({
    checkboxGroupInput(ns('selectVisualizations'), 'Select type(s) of visualization to save',
                       c("clustModel","clustPathway","comemberList","genePM","keggView","mdsModel"))
  })
  
  # comparison type
  output$compType = renderUI({
    selectInput(ns('compType'), 'Selelct comparison type:', 
                c("single", "multiple"),
                selected=as.character(DB$compType))
  })
  
  # select datasets for ACS-DE plot
  output$ACS_DEgroup = renderUI({
    checkboxGroupInput(ns('ACS_DEgroup'), "Select datasets for ACS-ADS-DE plot",DB$MergedStudyNames,selected = DB$MergedStudyNames)
  })
  
  # ACS-DE plots layout
  output$ACS_DEplots = renderUI({
    nCol = nRow = length(input$ACS_DEgroup)
    # Pcordi = cbind(rep(seq(1:length(input$ACS_DEgroup)),each = length(input$ACS_DEgroup)),
    #                rep(seq(1:length(input$ACS_DEgroup)),times = length(input$ACS_DEgroup)))
    # print(Pcordi)
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
    
    #listP = list()
    #for (n in 1:(2*numPath)) {
    #if((n %% 2) == 0){
    #listP[[n]] = DT::dataTableOutput(ns(paste0("clickedPathway",n)))
    #}else{
    #listP[[n]] = span(textOutput(ns(paste0("clickedPathwayName",n))),style="color:red")
    #}
    #}
    #listP
    studyPairs.index = combn(length(DB$MergedStudyNames),2)
    unlist(lapply(1:numPath, function(n){
      list(span(textOutput(ns(paste0("clickedPathwayName",n))),style="color:red"),
           DT::dataTableOutput(ns(paste0("clickedPathway",n))),
           #fluidRow(column(4,plotOutput(ns(paste0("mdsPathway",n)))),
           #column(4,plotOutput(ns(paste0("clustPathway",n)))),
           #column(4,plotOutput(ns(paste0("genePM",n))))),
           plotOutput(ns(paste0("mdsPathway",n)),width="auto",height="auto"),
           plotOutput(ns(paste0("clustPathway",n)),width="auto",height="auto"),
           plotOutput(ns(paste0("genePM",n)),width="auto",height="auto"),
           lapply(1:ncol(studyPairs.index), function(c){
             list(textOutput(ns(paste0("keggViewDat",n,c))),
                  plotOutput(ns(paste0("keggView",n,c)),width="auto",height="auto"))
           })
      )
    }),recursive = F)
    
    #lapply(1:numPath, function(n){
    #DT::dataTableOutput(ns(paste0("clickedPathway",n)))
    #})
    
  })
}

