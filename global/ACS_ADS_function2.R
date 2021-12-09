##################################
###### ACS scores ################
##################################

CSY <- function(d1,d2,index1){
  Sens <- calcSensC(d1[index1,],d2[index1,])
  Spec <- calcSpecC(d1[-index1,],d2[-index1,])
  
  CS <- Sens + Spec - 1
  if(is.nan(CS)) {
    CS <- 0
  }
  return(CS)
}

CSF <- function(d1,d2,index1,index2){
  Sens <- calcSensC(d1[index1,],d2[index1,])
  Prec <- calcPrecC(d1[index2,],d2[index2,])
  
  CS <- (2*Sens*Prec)/(Sens+Prec)
  if(is.nan(CS)) {
    CS <- 0
  }
  return(CS)
}

CSG <- function(d1,d2,index1,index2){
  Sens <- calcSensC(d1[index1,],d2[index1,])
  Spec <- calcSpecC(d1[-index1,],d2[-index1,])
  
  CS <- sqrt(Sens*Spec)
  if(is.nan(CS)) {
    CS <- 0
  }
  return(CS)
}


ECSY <- function(d1,d2){
  ESens <- calcESensC(d1,d2)
  ESpec <- calcESpecC(d1,d2)
  
  ECS <- ESens + ESpec - 1
  if(is.nan(ECS)) {
    ECS <- 0
  }
  return(ECS)
}

ECSF <- function(d1,d2){
  ESens <- calcESensC(d1,d2)
  EPrec <- calcEPrecC(d1,d2)
  
  ECS <- (2*ESens*EPrec)/(ESens+EPrec)
  if(is.nan(ECS)) {
    ECS <- 0
  }
  return(ECS)
}

ECSG <- function(d1,d2){
  ESens <- calcESensC(d1,d2)
  ESpec <- calcESpecC(d1,d2)
  
  ECS <- sqrt(ESens*ESpec)
  if(is.nan(ECS)) {
    ECS <- 0
  }
  return(ECS)
}

permCSY <- function(d1,d2){
  Sens <- calcSensC(d1,d2)
  Spec <- calcSpecC(d1,d2)
  
  permCS <- Sens + Spec - 1
  if(is.nan(permCS)) {
    permCS <- 0
  }
  return(permCS)
}

permCSF <- function(d1,d2){
  Sens <- calcSensC(d1,d2)
  Prec <- calcPrecC(d1,d2)
  
  permCS <- (2*Sens*Prec)/(Sens+Prec)
  if(is.nan(permCS)) {
    permCS <- 0
  }
  return(permCS)
}

permCSG <- function(d1,d2){
  Sens <- calcSensC(d1,d2)
  Spec <- calcSpecC(d1,d2)
  
  permCS <- sqrt(Sens*Spec)
  if(is.nan(permCS)) {
    permCS <- 0
  }
  return(permCS)
}


CS <- function(dat1,dat2,deIndex1,deIndex2,measure="Fmeasure"){
  ## same for both global and pathway
  if(measure=="youden"){
    CS <- CSY(dat1,dat2,deIndex1)
  } else if(measure=="Fmeasure"){
    CS <- CSF(dat1,dat2,deIndex1,deIndex2)
  } else if(measure=="geo.mean"){
    CS <- CSG(dat1,dat2,deIndex1)
  }
  return(CS)
}


ECS <- function(dat1,dat2,measure="Fmeasure"){
  ## same for both global and pathway
  if(measure=="youden"){
    ECS <- ECSY(dat1,dat2)
  } else if(measure=="Fmeasure"){
    ECS <- ECSF(dat1,dat2)
  } else if(measure=="geo.mean"){
    ECS <- ECSG(dat1,dat2)
  }
  return(ECS)
}

permCS <- function(dat1,dat2,measure="Fmeasure"){
  ## same for both global and pathway
  if(measure=="youden"){
    pemrCS <- permCSY(dat1,dat2)
  } else if(measure=="Fmeasure"){
    permCS <- permCSF(dat1,dat2)
  } else if(measure=="geo.mean"){
    permCS <- permCSG(dat1,dat2)
  }
  return(permCS)
}


##################################
###### ADS scores ################
##################################


DSY <- function(d1,d2,index1){
  Sens <- calcSensD(d1[index1,],d2[index1,])
  Spec <- calcSpecD(d1[-index1,],d2[-index1,])
  
  DS <- Sens + Spec - 1
  if(is.nan(DS)) {
    DS <- 0
  }
  return(DS)
}

DSF <- function(d1,d2,index1,index2){
  Sens <- calcSensD(d1[index1,],d2[index1,])
  Prec <- calcPrecD(d1[index2,],d2[index2,])
  
  DS <- (2*Sens*Prec)/(Sens+Prec)
  if(is.nan(DS)) {
    DS <- 0
  }
  return(DS)
}

DSG <- function(d1,d2,index1,index2){
  Sens <- calcSensD(d1[index1,],d2[index1,])
  Spec <- calcSpecD(d1[-index1,],d2[-index1,])
  
  DS <- sqrt(Sens*Spec)
  if(is.nan(DS)) {
    DS <- 0
  }
  return(DS)
}


EDSY <- function(d1,d2){
  ESens <- calcESensD(d1,d2)
  ESpec <- calcESpecD(d1,d2)
  
  EDS <- ESens + ESpec - 1
  if(is.nan(EDS)) {
    EDS <- 0
  }
  return(EDS)
}

EDSF <- function(d1,d2){
  ESens <- calcESensD(d1,d2)
  EPrec <- calcEPrecD(d1,d2)
  
  EDS <- (2*ESens*EPrec)/(ESens+EPrec)
  if(is.nan(EDS)) {
    EDS <- 0
  }
  return(EDS)
}

EDSG <- function(d1,d2){
  ESens <- calcESensD(d1,d2)
  ESpec <- calcESpecD(d1,d2)
  
  EDS <- sqrt(ESens*ESpec)
  if(is.nan(EDS)) {
    EDS <- 0
  }
  return(EDS)
}

permDSY <- function(d1,d2){
  Sens <- calcSensD(d1,d2)
  Spec <- calcSpecD(d1,d2)
  
  permDS <- Sens + Spec - 1
  if(is.nan(permDS)) {
    permDS <- 0
  }
  return(permDS)
}

permDSF <- function(d1,d2){
  Sens <- calcSensD(d1,d2)
  Prec <- calcPrecD(d1,d2)
  
  permDS <- (2*Sens*Prec)/(Sens+Prec)
  if(is.nan(permDS)) {
    permDS <- 0
  }
  return(permDS)
}

permDSG <- function(d1,d2){
  Sens <- calcSensD(d1,d2)
  Spec <- calcSpecD(d1,d2)
  
  permDS <- sqrt(Sens*Spec)
  if(is.nan(permDS)) {
    permDS <- 0
  }
  return(permDS)
}


DS <- function(dat1,dat2,deIndex1,deIndex2,measure="Fmeasure"){
  ## same for both global and pathway
  if(measure=="youden"){
    DS <- DSY(dat1,dat2,deIndex1)
  } else if(measure=="Fmeasure"){
    DS <- DSF(dat1,dat2,deIndex1,deIndex2)
  } else if(measure=="geo.mean"){
    DS <- DSG(dat1,dat2,deIndex1)
  }
  return(DS)
}


EDS <- function(dat1,dat2,measure="Fmeasure"){
  ## same for both global and pathway
  if(measure=="youden"){
    EDS <- EDSY(dat1,dat2)
  } else if(measure=="Fmeasure"){
    EDS <- EDSF(dat1,dat2)
  } else if(measure=="geo.mean"){
    EDS <- EDSG(dat1,dat2)
  }
  return(EDS)
}

permDS <- function(dat1,dat2,measure="Fmeasure"){
  ## same for both global and pathway
  if(measure=="youden"){
    pemrDS <- permDSY(dat1,dat2)
  } else if(measure=="Fmeasure"){
    permDS <- permDSF(dat1,dat2)
  } else if(measure=="geo.mean"){
    permDS <- permDSG(dat1,dat2)
  }
  return(permDS)
}



##################################
#### Global (expected value from marginal, 
#### permutate genes to get p-value) 
##################################


perm_global <- function(dat1,dat2,measure="Fmeasure",B){
  G <- nrow(dat1)
  
  #rawCS <- CS(dat1,dat2,deIndex1,deIndex2,measure)
  #expCS <- ECS(dat1,dat2,measure)
  #rawDS <- DS(dat1,dat2,deIndex1,deIndex2,measure)
  #expDS <- EDS(dat1,dat2,measure)
  
  out <- matrix(NA,B,4)
  colnames(out) <- c("permCS","permECS","permDS","permEDS")
  
  for(b in 1:B){
    dat1perm <- dat1[sample(1:G,G,replace = F),]
    dat2perm <- dat2[sample(1:G,G,replace = F),]
    out[b,"permCS"] <- permCS(dat1perm,dat2perm,measure)
    out[b,"permECS"] <- ECS(dat1perm,dat2perm,measure)
    out[b,"permDS"] <- permDS(dat1perm,dat2perm,measure)
    out[b,"permEDS"] <- EDS(dat1perm,dat2perm,measure)
  }  
  return(out)
}

ACS_global <- function(dat1,dat2,deIndex1,deIndex2,
                       measure="Fmeasure"){
  cs <- CS(dat1,dat2,deIndex1,deIndex2,measure)
  ecs <- ECS(dat1,dat2,measure)
  acs <- (cs - ecs)/(1-ecs)
  return(acs)
}    

pACS_global <- function(dat1,dat2,deIndex1,deIndex2,
                        measure="Fmeasure",acs,permOut){
  
  permcs <- permOut[,"permCS"]
  permecs <- permOut[,"permECS"]
  permacs <- (permcs - permecs)/(1-permecs)
  
  p_acs <- (sum(permacs>=acs) + 1)/(length(permacs)+1)
  return(p_acs)
}

ADS_global <- function(dat1,dat2,deIndex1,deIndex2,
                       measure="Fmeasure"){
  ds <- DS(dat1,dat2,deIndex1,deIndex2,measure)
  eds <- EDS(dat1,dat2,measure)
  ads <- (ds - eds)/(1-eds)
  return(ads)
}    

pADS_global <- function(dat1,dat2,deIndex1,deIndex2,
                        measure="Fmeasure",ads,permOut){
  permds <- permOut[,"permDS"]
  permeds <- permOut[,"permEDS"]
  permads <- (permds - permeds)/(1-permeds)
  
  p_ads <- (sum(permads>=ads) + 1)/(length(permads)+1)
  return(p_ads)
}


##################################
#### Pathway (expected value from global, 
#### permute genes to get p-value) 
##################################


margin_pathway <-  function(dat1,dat2,
                            select.pathway.list,
                            measure="Fmeasure"){
  select.pathways <- names(select.pathway.list)
  data_genes <- rownames(dat1)
  pathway.size <- sapply(select.pathway.list,function(x) {
    length(intersect(data_genes,x))})
  K <- length(select.pathways)
  G <- nrow(dat1)
  
  out <- matrix(NA,nrow=K,ncol=2)
  rownames(out) <- select.pathways
  colnames(out) <- c("ECS","EDS")
  
  R <- 20 ##fairly enough
  
  for(k in 1:K){
    #print(k)
    ecsk <- edsk <- rep(NA,R)
    pathsizek <- pathway.size[k]
    for(j in 1:R){
      index <- sample(1:G,pathsizek,replace=F)
      dat1.select <- dat1[index,]
      dat2.select <- dat2[index,]
      ecsk[j] <- ECS(dat1.select,dat2.select,measure)
      edsk[j] <- EDS(dat1.select,dat2.select,measure)
    }
    out[k,"ECS"] <- mean(ecsk)
    out[k,"EDS"] <- mean(edsk)
  }
  
  return(out)
  
}

perm_pathway <- function(dat1,dat2,
                         select.pathway.list,
                         measure="Fmeasure",B,parallel=F,n.cores=4){
  
  select.pathways <- names(select.pathway.list)
  data_genes <- rownames(dat1)
  pathway.size <- sapply(select.pathway.list,function(x) {
    length(intersect(data_genes,x))})
  K <- length(select.pathways)
  G <- nrow(dat1)
  
  #out <- array(1,dim=c(B,K,4),dimnames=
                 #list(1:B,select.pathways,
                      #c("permCS","permECS","permDS","permEDS")))
  
  #rawCS <- CS(dat1,dat2,deIndex1,deIndex2,measure)
  #expCS <- ECS(dat1,dat2,measure)
  #rawDS <- DS(dat1,dat2,deIndex1,deIndex2,measure)
  #expDS <- EDS(dat1,dat2,measure)
  
  out <- array(1,dim=c(B,K,2),dimnames=
                 list(1:B,select.pathways,c("permCS","permDS")))
  
  for(k in 1:K){
    #print(k)
    pathsizek <- pathway.size[k]
    if(parallel == T){
      permFunc = function(b){
        dat1perm <- dat1[sample(1:G,G,replace = F),]
        dat2perm <- dat2[sample(1:G,G,replace = F),]
        index <- sample(1:G,pathsizek,replace=F)
        dat1perm.select <- dat1perm[index,]
        dat2perm.select <- dat2perm[index,]
        permCS_res <- permCS(dat1perm.select,dat2perm.select,measure)
        permDS_res <- permDS(dat1perm.select,dat2perm.select,measure)
        return(list(permCS_res = permCS_res, permDS_res = permDS_res))
      }
      out.ls = mclapply(1:B, permFunc, mc.cores = n.cores)
      for(b in 1:B){
        out[b,k,"permCS"] <- out.ls[[b]]$permCS_res
        out[b,k,"permDS"] <- out.ls[[b]]$permDS_res
      }
    }else{
      for(b in 1:B){
        dat1perm <- dat1[sample(1:G,G,replace = F),]
        dat2perm <- dat2[sample(1:G,G,replace = F),]
        index <- sample(1:G,pathsizek,replace=F)
        dat1perm.select <- dat1perm[index,]
        dat2perm.select <- dat2perm[index,]
        
        out[b,k,"permCS"] <- permCS(dat1perm.select,dat2perm.select,measure)
        #out[b,k,"permECS"] <- ECS(dat1perm.select,dat2perm.select,measure)
        out[b,k,"permDS"] <- permDS(dat1perm.select,dat2perm.select,measure)
        #out[b,k,"permEDS"] <- EDS(dat1perm.select,dat2perm.select,measure)
        
      }  
    }
  }   
  return(out)
}

ACS_pathway <- function(dat1,dat2,deIndex1,deIndex2,
                        select.pathway.list,
                        measure="Fmeasure",marginOut){
  
  select.pathways <- names(select.pathway.list)
  data_genes <- rownames(dat1)
  pathway.size <- sapply(select.pathway.list,function(x) {
    length(intersect(data_genes,x))})
  K <- length(select.pathways)
  G <- nrow(dat1)
  
  acs <- rep(NA,K)
  names(acs) <- select.pathways
  
  for(k in 1:K){
    path_genek <- select.pathway.list[[k]]
    genek <- intersect(path_genek,data_genes)
    dat1_k <- dat1[genek,]
    dat2_k <- dat2[genek,]
    
    if(length(intersect(names(deIndex1), genek))<=3 ){
      deIndex1_k <- 1:nrow(dat1_k)
    } else {
      deIndex1_k <- which(rownames(dat1_k)%in%intersect(names(deIndex1), genek))
    } 
    
    if(length(intersect(names(deIndex2), genek))<=3 ){
      deIndex2_k <- 1:nrow(dat2_k)
    } else {
      deIndex2_k <- which(rownames(dat2_k)%in%intersect(names(deIndex2), genek))
    } 
    
    cs <- CS(dat1_k,dat2_k,deIndex1_k,deIndex2_k,measure)
    ecs <- marginOut[k,"ECS"]
    acs[k] <- (cs - ecs)/(1-ecs)
  }
  return(acs)
}    

pACS_pathway <- function(dat1,dat2,deIndex1,deIndex2,
                         select.pathway.list,
                         measure="Fmeasure",acs,permOut,marginOut){
  
  select.pathways <- names(select.pathway.list)
  K <- length(select.pathways)
  
  p_acs <- rep(NA,K)
  names(p_acs) <- select.pathways
  
  for(k in 1:K){
    
    permcs <- permOut[,k,"permCS"]
    ecs <- marginOut[k,"ECS"]
    #permecs <- permOut[,k,"permECS"]
    permacs <- (permcs - ecs)/(1-ecs)
    
    p_acs[k] <- (sum(permacs>=acs[k]) + 1)/(length(permacs)+1)
  }
  
  return(p_acs)
  
}  


ADS_pathway <- function(dat1,dat2,deIndex1,deIndex2,
                        select.pathway.list,
                        measure="Fmeasure",marginOut){
  
  select.pathways <- names(select.pathway.list)
  data_genes <- rownames(dat1)
  pathway.size <- sapply(select.pathway.list,function(x) {
    length(intersect(data_genes,x))})
  K <- length(select.pathways)
  G <- nrow(dat1)
  
  ads <- rep(NA,K)
  names(ads) <- select.pathways
  
  for(k in 1:K){
    path_genek <- select.pathway.list[[k]]
    genek <- intersect(path_genek,data_genes)
    dat1_k <- dat1[genek,]
    dat2_k <- dat2[genek,]
    
    if(length(intersect(names(deIndex1), genek))<=3 ){
      deIndex1_k <- 1:nrow(dat1_k)
    } else {
      deIndex1_k <- which(rownames(dat1_k)%in%intersect(names(deIndex1), genek))
    } 
    
    if(length(intersect(names(deIndex2), genek))<=3 ){
      deIndex2_k <- 1:nrow(dat2_k)
    } else {
      deIndex2_k <- which(rownames(dat2_k)%in%intersect(names(deIndex2), genek))
    } 
    
    ds <- DS(dat1_k,dat2_k,deIndex1_k,deIndex2_k,measure)
    eds <- marginOut[k,"EDS"]
    ads[k] <- (ds - eds)/(1-eds)
  }
  return(ads)
}    

pADS_pathway <- function(dat1,dat2,deIndex1,deIndex2,
                         select.pathway.list,
                         measure="Fmeasure",ads,permOut,marginOut){
  
  select.pathways <- names(select.pathway.list)
  K <- length(select.pathways)
  
  p_ads <- rep(NA,K)
  names(p_ads) <- select.pathways
  
  for(k in 1:K){
    
    permds <- permOut[,k,"permDS"]
    #permeds <- permOut[,k,"permEDS"]
    eds <- marginOut[k,"EDS"]
    permads <- (permds - eds)/(1-eds)
    
    p_ads[k] <- (sum(permads>=ads[k]) + 1)/(length(permads)+1)
    
  }
  
  return(p_ads)
  
}    

# ##################################
# #### MultiDataset function 
# ##################################
# multi_ACS_ADS_global <- function(mcmc.merge.list,dataset.names,
#                             measure="Fmeasure",B=100){
#   
#   names(mcmc.merge.list) <- dataset.names
#   M <- length(mcmc.merge.list)
#   P <- choose(M,2)
#   ACS <- ACSpvalue <- matrix(NA,M,M)
#   rownames(ACS) <- colnames(ACS) <- 
#     rownames(ACSpvalue) <- colnames(ACSpvalue) <- dataset.names
#   ADS <- ADSpvalue <- matrix(NA,M,M)
#   rownames(ADS) <- colnames(ADS) <- 
#     rownames(ADSpvalue) <- colnames(ADSpvalue) <- dataset.names
#   
#   diag(ACS) <- diag(ACSpvalue) <- diag(ADS) <- diag(ADSpvalue)  <- 1
#   
#   for(i in 1:(M-1)){
#     for(j in (i+1):M){
#       dat1 <- mcmc.merge.list[[i]]
#       dat2 <- mcmc.merge.list[[j]]
#       deIndex1 <- attr(dat1,"DEindex") 
#       deIndex2 <- attr(dat2,"DEindex") 
#       
#       permOut <- perm_global(dat1,dat2,measure="Fmeasure",B=B)
#       ACS[i,j] <- ACS[j,i] <- acs_global <- ACS_global(dat1,dat2,deIndex1,deIndex2,
#                                                       measure=measure)
#       ACSpvalue[i,j] <- ACSpvalue[j,i] <- pACS_global(dat1,dat2,deIndex1,deIndex2,
#                                                       measure=measure,acs=acs_global,permOut=permOut)
#       ADS[i,j] <- ADS[j,i] <- ads_global <- ADS_global(dat1,dat2,deIndex1,deIndex2,
#                                                        measure=measure)
#       ADSpvalue[i,j] <- ADSpvalue[j,i] <- pADS_global(dat1,dat2,deIndex1,deIndex2,
#                                                       measure=measure,ads=ads_global,permOut=permOut)
#       print(paste("pair: dataset ",i," and dataset ",j,sep=""))
#     }
#   }
#   dir.path <- "arsGlobal"
#   if (!file.exists(dir.path)) dir.create(dir.path)
#   write.csv(ACS,file=paste(paste(dir.path,"ACS_global_",sep="/"),M,".csv",sep=""))
#   write.csv(ACSpvalue,file=paste(paste(dir.path,"ACSpvalue_global_",sep="/"),M,".csv",sep=""))
#   write.csv(ADS,file=paste(paste(dir.path,"ADS_global_",sep="/"),M,".csv",sep=""))
#   write.csv(ADSpvalue,file=paste(paste(dir.path,"ADSpvalue_global_",sep="/"),M,".csv",sep=""))
#   
#   out <- list(ACS=ACS,ACSpvalue=ACSpvalue,ADS=ADS,ADSpvalue=ADSpvalue)
#   return(out)
# }
# 
# multi_ACS_ADS_pathway <- function(mcmc.merge.list,dataset.names,
#                              select.pathway.list,
#                              measure="Fmeasure",
#                              B=100,parallel=F,n.cores=4){
#   
#   names(mcmc.merge.list) <- dataset.names
#   M <- length(mcmc.merge.list)
#   P <- choose(M,2)
#   select.pathways <- names(select.pathway.list)
#   data_genes <- rownames(mcmc.merge.list[[1]])
#   pathway.size <- sapply(select.pathway.list,function(x) {
#     length(intersect(data_genes,x))})
#   K <- length(select.pathways)
#   
#   ACS <- ACSpvalue <- ADS <- ADSpvalue <- array(1,dim=c(K,M,M),dimnames=
#                                                 list(select.pathways,dataset.names,dataset.names))
#   for(i in 1:(M-1)){
#     for(j in (i+1):M){
#       print(paste("pair: dataset ",i," and dataset ",j,sep=""))
#       
#       dat1 <- mcmc.merge.list[[i]]
#       dat2 <- mcmc.merge.list[[j]]
#       deIndex1 <- attr(dat1,"DEindex") 
#       deIndex2 <- attr(dat2,"DEindex")         
#       names(deIndex1) <- rownames(dat1)[deIndex1]
#       names(deIndex2) <- rownames(dat2)[deIndex2]
#       data_genes <- rownames(dat1)
#       
#       print("Calculate ECS & EDS for each pathway...")
#       marginOut <- margin_pathway(dat1,dat2,select.pathway.list,
#                                   measure=measure)
#       print("Permutate genes for each pathway...")
#       permOut_pathway <- perm_pathway(dat1,dat2,select.pathway.list,
#                                       measure=measure,B=B,parallel=parallel,n.cores=n.cores)
#       ACS[,i,j] <- ACS[,j,i] <- acs_pathway <- ACS_pathway(dat1,dat2,deIndex1,deIndex2,
#                                                select.pathway.list,
#                                                measure=measure,marginOut=marginOut)
#       
#       ACSpvalue[,i,j] <- ACSpvalue[,j,i] <- pACS_pathway(dat1,dat2,deIndex1,deIndex2,
#                                                          select.pathway.list,
#                                                          measure=measure,acs=acs_pathway,
#                                                          permOut=permOut_pathway,
#                                                          marginOut=marginOut)
#       
#       ADS[,i,j] <- ADS[,j,i] <- ads_pathway <- ADS_pathway(dat1,dat2,deIndex1,deIndex2,
#                                                select.pathway.list,
#                                                measure=measure,marginOut=marginOut)
#       
#       ADSpvalue[,i,j] <- ADSpvalue[,j,i] <- pADS_pathway(dat1,dat2,deIndex1,deIndex2,
#                                                          select.pathway.list,
#                                                          measure=measure,ads=ads_pathway,
#                                                          permOut=permOut_pathway,
#                                                          marginOut=marginOut)
#     }
#   }
#   
#   ACS.mat <- ACSpvalue.mat <- ADS.mat <- ADSpvalue.mat <- data.frame(matrix(0,K,P))
#   rownames(ACS.mat) <- rownames(ACSpvalue.mat) <- rownames(ADS.mat) <- rownames(ADSpvalue.mat) <- select.pathways
#   
#   combinations <- combn(dataset.names,m=2)
#   colnames(ACS.mat) <- colnames(ACSpvalue.mat) <- colnames(ADS.mat) <- colnames(ADSpvalue.mat) <- apply(combinations,2,FUN=function(x) paste(x,collapse ="_"))
#   
#   
#   for(p in 1:P){
#     pairname <- colnames(ACS.mat)[p]
#     name1 <- strsplit(pairname,split="_",fixed=T)[[1]][1]
#     name2 <- strsplit(pairname,split="_",fixed=T)[[1]][2]
#     ACS.mat[,pairname] <- ACS[,name1,name2]
#     ACSpvalue.mat[,pairname] <- ACSpvalue[,name1,name2]
#     ADS.mat[,pairname] <- ADS[,name1,name2]
#     ADSpvalue.mat[,pairname] <- ADSpvalue[,name1,name2]
#     
#   }
#   
#   dir.path <- "arsPathway"
#   if (!file.exists(dir.path)) dir.create(dir.path)
#   write.csv(ACS.mat,file=paste(paste(dir.path,"ACS_pathway_",sep="/"),M,".csv",sep=""))
#   write.csv(ACSpvalue.mat,file=paste(paste(dir.path,"ACSpvalue_pathway_",sep="/"),M,".csv",sep=""))
#   write.csv(ADS.mat,file=paste(paste(dir.path,"ADS_pathway_",sep="/"),M,".csv",sep=""))
#   write.csv(ADSpvalue.mat,file=paste(paste(dir.path,"ADSpvalue_pathway_",sep="/"),M,".csv",sep=""))
#   
#   out <- list(ACS.mat=ACS.mat,ACSpvalue.mat=ACSpvalue.mat,ADS.mat=ADS.mat,ADSpvalue.mat=ADSpvalue.mat)
#   return(out)
# }
