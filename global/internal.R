##########################
###### Check functions ###
##########################

check.rawData <- function(data){
  G <- nrow(data)
  N <- ncol(data)
  if(!is.numeric(data)) {
    stop("expression data not numeric")
  }
  if(is.null(row.names(data))) {
    stop("gene symbol missing")
  }
  if(N <= 3) {
    stop("too few samples")
  }
  if(sum(duplicated(row.names(data)))>0){
    stop("duplicate gene symbols")
  }
}

check.groupData <- function(group){
  l <- nlevels(group)
  if( l != 2 ) {
    stop("not a two-class comparison")
  }
}

check.pData <- function(pData){
  G <- nrow(pData)
  if(!is.numeric(pData)) {
    stop("pvalue data not numeric")
  }
  if(is.null(row.names(pData))) {
    stop("gene symbol missing")
  }
  if(ncol(pData) <2) {
    stop("missing either p-value or effect size")
  }
}


check.compatibility <- function(data, group, case.label, ctrl.label){
  G <- nrow(data)
  N1 <- ncol(data)
  N2 <- length(group)
  if(N1 != N2) {
    stop("expression data and class label have unmatched sample size")
  }
  if(!all(group %in% c(case.label,ctrl.label))){
    stop("including class labels other than the case and control")
  }
  if(sum(group==case.label) <= 1 ||sum(group==ctrl.label) <= 1){
    stop("not enough samples in either case or control group")
  }
}

##########################
###### BayesP part ###
##########################

PtoZ <- function(p2, lfc) {
  sgn <- sign(lfc)
  z <- ifelse(sgn>0, qnorm(p2/2), qnorm(p2/2,lower.tail = F))
  return(z)
}

SelectGamma <- function(p){
  ## Gamma is DE proportion  = 1-pi0
  m <- length(p)
  lambda <- seq(0,0.95,by=0.01)
  pi0 <- sapply(lambda, function(x) sum(p>x)/(m*(1-x))  )
  
  # fit a natural cubic spline
  library(splines)
  dat <- data.frame(pi0=pi0, lambda=lambda)
  lfit <- lm(pi0 ~ ns(lambda, df = 3), data=dat)
  pi0hat <- predict(lfit, data.frame(lambda=1))
  gamma <- 1 - pi0hat
  return(gamma)
}


##########################
##Pathway enrich analysis###
##########################

gsa.fisher <- function(x, background, pathway) {
  ####x is the list of query genes
  ####backgroud is a list of background genes that query genes from
  ####pathway is a list of different pathway genes
  count_table<-matrix(0,2,2)
  x<-toupper(x)
  background<-toupper(background)
  index<-which(toupper(background) %in% toupper(x)==FALSE)
  background_non_gene_list<-background[index]
  x<-toupper(x)
  pathway<-lapply(pathway,function(x) intersect(toupper(background),toupper(x)))
  get.fisher <- function(path) {
    res <- NA
    ####in the gene list and in the pathway
    count_table[1,1]<-sum(x %in% path)
    #count_table[1,1]<-sum(is.na(charmatch(x,path))==0)
    ####in the gene list but not in the pathway
    count_table[1,2]<-length(x)-count_table[1,1]
    ####not in the gene list but in the pathway
    count_table[2,1]<-sum(background_non_gene_list%in% path)
    ####not in the gene list and not in the pathway
    count_table[2,2]<-length(background_non_gene_list)-count_table[2,1]
    matched_gene<-x[x %in% path]
    match_num<-length(matched_gene)
    overlap_info<-array(0,dim=4)
    names(overlap_info)<-c("DE in Geneset","DE not in Genese","NonDE in Geneset","NonDE out of Geneset")
    overlap_info[1]=count_table[1,1]
    overlap_info[2]=count_table[1,2]
    overlap_info[3]=count_table[2,1]
    overlap_info[4]=count_table[2,2]
    if(length(count_table)==4){
      res <- fisher.test(count_table, alternative="greater")$p}
    return(list(p_value=res,match_gene=matched_gene,match_num=match_num,
                fisher_table=overlap_info))
  }
  p_val<-rep(0,length(pathway))
  
  match_gene_list<-list(length(pathway))
  match_gene <- array(0,dim=length(pathway))
  num1<-array(0,dim=length(pathway))
  num2<-matrix(0,nrow=length(pathway),ncol=4)
  colnames(num2)<-c("DE in Geneset","DE not in Genese","NonDE in Geneset","NonDE out of Geneset")
  for(i in 1:length(pathway)){
    result<-get.fisher(pathway[[i]])
    p_val[i]<-result$p_value
    match_gene_list[[i]]<-result$match_gene
    match_gene[i]<-paste(match_gene_list[[i]],collapse="/")
    num1[i]<-result$match_num
    num2[i,]<-result$fisher_table
  }
  names(p_val) <- names(pathway)
  q_val <- p.adjust(p_val, "BH")
  
  summary<-data.frame(pvalue=p_val,
                      qvalue=q_val,
                      DE_in_Set=num2[,1],
                      DE_not_in_Set=num2[,2],
                      NonDE_in_Set=num2[,3],
                      NonDE_not_in_Set=num2[,4])
  
  a<-format(summary,digits=3)
  
  return(a)
}

fisher <- function(x){
  n <- length(x)
  y <- -2*log(x)
  Tf <- sum(y)
  return(1-pchisq(Tf,2*n))
}


##########################
###### SA functions ###
##########################

## energy function

E_tot <- function(delta.mat,a,delta.est){
  ## vector "a" of length n
  ## vector "delta_est" of length K+1: start from theta_0, then ordered from k=1 to K
  n <- nrow(delta.mat)
  K <- length(unique(a))
  #theta_0 <- delta.est[1]
  theta_0 <- 0
  E <- sum(sapply(1:n, function(x) {
    sum(sapply(1:n, function(y){
      if(a[x]==a[y]){
        (delta.mat[x,y] - delta.est[as.character(a[x])])^2
      } else{
        (delta.mat[x,y] - theta_0)^2
      }
    },simplify=T))
  }, simplify=T))
  
  return(E)
}

## Estimate of delta.est (the means)

Est_mean <- function(delta.mat,a){
  # the first element is always the off-diagonal parts
  n <- nrow(delta.mat)
  K <- length(unique(a))
  total <- rep(0,K+1)
  size <- rep(0,K+1)
  deltamean <- rep(0,K+1)
  names(deltamean) <- c(0,sort(unique(a)))
  for(k in 1:K){
    a_k <- sort(unique(a))[k]
    total[1+k] <- sum(delta.mat[a==a_k,a==a_k])
    size[1+k] <- sum(a==a_k)^2
    deltamean[1+k] <- total[1+k]/size[1+k]
  }
  deltamean[1] <- (sum(delta.mat) - sum(total[-1]))/(n*n - sum(size[-1]))
  return(deltamean)
}

## Trial = split or relocate

Split <- function(a) {
  n <- length(a)
  ua <- unique(a)
  if(length(ua)==n) {
    return(a)
  } else{
    #ua.pick <- sample(x=ua,size=1)
    #a[names(sample(x=which(a==ua.pick),size=1))] <- max(ua)+1
    a.pick <- sample(x=a,size=1)
    pick.ind <- sample(x=which(a==a.pick),size=1)
    a[pick.ind] <- max(a)+1
    return(a)
  }
}

Relocate <- function(a){
  n <- length(a)
  ua <- unique(a)
  if(length(ua)==1) {
    return(a)
  } else{
    pick.ind <- sample(x=1:n,size=1)
    #a.pick <- a[pick.ind]
    #a[pick.ind] <- sample(x=a[-which(a==a.pick)],size=1)
    a[pick.ind] <- sample(x=a[-pick.ind],size=1)
    return(a)
  }
}


##########################
###### Scatterness ###
##########################


scatter = function(dat,cluster.assign,sil_cut=0.1){
  sd_check = apply(dat, 1, sd)
  
  if(ncol(dat) == 1|!all(sd_check != 0) ){
    dist.dat = as.matrix(dist(dat))
  }else{
    dist.dat = 1 - cor(t(dat))
  }
  sil = silhouette(cluster.assign, dist=dist.dat,diss=T)
  
  new.dist = dist.dat
  cluster.assign2 = cluster.assign
  
  while(min(sil[,3]) < sil_cut){
    cluster.assign2<-cluster.assign2[-which(sil[,3] == min(sil[,3]))]
    for (i in unique(cluster.assign2)){
      if(length(which(cluster.assign2==i))==1){
        cluster.assign2<-cluster.assign2[-which(cluster.assign2==i)]
      }
    }
    temp<-cluster.assign2
    for(d in 1:length(cluster.assign2)){
      cluster.assign2[d]<-rank(unique(temp))[which(unique(temp)==temp[d])]
    }#rename cluster index, so it is integer from 1 to k
    
    new.dist<-new.dist[rownames(new.dist)%in%names(cluster.assign2),
                       colnames(new.dist)%in%names(cluster.assign2)]
    sil <- silhouette(cluster.assign2, dist=new.dist, diss=T)#recalculate silhoutte
  }
  
  scatter.index = which(!names(cluster.assign)%in%names(cluster.assign2))
  return(scatter.index)
}

##########################
###### Text mining #######
##########################

TextMine <- function(hashtb, pathways, pathway, result, scatter.index=NULL,permutation="all"){
  cat("Performing Text Mining Analysis...\n")
  hashtbAll = hashtb
  k <- length(unique(result))
  if(!is.null(scatter.index)){
    cluster = result
    cluster[scatter.index] = k+1
    hashtb = hashtb[hashtb [,3]%in%which(pathways %in% pathway),]
    tmk = list()
    nperm = 1000
    if (nrow(hashtb) == 0){
      for (i in 1:(k-1)){
        tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)
      }
    }
    else{
      for (i in 1:k){
        e = cluster[cluster == i]
        e = which(pathways %in% names(e))
        hashcl = hashtb[hashtb [,3]%in%e,]
        hashcl = hashcl[duplicated(hashcl[,2]) | duplicated(hashcl[,2], fromLast=TRUE),]
        if (nrow(hashcl) != 0){
          hashf = hashcl
          hashf[,2] = 1
          hashf = aggregate(hashf[,c("row","value")],by = hashf["phrase"],FUN=sum)
          
          rownames(hashf) = hashf[,"phrase"]
          hashf = hashf[,-1]
          colnames(hashf) = c("count","sum")
          hashap = hashcl
          hashap[,c(2,3,4)] = 0
          mperm = matrix(nrow = nrow(hashf),ncol = nperm)
          for (j in 1:nperm){
            if (permutation=="all"){
              subtb = hashtbAll[hashtbAll[,3]%in%sample(1:length(pathways),length(e)),]
            }else if(permutation=="enriched"){
              subtb = hashtb[hashtb[,3]%in%sample(unique(hashtb[,3]),length(e)),]
            }else{
              stop("Permutation should be 'all' or 'enriched' ")
            }
            subtb = rbind(subtb,hashap)
            subtb = subtb[subtb$phrase %in% hashap$phrase,]
            subtb[,2] = 1
            subtb = aggregate(subtb[,c("row","value")],by = subtb["phrase"],FUN=sum)
            rownames(subtb) = subtb[,"phrase"]
            subtb = subtb[,-1]
            colnames(subtb) = c("count","sum")
            mperm[,j] = subtb[,2]
          }
          hashf[,"p-value"] = apply(cbind(hashf[,2],mperm),1,
                                    function(x)((nperm + 2)-rank(x)[1])/(nperm + 1))
          hashf[,"q-vlaue"] = p.adjust(hashf[,"p-value"],method = "BH")
          tmk[[i]] = hashf[order(hashf[,3],-hashf[,2]),]
        }
        else {tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)}
      }
      tmk[[k+1]] = matrix(NA,nrow = 1,ncol = 4)
    }
  }else{
    hashtb = hashtb[hashtb [,2]%in%which(pathways %in% pathway),]
    tmk = list()
    nperm = 1000
    if (nrow(hashtb) == 0){
      for (i in 1:(k-1)){
        tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)
      }
    }
    else{
      for (i in 1:k){
        e = cluster[cluster == i]
        e = which(pathways %in% names(e))
        hashcl = hashtb[hashtb [,3]%in%e,]
        hashcl = hashcl[duplicated(hashcl[,2]) | duplicated(hashcl[,2], fromLast=TRUE),]
        if (nrow(hashcl) != 0){
          hashf = hashcl
          hashf[,2] = 1
          hashf = aggregate(hashf[,c("row","value")],by = hashf["phrase"],FUN=sum)
          
          rownames(hashf) = hashf[,"phrase"]
          hashf = hashf[,-1]
          colnames(hashf) = c("count","sum")
          hashap = hashcl
          hashap[,c(2,3,4)] = 0
          mperm = matrix(nrow = nrow(hashf),ncol = nperm)
          for (j in 1:nperm){
            if (permutation=="all"){
              subtb = hashtbAll[hashtbAll[,3]%in%sample(1:length(pathways),length(e)),]
            }else if(permutation=="enriched"){
              subtb = hashtb[hashtb[,3]%in%sample(unique(hashtb[,3]),length(e)),]
            }else{
              stop("Permutation should be 'all' or 'enriched' ")
            }
            subtb = rbind(subtb,hashap)
            subtb = subtb[subtb$phrase %in% hashap$phrase,]
            subtb[,2] = 1
            subtb = aggregate(subtb[,c("row","value")],by = subtb["phrase"],FUN=sum)
            rownames(subtb) = subtb[,"phrase"]
            subtb = subtb[,-1]
            colnames(subtb) = c("count","sum")
            
            mperm[,j] = subtb[,2]
          }
          hashf[,"p-value"] = apply(cbind(hashf[,2],mperm),1,
                                    function(x)((nperm + 2)-rank(x)[1])/(nperm + 1))
          hashf[,"q-vlaue"] = p.adjust(hashf[,"p-value"],method = "BH")
          tmk[[i]] = hashf[order(hashf[,3],-hashf[,2]),]
        }
        else {tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)}
      }
    }
  }
  return(tmk)
} # End of Text Mining


writeTextOut <- function(tm_filtered,k,pathway.summary) {
  cat("Cluster 1\n", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
  cat("Key words,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
  write.table(t(rownames(tm_filtered[[1]])[1:15]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F,
              append = T, row.names=F,col.names=F,na="")
  cat("q_value,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
  write.table(t(tm_filtered[[1]][1:15,4]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F,
              append = T, row.names=F,col.names=F,na="")
  cat("count,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
  write.table(t(tm_filtered[[1]][1:15,1]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F,
              append = T, row.names=F,col.names=F,na="")
  write.table(pathway.summary[[1]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T,
              append = T, row.names=F,col.names=F)
  for (i in 2:k){
    cat(paste("\nCluster ", i, "\n", sep = ""), file = paste("Clustering_Summary_K",k,".csv",sep=""), append = T)
    cat("Key words,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(rownames(tm_filtered[[i]])[1:15]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F,
                append = T, row.names=F,col.names=F,na="")
    cat("q_value,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(tm_filtered[[i]][1:15,4]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F,
                append = T, row.names=F,col.names=F,na="")
    cat("count,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(tm_filtered[[i]][1:15,1]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F,
                append = T, row.names=F,col.names=F,na="")
    write.table(pathway.summary[[i]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T,
                append = T, row.names=F,col.names=F)
    
  }
}

writeTextOut <- function(tm_filtered,k,pathway.summary,scatter.index=NULL) {
  if(is.null(dim(tm_filtered[[1]]))==TRUE|dim(tm_filtered[[1]])[1] == 0){
    print(paste("No phrase pass q-value threshold in cluster 1"))
    cat("Cluster 1\n", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(pathway.summary[[1]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T, append = T, row.names=F,col.names=F)
  } else {
    cat("Cluster 1\n", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    cat("Key words,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(as.character(rownames(tm_filtered[[1]])[1:15])), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("q_value,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(tm_filtered[[1]][1:15,4]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("count,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(tm_filtered[[1]][1:15,1]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    write.table(pathway.summary[[1]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T, append = T, row.names=F,col.names=F)
    
  }
  
  for (i in 2:k){
    if(is.null(dim(tm_filtered[[i]]))==TRUE|dim(tm_filtered[[i]])[1] == 0){
      print(paste("No phrase pass q-value threshold in cluster ",k,sep = ""))
      cat(paste("\nCluster ", i, "\n", sep = ""), file = paste("Clustering_Summary_K",k,".csv",sep=""), append = T)
      write.table(pathway.summary[[i]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T, append = T, row.names=F,col.names=F)
    } else {
      cat(paste("\nCluster ", i, "\n", sep = ""), file = paste("Clustering_Summary_K",k,".csv",sep=""), append = T)
      cat("Key words,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
      write.table(t(as.character(rownames(tm_filtered[[i]])[1:15])), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("q_value,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
      write.table(t(tm_filtered[[i]][1:15,4]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("count,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
      write.table(t(tm_filtered[[i]][1:15,1]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      write.table(pathway.summary[[i]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T, append = T, row.names=F,col.names=F)
    }
    
  }
  if(!is.null(scatter.index)){
    cat(paste("\nSingleton Term", "\n", sep = ""), file = paste("Clustering_Summary_K",k,".csv",sep=""), append = T)
    write.table(pathway.summary[[k+1]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T,
                append = T, row.names=F,col.names=F)
  }
}

textMine <- function(hashtb,pathways,cluster.assign,scatter.index=NULL,thres=0.05,permutation="all"){
  tmk <- TextMine(hashtb=hashtb, pathways= pathways,
                  pathway=names(cluster.assign), result=cluster.assign, scatter.index,permutation=permutation)
  C <- length(unique(cluster.assign))
  tm_filtered <- list()
  for (i in 1:C){
    tm_filtered[[i]] <- tmk[[i]][which((as.numeric(tmk[[i]][,4]) < thres)), ]
  }
  if(!is.null(scatter.index)){
    tm_filtered[[C+1]] = tmk[[C+1]]
    cluster = cluster.assign
    cluster[scatter.index] = C+1
    pathway.summary <- lapply(1:(C+1), function(x) names(which(cluster==x)))
    writeTextOut(tm_filtered,C,pathway.summary,scatter.index=scatter.index)
  }else{
    pathway.summary <- lapply(1:C, function(x) names(which(cluster.assign==x)))
    writeTextOut(tm_filtered,C,pathway.summary)
  }
  return(tm_filtered)
}

##########################
##### Comember Plot ######
##########################
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    if(high-low != 0){
      x <- (x - low)/(high - low)
    }else{
      if(is.matrix(x)){
        x = matrix(data=high, nrow = nrow(x), ncol = ncol(x))
      }else{
        x = rep(high, length(x))
      }
    }
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

##########################
######    parseXML     ###
##########################
parseRelation <- function(pathwayID, keggSpecies="hsa", binary = T, sep = "-") {
  # download xml file
  download.kegg(pathway.id = pathwayID, keggSpecies, kegg.dir = ".", file.type="xml")
  # generate relation matrix
  KEGG.pathID2name = lapply(KEGGREST::keggList("pathway",keggSpecies),function(x) strsplit(x," - ")[[1]][-length(strsplit(x," - ")[[1]])])
  names(KEGG.pathID2name) = gsub(paste0("path:",keggSpecies),"",names(KEGG.pathID2name))
  
  pathName = unlist(KEGG.pathID2name[pathwayID])
  
  xmlFile = paste0(getwd(), "/",keggSpecies, pathwayID,".xml")
  pathway = KEGGgraph::parseKGML(xmlFile)
  pathway = KEGGgraph::splitKEGGgroup(pathway)
  
  entries = KEGGgraph::nodes(pathway)
  types = sapply(entries, KEGGgraph::getType)
  relations = unique(KEGGgraph::edges(pathway)) ## to avoid duplicated edges
  relationNum = length(relations)
  entryNames = as.list(sapply(entries, KEGGgraph::getName))
  if(any(types == "group") || any(types=="map")){
    entryNames = entryNames[!(types %in% c("group","map"))]
  }
  entryIds = names(entryNames)
  entryNames = lapply(1:length(entryNames), function(i) paste(entryNames[[i]],collapse=sep))
  names(entryNames) = entryIds
  
  entryNames.unique = unique(entryNames)
  entryNum = length(entryNames.unique)
  
  relation.mat = matrix(0, entryNum, entryNum)
  rownames(relation.mat) = colnames(relation.mat) = entryNames.unique
  
  ## if no relation edge, just return
  if(relationNum == 0){
    print(paste0("There is no topological connected gene nodes in ", pathName))
    return(relation.mat)
  }
  
  entry1 = KEGGgraph::getEntryID(relations)[,1]
  entry2 = KEGGgraph::getEntryID(relations)[,2]
  for(i in 1:length(relations)){
    if(entry1[i] %in% names(entryNames) && entry2[i] %in% names(entryNames)){
      relation.mat[entryNames[[entry1[i]]],entryNames[[entry2[i]]]]=1
    }
    else{
      print(paste("connections not included:",entry1[i], entry2[i], sep=" "))
    }
  }
  
  file.remove(xmlFile)
  return(relation.mat)
}
##########################
###### KEGG module SA ####
##########################
SA_module_M = function(sp.mat, xmlG, M, nodes, B = 1000,
                       G.ini.list=NULL, reps_eachM = 100,topG_from_previous=10,
                       Tm0=10,mu=0.95,epsilon=1e-5,
                       N=1000,run=10000,seed=12345,sub.num=1){
  #Null distribution for M
  set.seed(seed)
  null.sp.dist = rep(NA,B)
  for(b in 1:B){
    permute.set <- sample(xmlG,M)
    permute.mat <- sp.mat[match(permute.set,row.names(sp.mat)),
                          match(permute.set,colnames(sp.mat))]
    null.sp.dist[b] <- mean(c(permute.mat[lower.tri(permute.mat)]))
  }
  null.sp.mean = mean(null.sp.dist)
  null.sp.median = median(null.sp.dist)
  if(is.null(G.ini.list)){
    p.mean.ls = c()
    G.module.ls = list()
    SP.ls = c()
    for (l in 1:reps_eachM) {
      ##Initialize
      if(length(nodes) == M){
        G.module = nodes
      }else{
        G.module = sample(nodes, M)
      }
      SPc = avgSP(G.module, sp.mat)
      r = 0
      count = 0
      Tm = Tm0
      while((length(nodes)>M) & (r < run) & (count < N) & (Tm >= epsilon)) {
        #pi = exp(-GPc/Tm) ## Boltzmann dist #may need a different Tm or -logP to be comparable?
        #print(SPc)
        ##New trial
        r = r+1
        a.node = sample(setdiff(nodes,G.module),sub.num)
        G.module_new = G.module
        G.module_new[sample(M,sub.num)] = a.node
        
        SPn = avgSP(G.module_new, sp.mat)
        
        if(SPn < SPc | SPc == Inf) {
          ##accept
          SPc = SPn;
          G.module = G.module_new;
        }else{
          count = count + 1;
          p = exp((SPc-SPn)/Tm)
          #print(p)
          r = min(1,p); ## acceptance prob.
          u <- runif(1);
          if(u>r) {
            ##not accept
            Tm <- Tm*mu
          } else {
            SPc = SPn;
            G.module = G.module_new;
          }
        }
      }
      G.module.ls[[l]] = G.module
      p.mean.ls[l] = (sum(null.sp.dist<= SPc) + 1)/(B+1)
      SP.ls[l] = SPc
    }
  }else{
    each.times = round(reps_eachM/topG_from_previous)
    case.index = expand.grid(1:length(G.ini.list),1:each.times)
    p.mean.ls = c()
    G.module.ls = list()
    SP.ls = c()
    for (l in 1:nrow(case.index)) {
      G.ini = G.ini.list[[case.index[l,1]]]
      ##Initialize
      if(length(nodes) == M){
        G.module = nodes
      }else{
        x = M-length(G.ini)
        G.module = c(G.ini,sample(setdiff(nodes,G.ini),x))
      }
      SPc = avgSP(G.module, sp.mat)
      r = 0
      count = 0
      Tm = Tm0
      while((length(nodes)>M) & (r < run) & (count < N) & (Tm >= epsilon)) {
        #pi = exp(-GPc/Tm) ## Boltzmann dist #may need a different Tm or -logP to be comparable?
        #print(SPc)
        ##New trial
        r = r+1
        a.node = sample(setdiff(nodes,G.module),sub.num)
        G.module_new = G.module
        G.module_new[sample(M,sub.num)] = a.node
        
        SPn = avgSP(G.module_new, sp.mat)
        
        if(SPn < SPc | SPc == Inf) {
          ##accept
          SPc = SPn;
          G.module = G.module_new;
        }else{
          count = count + 1;
          p = exp((SPc-SPn)/Tm)
          #print(p)
          r = min(1,p); ## acceptance prob.
          u <- runif(1);
          if(u>r) {
            ##not accept
            Tm <- Tm*mu
          } else {
            SPc = SPn;
            G.module = G.module_new;
          }
        }
      }
      G.module.ls[[l]] = G.module
      p.mean.ls[l] = (sum(null.sp.dist<= SPc) + 1)/(B+1)
      SP.ls[l] = SPc
    }
    
  }
  p.sd.ls = sqrt(p.mean.ls*(1-p.mean.ls)/B)
  r.p = rank(p.mean.ls,ties.method = "random")
  index = match(1:topG_from_previous,r.p)
  best.index = which(r.p == 1)
  
  minG = G.module.ls[[best.index]]
  sp = SP.ls[best.index]
  p.mean = p.mean.ls[best.index]
  p.sd = p.sd.ls[best.index]
  
  top.G = G.module.ls[index]
  top.pmean = p.mean.ls[index]
  top.psd = p.sd.ls[index]
  top.sp = SP.ls[index]
  
  return(list(minG = minG,sp = sp,p.mean = p.mean,p.sd = p.sd,
              top.G = top.G,top.sp = top.sp,top.pmean = top.pmean,top.psd = top.psd,
              null.sp.mean = null.sp.mean,null.sp.median = null.sp.median))
}

avgSP = function(G.module, sp.mat){
  m = length(G.module)
  set.mat = sp.mat[match(G.module,row.names(sp.mat)),
                   match(G.module,colnames(sp.mat))]
  G.sp = mean(c(set.mat[lower.tri(set.mat)]))
  return(G.sp)
}

KEGG_module_topology_plot = function(res_KEGG_module,which_to_draw = "all",filePath = getwd()){
  minG.ls = res_KEGG_module$minG.ls
  module.size = as.numeric(gsub("minG","",names(minG.ls)))
  module.type = res_KEGG_module$module.type
  
  mergePMmat = res_KEGG_module$mergePMmat
  KEGGspecies = res_KEGG_module$KEGGspecies
  KEGGpathwayID = res_KEGG_module$KEGGpathwayID
  KEGGpathwayID_spec = paste0(KEGGspecies,KEGGpathwayID)
  data.pair = res_KEGG_module$data.pair
  dat1.name = data.pair[1]
  dat2.name = data.pair[2]
  
  if(which_to_draw[[1]] == "all"){
    which_to_draw_index = 1:length(minG.ls)
  }else if(!is.numeric(which_to_draw)){
    stop("which_to_draw should be 'all' or a numeric vector")
  }else{
    which_to_draw_index = match(which_to_draw,module.size)
  }
  orig.path = getwd()
  setwd(filePath)
  for (j in 1:length(which_to_draw_index)) {
    index = which_to_draw_index[j]
    topologyG = minG.ls[[index]]$minG
    if(is.null(dim(topologyG))){
      signPM.mat = mergePMmat[topologyG,]
      row.names(signPM.mat) = sapply(row.names(signPM.mat), function(x) strsplit(x,"_")[[1]][1])
      
      res = pathview(gene.data = signPM.mat, pathway.id = KEGGpathwayID,
                     species = KEGGspecies, out.suffix = "", kegg.native = T,
                     key.pos = "bottomright", map.null=T,cex = 0.15)
      
      file.rename(paste(KEGGpathwayID_spec,"..multi.png",sep=""),
                  paste(KEGGpathwayID_spec,"_",dat1.name,"_",dat2.name,"_",names(minG.ls)[index],"_",module.type,".png",sep=""))
      file.remove(paste(KEGGpathwayID_spec,".xml",sep=""))
      file.remove(paste(KEGGpathwayID_spec,".png",sep=""))
      
      
    }else{
      for (i in 1:ncol(topologyG)) {
        topologyG0 = topologyG[,i]
        signPM.mat = mergePMmat[topologyG0,]
        row.names(signPM.mat) = sapply(row.names(signPM.mat), function(x) strsplit(x,"_")[[1]][1])
        
        res = pathview(gene.data = signPM.mat, pathway.id = KEGGpathwayID,
                       species = KEGGspecies, out.suffix = "", kegg.native = T,
                       key.pos = "bottomright", map.null=T,cex = 0.15)
        file.rename(paste(KEGGpathwayID_spec,"..multi.png",sep=""),
                    paste(KEGGpathwayID_spec,"_",dat1.name,"_",dat2.name,"_",names(minG.ls)[index],"_",i,"_",module.type,".png",sep=""))
        file.remove(paste(KEGGpathwayID_spec,".xml",sep=""))
        file.remove(paste(KEGGpathwayID_spec,".png",sep=""))
        
      }
    }
  }
  setwd(orig.path)
}

KEGG_module = function(mcmc.merge.list,dataset.names,
                       KEGGspecies="hsa",
                       KEGGpathwayID,
                       KEGG.dataGisTopologyG = FALSE,
                       KEGG.dataG2topologyG = NULL,
                       data.pair,
                       gene_type = c("discordant","concordant"),
                       DE_PM_cut = 0.2, minM = 4, maxM = NULL,
                       B = 1000, cores = 1,
                       search_method = c("Exhaustive","SA"),
                       reps_eachM = 1,
                       topG_from_previous=1,
                       Tm0=10,mu=0.95,epsilon=1e-5,N=3000,
                       Elbow_plot = T, filePath = getwd(),
                       seed = 12345, sep = "-"){
  if(minM < 2){
    stop("minM has to be larger than 1.")
  }
  
  dat1.name = data.pair[[1]]
  dat2.name = data.pair[[2]]
  dat1 = mcmc.merge.list[[match(dat1.name,dataset.names)]]
  dat2 = mcmc.merge.list[[match(dat2.name,dataset.names)]]
  signPM.mat = cbind(apply(dat1,1,mean),apply(dat2,1,mean))
  
  #match data names and gene names on KEGG topology
  if(KEGG.dataGisTopologyG == TRUE){
    topologyG = rownames(signPM.mat)
  }else if(!is.null(KEGG.dataG2topologyG)){
    topologyG = KEGG.dataG2topologyG[match(rownames(signPM.mat),KEGG.dataG2topologyG[,1]),2]
    na.index = which(is.na(topologyG))
    if(length(na.index != 0)){
      topologyG = topologyG[-na.index]
      signPM.mat = signPM.mat[-na.index,]
    }
    
  }else if(KEGGspecies == "hsa"){
    map.ls = as.list(org.Hs.eg.db::org.Hs.egALIAS2EG)
    topologyG = sapply(rownames(signPM.mat),function(g) map.ls[[g]][[1]])
    na.index = sapply(topologyG, is.null)
    if(sum(na.index) != 0){
      topologyG = topologyG[!na.index]
      signPM.mat = signPM.mat[!na.index,]
    }
    
  }else if(KEGGspecies == "mmu"){
    map.ls = as.list(org.Mm.eg.db::org.Mm.egALIAS2EG)
    topologyG = sapply(rownames(signPM.mat),function(g) map.ls[[g]][[1]])
    na.index = sapply(topologyG, is.null)
    if(sum(na.index) != 0){
      topologyG = topologyG[!na.index]
      signPM.mat = signPM.mat[!na.index,]
    }
    
  }else if(KEGGspecies == "rno"){
    map.ls = as.list(as.list(org.Rn.eg.db::org.Rn.egALIAS2EG))
    topologyG = sapply(rownames(signPM.mat),function(g) map.ls[[g]][[1]])
    na.index = sapply(topologyG, is.null)
    if(sum(na.index) != 0){
      topologyG = topologyG[!na.index]
      signPM.mat = signPM.mat[!na.index,]
    }
    
  }else if(KEGGspecies == "cel"){
    wormbase = biomaRt::useMart(biomart = "parasite_mart",
                                host = "https://parasite.wormbase.org",
                                port = 443)
    wormbase = useDataset(mart = wormbase, dataset = "wbps_gene")
    map.mat = getBM(attributes = c("entrezgene_name","wormbase_gseq"),
                    filters = "entrezgene_name",
                    values = rownames(signPM.mat),
                    mart = wormbase)
    topologyG = paste0("CELE_",map.mat[match(rownames(signPM.mat),map.mat[,"entrezgene_name"]),"wormbase_gseq"])
    
    na.index = which(topologyG == "CELE_NA")
    if(length(na.index != 0)){
      topologyG = topologyG[-na.index]
      signPM.mat = signPM.mat[-na.index,]
    }
  }else if(KEGGspecies == "dme"){
    map.ls0 = as.list(org.Dm.eg.db::org.Dm.egALIAS2EG)
    EntrezID = sapply(rownames(signPM.mat),function(g) map.ls0[[g]][[1]])
    map.ls = as.list(org.Dm.eg.db::org.Dm.egFLYBASECG)
    topologyG = paste0("Dmel_",sapply(EntrezID,function(g) ifelse(is.null(g), NA, map.ls[[g]][[1]])))
    
    na.index = which(topologyG == "Dmel_NA")
    if(length(na.index != 0)){
      topologyG = topologyG[-na.index]
      signPM.mat = signPM.mat[-na.index,]
    }
  }else{
    stop("Please provide mapping between data genes and topology genes when they are of different gene name types and species is not one of 'hsa', 'mmu','rno','cel'or'dme'.")
  }
  
  rownames(signPM.mat) = topologyG
  adjacent_mat = parseRelation(pathwayID = KEGGpathwayID, keggSpecies = KEGGspecies, sep = sep)
  xmlG = row.names(adjacent_mat)[grep(KEGGspecies,row.names(adjacent_mat))]
  xmlG = gsub(paste0(KEGGspecies,":"),"",xmlG)
  xmlG.ls = lapply(xmlG, function(x){
    strsplit(x,sep)[[1]]
  })
  
  row.names(adjacent_mat) = colnames(adjacent_mat) = gsub(paste0(KEGGspecies,":"),"",row.names(adjacent_mat))
  
  mergePMls = lapply(1:length(xmlG.ls), function(x){
    genes = xmlG.ls[[x]]
    cmG = intersect(topologyG,genes)
    if(length(cmG) !=0){
      sub.signPM.mat = matrix(signPM.mat[cmG,],ncol = 2)
      avgPM = apply(sub.signPM.mat, 2, mean)
      return(avgPM)
    }
  })
  names(mergePMls) = xmlG
  mergePMmat = do.call(rbind,mergePMls)
  
  #discordant/concordant genes definition
  #discordant/concordant genes definition
  if(all(abs(mergePMmat[,1])<=DE_PM_cut | abs(mergePMmat[,2])<=DE_PM_cut)){
    DE_PM_cut = -1
    print(paste0("All genes with DE strength greater than the cutoff value for posterior probability of DE are not connected. Removed the cutoff criterion to consider all ",gene_type," genes regardless of its DE stength."))
  }
  if(gene_type == "discordant"){
    
    if(sum(mergePMmat[,1]*mergePMmat[,2]<0) == 0){
      stop("No discordant genes are topologically connected. Probably dut to low DE strength in one study. Please check the genePM plot. ")
    }else{
      nodes = unique(rownames(mergePMmat)[which(mergePMmat[,1]*mergePMmat[,2]<0&abs(mergePMmat[,1])>DE_PM_cut&abs(mergePMmat[,2])>DE_PM_cut)])
    }
    
  }else if(gene_type == "concordant"){
    
    if(sum(mergePMmat[,1]*mergePMmat[,2]<0) == 0){
      stop("No concordant genes are topologically connected. Probably dut to low DE strength in one study. Please check the genePM plot. ")
    }else{
      nodes = unique(rownames(mergePMmat)[which(mergePMmat[,1]*mergePMmat[,2]>0&abs(mergePMmat[,1])>DE_PM_cut&abs(mergePMmat[,2])>DE_PM_cut)])
    }
  }else{
    stop("gene_type has to be 'discordant' or 'concordant'.")
  }
  
  undir_adj_mat = adjacent_mat
  for(i in 1:nrow(adjacent_mat)){
    for (j in 1:ncol(adjacent_mat)) {
      undir_adj_mat[i,j] = max(adjacent_mat[i,j],adjacent_mat[j,i])
      undir_adj_mat[j,i] = max(adjacent_mat[i,j],adjacent_mat[j,i])
    }
  }
  g = graph_from_adjacency_matrix(undir_adj_mat,mode="undirected")
  sp.mat <- shortest.paths(g)
  
  d = degree(g)
  sort.d = sort(d,decreasing = T)[nodes]
  
  sub.adj.mat = undir_adj_mat[match(nodes,row.names(undir_adj_mat)),
                              match(nodes,colnames(undir_adj_mat))]
  #dim(sub.adj.mat)
  #dim(undir_adj_mat)
  #sum(lower.tri(sub.adj.mat))/length(lower.tri(sub.adj.mat))
  #(sum(lower.tri(undir_adj_mat))-sum(lower.tri(sub.adj.mat)))/(length(lower.tri(undir_adj_mat))-length(lower.tri(sub.adj.mat)))
  
  if(is.null(maxM)){
    maxM = length(nodes)
  }else{
    maxM = min(length(nodes),maxM)
  }
  module.size = minM:maxM
  
  search_method = match.arg(search_method)
  if(search_method == "Exhaustive"){
    minG.ls = mclapply(1:length(module.size),function(i){
      m = module.size[i]
      m.combn = combn(x=nodes,m=m)
      
      ## observed ones:
      
      asp.m = rep(NA,ncol(m.combn))
      
      
      for(j in 1:ncol(m.combn)){
        node.set <- m.combn[,j]
        set.mat <- sp.mat[match(node.set,row.names(sp.mat)),
                          match(node.set,colnames(sp.mat))]
        asp.m[j] <- mean(c(set.mat[lower.tri(set.mat)]))
      }
      
      (minG = m.combn[,which(asp.m == min(asp.m))])
      
      set.seed(seed)
      asp.m.perm = rep(NA,B)
      
      for(b in 1:B){
        permute.set <- sample(xmlG,m)
        permute.mat <- sp.mat[match(permute.set,row.names(sp.mat)),
                              match(permute.set,colnames(sp.mat))]
        asp.m.perm[b] <- mean(c(permute.mat[lower.tri(permute.mat)]))
        
      }
      
      p.observed = (sum(asp.m.perm<= min(asp.m)) + 1)/(B+1)
      p.sd = sqrt(p.observed*(1-p.observed)/B)
      null.sp.mean = mean(asp.m.perm)
      null.sp.median = median(asp.m.perm)
      
      return(list(minG = minG, p.mean = p.observed, p.sd = p.sd,sp = min(asp.m),
                  null.sp.mean = null.sp.mean,
                  null.sp.median = null.sp.median))
    },mc.cores = cores)
  }else{
    minG.ls = list()
    for(i in 1:length(module.size)){
      M = module.size[i]
      print(M)
      if(M>min(module.size)){
        G.ini.list = minG.ls[[i-1]]$top.G
        minG.ls[[i]] = SA_module_M(sp.mat, xmlG, M, nodes, B = B,
                                   G.ini.list=G.ini.list, reps_eachM = reps_eachM,
                                   topG_from_previous=topG_from_previous,
                                   Tm0=Tm0,mu=mu,epsilon=epsilon,
                                   N=N,seed=seed)
      }else{
        minG.ls[[i]] = SA_module_M(sp.mat, xmlG, M, nodes, B = B,
                                   G.ini.list=NULL, reps_eachM = reps_eachM,
                                   topG_from_previous=topG_from_previous,
                                   Tm0=Tm0,mu=mu,epsilon=epsilon,
                                   N=N,seed=seed)
      }
    }
  }
  names(minG.ls) = paste0("minG",module.size)
  
  #Select best size
  p.mean = sapply(minG.ls, function(x) x[["p.mean"]])
  p.sd = sapply(minG.ls, function(x) x[["p.sd"]])
  #obs.mean.ratio = sapply(minG.ls, function(x) x[["sp"]])/sapply(minG.ls, function(x) x[["null.sp.mean"]])
  obs.median.ratio = sapply(minG.ls, function(x) x[["sp"]])/sapply(minG.ls, function(x) x[["null.sp.median"]])
  obs.median.ratio[is.na(obs.median.ratio)] = 0
  obs.sp = sapply(minG.ls, function(x) x[["sp"]])
  
  index = which.min(p.mean)
  p.cut = p.mean[index]+2*p.sd[index]
  finalSelect = names(minG.ls)[max(which(p.mean<p.cut & 1:length(p.mean)>=index))]
  
  KEGGpathwayID_spec = paste0(KEGGspecies,KEGGpathwayID)
  if(Elbow_plot == T){
    names(p.mean) = names(p.sd) = module.size
    pL = sapply(p.mean-2*p.sd, function(x) ifelse(x<0,0,x))
    df = data.frame(logp.observed = -log10(p.mean),
                    logp.max = -log10(pL),
                    logp.min = -log10(p.mean+2*p.sd),
                    size = module.size)
    png(paste0(filePath,"/",KEGGpathwayID_spec,"_",gene_type,"_",dat1.name,"_",dat2.name,"_",search_method,"_elbow_plot.png"))
    p = ggplot(df, aes(x=size, y=logp.observed)) +
      geom_errorbar(aes(ymin=logp.min, ymax=logp.max), width=.1) +
      geom_line() +
      geom_point() +
      ylim(0,3.1)+
      labs(title = paste0(KEGGpathwayID_spec,"_",gene_type,"_",dat1.name,"_",dat2.name),y = "-log10(p-value)")
    print(p)
    dev.off()
    
    df.ratio = data.frame(obs.sp, obs.median.ratio, size = module.size)
    
    # png(paste0(filePath,"/",KEGGpathwayID_spec,"_",gene_type,"_",dat1.name,"_",dat2.name,"_",search_method,"_avgSP.png"))
    # p = ggplot(df.ratio, aes(x=size, y=obs.sp)) +
    #   geom_line() +
    #   geom_point()+
    #   labs(title = paste0(KEGGpathwayID_spec,"_",gene_type,"_",dat1.name,"_",dat2.name),y = "module average shortest path value")
    # print(p)
    # dev.off()
    
    # png(paste0(filePath,"/",KEGGpathwayID_spec,"_",gene_type,"_",dat1.name,"_",dat2.name,"_",search_method,"avgSP_ratio_obs_to_median.png"))
    # p = ggplot(df.ratio, aes(x=size, y=obs.median.ratio)) +
    #   geom_line() +
    #   geom_point() +
    #   labs(title = paste0(KEGGpathwayID_spec,"_",gene_type,"_",dat1.name,"_",dat2.name),y = "average module shortest path/median(null average shortest path)")
    # print(p)
    # dev.off()
    
  }
  return(list(minG.ls=minG.ls,bestSize = finalSelect,
              mergePMmat = mergePMmat,
              KEGGspecies = KEGGspecies,
              KEGGpathwayID = KEGGpathwayID,
              data.pair = data.pair,
              module.type = gene_type))
}
