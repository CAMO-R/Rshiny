##' Analysis results for single pair: visualization outputs for each pathway
##' The \code{singleOutput} is function to generate visualization outputs 
##' for single pair: including heatmap of gene posterior mean, kegg 
##' pathway topology for each pathway
##' @title Analysis results for single pair: visualization outputs  
##' @param mcmc.merge.list: a list of merged MCMC output matrices.
##' @param select.pathway.list: a list of selected pathways (containing gene 
##' components).
##' @param dataset.names: a vector of dataset names.
##' @param output: two options: "genePM" (generating heatmap of gene 
##' posterior mean),"keggView" (generating kegg pathway topology, human KEGG only). 
##' cannot be empty.
##' @param kegg_pathname: KEGG pathway name list. For "keggView" only.
##' @param hs_gene_id: Human sapiens gene id. For "keggView" only.

##' @return stored output in created folders.
##' @export
##' @examples
##' \dontrun{
##' #mcmc.merge.list from the merge step
##' #select.pathway.list from the pathSelect step
##' dataset.names = c("hb","mb")
##' library(KEGG.db)
##' kegg_pathname <- unlist(as.list(KEGGPATHID2NAME))
##' library("org.Hs.eg.db")
##' hs_gene_id <- unlist(mget(x=rownames(mcmc.merge.list[[1]]),
##' envir=org.Hs.egALIAS2EG))
##' singleOutput(mcmc.merge.list,dataset.names,select.pathway.list,
##' output=c("genePM","keggView"))
##' }

singleOutput <- function(mcmc.merge.list,dataset.names,select.pathway.list,
                         output=c("genePM","keggView"),
                         kegg_pathname=NULL,hs_gene_id=NULL) {
  
### Single pair output: including per pathway gene PM heatmap, pathway topology (KEGG only, human only)

  #kegg_pathname, hs_gene_id are local input
  
   if(length(output)==0 || is.null(output)){
     stop("at least one type of output has to be chosen")
   }
  
   orig.path <- getwd()
   pathway.name <- names(select.pathway.list)
   K <- length(pathway.name)
   
   dat1 <- mcmc.merge.list[[1]]
   dat2 <- mcmc.merge.list[[2]]
   
   if("genePM" %in% output){
     dir.path <- "genePM"
     if (!file.exists(dir.path)) dir.create(dir.path)
     setwd(paste(orig.path,"/",dir.path,sep=""))
     for(k in 1:K){
       print(paste("genePM",k,sep=":"))
       pathk.name <- pathway.name[k]
       pathway.genes <- select.pathway.list[[k]]
       signPM.list <- list(apply(dat1,1,mean),apply(dat2,1,mean))
       names(signPM.list) <- dataset.names
       hm <- genePM(signPM.list, pathway.genes=pathway.genes, 
                          pathway.name=pathk.name)
     }
   }
   
   setwd(orig.path)
   
   if("keggView" %in% output){
     if(sum(grepl("KEGG",pathway.name))==0) {
       warning("No KEGG pathways") 
      } else{
       data("paths.hsa")
       dir.path <- "keggView"
       if (!file.exists(dir.path)) dir.create(dir.path)
       setwd(paste(orig.path,"/",dir.path,sep=""))
       
       kegg.pathway.name <- pathway.name[grep("KEGG",pathway.name)]
       K_KEGG <- length(kegg.pathway.name)

       for(k in 1:K_KEGG){
         print(paste("keggView",k,sep=":"))
         keggk.name <- kegg.pathway.name[k]
         overlap.genes <- intersect(rownames(dat1),select.pathway.list[[keggk.name]])
         signPM.mat <- cbind(apply(dat1[overlap.genes,],1,mean),
                             apply(dat2[overlap.genes,],1,mean))
         colnames(signPM.mat) <- dataset.names
         keggk.name1 <- gsub("KEGG ","",keggk.name) 
         pathwayID <- gsub("hsa","",names(paths.hsa)[which(paths.hsa==keggk.name1)])
         res <- keggView(mat=signPM.mat,pathwayID)
         if(grepl("/",keggk.name)){
           keggk.name <- gsub("/","-",keggk.name)
         }
         hsaName <- paste("hsa",pathwayID,sep="")
         file.rename(paste(hsaName,"..multi.png",sep=""), 
                     paste(keggk.name,".png",sep=""))
         file.remove(paste(hsaName,".xml",sep=""))
         file.remove(paste(hsaName,".png",sep=""))
       }
     } 
   } 
  setwd(orig.path)
  print("Single pair analysis completed.") 
}


# 
# mdsModel <- function(arsPath,model.name,pathway.name,sep) {
#   ## for each pathway, plot MDS for all models
#   ## arsP is a vector of choose(M,2) elements (named in paste(name1,name2,sep="")) 
#   ## model.name is a vector of model names
#   M <- length(model.name)
#   distF <- ARStransform(arsPath)
#   d <- matrix(NA,nrow=M,ncol=M)
#   rownames(d) <- colnames(d) <- model.name
#   
#   for(i in 1:(M-1)){
#     for(j in (i+1):M){
#       name1 <- rownames(d)[i]
#       name2 <- rownames(d)[j]
#       d[name1,name2] <- d[name2,name1] <- distF[paste(name1,name2,sep=sep)]
#     }
#   }
#   
#   diag(d) <- 0
#   dist <- as.dist(d,upper = TRUE, diag = TRUE)
#   fit <- sammon(d=dist, y= jitter(cmdscale(dist, 2)), k=2) # k is the number of dim
#   
#   x <- fit$points[,1]
#   y <- fit$points[,2]
#   xlimit <- ifelse(abs(min(x))>abs(max(x)),abs(min(x)),abs(max(x)))
#   ylimit <- ifelse(abs(min(y))>abs(max(y)),abs(min(y)),abs(max(y)))
#   
#   if(grepl("/",pathway.name)){
#     pathway.name <- gsub("/","-",pathway.name)
#   }
#   
#   color <- rainbow(M,s=0.5,v=1,alpha=1)
#   png(paste(pathway.name,".png",sep=""))
#   p<-ggplot() +
#     ggtitle(pathway.name) +
#     xlab("Coordinate 1") + ylab("Coordinate 2") + 
#     xlim(c(-xlimit-0.5,xlimit+0.5)) + ylim(c(-ylimit-0.5,ylimit+0.5)) + 
#     geom_point(aes(x, y), color = color  ,size=6) +
#     geom_text_repel(aes(x, y, label = rownames(d),fontface="bold"),size=8) + 
#     theme(plot.title = element_text(size = 15, hjust=0.5,face="bold"),
#           axis.text.x = element_text(size = 12),
#           axis.text.y = element_text(size = 12))
#   print(p)
#   dev.off()
# }
# 
# 
# SA_algo <- function(arsPvaluePath,model.name,sep,Tm=10,P=0.5,
#                        mu=0.9,epsilon=1e-5,N=1000,seed=12345){
#   ## Total possible configurations = choose(n,1) + choose(n,2) + ...
#   ## clustering models using SA algorithm
#   ## delta.mat is a matrix of pairwise -log10(arsPvalue) of M rows and M columns 
#   ## with model names
#   M <- length(model.name)
#   distF <- -log10(arsPvaluePath)
#   delta.mat <- matrix(NA,nrow=M,ncol=M)
#   rownames(delta.mat) <- colnames(delta.mat) <- model.name 
#   
#   for(i in 1:(M-1)){
#     for(j in (i+1):M){
#       name1 <- rownames(delta.mat)[i]
#       name2 <- rownames(delta.mat)[j]
#       delta.mat[name1,name2] <- delta.mat[name2,name1] <- distF[paste(name1,name2,sep=sep)]
#     }
#   }
#   diag(delta.mat) <- max(delta.mat,na.rm=T)
#   
#   n <- nrow(delta.mat) # =M
#   obs.name <- rownames(delta.mat)
#   
#   ## initialize
#   hc <- hclust(d=dist((max(delta.mat)-delta.mat)^2))
#   K <- 3
#   a <- cutree(hc, k = K)
#   names(a) <- obs.name
#   delta.est <- Est_mean(delta.mat,a) 
#   names(delta.est) <- c(0,1:K)
#   
#   count <- 0 #initial
#   Jc <- E_tot(delta.mat,a,delta.est)
#   pi <- exp(-Jc/Tm) ## Boltzmann dist
#   
#   while((count < N) && (Tm >= epsilon)) {
#     ##New trial
#     u <- runif(1)
#     if(u>0.5){
#       a_new <- Split(a)
#     } else{
#       a_new <- Relocate(a)
#     }
#     names(a_new) <- obs.name
#     K_new <- length(unique(a_new))
#     delta.est_new <- Est_mean(delta.mat,a_new)
#     Jn <- E_tot(delta.mat, a_new, delta.est_new)
#     pi_new <- exp(-Jn/Tm)
#     
#     if(Jn < Jc) { 
#       ##accept
#       Jc <- Jn;
#       a <- a_new;
#       K <- K_new;
#       delta.est <- delta.est_new;
#     } else {
#       count <- count + 1;
#       r <- min(1,pi_new/pi); ## acceptance prob.
#       u <- runif(1);
#       if(u>r) {
#         ##not accept
#         Tm <- Tm*mu
#       } else {
#         Jc <- Jn;
#         a <- a_new;
#         K <- K_new;
#         delta.est <- delta.est_new;
#       } 
#     }
#   }
#   return(a) #cluster assigned 
# }
# 
# clustModel <- function(arsPvaluePath,model.name, cluster.assign,pathway.name,sep){
#   ## delta.mat is a matrix of pairwise -log(arsPvalue) of M rows and M columns 
#   ## with model names 
#   ## cluster.assign: results from SA algorithm
#   ## pathway.name: pathway name
#   
#   M <- length(model.name)
#   distF <- -log10(arsPvaluePath)
#   delta.mat <- matrix(NA,nrow=M,ncol=M)
#   rownames(delta.mat) <- colnames(delta.mat) <- model.name 
#   
#   for(i in 1:(M-1)){
#     for(j in (i+1):M){
#       name1 <- rownames(delta.mat)[i]
#       name2 <- rownames(delta.mat)[j]
#       delta.mat[name1,name2] <- delta.mat[name2,name1] <- distF[paste(name1,name2,sep=sep)]
#     }
#   }
#   diag(delta.mat) <- max(delta.mat,na.rm=T)
#   
#   if(grepl("/",pathway.name)){
#     pathway.name <- gsub("/","-",pathway.name)
#   }
#   
#   png(paste(pathway.name,'.png',sep="_"))
#   hm <- heatmap.2(delta.mat[order(cluster.assign),order(cluster.assign)], 
#                  main=pathway.name,
#                  cexCol=1,cexRow=1,
#                  colsep=cumsum(table(cluster.assign)),
#                  rowsep=cumsum(table(cluster.assign)),
#                  sepwidth=c(0.05, 0.05),  # width of the borders
#                  sepcolor=c('white'),
#                  symbreaks=T,key=T, keysize=1,symkey=F, 
#                  dendrogram=c('none'),density.info="none", 
#                  trace="none",Rowv=F,Colv=F,
#                  srtCol=50, symm=F,
#                  col=greenred,breaks=seq(0,round(max(delta.mat)),by=0.01) )
#   dev.off()
#   return(hm)
# }
# 
# 
# clustModelOne <- function(arsPvaluePath,model.name, pathway.name,sep){
#   ## delta.mat is a matrix of pairwise -log(arsPvalue) of M rows and M columns 
#   ## with model names 
#   ## cluster.assign: results from SA algorithm
#   ## pathway.name: pathway name
#   
#   M <- length(model.name)
#   distF <- -log10(arsPvaluePath)
#   delta.mat <- matrix(NA,nrow=M,ncol=M)
#   rownames(delta.mat) <- colnames(delta.mat) <- model.name 
#   
#   for(i in 1:(M-1)){
#     for(j in (i+1):M){
#       name1 <- rownames(delta.mat)[i]
#       name2 <- rownames(delta.mat)[j]
#       delta.mat[name1,name2] <- delta.mat[name2,name1] <- distF[paste(name1,name2,sep=sep)]
#     }
#   }
#   diag(delta.mat) <- max(delta.mat,na.rm=T)
#   
#   if(grepl("/",pathway.name)){
#     pathway.name <- gsub("/","-",pathway.name)
#   }
#   
#   png(paste(pathway.name,'.png',sep="_"))
#   hm <- heatmap.2(delta.mat, 
#                   main=pathway.name,
#                   cexCol=1,cexRow=1,
#                   symbreaks=T,key=T, keysize=1,symkey=F, 
#                   dendrogram=c('none'),density.info="none", 
#                   trace="none",Rowv=F,Colv=F,
#                   srtCol=50, symm=F,
#                   col=greenred,breaks=seq(0,round(max(delta.mat)),by=0.01) )
#   dev.off()
#   return(hm)
# }
# 
# 
# genePM <- function(signPM.list, pathway.genes, pathway.name){
#   
#   M <- length(signPM.list)
#   model.name <- names(signPM.list)
#   std.genes <- intersect(names(signPM.list[[1]]),pathway.genes) 
#   G <- length(std.genes)
#   
#   mat <- matrix(0,nrow=G,ncol=M)
#   rownames(mat) <- std.genes
#   colnames(mat) <- model.name
#   
#   for (m in 1:M){
#     mat[std.genes,m] <- signPM.list[[m]][std.genes]
#   }
#   
#   if(grepl("/",pathway.name)){
#     pathway.name <- gsub("/","-",pathway.name)
#   }
#   
#   #pdf(paste(pathway.name,'.pdf',sep=""))
#   jpeg(paste(pathway.name,".jpeg",sep=""),quality = 100)
#   par(cex.main=1)
#   hm <- heatmap.2(mat, symm=F,main=pathway.name,
#                 cexCol=1,cexRow=0.6,
#                 colsep=c(1:ncol(mat)),
#                 rowsep=c(1:nrow(mat)),
#                 sepwidth=c(0.01, 0.2),  # width of the borders
#                 sepcolor=c('black'),scale='none',
#                 symbreaks=T,key=T, keysize=1,symkey=F, 
#                 dendrogram=c('row'),density.info="none", 
#                 trace="none",Rowv=T,Colv=F,
#                 col=bluered,breaks=seq(-1,1,by=0.001),
#                 srtCol=0,adjCol = c(NA,0.5))
#   dev.off()
#   return(hm)
# }    
# 
# 
# keggView <- function(mat,pathwayID){
#   ## human only 
#   #set your working dir, automatically save there
#   #mat is two column signed PM
#   #kegg.name <- gsub("KEGG ","",pathway.name) 
#   library(org.Hs.eg.db)
#   #database = c("org.Hs.eg.db")
#   #also need library("biomaRt") library("KEGG.db")
#   
#   #require(database)
#   
#   model.name <- colnames(mat)
#   std.genes <- rownames(mat)
#   
#   out <- unlist(mget(x=std.genes,envir=org.Hs.egALIAS2EG))
#   IDs = out[std.genes]
#   geneDat <- mat
#   rownames(geneDat) <- IDs
#   
#   pv.out <- pathview(gene.data = geneDat, pathway.id = pathwayID, 
#                      species = "hsa", out.suffix = "", kegg.native = T,
#                      key.pos = "bottomright", map.null=T,cex = 0.15)
#   
#   return(pv.out)
#   
# }   
