### ACS-DE & ADS-DE ###
#------------function-----------------
ARS_to_size <- function(ARSp,factor=8){
  if(ARSp > 0.05) {
    return(1) 
  } else {
    return(-log10(ARSp)*factor)
  }
}

ACS_ADS_DE <- function(ds1,ds2,DEevid1,DEevid2,ACSp,ADSp,cluster,
                       highlight.pathways = NULL,lb=0,ub=1,size.scale=4){
  
  P <- length(ACSp)
  ACS_size=sapply(ACSp,ARS_to_size)
  ADS_size=sapply(ADSp,ARS_to_size)
  
  if(!is.null(cluster)){
    color = c()
    RB = rainbow(length(unique(cluster)))
    for (i in 1:length(cluster)) {
      if(cluster[i] != "scatter"){
        color[i] = RB[as.numeric(cluster[i])]
      }else{
        color[i] = "grey50"
      }
    }
  }else{
    color = "black"
  }

  if(!is.null(highlight.pathways)){
    index = 1:P
    index[-highlight.pathways] = ""
  }else{
    index = rep("",P)
  }
  data_ACS <- data.frame(ds1_score=DEevid1,ds2_score=DEevid2,
                         ACS_size=ACS_size,index=index,
                         color_pos=color)
  data_ADS <- data.frame(ds1_score=DEevid1,ds2_score=DEevid2,
                         ADS_size=ADS_size,index=index,
                         color_neg=color)
  
  p_pos <-ggplot(data_ACS, aes(x=ds2_score, y=ds1_score,label=index)) +
    geom_point(size = data_ACS$ACS_size*size.scale, shape=16, color = data_ACS$color_pos)+
    geom_text(size=8*size.scale,parse=TRUE,color="black",hjust = -0.05,vjust=-0.05) + 
    theme_bw() +
    coord_fixed(ylim=c(lb,ub),xlim=c(lb,ub)) + 
    labs(x="",y="") +
    scale_x_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1)) + 
    scale_y_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1)) + 
    theme(#legend.title = element_blank(),
      axis.line = element_line(colour = "black"),
      #axis.line=element_blank(),
      axis.text.x = element_text(size = 60*size.scale,face = "bold"),
      axis.text.y = element_text(size = 60*size.scale,face = "bold"),
      panel.border = element_blank(),
      panel.grid.major = element_line(linetype = 'solid',#size = 2,
                                      colour = "white"),
      panel.grid.minor = element_line(linetype = 'solid',
                                      colour = "white"),
      panel.background = element_rect(fill = "#FFF1E1")) +
    annotate("text", x = (ub-0.15), y = lb, fontface=2,
             label = paste(ds2,sep=""),
             size=40*size.scale,colour="blue",hjust=0.6,vjust=0.1) +
    annotate("text", x = lb, y = (ub-0.15), fontface=2,
             label=paste(ds1,sep=""),
             size=40*size.scale,colour="blue",vjust=0,hjust=0.2)
  
  p_neg <-ggplot(data_ADS, aes(x=ds1_score, y=ds2_score,label=index)) + #label=index
    geom_point(size = data_ADS$ADS_size*size.scale, shape=16, color = data_ADS$color_neg)+
    geom_text(size=8*size.scale,parse=TRUE,color="black",hjust = -0.05,vjust=-0.05) + 
    theme_bw() +
    coord_fixed(ylim=c(lb,ub),xlim=c(lb,ub)) + 
    labs(x="",y="") +
    scale_x_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1)) + 
    scale_y_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1)) + 
    theme(#legend.title = element_blank(),
      axis.line = element_line(colour = "black"),
      #axis.line=element_blank(),
      axis.text.x = element_text(size = 60*size.scale,face = "bold"),
      axis.text.y = element_text(size = 60*size.scale,face = "bold"),
      panel.border = element_blank(),
      panel.grid.major = element_line(linetype = 'solid',#size = 2,
                                      colour = "white"),
      panel.grid.minor = element_line(linetype = 'solid',
                                      colour = "white"),
      panel.background = element_rect(fill = "#EBF5FF")) +
    annotate("text", x = (ub-0.15), y = lb, fontface=2,
             label = paste(ds1,sep=""),
             size=40*size.scale,colour="blue",hjust=0.6,vjust=0.1) +
    annotate("text", x = lb, y = (ub-0.15), fontface=2,
             label=paste(ds2,sep=""),
             size=40*size.scale,colour="blue",vjust=0,hjust=0.2)
  
  #ggsave(filename=paste(ds2,"_",ds1,"_ACS_figure",".pdf",sep=""),p_pos,
  #width = 10, height = 10)
  
  #ggsave(filename=paste(ds1,"_",ds2,"_ADS_figure",".pdf",sep=""),p_neg,
  #width = 10, height = 10)  
  
  plist <- list(p_pos,p_neg)
  return(plist)
} 


DEevid_ACS_plot_clicked <- function(ds1,ds2,DEevid1,DEevid2,ACSp,ADSp,clickedVec,
                                    highlight.pathways = NULL,lb=0,ub=1,size.scale=4){
  
  P <- length(ACSp)
  ACS_size=sapply(ACSp,ARS_to_size)
  ADS_size=sapply(ADSp,ARS_to_size)
  
  if(!is.null(highlight.pathways)){
    index = 1:P
    index[-highlight.pathways] = ""
  }else{
    index = rep("",P)
  }
  # index=1:P
  # index[ACS_size==1] <- ""
  
  # data_pos <- data.frame(ds1_score=ds1.pos,ds2_score=ds2.pos,
  #                        ACS_size=ACS_size,index=index,clicked=clickedVec)
  # 
  # data_neg <- data.frame(ds1_score=ds1.neg,ds2_score=ds2.neg,
  #                        ACS_size=ACS_size,index=index,clicked=clickedVec)
  data_ACS <- data.frame(ds1_score=DEevid1,ds2_score=DEevid2,
                         ACS_size=ACS_size,index=index,
                         color_pos=clickedVec)
  data_ADS <- data.frame(ds1_score=DEevid1,ds2_score=DEevid2,
                         ADS_size=ADS_size,index=index,
                         color_neg=clickedVec)
  
  
  p_pos <-ggplot(data_ACS, aes(x=ds2_score, y=ds1_score, label=index)) + #label=index
    #geom_text(size=10,parse=TRUE,color="black",hjust = -0.05,vjust=-0.05) +
    geom_point(shape=16,size=(data_ACS$ACS_size*size.scale),color=data_ACS$color_pos) +
    theme_bw() +
    coord_fixed(ylim=c(lb,ub),xlim=c(lb,ub)) +
    labs(x="",y="") +
    scale_x_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1),
                       labels=c("0","0.5","1")) +
    scale_y_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1),
                       labels=c("0","0.5","1")) +
    #scale_color_manual(values=c("red", "grey"))+
    theme(#legend.title = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(size=60*size.scale, face="bold"),
      axis.text.y = element_text(size=60*size.scale, face="bold"),
      panel.border = element_blank(),
      # plot.margin = unit(c(0,0,0,0),"cm"))
      panel.grid.major = element_line(linetype = 'solid',#size = 2,
                                      colour = "white"),
      panel.grid.minor = element_line(linetype = 'solid',
                                      colour = "white"),
      panel.background = element_rect(fill = "#FFF1E1"))
  
  p_pos <- p_pos + annotate("text", x = (ub-0.15), y = lb, fontface=2,
                            label = paste(ds2,sep=""),
                            size=40*size.scale,colour="blue",hjust=0.6,vjust=0.1)
  
  p_pos <- p_pos + annotate("text", x = lb, y = (ub-0.15), fontface=2,
                            label=paste(ds1,sep=""),
                            size=40*size.scale,colour="blue",vjust=0,hjust=0.2)
  
  
  p_neg <-ggplot(data_ADS, aes(x=ds1_score, y=ds2_score, label=index)) + #label=index
    #geom_text(size=10,parse=TRUE,color="black",hjust = -0.05,vjust=-0.05) +
    geom_point(shape=16,size=data_ADS$ADS_size*size.scale,color=data_ADS$color_neg) +
    geom_text(size=8*size.scale,parse=TRUE,color="black",hjust = -0.05,vjust=-0.05) +
    theme_bw() +
    coord_fixed(ylim=c(lb,ub),xlim=c(lb,ub)) +
    labs(x="",y="") +
    scale_x_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1)) +
    scale_y_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1)) +
    #scale_color_manual(values=c("red", "grey"))+
    theme(#legend.title = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(size=60*size.scale, face = "bold"),
      axis.text.y = element_text(size=60*size.scale, face = "bold"),
      #plot.margin = unit(c(0,0,0,0),"cm")
      panel.border = element_blank(),
      panel.grid.major = element_line(linetype = 'solid',#size = 2,
                                      colour = "white"),
      panel.grid.minor = element_line(linetype = 'solid',
                                      colour = "white"),
      panel.background = element_rect(fill = "#EBF5FF"))
  
  p_neg <- p_neg + annotate("text", x = (ub-0.15), y = lb, fontface=2,
                            label = paste(ds1,sep=""),
                            size=40*size.scale,colour="blue",hjust=0.6,vjust=0.1)
  p_neg <- p_neg + annotate("text", x = lb, y = (ub-0.15), fontface=2,
                            label=paste(ds2,sep=""),
                            size=40*size.scale,colour="blue",vjust=0,hjust=0.2)
  
  #ggsave(filename=paste("ACS_DE",ds1,"_",ds2,"_pos",".pdf",sep=""),p_pos,
  #width = 10, height = 10)
  
  #ggsave(filename=paste("ACS_DE",ds2,"_",ds1,"_neg",".pdf",sep=""),p_neg,
  #width = 10, height = 10)
  
  plist <- list(p_pos,p_neg)
  return(plist)
}

