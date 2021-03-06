#' Extrapolate objects from ClusterAnalysis
#' 
#' @description Extrapolate objects from ClusterAnalysis output list (see \code{\link{ClusterAnalysis}}).
#'
#' @param stability.list The list obtained from the ClusterAnalysis function. (see \code{\link{ClusterAnalysis}})
#' @param data CONNECTORList. (see \code{\link{DataImport}} or \code{\link{DataTruncation}})
#' @param G  The number of clusters.
#' 
#' @details 
#' \itemize{
#' \item{IndexesPlot.Extrapolation}{extrapolates from ClusterAnalysis output list the box plot fixing the h value.}
#' \item{ConsMatrix.Extrapolation}{extrapolates from ClusterAnalysis output list the Consensus Matrix for G fixed.}
#' \item{MostProbableClustering.Extrapolation}{extrapolates from ClusterAnalysis output list the most frequent clustering among the several runs obtained  for G fixed.}
#' }
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @name ExtrapolationNewFuncs
NULL
#> NULL
#' @import ggplot2 viridis reshape2 statmod
#' 
#' @rdname ExtrapolationNewFuncs
#' @export
IndexesPlot.ExtrapolationNew <- function(stability.list){
  Clusters.List<-stability.list$Clusters.List 
  G <- as.numeric(sub("G","",names(Clusters.List)))
  
  IndexesValues <- IndexeValues.generation(Clusters.List, G)
  Indexes.MostProb <- IndexesValues$Indexes.MostProb
  Indexes.Rep <- IndexesValues$Indexes.Rep
  
  Box.pl<- ggplot(data= Indexes.Rep)+
    facet_wrap(~Index,scales = "free")+
    geom_violin(aes(x=Cluster,y=V,fill=ClusterH,group=Cluster))+
    geom_line(data= Indexes.MostProb,aes(x=Cluster,y=V,col="Most probable"))+
    geom_jitter(aes(x=Cluster,y=V),color="black", width = .1, height = 0, alpha=0.5)+   
    scale_fill_viridis("",discrete = TRUE, alpha=0.6) +
    #geom_point(col="red",aes(x=Cluster,y=V,size=Freq/runs) )+
    labs(size="Counts freq.",col= "",x = "Number of Clusters",y="Indexes Values")+
    scale_x_continuous(breaks = G)+
    scale_color_manual(values = c("Most probable" = "blue") ) +
    theme(axis.text=element_text(size = 14, hjust = 0.5),
          axis.text.x=element_text(vjust=0.5, hjust=1),
          axis.title=element_text(size=16,face="bold"),
          axis.line = element_line(colour="black"),
          plot.title=element_text(size=30, face="bold", vjust=1, lineheight=0.6),
          legend.text=element_text(size=16),
          legend.position="bottom",
          legend.key=element_blank(),
          legend.title = element_text(size=16,face="bold"),
          legend.key.size = unit(.9, "cm"),
          legend.key.width = unit(.9,"cm"),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          plot.margin=unit(c(5,5,5,5),"mm"),
          strip.text = element_text(size = 20)) 
  
  return(list(Plot=Box.pl, IndexesValues = IndexesValues$DataIndexes) )
}

#' @rdname ExtrapolationNewFuncs
#' @export
ConsMatrix.ExtrapolationNew <- function(stability.list,data){
  Clusters.List<-stability.list$Clusters.List 
  runs <-stability.list$runs 

  G <- as.numeric(sub("G","",names(Clusters.List)))
  
  IndexesValues <- IndexeValues.generation(Clusters.List, G)
  Indexes.MostProb <- IndexesValues$Indexes.MostProb
  
  Freq.ConfigCl<-unique(Indexes.MostProb[,c("Cluster","Config")])
  
  ConsensusInfo<-lapply(1:length(G), function(Gind){
    ConsM.generation(Gind,Clusters.List,runs,data,Freq.ConfigCl)
    })
  names(ConsensusInfo)<-paste0("G",G)
  
  return(ConsensusInfo)
}

#' @rdname ExtrapolationNewFuncs
#' @export
MostProbableClustering.ExtrapolationNew <- function(stability.list, G){
  Clusters.List<-stability.list$Clusters.List 
  runs <-stability.list$runs 
  
  G.all <- as.numeric(sub("G","",names(Clusters.List)))
  
  IndexesValues <- IndexeValues.generation(Clusters.List, G.all)
  Indexes.MostProb <- IndexesValues$Indexes.MostProb
  
  Freq.ConfigCl<-unique(Indexes.MostProb[,c("Cluster","Config")])
  ############# the most probably clustering:
  IndexBestClustering <- Freq.ConfigCl[which(G.all == G),"Config"]
  MostProbableClustering<-Clusters.List[[which(G.all == G)]]$ClusterAll[[IndexBestClustering]]
  
  return(MostProbableClustering)
}

ConsM.generation<-function(Gind,ALL.runs,runs,data,Freq.ConfigCl) 
{
  FixedG.runs<-ALL.runs[[Gind]]$ClusterAll
  L1<- length(FixedG.runs)
  #Check some errors in the predictions:
  FixedG.runs.tmp <-lapply(1:L1,function(x){
    if(!is.character(FixedG.runs[[x]]$Error))
      FixedG.runs[[x]]
    else NA
  } )
  # if(length(which(is.na(FixedG.runs.tmp)))!=0){
  #   ErrorConfiguration <- list(FromFitting=ALL.runs[[Gind]]$ErrorConfigurationFit,
  #                              FromPrediction = ALL.runs[[which(is.na(FixedG.runs.tmp))]] )
  # }else{
  #   ErrorConfiguration <-list(FromFitting=ALL.runs[[Gind]]$ErrorConfigurationFit,
  #                             FromPrediction = NULL)
  # }
  ###
  
  L1<- length(FixedG.runs)
  ############# first we found the most probably clustering:
  IndexBestClustering <- Freq.ConfigCl[Gind,"Config"]
  BestClustering<-FixedG.runs[[IndexBestClustering]]
  
  ##########################################################
  #### Let build the consensus matrix
  
  consensus.runs<-lapply(1:L1, function(x){
    consensusM<-diag(0,length(data$LenCurv))
    fcm.G<-FixedG.runs[[x]]
    cl.vector<-fcm.G$FCM$cluster$cluster.member
    for(i in 1:length(data$LenCurv)) consensusM[i,which(cl.vector%in%cl.vector[i])]<-FixedG.runs[[x]]$ParamConfig.Freq
    consensusM
  })
  
  consensusM<-Reduce("+",consensus.runs)
  
  colnames(consensusM)<- data$LabCurv$ID
  row.names(consensusM)<- data$LabCurv$ID
  
  consensusM<-as.data.frame(consensusM)/runs
  
  ################
  ### ordering the samples to have all curves in the same cluster close with each otehrs.
  
  cl.memer<-BestClustering$FCM$cluster$cluster.m
  #  mat[do.call(order, as.data.frame(mat)),]
  
  fcm<-BestClustering
  curvename<-data$LabCurv$ID
  curve<-fcm$FCM$prediction$gpred
  rownames(curve)<-curvename
  
  grid<-fcm$FCM$fit$grid
  
  gauss.quad(10) -> gauss
  a <- min(grid)
  b <- max(grid)
  itempi <- (a+b)/2 + (b-a)/2*gauss$nodes
  
  match(itempi,grid) -> itimeindex 
  
  fxG <- (curve[,itimeindex] )^2
  int <- (b-a)/2 * rowSums( gauss$weights * fxG )
  dist.curve <- sqrt(int)
  names(which.min(dist.curve))->lowest.curve
  
  m.lowercurve<-matrix(curve[lowest.curve,itimeindex],nrow = length(curve[,1]),ncol = length(itimeindex),byrow = T)
  fxG <- (curve[,itimeindex]-m.lowercurve )^2
  int <- (b-a)/2 * rowSums( gauss$weights * fxG )
  dist.curve <- sqrt(int)
  ################## Sorting the names!!!
  
  #1) sorting depending by the l2 distance with the zero x-axis!
  curvename.ordered<-names(sort(dist.curve)) 
  
  #2) Sorting depending on the cluster membership
  names(cl.memer)<-curvename
  cl.memer<-sort(cl.memer)
  #3) Defining a dataframe in order to sort the curves depending by the cluster and distance!
  
  namefroml2<-1:length(curvename)
  names(namefroml2)<-curvename.ordered
  namefroml2<-namefroml2[names(cl.memer)] # ordering the curve ordered by the L2 distance
  
  curve.name.dtfr<-data.frame(namefroml2=namefroml2,cl.mem=cl.memer,1:length(curvename.ordered))
  
  # ind<-curve.name.dtfr[order(curve.name.dtfr$namefroml2,curve.name.dtfr$cl.mem,curve.name.dtfr$namefroml2),3]
  ind<-curve.name.dtfr[order(curve.name.dtfr$cl.mem,curve.name.dtfr$namefroml2),3]
  curvename.ordered<-names(cl.memer)[ind]
  consensusM <- consensusM[curvename.ordered,curvename.ordered]
  consensusM <- as.matrix(consensusM)
  consensusM[upper.tri(consensusM)] <- NA
  
  m<-melt( as.matrix(consensusM) , na.rm = TRUE)
  
  m$Var2<-factor(m$Var2,levels = curvename.ordered  )
  m$Var1<-factor(m$Var1,levels = rev(curvename.ordered )  )
  
  ############
  length.cl<-rev(as.numeric(table(cl.memer[ind])[unique(cl.memer[ind])]))
  coor<-c(0,cumsum(length.cl))
  coor.2<-sum(length.cl)-coor
  
  cluster.lines<-data.frame(x=c(coor[-length(coor)],coor[-length(coor)]),
                            y=c(coor.2[-1],coor.2[-length(coor.2)]),
                            xend=c(coor[-1],coor[-length(coor)]),
                            yend=c(coor.2[-1],coor.2[-1] ) ) +.5
  
  x.text<-(cluster.lines$xend-cluster.lines$x)/2 + cluster.lines$x
  y.text<-(cluster.lines$yend-cluster.lines$y)/2 + cluster.lines$y
  x.text <- x.text[1:length(length.cl)]
  y.text <- y.text[-(1:length(length.cl) )]
  
  lab = rev(BestClustering$FCM$cluster$cluster.names)[rev(BestClustering$FCM$cluster$cluster.names) %in% unique(names(BestClustering$FCM$cluster$cluster.member)) ]
  
  ### Freq in each cluster
  G <- unique(cl.memer)
  Freq.cl<-lapply(G,function(g){
    as.numeric(names(cl.memer[which(cl.memer == g)])) -> curve.cl.fixed
    indexes <- which(m$Var1 %in% curve.cl.fixed & m$Var2 %in% curve.cl.fixed & as.numeric(levels(m$Var1)[c(m$Var1)]) - as.numeric(levels(m$Var2)[c(m$Var2)])!=0)
    data.frame(Median= median(m[indexes,3]), Mean = round(mean(m[indexes,3]),2),Cluster = g)
  })
  Freq.cl<-do.call("rbind",Freq.cl)
  Freq.cl$Cluster<-BestClustering$FCM$cluster$cluster.names[Freq.cl$Cluster]
  Freq.cl[order(match(Freq.cl$Cluster, lab)),]
  MeanFreq<-mean(Freq.cl$Mean)
  Freq.cl$Mean[is.na(Freq.cl$Mean)] <- 1
  lab.new<-rev(sapply(1:length(G),function(i) paste(Freq.cl[i,c("Cluster","Mean")],collapse = ": ") ))
  ##
  
  ConsensusPlot <- ggplot() +
    geom_tile(data = m, 
              aes(x = Var1, y = Var2, fill = value)) +
    labs(x = "",y = "", fill = "Same cluster \ncounting") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_gradient2(midpoint = 0.5, low = "blue", 
                         mid = "yellow", high = "red") +
    geom_segment(data = cluster.lines, 
                 aes(x, y, xend = xend, yend = yend), size = 1.5, 
                 inherit.aes = F) +
    labs(title = "Consensus Matrix",
         subtitle = paste("Black line for the most probable clustering: ",MeanFreq) ) + 
    annotate(geom = "text", x = x.text, y = y.text, 
             label = lab.new , 
             size = 6)
  
  return(list(ConsensusMatrix = consensusM,
              ConsensusPlot = ConsensusPlot) )
} 

IndexeValues.generation <- function(Clusters.List, G){
  l.tight <- lapply(1:length(G),function(i){
    ClusterAll <- Clusters.List[[i]]$ClusterAll
    V.all<- lapply(1:(length(ClusterAll)),function(j){
      if(!is.character(ClusterAll[[j]]$Error) )
        data.frame(Config = j, V = ClusterAll[[j]]$Cl.Info$Tight.indexes, Freq= ClusterAll[[j]]$ParamConfig.Freq)
      else NA
    })
    V.all <- do.call("rbind",V.all)
    V.all$Cluster <-G[i]
    V.all$Index = "Tightness"
    V.all$ClusterH = paste("G:",G[i],"; h:",Clusters.List[[i]]$h.selected)
    return(V.all)
  })
  l.fdb <- lapply(1:length(G),function(i){
    ClusterAll <- Clusters.List[[i]]$ClusterAll
    V.all<- lapply(1:(length(ClusterAll)),function(j){
      if(!is.character(ClusterAll[[j]]$Error) )
        data.frame(Config = j, V = ClusterAll[[j]]$Cl.Info$Coefficents$fDB.index, Freq= ClusterAll[[j]]$ParamConfig.Freq)
      else NA
    })
    V.all <- do.call("rbind",V.all)
    V.all$Cluster <-G[i]
    V.all$Index = "fDB"
    V.all$ClusterH = paste("G:",G[i],"; h:",Clusters.List[[i]]$h.selected)
    return(V.all)
  })
  l.fdb1 <- lapply(1:length(G),function(i){
    ClusterAll <- Clusters.List[[i]]$ClusterAll
    V.all<- lapply(1:(length(ClusterAll)),function(j){
      if(!is.character(ClusterAll[[j]]$Error) )
        data.frame(Config = j, V = ClusterAll[[j]]$Cl.Info$Deriv.Coefficents$fDB.index,
                   Freq= ClusterAll[[j]]$ParamConfig.Freq)
      else NA
    })
    V.all <- do.call("rbind",V.all)
    V.all$Cluster <-G[i]
    V.all$Index = "Deriv.fDB"
    V.all$ClusterH = paste("G:",G[i],"; h:",Clusters.List[[i]]$h.selected)
    return(V.all)
  })
  l.fdb2 <- lapply(1:length(G),function(i){
    ClusterAll <- Clusters.List[[i]]$ClusterAll
    V.all<- lapply(1:(length(ClusterAll)),function(j){
      if(!is.character(ClusterAll[[j]]$Error) )
        data.frame(Config = j, V = ClusterAll[[j]]$Cl.Info$Deriv2.Coefficents$fDB.index, Freq= ClusterAll[[j]]$ParamConfig.Freq)
      else NA
    })
    V.all <- do.call("rbind",V.all)
    V.all$Cluster <-G[i]
    V.all$Index = "Deriv2.fDB"
    V.all$ClusterH = paste("G:",G[i],"; h:",Clusters.List[[i]]$h.selected)
    return(V.all)
  })
  ## Grouping all the indexes:
  Data <- list(Tight =do.call("rbind", l.tight),
               fDB =do.call("rbind", l.fdb),
               fDB1 =do.call("rbind", l.fdb1),
               fDB2 =do.call("rbind", l.fdb2))
  ##
  dt.fr <- rbind(do.call("rbind", l.tight), do.call("rbind",l.fdb))
  dt.fr.rep <- lapply(1:length(dt.fr[, 1]), function(i) {
    freq <- dt.fr[i, "Freq"]
    do.call("rbind", lapply(1:freq, function(j) dt.fr[i,]))
  })
  dt.fr.rep <- do.call("rbind", dt.fr.rep)
  dt.fr.max <- aggregate(dt.fr$Freq, by = list(Cluster = dt.fr$Cluster,
                                               Index = dt.fr$Index,
                                               ClusterH = dt.fr$ClusterH),
                         FUN = "max")
  colnames(dt.fr.max)[4] <- "Freq"
  dt.fr.max <- merge(dt.fr.max, dt.fr)
  dt.fr.max <- aggregate(dt.fr.max$V, by = list(Cluster = dt.fr.max$Cluster,
                                                Index = dt.fr.max$Index, ClusterH = dt.fr.max$ClusterH,
                                                Freq = dt.fr.max$Freq), FUN = "min")
  colnames(dt.fr.max)[5] <- "V"
  dt.fr.max <- merge(dt.fr.max, dt.fr)
  
  return(list(Indexes.MostProb = dt.fr.max,Indexes.Rep = dt.fr.rep,DataIndexes = Data))
}
