#' Extrapolate objects from ClusterAnalysis
#' 
#' @description Extrapolate objects from ClusterAnalysis output list (see \code{\link{ClusterAnalysis}}).
#'
#' @param stability.list The list obtained from the ClusterAnalysis function. (see \code{\link{ClusterAnalysis}})
#' @param G  The number of clusters.
#' @param q The quantiles used to calculate the time grid interval on which the distances are calculated. If NULL then the time grid outliers will be ignored through the distance calculation. If double (0<q<1) then the cutting is symmetrical w.r.t. the quantile setled (e.g., q = 0.25). If a vector, then the minimum value is used for the lower cutting and the maximum value for the upper cutting.
#' @details 
#' \itemize{
#' \item{IndexesPlot.Extrapolation}{extrapolates from ClusterAnalysis output list the box plot fixing the h value.}
#' \item{ConsMatrix.Extrapolation}{extrapolates from ClusterAnalysis output list the Consensus Matrix for G fixed.}
#' \item{MostProbableClustering.Extrapolation}{extrapolates from ClusterAnalysis output list the most frequent clustering among the several runs obtained  for G fixed.}
#' }
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @name ExtrapolationFuncs
NULL
#> NULL
#' @import ggplot2 viridis tidyr statmod dplyr
#' 
#' @rdname ExtrapolationFuncs
#' @export
IndexesPlot.Extrapolation <- function(stability.list,q=NULL){
  
  Clusters.List<-stability.list$Clusters.List 
  G <- as.numeric(sub("G","",names(Clusters.List)))
  
  Allerrors = which(sapply(1:length(Clusters.List),
                           function(x) is.character(Clusters.List[[x]]) ))
  if(length(Allerrors)>0)
  {
    errors = sapply(Allerrors, function(x) 
      paste("Cluster",G[x],"has the following errors: ", Clusters.List[[x]]) )
    
    print(errors)
    G = G[-Allerrors]
    Clusters.List = Clusters.List[-Allerrors]
  }
  if(length(G)==0)
  {
    return("All the CONNECTOR runs have errors!")
  }
  
  IndexesValues.list <- IndexesValues.calculation(Clusters.List, G,q)
  IndexesValues <- IndexesValues.extrap(IndexesValues.list,Clusters.List, G)
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

#' @rdname ExtrapolationFuncs
#' @export
ConsMatrix.Extrapolation <- function(stability.list,q=NULL){
  Clusters.List<-stability.list$Clusters.List 
  data <-stability.list$CONNECTORList
  runs <-stability.list$runs 
  
  G <- as.numeric(sub("G","",names(Clusters.List)))
  Allerrors = which(sapply(1:length(Clusters.List),
                           function(x) is.character(Clusters.List[[x]]) ))
  if(length(Allerrors)>0)
  {
    errors = sapply(Allerrors, function(x) 
      paste("Cluster",G[x],"has the following errors: ", Clusters.List[[x]]) )
    
    print(errors)
    G = G[-Allerrors]
    Clusters.List = Clusters.List[-Allerrors]
  }
  if(length(G)==0)
  {
    return("All the CONNECTOR runs have errors!")
  }
  
  IndexesValues.list <- IndexesValues.calculation(Clusters.List, G,q)
  IndexesValues <- IndexesValues.extrap(IndexesValues.list,Clusters.List, G)
  Indexes.MostProb <- IndexesValues$Indexes.MostProb
  
  Freq = Indexes.MostProb %>% 
    dplyr::select(Cluster,Config) %>%
    dplyr::distinct(Cluster, .keep_all = T)
  
  Freq.ConfigCl<- Freq[order(Freq$Cluster),]
  
  ConsensusInfo<-lapply(1:length(G), function(Gind){
    ConsM.generation(Gind,Clusters.List,runs,data,Freq.ConfigCl,q)
  })
  names(ConsensusInfo)<-paste0("G",G)
  
  return(ConsensusInfo)
}

#' @rdname ExtrapolationFuncs
#' @export
MostProbableClustering.Extrapolation <- function(stability.list, G,q=NULL){
  Clusters.List<-stability.list$Clusters.List 
  runs <-stability.list$runs 
  
  G.all <- as.numeric(sub("G","",names(Clusters.List)))
  Allerrors = which(sapply(1:length(Clusters.List),
                           function(x) is.character(Clusters.List[[x]]) ))
  if(length(Allerrors)>0)
  {
    errors = sapply(Allerrors, function(x) 
      paste("Cluster",G.all[x],"has the following errors: ", Clusters.List[[x]]) )
    
    print(errors)
    G.all = G.all[-Allerrors]
    Clusters.List = Clusters.List[-Allerrors]
  }
  if(length(G.all)==0)
  {
    return("All the CONNECTOR runs have errors!")
  }
  
  IndexesValues.list <- IndexesValues.calculation(Clusters.List, G.all,q)
  IndexesValues <- IndexesValues.extrap(IndexesValues.list,Clusters.List, G.all)
  Indexes.MostProb <- IndexesValues$Indexes.MostProb
  
  Freq.ConfigCl<-unique(Indexes.MostProb[,c("Cluster","Config")])
  ############# the most probably clustering:
  IndexBestClustering <- Freq.ConfigCl[which(Freq.ConfigCl$Cluster == G),"Config"]
  
  # If there are more configuration with the same freq. than we will selct the one with smaller fdb
  if(length(IndexBestClustering) >1){
    fdb = sapply(IndexBestClustering,function(ic){
      return(c(ic,IndexesValues.list[[paste0("G", G)]][[ic]]$Coefficents$fDB.index))
    })
    IndexBestClustering <- fdb[1,which.min(fdb[2,])]
  }
  ##
  MostProbableClustering<-Clusters.List[[paste0("G",G)]]$ClusterAll[[IndexBestClustering]]
  
  MostProbableClustering$Cl.Info <- IndexesValues.list[[paste0("G",G)]][[IndexBestClustering]]
  MostProbableClustering$CONNECTORList <- stability.list$CONNECTORList
  
  return(MostProbableClustering)
}

ConsM.generation<-function(Gind,ALL.runs,runs,data,Freq.ConfigCl,q) 
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
  
  # gauss.quad(10) -> gauss
  # a <- min(grid)
  # b <- max(grid)
  # itempi <- (a+b)/2 + (b-a)/2*gauss$nodes
  # 
  # match(itempi,grid) -> itimeindex 
  # 
  # itimeindex <- match(database$Time[database$ID == i],grid)
  # fxG <- (curve[,itimeindex] )^2
  # int <- (b-a)/2 * rowSums( gauss$weights * fxG )
  # dist.curve <- sqrt(int)
  
  dist.curve<-L2dist.curve20(clust = fcm$FCM$cluster$cluster.member,
                             fcm.curve = fcm$FCM$prediction,
                             database = data$Dataset,
                             q = q)
  
  names(dist.curve) <- curvename
  # names(which.min(dist.curve))->lowest.curve
  # 
  # m.lowercurve<-matrix(curve[lowest.curve,itimeindex],
  #                      nrow = length(curve[,1]),
  #                      ncol = length(itimeindex),byrow = T)
  # fxG <- (curve[,itimeindex]-m.lowercurve )^2
  # int <- (b-a)/2 * rowSums( gauss$weights * fxG )
  # dist.curve <- sqrt(int)
  ################## Sorting the names!!!
  
  #1) sorting depending by the l2 distance with the zero x-axis!
  curvename.ordered<-names(sort(dist.curve)) 
  
  #2) Sorting depending on the cluster membership
  
  cl.memer.tmp <- data.frame(curvename,names(cl.memer),cl.memer)
  
  cl.memer.tmp<-cl.memer.tmp[order(cl.memer.tmp$names.cl.memer.),]
  cl.memer <- cl.memer.tmp$cl.memer
  names(cl.memer) <- cl.memer.tmp$curvename
  
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
  
  m<-gather(as.data.frame(consensusM) , na.rm = TRUE)
  
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
  Freq.cl<-Freq.cl[order(match(Freq.cl$Cluster, lab)),]
  MeanFreq<-mean(Freq.cl$Mean,na.rm = T)
  Freq.cl$Mean[is.na(Freq.cl$Mean)] <- 1
  
  labText<-sapply(1:length(G),function(i) paste(Freq.cl[i,c("Cluster","Mean")],collapse = ": ") )
  
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
    labs(title = paste( "Consensus Matrix, G=", length(G) ),
         subtitle = paste("Black line for the most probable clustering: ",MeanFreq) ) + 
    annotate(geom = "text", x = x.text, y = y.text, 
             label = labText , 
             size = 6)
  
  return(list(ConsensusMatrix = consensusM,
              ConsensusPlot = ConsensusPlot) )
} 

IndexesValues.calculation <- function(Clusters.List, G,q){
  
  Indexes <- lapply(1:length(G),function(i){
    ## Check the same parameter configurations:
    Indexes.Uniq.Par<-length(Clusters.List[[i]]$ClusterAll)
    
    Cl.Info<- lapply(1:Indexes.Uniq.Par,function(j){
      if(is.character(Clusters.List[[i]]$ClusterAll[[j]]$Error) ){
        return(Clusters.List[[i]]$ClusterAll[[j]]$Error)
      }else{
        Clusters.List[[i]]$ClusterAll[[j]]$FCM$prediction -> fcm.prediction
        Clusters.List[[i]]$ClusterAll[[j]]$FCM$cluster$cl.info$class.pred -> cluster
        Clusters.List[[i]]$ClusterAll[[j]]$FCM$cluster$ClustCurve[,c("ID","Observation","Time")] -> database
        Clusters.List[[i]]$ClusterAll[[j]]$FCM$fit -> fcm.fit  
        ## Goodness coefficents calculation
        fcm.prediction$meancurves->meancurves
        distances <- L2dist.curve2mu(clust=cluster, fcm.curve = fcm.prediction, database = database, fcm.fit = fcm.fit, deriv = 0,q)
        distances.zero<-L2dist.mu20(clust=cluster,fcm.prediction,database = database,fcm.fit,deriv=0,q)
        Coefficents<-Rclust.coeff(clust=cluster, fcm.curve = fcm.prediction, database = database, fcm.fit = fcm.fit, deriv = 0,q)
        Deriv.Coefficents<-Rclust.coeff(clust=cluster, fcm.curve = fcm.prediction, database = database, fcm.fit = fcm.fit, deriv = 1,q)
        Deriv2.Coefficents<-Rclust.coeff(clust=cluster, fcm.curve = fcm.prediction, database = database, fcm.fit = fcm.fit, deriv = 2,q)
        
        return(list(Tight.indexes = sum(distances),
                    Coefficents=Coefficents,
                    Deriv.Coefficents=Deriv.Coefficents,
                    Deriv2.Coefficents=Deriv2.Coefficents))
      }
      
    })
    return(Cl.Info)
  })
  names(Indexes) <-paste0("G",G)
  
  return(Indexes)
}

IndexesValues.extrap <- function(IndexesValues.list,Clusters.List,G){
  
  l.tight <- lapply(1:length(G),function(i){
    ClusterAll <- IndexesValues.list[[i]]
    V.all<- lapply(1:(length(ClusterAll)),function(j){
      if(!is.character(ClusterAll[[j]]) )
        data.frame(Config = j, 
                   V = ClusterAll[[j]]$Tight.indexes, 
                   Freq= Clusters.List[[i]]$ClusterAll[[j]]$ParamConfig.Freq)
      # else NA
    })
    V.all <- do.call("rbind",V.all)
    V.all$Cluster <-G[i]
    V.all$Index = "Tightness"
    V.all$ClusterH = paste("G:",G[i],"; h:",Clusters.List[[i]]$h.selected)
    return(V.all)
  })
  l.fdb <- lapply(1:length(G),function(i){
    ClusterAll <- IndexesValues.list[[i]]
    V.all<- lapply(1:(length(ClusterAll)),function(j){
      if(!is.character(ClusterAll[[j]]) )
        data.frame(Config = j, 
                   V = ClusterAll[[j]]$Coefficents$fDB.index,
                   Freq= Clusters.List[[i]]$ClusterAll[[j]]$ParamConfig.Freq)
      # else NA
    })
    V.all <- do.call("rbind",V.all)
    V.all$Cluster <-G[i]
    V.all$Index = "fDB"
    V.all$ClusterH = paste("G:",G[i],"; h:",Clusters.List[[i]]$h.selected)
    return(V.all)
  })
  l.fdb1 <- lapply(1:length(G),function(i){
    ClusterAll <- IndexesValues.list[[i]]
    V.all<- lapply(1:(length(ClusterAll)),function(j){
      if(!is.character(ClusterAll[[j]]) )
        data.frame(Config = j, V = ClusterAll[[j]]$Deriv.Coefficents$fDB.index,
                   Freq= Clusters.List[[i]]$ClusterAll[[j]]$ParamConfig.Freq)
      #  else NA
    })
    V.all <- do.call("rbind",V.all)
    V.all$Cluster <-G[i]
    V.all$Index = "Deriv.fDB"
    V.all$ClusterH = paste("G:",G[i],"; h:",Clusters.List[[i]]$h.selected)
    return(V.all)
  })
  l.fdb2 <- lapply(1:length(G),function(i){
    ClusterAll <-IndexesValues.list[[i]]
    V.all<- lapply(1:(length(ClusterAll)),function(j){
      if(!is.character(ClusterAll[[j]]) )
        data.frame(Config = j, V = ClusterAll[[j]]$Deriv2.Coefficents$fDB.index, 
                   Freq= Clusters.List[[i]]$ClusterAll[[j]]$ParamConfig.Freq)
      #  else NA
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
  colnames(dt.fr.max)[5] <- "Obs"
  dt.fr.max <- merge(dt.fr.max, dt.fr)
  
  return(list(Indexes.MostProb = dt.fr.max,Indexes.Rep = dt.fr.rep,DataIndexes = Data))
}
