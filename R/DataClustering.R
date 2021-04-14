#' Cluster Stability
#' 
#'@description
#'
#'  Fits and clusters the data with respect to the Functional Clustering Model [Sugar and James]. The plots based on the Elbow Method and on the functional Davies Bouldin (fDB) indexes are generated to properly guide the choice of the number of clusters. As explained in [Sugar and James], to have a simple low-dimensional representation of the individual curves and to reduce the number of parameters to be estimated, h value must be equals or lower than \eqn{min(p,G-1)}.
#'  
#' @param data CONNECTORList. (see \code{\link{DataImport}} or \code{\link{DataTruncation}})
#' @param G  The vector/number of possible clusters.
#' @param p The dimension of the natural cubic spline basis. (see \code{\link{BasisDimension.Choice}})
#' @param h
#' @param runs Number of runs.
#' @param seed Seed for the kmeans function.
#' @param save If TRUE then the growth curves plot truncated at the "TruncTime" is saved into a pdf file.
#' @param path The folder path where the plot(s) will be saved. If it is missing, the plot is saved in the current working  directory.
#' @param Cores Number of cores to parallelize computations.
#' @param PercPCA=.85
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @return StabilityAnalysis returns a list of (i) lists, called ConsensusInfo, reporting for each G and h: the Consensus Matrix, either as a NxN matrix, where N is the number of samples, or plot, and the most probable clustering obtained from running several times the method; (ii) the box plots showing both the Elbow plot considering the total tightness and the box plots of the fDB indexes for each G; and finally, (iii) the seed. See \code{\link{BoxPlot.Extrapolation}} and \code{\link{MostProbableClustering.Extrapolation}}.
#' 
#' @details 
#'  Connector provides two different plots to properly guide the choice of the number of clusters:
#'  \itemize{
#'  \item{Elbow Plot:}{ a plot in which the total tightness against the number of clusters is plotted. A proper number of clusters can be inferred as large enough to let the total tightness drop down to relatively little values but as the smallest over which the total tightness does not decrease substantially (we look for the location of an "elbow" in the plot). }
#'  \item{Box Plot:}{ a plot in which for each number of clusters , G, the functional Davies-Buldin (fDB), the cluster separation measure index, is plotted as a boxplot. A proper number of clusters can be associated to the minimum fDB value. }
#'  }
#'  
#'    The proximities measures choosen is defined as follow
#'  \deqn{D_q(f,g) = \sqrt( \integral | f^{(q)}(s)-g^{(q)}(s) |^2 ds ), d=0,1,2}
#' where f and g are two curves and f^{(q)} and g^{(q)} are their q-th derivatives. Note that for q=0, the equation becomes the distance induced by the classical L^2-norm.
#' Hence, we can define the following indexes for obtaining a cluster separation measure:
#' \itemize{
#' \item{T:}{ the total tightness representing the dispersion measure given by
#' \deqn{ T = \sum_{k=1}^G \sum_{i=1}^n D_0(\hat{g}_i, \bar{g}^k)}; }
#'  \item{S_k:}{ 
#'  \deqn{ S_k = \sqrt{\frac{1}{G_k} \sum_{i=1}^{G_k} D_q^2(\hat{g}_i, \bar{g}^k);} }
#'  with G_k the number of curves in the k-th cluster;
#'  }
#'  \item{M_{hk}:}{ the distance between centroids (mean-curves) of h-th and k-th cluster 
#'  \deqn{M_{hk} =  D_q(\bar{g}^h, \bar{g}^k);}
#'  }
#'  \item{R_{hk}:}{  a measure of how good the clustering is,
#'  \deqn{R_{hk} =   \frac{S_h + S_k}{M_{hk}};}
#'  }
#'  \item{fDB_q:}{ functional Davies-Bouldin index, the cluster separation measure
#'  \deqn{fDB_q = \frac{1}{G} \sum_{k=1}^G \max_{h \neq k} {  \frac{S_h + S_k}{M_{hk}} } }
#'  }
#' }
#' 
#'
#' @seealso MostProbableClustering.Extrapolation, BoxPlot.Extrapolation, ConsMatrix.Extrapolation.
#' 
#' @import ggplot2 reshape2 RColorBrewer statmod parallel  Matrix splines statmod
#' @importFrom cowplot plot_grid
#' @export

ClusterAnalysis<-function(data,G,p,h=NULL,runs=50,seed=123,save=FALSE,path=NULL,Cores=1,PercPCA=.85)
{
  # library(statmod); library(parallel); library(Matrix); library(splines); library(statmod);library( reshape2)
  # source('~/GIT/R_packages_project/connector/R/FCM_fclust.R', echo=F)
  # source('~/GIT/R_packages_project/connector/R/ClusterCoefficents.R', echo=F)
  # Cores=2; seed=2404; tol = 0.001; maxit = 20; PercPCA=.85; runs=50; p=5
  #  data = CONNECTORList
  #  G=4

  Clusters.List<-list()
  for( g in 1:length(G) ){
    Clusters.List[[g]] = FCM.estimation(data = data,
                                        G = G[g],
                                        p = p,
                                        h=h,
                                        seed=seed,
                                        PercPCA= PercPCA,
                                        Cores=Cores,
                                        runs = runs
                                        )
  }
  
  names(Clusters.List)<-paste0("G",G)
  ###### box plot generation #####
  ##### Calculation of the tightness and DB (0,1,2) indexes
  l.tight <- lapply(1:length(G),function(i){
    ClusterAll <- Clusters.List[[i]]$ClusterAll
    V<- sapply(1:(length(ClusterAll)),function(j){
      if(!is.character(ClusterAll[[j]]$Error) )
        ClusterAll[[j]]$Cl.Info$Tight.indexes
      else NA
      })
    Freq<- sapply(1:(length(ClusterAll)),function(j){
      if(!is.character(ClusterAll[[j]]$Error))
        ClusterAll[[j]]$ParamConfig.Freq
      else NA
      })
    return(data.frame(V = na.omit(V), Cluster = G[i], Freq= na.omit(Freq), Index = "Elbow plot (Tightness)"))
  })
  l.fdb <- lapply(1:length(G),function(i){
    ClusterAll <- Clusters.List[[i]]$ClusterAll
    V<- sapply(1:(length(ClusterAll)),function(j){
      if(!is.character(ClusterAll[[j]]$Error)) 
        ClusterAll[[j]]$Cl.Info$Coefficents$fDB.index
      else NA
      })
    Freq<- sapply(1:(length(ClusterAll)),function(j){
      if(!is.character(ClusterAll[[j]]$Error))
        ClusterAll[[j]]$ParamConfig.Freq
      else NA
    })
    return(data.frame(V = na.omit(V), Cluster = G[i], Freq= na.omit(Freq), Index = "fDB indexes"))
  })
  
  dt.fr<-rbind(do.call("rbind",l.tight),do.call("rbind",l.fdb))
  dt.fr.max<-aggregate(dt.fr$Freq,by = list(Cluster=dt.fr$Cluster,Index=dt.fr$Index), FUN = "max")
  colnames(dt.fr.max)[3]<-"Freq"
  dt.fr.max<-merge(dt.fr.max,dt.fr)
  
  Box.pl<-ggplot(data= dt.fr)+
    facet_wrap(~Index,scales = "free")+
    geom_boxplot(aes(x=Cluster,y=V,group=Cluster))+
    geom_line(data= dt.fr.max,aes(x=Cluster,y=V,col="Most probable"))+
    geom_point(col="red",aes(x=Cluster,y=V,size=Freq/runs) )+
    scale_x_continuous(breaks = G)+
    scale_color_manual(values = c("Most probable" = "blue") )+
    labs(size="Counts freq.",col= "",x = "Number of Clusters",y="Indexes Values")+
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
  
  ConsensusInfo<-lapply(1:length(G), function(Gind){
    ConsM.generation(Gind,Clusters.List,runs,data)
    })
  names(ConsensusInfo)<-paste0("G",G)
    
  return( list(ConsensusInfo=ConsensusInfo,BoxPlots=Box.pl,seed=seed) )
}

FCM.estimation<-function(data,G,p=5,h=NULL,Cores=1,seed=2404,tol = 0.001, maxit = 20,PercPCA=.85,runs=50)
{
  nworkers <- detectCores()
  if(nworkers<Cores) Cores <- nworkers

  database<-data$Dataset
  data.funcit <-matrix(c(database$ID,database$Vol,database$Time),ncol=3,byrow=F)
  grid <- data$TimeGrid
  
  ############### Calculation and integration of the Gauss points into the timegrid 
  ############### for calculating the distances between curves
  ## we need the splines values calculated in the gauss$nodes in the time interval [a,b]
  
  gauss.quad(10) -> gauss
  
  a <- min(grid)
  b <- max(grid)
  itempi <- (a+b)/2 + (b-a)/2*gauss$nodes
  
  grid <- sort(c(grid,itempi))
  match(itempi,grid) -> itimeindex 
  
  gauss.info<-list(gauss=gauss,itimeindex=itimeindex,a=a,b=b)
  ##########################################
  ## Let define the input parameters for the FCM function
  
  points<-database$Vol
  ID<-database$ID
  times<-database$Time
  match(times,grid) -> timeindex
  
  if(is.null(h)){
    ##########################################
    ## First run with h = min(G-1,p)
    ALL.runs = Par.fitfclust(points=points,
                             ID=ID,
                             timeindex=timeindex,
                             p=p,
                             h=min(G-1,p),
                             G=G,
                             grid=grid,
                             tol=tol,
                             maxit=maxit,
                             Cores=Cores,
                             runs=runs)
    ALL.runs.tmp<-ALL.runs
    # deleting if there was some errors in the predictions:
    ALL.runs.tmp <-lapply(1:length(ALL.runs.tmp),function(x){
      if(!is.character(ALL.runs.tmp[[x]]$Error))
        ALL.runs.tmp[[x]]
      else NA
    } )
    if(length(which(is.na(ALL.runs.tmp)))!=0){
      ALL.runs.tmp <- ALL.runs.tmp[-which(is.na(ALL.runs.tmp))]
    }
    
    ###
    alpha.List<-lapply(1:length(ALL.runs.tmp) , function(x) ALL.runs.tmp[[x]]$parameters$alpha)
    percentagePCA.List<-lapply(1:length(alpha.List),function(i){
      alpha <- alpha.List[[i]]
      #alpha <- fcm.fit$parameters$alpha
    # Principal Components Analysis considering the alpha_k
      
      princomp(as.matrix(alpha)) -> pca
      # Number of principal components
      ncomp <- length(names(pca$sdev))
      # Principal components variances
      eigs <- pca$sdev^2
      # Percentage of variances explained by each component
      percentage <- eigs/sum(eigs)
      return(percentage)
    })
    
    percentagePCA <- do.call("rbind",percentagePCA.List)
    AllPerc<-apply(percentagePCA,1,cumsum)
    if(min(G-1,p)!=1){
      h <- sapply(1:length(AllPerc[1,]) , function(i) min(which(AllPerc[,i] >= PercPCA)) )
      h.selected = as.numeric(names(table(h)[1])) 
    }else{
      h.selected =1
    }
    
    h.out <- list(h=h.selected,
                Percenteges=AllPerc)
      }
    ##########################################
    if(h.selected != min(G-1,p)){
      ### Setting the new h for estimating the new param configuration
      ALL.runs = Par.fitfclust(points,ID,
                               timeindex,p,
                               h=h.selected,
                               G,grid,tol,maxit,
                               Cores=Cores,
                               runs=runs)
  }else{
    ## run with h passed by the user
    ALL.runs = Par.fitfclust(points=points,
                             ID=ID,
                             timeindex=timeindex,
                             p=p,
                             h=h,
                             G=G,
                             grid=grid,
                             tol=tol,
                             maxit=maxit,
                             Cores=Cores,
                             runs=runs)
    h.out <- h
  }
  
  ## Check the same parameter configurations:
  Indexes.Uniq.Par<-Unique.Param(ALL.runs)
  ## Obtain the cluster and all its info

  ClusterAll<-ClusterPrediction(ALL.runs,Indexes.Uniq.Par,data,gauss.info,G)
  ####################
  ClusterAll$h.selected<-h.out
  
  return(ClusterAll)
}

Par.fitfclust = function(points,ID,timeindex,p,h,G,grid,tol,maxit,Cores=1,runs=100,seed=2404)
{
  
  cl <- makeCluster(Cores)
  clusterSetRNGStream(cl, seed)
  clusterCall(cl, function(){ 
    library(Matrix)
    library(splines)
    library(statmod)
  })
  
  
  clusterExport(cl,list("fitfclust","fclustinit","fclustMstep","fclustEstep","fclustconst",
                        "points","ID","timeindex","G","p","h","grid","tol","maxit"),envir = environment() )
  
  ALL.runs<-parLapply(cl,1:runs, function(i){
    tryCatch({
      fitfclust(x=points,
                curve=ID,
                timeindex=timeindex,
                q=p,
                h=h,
                K=G,
                p=p,
                grid=grid,
                tol = tol,
                maxit = maxit)},
      error=function(e) {
        err<-paste("ERROR in fitfclust :",conditionMessage(e), "\n")
        err.list<-list(Error= err)
        print(err)
        return(err.list)
        })
    })
  
  stopCluster(cl)
  
  return(ALL.runs)
}

Unique.Param = function(List.runs.fitfclust)
{
  L1<- length(List.runs.fitfclust)
  # deleting if there was some errors in the predictions:
  List.runs.fitfclust <-lapply(1:L1,function(x){
    if(!is.character(List.runs.fitfclust[[x]]$Error))
      List.runs.fitfclust[[x]]
    else NA
  } )
  if(length(which(is.na(List.runs.fitfclust)))!=0){
    List.runs.fitfclust <- List.runs.fitfclust[-which(is.na(List.runs.fitfclust))]
  }
  ###
  L1=length(List.runs.fitfclust)
  Indexes.Param.list<-list()
  k=1
  seq.L2 = 1:L1
  while( length(seq.L2)!=0 ){
    i = seq.L2[1]
    for(j in seq.L2){
      
       EqParam = lapply(1:length(List.runs.fitfclust[[i]]$parameters), function(par){
         trunc(List.runs.fitfclust[[i]]$parameters[[par]],prec = 4) %in% trunc(List.runs.fitfclust[[j]]$parameters[[par]],prec = 4)
         } )
     
      table(unlist(EqParam))->t.EqP
      if(length(t.EqP)==1 && as.logical(names(t.EqP[1])))
      {
        if(length(Indexes.Param.list) < k) Indexes.Param.list[[k]] <- j
        else Indexes.Param.list[[k]] <- c(Indexes.Param.list[[k]],j)
      }
    }
    seq.L2 = seq.L2[-which(seq.L2%in%Indexes.Param.list[[k]])]
    #cat(seq.L2,"\n")
    k=k+1
  }
  
  return(Indexes.Param.list)
}

ClusterPrediction = function(List.runs.fitfclust,Indexes.Uniq.Par,data,gauss.info,G)
{
  L1<- length(List.runs.fitfclust)
  # deleting if there was some errors in the predictions:
  List.runs.fitfclust <-lapply(1:L1,function(x){
    if(!is.character(List.runs.fitfclust[[x]]$Error))
      List.runs.fitfclust[[x]]
    else NA
  } )
  if(length(which(is.na(List.runs.fitfclust)))!=0){
    ErrorConfigurationFit <- List.runs.fitfclust[which(is.na(List.runs.fitfclust))]
    List.runs.fitfclust <- List.runs.fitfclust[-which(is.na(List.runs.fitfclust))]
  }else{
    ErrorConfigurationFit<-NULL
  }
  ###
  ClusterAll<-lapply(1:length(Indexes.Uniq.Par),function(i){
    tryCatch({
      ## Select just one configuration among the same ones
      Indexes.Uniq.Par[[i]][1]->j
      List.runs.fitfclust[[j]]->fcm.fit
      ## Run the FCM prediction
      fclust.curvepred(fcm.fit) -> fcm.prediction
      fclust.pred(fcm.fit) -> fcm.PRED
      fcm.PRED$class.pred -> cluster
      if(length(unique(cluster))!=G){
        err<-list(Error= paste0("ERROR in prediction: number of clusters obtained is different from ",G),
                  Params= List.runs.fitfclust[[j]]) 
        return( err )
      }
       
      ##
      
      database<-data$Dataset
      ClustCurve <- data.frame(ID=database[,1],
                               Times=database[,3],
                               Vol=database[,2],
                               Cluster= rep(cluster,data$LenCurv))
      
      ## Goodness coefficents calculation
      fcm.prediction$meancurves->meancurves
      distances <- L2dist.curve2mu(clust=cluster, fcm.curve = fcm.prediction, gauss.info = gauss.info, fcm.fit = fcm.fit, deriv = 0)
      distances.zero<-L2dist.mu20(fcm.prediction,gauss.info,fcm.fit,deriv=0)
      Coefficents<-Rclust.coeff(clust=cluster, fcm.curve = fcm.prediction, gauss.info = gauss.info, fcm.fit = fcm.fit, deriv = 0)
      Deriv.Coefficents<-Rclust.coeff(clust=cluster, fcm.curve = fcm.prediction, gauss.info = gauss.info, fcm.fit = fcm.fit, deriv = 1)
      Deriv2.Coefficents<-Rclust.coeff(clust=cluster, fcm.curve = fcm.prediction, gauss.info = gauss.info, fcm.fit = fcm.fit, deriv = 2)
      ## Let name the cluster with A->Z from the lower mean curve to the higher.
      if( G !=1 )
      {
        Cl.order<-order(distances.zero)
      }else{
        Cl.order<-1
      }
      cl.names<-LETTERS[order(Cl.order)]
      names(cluster)<-cl.names[cluster]
      
      out.funcit<-list()
      output<-list()
      
      out.funcit$cluster <- list(ClustCurve=ClustCurve,
                                 cl.info=fcm.PRED,
                                 cluster.member=cluster,
                                 cluster.names=cl.names)
      out.funcit$fit <- fcm.fit
      out.funcit$prediction <- fcm.prediction
      out.funcit$TimeGrid<-data$TimeGrid
      
      output$FCM <- out.funcit
      output$Cl.Info<- list(Tight.indexes = sum(distances),
                            Coefficents=Coefficents,
                            Deriv.Coefficents=Deriv.Coefficents,
                            Deriv2.Coefficents=Deriv2.Coefficents)
      output$ParamConfig.Freq <- length(Indexes.Uniq.Par[[i]]) 
      return(output)
    }, error=function(e) {
      err<-paste("ERROR in prediction :",conditionMessage(e), "\n")
      err.list<-list(Error= err,
                Params= List.runs.fitfclust[[j]])
      print(err)
      return(err.list)
      }
    )
  })
  
  return(list(ClusterAll=ClusterAll,ErrorConfigurationFit=ErrorConfigurationFit) )
}

ConsM.generation<-function(Gind,ALL.runs,runs,data) 
{
  FixedG.runs<-ALL.runs[[Gind]]$ClusterAll
  L1<- length(FixedG.runs)
  # deleting if there was some errors in the predictions:
  FixedG.runs.tmp <-lapply(1:L1,function(x){
    if(!is.character(FixedG.runs[[x]]$Error))
      FixedG.runs[[x]]
    else NA
  } )
  if(length(which(is.na(FixedG.runs.tmp)))!=0){
    ErrorConfiguration <- list(FromFitting=ALL.runs[[Gind]]$ErrorConfigurationFit,
                               FromPrediction = FixedG.runs[which(is.na(FixedG.runs.tmp))] )
    FixedG.runs <- FixedG.runs[-which(is.na(FixedG.runs.tmp))]
  }else{
    ErrorConfiguration <-list(FromFitting=ALL.runs[[Gind]]$ErrorConfigurationFit,
                              FromPrediction = NULL)
  }
  ###
  
  L1<- length(FixedG.runs)
  ############# first we found the most probably clustering:
  ClustCounting<-sapply(1:L1,function(x){
      names(FixedG.runs[[x]]$FCM$cluster$cluster.member  )
    } )
  Freq<-sapply(1:L1,function(x){
      FixedG.runs[[x]]$ParamConfig.Freq 
    } ) 
  
  ClustString<-sapply(1:length(ClustCounting[1,]),function(x) paste ( table(ClustCounting[,x]) , collapse = "") )
  IndexBestClustering<-which.max(Freq)
  
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
  ErrorConfiguration
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
  lab.new<-rev(sapply(1:length(G),function(i) paste(Freq.cl[i,c("Cluster","Mean")],collapse = "; ") ))
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
         subtitle = paste("Black line for the most probable clustering: ",mean(Freq.cl$Mean)) ) + 
    annotate(geom = "text", x = x.text, y = y.text, 
             label = lab.new , 
             size = 6)
  
  return(list(ConsensusMatrix = consensusM,
              ConsensusPlot = ConsensusPlot,
              MostProbableClustering = BestClustering,
              ErrorConfiguration = ErrorConfiguration) )
} 







