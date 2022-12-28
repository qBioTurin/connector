#' Cluster Analysis
#' 
#'@description
#'
#'  Fits and clusters the data with respect to the Functional Clustering Model [Sugar and James]. Multiple runs of the algorithm are necessary since the algorithm is stochastic As explained in [Sugar and James], to have a simple low-dimensional representation of the individual curves and to reduce the number of parameters to be estimated, h value must be equals or lower than \eqn{min(p,G-1)}.
#'  
#' @param data CONNECTORList. (see \code{\link{DataImport}} or \code{\link{DataTruncation}})
#' @param G  The vector/number of possible clusters.
#' @param p The dimension of the natural cubic spline basis. (see \code{\link{BasisDimension.Choice}})
#' @param runs Number of runs.
#' @param seed Seed for the kmeans function.
#' @param save If TRUE then the growth curves plot truncated at the "TruncTime" is saved into a pdf file.
#' @param path The folder path where the plot(s) will be saved. If it is missing, the plot is saved in the current working  directory.
#' @param Cores Number of cores to parallelize computations.
#' @param PercPCA=.85
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @return StabilityAnalysis returns a list of (i) lists, called ConsensusInfo, reporting for each G and h: the Consensus Matrix, either as a NxN matrix, where N is the number of samples, or plot, and the most probable clustering obtained from running several times the method; (ii) the box plots showing both the Elbow plot considering the total tightness and the box plots of the fDB indexes for each G; and finally, (iii) the seed. See \code{\link{IndexesPlot.Extrapolation}} and \code{\link{MostProbableClustering.Extrapolation}}.
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
#' @import RColorBrewer statmod parallel Matrix splines 
#' @export

ClusterAnalysis<-function(data,G,p,h=NULL,runs=50,seed=2404,save=FALSE,path=NULL,Cores=1,PercPCA=.85,MinErrFreq= 0)
{
  ## Let define some parameters needed to FCM
  params <- list()
  database<-data$Dataset
  data.funcit <-matrix(c(database$ID,database$Observation,database$Time),ncol=3,byrow=F)
  params$grid <- grid <- data$TimeGrid
  params$points <- database$Observation
  params$ID <- database$ID
  params$timeindex <- match(database$Time,grid)
  #######
  h.gBefore <- NULL
  Clusters.List<-list()
  for( g in 1:length(G) ){
    FCM.Cluster <- FCM.estimation(data = data,
                                  G = G[g],
                                  p = p,
                                  h.user=h,
                                  seed=seed,
                                  PercPCA= PercPCA,
                                  Cores=Cores,
                                  runs = runs,
                                  params = params,
                                  #gauss.infoList = gauss.infoList,
                                  h.gBefore = h.gBefore
    )
    if(is.character(FCM.Cluster)){
      # If NULL it means that it was not found any errors during the model runs. Otherwise the max h value without errors is used thereafter.
      Clusters.List[[g]] = FCM.Cluster
    }else{
      # If NULL it means that it was not found any errors during the model runs. Otherwise the max h value without errors is used thereafter.
      h.gBefore <- FCM.Cluster$h.gBefore
      Clusters.List[[g]] = FCM.Cluster$ClusterAll
    }
  }
  
  names(Clusters.List)<-paste0("G",G)
  
  return(list(Clusters.List=Clusters.List, CONNECTORList = data, seed=seed,runs = runs))
}

FCM.estimation<-function(data,G,params,gauss.infoList=NULL,h.gBefore,p=5,h.user=NULL,Cores=1,seed=2404,tol = 0.001, maxit = 20,PercPCA=.85,runs=50,MinErrFreq= 0)
{
  nworkers <- detectCores()
  if(nworkers<Cores) Cores <- nworkers
  
  if(is.null(h.user)){
    ##########################################
    ## First run with h = min(G-1,p) or with the h from the smaller G (the upper bound of the errors runs was found)
    if(is.null(h.gBefore)){
      h.selected = min(G-1,p)
    }else{
      h.selected = h.gBefore
    }
    
    h.found = F
    tentative = 1
    while(!h.found)
    {
      ALL.runs = Par.fitfclust(points=params$points,
                               ID=params$ID,
                               timeindex=params$timeindex,
                               p=p,
                               h=h.selected,
                               G=G,
                               grid=params$grid,
                               tol=tol,
                               maxit=maxit,
                               Cores=Cores,
                               runs=runs)
      
      ## Check the same parameters configurations:
      Indexes.Uniq.Par<-Unique.Param(ALL.runs)
      if(is.character(Indexes.Uniq.Par)){
        return(Indexes.Uniq.Par)
      }
      
      ## Let calculate the clustering and the fDB indexes
      ClusterAll<-ClusterPrediction(ALL.runs,Indexes.Uniq.Par,data,gauss.infoList,G)
      # l.fdb <- fdb.calc(ClusterAll$ClusterAll)
      # rep(l.fdb$V,l.fdb$Freq) -> fDBindexes
      
      ## Check the number of errors:
      ALL.runs.tmp <-lapply(1:length(ALL.runs),function(x){
        if(!is.character(ALL.runs[[x]]$Error))
          ALL.runs[[x]]
        else NA
      } )
      ClusterAll.tmp <-lapply(1:length(ClusterAll$ClusterAll),function(x){
        if(!is.character(ClusterAll$ClusterAll[[x]]$Error))
          ClusterAll$ClusterAll[[x]]
        else NA
      } )
      
      if(length( which(is.na(ClusterAll.tmp)) )!=0 ) {
        if (length(which(is.na(ClusterAll.tmp))) == length(ClusterAll.tmp) ) {
          errors = unique(sapply(1:length(ClusterAll$ClusterAll), function(x) ClusterAll$ClusterAll[[x]]$Error))
          return(paste("All runs have errors, that are the following:\n ",errors) )
        }else{
          Nerr.Clus <- sum(sapply( which(is.na(ClusterAll.tmp)) ,
                                   function(i) length(Indexes.Uniq.Par[[i]] )))
        }
      }else{
        Nerr.Clus <- 0
      }
      Nerr.Fitting <- length(which(is.na(ALL.runs.tmp)))
      Err.Freq <- (Nerr.Fitting+Nerr.Clus)/runs
      
      if(MinErrFreq > 1 ){
        warning("MinErrFreq paramter should belong to [0,1]. The default is used.")
        MinErrFreq <- 0
        
      }
      
      # h reached the minimum value or the errors freq is higher than the set MinErrFreq
      if(Err.Freq <= MinErrFreq | h.selected == 1 ){
        h.found = T
        h.out = h.selected
        if(tentative > 1){
          h.gBefore <- h.selected
        }
      }else{
        h.selected = h.selected - 1
        tentative = tentative + 1
      }
      
    }
    
  }else{
    ## run with h passed by the user
    ALL.runs = Par.fitfclust(points=params$points,
                             ID=params$ID,
                             timeindex=params$timeindex,
                             p=p,
                             h=h.user,
                             G=G,
                             grid=params$grid,
                             tol=tol,
                             maxit=maxit,
                             Cores=Cores,
                             runs=runs)
    h.out <- h.user 
    ## Check the same parameter configurations:
    Indexes.Uniq.Par<-Unique.Param(ALL.runs)
    ## Obtain the cluster and all its info
    ClusterAll<-ClusterPrediction(ALL.runs,Indexes.Uniq.Par,data,gauss.infoList,G)
  }
  
  ####################
  ClusterAll$h.selected<-h.out
  
  return(list(h.gBefore = h.gBefore, ClusterAll = ClusterAll) )
}

Par.fitfclust = function(points,ID,timeindex,p,h,G,grid,tol,maxit,Cores=1,runs=100,seed=2404)
{
  
  if(Cores == 1){
    ALL.runs<-lapply(1:runs, function(i,points,ID,timeindex,G,p,h,grid,tol,maxit){
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
          #print(err)
          return(err.list)
        })
    },points,ID,timeindex,G,p,h,grid,tol,maxit)
  }else{

  type <- if (exists("mcfork", mode="function")) "FORK" else "PSOCK"
  cl <- makeCluster(Cores, type = type)
  clusterSetRNGStream(cl, seed)
  clusterCall(cl, function(){ 
    library(Matrix)
    library(splines)
    library(statmod)
  })
  
  
  clusterExport(cl,list("fitfclust","fclustinit","fclustMstep","fclustEstep","fclustconst",
                        "points","ID","timeindex","G","p","h","grid","tol","maxit"),
                envir = environment() )
  
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
        #print(err)
        return(err.list)
      })
  })
  
  stopCluster(cl)
  
  }
  
  return(ALL.runs)
}

Unique.Param = function(List.runs.fitfclust)
{
  List.runs.fitfclust.tmp = List.runs.fitfclust
  L1 <- length(List.runs.fitfclust)
  List.runs.fitfclust <- lapply(1:L1, function(x) {
    if (!is.character(List.runs.fitfclust[[x]]$Error)) 
      List.runs.fitfclust[[x]]
    else NA
  })
  if (length(which(is.na(List.runs.fitfclust))) == length(List.runs.fitfclust)) {
    errors = unlist(unique(List.runs.fitfclust.tmp))
    return(paste("All runs have errors, that are the following:\n ",errors) )
  }
  
  # if(length(which(is.na(List.runs.fitfclust)))!=0){
  #   List.runs.fitfclust <- List.runs.fitfclust[-which(is.na(List.runs.fitfclust))]
  # }
  ###
  L1=length(List.runs.fitfclust)
  Indexes.Param.list<-list()
  k=1
  seq.L2 = 1:L1
  
  if(length(which(is.na(List.runs.fitfclust)))>0){
    seq.L2 <- seq.L2[-which(is.na(List.runs.fitfclust))]
  }
  
  while( length(seq.L2)!=0 ){
    i = seq.L2[1]
    for(j in seq.L2){
      EqParam = lapply(1:length(List.runs.fitfclust[[i]]$parameters), function(par){
        trunc(List.runs.fitfclust[[i]]$parameters[[par]],prec = 3) %in% trunc(List.runs.fitfclust[[j]]$parameters[[par]],prec = 4)
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

ClusterPrediction = function(List.runs.fitfclust,Indexes.Uniq.Par,data,gauss.infoList,G)
{
  L1<- length(List.runs.fitfclust)
  # deleting if there was some errors in the predictions:
  List.runs.fitfclust <-lapply(1:L1,function(x){
    if(!is.character(List.runs.fitfclust[[x]]$Error))
      List.runs.fitfclust[[x]]
    else NA
  } )
  # if(length(which(is.na(List.runs.fitfclust)))!=0){
  #   ErrorConfigurationFit <- List.runs.fitfclust[which(is.na(List.runs.fitfclust))]
  #   List.runs.fitfclust <- List.runs.fitfclust[-which(is.na(List.runs.fitfclust))]
  # }else{
  #   ErrorConfigurationFit<-NULL
  # }
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
                               Time=database[,3],
                               Observation=database[,2],
                               Cluster= rep(cluster,data$LenCurv))
      
      ## Goodness coefficents calculation
      fcm.prediction$meancurves->meancurves
      #distances <- L2dist.curve2mu(clust=cluster, fcm.curve = fcm.prediction, database = database, fcm.fit = fcm.fit, deriv = 0)
      distances.zero<-L2dist.mu20(clust=cluster,fcm.prediction,database = database,fcm.fit,deriv=0,q=NULL)
      #Coefficents<-Rclust.coeff(clust=cluster, fcm.curve = fcm.prediction, database = database, fcm.fit = fcm.fit, deriv = 0)
      #Deriv.Coefficents<-Rclust.coeff(clust=cluster, fcm.curve = fcm.prediction, database = database, fcm.fit = fcm.fit, deriv = 1)
      #Deriv2.Coefficents<-Rclust.coeff(clust=cluster, fcm.curve = fcm.prediction, database = database, fcm.fit = fcm.fit, deriv = 2)
      
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
      # output$Cl.Info<- list(Tight.indexes = sum(distances),
      #                       Coefficents=Coefficents,
      #                       Deriv.Coefficents=Deriv.Coefficents,
      #                       Deriv2.Coefficents=Deriv2.Coefficents)
      output$ParamConfig.Freq <- length(Indexes.Uniq.Par[[i]]) 
      return(output)
    }, error=function(e) {
      err<-paste("ERROR in prediction :",conditionMessage(e), "\n")
      err.list<-list(Error= err,
                     Params= List.runs.fitfclust[[j]])
      #print(err)
      return(err.list)
    }
    )
  })
  
  return(list(ClusterAll=ClusterAll) )
}

fdb.calc <-function(ClusterAll){
  
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
  
  return(data.frame(V = na.omit(V),
                    Freq= na.omit(Freq))
  )
}



