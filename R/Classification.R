#' Classification New Curves
#'
#' @description
#'
#' @param newdata Dataframe of three columns storing the new longitudinal data to classify into the clusters stored in *clusterdata*. The first columns, called "\emph{ID}", stores the identification number of each sample. The second column, labeled "\emph{Observation}", contains  the data observation over the time and the respective time point is reported in the third column, labeled  "\emph{Time}"
#' @param clusterdata The list obtained from extrapolating the most probable clustering from the StabilityAnalysis function output. (see \code{\link{StabilityAnalysis}} and \code{\link{clusterdata.Extrapolation}}). 
#' @param entropyCutoff  Entropy filter value for the definition of the UNCLASSIFIED cluster. The default is 1, and it works together with  probCutoff parameter.
#' @param probCutoff  Probability filter value for the definition of the UNCLASSIFIED cluster. The default is 0.6, and it works together with  entropyCutoff parameter.
#' @param Cores Number of cores to parallelize computations.

#' @description Filter formula for the definion of the UNCLASSIFIED cluster 
#'              $ Entropy >= entropyCutoff && MajorProb<=probCutoff$ then $"Unclassified"$
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @return 
#' 
#' @examples
#'
#' @import dplyr ggplot2 tidyr mvtnorm parallel
#' @export
#' 
#' 
ClassificationNewCurves<-function(newdata, clusterdata, entropyCutoff =1,probCutoff = 0.6, Cores=1)
{
  nworkers <- detectCores()
  if(nworkers<Cores) Cores <- nworkers
  grid <- clusterdata$FCM$TimeGrid
  
  ### delete points that are outside of the old grid ###
  newdata = newdata %>% filter(Time <= max(grid), Time >= min(grid))
  
  ### Lets obtain the new S of the new curves exploiting the S from FCM
  ## The old S is not the FUllS stored in clusterdata since this one is the U from the svd(old S)
  
  par <- clusterdata$FCM$fit$parameters
  p = unique(dim(par$Gamma))
  Nclusters =  length(clusterdata$FCM$cluster$cluster.names)
  
  Lambda <- par$Lambda
  alpha <- par$alpha
  lambda.zero <- as.vector(par$lambda.zero)
  Lambda.alpha <- lambda.zero + Lambda %*% t(alpha)
  Gamma <- par$Gamma
  
  cbind(1,splines::ns(grid, df = (p - 1))) -> Sold
  svdSold = svd(Sold)
  U = svdSold$u
  V = svdSold$v
  D = svdSold$d
  Dinv = 1/D *diag(nrow = p)
  
  # Now we have to find the inverse of the covariance and meancurves
  
  Lambda.alpha.new = V %*% Dinv %*% Lambda.alpha
  Gamma.new = V %*% Dinv %*% Gamma %*% Dinv %*% t(V)
  
  # Lets calculate the new S of the new curves
  
  newGrid = sort(unique(newdata$Time))
  Snew = matrix(1,ncol = p , nrow = length(newGrid))
  Snew[,2:(dim(Sold)[2])] = sapply(2: (dim(Sold)[2]), 
                                   function(i) 
                                     stats::spline(x = grid, y = Sold[,i], xout = newGrid )$y
  )
  row.names(Snew) = newGrid
  
  
  ###
  
  cl <- makeCluster(Cores)
  clusterCall(cl, function(){ 
    library(dplyr)
    library(ggplot2)
    library(mvtnorm)
    library(tidyr)
  })
  clusterExport(cl,list("clusterdata","newdata","Lambda.alpha.new","Snew", 
                        "Gamma.new","ClassificationSingleCurve"),envir = environment() )
  IDcurves = unique(newdata$ID)
  
  ALL.runs<-parLapply(cl,IDcurves, function(i){
    tryCatch({
      ClassificationSingleCurve(clusterdata,
                                newdata %>% filter(ID == i),
                                Snew = Snew,
                                Gamma.new = Gamma.new,
                                Lambda.alpha.new = Lambda.alpha.new)
    },
    error=function(e) {
      err<-paste("ERROR:",conditionMessage(e), "\n")
      err.list<-list(Error= err)
      #print(err)
      return(err.list)
    })
  })
  
  stopCluster(cl)
  
  names(ALL.runs) = paste0("ID_",IDcurves)
  df = as.data.frame(t(sapply(ALL.runs,"[[",3)),row.names = F)
  df$ID = IDcurves
  df = df[, c("ID",sort(colnames(df%>%select(-ID))) )]
  
  # Entropy calculation
  
  df_Entrop = df %>%
    tidyr::gather(-ID,key = "Cluster",value = "Prob") %>%
    group_by(ID) %>%
    mutate(Entropy = -sum(Prob*log2(Prob)),
           MajorProb = max(Prob) )%>%
    mutate(ClusterOld = Cluster,
           Cluster = ifelse( !is.na(Entropy) & (Entropy<entropyCutoff | MajorProb>probCutoff) , Cluster[which(Prob == MajorProb)], "Unclassified") ) %>%
    ungroup() %>%
    tidyr::spread(key = "ClusterOld",value = "Prob")
  
  return(list(ClassMatrix = df, ClassMatrix_entropy = df_Entrop, ListClassID =  ALL.runs ) )
}

ClassificationSingleCurve = function(clusterdata, newdata_sing, Snew, Gamma.new, Lambda.alpha.new){
  par <- clusterdata$FCM$fit$parameters
  p = unique(dim(par$Gamma))
  Nclusters =  length(clusterdata$FCM$cluster$cluster.names)
  
  # newdata_sing = newdata %>% filter(ID == 9)
  
  
  ##
  Pcl = lapply(1:Nclusters, function(i,Snew, Gamma.new, Lambda.alpha.new){
    ### Truncate the curve w.r.t the maximum time among the curves belonging the corresponding cluster.
    MaxTimeCluster = clusterdata$FCM$cluster$ClustCurve %>% filter(Cluster == i) %>% summarise(T_max = max(Time))
    newdata_sing_Cl = newdata_sing %>% dplyr::filter(Time <= MaxTimeCluster$T_max)
    
    if(length(newdata_sing_Cl$Time) < 2)
      return(data.frame(log_pi = 0, pi = 0, cluster = i) )
    
    ## defining the new spline basis matrix S
    Sx  = Snew[row.names(Snew) %in% newdata_sing_Cl$Time,]
    x = newdata_sing_Cl$Observation
    n = length(x)
    Sigma <- par$sigma * diag(n) + Sx %*% Gamma.new %*% t(Sx)
    
    mu_i =  Sx %*% Lambda.alpha.new[,i] #(par$lambda.zero + par$Lambda %*% alphai )
    x_mu = (x - mu_i)
    
    pi = mvtnorm::dmvnorm(x = x,
                          mean = mu_i,
                          sigma = Sigma )
    
    log_pi = mvtnorm::dmvnorm(x = x,
                              mean = mu_i,
                              sigma = Sigma, log = T ) + log(par$pi[i])
    
    return(data.frame(log_pi = log_pi, pi = pi, cluster = i) )
  }, Snew =Snew, Gamma.new = Gamma.new, Lambda.alpha.new = Lambda.alpha.new)
  
  Pcl =  do.call(rbind,Pcl)
  Pcl$Cluster = clusterdata$FCM$cluster$cluster.names[Pcl[,"cluster"]]
  
  ## calculate the probs to belong in the clusters
  
  Pcl$class = Pcl$pi * par$pi / sum(Pcl$pi * par$pi)
  
  newdata_sing$Cluster = Pcl$Cluster[which.max(Pcl$log_pi)]
  Pclass = Pcl$class
  names(Pclass) = Pcl$Cluster
  
  ### Plotting with the new
  
  ClustCurve = clusterdata$FCM$cluster$ClustCurve
  ClustCurve$Cluster = clusterdata$FCM$cluster$cluster.names[ClustCurve$Cluster]
  ClustCurve = ClustCurve %>% group_by(Cluster) %>% mutate(Tmax = max(Time)) %>% ungroup()
  
  MeanCurve = as.data.frame(clusterdata$FCM$prediction$meancurves)
  colnames(MeanCurve) = clusterdata$FCM$cluster$cluster.names
  MeanCurve$Time = clusterdata$CONNECTORList$TimeGrid
  MeanCurve = MeanCurve %>% tidyr::gather(-Time, value = "Mean",key = "Cluster")
  
  MeanCurve = merge(MeanCurve,ClustCurve %>% select(Cluster,Tmax) %>% distinct()) %>% group_by(Cluster) %>% filter(Time <=Tmax) %>% ungroup()
  #### plot generation
  
  ProbAnnot = data.frame(Cluster = Pcl$Cluster,
                         Prob = round(Pcl$class,digits = 2),
                         x = rep(-Inf,Nclusters),
                         y = rep(Inf,Nclusters))
  
  pl = ggplot()+
    geom_line(data = MeanCurve,aes(x = Time, y = Mean,linetype=Cluster)) +
    geom_line(data = ClustCurve,aes(x = Time, y = Observation, group = ID), col = "grey",alpha = .4) +
    facet_grid(~Cluster)+
    geom_line(data = newdata_sing,aes(x = Time, y = Observation),col = "red") + 
    geom_text(data = ProbAnnot, aes(x = x, y = y, label = paste0("Prob: ",Prob)),hjust = 0, vjust = 1)
  
  return(list(plot = pl, weight = Pcl, prob = Pclass))
}



