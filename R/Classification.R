#' Classification New Curves
#'
#' @description
#'
#' @param newdata 
#' @param clusterdata The list obtained from extrapolating the most probable clustering from the StabilityAnalysis function output. (see \code{\link{StabilityAnalysis}} and \code{\link{clusterdata.Extrapolation}}). 
#' @param Cores Number of cores to parallelize computations.

#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @return 
#' 
#' @examples
#'
#' @import dplyr ggplot2 tidyr mvtnorm parallel
#' @export
#' 
ClassificationNewCurves<-function(newdata, clusterdata, Cores=1)
{
  nworkers <- detectCores()
  if(nworkers<Cores) Cores <- nworkers
  
  cl <- makeCluster(Cores)
  clusterCall(cl, function(){ 
    library(dplyr)
    library(ggplot2)
    library(mvtnorm)
    library(tidyr)
  })
  clusterExport(cl,list("clusterdata","newdata","ClassificationSingleCurve"),envir = environment() )
  IDcurves = unique(newdata$ID)
  
  ALL.runs<-parLapply(cl,IDcurves, function(i){
    tryCatch({
      ClassificationSingleCurve(clusterdata,
                                newdata %>% filter(ID == i),
                                FullS = T)
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
  
  return(list(ClassMatrix = df, ListClassID =  ALL.runs ) )
}


ClassificationSingleCurve = function(clusterdata,newdata,FullS =F){
  par <- clusterdata$FCM$fit$parameters
  p = unique(dim(par$Gamma))
  Nclusters =  length(clusterdata$FCM$cluster$cluster.names)
  
  ## defining the new spline basis matrix S
  Fullgrid = unique(clusterdata$CONNECTORList$Dataset$Time)
  if(FullS){
    grid = unique(c(Fullgrid,newdata$Time) )
  }else{
    grid = unique(c(min(Fullgrid),newdata$Time,max(Fullgrid)) )
  }
  grid = sort(grid)
  Sfull <- cbind(1, splines::ns(grid, df = (p - 1)))
  S <- svd(Sfull)$u
  Sij = S[which(grid %in% newdata$Time),]
  xij = newdata$Observation
  n = length(xij)
  ##
  Pcl = sapply(1:Nclusters, function(i){
    alphai = matrix(par$alpha[i,],ncol = 1)
    Sigma <- par$sigma * diag(n) + Sij %*% par$Gamma %*% t(Sij)
    mu_i =  Sij %*% (par$lambda.zero + par$Lambda %*% alphai )
    x_mu = (xij - mu_i)
    if(par$pi[i] != 0)
      log_pi = - 1/2 * t(x_mu) %*% solve(Sigma) %*% x_mu + log(par$pi[i])
    else 
      log_pi = 0
    
    return(log_pi)
  })
  
  names(Pcl) = clusterdata$FCM$cluster$cluster.names
  newdata$Cluster = clusterdata$FCM$cluster$cluster.names[which.max(Pcl)]
  
  ## calculate the probs to belong in the clusters
  
  
  ### Plotting with the new
  
  ClustCurve = clusterdata$FCM$cluster$ClustCurve
  ClustCurve$Cluster = clusterdata$FCM$cluster$cluster.names[ClustCurve$Cluster]
  MeanCurve = as.data.frame(clusterdata$FCM$prediction$meancurves)
  colnames(MeanCurve) = clusterdata$FCM$cluster$cluster.names
  MeanCurve$Time = clusterdata$CONNECTORList$TimeGrid
  MeanCurve = MeanCurve %>% tidyr::gather(-Time, value = "Mean",key = "Cluster")
  
  
  #### Prob calculation
  
  Prob = sapply(1:Nclusters, function(i){
    alphai = matrix(par$alpha[i,],ncol = 1)
    Sigma <- par$sigma * diag(n) + Sij %*% par$Gamma %*% t(Sij)
    mu_i =  Sij %*% (par$lambda.zero + par$Lambda %*% alphai)
    pi = mvtnorm::dmvnorm(x = xij,
                          mean = mu_i,
                          sigma = Sigma )
    return(pi)
  })
  
  Pclass = Prob * par$pi / sum(Prob * par$pi)
  
  names(Pclass) = clusterdata$FCM$cluster$cluster.names
  
  #### plot generation
  
  ProbAnnot = data.frame(Cluster = clusterdata$FCM$cluster$cluster.names,
                         Prob = round(Pclass,digits = 2),
                         x = rep(-Inf,Nclusters),
                         y = rep(Inf,Nclusters))
  
  pl = ggplot()+
    geom_line(data = MeanCurve,aes(x = Time, y = Mean,linetype=Cluster),alpha = .4) +
    geom_line(data = ClustCurve,aes(x = Time, y = Observation, group = ID, col = as.factor(ID)),alpha = .4) +
    facet_grid(~Cluster)+
    geom_line(data = newdata,aes(x = Time, y = Observation),col = "red") + 
    geom_text(data = ProbAnnot, aes(x = x, y = y, label = paste0("Prob: ",Prob)),hjust = 0, vjust = 1)
  
  return(list(plot = pl, weight = Pcl,prob = Pclass))
}



