#' Fitting and Clustering
#'
#'
#' @description
#' Fitting and Clustering the cancer growth data with respect to the following growth functions: Malthus, Gompertz and Logistic. The cancer growth fitted functions are then clustered.
#'
#' @param data CONNECTORList.
#' @param k  Number of clusters.
#' @param model Model name, i.e. Malthus, Gompertz and Logistic.
#' @param fitting.method 
#' @param lower
#' @param upper
#' @param seed
#' @return A list containing per each curve (i) the parameters of the fitting, (ii) the centers, (iii) the cluster of affinity and (iv) the mean values.
#' @examples
#'
#' @import GenSA DEoptim optimr
#' @export
FittingAndClustering <- function(data,k,model, fitting.method="optimr",init=NULL,lower,upper,seed=NULL)
{

  data->dati
  time <- dati$Dataset$Time
  t <- sort(unique(time))
  x <- dati$Dataset$Vol
  curve <- dati$Dataset$ID
  number_curves<-max(curve)
  
  data.fit <- matrix(c(curve,x,time),ncol=3,byrow=F)
  meancurves <- matrix(numeric(length(t)*k),ncol=k)
  
  
if((length(lower)!=3 || length(upper)!=3 ))
{
    if(model=="Malthus" & (length(lower)!=2 || length(upper)!=2 ) ) 
    {
      stop("Considering the Malthus model, the length of buonds should be 2.")
    }
    else{
      stop("Considering the Logistic or Gompertz model, the length of buonds should be 3.")
    }
}
 
if(fitting.method!="DEoptim")
{
  if((length(init)!=3 ))
  {
    if(model=="Malthus" & (length(init)!=2 ) ) 
      {
        stop("Considering the Malthus model, the length of the initial vector should be 2.")
      }else{
  stop("Considering the Logistic or Gompertz model, the length of the initial vector should be 3.")
}
}
}


    
  
  if(!is.null(seed))
  {
    set.seed(seed)
  }
  
  
  # fitgompertz
  if(model=="Gompertz")
  {
    lsfunction <- function(x,t,y)
    {
      sum((x[1]*exp(x[2]/x[3]*(1-exp(-x[3]*t))) - y )^2)
      # sum(((x[1]*exp(-exp( (x[2]*exp(1)/x[1])*(x[3]-t) + 1) ) )-y)^2)
    }
    
  }
  
  # fitlogistic
  if(model=="Logistic")
  {
    lsfunction <- function(x,t,y)
    { 
      sum( ( ( x[1]*x[2]/( x[1] + (x[2]-x[1])*exp(- x[3]*t) ) ) - y )^2)
      # sum(((x[1]/(1+exp((4*x[2]/x[1])*(x[3]-t)+2)))-y)^2)
    }
    
  }
  
  # fitMalthus
  if(model=="Malthus")
  {
    
    lsfunction <- function(x,t,y)
    {
      sum(((x[1]*exp(x[2]*t))-y)^2)
    } 
    
  }
  
  
  ######################  FITTING THE DATA
  
  if(fitting.method == "GenSA"){
    
    fn<-function(i){
      
      m <- which(data.fit[,1] == i)
      tdata <- time[m]
      Voldata <- x[m]
      GenSA(unlist(init), lsfunction,lower=lower,upper=upper,t=tdata,y=Voldata)$par
    }
    
  }
  
  if(fitting.method == "optimr"){
    
    fn<-function(i){
      
      m <- which(data.fit[,1] == i)
      tdata <- time[m]
      Voldata <- x[m]
      # optim(init,lsfunction,method="L-BFGS-B",lower=lower,t=tdata,y=Voldata)$par
      params<-opm(unlist(init),lsfunction,lower=lower,upper=upper,method = "ALL",t=tdata,y=Voldata)
      ##### let take the params combination with minor value of SSE
      params<-summary(params,order=value)
      ##### just the convent method
      covergent<-params[params$convergence==0,]
      
      if(empty(covergent)) {
        warning("Using optmir, no methods are successful convergent!")
      }else{   
        covergent[1,1:length(init)] 
      }
      
    }
    
  }
  

  if(fitting.method == "DEoptim"){
    
    fn<-function(i){
      
      m <- which(data.fit[,1] == i)
      tdata <- time[m]
      Voldata <- x[m]
      DEoptim( lsfunction, lower=lower,upper=upper,t=tdata,y=Voldata,control = DEoptim.control(trace = F, itermax = 2000))$optim$bestmem
    }
    
  }
  
  parameters<- t(sapply(1:number_curves,fn))
  
  x<-parameters[1,]
  # set.seed(2404)
  group <- kmeans(parameters,k,iter.max = 50)
  cluster <- group$cluster
  center <- group$centers
  curves <- data.frame(ID=data$Dataset$ID, Times=data$Dataset$Time, Vol=data$Dataset$Vol, Cluster= rep(cluster,data$LenCurv))
  
  
  seed <- .Random.seed
  ###### meancurves calculating ######
  
  for(i in c(1:k))
  {
    center[i,]->x
    if(model=="Gompertz")
    {
      # meancurves[,i] <- (x[1]*exp(-exp((x[2]*exp(1)/x[1])*(x[3]-t)+1)))
      meancurves[,i] <-x[1]*exp(x[2]/x[3]*(1-exp(-x[3]*t)))
    }
    if(model=="Logistic")
    {
      # meancurves[,i] <-(x[1]/(1+exp((4*x[2]/x[1])*(x[3]-t)+2)))
      meancurves[,i] <-x[1]*x[2]/( x[1] + (x[2]-x[1])*exp(-x[3]*t) )
    }
    if(model=="Malthus")
    {
      meancurves[,i] <- x[1]*exp(x[2]*t)
    }
  }
  
 OUT<-list(parameters=parameters,cluster=cluster,center=center,meancurves=meancurves,Summary=curves,seed=seed)

  return(OUT)
}
