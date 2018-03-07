#' Clustering
#'
#' @description
#' Fitting the cancer growth data with respect to the following growth function: Malthus, Gompertz and Logistic. The cancer growth fitted functions are then clustered.
#'
#' @param data CONNECTORList.
#' @param k  Number of clusters.
#' @param model Model name, i.e. Malthus, Gompertz and Logistic.
#' @return A list containing per each curve (i) the parameters of the fitting, (ii) the centers, (iii) the cluster of affinity and (iv) the mean values.
#' @examples
#'
#' @import grofit
#' @export
clustering<- function(data,k,model)
{

  data->dati
  time <- dati$Dataset$Time
  t <- sort(unique(time))
  x <- dati$Dataset$Vol
  curve <- dati$Dataset$ID
  number_curves<-max(curve)

  data.fit <- matrix(c(curve,x,time),ncol=3,byrow=F)
  parameters <- matrix(0,ncol=3,nrow=number_curves)
  meancurves <- matrix(numeric(length(t)*k),ncol=k)

  if(model=="Malthus")
  {
    parameters <- matrix(0,ncol=2,nrow=number_curves)
  }

  for (i in c(1:number_curves))
  {
    m <- which(data.fit[,1] == i)
    tdata <- time[m]
    Voldata <- x[m]


    # fitgompertz
    if(model=="Gompertz")
    {
      lsgompertz <- function(x,t,y)
      {
        sum(((x[1]*exp(-exp((x[2]*exp(1)/x[1])*(x[3]-t)+1)))-y)^2)
      }
      initgompertz(tdata,Voldata,1,1,1) -> initgomp
      initgomp<- list(A=initgomp$A,mu=initgomp$mu,lambda=initgomp$lambda)
      parameters[i,]<- optim(initgomp,lsgompertz,method="L-BFGS-B",lower=c(0,0,0),t=tdata,y=Voldata)$par
    }

    # fitlogistic
    if(model=="Logistic")
    {
      lslogistic <- function(x,t,y)
      {
        sum(((x[1]/(1+exp((4*x[2]/x[1])*(x[3]-t)+2)))-y)^2)
      }
      initlogistic(tdata,Voldata,1,1,1) -> initlog
      initlog<- list(A=initlog$A,mu=initlog$mu,lambda=initlog$lambda)
      parameters[i,]<- optim(initlog,lslogistic,method="L-BFGS-B",lower=c(0,0,0),t=tdata,y=Voldata)$par

    }

    # fitMalthus
    if(model=="Malthus")
    {
      malthus<-function(x,t)
      {
        x[1]*exp(x[2]*t)
      }
      lsmalthus <- function(x,t,y)
      {

        sum(((x[1]*exp(x[2]*t))-y)^2)
      }
      par <- matrix(numeric(200),ncol=2)
      res <- numeric(100)
      for(j in 1:100)
      {
        fitmal <- optim(abs(rnorm(2)),lsmalthus,t=tdata,y=Voldata)
        par[j,] <- fitmal$par
        res[j] <- sum(abs(malthus(par[j,],tdata)-Voldata)^2)
      }
      pos <- which(res==min(res))
      parameters[i,]<- par[pos,]
    }
  }

  set.seed(2404)
  group <- kmeans(parameters,k)
  cluster <- group$cluster
  center <- group$centers

###### meancurves calculating ######

  for(i in c(1:k))
  {
    center[i,]->x
    if(model=="Gompertz")
      {
        meancurves[,i] <- (x[1]*exp(-exp((x[2]*exp(1)/x[1])*(x[3]-t)+1)))
      }
    if(model=="Logistic")
      {
        meancurves[,i] <- (x[1]/(1+exp((4*x[2]/x[1])*(x[3]-t)+2)))
      }
    if(model=="Malthus")
      {
        meancurves[,i] <- x[1]*exp(x[2]*t)
      }
  }

  out<-list(parameters=parameters,cluster=cluster,center=center,meancurves=meancurves)
  return(out)
}
