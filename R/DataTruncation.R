#' DataTruncation
#'
#' @description
#'
#' Truncates the functional data (the time series) at a specific time point chosen by the user.
#'
#' @param data CONNECTORList. (see \code{\link{DataImport}})
#' @param feature The column name reported in the AnnotationFile containing the feature  to be investigated.
#' @param truncTime  A two dimension vector of integers corresponding to the time points where the curves will be truncated. If an integer number is passed, than it will be considered as the upper time point by default.
#' @param labels  Vector containing the text for the title of axis names and plot title.
#' @param save If TRUE then the plot of the time series truncated at the "TruncTime" is saved into a pdf file.
#' @param path The folder path where the plot(s) will be saved. If it is missing, the plot is saved in the current working  directory.
#' 
#' @return  DataTruncation returns the line plot of time series with a vertical line at the truncation time and the CONNECTORList updated with the following elements: (i) a data frame with three variables (ID curves, observation and time values truncated at the chosen time), (ii) a vector collecting the number of truncated observations collected per sample, a data frame matching curses  with the chosen feature, (iv) the vector storing all the truncated time points of the samples (i.e. truncated time grid).
#' Furthermore, it prints an updated summary of the input data, updating the minimum and the maximum curve length after the data truncation.
#'
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#' 
#' @examples
#'
#'TimeSeriesFile<-"data/475dataset.xls"
#'AnnotationFile <-"data/475info.txt"
#'
#'CONNECTORList <- DataImport(TimeSeriesFile,AnnotationFile)
#'
#'CONNECTORList<- DataTruncation(CONNECTORList,"Progeny",truncTime=50,labels = c("time","volume","Tumor Growth"))
#'CONNECTORList<- DataTruncation(CONNECTORList,"Progeny",truncTime=c(20,50),labels = c("time","volume","Tumor Growth"))
#'
#' @import ggplot2
#' @export
DataTruncation <- function(data,feature,truncTime=NULL,labels=NULL,save=FALSE,path=NULL)
{
  lencurve.before<-data$LenCurv
  
  ### Variables initialization
  growth.curve.ls <- PlotTimeSeries(data,feature,labels = labels)
  ### Plot curves with truncation time
  growth.curve.tr <- growth.curve.ls$PlotTimeSeries_plot
  
  if(! is.null(truncTime))
    growth.curve.tr <- growth.curve.tr + geom_vline(xintercept=truncTime, color="black", size=1)
  
  if(save==TRUE)
  {
    if(is.null(path))
    {
      path <- getwd()
    }
    ggsave(filename="DataTruncation.pdf",plot =growth.curve.tr,width=29, height = 20, units = "cm",scale = 1,path=path )
  }
  
  # Data truncation
  if(!is.null(truncTime)) 
    data.tr <- DataTrunc(data,truncTime=truncTime)
  else 
    data.tr <- data
  
  if(is.double(data.tr) )
  {
    stop("Considering this truncation time some curves with just one point are present.\n  A larger truncation time must be chosen.\n")
  }
  
  data.tr[["ColFeature"]]<-growth.curve.ls$data$FeatureColour
  data.tr[["PlotTimeSeries_plot"]]<-growth.curve.tr
  
  cat("############################################################## \n######## Summary of the trunc. data ############\n")
  cat("\n Number of curves:",length(data.tr$LenCurv),";\n Min curve length: ",min(data.tr$LenCurv),"; Max curve length: ",max(data.tr$LenCurv),";\n")
  
  Curve.cutted.ind<-which(!data.tr$LenCurv==lencurve.before)
  points.deleted<-lencurve.before[Curve.cutted.ind]-data.tr$LenCurv[Curve.cutted.ind]
  
  cat("\n Number of truncated curves:",length(Curve.cutted.ind),";\n Min points deleted: ",min(points.deleted),"; Max points deleted: ",max(points.deleted),";\n")
  
  cat("############################################################## \n")
  
  
  return(data.tr)
}

#######################à

DataTrunc <- function(data,truncTime=NULL)
{
  #### Truncates the growth data at a specific time point chosen by the user.
  # Variables inizialization
  dataset <- data$Dataset
  sample.size <- length(unique(dataset[,1]))
  lencurv.tr <- numeric(sample.size)
  
  # Data truncation
  
  if(!is.null(truncTime))
  {
    max.time<-max(dataset[,3])
    min.time<-min(dataset[,3])
    
    if(length(truncTime)>1) 
    {
      max.truncTime <- max(truncTime)
      min.truncTime <- min(truncTime)
    }else{
      min.truncTime<- min.time
      max.truncTime <- truncTime
    }
    
    if(max.time<max.truncTime)  warning("Max truncation time greater than maximum time in the dataset.")
    if(min.time>min.truncTime)  warning("Min truncation time smaller than minimum time in the dataset.")
    
    dataset.tr <- dataset[dataset[,3]<=max.truncTime & dataset[,3]>=min.truncTime ,]
    
    for (i in 1:sample.size)  lencurv.tr[i]<-length(dataset.tr[,1][dataset.tr[,1]==i])
    timegrid.tr <- data$TimeGrid[data$TimeGrid<=max.truncTime & data$TimeGrid>=min.truncTime ]
    data.tr=list(Dataset=dataset.tr,LenCurv=lencurv.tr,LabCurv=data$LabCurv,TimeGrid=timegrid.tr)
  }
  else data.tr <- data
  
  if( length(which(data.tr$LenCurv < 2))!=0 ) 
  {
    warning("Curves with less than 2 points are now present!!!")
  }
  
  return(data.tr)
  ##### The CONNECTORList updated with the following arguments: a data frame with three variables (ID curves, observation and time values truncated at the chosen time), a vector collecting the number of truncated observations collected per sample, a data frame with curves labeled according to target file feature chosen and a vector for overall truncated time grid.
}
