#' Data Import
#'
#'@description
#' Reads the files and creates a list storing all the information.
#'
#'
#' @param TimeSeriesFile The name of the excel file storing the  growth evolution data. See \sQuote{Details}.
#'
#' @param AnnotationFile The name of a txt/csv file storing  the annotation data. See \sQuote{Details}.
#'
#' @return   DataImport returns a list, called ConnectorList, with four arguments: (i) Dataset: data frame with three columns (i.e. ID, data and time values) encoding the growth data for each sample, (ii) LenCurv: the vector reporting the number of observations collected per sample  (iii) LabCurv: the data frame matching the samples with their annotations, (iv) TimeGrid: the vector storing all the sample time points (i.e. time grid). Furthermore, it prints a brief summary of the input data, i.e. the total number of curves (samples), the minimum and the maximum curve length.
#' 
#' @details Two files are requested to run the data analysis exploiting the CONNECTOR package:
#' \itemize{
#' \item the excel file, namely \emph{TimeSeriesFile}, reporting the growth evolution data,
#' \item the csv file, namely \emph{AnnotationFile}, containing the annotation information associated with the samples.
#' }
#' Hence, the growth data associated with an experiment must be stored into TimeSeriesFile as a table with two columns for each sample.
#' The first column, labeled  \emph{Time}, contains the time points of a sample. The second column, labeled by the sample name, contains  the data observation over the time.
#' 
#' Instead, the second file (i.e. AnnotationFile)  stores  the annotated information associated with the samples as a table in csv format so that  number of rows is equals to the total number of samples.
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#' 
#' @examples
#'
#' TimeSeriesFile <- "data/475dataset.xls"
#' AnnotationFile <- "data/475info.txt"
#'
#' CONNECTORList <- DataImport(TimeSeriesFile,AnnotationFile)
#'
#' @import readxl dplyr
#' @export
DataImport <- function(TimeSeriesFile,AnnotationFile) {
  ###Read Data File
  dataset <- read_excel(TimeSeriesFile,col_names=T)
  if( all(dataset[1,]==rep(0,length(dataset[1,]))) )
  {
    warning("The first line of the TimeSeriesFile is zero.",immediate. = T)
  }
  
  ### Read Target File
  labcurv  <- read.csv(file=AnnotationFile,header=TRUE)
  colnames(labcurv)= c("IDSample",colnames(labcurv)[-(1)] )
  labcurv$ID = 1:length(labcurv$IDSample)
  labcurveID <-  labcurv$ID
  names(labcurveID) <- labcurv$IDSample
  
  ###Check the column names
  c_names<-colnames(dataset[2*(1:(length(dataset[1,])/2))])
  if(length(c_names)!=(length(labcurv$ID)))
  {
    warning("Number of columns in the excel file is different from the number of samples stored in the target file.")
    return()
  }else{
    
    if(all(c_names==labcurv$IDSample)==FALSE)
    {
      warning("SampleNames in the target file do not correspond to the names in the excel file.")
      return()
    }
  }
  
  ### Inizialize :
  dataset<-as.matrix(dataset)
  
  ### vector for curves lenghts
  nvar       <- dim(dataset)[2]
  nobs       <- dim(dataset)[1]
  
  samplesize <- nvar/2
  lencurv    <- numeric(samplesize)
  ### vectors for times, observation and ID curves values
  
  TimeIndex  <- seq(1,nvar,2)
  ObsIndex   <- seq(2,nvar,2)
  
  lencurv.time<-sapply(TimeIndex, function(x){length(dataset[!is.na(dataset[,x]),x])})
  lencurv.obs<-sapply(ObsIndex, function(x){length(dataset[!is.na(dataset[,x]),x])})
  
  if(!all(lencurv.time==lencurv.obs)) # check if times are without observation
  {
    
    ind.diff<-which(!lencurv.time==lencurv.obs)
    lencurv.obs.diff<-lencurv.obs[ind.diff]
    lencurv.time.diff<-lencurv.time[ind.diff]
    
    for(i in 1:length(ind.diff) )
    {
      na.obs<-is.na(dataset[,ind.diff[i]*2 ] )
      na.time<-is.na(dataset[,ind.diff[i]*2-1 ] )
      i.different<-which(!(na.time==na.obs))
      dataset[i.different,ind.diff[i]*2 ]<-NA
      dataset[i.different,ind.diff[i]*2-1 ]<-NA
    }
    
    warning(paste(length(ind.diff),"time(s) without the corresponding observation (or viceversa) is(are) present.\n These cases will be deleted.\n") )
  }
  
  lencurv<-sapply(ObsIndex, function(x){length(dataset[!is.na(dataset[,x]),x])})
  
  TimeValue  <- as.vector(dataset[,TimeIndex])
  TimeValue<- as.double(TimeValue[!is.na(TimeValue)])
  
  ObsValue  <- as.vector(dataset[,ObsIndex])
  ObsValue<-  as.double(ObsValue[!is.na(ObsValue)] )
  
  IDsample<- rep(c_names,times=lencurv)
  timegrid <- sort(unique(TimeValue))
  
  ### ID, observation and time data frame
  
  dataset <- data.frame(ID=IDsample,Observation=ObsValue,Time=TimeValue)
  dataset$ID = labcurveID[dataset$ID] 
  
  dataset = dataset %>% 
    arrange(ID,Time)
  
  check.times<-aggregate(dataset$Time, by = list(dataset$ID),function(x){
    table(x)->number.times
    if(length(unique(number.times))>1) return("Multiple Times")
    else return("Perfect")
  })
  
  if(length(unique(check.times$x))>1){
    warning("Samples with multiple observations at the same time point are present. One observation per time point is required")
    return()
  }
  
  alldata <- list(Dataset=dataset,LenCurv=lencurv,LabCurv=labcurv,TimeGrid=timegrid)
  
  cat("############################### \n######## Summary ##############\n")
  cat("\n Number of curves:",samplesize,";\n Min curve length: ",min(lencurv),"; Max curve length: ",max(lencurv),".\n")
  cat("############################### \n")
  
  return(alldata)
}
