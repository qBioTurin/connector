#' DataTrunc
#' @description
#'
#' Truncates the cancer growth data at a specific time point chosen by the user.
#'
#' @param data CONNECTORList.
#' @param truncTime  An integer number corresponding to the time at which truncate the curves.
#' @return  The CONNECTORList updated with the following arguments: a data frame with three variables (ID curves, volume and time values truncated at the chosen time), a vector collecting the number of truncated observations collected per sample, a data frame with curves labeled according to target file feature chosen and a vector for overall truncated time grid.
#' @example
#'
#'GrowDataFile<-"data/1864dataset.xls"
#'AnnotationFile <-"data/1864info.txt"
#'
#'CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#'
#'CONNECTORList<- DataTrunc(CONNECTORList,truncTime=60)
#'
#' @export
DataTrunc <- function(data,truncTime=NULL)
{
   # Variables inizialization
   dataset <- data$Dataset
   sample.size <- max(unique(dataset[,1]))
   lencurv.tr <- numeric(sample.size)

   # Data truncation

  if(!is.null(truncTime))
  {
    max.time<-max(dataset[,3])
    if(max.time<truncTime)  warning("Truncation time greater than maximum time in the dataset.")
   dataset.tr <- dataset[dataset[,3]<=truncTime,]
   for (i in 1:sample.size)  lencurv.tr[i]<-length(dataset.tr[,1][dataset.tr[,1]==i])
   timegrid.tr <- data$TimeGrid[data$TimeGrid<=truncTime]
   data.tr=list(Dataset=dataset.tr,LenCurv=lencurv.tr,LabCurv=data$LabCurv,TimeGrid=timegrid.tr)
  }
  else data.tr <- data

  return(data.tr)
}
