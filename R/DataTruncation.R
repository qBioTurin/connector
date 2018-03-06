#' DataTruncation
#'
#' @description
#'
#' Truncates the cancer growth data at a specific time point chosen by the user.
#'
#' @param data CONNECTORList.
#' @param feature The column name reported in the AnnotationFile containing the feature interesting for the user to be investigated.
#' @param truncTime  An integer number corresponding to the time at which truncate the curves.
#' @param labels   Vector containing the text for the title of axis and plot title.
#' @param save If TRUE the plot is saved in a pdf file.
#' @param path Path to save plot to (combined with file name). If it is missing, the plot is saved in the working directory.
#' @return DataTruncation returns the line plot of growth curves with a vertical line at the truncation time and the CONNECTORList updated with the following arguments: a data frame with three variables (ID curves, volume and time values truncated at the chosen time), a vector collecting the number of truncated observations collected per sample, a data frame with curves labeled according to target file feature chosen and a vector for overall truncated time grid.
#'
#' @example
#'
#'GrowDataFile<-"data/1864dataset.xls"
#'AnnotationFile <-"data/1864info.txt"
#'
#'CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#'
#'CONNECTORList<- DataTruncation(CONNECTORList,"Progeny",truncTime=60,labels = c("time","volume","Tumor Growth"))
#'
#' @import ggplot2
#' @export
DataTruncation <- function(data,feature,truncTime=NULL,labels=NULL,save=FALSE,path=NULL)
{

### Variables initialization
growth.curve.ls <- GrowthCurve(data,feature,labels = labels)
### Plot growth curves with truncation time
growth.curve.tr <- growth.curve.ls$GrowthCurve_plot
if(! is.null(truncTime)) growth.curve.tr <- growth.curve.tr + geom_vline(xintercept=truncTime, color="black", size=1)

if(save==TRUE)
{
 if(is.null(path))
 {
  path <- getwd()
 }
 ggsave(filename="DataTruncation.pdf",plot =growth.curve.tr,width=29, height = 20, units = "cm",scale = 1,path=path )
}

### Truncated dataset
dataset <- data$Dataset
sample.size <- max(unique(dataset[,1]))
lencurv.tr <- numeric(sample.size)

# Data truncation
if(!is.null(truncTime)) data.tr <- DataTrunc(data,truncTime=truncTime)
else data.tr <- data

data.tr$LabCurv <- data.tr$LabCurv[c("ID",feature)]

plot(growth.curve.tr)

return(data.tr)
}
