#' Data Visualitation
#'
#' Visualization of the time grid and the line plot of cancer growth data.
#'
#'
#' @param data CONNECTORList.
#' @param feature The column name reported in the AnnotationFile containing the feature interesting for the user to be investigated.
#' @param labels Vector containing the text for the title of axis and plot title.
#' @param save If TRUE the plot is saved in a pdf file.
#' @param path Path to save plot to (combined with file name). If it is missing, the plot is saved in the working directory.
#' @return A plot with the density time grid and the line plot of cancer growth data as a ggplot object.
#' @example
#'
#'GrowDataFile<-"data/1864dataset.xls"
#'AnnotationFile <-"data/1864info.txt"
#'
#'CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#'
#'DataVisualization(CONNECTORList,"Progeny",labels = c("time","volume","Tumor Growth"))
#'
#' @import ggplot2 cowplot
#' @export
DataVisualization <- function(data,feature,labels=NULL,save=FALSE,path=NULL)
{

 ### Variables initialization
 growth.curves <- GrowthCurve(data,feature,labels=labels)
 plot1 <- growth.curves$GrowthCurve_plot
 plot2 <- TimeGridDensity(data)

 plots <- plot_grid(plotlist=list(plot1 , plot2))

 if(save==TRUE)
 {
     ggsave(filename="DataVisualization.pdf",plot =plots,width=29, height = 20, units = "cm",scale = 1,path=path )
 }
 return(plots)
}

