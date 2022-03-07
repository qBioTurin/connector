#' Data Visualitation
#'
#' Computes the time grid and the line plot of growth data.
#'
#'
#' @param data CONNECTORList. (see \code{\link{DataImport}})
#' @param feature The column name reported in the TimeSeriesFile  containing the feature to be investigated.
#' @param labels The vector containing the text for the axis names and plot title.
#' @param save If TRUE then the plot is saved into a pdf file.
#' @param path The folder path  where the plot will be saved. If it is missing, the plot is saved in the current working  directory.
#' 
#' @return  Data Visualization returns a plot with the density time grid and the line plot of growth data as a ggplot object.
#'  In details, a point $p_{x,y}$ of the time grid density  is defined by a pair of coordinates $p_{x,y}=\left( x,y\right) \ $ and by a colour. $p_{x,y}$ is defined if only if exists at least one sample with two observations at time $x\ $ and $y$.
#'  The colour associates with it encodes the frequency of samples in which $p_{x,y}$ is present.
#' 
#' @examples
#'
#'TimeSeriesFile<-"data/745dataset.xls"
#'AnnotationFile <-"data/745info.txt"
#'
#'CONNECTORList <- DataImport(TimeSeriesFile,AnnotationFile)
#'
#'DataVisualization(CONNECTORList,"Progeny",labels = c("time","volume","Tumor Growth"))
#'
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#' 
#' @seealso  \code{\link{PlotTimeSeries}}, code{\link{TimeGridDensity}}.
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @export
DataVisualization <- function(data,feature,labels=NULL,save=FALSE,path=NULL)
{
  
  ### Variables initialization
  TimeSeries.curves <- PlotTimeSeries(data,feature,labels=labels)
  plot1 <- TimeSeries.curves$PlotTimeSeries_plot
  plot2 <- TimeGridDensity(data)
  
  plots <- plot_grid(plotlist=list(plot1 , plot2$TimeGrid_plot))
  
  if(save==TRUE)
  {
    ggsave(filename="DataVisualization.pdf",plot =plots,width=29, height = 20, units = "cm",scale = 1,path=path )
  }
  return(plots)
}

