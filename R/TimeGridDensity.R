#' Time Grid Density
#'
#'@description
#'Plots the grid of the time points at which the growth data are collected.
#'Each point of the grid is colored with respect to its frequency in the input data.
#'
#' @param  data CONNECTORList.
#' @param  save If TRUE it is saved in a pdf the density time grid plot.
#' @param  path	Path to save plot to (combined with file name). If it is missing, the plot is saved in the working  directory.
#' @return The time grid density plot as a ggplot object.
#' @examples
#' GrowDataFile<-"data/1864dataset.xls"
#' AnnotationFile <-"data/1864info.txt"
#'
#' CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#' TimeGridDensity(CONNECTORList)
#'
#' @import ggplot2
#' @export
TimeGridDensity <- function(data,save=FALSE,path=NULL)
{

 ### Variables initialization
 TimeMeasure <-data$Dataset[,c(1,3)]
 SampleSize <- length(unique(data$Dataset[,1]))
 LenCurve <- data$LenCurv
 PointsCoord<-matrix(0,nrow =LenCurve%*%LenCurve,ncol = 2)
 k<-0

 # Generate time grid
  for(i in 1:SampleSize)
  {
    l<-length(TimeMeasure[TimeMeasure[,1]==i,1])
    PointsCoord[(k+1):(k+l^2),1]<-rep(TimeMeasure[TimeMeasure[,1]==i,2],l)
    PointsCoord[(k+1):(k+l^2),2]<-rep(TimeMeasure[TimeMeasure[,1]==i,2],each=l)
    k<-k+l^2
  }

  ### Generate data frame for ggplot
  Time1 <- PointsCoord[,1]
  Time2 <- PointsCoord[,2]
  df <- data.frame(Time1=Time1,Time2=Time2,
  d  <- densCols(Time1, Time2, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

  ### Plot density grid
  TimeGrid_plot <- ggplot(df) +
                   geom_point(aes(Time1, Time2,col=d), size = 4) +
	               coord_fixed(ratio = 1) +
                   scale_color_identity() +
                   theme_bw() +
	               labs(title="Time grid",x="Time", y = "Time")+
                   theme(plot.title = element_text(hjust = 0.5),title =element_text(size=12, face='bold'))

  if(save==TRUE)
  {
    if(is.null(path))
	{
	path <- getwd()
    ggsave(filename="TimeGrid.pdf",plot =TimeGrid_plot,width=29, height = 20, units = "cm",scale = 1,path=path )
    }
  }
  return(TimeGrid_plot)
}
