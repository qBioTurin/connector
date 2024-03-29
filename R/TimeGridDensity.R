#' Time Grid Density
#'
#'@description
#'Plots the grid of the time points at which the growth data are collected. Each point of the grid is colored with respect to its frequency in the input data.
#'
#' @param  data CONNECTORList. (see \code{\link{DataImport}})
#' @param  save If TRUE then the density time grid plot is saved into a pdf file.
#' @param  path	 The folder path where the plot will be saved. If it is missing, the plot is saved in the current working  directory.
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @return Time Grid Density returns the time grid density plot as a ggplot object, in which a point $p_{x,y}$ is defined by a pair of coordinates $p_{x,y}=\left( x,y\right) \ $ and by a colour. $p_{x,y}$ is defined if only if exists at least one sample with two observations at time $x\ $ and $y$.
#'  The colour associates with it encodes the frequency of samples in which $p_{x,y}$ is present.
#' 
#' @examples
#' GrowDataFile<-"data/745dataset.xls"
#' AnnotationFile <-"data/745info.txt"
#'
#' CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#' TimeGridDensity(CONNECTORList)
#'
#' @import ggplot2 grid dplyr
#' @importFrom ggplotify as.ggplot 
#' @importFrom gridExtra arrangeGrob
#' @export
TimeGridDensity <- function(data,save=FALSE,path=NULL)
{

 ### Variables initialization
 TimeMeasure <-data$Dataset[,c(1,3)]
 SampleSize <- length(unique(data$Dataset[,1]))
 SampleID <- unique(data$Dataset[,1])
 LenCurve <- data$LenCurv
 PointsCoord<-matrix(0,nrow =LenCurve%*%LenCurve,ncol = 2)
 k<-0

 # Generate time grid
  for(i in SampleID)
  {
    l<-length(TimeMeasure[TimeMeasure[,1]==i,1])
    PointsCoord[(k+1):(k+l^2),1]<-rep(TimeMeasure[TimeMeasure[,1]==i,2],l)
    PointsCoord[(k+1):(k+l^2),2]<-rep(TimeMeasure[TimeMeasure[,1]==i,2],each=l)
    k<-k+l^2
  }

  ### Generate data frame for ggplot
  Time1 <- PointsCoord[,1]
  Time2 <- PointsCoord[,2]
  
  df <- data.frame(Time1=Time1,Time2=Time2)
  df1 <- df %>%
    group_by(Time1, Time2, .drop = F)  %>% 
    tally %>%
    ungroup() %>%
    as.data.frame()
  

  ### Plot density grid
#   TimeGrid_plot <- ggplot(df,aes(Time1, Time2)) +
#                     stat_density_2d( aes(fill=stat(nlevel)), geom = "polygon" )+
#                     scale_fill_gradientn(colours = c("#baffc9","#FF0000"),
#                                          name="Freq. of \nobservations")+
#                   theme_bw() +
# 	                labs(title="Time grid",x="Time", y = "Time")+
#                   theme(plot.title = element_text(hjust = 0.5),
#                         title =element_text(size=12, face='bold'))

  
  TimeGrid_plot <-  ggplot() + 
    geom_point(data=df1,aes(Time1, Time2,col=n/max(n)))+
    stat_density_2d(data=df,aes(Time1, Time2,fill=stat(nlevel)), geom = "polygon" )+
    scale_fill_gradientn(colours = c("#baffc9","#FF0000"),
                         name="Freq. of \nobservations")+
    scale_color_gradientn(colours = c("#baffc9","#FF0000"),
                         name="Freq. of \nobservations")+
    theme_bw() +
    labs(title="Grid of time points",x="Time", y = "Time")+
    theme(plot.title = element_text(hjust = 0.5),
          title =element_text(size=12, face='bold'))+
    guides(color = "none")
  
  BarPlot<-ggplot(data=data.frame(x=data$Dataset$Time), aes(x) )+
    geom_boxplot(width=0.4)+ coord_flip()+
    #labs(x="",y="Number of observations per time point")+
    labs(x="",y="Time points distribution")+
    theme(legend.text=element_text(size=14),
          legend.position=c(.9,1.1),
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.key.size = unit(.9, "cm"),
          legend.key.width = unit(.9,"cm"),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          plot.margin=unit(c(10,5,5,5),"mm"), 
          strip.background = element_blank(),
          strip.text.x = element_blank() ,
          title =element_text(size=12, face='bold'))

  Plot.new <-arrangeGrob(cbind(ggplotGrob(TimeGrid_plot+
                                             theme(legend.position="left") ) ,
                                           ggplotGrob(BarPlot),size = "last")) 
  
  Plot.new <-ggplotify::as.ggplot(Plot.new)
  
  if(save==TRUE)
  {
    if(is.null(path))
	{
	path <- getwd()
    ggsave(filename="TimeGrid.pdf",plot =TimeGrid_plot,width=29, height = 20, units = "cm",scale = 1,path=path )
    }
  }
  return(list(TimeGrid_plot=TimeGrid_plot,TimeGrid_Boxplot=Plot.new))
}
