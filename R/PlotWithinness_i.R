#' i-th Cluster Withinness Plot
#'
#' Plots the i-th cluster withinness measure circles.
#'
#' @param ClustCurve Data frame with five arguments : time, growth values, ID, cluster membership and feature values for each sample.
#' @param MeanCurves Matrix with the meancurve values on the columns according to different clusters.
#' @param i Value for the cluster involved in withinness computation.
#' @param centroids If TRUE the WithCluster_MeanDist function is used, otherwise withinness is calculated using WithCluster_CurvDist.
#' @param shift Vector representing the distances between the first withinness circle center and the other centers. 
#' @return Returns a plot representing the i-th cluster withinness measures circles.
#' 
#' @import ggforce
#' @export
#'
PlotWithinness.i <- function(ClustCurve,MeanCurves,i,shift=0)
{
  dataplot <- DataFrameWithinness.i(ClustCurve,MeanCurves,i,shift=shift)
  ### Data frame for plot
  circles <- dataplot$circles
  WithDist <- dataplot$WithDist
  feature <- WithDist[,3]
  feature.name <- colnames(WithDist)[3]
  nfeature <- length(unique(ClustCurve[,feature.name]))[1]
  feature.palette <- rainbow(nfeature)
  ### Plot distances
  plot.i <- ggplot() + geom_circle(aes(x0=x0, y0=y0, r=r), data=circles,size=1)
  plot.i <- plot.i  + geom_point(data=WithDist,aes(x=x1,y=y1,col=feature),shape=(i-1)) + labs(title=paste("Cluster",ClustSymbol[i],"withinness"),x="distance",y="distance")
  plot.i <- plot.i + scale_colour_manual(values = feature.palette[unique(feature)],name=feature.name) + xlab("distance x") + ylab("distance y")
  return(plot.i)
 }
