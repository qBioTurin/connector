#' Withinness data frame
#'
#' Creates a dataframe for i-th cluster withinness measures.
#'
#' @param ClustCurve Data frame with five arguments : time, growth values, ID, cluster membership and feature values for each sample.
#' @param MeanCurves Matrix with the meancurve values on the columns according to different clusters.
#' @param i Value included from 0 to K (number of clusters) representing the cluster involved in withinness computation.
#' @param ClustSymbol Vector of the cluster symbols.
#' @param shift Vector representing the distances between the first withinness circle center and the other centers.
#' @return Returns a list composed by two data frames: (i) "circles" containing coordinates and radius for i-th cluster withinness circles measures, (ii) "WithDist" containing the i-th cluster curves withinness distance from i-th cluster centroid/meancurve.
#' 
#' @export
DataFrameWithinness.i <- function(ClustCurve,MeanCurves,i,ClustSymbol,shift=0)
{
  K <- length(unique(ClustCurve[,4]))

  ### Centroid withinness distance
  Withinness.i <- WithCluster_MeanDist(ClustCurve,MeanCurves,i)

  ### Mean and standard deviation distance
  mean.dist <- mean(Withinness.i)
  st.dev <- sd(Withinness.i)

  ### Data frame for ggplot
  circles.i <- data.frame(
  x0 <- rep(0+shift,3),
  y0 <- rep(0,3),
  r <- rep(mean.dist, 3)+c(st.dev,0,-st.dev),
  distance <- c("sd","mean","sd"),
  Cluster <- factor(ClustSymbol[i],levels=ClustSymbol)
  )
  colnames(circles.i) <- c("x0","y0","r","distance","Cluster")

  WithDist.i <- data.frame(
  x1 <- Withinness.i + shift,
  y1 <- numeric(length(Withinness.i)),
  Cluster <- factor((i-1),levels=c(0:(K-1)))
  )
  colnames(WithDist.i) <- c("x1","y1","Cluster")
 return(DataFrame.i=list(circles=circles.i,WithDist=WithDist.i))
 }



