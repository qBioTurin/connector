#' Centroid betweenness
#'
#' Calculates the Hausdorff distance between the i-th cluster meancurve and the other clusters meancurves.
#'
#' @param ClustCurve Data frame with five arguments : time, growth values, ID, cluster membership and feature values for each sample.
#' @param MeanCurves Matrix with the meancurve values on the columns according to different clusters.
#' @param i Value included from 0 to K (number of clusters) representing the cluster involved in withinness computation.
#' @return Returns a numeric vector for i-th cluster meancurve hausdorff distance from other clusters meancurves.
#' @examples
#' @export
BetweenCluster_MeanDist <- function(ClustCurve,MeanCurves,i)
{
  K <- length(unique(ClustCurve[,4]))
  ClustSymbol <- cluster.symbol(K)
  TimeGrid <- sort(unique(ClustCurve[,2]))
  ### i-th cluster curves data
  ClustCurve.i <- ClustCurve[ClustCurve[,4]==i,]
  ### i-th cluster max obs time
  tmax <- max(ClustCurve.i[,2])
  ### i-th cluster meancurve truncated at tmax
  MeanCurve.i <- MeanCurves[TimeGrid <= tmax,i]
  ### i-th cluster obs time grid
  TimeGrid.i <- TimeGrid[TimeGrid <= tmax]
  ### i-th cluster meancurve
  A <- cbind(TimeGrid.i,MeanCurve.i)
  ### other clusters
  other.cluster <- sort(unique(ClustCurve[-which(ClustCurve[,4]==i),4]))

  ### Betweenness centroid dist.curveance
  betweencentroid.i <- matrix(numeric(length(other.cluster)),nrow=1)
  count <- 1
  for (j in other.cluster)
   {
    ClustCurve.j <- ClustCurve[ClustCurve[,4]==j,]
    tmax.j <- max(ClustCurve.j[,2])
    MeanCurve.j <- MeanCurves[TimeGrid <= tmax.j,j]
    TimeGrid.j <- TimeGrid[TimeGrid <= tmax.j]
    B <- cbind(TimeGrid.j,MeanCurve.j)
    betweencentroid.i[count] <- hausdorff(A,B)
	count <- count +1
   }
   between <- matrix(numeric(K),nrow=1)
   between[,-i] <- betweencentroid.i
   colnames(between) <- paste(ClustSymbol,"dist")
  return(BetweenCentr=between)
}
