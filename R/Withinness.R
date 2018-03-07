#' Overall Clusters withinness
#'
#' Compute the withinness across all clusters obtained.
#'
#' @param ClustCurve A data frame with five arguments: time, growth values, ID, cluster membership and feature values for each sample. 
#' @param MeanCurves Matrix with the meancurve values on the columns according to different clusters.
#' @param centroids If TRUE the WithCluster_MeanDist function is used, otherwise withinness is calculated using WithCluster_CurvDist.
#' @return Returns a matrix with four columns: mean, standard deviation, minimum and maximum withinness distance across the clusters.
#' @examples
#' @export
Withinness <- function(ClustCurve,MeanCurves,centroids=TRUE)
{
  K <- length(unique(ClustCurve[,4]))
  ClustSymbol <- cluster.symbol(K)
  ### Withinness matrix
  withinness <- matrix(numeric(K*4),ncol=K)
  rownames(withinness) <- c("mean","sd","min","max")
  colnames(withinness) <- paste("Cluster ",ClustSymbol,sep="")

  for (i in 1:K)
  {
    if(centroids==FALSE) ### i-th cluster curves distance
	{
	  within.i <- WithCluster_CurvDist(ClustCurve,i)
	}
    else                 ### i-th cluster curves and meancurve distance
	  within.i <- WithCluster_MeanDist(ClustCurve,MeanCurves,i)

    ### i-th cluster mean distance
    MeanDist <- mean(within.i)
	### i-th cluster standard deviation distance
    if(length(within.i)==1) { StDev <- 0 }
    else                        StDev <- sd(within.i)
	MinDist <- min(within.i)
	MaxDist <- max(within.i)
	### i-th cluster withinness data
    withinness[,i] <- cbind(MeanDist,StDev,MinDist,MaxDist)
  }
 return(withinness)
}
