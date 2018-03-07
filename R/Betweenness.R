#' Betweenness
#'
#' Computes the Betweenness across all clusters.
#'
#' @param ClustCurve Data frame with five arguments : time, growth values, ID, cluster membership and feature values for each sample.
#' @param MeanCurves Matrix with the meancurve values on the columns according to different clusters.
#' @return Returns a list with three arguments: (i) Betweenness, a matrix composed by two columns, the minimum cluster distance and the cluster name that achieves the minimum, (ii) CentroidDist, a K x K matrix (where K is the cluster number) containing centroids cluster distance among each cluster couple, (iii) Classes, a vector containing the cluster name membership for each dataset sample.
#' 
#' @export
Betweenness <- function(ClustCurve,MeanCurves)
{
  K <- length(unique(ClustCurve[,4]))
  ClustSymbol <- cluster.symbol(K)
  ClassCurve <- unique(ClustCurve[,c(1,4)])[,2]
  classes <- ClustSymbol[ClassCurve]
  centroid.dist <- matrix(numeric(K*K),nrow=K)
  between <- matrix(numeric(K*2),ncol=2)
  colnames(between) <- c("Min distance","Nearest cluster")
  rownames(between) <- ClustSymbol
  rownames(centroid.dist) <- ClustSymbol
  colnames(centroid.dist) <- rownames(centroid.dist)

  for(i in 1:K)
    {
	centroid.dist[i,] <-  BetweenCluster_MeanDist(ClustCurve,MeanCurves,i)
	ClustSymbol.i <- ClustSymbol[-i]
	min.dist <- min(centroid.dist[i,][centroid.dist[i,]!=0])
	argmin.dist <- ClustSymbol.i[which.min(centroid.dist[i,][centroid.dist[i,]!=0])]
	between[i,] <- c(min.dist,argmin.dist)
	}
 return(list(Betweenness=between,CentroidDist=centroid.dist,Classes=classes))
}
