#' Curve-Curve withinness 
#'
#' @description
#'
#' Computes the Hausdorff distance among curves belonging to the i-th cluster.
#'
#' @param ClustCurve Data frame with five arguments : time, growth values, ID, cluster membership and feature values for each sample.
#' @param i Value included from 0 to K (number of clusters) representing the cluster involved in withinness computation.
#' @return Returns a vector reporting the hausdorff distance among the  curves belonging to the i-th cluster.
#' @examples
#' @export
WithCluster_CurvDist <- function(ClustCurve,i)
{

  ### i-th cluster curves data
  ClustCurve.i <- ClustCurve[ClustCurve[,4]==i,]
  ### i-th cluster curves ID
  index <- sort(unique(ClustCurve.i[,1]))
  ### i-th cluster magnitudo
  ni <- length(index)

  if(ni!=1)
   {
   within.i <- numeric((1/2)*(factorial(ni)/(factorial(ni-2))))
   count <- 0
   for (k in 1:(ni-1))
    {
      A <- ClustCurve.i[ClustCurve.i[,1]==index[k],2:3]
      for (j in (k+1):ni)
       {
         B <- ClustCurve.i[ClustCurve.i[,1]==index[j],2:3]
         count <- count +1
         within.i[count] <- hausdorff(A,B)
       }
    }
   }
  else within.i <- 0
  ### i-th cluster withinness
  return(within.i)
}

