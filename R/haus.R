#' Hausdorff distance
#'
#' Calculates the Hausdorff distance between two sets.
#'
#' @param A A set.
#' @param B A set.
#' @return Returnst the Hausdorff distance between A and B.
#' 
#' @export

hausdorff <- function(A,B)
{
  H <- max(distance(A,B),distance(B,A))
  return(H)
}

eucl.dist <- function(x) {sqrt(sum(x^2))}

distance <- function(A,B)
{
   n <- nrow(A)
   m <- nrow(B)
   a <- numeric(n)
   b <- numeric(m)

   for (i in 1:n)
      {
        for (j in 1:m)
        {
          b[j] <- eucl.dist(A[i,]-B[j,])
        }
          a[i] <- min(b)
      }
   return(max(a))
}
