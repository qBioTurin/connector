% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExtrapolationFunctions.R
\name{ExtrapolationFuncs}
\alias{ExtrapolationFuncs}
\alias{IndexesPlot.Extrapolation}
\alias{ConsMatrix.Extrapolation}
\alias{MostProbableClustering.Extrapolation}
\title{Extrapolate objects from ClusterAnalysis}
\usage{
IndexesPlot.Extrapolation(stability.list, q = NULL)

ConsMatrix.Extrapolation(stability.list, q = NULL)

MostProbableClustering.Extrapolation(stability.list, G, q = NULL)
}
\arguments{
\item{stability.list}{The list obtained from the ClusterAnalysis function. (see \code{\link{ClusterAnalysis}})}

\item{q}{The quantiles used to calculate the time grid interval on which the distances are calculated. If NULL then the time grid outliers will be ignored through the distance calculation. If double (0<q<1) then the cutting is symmetrical w.r.t. the quantile setled (e.g., q = 0.25). If a vector, then the minimum value is used for the lower cutting and the maximum value for the upper cutting.}

\item{G}{The number of clusters.}
}
\description{
Extrapolate objects from ClusterAnalysis output list (see \code{\link{ClusterAnalysis}}).
}
\details{
\itemize{
\item{IndexesPlot.Extrapolation}{extrapolates from ClusterAnalysis output list the box plot fixing the h value.}
\item{ConsMatrix.Extrapolation}{extrapolates from ClusterAnalysis output list the Stability Matrix for G fixed.}
\item{MostProbableClustering.Extrapolation}{extrapolates from ClusterAnalysis output list the most frequent clustering among the several runs obtained  for G fixed.}
}
}
\author{
Cordero Francesca, Pernice Simone, Sirovich Roberta
}
