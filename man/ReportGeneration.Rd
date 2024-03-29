% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReportGeneration.R
\name{ReportGeneration}
\alias{ReportGeneration}
\title{Report generation}
\usage{
ReportGeneration(
  data,
  stability.list,
  G,
  feature = "ID",
  namefile = "Report",
  path = NULL,
  truncTime = NULL,
  p = NULL,
  pRange = NULL
)
}
\arguments{
\item{data}{CONNECTORList. (see \code{\link{DataImport}})}

\item{stability.list}{The list obtained from the ClusterAnalysis function. (see \code{\link{ClusterAnalysis}})}

\item{G}{Vector of integers representing the number of clusters.}

\item{feature}{The column name reported in the AnnotationFile  containing the feature to be investigated.}

\item{namefile}{Report name.}

\item{path}{The folder path  where the report will be saved. If it is missing, the report is saved in the current working  directory.}

\item{truncTime}{A two dimension vector of integers corresponding to the time points where the curves will be truncated. If an integer number is passed, than it will be considered as the upper time point by default.}

\item{p}{The dimension of the natural cubic spline basis. (see \code{\link{BasisDimension.Choice}})}

\item{pRange}{The vector of the dimension of the natural cubic spline basis, by which the cross-validated loglikelihood is calculated and plotted. It could be the path of the Rds file storing the output of the \code{\link{BasisDimension.Choice}} function. If NULL, the p selection step is omitted.}
}
\value{

}
\description{
Report generation
}
\examples{

}
\seealso{
\code{\link{GrowthCurve}}, code{\link{TimeGridDensity}}.
}
\author{
Cordero Francesca, Pernice Simone, Sirovich Roberta
}
