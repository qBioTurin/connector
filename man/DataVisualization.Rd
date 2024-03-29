% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataVisualization.R
\name{DataVisualization}
\alias{DataVisualization}
\title{Data Visualitation}
\usage{
DataVisualization(data, feature, labels = NULL, save = FALSE, path = NULL)
}
\arguments{
\item{data}{CONNECTORList. (see \code{\link{DataImport}})}

\item{feature}{The column name reported in the TimeSeriesFile  containing the feature to be investigated.}

\item{labels}{The vector containing the text for the axis names and plot title.}

\item{save}{If TRUE then the plot is saved into a pdf file.}

\item{path}{The folder path  where the plot will be saved. If it is missing, the plot is saved in the current working  directory.}
}
\value{
Data Visualization returns a plot with the density time grid and the line plot of growth data as a ggplot object.
In details, a point $p_{x,y}$ of the time grid density  is defined by a pair of coordinates $p_{x,y}=\left( x,y\right) \ $ and by a colour. $p_{x,y}$ is defined if only if exists at least one sample with two observations at time $x\ $ and $y$.
The colour associates with it encodes the frequency of samples in which $p_{x,y}$ is present.
}
\description{
Computes the time grid and the line plot of growth data.
}
\examples{

TimeSeriesFile<-"data/745dataset.xls"
AnnotationFile <-"data/745info.txt"

CONNECTORList <- DataImport(TimeSeriesFile,AnnotationFile)

DataVisualization(CONNECTORList,"Progeny",labels = c("time","volume","Tumor Growth"))

}
\seealso{
\code{\link{PlotTimeSeries}}, code{\link{TimeGridDensity}}.
}
\author{
Cordero Francesca, Pernice Simone, Sirovich Roberta
}
