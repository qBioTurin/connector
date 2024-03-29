% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataClustering.R
\name{ClusterAnalysis}
\alias{ClusterAnalysis}
\title{Cluster Analysis}
\usage{
ClusterAnalysis(
  data,
  G,
  p,
  h = NULL,
  runs = 50,
  seed = 2404,
  save = FALSE,
  path = NULL,
  Cores = 1,
  PercPCA = 0.85,
  MinErrFreq = 0,
  pert = 0.01
)
}
\arguments{
\item{data}{CONNECTORList. (see \code{\link{DataImport}} or \code{\link{DataTruncation}})}

\item{G}{The vector/number of possible clusters.}

\item{p}{The dimension of the natural cubic spline basis. (see \code{\link{BasisDimension.Choice}})}

\item{runs}{Number of runs.}

\item{seed}{Seed for the kmeans function.}

\item{save}{If TRUE then the growth curves plot truncated at the "TruncTime" is saved into a pdf file.}

\item{path}{The folder path where the plot(s) will be saved. If it is missing, the plot is saved in the current working  directory.}

\item{Cores}{Number of cores to parallelize computations.}

\item{pert}{....}
}
\value{
StabilityAnalysis returns a list of (i) lists, called ConsensusInfo, reporting for each G and h: the Consensus Matrix, either as a NxN matrix, where N is the number of samples, or plot, and the most probable clustering obtained from running several times the method; (ii) the box plots showing both the Elbow plot considering the total tightness and the box plots of the fDB indexes for each G; and finally, (iii) the seed. See \code{\link{IndexesPlot.Extrapolation}} and \code{\link{MostProbableClustering.Extrapolation}}.
}
\description{
Fits and clusters the data with respect to the Functional Clustering Model \link{Sugar and James}. Multiple runs of the algorithm are necessary since the algorithm is stochastic As explained in \link{Sugar and James}, to have a simple low-dimensional representation of the individual curves and to reduce the number of parameters to be estimated, h value must be equals or lower than \eqn{min(p,G-1)}.
}
\details{
Connector provides two different plots to properly guide the choice of the number of clusters:
\itemize{
\item{Elbow Plot:}{ a plot in which the total tightness against the number of clusters is plotted. A proper number of clusters can be inferred as large enough to let the total tightness drop down to relatively little values but as the smallest over which the total tightness does not decrease substantially (we look for the location of an "elbow" in the plot). }
\item{Box Plot:}{ a plot in which for each number of clusters , G, the functional Davies-Buldin (fDB), the cluster separation measure index, is plotted as a boxplot. A proper number of clusters can be associated to the minimum fDB value. }
}

The proximities measures choosen is defined as follow
\deqn{D_q(f,g) = \sqrt( \integral | f^{(q)}(s)-g^{(q)}(s) |^2 ds ), d=0,1,2}
where f and g are two curves and f^{(q)} and g^{(q)} are their q-th derivatives. Note that for q=0, the equation becomes the distance induced by the classical L^2-norm.
Hence, we can define the following indexes for obtaining a cluster separation measure:
\itemize{
\item{T:}{ the total tightness representing the dispersion measure given by
\deqn{ T = \sum_{k=1}^G \sum_{i=1}^n D_0(\hat{g}_i, \bar{g}^k)}; }
\item{S_k:}{
\deqn{ S_k = \sqrt{\frac{1}{G_k} \sum_{i=1}^{G_k} D_q^2(\hat{g}_i, \bar{g}^k);} }
with G_k the number of curves in the k-th cluster;
}
\item{M_{hk}:}{ the distance between centroids (mean-curves) of h-th and k-th cluster
\deqn{M_{hk} =  D_q(\bar{g}^h, \bar{g}^k);}
}
\item{R_{hk}:}{  a measure of how good the clustering is,
\deqn{R_{hk} =   \frac{S_h + S_k}{M_{hk}};}
}
\item{fDB_q:}{ functional Davies-Bouldin index, the cluster separation measure
\deqn{fDB_q = \frac{1}{G} \sum_{k=1}^G \max_{h \neq k} {  \frac{S_h + S_k}{M_{hk}} } }
}
}
}
\seealso{
MostProbableClustering.Extrapolation, BoxPlot.Extrapolation, ConsMatrix.Extrapolation.
}
\author{
Cordero Francesca, Pernice Simone, Sirovich Roberta
}
