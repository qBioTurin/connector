

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#installation" id="toc-installation">Installation</a>
    -   <a href="#quick-start" id="toc-quick-start">Quick Start</a>
    -   <a href="#docker" id="toc-docker">Docker</a>
-   <a href="#case-study-on-tumor-growth-dataset."
    id="toc-case-study-on-tumor-growth-dataset.">Case study on tumor growth
    dataset.</a>
    -   <a href="#data-importing-by-files" id="toc-data-importing-by-files">Data
        Importing by files</a>
    -   <a href="#data-importing-by-data-frame"
        id="toc-data-importing-by-data-frame">Data importing by data frame</a>
    -   <a href="#data-visualization" id="toc-data-visualization">Data
        Visualization</a>
    -   <a href="#model-selection-tools" id="toc-model-selection-tools">Model
        Selection Tools</a>
        -   <a href="#the-spline-basis-dimension---p"
            id="toc-the-spline-basis-dimension---p">The spline basis dimension -
            p</a>
        -   <a href="#the-number-of-clusters---g"
            id="toc-the-number-of-clusters---g">The number of clusters - G</a>
    -   <a href="#output" id="toc-output">Output</a>
        -   <a href="#connector-clusters" id="toc-connector-clusters">CONNECTOR
            clusters</a>
        -   <a href="#discrimination-plot"
            id="toc-discrimination-plot">Discrimination Plot</a>
        -   <a href="#discrimination-function"
            id="toc-discrimination-function">Discrimination Function</a>
        -   <a href="#estimated-curve" id="toc-estimated-curve">Estimated Curve</a>
-   <a href="#references" id="toc-references">References</a>

# Introduction

The transition from the evaluation of a single time point to the
examination of the entire dynamic evolution of a system is possible only
in the presence of the proper framework. Here we introduce CONNECTOR, a
data-driven framework able to analyze and inspect longitudinal
high-dimensional data in a straightforward and revealing way. CONNECTOR
is a tool for the unsupervised analysis of longitudinal data, that is,
it can process any sample consisting of measurements collected
sequentially over time. CONNECTOR is built on the model-based approach
for clustering functional data presented in (James and Sugar 2003),
which is particularly effective when observations are sparse and
irregularly spaced, as growth curves usually are.

The CONNECTOR framework is based on the following steps:

**(1)** The **pre-processing step** consists of a visualization of the
longitudinal data by a line plot and a time grid heat map helping in the
inspection of the sparsity of the time points.

**(2)** Then, by means of the FDA analysis, the sampled curves are
processed by CONNECTOR with a functional clustering algorithm based on a
mixed effect model. This step requires a **model selection phase**, in
which the dimension of the spline basis vector measures and number of
clusters are computed that help the user to properly set the two free
parameters of the model. The selection of the optimal set of parameters
is supported by the generation of several plots.

**(3)** Once the model selection phase is completed, the **output** of
CONNECTOR is composed of several graphical visualizations to easily mine
the results.

<!-- For a more detailed review on tools and methods developed by the authors on parameter settings, please see [our NM paper].
The main idea is to model individuals curves $g_i$ with basis functions
\begin{equation}\label{eq:gspline}
g_i(t) = \mathbf{s}(t)^T  \boldsymbol{\eta}_i,
\end{equation}
where $\mathbf{s}(t)$ is a $p-$dimensional spline basis vector and $\boldsymbol{\eta}_i$ is a vector of spline coefficients. The $\boldsymbol{\eta}_i$'s are treated with a random-effects model rather than considering them as parameters and fitting a separate spline curve for each individual. The model is fitted through an EM procedure.
Particular attention will be devoted to model selection. In particular, one must choose the number of clusters to fit, and the dimension of the spline basis. `connector` provides a complete toolkit to lead the user through such decisions, one of the most difficult problems in cluster analysis.  -->

# Installation

The `connector` package is hosted on the GitHub platform. The simplest
way to obtain `connector` is to install it using `devtools`. Type the
following commands in R console:

    # Install
    install.packages("devtools",repos = "http://cran.us.r-project.org")
    library(devtools)
    install_github("qBioTurin/connector", ref="master",dependencies=TRUE)
    # Load
    library(connector)

Users may change the `repos` options depending on their locations and
preferences. For more details, see `help(install.packages)`.

The following packages are required before installing *CONNECTOR*:

    install.packages(c('cowplot', 'fda', 'flexclust','ggplot2',
                       'MASS', 'Matrix', 'plyr','ggplotify',
                       'RColorBrewer', 'readxl',
                       'reshape2', 'splines', 'statmod',
                       'sfsmisc','shinyWidgets', 'viridis',
                       'dashboardthemes','shinybusy',
                       'shinydashboard','shinyjs','tidyr'))

## Quick Start

A demo is available to provide the list of commands to perform the
analysis of the example treated in this tutorial. To run the example
just type:

    demo("MainCommandsList", package = "connector")

## Docker

To use Docker it is necessary install Docker on your machine. For more
info see this document
[https://docs.docker.com/engine/installation/](Docker%20installation).
Ensure your user has the rights to run docker (no sudo). To create the
docker group and add your user:

Create the docker group.

      $ sudo groupadd docker

Add your user to the docker group.

      $ sudo usermod -aG docker $USER

Log out and log back in, your group membership is re-evaluated.

To run CONNECTOR from docker Permalink without installing CONNECTOR it
is necessary download the image from dockerhub

      $ docker pull qbioturin/connector:latest

Run the RStudio Docker image with CONNECTOR installed and work with
RStudio in the browser

      $ docker run --cidfile=dockerID -e DISABLE_AUTH=true -p 8787:8787 -d qbioturin/connector:latest
      $ open http://localhost:8787/

From CONNECTOR, it is possible to download and run the RStudio Docker
image by using the CONNECTOR functions called *downloadContainers* and
*docker.run* as follows:

     $ downloadContainers()
     $ docker.run()

# Case study on tumor growth dataset.

To illustrate how `connector` works we use patient-derived xenografts
(PDXs) lines from 21 models generated from chemotherapy-naive high-grade
serous epithelial ovarian cancer.

## Data Importing by files

The analysis starts from two distinct files.

1.  An excel file reporting the discretely sampled curve data. The
    longitudinal data must be reported in terms of time *t* and y-values
    *y*. Each sample is described by two columns: the first column,
    named *time*, includes the lags; while the second column, named with
    the *ID sample*, contains the list of *y* values. Note that, each
    sample must have different *ID sample*. Hence, the 21 tumor growth
    curves are collected in a file composed of 48 columns. See Figure .

2.  A csv (or txt) file contains the annotations associated with the
    sampled curves. The first column reports the **IDSample**, the other
    columns report the features relevant for the analysis, one per
    column. Note that the column *ID sample* must contain the same ID
    names which appear in the excel file of the sampled curves. See
    Figure .

<figure>
<img src="./Fig/excel_file.png" style="width:90.0%"
alt="First lines of the excel file of the sampled curve data." />
<figcaption aria-hidden="true">First lines of the excel file of the
sampled curve data.</figcaption>
</figure>

<figure>
<img src="./Fig/txt_file.png" style="width:70.0%"
alt="First lines of the csv file of the annotated features." />
<figcaption aria-hidden="true">First lines of the csv file of the
annotated features.</figcaption>
</figure>

The two files are imported by the `DataImport` function, which arguments
are the file names. In this example `TimeSeriesFile` and
`AnnotationFile`.

    # find the full path names of the example files
    TimeSeriesFile<-system.file("data", "475dataset.xlsx", package = "connector")
    AnnotationFile <-system.file("data", "475info.txt", package = "connector")
    # import the samples
    CONNECTORList<-DataImport(TimeSeriesFile = TimeSeriesFile,
                              AnnotationFile = AnnotationFile)
    ## ############################### 
    ## ######## Summary ##############
    ## 
    ##  Number of curves: 21 ;
    ##  Min curve length:  7 ; Max curve length:  18 .
    ## ###############################

A list of four objects is created:

    # show the CONNECTORList structure
    str(CONNECTORList)
    ## List of 4
    ##  $ Dataset :'data.frame':    265 obs. of  3 variables:
    ##   ..$ ID         : int [1:265] 1 1 1 1 1 1 1 1 1 1 ...
    ##   ..$ Observation: num [1:265] 63.5 0 78.1 137.9 188.4 ...
    ##   ..$ Time       : num [1:265] 9 16 23 31 36 44 51 57 64 72 ...
    ##  $ LenCurv : int [1:21] 15 15 15 14 14 7 7 11 7 11 ...
    ##  $ LabCurv :'data.frame':    21 obs. of  5 variables:
    ##   ..$ IDSample    : chr [1:21] "475_P1b" "475_P1bP2a" "475_P1bP2b" "475_P1bP2aP3a" ...
    ##   ..$ Progeny     : chr [1:21] "  P1" "  P2" "  P2" "  P3" ...
    ##   ..$ Source      : chr [1:21] "  P0" "  475_P1b" "  475_P1b" "  475_P1bP2a" ...
    ##   ..$ Real.Progeny: chr [1:21] "  P3" "  P4" "  P4" "  P5" ...
    ##   ..$ ID          : int [1:21] 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ TimeGrid: num [1:66] 3 5 6 8 9 10 12 13 14 15 ...

The components of the `CONNECTORList` are:

1.  `$Dataset`, a data frame with three variables: `ID` of the curve,
    `Observation` the *y* values and `Time` the time lags;

2.  `$LenCurv`, the vector reporting the number of observations per
    sample;

3.  `$LabCurv`, a data frame matching the sample with the corresponding
    annotated features. The variables are extracted and named from the
    `AnnotationFile`;

4.  `$TimeGrid`, the vector containing the time grid points.

## Data importing by data frame

The data can be impoted by two data frames used to generate
`CONNECTORList`:

1.  *TimeSeriesFrame*: dataframe of three columns storing the time
    series. The first columns, called “”, stores the sample
    identification. The second column, labeled “”, contains the
    observations over the time, and the third column, called “, reports
    the time point.

2.  *AnnotationFrame*: dataframe with a number of rows equal to the
    number of samples in the *TimeSeriesFrame*. The first column must be
    called “” reporting the sample identifiers. The other columns
    correspond to the features associate with the respective sample (one
    column per sample). If the dataframe is composed of one column, in
    the `CONNECTORList` only the feature will be reported.

Hence, `CONNECTORList` can be generated by exploiting the
*DataFrameImport* function.

    # show the CONNECTORList structure
    CONNECTORList <- DataFrameImport(GrowDataFrame = GrowDataFrame,
                                     AnnotationFrame = AnnotationFrame)

## Data Visualization

The `PlotTimeSeries` function generates the plot of the sampled curves,
coloured by the selected feature from the `AnnotationFile`. Figure
reports the line plot generated by the following code:

    CurvesPlot<-PlotTimeSeries(data = CONNECTORList,
                               feature = "Progeny")
    CurvesPlot

![Sampled curves, coloured by progeny feature.
](/private/var/folders/81/9j7tyzks7c34104pkjdybty40000gn/T/RtmpyJx8w5/preview-1559d662b3ac0.dir/CONNECTORguide_files/figure-markdown_strict/unnamed-chunk-14-1.png)

As we can see from , only few samples have observations at the last time
points. The `DataVisualization` function plots the time grid heat map
helping in the inspection of the sparsity of the time points. The
`DataTruncation` function have been developed to truncate the time
series at specific time value.

    # Growth curves and time grid visualization 
    Datavisual<-DataVisualization(data = CONNECTORList,
                                  feature = "Progeny", 
                                  labels = c("Time","Volume","Tumor Growth"))
    Datavisual

<img src="/private/var/folders/81/9j7tyzks7c34104pkjdybty40000gn/T/RtmpyJx8w5/preview-1559d662b3ac0.dir/CONNECTORguide_files/figure-markdown_strict/unnamed-chunk-15-1.png" alt="Sampled curves, coloured by progeny feature (left panel) and time grid density (right panel). \label{fig:gridDensity}"  />
<p class="caption">
Sampled curves, coloured by progeny feature (left panel) and time grid
density (right panel).
</p>

In Figure the time grid density is showed. Each point
*p*<sub>*x*, *y*</sub> is defined by a pair of coordinates
*p*<sub>*x*, *y*</sub> = (*x*,*y*)  and by a colour.
*p*<sub>*x*, *y*</sub> is defined if only if exists at least one sample
with two observations at time *x*  and *y*. The colour encodes the
number of samples in which *p*<sub>*x*, *y*</sub> is reported. In our
example, according to Figure we decide to truncate the observations at
70 days.

    # data truncation
    trCONNECTORList<-DataTruncation(data = CONNECTORList,
                                    feature="Progeny",
                                    truncTime = 70,
                                    labels = c("Time","Volume","Tumor Growth"))
    ## ############################################################## 
    ## ######## Summary of the trunc. data ############
    ## 
    ##  Number of curves: 21 ;
    ##  Min curve length:  7 ; Max curve length:  16 ;
    ## 
    ##  Number of truncated curves: 6 ;
    ##  Min points deleted:  2 ; Max points deleted:  6 ;
    ## ##############################################################
    # the trCONNECTORList structure
    str(trCONNECTORList, max.level = 1)
    ## List of 6
    ##  $ Dataset            :'data.frame': 241 obs. of  3 variables:
    ##  $ LenCurv            : num [1:21] 9 9 9 14 14 7 7 9 7 9 ...
    ##  $ LabCurv            :'data.frame': 21 obs. of  5 variables:
    ##  $ TimeGrid           : num [1:50] 3 5 6 8 9 10 12 13 14 15 ...
    ##  $ ColFeature         : chr [1:5] "#FF0000" "#CCFF00" "#00FF66" "#0066FF" ...
    ##  $ PlotTimeSeries_plot:List of 9
    ##   ..- attr(*, "class")= chr [1:2] "gg" "ggplot"

The output of the function `DataTruncation` is a list of six objects
reporting th truncated versions of the `DataImport` function. Moreover,
the plot stored in `$PlotTimeSeries_plot` shows the sampled curves plot
with a vertical line indicating the truncation time, see Figure .

    # plot
    trCONNECTORList$PlotTimeSeries_plot

![Sampled curves, coloured by progeny feature and truncation lag
(vertical solid black line).
](/private/var/folders/81/9j7tyzks7c34104pkjdybty40000gn/T/RtmpyJx8w5/preview-1559d662b3ac0.dir/CONNECTORguide_files/figure-markdown_strict/unnamed-chunk-18-1.png)

## Model Selection Tools

Before running the fitting and clustering approach, we have to properly
chose the two free parameters:

1.  the spline basis dimension, *p*;

2.  the number of clusters, *G*.

We developed several functions to enable the user to properly set the
free parameters.

### The spline basis dimension - p

The dimension of the spline basis can be chose by exploiting the
`BasisDimension.Choice` function, by taking the
*p* ∈ \[*p*<sub>*m**i**n*</sub>,*p*<sub>*m**a**x*</sub>\] value
corresponding to the largest cross-validated likelihood, as proposed in
(James, Hastie, and Sugar 2000), where the interval of values
\[*p*<sub>*m**i**n*</sub>,*p*<sub>*m**a**x*</sub>\] is given by the
user. In particular, a ten-fold cross-validation approach is
implemented: firstly the data are split into 10 equal-sized parts,
secondly the model is fitted considering 9 parts and the computation of
the log-likelihood on the excluded part is performed.

In details two plots are returned:

1.  The `$CrossLogLikePlot` return the plot of the mean tested
    log-likelihoods versus the dimension of the basis, see Figure . Each
    gray dashed line corresponds to the cross-log-likelihood values
    obtained on different test/learning sets and the solid black line is
    their mean value. The resulting plot could be used as a guide to
    choose the largest cross-validated likelihood. Specifically, the
    optimal value of `p` is generally the smallest ensuring the larger
    values of the mean cross-log-likelihood function.

2.  The `$KnotsPlot` shows the time series curves (top) together with
    the knots distribution over the time grid (bottom). The knots divide
    the time domain into contiguous intervals, and the curves are fitted
    with separate polynomials in each interval. Hence, this plot allows
    the user to visualize whether the knots properly split the time
    domain considering the distribution of the sampled observations. The
    user should choose *p* such that the number of cubic polynomials
    (i.e., *p* − 1) defined in each interval is able to follow the
    curves dynamics.

<!-- -->

    # ten-fold cross-validation 
    CrossLogLike<-BasisDimension.Choice(data = trCONNECTORList,
                                        p = 2:6 )
    CrossLogLike$CrossLogLikePlot

<img src="/private/var/folders/81/9j7tyzks7c34104pkjdybty40000gn/T/RtmpyJx8w5/preview-1559d662b3ac0.dir/CONNECTORguide_files/figure-markdown_strict/unnamed-chunk-19-1.png" alt="Cross-validated loglikelihood functions. \label{fig:crossloglike}"  />
<p class="caption">
Cross-validated loglikelihood functions.
</p>

    CrossLogLike$KnotsPlot

<img src="/private/var/folders/81/9j7tyzks7c34104pkjdybty40000gn/T/RtmpyJx8w5/preview-1559d662b3ac0.dir/CONNECTORguide_files/figure-markdown_strict/unnamed-chunk-20-1.png" alt="Knots ditribution. \label{fig:Knotscrossloglike}" width="60%" />
<p class="caption">
Knots ditribution.
</p>

    # set p
    p <- 3

In our example, the optimal value of *p* is 3.

### The number of clusters - G

`connector` provides two different plots to properly guide to set the
number of clusters. Two measures of proximity are introduced: the *total
tightness* *T* and the *functional Davied-Bouldin index* *f**D**B*. Both
measures rely on the family of semi-metrics between curves defined as
where *f* and *g* are two curves and *f*<sup>(*q*)</sup> and
*g*<sup>(*q*)</sup> are their *q*th derivatives. Note that for *q* = 0,
eq. is the distance induced by the classical *L*<sup>2</sup>−norm. It
turns out that *D*<sub>*q*</sub> can be reliably calculated in our
setting, since we are interested in proximity measures between curves
and center-curves for each cluster (tightness of the cluster), as well
as center-curve and center-curve of different clusters (separateness of
clusters). *D*<sub>*q*</sub> is calculated taking advantage of the
spline representation of the estimated curves and mean curves, see eq. .

The *total tightness* *T* is the dispersion measure defined as where
*ĝ*<sub>*i*</sub> is the estimated *i*–th curve given in eq. and
*ḡ*<sup>*k*</sup> is the center of *k*–th cluster given in eq. .

As the number of clusters increases, the total tightness decreases to
zero, the value which is attained when the number of fitted clusters is
equal to the number of sampled curves. In this case, any *k*th cluster
mean curve coincides with an estimated curve and
*D*<sub>0</sub>(*ĝ*<sub>*i*</sub>,*ḡ*<sup>*k*</sup>) = 0 for any *i* and
*k*.

A proper number of clusters can be inferred as large enough to let the
total tightness drop down to little values over which it does not
decrease substantially. Hence, we look for the location of an *elbow* in
the plot of the total tightness against the number of clusters. The
second index, which is a cluster separation measure, is called
*functional David Bouldin* (fDB) index. Specifically, we defined it as
follows where, for each cluster *k*<sub>1</sub> and *k*<sub>2</sub> with
*G*<sub>*k*</sub> the number of curves in the *k*th cluster. The
significance of eq. can be understood as the average of the blend
measures of each cluster from its most overlapping cluster. The optimal
choice of clusters, then, will be that which minimizes this average
blend. Furthermore, the random initialization of the k-means algorithm
to get initial cluster memberships into the FCM algorithm and the
stochasticity characterizing the method lead to some variability among
different runs. For these reasons, multiple runs are necessary to
identify the most frequent clustering fixed a number of clusters (*G*).

To effectively take advantage of those two measures, `connector`
supplies the function `ClusterAnalysis` which repeats the clustering
procedure a number of times equal to the parameter `runs` and for each
of the number of clusters given in the parameter `G`. The output of the
function is a list of three objects (i) `$Clusters.List`: the list of
all the clustering divisions obtained varying among the input *G*
values, (ii) `$seed`: the seed sets before running the method, and (iii)
`$runs`: the number of runs.

Specifically, the object storing the clustering obtained with *G* = 2
(i.e., `ClusteringList$Clusters.List$G2`) is a list of three elements:

1.  `$ClusterAll`: the list of all the possible FCM parameters values
    (see Sec. *Details on the functional clustering model*) that can be
    obtained through each run. In details, the item `$ParamConfig.Freq`
    reports the number of times that the respectively parameters
    configuration (stored in `$FCM`) is found. Let us note, that the sum
    might be not equal to the number of runs since errors may occur;

2.  `$ErrorConfigurationFit`: list of errors that could be obtained by
    running the FCM method with extreme parameters configurations;

3.  `$h.selected`: the dimension of the mean space, denoted as *h*,
    selected to perform the analysis. In details, this parameter gives a
    further parametrization of the mean curves allowing a
    lower-dimensional representation of the curves with means in a
    restricted subspace. Its value is choose such that the inequality
    *h* ≤ min (*G*−1, *p*) holds. Better clusterings are obtained with
    higher values of *h*, but this may lead to many failing runs. In
    this cases *h* values is decreased until the number of successful
    runs are higher than a specific constraint which can be defined by
    the user (by default is 100% of successful runs).

<!-- -->

    ClusteringList <-ClusterAnalysis(data = trCONNECTORList,
                                     G = 2:5,
                                     p = p,
                                     runs = 100)

    # the output structure
    str(ClusteringList, max.level = 2, vec.len=1)
    ## List of 4
    ##  $ Clusters.List:List of 4
    ##   ..$ G2:List of 2
    ##   ..$ G3:List of 2
    ##   ..$ G4:List of 2
    ##   ..$ G5:List of 2
    ##  $ CONNECTORList:List of 6
    ##   ..$ Dataset            :'data.frame':  241 obs. of  3 variables:
    ##   ..$ LenCurv            : num [1:21] 9 9 ...
    ##   ..$ LabCurv            :'data.frame':  21 obs. of  5 variables:
    ##   ..$ TimeGrid           : num [1:50] 3 5 ...
    ##   ..$ ColFeature         : chr [1:5] "#FF0000" ...
    ##   ..$ PlotTimeSeries_plot:List of 9
    ##   .. ..- attr(*, "class")= chr [1:2] "gg" ...
    ##  $ seed         : num 2404
    ##  $ runs         : num 100
    str(ClusteringList$Clusters.List$G2, max.level = 3, vec.len=1)
    ## List of 2
    ##  $ ClusterAll:List of 2
    ##   ..$ :List of 2
    ##   .. ..$ FCM             :List of 4
    ##   .. ..$ ParamConfig.Freq: int 43
    ##   ..$ :List of 2
    ##   .. ..$ FCM             :List of 4
    ##   .. ..$ ParamConfig.Freq: int 57
    ##  $ h.selected: num 1

The function `IndexesPlot.Extrapolation` generated two plots reported in
Figure . In details, the tightness and fDB indexes are plotted for each
value of *G*. Each value of *G* is associated with a violin plot created
by the distribution of the index values collected from the runs
performed. The blue line shows the most probable configuration which
will be exploited hereafter.

    IndexesPlot.Extrapolation(ClusteringList)-> indexes
    indexes$Plot

<img src="/private/var/folders/81/9j7tyzks7c34104pkjdybty40000gn/T/RtmpyJx8w5/preview-1559d662b3ac0.dir/CONNECTORguide_files/figure-markdown_strict/unnamed-chunk-24-1.png" alt="\label{fig:TandfDB} Violin Plots of the {\it total tightness} $T$ calculated on each run and for different number of clusters $G$ (right panel). Violin Plots of the {\it functional DB index} fDB calculated on each run and for different number of clusters $G$ (left panel)."  />
<p class="caption">
Violin Plots of the {} *T* calculated on each run and for different
number of clusters *G* (right panel). Violin Plots of the {} fDB
calculated on each run and for different number of clusters *G* (left
panel).
</p>

<!-- The indexes show that $G=3$ and $G=4$ may be good choices for the parameter. Specifically, the right panel shows that $G=3$ may be good choice for the parameter, while the fDB indexes plotted in the left panel lead the choice to $G=4$, a number of clusters that explicitly minimizes the fDB index. -->

*G* = 4 corresponds to the optimal choice for our example, since it
represents the elbow for the tightness indexes (right panel) and
explicitly minimizes the fDB indexes (left panel).

The variability of the two measures among runs, exhibited in Figure , is
related to the random initialization of the k-means algorithm to get
initial cluster memberships from points. The stability of the clustering
procedure can be visualized through the stability matrix extrapolated by
the function `ConsMatrix.Extrapolation`, as shown in Figure . The plot
informs about the stability of the final clustering across different
runs. Indeed, each cell of the matrix is coloured proportionally to the
frequency of the two corresponding curves belonging to the same cluster
across different runs. Hence the larger the frequencies are (which
corresponds to warmer colours of the cells in the plot), the more stable
is the final clustering. For each cluster the average stability value is
reported.

    ConsMatrix<-ConsMatrix.Extrapolation(stability.list = ClusteringList)
    str(ConsMatrix, max.level = 2, vec.len=1)
    ## List of 4
    ##  $ G2:List of 2
    ##   ..$ ConsensusMatrix: num [1:21, 1:21] 1 1 ...
    ##   .. ..- attr(*, "dimnames")=List of 2
    ##   ..$ ConsensusPlot  :List of 9
    ##   .. ..- attr(*, "class")= chr [1:2] "gg" ...
    ##  $ G3:List of 2
    ##   ..$ ConsensusMatrix: num [1:21, 1:21] 1 1 ...
    ##   .. ..- attr(*, "dimnames")=List of 2
    ##   ..$ ConsensusPlot  :List of 9
    ##   .. ..- attr(*, "class")= chr [1:2] "gg" ...
    ##  $ G4:List of 2
    ##   ..$ ConsensusMatrix: num [1:21, 1:21] 1 1 ...
    ##   .. ..- attr(*, "dimnames")=List of 2
    ##   ..$ ConsensusPlot  :List of 9
    ##   .. ..- attr(*, "class")= chr [1:2] "gg" ...
    ##  $ G5:List of 2
    ##   ..$ ConsensusMatrix: num [1:21, 1:21] 1 0.81 ...
    ##   .. ..- attr(*, "dimnames")=List of 2
    ##   ..$ ConsensusPlot  :List of 9
    ##   .. ..- attr(*, "class")= chr [1:2] "gg" ...
    ConsMatrix$G4$ConsensusPlot

<img src="/private/var/folders/81/9j7tyzks7c34104pkjdybty40000gn/T/RtmpyJx8w5/preview-1559d662b3ac0.dir/CONNECTORguide_files/figure-markdown_strict/unnamed-chunk-25-1.png" alt="\label{fig:ConsMatg4} Stability Matrix for G = 4." width="60%" />
<p class="caption">
Stability Matrix for G = 4.
</p>

<!-- Hence, given the higher stability of the cluster membership among the curves, we are choosing $G=4$ for the illustrated example. -->

Once the free parameters are all set, the function
`MostProbableClustering.Extrapolation` can be used to fix the most
probable clustering with given dimension of the spline basis *p*, and
number of clusters *G* and save the result in a dedicated object.

    CONNECTORList.FCM.opt<-MostProbableClustering.Extrapolation(
      stability.list = ClusteringList,
      G = G )

## Output

The output of `connector` consists of the line plot of the samples
grouped by the `connector` cluster, the discriminant plot, the
discriminant function, and the estimation of the entire curve for each
single subject.

### CONNECTOR clusters

`connector` provides the function `ClusterWithMeanCurve` which plots the
sampled curves grouped by cluster membership, together with the mean
curve for each cluster, see Figure . The function prints as well the
values for *S*<sub>*h*</sub>, *M*<sub>*h**k*</sub> and fDB given in
equation .

    FCMplots<- ClusterWithMeanCurve(clusterdata = CONNECTORList.FCM.opt,
                                    feature = "Progeny",
                                    labels = c("Time","Volume"),
                                    title = "FCM model")
    ## 
    ## ######## S indexes ############
    ## 
    ## 
    ## |          |        S|       S_1|      S_2|
    ## |:---------|--------:|---------:|--------:|
    ## |Cluster B | 1702.753| 114.57314| 5.128564|
    ## |Cluster A |  654.418|  29.34886| 1.021555|
    ## |Cluster C | 2791.051| 136.25605| 6.430096|
    ## 
    ## ##############################################################
    ## ############################################################## 
    ## 
    ##         ######## M indexes ############
    ## 
    ## 
    ## |          | Cluster B| Cluster A| Cluster C|
    ## |:---------|---------:|---------:|---------:|
    ## |Cluster B |     0.000|  2851.265|  5731.600|
    ## |Cluster A |  2851.265|     0.000|  8831.041|
    ## |Cluster C |  5731.600|  8831.041|     0.000|
    ## 
    ## ##############################################################
    ## ######## R indexes ############
    ## 
    ## 
    ## |          |         R|       R_1|      R_2|
    ## |:---------|---------:|---------:|--------:|
    ## |Cluster B | 0.8267106| 1.0914505| 1.530825|
    ## |Cluster A | 0.8267106| 1.0914505| 1.530825|
    ## |Cluster C | 0.7840400| 0.8803222| 1.255647|
    ## 
    ## ##############################################################
    ## ######## fDB indexes ############
    ## 
    ## 
    ## |       fDB|    fDB_1|    fDB_2|
    ## |---------:|--------:|--------:|
    ## | 0.8124871| 1.021074| 1.439099|
    ## 
    ## ##############################################################

<img src="/private/var/folders/81/9j7tyzks7c34104pkjdybty40000gn/T/RtmpyJx8w5/preview-1559d662b3ac0.dir/CONNECTORguide_files/figure-markdown_strict/unnamed-chunk-27-1.png" alt="\label{fig:clusters} Sampled curves grouped by cluster membership." width="60%" />
<p class="caption">
Sampled curves grouped by cluster membership.
</p>

Finally, to inspect the composition of the clusters, the function
`CountingSamples` reports the number and the name of samples in each
cluster according to the feature selected by the user.

    NumberSamples<-CountingSamples(clusterdata = CONNECTORList.FCM.opt,
                                   feature = "Progeny")
    str(NumberSamples, max.level = 2)
    ## List of 2
    ##  $ Counting    :'data.frame':    15 obs. of  3 variables:
    ##   ..$ Cluster: chr [1:15] "A" "A" "A" "A" ...
    ##   ..$ Progeny: Factor w/ 5 levels "  P1","  P2",..: 1 2 3 4 5 1 2 3 4 5 ...
    ##   ..$ n      : int [1:15] 1 2 0 0 0 0 0 2 2 2 ...
    ##  $ ClusterNames:'data.frame':    21 obs. of  2 variables:
    ##   ..$ ID     : int [1:21] 1 2 3 4 5 6 7 8 9 10 ...
    ##   ..$ Cluster: chr [1:21] "A" "A" "A" "B" ...

### Discrimination Plot

A detailed visualization of the clusterization of the sample curves can
be obtained by means of the `DiscriminantPlot` function.  
The low-dimensional plots of curve datasets, enabling a visual
assessment of clustering. The curves can be projected the into the lower
dimensional space of the mean space, in this manner the curve can be
plot as points (with coordinates the functional linear discriminant
components). In details, when *h* is equal to 1 or 2, than the
`DiscriminantPlot` function return the plots of the curves projected
onto the *h* dimensional mean space:

1.  when *h* = 1, the functional linear discriminant *α*\_1 is plotted
    versus its variability;
2.  when *h* = 2, the functional linear discriminant *α*\_1 is plotted
    versus the functional linear discriminant *α*\_2.

If *h* &gt; 2, the principal component analysis is exploited to
extrapolate the components of the *h*-space with more explained
variability (which is reported in brackets), that are then plotted
together.

In the case study here described, we get the plots in Figure which is
colored by cluster membership and in Figure which is colored by the
`"Progeny"` feature.

    DiscrPlt<-DiscriminantPlot(clusterdata = CONNECTORList.FCM.opt,
                               feature = "Progeny")
    DiscrPlt$ColCluster

<img src="/private/var/folders/81/9j7tyzks7c34104pkjdybty40000gn/T/RtmpyJx8w5/preview-1559d662b3ac0.dir/CONNECTORguide_files/figure-markdown_strict/unnamed-chunk-29-1.png" alt="\label{fig:DiscrPlotCL} Curves projected onto the 3-dimensional mean space: symbols are coloured by cluster membership." width="60%" />
<p class="caption">
Curves projected onto the 3-dimensional mean space: symbols are coloured
by cluster membership.
</p>

    DiscrPlt$ColFeature

<img src="/private/var/folders/81/9j7tyzks7c34104pkjdybty40000gn/T/RtmpyJx8w5/preview-1559d662b3ac0.dir/CONNECTORguide_files/figure-markdown_strict/unnamed-chunk-29-2.png" alt="\label{fig:DiscrPlotF} Curves projected onto the 3-dimensional mean space: symbols are coloured by progeny." width="60%" />
<p class="caption">
Curves projected onto the 3-dimensional mean space: symbols are coloured
by progeny.
</p>

### Discrimination Function

There will be *h* discriminant functions and each curve shows the times
with higher discriminatory power, which are the times corresponding to
largest absolute (positive or negative) values on the y-axis, see Figure
.

    MaxDiscrPlots<-MaximumDiscriminationFunction(clusterdata = CONNECTORList.FCM.opt)
    MaxDiscrPlots[[1]]

<img src="/private/var/folders/81/9j7tyzks7c34104pkjdybty40000gn/T/RtmpyJx8w5/preview-1559d662b3ac0.dir/CONNECTORguide_files/figure-markdown_strict/unnamed-chunk-30-1.png" alt="\label{fig:discrimination} Discriminant curve." width="60%" />
<p class="caption">
Discriminant curve.
</p>

Large absolute values correspond to large weights and hence large
discrimination between clusters. From the Figure it is possible assess
that earlier measurements are crucial in determining cluster assignment
as well as latter ones.

### Estimated Curve

The functional clustering procedure predict unobserved portions of the
true curves for each subject. The estimated curves are returned by
`connector` and plotted with confidence intervals as well that is
possible to visualize using the `Spline.plots` function.

    PlotSpline = Spline.plots(FCMplots)
    PlotSpline$`1`

<img src="/private/var/folders/81/9j7tyzks7c34104pkjdybty40000gn/T/RtmpyJx8w5/preview-1559d662b3ac0.dir/CONNECTORguide_files/figure-markdown_strict/unnamed-chunk-31-1.png" alt="\label{fig:Spline} Fitting of the Sample with ID = 1: in blue the observations characterinzing the curve, in red the estimated spline from the fitting, and in black the mean curve of the associated cluster. The grey area represents the confidence interval." width="60%" />
<p class="caption">
Fitting of the Sample with ID = 1: in blue the observations
characterinzing the curve, in red the estimated spline from the fitting,
and in black the mean curve of the associated cluster. The grey area
represents the confidence interval.
</p>

<!-- # Details on the functional clustering model  -->
<!-- The curves, $g_i(t)$ for each $i$th selected individual, are supposed to be observed with measurement errors and only at few discrete time points. Hence the vector $\mathbf{Y}_i$ of observed values at times $t_{i_1}, \dots , t_{i_{n_i}}$ is given as -->
<!-- \begin{equation*} -->
<!--    \mathbf{Y}_i = \mathbf{g}_i + \boldsymbol{\varepsilon}_i, -->
<!-- \end{equation*} -->
<!-- where $\mathbf{g}_i$  and $\boldsymbol{\varepsilon}_i$ are the vectors of true values and measurement errors at time grid, respectively. As there are only finite number of observations, individual curves are modeled using basis functions, in particular cubic splines. Let -->
<!-- \begin{equation}\label{eq:gspline} -->
<!--    g_i(t) = \mathbf{s}(t)^T  \boldsymbol{\eta}_i, -->
<!-- \end{equation} -->
<!-- where $\mathbf{s}(t)$ is a $p-$dimensional spline basis vector and $\boldsymbol{\eta}_i$ is a vector of spline coefficients. The $\boldsymbol{\eta}_i$'s are treated with a random-effects model rather than considering them as parameters and fitting a separate spline curve for each individual. Cluster means are furthermore rewritten as -->
<!-- \begin{equation*} -->
<!--    \boldsymbol{\mu}_{k} = \boldsymbol{\lambda}_0 + \Lambda \boldsymbol{\alpha}_k, -->
<!-- \end{equation*} -->
<!-- where $\boldsymbol{\lambda}_0$ and $\boldsymbol{\alpha}_k$ are $p-$ and $h-$ dimensional vectors, $\Lambda$ is a $(p,h)$ matrix and $h \leq \min(p,G-1)$, where $G$ denote the true number of clusters. This parametrization allows a lower-dimensional representation of the curves with means in a restricted subspace (for $h < G-1$). -->
<!-- With this formulation, the functional clustering model can be written as -->
<!-- \begin{eqnarray} -->
<!--    \mathbf{Y}_i =S_i  \cdot ( \boldsymbol{\lambda}_0 + \Lambda \boldsymbol{\alpha}_{\mathbf{z}_i} +  \boldsymbol{\gamma}_i) +  \boldsymbol{\varepsilon}_i, \quad i=1, \dots, n,\nonumber\\ -->
<!--    \boldsymbol{\varepsilon}_i \sim  \mathcal{N} (\mathbf{0},R), \quad  \boldsymbol{\gamma}_i \sim   \mathcal{N} (\mathbf{0},\Gamma), \qquad \qquad -->
<!-- \end{eqnarray} -->
<!-- where $S_i = (\mathbf{s}(t_{i_1}),\dots,\mathbf{s}(t_{i_{n_i}}))^T$ is the spline basis matrix for the $i-$th curve.  -->
<!-- The model is fitted following [@james2003clustering] and all the estimated parameters and the predicted cluster membership are returned. Notice that the $k$th cluster mean curve can be retrieved as -->
<!-- \begin{equation}\label{eq:splinemeancurve} -->
<!--    \bar{g}^k(t) = \mathbf{s}(t)^T ( \hat{\boldsymbol{\lambda}}_0 + \hat{\Lambda} \hat{ \boldsymbol{\alpha}}_k). -->
<!-- \end{equation} -->
<!-- Moreover, the functional clustering procedure can accurately predict unobserved portions of the curves $g_i(t)$ by means of the natural estimate -->
<!-- \begin{equation}\label{eq:gsplinepred} -->
<!--    \hat{g}_i(t) =  \mathbf{s}(t)^T  \hat{\boldsymbol{\eta}}_i, -->
<!-- \end{equation} -->
<!-- where $\hat{\boldsymbol{\eta}}_i$ is a prediction for $\boldsymbol{\eta}_i$ which is proven to be optimally computed as $\mathbb{E}( \boldsymbol{\eta}_i \;|\; \mathbf{Y}_i)$ and explicitly given in [@james2003clustering], eq. (17).  -->

# References

James, Gareth M, Trevor J Hastie, and Catherine A Sugar. 2000.
“Principal Component Models for Sparse Functional Data.” *Biometrika* 87
(3): 587–602.

James, Gareth M, and Catherine A Sugar. 2003. “Clustering for Sparsely
Sampled Functional Data.” *Journal of the American Statistical
Association* 98 (462): 397–408.
