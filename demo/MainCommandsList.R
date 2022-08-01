library(connector)

### Data files
TimeSeriesFile<-system.file("data", "475treatedDataset.xlsx", package = "connector")
AnnotationFile <-system.file("data", "475treatedInfo.txt", package = "connector")

### Merge curves and target file
CONNECTORList <- DataImport(TimeSeriesFile,AnnotationFile)

### Visualization
# To plot just the longitudinal data
GrowPlot<- PlotTimeSeries(CONNECTORList,"Progeny")
GrowPlot$PlotTimeSeries_plot
# To plot just the Time Grid 
Timegrid <- TimeGridDensity(CONNECTORList)
Timegrid$TimeGrid_plot
# To visualize both the plots together
Datavisual<-DataVisualization(CONNECTORList,
                              feature="Progeny",
                              labels = c("Time","Volume","Tumor Growth"))
Datavisual

### Truncation
CONNECTORList.trunc<- DataTruncation(CONNECTORList,
                                     feature="Progeny",
                                     truncTime = 70,
                                     labels = c("Time","Volume","Tumor Growth"))

### Calculation of p
CrossLogLike<-BasisDimension.Choice(CONNECTORList.trunc,2:6,Cores = 2)

CrossLogLike$CrossLogLikePlot
CrossLogLike$KnotsPlot

# p is 
p<-3

### Cluster Analysis to set G
S.cl <-ClusterAnalysis(CONNECTORList.trunc,
                       G=2:5,
                       p=p,
                       runs=50,
                       Cores=1)

IndexesPlot.Extrapolation(S.cl)-> indexes
indexes$Plot

ConsMatrix.Extrapolation(S.cl)-> ConsInfo
ConsInfo$G3$ConsensusPlot
ConsInfo$G4$ConsensusPlot

MostProbableClustering.Extrapolation(S.cl,4) ->MostProbableClustering

FCMplots<- ClusterWithMeanCurve(clusterdata = MostProbableClustering,
                               feature = "Progeny",
                               labels = c("Time","Volume"),
                               title= ("FCM model"))


PlotSpline = Spline.plots(FCMplots)
PlotSpline$`1`

### Discriminant Plot (goodness of the cluster)
DiscriminantPlot(clusterdata = MostProbableClustering,
                 feature = "Progeny")

### Counting samples distribution into the clusters

NumberSamples<-CountingSamples(clusterdata=MostProbableClustering,
                               feature = "Progeny")
NumberSamples$Counting
NumberSamples$ClusterNames

######
# Advanced Analysis
#######
### Plotting discriminant functions
######
MaxDiscrPlots<-MaximumDiscriminationFunction(clusterdata = MostProbableClustering)

MaxDiscrPlots[[1]]
