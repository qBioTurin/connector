library(connector)

### Data files
GrowDataFile<-system.file("data", "475dataset.xlsx", package = "connector")
AnnotationFile <-system.file("data", "475info.txt", package = "connector")

### Merge curves and target file
CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)

### Visualization

GrowPlot<-GrowthCurve(CONNECTORList,"Progeny")
Timegrid <- TimeGridDensity(CONNECTORList)

Datavisual<-DataVisualization(CONNECTORList,
                              feature="Progeny",
                              labels = c("Time","Volume","Tumor Growth"))

Datavisual


### Truncation
CONNECTORList.trunc<- DataTruncation(CONNECTORList,feature="Progeny",70,
                                     labels = c("Time","Volume","Tumor Growth"))

### Calculation of p
CrossLogLike<-BasisDimension.Choice(CONNECTORList.trunc,2:6,Cores = 2)
save(CrossLogLike,file = "CrossLL.Rds")

CrossLogLike$CrossLogLikePlot
CrossLogLike$KnotsPlot

# p is 
p<-3

#### New part:
S.cl <-ClusterAnalysis(CONNECTORList.trunc,G=2:5,
                       p=p,
                       runs=100,
                       Cores=2)

IndexesPlot.Extrapolation(S.cl)-> indexes
indexes$Plot 

ConsMatrix.Extrapolation(S.cl)-> ConsInfo
MostProbableClustering.Extrapolation(S.cl,4) ->MostProbableClustering

ConsInfo$G4$ConsensusPlot
ConsInfo$G5$ConsensusPlot


FMplots<- ClusterWithMeanCurve(clusterdata = MostProbableClustering,
                               feature = "Progeny",
                               labels = c("Time","Volume"),
                               title= (" FCM model"))


# Spline.plots(FCMplots,All=T,path="~/Desktop/")

### Disriminant Plot (goodness of the cluster) just for h = 1 or 2
DiscriminantPlot(clusterdata = MostProbableClustering,
                 feature = "Progeny")

### Counting samples distribution into the clusters

NumberSamples<-CountingSamples(clusterdata=MostProbableClustering,
                               feature = "Progeny")


######
# Advanced Analysis
#######
### Plotting discriminant functions
######
MaxDiscrPlots<-MaximumDiscriminationFunction(clusterdata = MostProbableClustering)

MaxDiscrPlots
