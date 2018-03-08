library(Prova1)

### Data files
GrowDataFile<-"data/1864dataset.xls"
AnnotationFile <-"data/1864info.txt"

### Merge curves and target file
CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)

### Visualization
DataVisualization(CONNECTORList,feature="Progeny",labels = c("time","volume","Tumor Growth"))

### Truncation

CONNECTORList<- DataTruncation(CONNECTORList,feature="Progeny",truncTime=60,labels = c("time","volume","Tumor Growth"))

### Calculation of h

pca <- PCA.Analysis(CONNECTORList$Dataset)
pca$plot

### Calculation of k

CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:6),h=2)
CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:6),PCAperc = pca$perc)

CONNECTORList.FCM.k4.h2<- CONNECTORList.FCM$FCM_all$`k= 4`$`h= 2`

### Fitting and clustering

FCMplots <- ClusterWithMeanCurve(clusterdata = CONNECTORList.FCM.k4.h2,data= CONNECTORList, k = 4, model = "FCM",feature = "Progeny",labels = c("Time","Volume"))

MalthusPlots<- ClusterWithMeanCurve(data = CONNECTORList,k = 4,model="Malthus",feature = "Progeny")

### Fitting and clustering considering all the models

CONNECTORList.models <- FittingAndClustering(data = CONNECTORList, clusterdata = CONNECTORList.FCM, h = 2, k=4, feature = "Progeny", labels = c("time","volume"))

CONNECTORList.models <- FittingAndClustering(data = CONNECTORList, clusterdata = CONNECTORList.FCM.k4.h2, feature = "Progeny", labels = c("time","volume"))

### Withinness and betweenness plot for Malthus and FCM

Malthus.ClustCurve <-  CONNECTORList.models$Malthus$Information$ClustCurve
Malthus.MeanCurves <-  CONNECTORList.models$Malthus$Information$meancurves

PlotWithinnessBetweenness(Malthus.ClustCurve,Malthus.MeanCurves,Title = "Malthus Cluster betweenness and withinness")

FCM.ClustCurve <-  CONNECTORList.models$FCM$Information$ClustCurve
FCM.MeanCurves <-  CONNECTORList.models$FCM$Information$meancurves

PlotWithinnessBetweenness(FCM.ClustCurve,FCM.MeanCurves,Title = "FCM Cluster betweenness and withinness")

### Counting the samples w.r.t. all the models

NumberSamples<-CountingSamples(CONNECTORList.models,feature = "Progeny")
