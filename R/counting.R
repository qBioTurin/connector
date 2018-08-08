#' Counting samples
#'
#' @description
#' For each cluster it is counted the number of samples distinguished by the feature.
#'
#' @param Clusterdata 
#' @param data
#' @param feature The column name reported in the AnnotationFile containing the feature interesting for the user to be investigated. If NULL all the features will be considered.
#' @return Returns a list storing for each model analyzed a matrix composed by three columns: (i) cluster, (ii) the feature name, and (iii) the number of samples. 
#' @examples
#'
#' ### Data files
#' GrowDataFile<-"data/1864dataset.xls"
#' AnnotationFile <-"data/1864info.txt"
#'
#' ### Merge curves and target file
#' CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#' CONNECTORList<- DataTruncation(CONNECTORList,feature="Progeny",truncTime=60,save=TRUE,path="~/Desktop/ImagesPerFrancesca/",labels = c("time","volume","Tumor Growth"
#'
#' CONNECTORList.FCM <- ClusterChoice(CONNECTORList,K=c(2:6),h=2)
#'
#' ### Using the CONNECTORList.FCM and the values: h=2 and k=4.
#'
#' CONNECTORList.models <- FittingAndClustering( data= CONNECTORList, clusterdata = CONNECTORList.FCM, h = 2, k=4, feature = "Progeny", labels = c("time","volume"))
#'
#' ###  Counting the samples for all the models.
#'
#' CountingSamples(CONNECTORList.models)
#'
#' ###  Counting the samples just for the FCM.
#'
#' CountingSamples(CONNECTORList.models, "FCM")
#'
#' @import plyr
#' @export
#' 
CountingSamples<-function(clusterdata,data,feature="ID")
{
  
  
  if(isS4(clusterdata))
  {
    
    ClustCurve <-clusterdata@models$fitfclust@fit$ClustCurves
    
  }else{
    
    ClustCurve <- clusterdata$Summary
    
  }
  
  ClustCurve <- data.frame(merge(ClustCurve,data$LabCurv,by="ID"))
  
  a<-count(ClustCurve, c("ID", "Cluster",feature))[,-length(count(ClustCurve, c("ID", "Cluster",feature)))]
  Counting<-count(a,c( "Cluster",feature))
  
  return(Counting)
}

