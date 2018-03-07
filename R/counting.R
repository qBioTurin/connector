#' Counting samples
#'
#' @description
#' In each cluster it is counted the number of samples belonging to the feature chosen by the user.
#'
#' @param clusterdata Output of the ClusterWithMeanCurve or FittingAndClustering function.
#' @param Model Vector of the model names, i.e. FCM, Malthus, Logistic or Gompertz. If NULL then all the models will be considered.
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
CountingSamples<-function(clusterdata,Model=NULL)
{

 if(is.null(Model))
 {
   Model<-c("FCM","Malthus","Gompertz","Logistic")
 }

 Counting<-list()

 for( i in Model)
 {
   # Possibility to consider the ClusterWithMeanCurve output

   if(is.null(clusterdata[[i]])) {
     if(is.null(cluster$Information))
     {
       warning("File in input is different from the ClusterWithMeanCurve or FittingAndClustering output")

       break
     }
     ClustCurve<-clusterdata$Information$ClustCurve
   }else{ ClustCurve<-clusterdata[[i]]$Information$ClustCurve}


   feature<-tail(colnames(ClustCurve),1)
   a<-count(ClustCurve, c("ID", "Cluster",feature))[,-4]
   Counting[[paste(i)]]<-count(a,c( "Cluster",feature))
 }

 return(Counting)
}
