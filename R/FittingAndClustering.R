#' Fitting and Clustering
#'
#'@description
#' Fitting and clustering the data with respect to all FCM, Malthus, Gompertz and Logistic models.
#'
#' @param data CONNECTORList.
#' @param clusterdata Object, or a list of objects belonging to the class funcyOutList.
#' @param k Number of clusters.
#' @param h Dimension of the cluster mean space.
#' @param feature  The column name reported in the AnnotationFile containing the feature interesting for the user to be investigated.
#' @param save If TRUE the plots of the growth curves are saved in different pdf files, one per model.
#' @param path Path to save plot to (combined with file name). If it is missing, the plot is saved in the working  directory.
#' @param labels  Vector containing the text for the title of axis.
#' @return  Returns a list with five members. In details, a list for each model among FCM, Malthus, Gompertz and Logistic model that stores (i) the cluster mean curves plot, (ii) the growth curves plots partitioned according to cluster membership, (iii) all the information related to the clustering, i.e. the mean curve values, the samples cluster membership and a dataframe reporting the for each sample the corresponding ID, time and observation values, feature and belonging cluster. The last component is the cluster mean curves plot for all four models.
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
#' ### Using directly the list corresponding to the values h=2 and k=4 
#'
#' CONNECTORList.FCM.k4.h2<- CONNECTORList.FCM$FCM_all$`k= 4`$`h=2`
#' CONNECTORList.models <- FittingAndClustering( data= CONNECTORList, clusterdata = CONNECTORList.FCM.k4.h2, feature = "Progeny", labels = c("time","volume"))
#'
#' @import ggplot2 cowplot
#' @export
FittingAndClustering<-function(data,clusterdata,h=NULL,k=NULL,feature,save=FALSE,path=NULL,labels=NULL)
{
  out<-list()

  if( !is.null(clusterdata$FCM_all) )
  {
    if(!is.null(h) & !is.null(k) )
    {
      FCM_all<-clusterdata$FCM_all
      out.funcit <-FCM_all[[paste("k=",k)]][[paste("h=",h)]]
    }else{
      warning("Insert a value for h or/and k.")
    }
  }else{
    out.funcit<-clusterdata
  }

  models<-c("FCM","Malthus","Gompertz","Logistic")
  for(i in models)
  {
    out[[paste(i)]]<-ClusterWithMeanCurve(out.funcit,data,k,i,feature = feature,labels=labels)

  }

  mcurves<-list(out$FCM$plotMeanCurve,out$Malthus$plotMeanCurve,out$Gompertz$plotMeanCurve,out$Logistic$plotMeanCurve)
  out$MeanCurves<-plot_grid(plotlist = mcurves)
  if(save==TRUE)
  {
     if(is.null(path)) path <- getwd()
     ggsave(filename="MeanCurves.pdf",plot =out$MeanCurves,width=29, height = 20, units = "cm",scale = 1,path = path )
     ggsave(filename = "FCMCluster.pdf",plot=out$FCM$plotsCluster$ALL,width=29, height = 20, units = "cm",scale = 1,path = path)
     ggsave(filename = "MalthusCluster.pdf",plot=out$Malthus$plotsCluster$ALL,width=29, height = 20, units = "cm",scale = 1,path = path)
     ggsave(filename = "GompertzCluster.pdf",plot=out$Gompertz$plotsCluster$ALL,width=29, height = 20, units = "cm",scale = 1,path = path)
     ggsave(filename = "LogisticCluster.pdf",plot=out$Logistic$plotsCluster$ALL,width=29, height = 20, units = "cm",scale = 1,path = path)
  }

  return(out)
}
