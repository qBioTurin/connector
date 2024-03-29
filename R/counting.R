#' Counting samples
#'
#' @description
#' Counts for each combination of cluster and feature the number of samples and their name are returned.
#'
#' @param clusterdata The list obtained from extrapolating the most probable clustering from the StabilityAnalysis function output. (see \code{\link{StabilityAnalysis}} and \code{\link{MostProbableClustering.Extrapolation}}). 
#' @param feature The column name reported in the AnnotationFile containing the feature  to be investigated.
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @return CountingSample returns a list containing  1) a matrix composed by three columns: (i) cluster, (ii) the feature name, and (iii) the number of samples, and 2) a dataframe storing the sample names reported in the AnnotationFile and the corresponding  cluster membership.
#' 
#' 
#' @examples
#'
#' ### Data files
#' GrowDataFile<-"data/1864dataset.xls"
#' AnnotationFile <-"data/1864info.txt"
#' 
#' ### Merge curves and target file
#' CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#' 
#'### Truncation
#'
#'CONNECTORList<- DataTruncation(CONNECTORList,feature="Progeny",60,labels = c("time","volume","Tumor Growth"))
#'
#' CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:6),h=2)
#' CONNECTORList.FCM.k4.h2<- CONNECTORList.FCM$FCM_all$`k= 4`$`h= 2`
#'
#' NumberSamples<-CountingSamples(clusterdata=CONNECTORList.FCM.k4.h2,CONNECTORList,feature = "Progeny")
#' @import dplyr 
#' @export
#' 
CountingSamples<-function(clusterdata,feature="ID")
{
  data = clusterdata$CONNECTORList
  if (is.null(clusterdata$FCM) & is.null(clusterdata$cluster)) {
    warning("In input is needed the FCM or StabilityAnalysis file. ", 
            immediate. = T)
  }
  else {
    if (!is.null(clusterdata$FCM)) 
      clusterdata <- clusterdata$FCM
  }
  ClustCurve <- clusterdata$cluster$ClustCurve
  ClustCurve <- merge(ClustCurve %>% select(ID,Cluster) %>% distinct() ,
                      data$LabCurv,
                      by = "ID")
  
  ClustCurve$Cluster <- clusterdata$cluster$cluster.names[ClustCurve$Cluster]
  countRes = ClustCurve[,c("ID", "Cluster", feature)]
  
  NAs = which(is.na(countRes[,feature]))
  if(length(NAs) > 0){
    countRes[NAs,feature] = "NA"
  }
  countRes[,feature] = as.factor(countRes[,feature])
  
  Counting <-  countRes %>%
    group_by_at(.vars = vars(c("Cluster", feature)), .drop = F)  %>% 
    tally %>%
    ungroup() %>%
    as.data.frame()
  
  return(list(Counting=Counting,
              ClusterNames=  ClustCurve[,c("ID", "Cluster")] ))
}

