#' Withinness and Betweenness Visualization
#'
#' Visualization of the plots reporting the withinness and betweenness clustering measures.
#'
#' @param ClustCurve Data frame with five arguments : time, growth data, sample ID, cluster membership and feature values for each sample.
#' @param MeanCurves Matrix with the meancurve values on the columns according to different clusters.
#' @param Title Title to associate to the withinness and Betweenness plot.
#' @param save If TRUE the Betweenness and Withinness plot is saved in a pdf file.
#' @param path Path to save plot to (combined with file name). If it is missing, the plot is saved in the working  directory.
#' @return Returns the Withinness and Betweenness plot.
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
#' CONNECTORList.models <- FittingAndClustering( data= CONNECTORList, FCM_all = CONNECTORList.FCM, h = 2, k=4, feature = "Progeny", labels = c("time","volume"))
#' 
#' ### Withinness and betweenness plot for Matlhus model and FCM
#' 
#' Malthus.ClustCurve <-  CONNECTORList.models$Malthus$Information$ClustCurve
#' Malthus.MeanCurves <-  CONNECTORList.models$Malthus$Information$meancurves
#' PlotWithinnessBetweenness(Malthus.ClustCurve,Malthus.MeanCurves,Title = "Malthus Cluster betweenness and withinness")
#' 
#' FCM.ClustCurve <-  CONNECTORList.models$FCM$Information$ClustCurve
#' FCM.MeanCurves <-  CONNECTORList.models$FCM$Information$meancurves
#' 
#' PlotWithinnessBetweenness(FCM.ClustCurve,FCM.MeanCurves,Title = "FCM Cluster betweenness and withinness")
#' 
#' @import ggplot2 ggforce
#' @export
PlotWithinnessBetweenness <- function(ClustCurve,MeanCurves,Title=NULL,save=TRUE,path=NULL)
{

  K <- length(unique(ClustCurve[,4]))
  ClustSymbol<-cluster.symbol(K)
  nfeature <- length(unique(ClustCurve[,5]))[1]
  cluster.palette <- rainbow(K)
  #feature.palette <- rainbow(nfeature+3)
  feature.lev <- sort(unique(ClustCurve[,5]))

  ## calculate the max, min, mean and sd distance in the k-clusters;
  within.all <- Withinness(ClustCurve,MeanCurves,centroids=TRUE)
  ## calculate the nearest cluster
  between.all <- Betweenness(ClustCurve,MeanCurves)$Betweenness

  if(length(which(table(between.all[,2])==1))== 0 )
      {first<-names(which(table(c(between.all[,2],rownames(between.all)))==1)[1])}
  else{ first<-names(which(table(between.all[,2])==1)[1])}
  ## sort the cluster and save the distance among them
  shift<-sort(Betweenness(ClustCurve,MeanCurves)$CentroidDist[paste(first),])

  ClustSymbol.sorted <- names(shift)
  i<-match(ClustSymbol.sorted, ClustSymbol)

  cluster.magnitudo <- apply(table(unique(ClustCurve[,c(4,1)])),1,sum)[i]
  circles <- data.frame(
  x0 = numeric(3*K)     ,
  y0 = numeric(3*K)     ,
  r = numeric(3*K)      ,
  distance = numeric(3*K),
  Cluster = numeric(3*K)
  )
  colnames(circles) <- c("x0","y0","r","distance","Cluster")

  n <- length(unique(ClustCurve[,1]))
  WithDist <- data.frame(
  x1 = numeric(n),
  y1 = numeric(n),
  Cluster = numeric(n),
  feature = numeric(n)
	)
  colnames(WithDist) <- c("x1","y1","Cluster","feature")
  counter <- 1
  index <- matrix(seq(1,3*K),nrow=K,byrow=TRUE)

  for(k in 1:K)
  {
  ### Data frames
    dataplot <- DataFrameWithinness.i(ClustCurve,MeanCurves,i[k],ClustSymbol,shift=shift[k])
    circles[index[k,],] <- dataplot$circles
    WithDist[counter:(cumsum(cluster.magnitudo)[k]),] <- dataplot$WithDist
    counter <- cumsum(cluster.magnitudo)[k] + 1
  }

  circles$distance <- factor(circles$distance)
  circles$Cluster <- factor(circles$Cluster)

  WithDist$Cluster <- factor(WithDist$Cluster)
  WithDist$feature <- factor(WithDist$feature)

  plots <- ggplot() + geom_circle(aes(x0=x0, y0=y0, r=r,linetype=distance,color=Cluster), data=circles,size=1)
  plots <- plots + scale_shape_manual(values=c(0:(K-1)))+
                    geom_point(data=WithDist,aes(x=x1,y=y1,shape=Cluster),size=4)+
                    scale_linetype_manual("",values = c( "2"= "dashed","1"="solid"),
                                           labels=c("Mean ","Mean+- sd"))

  if(is.null(Title)) Title<-"Cluster betweenness and withinness"

  plots <- plots  + geom_text(aes(x=x0, y=y0,label=Cluster), data=circles) +
                    #scale_colour_manual(values = cluster.palette,name="Cluster Colors") +
                    labs(title=Title,x="distance",y="distance")+
                    theme(plot.title = element_text(hjust = 0.5),                                                                      axis.line = element_line(colour = "black"))


  if(save==TRUE)
  {
   if(is.null(path)) path <- getwd()
   ggsave(filename="Betweenness&Withinness.pdf",plot =plots,width=29, height = 20, units = "cm",scale = 1,path=path)
  }

  return(plots)
}

