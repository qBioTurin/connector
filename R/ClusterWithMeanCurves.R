#' Mean Cluster Curves
#'
#'@description
#' Visualize the plot of the mean cluster curves given one model (i.e. FCM, Malthus, Gompertz or Logistic).
#' ClusterWithMeanCurve  generates k plots, where k is the number of clusters: in each plot the curves are colored according to the feature chosen by the user.
#'
#' @param clusterdata Object belonging to the class funcyOutList.
#' @param data CONNECTORList.
#' @param k Number of clusters.
#' @param model  Model name, i.e. FCM, Malthus, Gompertz or Logistic.
#' @param feature The column name reported in the AnnotationFile containing the feature interesting for the user to be investigated.
#' @param labels  Vector containing the text for the title of axis.
#' @return The plot reporting the mean curves of all clusters, for each cluster the plot of the growth curves and a list of data reporting the mean curve values and for each sample the ID of the belonging cluster.
#'
#' @examples
#'
#'GrowDataFile<-"data/1864dataset.xls"
#'AnnotationFile <-"data/1864info.txt"
#'
#'CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#'
#'CONNECTORList<- DataTruncation(CONNECTORList,"Progeny",truncTime=60,labels = c("time","volume","Tumor Growth"))
#'
#'# Considering the FCM
#'
#' CONNECTORList.FCM.k4.h2<- CONNECTORList.FCM$FCM_all$`k= 4`$`h= 2`
#'
#' FCMplots <- ClusterWithMeanCurve(CONNECTORList.FCM.k4.h2,CONNECTORList,k = 4,model="FCM",feature = "Progeny",labels = c("Time","Volume"))
#'
#'# Considering the Malthud model
#'
#' MalthusPlots<- ClusterWithMeanCurve(database = CONNECTORList,k = 4,model="Malthus",feature = "Progeny")
#'
#' @import ggplot2 cowplot
#' @export
ClusterWithMeanCurve<-function(clusterdata=NULL,data,k,model,feature,labels=NULL)
{
  if(is.null(clusterdata) & model=="FCM")
  {
    warning("The model FCM needs an object of the class funcyOutList." )
    break
  }

  if(is.null(labels))
  {
    axis.x<-""
    axis.y<-""
    title<-""

  }else{
    axis.x<-labels[1]
    axis.y<-labels[2]
    title<-labels[3]
  }

  symbols<-cluster.symbol(k)
  Information<-list()
  time <- sort(unique(data$Dataset$Time))

  if(model=="FCM")
  {
    Cluster(clusterdata)->classes->Information$classes
    out.fit<-clusterdata@models$fitfclust@fit

    # Check if it is regular
    if(clusterdata@reg==1)
    {
      fitfclust.curvepred(out.fit)$meancurves->meancurves->Information$meancurves
    } else{
      fitfclust.curvepredIrreg(out.fit)$meancurves->meancurves->Information$meancurves
    }
  }
  else{

    clustering(data,k,model) ->classification
    classification$cluster ->classes-> Information$classes
    classification$center -> Information$center
    classification$meancurves->meancurves->Information$meancurves
  }
  classificate <- rep(classes,data$LenCurv)
  curves <- data.frame(Times=data$Dataset$Time,Vol=data$Dataset$Vol,ID=data$Dataset$ID,Cluster=classificate, Info=rep(t(data$LabCurv[feature]),data$LenCurv))
  Information$ClustCurve <- data.frame(merge(curves[,1:4],data$LabCurv,by="ID"))

  # cut the meancurves at the growth curves' maximum time
  time1<-sort(unique(data$Dataset$Time))
  meancurves_truncated<-c()
  time3<-c()
  cluster<-c()
  for(clust in 1:k)
  {
    time2<-sort(unique(curves[curves$Cluster==clust,]$Times))
    m<-meancurves[,clust][time1<=max(time2)]
    time3<-c(time3,time1[time1<=max(time2)])
    meancurves_truncated<-c(meancurves_truncated,m)
    cluster<-c(cluster,rep(clust,length(time1[time1<=max(time2)])))
  }



  plot_data<-data.frame(time=time3,means=meancurves_truncated,clusters=cluster)
  PlotMeanCurveFCM<-ggplot()+
                    geom_line(data=plot_data, aes(x=time,y=means,group=clusters) )+
                    labs(title=paste(model," cluster mean curves"), x=axis.x, y = axis.y)+
                    theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())

    col<-as.character(unique(curves$Info))
    col1<-rainbow(length(col))
    plots<-list()
    ymax<-max(curves$Vol)
    for(i in 1:k)
    {
      plots[[paste(symbols[i],"Cluster")]]<-ggplot()+
        geom_line(data=plot_data[plot_data$clusters==i,], aes(x=time,y=means),size =1.3 )+
        labs(title=paste(model,"",symbols[i],"Cluster"), x=axis.x, y = axis.y)+
        geom_line(data = curves[curves$Cluster==i,],aes(x=Times,y=Vol,group=ID,color=factor(Info)))+
        scale_colour_manual(values = col1,limits=col,breaks=col,name=feature)+
        theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())+
        ylim(0,ymax)
    }
     plots[["ALL"]]<-plot_grid(plotlist = plots)
     plots$ALL

 return(list(plotMeanCurve=PlotMeanCurveFCM,plotsCluster=plots,Information=Information))
}


