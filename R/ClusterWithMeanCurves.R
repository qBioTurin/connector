#' Mean Cluster Curves
#'
#'@description
#' Visualize the plot of the mean cluster curves given one model (i.e. FCM, Malthus, Gompertz or Logistic).
#' ClusterWithMeanCurve  generates k plots, where k is the number of clusters: in each plot the curves are colored according to the feature chosen by the user.
#'
#' @param clusterdata.FCM Object belonging to the class funcyOutList.
#' @param clusterdata.Models Object, or a list of objects storing the parameters, clusters, centers, meancurves and seed derived from fitting and clustering the data using Malthus, Gompertz or/and Logistic model.
#' @param data CONNECTORList.
#' @param model  Model name, i.e. FCM, Malthus, Gompertz or Logistic. The default is "ALL" that runs all the models but it needs in input both clusterdata.FCM and clusterdata.Models .
#' @param feature The column name reported in the AnnotationFile containing the feature interesting for the user to be investigated.
#' @param labels  Vector containing the text for the title of axis.
#' @param save
#' @return The plot reporting the mean curves of all clusters, for each cluster the plot of the growth curves and a list of data reporting the mean curve values and for each sample the ID of the belonging cluster.
#'
#' @examples
#'
#'GrowDataFile<-"data/1864dataset.xls"
#'AnnotationFile <-"data/1864info.txt"
#'
#'CONNECTORList <- DataImport(GrowDataFile,AnnotationFile).............
#'
#' @import ggplot2 cowplot
#' @export
ClusterWithMeanCurve<-function(clusterdata, data, feature ,title="", labels=c("","") ,path=NULL ,save=FALSE)
{
axis.x<-labels[1]
axis.y<-labels[2]
  
    if(isS4(clusterdata))
    {
      k<-clusterdata@k
      
      Cluster(clusterdata)->classes
      out.fit<-clusterdata@models$fitfclust@fit
      # Check if it is regular
      if(clusterdata@reg==1)
      {
        fitfclust.curvepred(out.fit)$meancurves->meancurves
      } else{
        fitfclust.curvepredIrreg(out.fit)$meancurves->meancurves
      }
    }else{
      
      clusterdata -> classification
      classification$cluster ->classes
      classification$meancurves->meancurves
      length(meancurves[1,])->k
    }
    
    symbols<-cluster.symbol(k)
    
    ##################
    
    classificate <- rep(classes,data$LenCurv)
    curves <- data.frame(Times=data$Dataset$Time,Vol=data$Dataset$Vol,ID=data$Dataset$ID,Cluster=classificate, Info=rep(t(data$LabCurv[feature]),data$LenCurv))
  
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
    PlotMeanCurve<-ggplot()+
                      geom_line(data=plot_data, aes(x=time,y=means,group=clusters) )+
                      labs(title=title, x=axis.x, y = axis.y)+
                      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())
  
      col<-as.character(unique(curves$Info))
      col1<-rainbow(length(col))
      plots<-list()
      ymax<-max(curves$Vol)
      for(i in 1:k)
      {
        plots[[paste(symbols[i],"Cluster")]]<-ggplot()+
          geom_line(data=plot_data[plot_data$clusters==i,], aes(x=time,y=means),size = 1.2 )+
          labs(title=paste(title,"",symbols[i],"Cluster"), x=axis.x, y = axis.y)+
          geom_line(data = curves[curves$Cluster==i,],aes(x=Times,y=Vol,group=ID,color=factor(Info)))+
          scale_colour_manual(values = col1,limits=col,breaks=col,name=feature)+
          theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())+
          ylim(0,ymax)
      }
       plots[["ALL"]]<-plot_grid(plotlist = plots)
       
       out<-list(plotMeanCurve=PlotMeanCurve,plotsCluster=plots)

  if(save)
  {
    if(is.null(path)) path <- getwd()
    for(i in 1:k)
    {
    ggsave(filename = paste(symbols[i],"Cluster.pdf",sep=""),plot=out$plotsCluster[[paste(symbols[i],"Cluster")]],width=29, height = 20, units = "cm",scale = 1,path = path)
    }
    ggsave(filename = "ALLCluster.pdf",plot=out$plotsCluster$ALL,width=29, height = 20, units = "cm",scale = 1,path = path)
    ggsave(filename = "MeanCurve.pdf",plot=out$plotMeanCurve,width=29, height = 20, units = "cm",scale = 1,path = path)
}

 return(out)
}


