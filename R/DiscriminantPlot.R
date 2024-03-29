#' Discriminant Plot
#'
#'@description
#' Gives a simple visualization regarding in which cluster each sample belongs.
#'
#' @param clusterdata The list obtained from extrapolating the most probable clustering from the StabilityAnalysis function output. (see \code{\link{StabilityAnalysis}} and \code{\link{MostProbableClustering.Extrapolation}}). 
#' @param data CONNECTORList. (see \code{\link{DataImport}} or \code{\link{DataTruncation}})
#' @param save If TRUE then the discriminant plots colored depending on the feature and the cluster membership are saved. 
#' @param path The folder path where the plot(s) will be saved. If it is missing, the plot is saved in the current working  directory.
#' 
#'  @return
#'  DiscriminantPlot visualizes two versions of the same discriminant linear plot, showing in a h-dimensional space the cluster membership of each curve.
#'  In the first case the symbols identifying the elements of different clusters are colored depending on the cluster membership,
#'  and in the second plot the symbols are colored depending the selected feature.
#'  If h > 2 the principal component analysis is exploited to extrapolate the components of the h-space with more explained variability.
#' 
#' @references
#' Gareth M. James and Catherine A. Sugar, (2003). Clustering for Sparsely Sampled Functional Data. Journal of the American Statistical Association.
#'
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'
#' @import  ggplot2 fda
#' @export
#' 

DiscriminantPlot<-function(clusterdata,feature,save=FALSE,path=NULL)
{ 
  data = clusterdata$CONNECTORList
  G<- length(clusterdata$FCM$prediction$meancurves[1,])
  h<- length(clusterdata$FCM$fit$parameters$Lambda[1,])
  discrplot<-list()
  
  outpred<-clusterdata$FCM$cluster$cl.info
  
  projectedcurve=outpred$alpha.hat
  projectedclustcenters=clusterdata$FCM$fit$parameters$alpha
  outpred$Calpha -> stdevalpha
  outpred$class.pred -> classes
  outpred$distance -> distanze
  Feature <- data$LabCurv[paste(feature)]
  
  col<- as.character(unlist(unique(Feature)) )
  colFetaure <- rainbow(dim(unique(data$Lab[feature]))[1])
  
  symbols<-clusterdata$FCM$cluster$cluster.names
  
  if(h==1)
  {     
    DataFrame<-data.frame(PrjCurv=projectedcurve,stdaplh=stdevalpha,Feature=unlist(Feature),Cluster=symbols[classes] )
    
    cl.names<-1:G-1
    names(cl.names)<-symbols
    
    discrplot[["ColCluster"]]<-ggplot(data = DataFrame)+
      geom_point(aes(x=PrjCurv,y=stdaplh,colour=Cluster,pch=Cluster),size=4)+
      geom_vline(xintercept = projectedclustcenters)+
      scale_shape_manual("Cluster",values=cl.names)+ 
      labs(title="Discriminant plot")+
      xlab("Alpha")+ylab('Standard Deviation')+
      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())
    
    
    discrplot[["ColFeature"]]<-ggplot(data = DataFrame)+
      geom_point(aes(x=PrjCurv,y=stdaplh,colour=as.factor(Feature),shape=Cluster),size=4)+
      geom_vline(xintercept = projectedclustcenters)+
      xlab("Alpha")+ylab('Standard Deviation')+
      labs(title="Discriminant plot")+
      scale_colour_manual(values = colFetaure,limits=col,breaks=col,name=feature )+
      scale_shape_manual("Cluster",values=cl.names) +
      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())
  }
  else if(h==2){
    DataFrameSamples<-data.frame(PrjCurv1=projectedcurve[,1],PrjCurv2=projectedcurve[,2],Feature=unlist(Feature),Cluster=symbols[classes] )
    DataFrameCluster<-data.frame(projectedclustcenters1=projectedclustcenters[,1],projectedclustcenters2=projectedclustcenters[,2],Cluster=symbols,Center=paste("c.",symbols,sep=""))
    
    cl.names<-1:G-1
    names(cl.names)<-symbols
    
    discrplot[["ColCluster"]]<-ggplot()+
      geom_point(data = DataFrameSamples,aes(x=PrjCurv1,y=PrjCurv2,colour=Cluster,pch=Cluster),size=4)+
      geom_text(data = DataFrameCluster,aes(x=projectedclustcenters1,y=projectedclustcenters2,label=paste("c",symbols,sep=""),colour=Cluster),size=5,show.legend = F)+
      scale_shape_manual("Cluster",values=cl.names,labels=paste(names(cl.names)," (c",symbols,")",sep=""),breaks=names(cl.names)) +
      scale_colour_discrete("Cluster",limits=names(cl.names),breaks=names(cl.names),labels=paste(names(cl.names)," (c",symbols,")",sep=""))+
      xlab("Alpha 1")+ylab('Alpha 2')+
      labs(title="Discriminant plot")+
      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())
    
    discrplot[["ColFeature"]]<-ggplot()+geom_point(data = DataFrameSamples,aes(x=PrjCurv1,y=PrjCurv2,colour=as.factor(Feature),pch=Cluster),size=4)+
      geom_text(data = DataFrameCluster,aes(x=projectedclustcenters1,y=projectedclustcenters2,label=Center),size=5,show.legend = F)+
      xlab("Alpha 1")+ylab('Alpha 2')+
      labs(title="Discriminant plot")+
      scale_colour_manual(values = colFetaure,limits=col,breaks=col,name=feature )+
      scale_shape_manual("Cluster",values=cl.names,labels=paste(names(cl.names),", c.",symbols,sep=""),breaks=names(cl.names)) +
      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())
    
  }else{
    colnames(projectedcurve) = paste0("Alpha",1:length(projectedcurve[1,]))
    prcomp(projectedcurve) -> pca
    # Principal components variances
    eigs <- pca$sdev^2
    # Percentage of variances explained by each component
    percentage <- eigs/sum(eigs)*100
    
    # all.equal(pca$x,
    #           (projectedcurve -
    #              matrix(pca$center,ncol = 3,nrow=length(projectedcurve[,1]),byrow = T)) %*%
    #             pca$rotation)
    
    projectedclustcentersPCA =
              (projectedclustcenters -
                 matrix(pca$center,ncol = h,nrow=length(projectedclustcenters[,1]),byrow = T)) %*%
                pca$rotation
    
    DataFrameSamples<-data.frame(PrjCurv1=pca$x[,"PC1"],
                                 PrjCurv2=pca$x[,"PC2"],
                                 Feature=unlist(Feature),
                                 Cluster=symbols[classes] )
    DataFrameCluster<-data.frame(projectedclustcenters1=projectedclustcentersPCA[,"PC1"],
                                 projectedclustcenters2=projectedclustcentersPCA[,"PC2"],
                                 Cluster=symbols,
                                 Center=paste("c.",symbols,sep=""))
    
    cl.names<-1:G-1
    names(cl.names)<-symbols
    
    discrplot[["ColCluster"]]<-ggplot()+
      geom_point(data = DataFrameSamples,aes(x=PrjCurv1,y=PrjCurv2,colour=Cluster,pch=Cluster),size=4)+
      geom_text(data = DataFrameCluster,aes(x=projectedclustcenters1,y=projectedclustcenters2,label=paste("c",symbols,sep=""),colour=Cluster),size=5,show.legend = F)+
      scale_shape_manual("Cluster",values=cl.names,labels=paste(names(cl.names)," (c",symbols,")",sep=""),breaks=names(cl.names)) +
      scale_colour_discrete("Cluster",limits=names(cl.names),breaks=names(cl.names),labels=paste(names(cl.names)," (c",symbols,")",sep=""))+
      xlab(paste0("PC1 (" , round(percentage[1]), "% )" ) ) +
      ylab(paste0("PC2 (" , round(percentage[2]), "% )" ) ) +
      labs(title="Discriminant plot")+
      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())
    
    discrplot[["ColFeature"]]<-ggplot()+geom_point(data = DataFrameSamples,aes(x=PrjCurv1,y=PrjCurv2,colour=as.factor(Feature),pch=Cluster),size=4)+
      geom_text(data = DataFrameCluster,aes(x=projectedclustcenters1,y=projectedclustcenters2,label=Center),size=5,show.legend = F)+
      xlab(paste0("PC1 (" , round(percentage[1]), "% )" ) ) +
      ylab(paste0("PC2 (" , round(percentage[2]), "% )" ) ) +      labs(title="Discriminant plot")+
      scale_colour_manual(values = colFetaure,limits=col,breaks=col,name=feature )+
      scale_shape_manual("Cluster",values=cl.names,labels=paste(names(cl.names),", c.",symbols,sep=""),breaks=names(cl.names)) +
      theme(plot.title = element_text(hjust = 0.5),axis.line = element_line(colour = "black"),panel.background = element_blank())
    
  }
  
  if(save==TRUE)
  {
    if(is.null(path))
    {
      path <- getwd()
    }
    ggsave(filename="DiscrPlot_ColFeature.pdf",plot =discrplot[["ColFeature"]],width=29, height = 20, units = "cm",scale = 1,path=path )
    ggsave(filename="DiscrPlot_ColClust.pdf",plot =discrplot[["ColCluster"]],width=29, height = 20, units = "cm",scale = 1,path=path )
  }
  
  return(discrplot)
}
