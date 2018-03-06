#' Cluster Choice
#'
#'@description
#'
#'Calculates the BIC and AIC values considering k number of clusters and h dimension of the cluster mean space, and generates the plot based on the Elbow Method.
#'
#' @param data CONNECTORList.
#' @param k Number of clusters.
#' @param h Dimension of the cluster mean space. If NULL, ClusterChoice set the h value equals to the number of PCA components needed to explain the 95 perc. of variability of the data.
#' @param PCAperc The PCA percentages calculated with the function PCAbarplot, if NULL then it must be inserted in input a value for h.
#'  @return
#' The matrices of the AIC and BIC values, a list of FCMList objects belonging to class funcyOutList for each h and k and  the plot generated using the Elbow Method.
#' @examples
#'
#'GrowDataFile<-"data/1864dataset.xls"
#'AnnotationFile <-"data/1864info.txt"
#'
#'CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#'
#'CONNECTORList<- DataTruncation(CONNECTORList,"Progeny",truncTime=60,labels = c("time","volume","Tumor Growth"))
#'
#' # Specifying the h value
#'CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:6),h=2)
#'
#' # Using the PCA percentaes instead of the h value.
#'
#' perc <- PCA.Analysis(CONNECTORList$Dataset)$perc
#'
#' CONNECTORList.FCM <- ClusterChoice(CONNECTORList,k=c(2:6),PCAperc=perc)
#'
#' @import funcy ggplot2
#' @export
ClusterChoice<-function(data,k,h=NULL,PCAperc=NULL)
{
  K<-k
  database<-data$Dataset

  if(is.null(h))
  {
    if(is.null(PCAperc)) print("Choose a value for h or insert the PCA percentages.")
      else{ h<-min(which(cumsum(PCAperc)>=95)) }
  }

  H<-1:h

  output_k<-list()
  row_names <-c(paste("k=",K))
  col_names<-c(paste("h=",H))
  matrix_AIC<-matrix(0,nrow = length(K),ncol = length(H),dimnames=list(row_names,col_names))
  matrix_BIC<-matrix(0,nrow = length(K),ncol = length(H),dimnames=list(row_names,col_names))


  data.funcit <-matrix(c(database$ID,database$Vol,database$Time),ncol=3,byrow=F)

  Tot.within<-matrix(0,nrow = length(K),ncol = length(H),dimnames=list(row_names,col_names))

  # return a list of K lists, in which is is stored the output for all h
  # We also create two matrixes with the BIC and AIC values
  for(k in K)
  {
    output_h<-list()
    for(h in H)
    {
      mycontfclust = new("funcyCtrl",baseType="splines",dimBase=5,init="kmeans",nrep=10,redDim=h)
      out.funcit<- funcit(data.funcit,seed=2404,k,methods="fitfclust",funcyCtrl=mycontfclust,save.data=TRUE)
      output_h[[paste("h=",h)]] <- out.funcit

      minDist<-out.funcit@models$fitfclust@cldist[,1]
      Tot.within[which(K==k),which(H==h)]<-sum(minDist)

      matrix_BIC[which(K==k),which(H==h)]<-output_h[[paste("h=",h)]]@models$fitfclust@BIC
      matrix_AIC[which(K==k),which(H==h)]<-output_h[[paste("h=",h)]]@models$fitfclust@AIC
    }
    output_k[[paste("k=",k)]]<-output_h
  }

  Tot.within<-data.frame(dist=c(Tot.within),K=rep(K,length(H)),H=factor(rep(H,each=length(K))))
  ElbowMethod<-ggplot(data=Tot.within,aes(x=K))+ geom_point(aes(y=dist,col=H))+
                                    geom_line(aes(y=dist,col=H))+
                                    labs(title="Elbow method",x="Cluster",y="total within-cluster")+
                                    theme(text = element_text(size=20))




  return(list(FCM_all=output_k,matrix_BIC=matrix_BIC,matrix_AIC=matrix_AIC,ElbowMethod=ElbowMethod))
}
