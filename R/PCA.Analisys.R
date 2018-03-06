#' PCA Analysis
#'
#'@description
#'
#' Generates a bar plot reporting the percentage values associated to each component identified by the Principal Component Analysis (PCA).
#'
#' @param data CONNECTORList.
#' @param dimBase Dimension of the basis.
#' @param save If TRUE the bar plot of the PCA components is saved in a pdf file.
#' @param path  Path to save plot to (combined with file name). If it is missing, the plot is saved in the working  directory.
#' @return The bar plot of the PCA components and the vector of percentage values of the PCA components.
#' @examples
#'
#' GrowDataFile<-"data/1864dataset.xls"
#' AnnotationFile <-"data/1864info.txt"
#'
#' CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#'
#' CONNECTORList<- DataTruncation(CONNECTORList,"Progeny",truncTime=60,labels = c("time","volume","Tumor Growth"))
#'
#' PCA<-PCA.Analysis(CONNECTORList$Dataset)
#' PCA$plot
#'
#' @import ggplot2
#' @export
PCA.Analysis <- function(data,dimBase=5,save=FALSE,path=NULL)
{
  TimeGrid <- c(1:max(data[,3]))

  data<-as.matrix(data)
  # curves splines basis coefficients
  res <- makeCoeffs(data=data, reg=FALSE, dimBase=dimBase,
                     grid=TimeGrid, pert=0.01)

  # Principal Components Analysis

  princomp(as.matrix(res$coeffs)) -> pca
  # Number of principal components
  ncomp <- length(names(pca$sdev))
  # Principal components variances
  eigs <- pca$sdev^2
  # Percentage of variances explained by each component
  percentage <- eigs/sum(eigs)*100

  # PCA bar plot
  PCA_barplot<-ggplot(data=data.frame(comp=paste("Comp.",1:ncomp),Variances=eigs,perc=paste(signif(percentage,4),"%",sep="")), aes(x=comp, y=Variances)) +
    geom_bar(stat="identity", fill="steelblue")+
    geom_text(aes(label=perc), vjust=-.3,  size=3.5)+
    labs(title="PCA barplot", x="Components", y = "Variances")+
    theme(plot.title = element_text(hjust = 0.5))

  if(save==TRUE)
  {
    ggsave(filename="PCA_Analysis.pdf",plot =PCA_barplot,width=29, height = 20, units = "cm",scale = 1,path = path)
  }
  return(list(plot=PCA_barplot,perc=percentage))
  }

