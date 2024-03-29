#' Return the plots of the spline exploited to fit and cluster the samples.
#' 
#' @description 
#' Returns the Cubic Spline plots from the ClusterWithMeanCurve output list (see \code{\link{ClusterWithMeanCurve}}). Each plot shows (i) in blue the sample curve, (ii) in red the cubic spline estimated from the FCM, (iii) in black the correspondive cluster mean curve, and finally (iv) the grey area represents the confidence interval.
#'
#' @param FCM.plots CONNECTORList obtained from the ClusterWithMeanCurve function. (see \code{\link{ClusterWithMeanCurve}})
#' @param All If TRUE then all the plots, for each sample curve, are saved.
#' @param SampleNumber The number/vector of the sample ID(s) referring to the curve(s) that will be saved.
#' @param path The folder path where the plot(s) will be saved. If it is missing, the plot is saved in the current working  directory.
#' @param save If TRUE then the spline plots are saved (a pdf file will be generated for each curve). 
#'
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @import ggplot2
#' @export
#' 
Spline.plots <- function(FCM.plots, All=TRUE, SampleNumber = NULL, save = FALSE, path = NULL){
  ListFitting = list()
  n.samples<-length(FCM.plots$spline.plots)
  
  if(All){ 
    samplesIndexes<- 1:n.samples
  }else if(!is.null(SampleNumber)) {
    # select only the sample number between 1 and n.samples
    samplesIndexes<-SampleNumber[which( SampleNumber %in% 1:n.samples == TRUE)]
    if(length(samplesIndexes)==0) {
      warning(paste0("The number of the sample must be between 0 and ",n.samples,", that is the number of samples. It will be saved the first plot."),immediate. = TRUE)
      samplesIndexes<-1
    }
  }else{
    warning(paste0("Select a number of the sample between 0 and ",n.samples,", or put the parameter All = True. Since neither the parameter All or SampleNumber were selected, it will be saved only the first curve."),immediate. = TRUE)
    samplesIndexes<-1
  }
  
  if(is.null(path)) path <- getwd()
  
  for(i in samplesIndexes)
  {
    ListFitting[[paste(i)]] = FCM.plots$spline.plots[[paste("Sample ",i)]]
    
    if(save)
      ggsave(filename = paste("Spline",i,"sample.pdf",sep="_"),plot=FCM.plots$spline.plots[[paste("Sample ",i)]],width=29, height = 20, units = "cm",scale = 1,path = path)
  }
  return(ListFitting)
}