#' Report generation
#'
#' @param data CONNECTORList. (see \code{\link{DataImport}})
#' @param stability.list The list obtained from the ClusterAnalysis function. (see \code{\link{ClusterAnalysis}})
#' @param feature The column name reported in the AnnotationFile  containing the feature to be investigated.
#' @param truncTime  A two dimension vector of integers corresponding to the time points where the curves will be truncated. If an integer number is passed, than it will be considered as the upper time point by default.
#' @param pRange The vector of the dimension of the natural cubic spline basis, by which the cross-validated loglikelihood is calculated and plotted. It could be the path of the Rds file storing the output of the \code{\link{BasisDimension.Choice}} function. If NULL, the p selection step is omitted.
#' @param p The dimension of the natural cubic spline basis. (see \code{\link{BasisDimension.Choice}})
#' @param G  Vector of integers representing the number of clusters.
#' @param path The folder path  where the report will be saved. If it is missing, the report is saved in the current working  directory.
#' @param namefile Report name.
#' 
#' @return  
#' @examples
#'
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#' 
#' @seealso  \code{\link{GrowthCurve}}, code{\link{TimeGridDensity}}.
#' @import markdown parallel
#' @importFrom knitr kable
#' @export

ReportGeneration <- function(data, stability.list, G, feature = "ID", namefile = "Report", path = NULL,truncTime = NULL, p=NULL,pRange=NULL)
{
  # Copy the report file to a temporary directory before processing it, in
  # case we don't have write permissions to the current working dir (which
  # can happen when deployed).
  report.path <- system.file("Shiny","report.Rmd", package = "connector")
  tempReport <- file.path(tempdir(), "report.Rmd")
  file.copy(report.path, tempReport, overwrite = TRUE)
  # Set up parameters to pass to Rmd document
  
  if(is.null(path)){
    path = getwd()
  }
  
  infoReport <-  list(data = data,
                      clusterdata = stability.list,
                      G = G,
                      p = p,
                      feature = feature,
                      CrossLL = pRange,
                      truncTime = truncTime)
  
  # infoReport <-  list(data = CONNECTORList,
  #                     clusterdata = S.cl,
  #                     p = p,
  #                     G = 4,
  #                     feature = "Progeny",
  #                     CrossLL = 3:7,
  #                     truncTime = 70)
  # Knit the document, passing in the `params` list, and eval it in a
  # child of the global environment (this isolates the code in the document
  # from the code in this app).
  rmarkdown::render(output_format = "html_document",
                    output_file = paste0(namefile,".html"),
                    output_dir = path,
                    input =  report.path,
                    params = list(infoReport = infoReport),
                    envir = new.env(parent = globalenv())
  )
  
}
 

 
