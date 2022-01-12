#' Report generation
#'
#' @param data CONNECTORList. (see \code{\link{DataImport}})
#' @param clusterdata The list obtained from extrapolating the most probable clustering from the StabilityAnalysis function output. (see \code{\link{StabilityAnalysis}} and \code{\link{MostProbableClustering.Extrapolation}}).
#' @param feature The column name reported in the AnnotationFile  containing the feature to be investigated.
#' @param truncTime
#' @param p The vector of the dimension of the natural cubic spline basis, by which the cross-validated loglikelihood is calculated and plotted. It could be the path of the Rds file storing the output of the \code{\link{BasisDimension.Choice}} function. If NULL, the p selection step is omitted.
#' @param pRange
#' @param path The folder path  where the plot will be saved. If it is missing, the plot is saved in the current working  directory.
#' @param file
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

ReportGeneration <- function(data, clusterdata, G, feature = "ID", file = "Report.pdf",truncTime = NULL, p=NULL,pRange=NULL)
{
  # Copy the report file to a temporary directory before processing it, in
  # case we don't have write permissions to the current working dir (which
  # can happen when deployed).
  report.path <- system.file("Shiny","report.Rmd", package = "connector")
  tempReport <- file.path(tempdir(), "report.Rmd")
  file.copy(report.path, tempReport, overwrite = TRUE)
  # Set up parameters to pass to Rmd document
  
  infoReport <-  list(data = data,
                      clusterdata = clusterdata,
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
                    output_file = "Report.html",output_dir = "inst/Shiny/",
                    input = "inst/Shiny/report.Rmd",
                    params = list(infoReport = infoReport),
                    envir = new.env(parent = globalenv())
  )
  
}
 

 
