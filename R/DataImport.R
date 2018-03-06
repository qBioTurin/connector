#' Data Import
#'
#'@description
#' Reads the files in a table format and create a list storing all the information.
#'
#'
#' @param file1 the name of the excel file which the data are to be read from. The growth data associated with an experiment must be reported using a pair of columns. The first column must contains the time points, suggested column header "Time". The second column must contains the data volume, the column header is the sample name.
#'
#' @param file2 the name of a cvs file which the annontation data are reported. Each row of the file containes: the identifier (ID) of the sample, e.g. an integer number, the sample name used in file1 and a list of features associated with the sample.
#'
#' @return List with four arguments: (i) data frame reporting three variables (ID, data and time values), (ii) the vector reporting the number of observations collected per sample, (iii) the data frame with curves labeled according to target file features and (iv) the vector for overall time grid.
#'
#' @examples
#'
#' GrowDataFile<-"data/1864dataset.xls"
#' AnnotationFile <-"data/1864info.txt"
#'
#' CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#'
#' @import readxl
#' @export
DataImport <- function(file1,file2) {
 ###Read Data File
  dataset <- read_excel(file1,col_names=T)
 ### Read Target File
  labcurv  <- read.csv(file=file2,header=TRUE)

 ###Check the column names
  c_names<-colnames(dataset[2*(1:(length(dataset[1,])/2))])
  if(length(c_names)!=(length(labcurv$ID)))
  {
    warning("Number of columns in the excel file is different from the number of samples stored in the target file.")

  }else{

    if(all(c_names==labcurv$SampleName)==FALSE)
    {
      warning("SampleNames in the target file do not correspond to the names in the excel file.")
    }
  }



  ### Inizialize :
  ### vector for curves lenghts
  nvar       <- dim(dataset)[2]
  nobs       <- dim(dataset)[1]
  tot        <- nobs*nvar
  samplesize <- dim(dataset)[2]/2
  lencurv    <- numeric(samplesize)
  ### vectors for times, volume and ID curves values
  tot        <- length(dataset)
  TimeIndex  <- seq(1,nvar,2)
  VolIndex   <- seq(2,nvar,2)
  TimeValue  <- as.matrix(dataset[TimeIndex])
  VolValue   <- as.matrix(dataset[VolIndex])

  times      <- numeric(tot)
  vol        <- numeric(tot)
  ID         <- numeric(tot)

  ### Organize times, volume and ID curves values, removing NA
  for (cappa in 1:samplesize)
  {
    tempv <- as.double(VolValue[,cappa])
    tempv <- tempv[!is.na(tempv)]
    lencurv[cappa] <- length(tempv)
    tempt <-as.double(TimeValue[,cappa])
    tempt <-tempt[!is.na(tempt)]

    if (cappa == 1)
    {
      vol[1:lencurv[cappa]]    <- tempv
      times[1:lencurv[cappa]]  <- tempt

    }

    else
    {
      lcum <- cumsum(lencurv)
      vol[(lcum[cappa-1]+1):lcum[cappa]] <- tempv
      times[(lcum[cappa-1]+1):lcum[cappa]] <- tempt
    }
  }

  ndata    <- sum(lencurv)
  vol      <- vol[1:ndata]
  times    <- times[1:ndata]
  ID       <- rep(labcurv$ID,times=lencurv)
  timegrid <- 1:max(times)

  ### ID, volume and time data frame

  dataset <- data.frame(ID=ID,Vol=vol,Time=times)
  alldata <- list(Dataset=dataset,LenCurv=lencurv,LabCurv=labcurv,TimeGrid=timegrid)


  return(alldata)
}
