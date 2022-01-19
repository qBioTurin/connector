#' Choice of the B-spline dimension 
#'
#'@description
#'	Generates a line plot reporting the cross-validated loglikelihood value for each number of knots. In details, for each number of knots 10\% of the curves from the whole data are removed and treated as a test set, then the remaing curves are fitted using the FCM and the loglikelihood on the test set is calculated. The process is then repeated nine more times.
#'  
#'
#' @param data CONNECTORList. (see \code{\link{DataImport}})
#' @param p The vector of the dimension of the natural cubic spline basis.
#' @param save If TRUE then the plot is saved into a pdf file.
#' @param path The folder path  where the plot will be saved. If it is missing, the plot is saved in the current working  directory.
#' @param Cores Number of cores to parallelize computations.
#' 
#' @return
#' DimensionBasis.Choice returns line plot of the cross-validated loglikelihood for each value of p, in grey the result of all ten repetitions of the likelihood calculation and in black the mean of them.
#'
#' 
#'  
#' @references
#' Gareth M. James and Catherine A. Sugar, (2000). Principal component models for sparse functional data.
#' 
#' @examples
#'
#'GrowDataFile<-"data/745dataset.xls"
#'AnnotationFile <-"data/745info.txt"
#'
#'CONNECTORList <- DataImport(GrowDataFile,AnnotationFile)
#'
#'CONNECTORList<- DataTruncation(CONNECTORList,"Progeny",truncTime=50,labels = c("time","volume","Tumor Growth"))
#'
#'CrossLogLike<-BasisDimension.Choice(CONNECTORList,2:6)
#'CrossLogLike$CrossLogLikePlot
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'
#' @import ggplot2 fda parallel grid
#' @importFrom ggplotify  as.ggplot
#' @importFrom plyr ldply
#' @importFrom MASS ginv
#' @export
#' 
BasisDimension.Choice<-function(data,p,save=FALSE,path=NULL,Cores = 1)
{
  samples <- unique(data$Dataset$ID)
  p.values<-p
  crossvalid<-list()
  n_sample<-length(data$LenCurv)
  perc <- as.integer(n_sample*0.1)
  if(perc<1) perc<-1
  
  grid<-data$TimeGrid
  length(grid)->Lgrid
  
  pmax<-min(data$LenCurv)
  
  # if(max(p)>pmax)
  # {
  #   warning(paste("The maximum value of p should be:", pmax,"\n") )
  # }
  
  if(max(p)>Lgrid)
  {
    p.values<-min(p):Lgrid
    warning(paste("The maximum value of p should not be larger than the length of the grid. Thus the range of p is fixed from",min(p),"to",Lgrid,".\n") )
  }
  
  nworkers <- detectCores()
  if(nworkers<Cores) Cores <- nworkers
  if(Cores > 10) Cores <- 10
  
  cl <- makeCluster(Cores)
  crossvalid<-parLapply(cl,1:10, function(step){
    tryCatch({
    SampleTestSet<-sample(samples,perc)
    SampleTestSet<-SampleTestSet[order(SampleTestSet)]
    
    TestSet<-data$Dataset[data$Dataset$ID%in%SampleTestSet,]
    TestSet <-data.frame(ID=rep(1:(length(SampleTestSet)),data$LenCurv[which(samples %in% SampleTestSet)]),
                         Vol=TestSet$Vol[],
                         Time=TestSet$Time)
    
    TrainingSet<-data$Dataset[-which(data$Dataset$ID%in%SampleTestSet),]
    
    #### Let define the FCM input matrix of ID sample, Points and time indexes per sample into the grid.
    
    data.funcit <-matrix(c(rep(1:(n_sample-length(SampleTestSet)),data$LenCurv[-which(samples %in% SampleTestSet)]),TrainingSet$Vol,match(TrainingSet$Time,grid)),
                         ncol=3,byrow=F)
    
    CrossLikelihood<-sapply(p.values, CalcLikelihood, data.funcit,TestSet,grid)
    
    return(data.frame(p=p.values,CrossLikelihood=CrossLikelihood,sim=step))
    }, error=function(e) {
      err<-paste("ERROR in prediction :",conditionMessage(e), "\n")
      #print(err)
      return(NULL)
    }
    )
  })
  stopCluster(cl)
  crossvalid <- crossvalid[!sapply(crossvalid,is.null)] 
  
  Knots.list<-lapply(p.values,function(p){
      Spline<-ns(grid, df = (p - 1))
      df<-data.frame(Knots= c(attr(Spline,"Boundary.knots"),attr(Spline,"knots")),
                     Name=c(rep("Boundary knots",2),rep("Knots",length(attr(Spline,"knots"))) ) )
      df$p = paste("p = ",p)
      return(df)
  })
  Knots.df <- do.call("rbind",Knots.list)
  Knots.df$p.num <- as.numeric(sub("p = ","",Knots.df$p))

  Knots.Plot<-ggplot(Knots.df) +
    scale_x_continuous(breaks = as.integer(seq(min(grid),
                                               max(grid),
                                               length.out = length(grid)/5))) +
    geom_hline(aes(yintercept = p.num), linetype = "dashed", color = "grey") +
    geom_point(aes(x = Knots, y = p.num, col = Name), shape = 3) +
    scale_color_manual(values = c(`Boundary knots` = "Orange",Knots = "blue")) +
    geom_boxplot(data = data.frame(y = max(Knots.df$p.num) +1 ,  x = grid),
                 aes(x, y),
                 width = 0.4, col = "black") +
    scale_y_continuous(breaks = c(unique(Knots.df$p.num),max(Knots.df$p.num) +1),
                       labels = c(paste("p = ", unique(Knots.df$p.num)),
                                  "Time grid \n distribution"))+
    theme(axis.text=element_text(size = 15, hjust = 0.5),
          axis.text.x=element_text(angle=+90),
          axis.title=element_text(size=18,face="bold"),
          axis.line = element_line(colour="black"),
          plot.title=element_text(size=20, face="bold", vjust=1, lineheight=0.6),
          legend.text=element_text(size=14),
          legend.position="bottom",
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.key.size = unit(.9, "cm"),
          legend.key.width = unit(.9,"cm"),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          plot.margin=unit(c(10,5,5,5),"mm"), 
          strip.background = element_blank(),
          strip.text.x = element_blank() )+
    labs(x = "Time grid", y = "" )
  
  dataplot <- data$Dataset
  GrowthCurve <- ggplot(data= dataplot, aes(x=Time, y=Vol,group=ID) )+
    scale_x_continuous(breaks = as.integer(seq(min(grid),max(grid),length.out = length(grid)/5) ) ) +
    geom_line() +
    geom_point()+
    labs(x="",y="")+
    theme(plot.title = element_text(hjust = 0.5),title =element_text(size=10, face='bold'))
 
  Knots.Plot.new <-arrangeGrob(rbind(ggplotGrob(GrowthCurve), ggplotGrob(Knots.Plot), size = "last"))
  Knots.Plot.new <-ggplotify::as.ggplot(Knots.Plot.new)

  ALLcrossvalid<-ldply(crossvalid, rbind)
  
  mean<-sapply(p.values, function(x){ mean(ALLcrossvalid$CrossLikelihood[ALLcrossvalid$p==x]) } )
  meandata<-data.frame(p=p.values,mean=mean)
  
  ValidationPlot<-ggplot()+
                  geom_line(data=meandata,aes(x=p,y=mean,linetype="mean",col="mean"),size=1.2)+
                  geom_line(data=ALLcrossvalid,aes(x=p,y=CrossLikelihood,group=sim,linetype="test",col="test"),size=1.1)+
                  geom_point(data=meandata,aes(x=p,y=mean),size=2 )+
                  labs(title="Cross-LogLikelihood Plot",y= "LogLikelihood ", x= "Natural Spline Dimension")+
                  scale_color_manual("",breaks=c("mean","test"),values = c("black","grey"))+
                  scale_linetype_manual("",breaks=c("mean","test"),values = c("solid","dashed"))+
                  theme(legend.title=element_blank(),
                        plot.title = element_text(hjust = 0.5))
  
  if(save==TRUE)
  {
    ggsave(filename="CrossLogLikePlot_Estimation_p.pdf",plot = ValidationPlot,width=29, height = 20, units = "cm",scale = 1,path=path )
  }
  
  
  return(list(CrossLogLikePlot=ValidationPlot,KnotsPlot = Knots.Plot.new))
}


CalcLikelihood<-function(p,data.funcit,TestSet,grid){
 
  perc<-max(TestSet[,1])
  
  #grid<-sort(unique(data.funcit[,3] ) )
  points<-data.funcit[,2]
  ID<-data.funcit[,1]
  timeindex<-data.funcit[,3]
  
  fcm.fit <- fitfclust(x=points,curve=ID,timeindex=timeindex,q=p,h=1,K=1,p=p,grid=grid,seed=2404)
  
  Gamma<-as.matrix(fcm.fit$par$Gamma)
  sigma<-fcm.fit$par$sigma
  Lambda<- fcm.fit$par$Lambda
  alpha<- fcm.fit$par$alpha
  mu<- fcm.fit$par$lambda.zero+Lambda*c(alpha)
  
  ##### make basis considering the test set
  
  # TimeGrid<- unique(sort(TestSet[,3]))
  # 
  # tempBase <-cbind(1, ns(TimeGrid, df = (p - 1)))
  # base <- svd(tempBase)$u
  
  base<-fcm.fit$FullS
  
  
  Likelihood<-function(x,base,TestSet)
  {
    data.temp<-TestSet[TestSet$ID==x,]
    time.temp<-data.temp$Time
    base.i<-base[grid%in%time.temp,]
    Yi<-data.temp$Vol
    n.i<-length(time.temp)
    Vi<-base.i%*%Gamma%*%t(base.i)+sigma^2*diag(n.i)
    mui<-base.i%*%mu
    invVi<-ginv(Vi)
    -n.i/2*log(det(Vi)) - n.i/2*log(2*pi) - 1/2*t(Yi-mui)%*%invVi%*%(Yi-mui)
  }
  
  Li<-sapply(1:perc,Likelihood,base,TestSet)
  Likelihood<-sum(Li)
  return( Likelihood )
}
