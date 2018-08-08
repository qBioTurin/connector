PlotWithinness.i <- function(ClustCurve,MeanCurves,i,shift=0)
{
  ### Returns a plot representing the i-th cluster withinness measures circles.
  
  dataplot <- DataFrameWithinness.i(ClustCurve,MeanCurves,i,shift=shift)
  ### Data frame for plot
  circles <- dataplot$circles
  WithDist <- dataplot$WithDist
  feature <- WithDist[,3]
  feature.name <- colnames(WithDist)[3]
  nfeature <- length(unique(ClustCurve[,feature.name]))[1]
  feature.palette <- rainbow(nfeature)
  ### Plot distances
  plot.i <- ggplot() + geom_circle(aes(x0=x0, y0=y0, r=r), data=circles,size=1)
  plot.i <- plot.i  + geom_point(data=WithDist,aes(x=x1,y=y1,col=feature),shape=(i-1)) + labs(title=paste("Cluster",ClustSymbol[i],"withinness"),x="distance",y="distance")
  plot.i <- plot.i + scale_colour_manual(values = feature.palette[unique(feature)],name=feature.name) + xlab("distance x") + ylab("distance y")
  return(plot.i)
}

#######################

BetweenCluster_CurvDist <- function(ClustCurve,i)
{
  
  K <- length(unique(ClustCurve[,4]))
  ClustSymbol <- cluster.symbol(K)
  ### i-th cluster info
  ### curves data
  ClustCurve.i <- ClustCurve[ClustCurve[,4]==i,]
  ### curves ID
  index.i <- sort(unique(ClustCurve.i[,1]))
  ### magnitudo
  ni <- length(index.i)
  
  ### other clusters info
  ### curves data
  others <- ClustCurve[-which(ClustCurve[,4]==i),]
  ### curves cluster and ID
  clusterID.others <- unique(others[,c(4,1)])
  clusterID.others[,1] <- ClustSymbol[clusterID.others[,1]]
  ### curves ID
  index.others <- sort(clusterID.others[,2])
  ### other clusters magnitudo
  n.others <- length(index.others)
  ### i-th cluster distance matrix
  dist.curve <- matrix(numeric(ni*n.others),ncol=ni)
  colnames(dist.curve) <- paste("curve",index.i,"distance")
  for(j in 1:ni)
  {
    A <- ClustCurve.i[ClustCurve.i[,1]==index.i[j],2:3]
    for(k in 1: n.others)
    {
      B <- others[others[,1]==index.others[k],2:3]
      dist.curve[k,j] <- hausdorff(A,B)
    }
  }
  
  dist.info <- cbind(clusterID.others,dist.curve)
  other.ncluster <- K-1
  other.cluster <- unique(clusterID.others[,1])
  
  ### Min between dist.curveance
  betw.m <- numeric(other.ncluster)
  ### Max between dist.curveance
  betw.M <- betw.m
  
  ### The nearest other cluster curve to i-th cluster
  near.curve <- betw.m
  ### The farest other cluster curve to i-th cluster
  far.curve <- betw.m
  
  ### Betweenness
  for (p in 1: other.ncluster)
  {
    betw.m[p] <- min(dist.info[dist.info[,1]==sort(other.cluster)[p],-c(1:2)])
    betw.M[p] <- max(dist.info[dist.info[,1]==sort(other.cluster)[p],-c(1:2)])
  }
  
  ### Nearest and farest curve to i-th cluster
  for (p in 1: other.ncluster)
  {
    near.curve[p] <- index.others[which(dist.curve==betw.m[p],arr.ind=T)[1]]
    far.curve[p] <- index.others[which(dist.curve==betw.M[p],arr.ind=T)[1]]
  }
  between.i <- cbind(betw.m,betw.M)
  rownames(between.i) <- paste("cluster",ClustSymbol[-i])
  colnames(between.i) <- c("betw min","betw max")
  return(list(CurveDistance=dist.info,Between=between.i,NearCurveID=near.curve,FarCurveID=far.curve))
  
  ## Returns a list with four arguments: (i) CurveDistance, a matrix composed by the hausdorff distance among the samples belonging to the i-th cluster, (ii) Between, a matrix for minimum and maximum Hausdorff distance of i-th cluster curves from the other clusters, (iii) NearCurveID and (iv) FarCurveID are the nearest and the farthest samples from the i-th cluster.
}

#######################

BetweenCluster_MeanDist <- function(ClustCurve,MeanCurves,i)
{
  ### Calculates the Hausdorff distance between the i-th cluster meancurve and the other clusters meancurves.
  
  K <- length(unique(ClustCurve[,4]))
  ClustSymbol <- cluster.symbol(K)
  TimeGrid <- sort(unique(ClustCurve[,2]))
  ### i-th cluster curves data
  ClustCurve.i <- ClustCurve[ClustCurve[,4]==i,]
  ### i-th cluster max obs time
  tmax <- max(ClustCurve.i[,2])
  ### i-th cluster meancurve truncated at tmax
  MeanCurve.i <- MeanCurves[TimeGrid <= tmax,i]
  ### i-th cluster obs time grid
  TimeGrid.i <- TimeGrid[TimeGrid <= tmax]
  ### i-th cluster meancurve
  A <- cbind(TimeGrid.i,MeanCurve.i)
  ### other clusters
  other.cluster <- sort(unique(ClustCurve[-which(ClustCurve[,4]==i),4]))
  
  ### Betweenness centroid dist.curveance
  betweencentroid.i <- matrix(numeric(length(other.cluster)),nrow=1)
  count <- 1
  for (j in other.cluster)
  {
    ClustCurve.j <- ClustCurve[ClustCurve[,4]==j,]
    tmax.j <- max(ClustCurve.j[,2])
    MeanCurve.j <- MeanCurves[TimeGrid <= tmax.j,j]
    TimeGrid.j <- TimeGrid[TimeGrid <= tmax.j]
    B <- cbind(TimeGrid.j,MeanCurve.j)
    betweencentroid.i[count] <- hausdorff(A,B)
    count <- count +1
  }
  between <- matrix(numeric(K),nrow=1)
  between[,-i] <- betweencentroid.i
  colnames(between) <- paste(ClustSymbol,"dist")
  return(BetweenCentr=between)
  ### Returns a numeric vector for i-th cluster meancurve hausdorff distance from other clusters meancurves.
}

#######################

Betweenness <- function(ClustCurve,MeanCurves)
{
  ### Computes the Betweenness across all clusters.
  
  K <- length(unique(ClustCurve[,4]))
  ClustSymbol <- cluster.symbol(K)
  ClassCurve <- unique(ClustCurve[,c(1,4)])[,2]
  classes <- ClustSymbol[ClassCurve]
  centroid.dist <- matrix(numeric(K*K),nrow=K)
  between <- matrix(numeric(K*2),ncol=2)
  colnames(between) <- c("Min distance","Nearest cluster")
  rownames(between) <- ClustSymbol
  rownames(centroid.dist) <- ClustSymbol
  colnames(centroid.dist) <- rownames(centroid.dist)
  
  for(i in 1:K)
  {
    centroid.dist[i,] <-  BetweenCluster_MeanDist(ClustCurve,MeanCurves,i)
    ClustSymbol.i <- ClustSymbol[-i]
    min.dist <- min(centroid.dist[i,][centroid.dist[i,]!=0])
    argmin.dist <- ClustSymbol.i[which.min(centroid.dist[i,][centroid.dist[i,]!=0])]
    between[i,] <- c(min.dist,argmin.dist)
  }
  return(list(Betweenness=between,CentroidDist=centroid.dist,Classes=classes))
  
  ### Returns a list with three arguments: (i) Betweenness, a matrix composed by two columns, the minimum cluster distance and the cluster name that achieves the minimum, (ii) CentroidDist, a K x K matrix (where K is the cluster number) containing centroids cluster distance among each cluster couple, (iii) Classes, a vector containing the cluster name membership for each dataset sample.
}

##########

DataFrameWithinness.i <- function(ClustCurve,MeanCurves,i,ClustSymbol,shift=0)
{
  ### Creates a dataframe for i-th cluster withinness measures.
  
  K <- length(unique(ClustCurve[,4]))
  
  ### Centroid withinness distance
  Withinness.i <- WithCluster_MeanDist(ClustCurve,MeanCurves,i)
  
  ### Mean and standard deviation distance
  mean.dist <- mean(Withinness.i)
  st.dev <- sd(Withinness.i)
  
  ### Data frame for ggplot
  circles.i <- data.frame(
    x0 <- rep(0+shift,3),
    y0 <- rep(0,3),
    r <- rep(mean.dist, 3)+c(st.dev,0,-st.dev),
    distance <- c("sd","mean","sd"),
    Cluster <- factor(ClustSymbol[i],levels=ClustSymbol)
  )
  colnames(circles.i) <- c("x0","y0","r","distance","Cluster")
  
  WithDist.i <- data.frame(
    x1 <- Withinness.i + shift,
    y1 <- numeric(length(Withinness.i)),
    Cluster <- factor((i-1),levels=c(0:(K-1)))
  )
  colnames(WithDist.i) <- c("x1","y1","Cluster")
  return(DataFrame.i=list(circles=circles.i,WithDist=WithDist.i))
  
  ### Returns a list composed by two data frames: (i) "circles" containing coordinates and radius for i-th cluster withinness circles measures, (ii) "WithDist" containing the i-th cluster curves withinness distance from i-th cluster centroid/meancurve.
}

####################

WithCluster_MeanDist <- function(ClustCurve,MeanCurves,i)
{
  #' Computes the Hausdorff distance between the curves belonging to the i-th cluster and the corresponding mean curve.
  
  TimeGrid <- sort(unique(ClustCurve[,2]))
  ### i-th cluster curves data
  ClustCurve.i <- ClustCurve[ClustCurve[,4]==i,]
  ### i-th cluster curves ID
  index <- sort(unique(ClustCurve.i[,1]))
  ### i-th cluster magnitudo
  ni <- length(index)
  withmean.i <- numeric(ni)
  ### i-th cluster max obs time
  tmax <- max(ClustCurve.i[,2])
  ### i-th cluster meancurve truncated at tmax
  MeanCurve.i <- MeanCurves[,i][which(TimeGrid <= tmax)]
  ### i-th cluster obs time grid
  TimeGrid.i <- TimeGrid[TimeGrid <= tmax]
  ### i-th cluster meancurve
  A <- matrix(c(TimeGrid.i,MeanCurve.i),ncol=2)
  for (k in 1:ni)
  {
    ### i-th cluster curves
    B <- ClustCurve.i[ClustCurve.i[,1]==index[k],2:3]
    ### Hausdorff distance between i-th cluster meancurve and k-th curve respectively
    withmean.i[k] <- hausdorff(A,B)
  }
  return(withmean.i)
  ### Returns a vector reporting the hausdorff distance of the i-th cluster curves from the corresponding cluster meancurve.
}

###############

WithCluster_CurvDist <- function(ClustCurve,i)
{
  ###### Computes the Hausdorff distance among curves belonging to the i-th cluster.
  
  ### i-th cluster curves data
  ClustCurve.i <- ClustCurve[ClustCurve[,4]==i,]
  ### i-th cluster curves ID
  index <- sort(unique(ClustCurve.i[,1]))
  ### i-th cluster magnitudo
  ni <- length(index)
  
  if(ni!=1)
  {
    within.i <- numeric((1/2)*(factorial(ni)/(factorial(ni-2))))
    count <- 0
    for (k in 1:(ni-1))
    {
      A <- ClustCurve.i[ClustCurve.i[,1]==index[k],2:3]
      for (j in (k+1):ni)
      {
        B <- ClustCurve.i[ClustCurve.i[,1]==index[j],2:3]
        count <- count +1
        within.i[count] <- hausdorff(A,B)
      }
    }
  }
  else within.i <- 0
  ### i-th cluster withinness
  
  # Returns a vector reporting the hausdorff distance among the  curves belonging to the i-th cluster.
  return(within.i)
}

##################

Withinness <- function(ClustCurve,MeanCurves,centroids=TRUE)
{
  ### Compute the withinness across all clusters obtained.
  
  K <- length(unique(ClustCurve[,4]))
  ClustSymbol <- cluster.symbol(K)
  ### Withinness matrix
  withinness <- matrix(numeric(K*4),ncol=K)
  rownames(withinness) <- c("mean","sd","min","max")
  colnames(withinness) <- paste("Cluster ",ClustSymbol,sep="")
  
  for (i in 1:K)
  {
    if(centroids==FALSE) ### i-th cluster curves distance
    {
      within.i <- WithCluster_CurvDist(ClustCurve,i)
    }
    else                 ### i-th cluster curves and meancurve distance
      within.i <- WithCluster_MeanDist(ClustCurve,MeanCurves,i)
    
    ### i-th cluster mean distance
    MeanDist <- mean(within.i)
    ### i-th cluster standard deviation distance
    if(length(within.i)==1) { StDev <- 0 }
    else                        StDev <- sd(within.i)
    MinDist <- min(within.i)
    MaxDist <- max(within.i)
    ### i-th cluster withinness data
    withinness[,i] <- cbind(MeanDist,StDev,MinDist,MaxDist)
  }
  #### Returns a matrix with four columns: mean, standard deviation, minimum and maximum withinness distance across the clusters.
  return(withinness)
}

##################

