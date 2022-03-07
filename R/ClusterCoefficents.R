#' Functional Davies–Bouldin (fDB) index and L2 (or with derivatives) distances
#' 
#' @description L2dist.curve2mu calculates the L2 (or with derivatives) distance between all the spline estimated for each sample with respect to the corresponding cluster center curve. L2dist.mu2mu calculates the L2 distance between the cluster centers curves. L2dist.mu20 calculates the L2 distance between the cluster centers curves and zero. Sclust.coeff calculates the S coefficents for each cluster. Rclust.coeff calculates the R coefficents for each cluster.
#' 
#' @param clust The vector reporting the cluster membership for each sample.
#' @param fcm.curve The list obtained applying the function fclust.curvepred to the fitfclust output, the FCM algorithm, more details in [Sugar and James].
#' @param database ...
#' @param fcm.fit The fitfclust output, the FCM algorithm, more details in [Sugar and James].
#' @param deriv The derivative values (0, 1 or 2).
#' 
#' @author Cordero Francesca, Pernice Simone, Sirovich Roberta
#'  
#' @return 
#'     The following indexes are calculated for obtaining a cluster separation measure:
#' \itemize{
#' \item{T:}{ the total tightness representing the dispersion measure given by
#' \deqn{ T = \sum_{k=1}^G \sum_{i=1}^n D_0(\hat{g}_i, \bar{g}^k); } (see \code{\link{L2dist.curve2mu}}) }
#'  \item{S_k:}{ 
#'  \deqn{ S_k = \sqrt{\frac{1}{G_k} \sum_{i=1}^{G_k} D_q^2(\hat{g}_i, \bar{g}^k);} }
#'  with G_k the number of curves in the k-th cluster; (see \code{\link{Sclust.coeff}})
#'  }
#'  \item{M_{hk}:}{ the distance between centroids (mean-curves) of h-th and k-th cluster 
#'  \deqn{M_{hk} =  D_q(\bar{g}^h, \bar{g}^k);  } (see\code{\link{L2dist.mu2mu}})
#'  }
#'  \item{R_{hk}:}{  a measure of how good the clustering is,
#'  \deqn{R_{hk} = \frac{S_h + S_k}{M_{hk}}; } (see\code{\link{Rclust.coeff}} )}
#'  \item{fDB_q:}{ functional Davies-Bouldin index, the cluster separation measure
#'  \deqn{fDB_q = \frac{1}{G} \sum_{k=1}^G \max_{h \neq k} {  \frac{S_h + S_k}{M_{hk}} }. }(see\code{\link{Rclust.coeff}}) }
#' }
#' Where the proximities measures choosen is defined as follow
#'  \deqn{D_q(f,g) = \sqrt( \integral | f^{(q)}(s)-g^{(q)}(s) |^2 ds ), d=0,1,2}
#' whit f and g are two curves and f^{(q)} and g^{(q)} are their q-th derivatives. Note that for q=0, the equation becomes the distance induced by the classical L^2-norm.
#' 
#'  @references
#' Gareth M. James and Catherine A. Sugar, (2003). Clustering for Sparsely Sampled Functional Data. Journal of the American Statistical Association.
#' 
#' @import sfsmisc dplyr 
#' @name DBindexL2dist
NULL
#> NULL

#' @rdname DBindexL2dist
L2dist.curve2mu <- function(clust,fcm.curve,database,fcm.fit=NULL,deriv=0,q){
  n.curves<-length(fcm.curve$gpred[,1])
  dist.curve2mu <-rep(0,n.curves)
  # Definition of the grid depending on the cluster membership
  cluster.grid<-ClusterGrid.truncation(clust,database,q)
  grid <-sort(unique(database$Time))
  G<-length(fcm.curve$meancurves[1,])
  
  Database.cl<-merge(data.frame(ID = unique(database$ID),Cluster=clust),
                     database,
                     by="ID")
  
  grid.cl <- lapply(1:length(cluster.grid),function(i)
    data.frame(Cluster =i,Time= cluster.grid[[paste0("G",i)]])
  )
  grid.cl <- do.call("rbind",grid.cl)
  Database.subset <- merge(Database.cl,grid.cl,by = c("Cluster","Time"),all.y= T)
  ShortCurves = Database.subset %>%
    dplyr::group_by(ID) %>%
    dplyr::count(ID) %>% filter(n < 3)
  
  if(length(ShortCurves$ID) > 0 ) 
    stop(paste0("The quantile values must be smaller in order to calculate correcty the fdb indexes. The following curves with the IDs:", 
                paste(ShortCurves$ID,collapse = ", "), " are too sparse and during the cutting procedure for calculating the indexes, these curves are too short." ))
  
  if(deriv==0)
  {
    dmeancurves <- fcm.curve$meancurves
    dgpred <- fcm.curve$gpred
  }
  else{
    if(is.null(fcm.fit)) warning("The fcm.fit is needed to calculate the derivatives!! ")
    dspl <-basis.derivation(fcm.fit,deriv)
    u.dspl<-dspl$u.dspl
    dmeancurves<-dspl$dmeancurves
    etapred <- fcm.curve$etapred
    
    matrix(0,n.curves,nrow(u.dspl)) -> dgpred
    dgpred <- t(sapply(1:n.curves ,function(ind) as.vector(u.dspl %*% etapred[ind,])))
  }
  
  meancurves <- data.frame(Time = rep(grid,G),
                           MeanCurve = c(dmeancurves),
                           Cluster = rep(1:G,each=length(grid)) )
  
  gpred <- data.frame(Time = rep(grid,n.curves),
                      gpred = c(t(dgpred)),
                      ID = rep(unique(database$ID),each=length(grid)), 
                      Cluster = rep(clust,each=length(grid)) )
  
  gpred.subset <- merge(gpred,Database.subset,by = c("Time","ID","Cluster"),all.y= T)
  Dataset.all<- merge(meancurves,gpred.subset,by = c("Time","Cluster"),all.y= T)
  
  fxk <- Dataset.all %>%
    group_by(ID , Time)
  
  fxk <- fxk %>% summarise(Diff = (gpred - MeanCurve)^2, Time = Time, ID = ID)
  
  
  dist.curve2mu <- fxk %>%
    group_by(ID) %>%
    group_map( ~ integrate.xy(x = .x$Time ,fx = .x$Diff ,use.spline = F))
  
  return(sqrt(unlist(dist.curve2mu)))
}

#' @rdname DBindexL2dist
L2dist.mu2mu <- function(clust,fcm.curve,database,fcm.fit=NULL,deriv=0,q){
  n.curves<-length(fcm.curve$gpred[,1])
  # gauss.info<-gauss.infoList[[n.curves+1]]
  # itimeindex<-gauss.info$itimeindex
  #  
  # weights <- gauss.info$gauss$weights 
  # a<-gauss.info$a
  # b<-gauss.info$b
  # Definition of the grid depending on the cluster membership
  cluster.grid<-ClusterGrid.truncation(clust,database,q) 
  grid <-sort(unique(database$Time))
  
  fcm.curve$meancurves->meancurves
  k<-length(meancurves[1,])
  dist.mu2mu<-matrix(0,ncol = k,nrow = k)
  
 #clust.ind <- combn(1:k,2)
  clust.ind <- expand.grid(1:k,1:k)
  
  if(deriv==0)
  {
    fxk <- (meancurves[,clust.ind$Var1]-meancurves[,clust.ind$Var2])^2
  }else{
    dspl <-basis.derivation(fcm.fit,deriv)
    u.dspl<-dspl$u.dspl
    dmeancurves<-dspl$dmeancurves
    
    fxk <- (dmeancurves[,clust.ind$Var1]-dmeancurves[,clust.ind$Var2])^2

  }
  # 
  # if(length(as.matrix(fxk)[1,] ) == 1 ){
  #   int <- (b-a)/2 * sum(weights * fxk)
  # }else int <- (b-a)/2 * colSums(weights * fxk)
  # 
  # dist.mu2mu[,] <- sqrt(int)
  int <- sapply( 1:length(fxk[1,]), function(i){
    cl.1 <- clust.ind[i,1]
    cl.2 <- clust.ind[i,2]
    # to claculate the diffeernce bertween the mean curves we consider the union of the two time intervals associated with each mean curve
    t.min<-min(min(cluster.grid[[cl.1]]),min(cluster.grid[[cl.2]]))
    t.max<-max(max(cluster.grid[[cl.1]]),max(cluster.grid[[cl.2]]))
    itime<- which(grid >= t.min & grid <= t.max )
    integrate.xy(grid[itime],fxk[itime,i],use.spline = F)
  } )
  dist.mu2mu[,] <- sqrt(int)
  return(dist.mu2mu)
}

#' @rdname DBindexL2dist
L2dist.mu20 <- function(clust,fcm.curve,database,fcm.fit=NULL,deriv=0,q){
  n.curves<-length(fcm.curve$gpred[,1])
  # gauss.info<-gauss.infoList[[n.curves+1]]
  # 
  # itimeindex <- gauss.info$itimeindex
  # weights <- gauss.info$gauss$weights 
  # a<-gauss.info$a
  # b<-gauss.info$b
  # 
  fcm.curve$meancurves->meancurves
  cluster.grid<-ClusterGrid.truncation(clust,database,q) 
  grid <-sort(unique(database$Time))
  k<-length(meancurves[1,])
  
  if(deriv==0)
  {
    fxk <- (meancurves)^2
  }else{
    dspl <-basis.derivation(fcm.fit,deriv)
    u.dspl<-dspl$u.dspl
    dmeancurves<-dspl$dmeancurves
    
    fxk <- (dmeancurves)^2
    
  }
  
  # if(length(as.matrix(fxk)[1,] ) == 1 ){
  #   int <- (b-a)/2 * sum(weights * fxk)
  # }else int <- (b-a)/2 * colSums(weights * fxk)
  # 
  int <- sapply( 1:length(fxk[1,]), function(i){
    cl.grid <- cluster.grid[[i]]
    itime <- which(grid %in% cl.grid)
    integrate.xy(grid[itime],fxk[itime,i],use.spline = F)
  })
  
  dist.mu20 <- sqrt(int)
  
  return(dist.mu20)
}

#' @rdname DBindexL2dist
L2dist.curve20 <- function(clust,fcm.curve,database,fcm.fit=NULL,deriv=0,q){
  n.curves<-length(fcm.curve$gpred[,1])
  dist.curve2mu <-rep(0,n.curves)
  grid <-sort(unique(database$Time))
  cluster.grid<-ClusterGrid.truncation(clust,database,q)
  grid <-sort(unique(database$Time))
  G<-length(fcm.curve$meancurves[1,])
  
  Database.cl<-merge(data.frame(ID = unique(database$ID),Cluster=clust),
                     database,
                     by="ID")
  if("Time" %in% colnames(Database.cl))
    colnames(Database.cl)[colnames(Database.cl) == "Time"] <- "Time"
  grid.cl <- lapply(1:length(cluster.grid),function(i)
    data.frame(Cluster =i,Time= cluster.grid[[paste0("G",i)]])
  )
  grid.cl <- do.call("rbind",grid.cl)
  
  Database.subset <- merge(Database.cl,grid.cl,by = c("Cluster","Time"),all.y= T)
  
  if(deriv==0)
  {
    dgpred <- fcm.curve$gpred
  }
  else{
    if(is.null(fcm.fit)) warning("The fcm.fit is needed to calculate the derivatives!! ")
    dspl <-basis.derivation(fcm.fit,deriv)
    u.dspl<-dspl$u.dspl
    etapred <- fcm.curve$etapred
    
    matrix(0,n.curves,nrow(u.dspl)) -> dgpred
    dgpred <- t(sapply(1:n.curves ,function(ind) as.vector(u.dspl %*% etapred[ind,])))
  }

  gpred <- data.frame(Time = rep(grid,n.curves),
                      gpred = c(t(dgpred)),
                      ID = rep(unique(database$ID),each=length(grid)), 
                      Cluster = rep(clust,each=length(grid)) )
  
  Dataset.all <- merge(gpred,Database.subset,by = c("Time","ID","Cluster"),all.y= T)

  fxk <- Dataset.all %>%
    group_by(ID , Time)
  
  fxk <- fxk %>% summarise(Diff = (gpred)^2, Time = Time, ID = ID)
  
  dist.curve2mu <- fxk %>%
    group_by(ID) %>%
    group_map( ~ integrate.xy(x = .x$Time ,fx = .x$Diff ,use.spline = F))
  
  return(sqrt(unlist(dist.curve2mu)))
}

#' @rdname DBindexL2dist
Sclust.coeff <- function(clust,fcm.curve,database,fcm.fit=NULL,deriv=0,q){

  distances <- L2dist.curve2mu(clust,fcm.curve,database,fcm.fit,deriv,q)
  k<-length(fcm.curve$meancurves[1,])

  out <- sapply(1:k,
   function(x){
    curve <- which(clust == x)
    if(length(curve)==0) warning("Empty clusters imply a NaN coefficents",immediate. = TRUE)
    sqrt(1/length(curve)*sum(distances[curve]^2))
   }
 )
  
  return(out)
}

#' @rdname DBindexL2dist
Rclust.coeff <- function(clust,fcm.curve,database,fcm.fit=NULL,deriv=0,q){
  k<-length(fcm.curve$meancurves[1,])
  
  emme <- L2dist.mu2mu(clust,fcm.curve,database,fcm.fit,deriv,q)
  cl   <- L2dist.mu20(clust,fcm.curve,database,fcm.fit,deriv,q)
  ######### Let name the cluster with A->Z from the lower mean curve to the higher.
 
  if( k !=1 )
  {
    Cl.order<-order(cl)
  }else{
    Cl.order<-1
  }
  
  symbols<-LETTERS[order(Cl.order)]
  row.names(emme)<-symbols
  colnames(emme)<-symbols
  
  esse <- Sclust.coeff(clust,fcm.curve,database,fcm.fit,deriv,q) 
  names(esse)<-symbols
  
  erre <- matrix(0,nrow=k,ncol=k,dimnames = list(symbols,symbols))

  S.matrix<-matrix(rep(esse,length(esse)),ncol=length(esse))
  sum.2S <- S.matrix + t(S.matrix) 
  erre <- sum.2S / emme
  diag(erre)<-0
  errei <- apply(erre,1,max)
  errebar <- mean(errei)
  
  return(list(erre=erre,errei=errei,fDB.index=errebar,emme=emme,esse=esse))
}


###############
# basis.derivation calculates the derivatives of the spline basis considering the matrix decomposition.
###############

basis.derivation<- function(fcm.fit,deriv)
{
  Lambda <- fcm.fit$par$Lambda
  alpha <- fcm.fit$par$alpha
  lambda.zero <- as.vector(fcm.fit$par$lambda.zero)
  Lambda.alpha <- lambda.zero + Lambda %*% t(alpha)
  q<-length(fcm.fit$FullS[1,])
  
  grid <- fcm.fit$grid
  
  cbind(1,ns(grid, df = (q - 1))) -> spl
  singular.spl <- svd(spl)
  cbind(0,ns.deriv(grid, df = (q-1), deriv = deriv)) -> dspl
  utrasf <- function(dspl){
    out <- matrix(0, nrow = dim(dspl)[1], ncol = dim(dspl)[2])
    for (colonna in 1:dim(dspl)[2]){
      out[,colonna] <- 1/singular.spl$d[colonna] * dspl %*% singular.spl$v[,colonna]
    }
    return(out)
  }
  u.dspl <- utrasf(dspl)
  dmeancurves <- u.dspl %*% Lambda.alpha
  return(list(u.dspl=u.dspl,dmeancurves=dmeancurves) )
}

#################
## ns.deriv is the ns function modified in order to make the function splineDesign calculates the derivatives!
#################

ns.deriv <- function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x), deriv) {
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax)) 
    x <- x[!nax]
  outside <- if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    (ol <- x < Boundary.knots[1L]) | (or <- x > Boundary.knots[2L])
  }
  else {
    if (length(x) == 1L) 
      Boundary.knots <- x * c(7, 9)/8
    FALSE
  }
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - 1L - intercept
    if (nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d", 
                       1L + intercept), domain = NA)
    }
    knots <- if (nIknots > 0L) {
      knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
                         2L)[-c(1L, nIknots + 2L)]
      quantile(x[!outside], knots)
    }
  }
  else nIknots <- length(knots)
  
  Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
  if (any(outside)) {
    warning(gettextf("I cannot derive"), domain = NA)
    #     basis <- array(0, c(length(x), nIknots + 4L))
    #     if (any(ol)) {
    #         k.pivot <- Boundary.knots[1L]
    #         xl <- cbind(1, x[ol] - k.pivot)
    #         tt <- splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0, 
    #             1))
    #         basis[ol, ] <- xl %*% tt
    #     }
    #     if (any(or)) {
    #         k.pivot <- Boundary.knots[2L]
    #         xr <- cbind(1, x[or] - k.pivot)
    #         tt <- splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0, 
    #             1))
    #         basis[or, ] <- xr %*% tt
    #     }
    #     if (any(inside <- !outside)) 
    #         basis[inside, ] <- splineDesign(Aknots, x[inside], 
    #             4)
  }
  else basis <- splineDesign(Aknots, x, ord = 4L, derivs = deriv)
  const <- splineDesign(Aknots, Boundary.knots, ord = 4L, derivs = c(2L, 
                                                                     2L))
  if (!intercept) {
    const <- const[, -1, drop = FALSE]
    basis <- basis[, -1, drop = FALSE]
  }
  qr.const <- qr(t(const))
  basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L), 
                                                     drop = FALSE])
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = 3L, knots = if (is.null(knots)) numeric() else knots, 
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("ns", "basis", "matrix")
  basis
}

####

ClusterGrid.truncation<-function(clust,database,q){
  grid.list<-lapply(unique(clust),function(j){
    id.curves <- which(clust == j)
    dbSubset= database[database$ID %in% id.curves, ]
    grid<- dbSubset$Time
    if(is.null(q)){ # Delete the outliers
      lw<- max(min(grid), as.numeric(quantile(grid, 0.25)) - (IQR(grid)*1.5)) #lower whisker
      up<- min(max(grid), as.numeric(quantile(grid, 0.75)) + (IQR(grid)*1.5)) # upper whisker
    }else{
      if(length(q)==1) q = c(q,1-q)
      q.up = max(q)
      q.lw = min(q)
      if(q.lw > 0.45) q.lw <- 0.45
      if(q.up < 0.55) q.up <- 0.55
      lw<- as.numeric(quantile(grid, q.lw))
      up<- as.numeric(quantile(grid, q.up))
    }
    return(sort(unique(grid[grid >= lw & grid <= up])))
    
  })
  names(grid.list) <- paste0("G",unique(clust))
  return(grid.list)
}

