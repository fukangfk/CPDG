# The estimation algorithm of Zhang et al.(2024)

# calculate \hat{M} and its eigen-decomposition 
# INPUT:
#  @Y: the observation matrix (N-by-p)
#  @autoupper: maximum of lag k
# OUTPUT:
#  @dWsum: the matrix \hat{M}
#  @value: eigenvalue of \hat{M}
#  @vector: eigenvector of \hat{M}
autoeigensum <- function(Y,autokupper){
  autoN <- nrow(Y); autop <- ncol(Y)
  dWsum <- matrix(rep(0,autop*autop),autop,autop)
  for(autok in 0:autokupper){
    dW <- cov(Y[(autok+1):autoN,],Y[1:(autoN-autok),])
    dWW <- dW%*%t(dW)
    dWsum <- dWsum+dWW
  }
  eiWW <- eigen(dWsum)
  return(list(dWsum=dWsum, value=eiWW$value,
              vector=eiWW$vector))
}

# calculate the auto-covariance matrix \hat{\Sigma}_{y}(k) and its eigen-decomposition 
# INPUT:
#  @Y: the observation matrix (N-by-p)
#  @autok: lag k
# OUTPUT:
#  @dWW: the matrix \hat{\Sigma}_{y}(k)
#  @value: eigenvalue of \hat{\Sigma}_{y}(k)
#  @vector: eigenvector of \hat{\Sigma}_{y}(k)
autoeigen <- function(Y,autok){
  autoN <- nrow(Y); autop <- ncol(Y)
  dW <- cov(Y[(autok+1):autoN,],Y[1:(autoN-autok),])
  dWW <- dW%*%t(dW)
  eiWW <- eigen(dWW)
  return(list(dWW=dWW, value=eiWW$value,
              vector=eiWW$vector))
}

# estimation of the number of factors
# INPUT:
#  @Y: the observation matrix (N-by-p)
#  @autokupper: maximum of lag k
#  @pmax: prior upper bound of total number of factors (pmax = p/4 default)
# OUTPUT:
#  @R: the eigenvalue-ratio 
#  @r0: the number of common factors
#  @r1: the total number of factors
estifactor <- function(Y,autokupper,pmax=NULL){
  dYN <- nrow(Y); dYp <- nrow(Y)
  if(is.null(pmax)){pmax <- ceiling(dYp/4)}
  aa0 <- autoeigen(Y,0)$value
  for(autok in 1:autokupper){
    aa <- autoeigen(Y,autok)$value
    aa0 <- aa0+aa
  }
  aa1 <- aa0[1:pmax]/aa0[2:(pmax+1)]
  bb <- searchpeak(aa1)
  return(list(R=aa1, r0=bb$r0, r1=bb$r1))
}

# main function
# estimation of the factor and loading
# INPUT:
#  @Y: the observation matrix (N-by-p)
#  @autokupper: maximum of lag k
#  @rvec: the number of common factors and total number of factors (2-dim vector)
# OUTPUT:
#  @R: the eigenvalue-ratio 
#  @r0: the number of common factors
#  @r1: the total number of factors
#  @A: estimation for the loadings for common factor
#  @PA: the projection for A
#  @B: estimation for the loadings cluster-specific factor
#  @PB: the projection for B
LFMest <- function(Y,autokupper,rvec=NULL){
  n <- nrow(Y); p <- ncol(Y)
  pmax <- ceiling(p/4)
  
  # step 1: estimation for the number of factors
  if(is.null(rvec)){
    r.result <- estifactor(Y,autokupper)
    r0 <- r.result$r0; r1 <- r.result$r1
  } else {
    r0 <- rvec[1]; r1 <- rvec[2]
  }
  
  
  # step 2: estimation for the loadings of common factor
  Mm <- autoeigensum(Y,autokupper)
  Ahat <- Mm$vector[,1:r0]
  Ahatproj <- Ahat%*%t(Ahat)
  
  # step 3: estimation for the loadings of cluster-specific factor
  dYY <- Y-Y%*%Ahat%*%t(Ahat)
  Bhat <- autoeigensum(dYY,autokupper)$vector[,1:(r1-r0)]
  Bhatproj <- Bhat%*%t(Bhat)
  
  return(list(r0=r0, A=Ahat, PA=Ahatproj,
              r1=r1, B=Bhat, PB=Bhatproj))
}

# Identification for the components not belonging to any clusters
# INPUT:
#  @B: the cluster-specific loading matrix (p-by-r)
#  @thresholdtype: the type of the threshold (1,2, or 3)
# OUTPUT:
#  @cluster0: a vector recording the components 
#             belonging to one group (0) or no group(-1)
#  @ind.cluster0: a vector recording the locations of the components 
#                 belonging to no group
#  @ind.cluster1: a vector recording the locations of the components 
#                 belonging to one group
#  @numcluster0: the number of components belonging to no group
#  @numcluster1: the number of components belonging to one group
nocluster <- function(B,thresholdtype){
  dp <- nrow(B); dr <- ncol(B)
  switch (thresholdtype,
    threshold <- dr/(dp*log(dp)*log(dp)),
    threshold <- dr/(dp*log(dp)),
    threshold <- dr/(dp*log(log(dp)))
  )
  
  bi <- apply(B, 1, function(x){norm(x,"2")^2})
  
  cluster0 <- rep(-1,dp)
  for(i in 1:dp) {
    if(bi[i] > threshold) {cluster0[i] <- 0}
  }
  
  ind.cluster0 <- which(cluster0 == -1)
  ind.cluster1 <- which(cluster0 == 0)
  numcluster0 <- length(ind.cluster0)
  numcluster1 <- length(ind.cluster1)
  
  list(cluster0=cluster0,
       ind.cluster0=ind.cluster0,ind.cluster1=ind.cluster1,
       numcluster0=numcluster0,numcluster1=numcluster1)
}

# Estimate the number of group
# INPUT:
#  @B: the cluster-specific loading matrix (p-by-r)
#  @threshold: the threshold (1-1/log(dgpN))
# OUTPUT:
#  @numgroup: the estimator for the number of groups
estMf <- function(B,threshold){
  aBB <- abs(B%*%t(B))
  eig.aBB <- eigen(aBB)
  ev <- eig.aBB$value
  numgroup <- length(ev[ev>threshold])
  return(numgroup)
}

# K-means Clustering
# INPUT:
#  @B0: the cluster-specific loading matrix (p0-by-r)
#         by the rows with their indices that belongs to no group
#  @eM: the number of groups
# OUTPUT:
#  @kmz: group-member vector
clusterppp <- function(B0,eM){
  dp <- nrow(B0); dr <- ncol(B0)
  #B1 <- matrix(0,dp,dp)
  #for(i in 1:dp) {
  #  for(j in 1:dp) {
  #    B1[i,j] <- abs(sum(B0[i,]*B0[j,]))/
  #      sqrt(sum(B0[i,]^2))/sqrt(sum(B0[j,]^2))
  #  }
  #}
  BB0 <- B0%*%t(B0)
  B1 <- diag((diag(BB0))^(-1/2))%*%abs(BB0)%*%diag((diag(BB0))^(-1/2))
  #kmz <- kmeans(B1, centers=eM)$cluster

  # a method for choosing the initial centers
  bbb <- rep(0,dp)
  for(i in 1:dp) {bbb[i] <- sum(B1[i,]^2)}
  Inival <- rep(0,eM)
  Inival[1] <- which.max(bbb)
  for(i in 2:eM){
    bbb <- rep(0,dp)
    ccc <- matrix(0,dp,i-1)
    for(jj in 1:dp){
      for(ii in 1:(i-1)){
        ccc[jj,ii] <- sum((B1[jj,]-B1[Inival[ii],])^2)
      }
      bbb[jj] <- min(ccc[jj,]) 
    }
    Inival[i] <- which.max(bbb)
  }
  kmz <- kmeans(B1,B1[Inival,])$cluster
  return(kmz)
}

# Calculate cluster error
# INPUT:
#  @Emv: a vector of estimated cluster membership
#  @Tmv: the true membership vector
clusterror <- function(Emv, Tmv) {
  n <- length(Emv)
  K <- length(as.matrix(table(Tmv)))
  Kpermn <- t(nptest::permn(K))
  L <- matrix(0, nrow = nrow(Kpermn), ncol = n)
  L[1:nrow(Kpermn), 1:n] <- Kpermn[1:nrow(Kpermn),Emv[1:n]]
  L <- unique(L) # delete duplicate rows
  MC <- NULL
  for (i in 1:nrow(L)) {
    MC[i] <- sum(as.numeric(Tmv != L[i,]))
  }
  return(min(MC)/n)
}

# Run the simulation
# INPUT:
#  @dN: the number of observation (scalar)
#  @dPm: the dimensional of the segment in each group (scalar)
#  @dRp: the ratio between the number of components belonging to no group 
#          and the number of components in one group
#  @dR0: the number of the factor in each group (scalar)
#  @dM: the number of group (scalar)
#  @iffac: whether to use real number of factors
# OUTPUT:
#  @PAerror: the errors of the estimators for A
#  @PBerror: the errors of the estimators for B
#  @E1error: the error E_1
#  @E2error: the error E_2
#  @Clerror: cluster error
#  @eM: the estimated number of clusters
simf <- function(dN,dPm,dRp,dR0,dM,iffac=TRUE){
  # generate data
  res.Y <- dgpf(dN,dPm,dRp,dR0,dM)
  Y <- res.Y$dY; PA <- res.Y$dPA; PAB <- res.Y$dPAB
  # Step1: estimation for the number of factors
  if(iffac){
    rvec <- c(dR0,dR0*(dM+1))
  } else {
    num.factor <- estifactor(Y,5)
    rvec <- c(num.factor$r0,num.factor$r1)
  }
  
  # Step 2&3: estimation for the loadings for common and cluster-specific factor
  res.AB <- LFMest(Y,5,rvec)
  Ahat <- res.AB$A; Bhat <- res.AB$B
  PAhat <- res.AB$PA; PBhat <- res.AB$PB
  
  PAerror <- comperror(PA,PAhat)
  PBerror <- comperror(PAB,PBhat)
  
  # Step 4: Identification for the components not belonging to any clusters
  E1error <- E2error <- NULL
  for (k in 1:3) {
    res.nogroup <- nocluster(Bhat,k)
    id.group <- res.nogroup$ind.cluster1
    id.nogroup <- res.nogroup$ind.cluster0
    num.group <- res.nogroup$numcluster1
    sss <- sum(id.group>dPm*(dM))
    E1error[k] <- (dPm*dM-num.group+sss)/(dPm*dM)
    E2error[k] <- sss/(dPm*dRp)
    if(k == 2){
      id.nogroup2 <- id.nogroup
      id.group2 <- id.group
      num.group2 <- num.group
    }
  }
  
  # Step 5: K-means clustering.
  tau <- 1-1/log(dN)
  dLabel <- c(rep(1:dM,each = dPm),rep(0,dPm*dRp))
  Mest <- estMf(Bhat,tau)
  B0hat <- Bhat[-id.nogroup2,]
  grouphat <- clusterppp(B0hat,Mest)
  labelhat <- numeric(dPm*(dM+dRp))
  for (i in 1:num.group2) {labelhat[id.group2[i]] <- grouphat[i]}
  id <- intersect(id.group2,which(dLabel!=0))
  Clerror <- clusterror(labelhat[id],dLabel[id])
  
  return(list(PAerror=PAerror,PBerror=PBerror,
              E1error=E1error,E2error=E2error,
              Clerror=Clerror,eM=Mest))
}
