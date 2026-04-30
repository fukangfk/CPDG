rm(list = ls())
library(Matrix)

# Estimation the number of the global and group pseudo-factors
# INPUT:
#  @Y: the observation matrix (T-by-p)
# OUTPUT:
#  @R: the eigenvalue-ratio statistics
#  @rvechat: the sequence of the number of factors (2-dim vector)
estinumfactor <- function(Y){
  dT <- nrow(Y); dp <- ncol(Y)
  
  pmax <- ceiling(dp/4)
  S.AutoCov <- t(Y)%*%Y/dT
  aa0 <- eigen(S.AutoCov)$value
  aa1 <- aa0[1:pmax]/aa0[2:(pmax+1)]
  bb <- searchpeak(aa1)
  r0 <- bb$r0; r1 <- bb$r1
  
  R <- cbind(0,aa1)
  rvechat <- c(r0,r1-r0)
  
  return(list(R=R, rvechat=rvechat))
}

# calculate the locations of two largest peaks in the vector
# INPUT:
#  @x: a numeric sequence
# OUTPUT:
#  @r0: the index corresponding to the largest value in x
#  @r1: the index corresponding to the second largest value in x
searchpeak <- function(x){
  Ratio00 <- c(1,x)
  Rmax <- length(x)
  Ratio00a <- splus2R::peaks(Ratio00)
  Ratio00p0 <- Ratio00[Ratio00a==TRUE]
  peaknum <- length(Ratio00p0)
  Ratio00peakloc <- rep(0,peaknum)
  j <- 0
  for(i in 1:(Rmax+1)){
    if(Ratio00a[i] == TRUE){
      j <- j+1
      Ratio00peakloc[j] <- i
    }
  }
  r0i <- which.max(Ratio00p0)
  Ratio00p0a <- Ratio00p0[-r0i]
  r1j <- which.max(Ratio00p0a)
  if(r0i == r1j){
    r0esti <- Ratio00peakloc[r0i]-1
    r1esti <- Ratio00peakloc[r1j+1]-1
    
  }
  if(r0i>r1j){
    r0esti <- Ratio00peakloc[r1j]-1
    r1esti <- Ratio00peakloc[r0i]-1
  }
  if(r0i<r1j){
    r0esti <- Ratio00peakloc[r0i]-1
    r1esti <- Ratio00peakloc[r1j+1]-1
  }
  list(r0=r0esti,r1=r1esti)
}

# calculate the eigen-space with top r eigen-vectors of matrix A
# INPUT:
#  @A: a numeric matrix (n-by-n)
#  @r: a positive number
# OUTPUT:
#  @A: a matrix with top r eigen-vectors of matrix A  (n-by-r)
eigenspace <- function(A,r){
  eig.A <- eigen(A)
  M <- eig.A$vectors[,1:r]
  return(M)
}

# PCA estimator for the global and group pseudo-factors and loading
# INPUT:
#  @Y: the observation matrix (T-by-p)
#  @rvec: the sequence of the number of factors (2-dim vector)
# OUTPUT:
#  @Ahat: estimated global loading matrix (p-by-r0)
#  @Ghat: estimated global pseudo-factor matrix (T-by-r0)
#  @Bhat: estimated group factor loading matrix (p-by-r)
#  @Fhat: estimated group pseudo-factor matrix (T-by-r)
gwest <- function(Y, rvec = NULL){
  dT <- nrow(Y); dp <- ncol(Y)
  if(is.null(rvec)){
    R.result <- estinumfactor(Y,glabel)$rvechat
    r0 <- R.result[1]; r1 <- R.result[-1]
  } else {
    r0 <- rvec[1]; r1 <- rvec[-1]
  }
  
  M1 <- ((Y)%*%t(Y))/(dT*dp)
  Ghat <- sqrt(dT)*eigenspace(M1,r0) # the estimated global pseudo-factors
  Ahat <- t(Y)%*%Ghat/dT # the estimated global loading
  

  dYY <- Y-Ghat%*%t(Ghat)%*%Y/dT
  M2 <- ((dYY)%*%t(dYY))/(dT*dp)
  Fhat <- sqrt(dT)*eigenspace(M2,r1) # the estimated group pseudo-factors
  Bhat <- t(dYY)%*%Fhat/dT # the estimated group loading
  
  return(list(Ghat = Ghat, Ahat = Ahat, Fhat = Fhat, Bhat = Bhat))
}

# calculate the value of quasi-likelihood given the pre-specific break-pair
# INPUT:
#  @Y: the observation matrix (T-by-p)
#  @rvec: the sequence of the number of factors (2-dim vector)
#  @Ghat: estimated global pseudo-factors (T-by-r0)
#  @Fhat: estimated group pseudo-factors (T-by-r)
#  @Ks: the pre-specific break-pair (2-dim vector)
# OUTPUT:
#  @value1: the value of global QL function
#  @value2: the value of group QL function
QMLEV <- function(Y,rvec,Ghat,Fhat,Ks){
  N <- nrow(Y); p <- ncol(Y)
  r0 <- rvec[1]; r1 <- sum(rvec)
  L <- length(Ks)/2
  KKSa <- c(0,Ks[1],N); KKSb <- c(0,Ks[2],N)
  
  SigmaG <- array(0,dim = c(r0,r0,L+1))
  SigmaF <- array(0,dim = c(r1-r0,r1-r0,L+1))
  
  value1 <- value2 <- 0
  for (ll in 1:(L+1)) {
    la1 <- KKSa[ll]+1; la2 <- KKSa[ll+1]
    lb1 <- KKSb[ll]+1; lb2 <- KKSb[ll+1]
    Deltaa <- (la2-la1); Deltab <- (lb2-lb1)
    
    SigmaG[,,ll] <- (t(Ghat[la1:la2,])%*%Ghat[la1:la2,])/Deltaa
    SigmaF[,,ll] <- (t(Fhat[lb1:lb2,])%*%Fhat[lb1:lb2,])/Deltab
    value1 <- value1 - Deltaa*log(det(as.matrix(SigmaG[,,ll])))
    value2 <- value2 - Deltab*log(det(as.matrix(SigmaF[,,ll])))
  }
  return(list(value1=value1,value2=value2))
}

# estimation for global and group breakpoints
# INPUT:
#  @Y: the observation matrix (T-by-p)
#  @rvec: the sequence of the number of factors (2-dim vector)
# OUTPUT:
#  @Khat1: the estimated global breakpoint
#  @Khat2: the estimated group breakpoint
QMLE.CP <- function(Y,rvec = NULL){
  N <- nrow(Y); p <- ncol(Y)
  
  if(is.null(rvec)){
    rvec <- estinumfactor(Y,glabel)$rvechat
  }
  
  Factor.result <- gwest(Y,rvec)
  Ghat <- as.matrix(Factor.result$Ghat)
  Fhat <- as.matrix(Factor.result$Fhat)
  
  li.value1 <- li.value2 <- rep(-Inf,N)
  for (tt in (N*0.25):(N*0.75)) {
    Ks <- c(tt,tt)
    res.K <- QMLEV(Y,rvec,Ghat,Fhat,Ks)
    li.value1[tt] <- res.K$value1
    li.value2[tt] <- res.K$value2
  }
  
  
  # li.value[is.na(li.value)] <- -Inf
  Khat1 <- which(li.value1==max(li.value1), arr.ind = TRUE)
  Khat2 <- which(li.value2==max(li.value2), arr.ind = TRUE)
  return(list(Khat=c(Khat1,Khat2)))
}
