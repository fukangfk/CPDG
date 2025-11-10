rm(list = ls())
library(Matrix)
# Data generated procedure function
# generate date from model (1)
# INPUT:
#  @dgpN: the number of observation (scalar)
#  @dgpPd: the dimensional of the segment in each group (M-by-1)
#  @dgpR: the number of the pre- and post-break factors ((M+1)-by-2)
#  @dgpM: the number of group (scalar)
#  @dgpK: the global and group breakpoints (2-dim vector)
#  @dgpdelta: the strength of weak factor (scalar (0,1))
# OUTPUT:
#  @dY: the observation matrix (T-by-p)
#  @glabel: the clustering label of the observation (p-dim vector)
#  @dA: global factor loading array (p-by-r0-by-2)
#  @dX: global factor matrix (T-by-r0)
#  @dB: group factor loading array (p-by-r-by-2)
#  @dZ: group factor matrix (T-by-r)
dgpf <- function(dgpN, dgpPd, dgpR, dgpM, dgpK, dgpdelta){
  # (r_0,r_1,...r_M)
  dgpRR <- apply(dgpR, 1, max)
  # r_0; (r_1,...r_M)
  dgpR0 <- dgpRR[1]; dgpRd <- dgpRR[-1]
  # p
  dgpP <- sum(dgpPd)
  # r:=r_1+...+r_M
  dgpRs <- sum(dgpRd)
  # single breakpoint
  dgpL <- length(dgpK)/2
  # k_0; k_0^g
  dgpKa <- dgpK[1]; dgpKb <- dgpK[2]
  
  # Y -- observation term (T x p)
  dgpY <- matrix(rep(0,dgpN*dgpP), ncol = dgpP)
  # E -- idiosyncratic term (T x p)
  dgpE <- matrix(rep(0, dgpN*dgpP), ncol = dgpP) 
  # X -- global factors (T x r_0)
  dgpX <- matrix(rep(0, dgpN*dgpR0), ncol = dgpR0)
  # Z -- group factors (T x r)
  dgpZ <- matrix(rep(0, dgpN*dgpRs), ncol = dgpRs)
  # (B_1^*,B_2^*) -- group pre- & post-break loading (p x r)
  dgpB <- array(rep(0,dgpP*dgpRs*(dgpL+1)), dim = c(dgpP,dgpRs,(dgpL+1)))
  # (A_1, A_2) -- global pre- & post-break loading (p x r_0)
  dgpA <- array(rep(0,dgpP*dgpR0*(dgpL+1)), dim = c(dgpP,dgpR0,(dgpL+1)))
  
  # grouped label
  glabel <- rep(1:dgpM, times = dgpPd)
  
  
  ZeroOne <- c(-1, 1)
  # generate X
  for(j in 1:dgpR0){ 
    a <- runif(1, 0.4, 0.95)*sample(ZeroOne,1)
    dgpX[,j] <- arima.sim(model=list(ar=a), n=dgpN)
    dgpX[,j] <- dgpX[,j]/sd(dgpX[,j])*runif(1,1,2)
  }
  # generate Z
  for(j in 1:dgpRs){ 
    a <- runif(1, 0.4, 0.95)*sample(ZeroOne,1)
    dgpZ[,j] <- arima.sim(model=list(ma=a), n=dgpN)
    dgpZ[,j] <- dgpZ[,j]/sd(dgpZ[,j])*runif(1,1,2)
  }
  
  # generate E
  for(j in 1:dgpP){ 
    a <- runif(1, 0.4, 0.95)*sample(ZeroOne,1)
    dgpE[,j] <- arima.sim(model=list(ma=a), n=dgpN, innov=rnorm(dgpN, 0, 0.5))
  }
  
  # (0,k_0,T); (0,k_0^g,T)
  dgpKKa <- c(0,dgpKa,dgpN); dgpKKb <- c(0,dgpKb,dgpN)
  
  # global and group signals
  dgpxi0 <- dgpeta0 <- matrix(rep(0, dgpN*dgpP), ncol = dgpP) 
  for (l in 1:(dgpL+1)) {
    # generate A_l
    dgpA[,,l] <- cbind(matrix(runif(dgpP*dgpR[1,l], -1, 1), ncol = dgpR[1,l]),
                       matrix(0,nrow = dgpP,ncol = dgpRR[1]-dgpR[1,l]))
    # generate B_l^*
    i1 <- 1; i2 <- dgpPd[1]
    j1 <- 1; j2 <- dgpRd[1]
    for(j in 1:dgpM) { 
      if(j > 1){
        i1 <- i1+dgpPd[j-1]; i2 <- i2+dgpPd[j]
        j1 <- j1+dgpRd[j-1]; j2 <- j2+dgpRd[j]
      }
      dgpB[i1:i2,j1:j2,l] <- cbind(dgpPd[j]^(-dgpdelta/2)*matrix(runif(dgpPd[j]*dgpR[j+1,l], -1, 1), ncol=dgpR[j+1,l]),
                                   matrix(0,nrow = dgpPd[j],ncol = dgpRR[j+1]-dgpR[j+1,l]))
    }
    
    la1 <- dgpKKa[l]+1; la2 <- dgpKKa[l+1]
    dgpxi0[la1:la2,] <- (dgpX[la1:la2,])%*%t(dgpA[,,l])
    lb1 <- dgpKKb[l]+1; lb2 <- dgpKKb[l+1]
    dgpeta0[lb1:lb2,] <- (dgpZ[lb1:lb2,])%*%t(dgpB[,,l])
  }
  
  dgpY0 <- dgpxi0 + dgpeta0 + dgpE
  dgpdY0 <- matrix(0, nrow = nrow(dgpY0), ncol = ncol(dgpY0))
  for (i in 1:ncol(dgpY0)) {
    dgpdY0[,i] <- dgpY0[,i]-mean(dgpY0[,i])
  }
  
  return(list(dY=dgpdY0, glabel=glabel,
              dA=dgpA, dX=dgpX, 
              dB=dgpB, dZ=dgpZ))
}

# Data generated procedure function (loading rotation change)
# generate date from model (1)
# INPUT:
#  @dgpN: the number of observation (scalar)
#  @dgpPd: the dimensional of the segment in each group (M-by-1)
#  @dgpR: the number of the pre- and post-break factors ((M+1)-by-2)
#  @dgpM: the number of group (scalar)
#  @dgpK: the global and group breakpoints (2-dim vector)
#  @dgpdelta: the strength of weak factor (scalar (0,1))
#  @dgpC: the rotation for the global loading
#  @dgpD: the rotation for the group loading
# OUTPUT:
#  @dY: the observation matrix (T-by-p)
#  @glabel: the clustering label of the observation vector (p-dim vector)
#  @dA: global factor loading array (p-by-r0-by-2)
#  @dX: global factor matrix (T-by-r0)
#  @dB: group factor loading array (p-by-r-by-2)
#  @dZ: group factor matrix (T-by-r)
dgpf2 <- function(dgpN, dgpPd, dgpR, dgpM, dgpK, dgpdelta, dgpC, dgpD){
  # (r_0,r_1,...r_M)
  dgpRR <- apply(dgpR, 1, max)
  # r_0; (r_1,...r_M)
  dgpR0 <- dgpRR[1]; dgpRd <- dgpRR[-1]
  # p
  dgpP <- sum(dgpPd)
  # r:=r_1+...+r_M
  dgpRs <- sum(dgpRd)
  # single breakpoint
  dgpL <- length(dgpK)/2
  # k_0; k_0^g
  dgpKa <- dgpK[1]; dgpKb <- dgpK[2]
  
  # Y -- observation term (T x p)
  dgpY <- matrix(rep(0,dgpN*dgpP), ncol = dgpP)
  # E -- idiosyncratic term (T x p)
  dgpE <- matrix(rep(0, dgpN*dgpP), ncol = dgpP) 
  # X -- global factors (T x r_0)
  dgpX <- matrix(rep(0, dgpN*dgpR0), ncol = dgpR0)
  # Z -- group factors (T x r)
  dgpZ <- matrix(rep(0, dgpN*dgpRs), ncol = dgpRs)
  # (B_1^*,B_2^*) -- group pre- & post-break loading (p x r)
  dgpB <- array(rep(0,dgpP*dgpRs*(dgpL+1)), dim = c(dgpP,dgpRs,(dgpL+1)))
  # (A_1, A_2) -- global pre- & post-break loading (p x r_0)
  dgpA <- array(rep(0,dgpP*dgpR0*(dgpL+1)), dim = c(dgpP,dgpR0,(dgpL+1)))
  
  
  # grouped label
  glabel <- rep(1:dgpM, times = dgpP)
  
  ZeroOne <- c(-1, 1)
  # generate X
  for(j in 1:dgpR0){ 
    a <- runif(1, 0.4, 0.95)*sample(ZeroOne,1)
    dgpX[,j] <- arima.sim(model=list(ar=a), n=dgpN)
    dgpX[,j] <- dgpX[,j]/sd(dgpX[,j])*runif(1,1,2)
  }
  # generate Z
  for(j in 1:dgpRs){ 
    a <- runif(1, 0.4, 0.95)*sample(ZeroOne,1)
    dgpZ[,j] <- arima.sim(model=list(ma=a), n=dgpN)
    dgpZ[,j] <- dgpZ[,j]/sd(dgpZ[,j])*runif(1,1,2)
  }
  
  # generate E
  for(j in 1:dgpP){ 
    a <- runif(1, 0.4, 0.95)*sample(ZeroOne,1)
    dgpE[,j] <- arima.sim(model=list(ma=a), n=dgpN, innov=rnorm(dgpN, 0, 0.5))
  }
  
  # (0,k_0,T); (0,k_0^g,T)
  dgpKKa <- c(0,dgpKa,dgpN); dgpKKb <- c(0,dgpKb,dgpN)
  
  # global and group signals
  dgpxi0 <- dgpeta0 <- matrix(rep(0, dgpN*dgpP), ncol = dgpP) 
  
  for (l in 1:(dgpL+1)) {
    if(l == 1){
      # generate A_1
      dgpA[,,l] <- cbind(matrix(runif(dgpP*dgpR[1,l], -1, 1), ncol = dgpR[1,l]),
                         matrix(0,nrow = dgpP,ncol = dgpRR[1]-dgpR[1,l]))
      # generate B_1^*
      i1 <- 1; i2 <- dgpPd[1]
      j1 <- 1; j2 <- dgpRd[1]
      for(j in 1:dgpM) { 
        if(j > 1){
          i1 <- i1+dgpPd[j-1]; i2 <- i2+dgpPd[j]
          j1 <- j1+dgpRd[j-1]; j2 <- j2+dgpRd[j]
        }
        dgpB[i1:i2,j1:j2,l] <- cbind(dgpPd[j]^(-dgpdelta/2)*matrix(runif(dgpPd[j]*dgpR[j+1,l], -1, 1), ncol=dgpR[j+1,l]),
                                     matrix(0,nrow = dgpPd[j],ncol = dgpRR[j+1]-dgpR[j+1,l]))
      }
    } else {
      # generate A_2
      dgpA[,,l] <- dgpA[,,1] %*% dgpC
      # generate B_2^*
      dgpB[,,l] <- dgpB[,,1] %*% as.matrix(bdiag(dgpD,dgpD,dgpD))
    }
    
    la1 <- dgpKKa[l]+1; la2 <- dgpKKa[l+1]
    dgpxi0[la1:la2,] <- (dgpX[la1:la2,])%*%t(dgpA[,,l])
    lb1 <- dgpKKb[l]+1; lb2 <- dgpKKb[l+1]
    dgpeta0[lb1:lb2,] <- (dgpZ[lb1:lb2,])%*%t(dgpB[,,l])
  }
  
  dgpY0 <- dgpxi0 + dgpeta0 + dgpE
  dgpdY0 <- matrix(0, nrow = nrow(dgpY0), ncol = ncol(dgpY0))
  for (i in 1:ncol(dgpY0)) {
    dgpdY0[,i] <- dgpY0[,i]-mean(dgpY0[,i])
  }
  
  return(list(dY=dgpdY0, glabel=glabel, 
              dA=dgpA, dX=dgpX, 
              dB=dgpB, dZ=dgpZ))
}

# estimation the number of the global and group pseudo-factors
# INPUT:
#  @Y: the observation matrix (T-by-p)
#  @glabel: the clustering label of the observation vector
# OUTPUT:
#  @R: the eigenvalue-ratio statistics
#  @rvechat: the sequence of the number of factors ((M+1)-dim vector)
#  @rtotal: the number of total factor by removing 
#           the observation of m-th group ((M+1)-dim vector)
estinumfactor <- function(Y,glabel){
  dT <- nrow(Y); dp <- ncol(Y)
  dM <- max(glabel) 
  dpm <- as.vector(table(glabel))
  
  rmm <- numeric(dM)
  pmax <- ceiling(dp/4)
  S.AutoCov <- t(Y)%*%Y/dT
  aa0 <- eigen(S.AutoCov)$value
  aa1 <- aa0[1:pmax]/aa0[2:(pmax+1)]
  bb <- searchpeak(aa1)
  r0 <- bb$r0; r1 <- bb$r1
  
  R <- cbind(0,aa1)
  for (tt in 1:dM) {
    index <- which(glabel == tt)
    pmax <- ceiling(sum(dpm[-tt])/4)
    S.AutoCov <- t(Y[,-c(index)])%*%Y[,-c(index)]/dT
    aa0 <- eigen(S.AutoCov)$value
    aa1 <- aa0[1:pmax]/aa0[2:(pmax+1)]
    bb <- searchpeak(aa1)
    rmm[tt] <- bb$r1
    aaa <- cbind(tt,aa1)
    R <- rbind(R,aaa)
  }
  
  rvechat <- rep(0,dM+1)
  rvechat[1] <- r0
  rvechat[2:(dM+1)] <- r1 - rmm
  rtotal <- c(r1,rmm)
  
  return(list(R=R, rvechat=rvechat, rtotal=rtotal))
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
#  @glabel: the clustering label of the observation vector (p-dim vector)
#  @rvec: the sequence of the number of factors ((M+1)-dim vector)
# OUTPUT:
#  @Ahat: estimated global loading matrix (p-by-r0)
#  @Ghat: estimated global pseudo-factor matrix (T-by-r0)
#  @Bhat: estimated group factor loading matrix (p-by-r)
#  @Fhat: estimated group pseudo-factor matrix (T-by-r)
gwest <- function(Y, glabel, rvec = NULL){
  dT <- nrow(Y); dp <- ncol(Y)
  if(is.null(rvec)){
    R.result <- estinumfactor(Y,glabel)$rvechat
    r0 <- R.result[1]; r1 <- R.result[-1]
  } else {
    r0 <- rvec[1]; r1 <- rvec[-1]
  }
  dM <- max(glabel)
  dpm <- as.vector(table(glabel)) 
  rbar <- sum(r1)
  
  
  M1 <- ((Y)%*%t(Y))/(dT*dp)
  Ghat <- sqrt(dT)*eigenspace(M1,r0) # the estimated global pseudo-factors
  Ahat <- t(Y)%*%Ghat/dT # the estimated global loading
  
  
  dYY <- Y-Ghat%*%t(Ghat)%*%Y/dT
  Bhatm <- vector(mode = "list", length = dM)
  PBhatm <- vector(mode = "list", length = dM)
  Fhatm <- vector(mode = "list", length = dM)
  Fhat <- matrix(0,nrow = dT, ncol = rbar)
  for (m in 1:dM) {
    index <- which(glabel == m)
    dYYm <- dYY[,index]
    M2m <- ((dYYm)%*%t(dYYm))/(dT*dpm[m])
    # the estimated m-th group pseudo-factors
    Fhatm[[m]] <- sqrt(dT)*eigenspace(M2m,r1[m])
    # the estimated m-th loading
    Bhatm[[m]] <- t(dYYm)%*%Fhatm[[m]]/dT
  }
  # the estimated group pseudo-factors
  Fhat <- as.matrix(do.call(cbind, Fhatm))
  # the estimated group loading
  Bhat <- as.matrix(Matrix(bdiag(Bhatm),sparse = FALSE))
  
  
  return(list(Ghat = Ghat, Ahat = Ahat, Fhat = Fhat, Bhat = Bhat))
}

# calculate the value of quasi-likelihood given the pre-specific break-pair
# INPUT:
#  @Y: the observation matrix (T-by-p)
#  @rvec: the sequence of the number of factors ((M+1)-dim vector)
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
#  @glabel: the clustering label of the observation vector (p-dim vector)
#  @rvec: the sequence of the number of factors ((M+1)-dim vector)
# OUTPUT:
#  @Khat1: the estimated global breakpoint
#  @Khat2: the estimated group breakpoint
QMLE.CP <- function(Y,glabel,rvec = NULL){
  N <- nrow(Y); p <- ncol(Y); dM <- max(glabel)
  
  if(is.null(rvec)){
    rvec <- estinumfactor(Y,glabel)$rvechat
  }
  
  Factor.result <- gwest(Y,glabel,rvec)
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
