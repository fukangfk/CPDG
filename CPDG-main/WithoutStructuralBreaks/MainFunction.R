library(Matrix)
# Data generated procedure function
# generate date from the grouped factor model
# INPUT:
#  @dgpN: the number of observation (scalar)
#  @dgpPm: the dimensional of the segment in each group (scalar)
#  @dgpR0: the number of the factor in each group (scalar)
#  @dgpM: the number of group (scalar)
#  @dgpDelta: the strength of weak factor (scalar (0,1))
# OUTPUT:
#  @dY: the observation matrix (T-by-p)
#  @dA: global factor loading matrix (p-by-r0)
#  @dG: global factor matrix (T-by-r0)
#  @dPA: the projection matrix of global factor loading matrix (p-by-p)
#  @dB: group factor loading matrix (p-by-r)
#  @dF: group factor matrix (T-by-r)
#  @dPB: the projection matrix of group factor loading matrix (p-by-p)
#  @dPAB: the matrix (I-A'A)B
#  @glabel: the clustering label of the observation vector (p-dim vector)
dgpf <- function(dgpN, dgpPm, dgpR0, dgpM, dgpdelta){
  # (p_1,...,p_M)
  dgpPd <- rep(dgpPm, dgpM)
  # (r_0,r_1,...,r_M)
  dgpRd <- rep(dgpR0, dgpM)
  # p
  dgpP <- dgpPm*dgpM
  # r:=r_1+...+r_M
  dgpR <- dgpR0*dgpM
  
  # Y -- observation term (T x p)
  dgpY <- matrix(rep(0,dgpN*dgpP), ncol = dgpP)
  # E -- idiosyncratic term (T x p)
  dgpE <- matrix(rep(0, dgpN*dgpP), ncol = dgpP) 
  # G  -- global factors (T x r_0)
  dgpG <- matrix(rep(0, dgpN*dgpR0), ncol = dgpR0)
  # F -- group factors (T x r)
  dgpF <- matrix(rep(0, dgpN*dgpR), ncol = dgpR)
  # B -- group loading (p x r)
  dgpB <- matrix(rep(0,dgpP*dgpR), ncol = dgpR)
  # A -- global loading (p x r_0)
  dgpA <- matrix(rep(0,dgpP*dgpR0), ncol = dgpR0)
  
  
  # the clustering label of the observation vector
  glabel <- rep(1:dgpM, each = dgpPm)
  
  ZeroOne <- c(-1, 1)
  # generate G
  for(j in 1:dgpR0){ 
    a <- runif(1,0.4,0.95)*sample(ZeroOne,1)
    dgpG[,j] <- arima.sim(model=list(ar=a), n=dgpN)
    dgpG[,j] <- dgpG[,j]/sd(dgpG[,j])*runif(1,1,2)
  }
  # generate F
  for(j in 1:dgpR){ 
    a <- runif(1,0.4,0.95)*sample(ZeroOne,1)
    dgpF[,j] <- arima.sim(model=list(ma=a), n=dgpN)
    dgpF[,j] <- dgpF[,j]/sd(dgpF[,j])*runif(1,1,2)
  }
  
  # generate E
  for(j in 1:dgpP){ 
    a <- runif(1,0.4,0.95)*sample(ZeroOne,1)
    dgpE[,j] <- arima.sim(model=list(ma=a), n=dgpN, innov=rnorm(dgpN, 0, 0.5))
  }
  
  # generate A
  dgpA <- matrix(runif(dgpP*dgpR0,-1,1), ncol = dgpR0)
  # generate B
  i1 <- 1; i2 <- dgpPd[1]
  j1 <- 1; j2 <- dgpRd[1]
  for(j in 1:dgpM) { 
    if(j > 1){
      i1 <- i1+dgpPd[j-1]; i2 <- i2+dgpPd[j]
      j1 <- j1+dgpRd[j-1]; j2 <- j2+dgpRd[j]
    }
    dgpB[i1:i2,j1:j2] <- dgpPm^(-dgpdelta/2)*matrix(runif(dgpPd[j]*dgpRd[j],-1,1), ncol=j2-j1+1)
  }
  
  dgpY0 <- (dgpG)%*%t(dgpA) + (dgpF)%*%t(dgpB) + dgpE
  dgpdY0 <- matrix(0, nrow = nrow(dgpY0), ncol = ncol(dgpY0))
  for (i in 1:ncol(dgpY0)) {
    dgpdY0[,i] <- dgpY0[,i]-mean(dgpY0[,i])
  }
  
  dgpAproj <- dgpA%*%(solve(t(dgpA)%*%dgpA))%*%t(dgpA)
  dgpBproj <- dgpB%*%(solve(t(dgpB)%*%dgpB))%*%t(dgpB)
  BAB0 <- dgpB-dgpAproj%*%dgpB
  dgpABproj <- BAB0%*%(solve(t(BAB0)%*%BAB0))%*%t(BAB0)
  
  
  return(list(dY=dgpdY0, dA=dgpA, dG=dgpG, dPA=dgpAproj,
              dB=dgpB, dF=dgpF, dPB=dgpBproj, dPAB=dgpABproj,
              glabel=glabel))
}

# calculate the error in sense of Frobenius or spectral norm
# INPUT:
#  @A1: a matrix or vector
#  @A2: a matrix or vector (the dimension of A2 is identical that of A1)
#  @type: the type of norm: Frobenius norm (default) or spectral norm
# OUTPUT:
#  @error: ||A1-A2||_F or ||A1-A2||_2
compError <- function(A1, A2, type = c("F","2")){
  E <- A1-A2
  if(is.null(type)){
    error <- norm(E, type = "F")
  }
  error <- norm(E, type = type)
  return(error)
}

# a multivariate measure of fit between G and Ghat
# INPUT:
#  @G: the true matrix
#  @Ghat: the estimated matrix
# OUTPUT:
#  @measure: result
traceratio <- function(G,Ghat){
  top <- t(G)%*%Ghat%*%(solve(t(Ghat)%*%Ghat))%*%t(Ghat)%*%G
  bottom <- t(G)%*%G
  measure <- psych::tr(top)/psych::tr(bottom)
  return(measure)
}

# estimation the number of the strong and weak factors
# INPUT:
#  @Y: the observation matrix (T-by-p)
#  @glabel: the clustering label of the observation vector
# OUTPUT:
#  @R: the eigenvalue-ratio 
#  @rvechat: the sequence of the number of factors ((M+1)-dim vector)
#  @rtotal: the total number of factor by deleting 
#           the observation of group m ((M+1)-dim vector)
estinumfactor <- function(Y,glabel){
  dT <- nrow(Y); dp <- ncol(Y)
  dM <- max(glabel); dpm <- as.vector(table(glabel))
  
  
  rmm <- numeric(dM)
  pmax <- ceiling(dp/4)
  S.AutoCov <- t(Y)%*%Y/(dT)
  aa0 <- eigen(S.AutoCov)$value
  aa1 <- aa0[1:pmax]/aa0[2:(pmax+1)]
  bb <- searchpeak(aa1)
  r0 <- bb$r0; r1 <- bb$r1
  
  R <- cbind(0,aa1)
  for (tt in 1:dM) {
    index <- which(glabel == tt)
    pmax <- ceiling(sum(dpm[-tt])/4)
    S.AutoCov <- t(Y[,-c(index)])%*%Y[,-c(index)]/(dT)
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
#  @M: a matrix with top r eigen-vectors of matrix A  (n-by-r)
eigenspace <- function(A,r){
  eig.A <- eigen(A)
  M <- eig.A$vectors[,1:r]
  return(M)
}

# PCA estimator for the global and group-wise factors and loading matrices
# INPUT:
#  @Y: the observation matrix (T-by-p)
#  @glabel: the clustering label of the observation vector (p-dim vector)
#  @rvec: the sequence of the number of factors ((M+1)-dim vector)
# OUTPUT:
#  @Ahat: estimated global factor loading matrix (p-by-r0)
#  @Ghat: estimated global factor matrix (T-by-r0)
#  @PAhat: estimated projection matrix of global factor loading matrix (p-by-p)
#  @Bhat: estimated group factor loading matrix (p-by-r)
#  @Fhat: estimated group factor matrix (T-by-r)
#  @PBhat: estimated projection matrix of group factor loading matrix (p-by-p)
gwest <- function(Y, glabel, rvec = NULL){
  dT <- nrow(Y); dp <- ncol(Y)
  if(is.null(rvec)){
    R.result <- estinumfactor(Y,glabel)$rvechat
    r0 <- R.result[1]; r1 <- R.result[-1]
  } else {
    r0 <- rvec[1] # the number of strong factors
    r1 <- rvec[-1] # the sequence of the number of weak factors
  }
  dM <- max(glabel) # the number of the group
  dpm <- as.vector(table(glabel)) # the length of each segment
  rbar <- sum(r1) # the total number of weak factors
  
  
  # estimate the global factor
  M1 <- ((Y)%*%t(Y))/(dT*dp)
  Ghat <- sqrt(dT)*eigenspace(M1,r0) # factor
  Ahat <- t(Y)%*%Ghat/dT # loading
  PAhat <- Ahat%*%(solve(t(Ahat)%*%Ahat))%*%t(Ahat) # projection matrix of loading
  
  # estimate the group-wise factor
  dYY <- Y-Ghat%*%t(Ghat)%*%Y/dT
  Bhatm <- vector(mode = "list", length = dM)
  PBhatm <- vector(mode = "list", length = dM)
  Fhatm <- vector(mode = "list", length = dM)
  Fhat <- matrix(0,nrow = dT, ncol = rbar)
  for (m in 1:dM) {
    index <- which(glabel == m)
    dYYm <- dYY[,index]
    M2m <- ((dYYm)%*%t(dYYm))/(dT*dpm[m])
    Fhatm[[m]] <- sqrt(dT)*eigenspace(M2m,r1[m])
    Bhatm[[m]] <- t(dYYm)%*%Fhatm[[m]]/dT
    PBhatm[[m]] <- Bhatm[[m]]%*%solve(t(Bhatm[[m]])%*%Bhatm[[m]])%*%t(Bhatm[[m]])
  }
  Bhat <- as.matrix(Matrix(bdiag(Bhatm),sparse = FALSE)) # loading
  Fhat <- as.matrix(do.call(cbind, Fhatm)) # factor
  PBhat <- as.matrix(Matrix(bdiag(PBhatm),sparse = FALSE)) # projection matrix of loading
  
  return(list(Ghat = Ghat, Ahat = Ahat, PAhat = PAhat,
              Fhat = Fhat, Bhat = Bhat, PBhat = PBhat))
}

