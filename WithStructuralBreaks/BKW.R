# the Code for Baltagi et al.(2017, 2021)
# Data generated procedure function
# generate date from the general vector factor model with structural breaks
# INPUT:
#  @dgpN: the number of observation (scalar)
#  @dgpP: the dimensional of observation (scalar)
#  @dgpR: the number of the factor (scalar)
#  @dgpK: a sequence of breakpoints (L-dim vector)
# OUTPUT:
#  @dY: the observation matrix (T-by-p)
#  @dA: factor loading array (p-by-r0-by-2)
#  @dX: factor matrix (T-by-r)
#  @dPA: the projection matrix of factor loading matrix (p-by-p)
dgpf0 <- function(dgpN, dgpP, dgpR, dgpK){
  dgpL <- length(dgpK)
  
  # Y -- observation term (T x p)
  dgpY <- matrix(rep(0,dgpN*dgpP), ncol = dgpP)
  # E -- idiosyncratic term (T x p)
  dgpE <- matrix(rep(0, dgpN*dgpP), ncol = dgpP) 
  # X -- common factors (T x r)
  dgpX <- matrix(rep(0, dgpN*dgpR), ncol = dgpR)
  # (A_1,A_2) -- pre- and post-break loading (p x r)
  dgpA <- array(rep(0,dgpP*dgpR*(dgpL+1)), dim = c(dgpP,dgpR,(dgpL+1)))
  # (PA_1,PA_2) -- projection matrix for A_1 and A_2 
  dgpAproj <- array(rep(0,dgpP*dgpP*(dgpL+1)), dim = c(dgpP,dgpP,(dgpL+1)))
  
  ZeroOne <- c(-1, 1)
  # generate X
  for(j in 1:dgpR){ 
    a <- runif(1, 0.4, 0.95)*sample(ZeroOne,1)
    dgpX[,j] <- arima.sim(model=list(ar=a), n=dgpN)
    dgpX[,j] <- dgpX[,j]/sd(dgpX[,j])*runif(1,1,2)
  }
  
  # generate E
  for(j in 1:dgpP){ 
    a <- runif(1, 0.4, 0.95)*sample(ZeroOne,1)
    dgpE[,j] <- arima.sim(model=list(ma=a), n=dgpN, innov=rnorm(dgpN, 0, 0.5))
  }
  
  # (0,k_0,T)
  dgpKK <- c(0,dgpK,dgpN)
  
  for (l in 1:(dgpL+1)) {
    # generate A_l
    dgpA[,,l] <- matrix(runif(dgpP*dgpR, -1, 1), ncol = dgpR)
    
    l1 <- dgpKK[l]+1; l2 <- dgpKK[l+1]
    dgpY0 <- (dgpX[l1:l2,])%*%t(dgpA[,,l]) + dgpE[l1:l2,]
    dgpdY0 <- matrix(0, nrow = nrow(dgpY0), ncol = ncol(dgpY0))
    for (i in 1:ncol(dgpY0)) {
      dgpdY0[,i] <- dgpY0[,i]-mean(dgpY0[,i])
    }
    dgpY[l1:l2,] <- dgpdY0
    
    # generate PA_l
    dgpAproj[,,l] <- dgpA[,,l]%*%(solve(t(dgpA[,,l])%*%dgpA[,,l]))%*%t(dgpA[,,l])
  }
  
  return(list(dY=dgpdY0, dA=dgpA, dX=dgpX, dPA=dgpAproj))
}

# estimation for factor and loading given the number of factor
# INPUT:
#  @Y: the observation matrix (T-by-p)
#  @r: the number of factors (scalar)
# OUTPUT:
#  @Ahat: estimated factor loading matrix (p-by-r)
#  @Ghat: estimated factor matrix (T-by-r)
estifactor0 <- function(Y,r){
  dN <- nrow(Y); dp <- ncol(Y)
  
  # estimate the factor and loading
  M1 <- (Y%*%t(Y))/(dN*dp)
  Ghat <- sqrt(dN)*eigenspace(M1,r)
  Ahat <- t(Y)%*%Ghat/dN
  
  return(list(Ghat=Ghat,Ahat=Ahat))
}

# estimation for the number of factor
# INPUT:
#  @Y: the observation matrix (T-by-p)
#  @kmax: a prior upper bound of r
estinumfactor0 <- function(Y, kmax, standardize = FALSE){
  dN <- nrow(Y); dp <- ncol(Y)
  
  if (standardize){Y <- scale(Y, center = TRUE, scale = TRUE)}
  v <- rep(0, kmax)
  kfactor <- 1:kmax
  bNT <- (dN+dp)/(dN*dp)
  cNT <- min(dN,dp)
  bev <- eigen((Y%*%t(Y))/dN)$values
  ER <- bev[1:kmax]/bev[2:(kmax+1)]
  
  for(k in 1:kmax){
    v[k] <- sum(bev[(k+1):dN])
  }
  PC1 <- v-v[kmax]*bNT*log(bNT)*kfactor
  PC2 <- v+v[kmax]*bNT*log(cNT)*kfactor
  PC3 <- v+v[kmax]*log(cNT)/cNT*kfactor
  
  IC1 <- log(v)-bNT*log(bNT)*kfactor
  IC2 <- log(v)+bNT*log(cNT)*kfactor
  IC3 <- log(v)+log(cNT)/cNT*kfactor
  
  AIC1 <- v+v[kmax]*(2/dN)*kfactor
  AIC2 <- v+v[kmax]*(2/dp)*kfactor
  AIC3 <- v+v[kmax]*(2*(dp+dN-kfactor)/(dp*dN))*kfactor
  
  BIC1 <- v+v[kmax]*(log(dp)/dp)*kfactor
  BIC2 <- v+v[kmax]*(log(dN)/dN)*kfactor
  BIC3 <- v+v[kmax]*((dp+dN-kfactor)*log(dp*dN)/(dp*dN))*kfactor
  
  return(list(ER=which.max(ER),
              PC1=which.min(PC1),
              PC2=which.min(PC2),
              PC3=which.min(PC3),
              IC1=which.min(IC1),
              IC2=which.min(IC2),
              IC3=which.min(IC3),
              AIC1=which.min(AIC1),
              AIC2=which.min(AIC2),
              AIC3=which.min(AIC3),
              BIC1=which.min(BIC1),
              BIC2=which.min(BIC2),
              BIC3=which.min(BIC3)))
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

# calculate the LS values given the pre-specific breakpoints
# INPUT:
#  @Y: the observation matrix (T-by-p)
#  @r: the number of factors (scalar)
#  @plfG: estimated factors (T-by-r)
#  @Ks: pre-specific breakpoints (L-dim vector)
# OUTPUT:
#  @value: the value of LS function
BKWLS <- function(Y,r,plfG,Ks){
  N <- nrow(Y); p <- ncol(Y)
  
  L <- length(Ks)
  KKS <- c(0,Ks,N)
  SigmaG <- array(0,dim = c(r,r,L+1))
  
  value <- 0
  for (ll in 1:(L+1)) {
    l1 <- KKS[ll]+1; l2 <- KKS[ll+1]
    Delta <- (l2-l1)
    SigmaG[,,ll] <- (t(plfG[l1:l2,])%*%plfG[l1:l2,])/Delta
    for (ii in l1:l2) {
      aaa <- plfG[ii,]%o%plfG[ii,]-SigmaG[,,ll]
      value <- value + sum(aaa^2)
    }
  }
  return(value)
}

# estimation for breakpoints
# INPUT:
#  @Y: the observation matrix (T-by-p)
#  @L: the number of break dates (1 or 2)
# OUTPUT:
#  @Khat: the estimated breakpoints
BKWLS.CP <- function(Y,L){
  N <- nrow(Y); p <- ncol(Y)
  R.result <- estinumfactor0(Y,ceiling(p/4))
  r <- R.result$IC1
  Factor.result <- estifactor0(Y,r)
  Ghat <- as.matrix(Factor.result$Ghat)
  
  if(L==1){
    li.value <- rep(Inf,N)
    for (tt in (N*0.25):(N*0.75)) {
      Ks <- c(tt)
      li.value[tt] <- BKWLS(Y,r,Ghat,Ks)
    }
  }
  
  if(L==2){
    li.value <- matrix(rep(Inf,N*N),ncol = N)
    for (tt in (N*0.25):(N*0.75-(0.15*N))) {
      for (pp in (tt+(0.15*N)):(N*0.75)) {
        Ks <- c(tt,pp)
        li.value[tt,pp] <- BKWLS(Y,r,Ghat,Ks)
      }
    }
  }
  
  li.value[is.na(li.value)] <- Inf
  Khat <- which(li.value==min(li.value), arr.ind = TRUE)
  return(list(Khat=as.vector(Khat)))
}


