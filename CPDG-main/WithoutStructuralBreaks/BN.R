StandardPCA <- function(Y,r){
  dT <- nrow(Y); dp <- ncol(Y)
  
  M1 <- ((Y)%*%t(Y))/(dT*dp)
  Ghat <- sqrt(dT)*eigenspace(M1,r) # factor
  Ahat <- t(Y)%*%Ghat/dT # loading
  PGhat <- Ghat%*%(solve(t(Ghat)%*%Ghat))%*%t(Ghat)
  PAhat <- Ahat%*%(solve(t(Ahat)%*%Ahat))%*%t(Ahat)
  Chat <- Ghat%*%t(Ahat)
  
  return(list(Ghat = Ghat,Ahat = Ahat,
              PAhat = PAhat,PGhat = PGhat,
              Chat = Chat))
}