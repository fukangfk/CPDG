rm(list = ls())
source('MainFunction.R')
source('ZPYZ.R')
source('BN.R')
############################################################
# estimation for the number of factors
############################################################
dgpN <- 400 # the sample size
dgpM <- 3 # the number of groups
dgpR0 <- 2 # the number of factor in each group
dgpPm <- 100 # the length of segments in each group
dgpdelta <- 0.5 # the strength of weak factor
dgpP <- dgpPm*dgpM # the dimension of observation
glabel <- rep(1:dgpM, each = dgpPm) # the clustering label
rvec <- rep(dgpR0,1+dgpM) # the sequence of the number of factor
tot.sim <- 100
RhatM1 <- matrix(0,nrow = tot.sim,ncol = dgpM+1)
for (ttt in 1:tot.sim) {
  res.Y <- dgpf(dgpN, dgpPm, dgpR0, dgpM, dgpdelta)
  Y <- res.Y$dY
  RhatM1[ttt,] <- estinumfactor(Y,glabel)$rvechat
  if(ttt %% 10 ==  0 ){cat("...")}
  if(ttt %% 100 == 0 ){cat(" ", ttt,"\n")}
}
colMeans(RhatM1 == matrix(rep(rvec,tot.sim),nrow = tot.sim, byrow = TRUE))


############################################################
# estimation for the factors and loadings
############################################################
dgpN <- 300 # the sample size
dgpM <- 3 # the number of groups
dgpR0 <- 2 # the number of factor in each group
dgpPm <- 100 # the length of segments in each group
dgpP <- dgpPm*dgpM # the dimension of observation
dgpdelta <- 0.1 # the strength of weak factor
glabel <- rep(1:dgpM, each = dgpPm) # the clustering label
rvec <- rep(dgpR0,1+dgpM) # the sequence of the number of factor
tot.sim <- 100
PAerror1 <- PGerror1 <- PBerror1 <- PFerror1 <- NULL
PAerror2 <- PBerror2 <- NULL
PGFerror <- NULL
for (ttt in 1:tot.sim) {
  res.Y <- dgpf(dgpN, dgpPm, dgpR0, dgpM, dgpdelta)
  Y <- res.Y$dY
  dA <- res.Y$dA; dG <- res.Y$dG; dPA <- res.Y$dPA
  dB <- res.Y$dB; dF <- res.Y$dF; dPB <- res.Y$dPB
  dPAB <- res.Y$dPAB; dGF <- cbind(dG,dF)
  
  dPG <- dG%*%(solve(t(dG)%*%dG))%*%t(dG)
  dPF <- dF%*%(solve(t(dF)%*%dF))%*%t(dF)
  dPGF <- dGF%*%(solve(t(dGF)%*%dGF))%*%t(dGF)
  
  # our estimation
  result1 <- gwest(Y,glabel,rvec)
  Ahat1 <- result1$Ahat; Ghat1 <- result1$Ghat; PAhat1 <- result1$PAhat
  Bhat1 <- result1$Bhat; Fhat1 <- result1$Fhat; PBhat1 <- result1$PBhat
  
  PGhat <- Ghat1%*%(solve(t(Ghat1)%*%Ghat1))%*%t(Ghat1)
  PFhat <- Fhat1%*%(solve(t(Fhat1)%*%Fhat1))%*%t(Fhat1)
  
  PAerror1[ttt] <- compError(dPA,PAhat1,type = "F")
  PGerror1[ttt] <- compError(dPG,PGhat,type = "F")
  PBerror1[ttt] <- compError(dPB,PBhat1,type = "F")
  PFerror1[ttt] <- compError(dPF,PFhat,type = "F")

  ################################################
  # estimation of ZhangPanYaoZhou
  result2 <- LFMest(Y,5,c(dgpR0,dgpR0*(dgpM+1)))
  Ahat2 <- result2$A; PAhat2 <- result2$PA
  Bhat2 <- result1$B; PBhat2 <- result2$PB
  
  PAerror2[ttt] <- compError(dPA,PAhat2,type = "F")
  PBerror2[ttt] <- compError(dPAB,PBhat2,type = "F")
  
  ################################################
  # estimation of BaiNg
  result3 <- StandardPCA(Y,sum(rvec))
  PGFhat  <- result3$PGhat
  
  PGFerror[ttt] <- compError(dPGF,PGFhat,type = "F")

  
  if(ttt %% 10 ==  0 ){cat("...")}
  if(ttt %% 100 == 0 ){cat(" ", ttt,"\n")}
}
cbind(mean(PAerror1),sd(PAerror1))
cbind(mean(PBerror1),sd(PBerror1))
cbind(mean(PGerror1),sd(PGerror1))
cbind(mean(PFerror1),sd(PFerror1))

cbind(mean(PAerror2),sd(PAerror2))
cbind(mean(PBerror2),sd(PBerror2))

cbind(mean(PGFerror),sd(PGFerror))




