rm(list = ls())
source('MainFunction.R')
source('ZPYZ.R')
source('BN.R')
############################################################
# estimation for the number of factors
############################################################
dgpNs <- c(200,400) # the sample size
dgpM <- 3 # the number of groups
dgpR0 <- 2 # the number of factor in each group
dgpPs <- c(50,75,100) # the length of segments in each group
dgpdeltas <- c(0.3,0.7) # the strength of weak factor
tot.sim <- 1000
result <- NULL
for (nnn in 1:length(dgpNs)) {
  dgpN <- dgpNs[nnn]
  for (ppp in 1:length(dgpPs)) {
    dgpPm <- dgpPs[ppp]
    rvec <- rep(dgpR0,1+dgpM)
    for (ddd in 1:length(dgpdeltas)) {
      dgpdelta <- dgpdeltas[ddd]
      RRR <- matrix(rep(0,tot.sim*(dgpM+1)),nrow = tot.sim)
      for (ttt in 1:tot.sim) {
        res.Y <- dgpf(dgpN, dgpPm, dgpR0, dgpM, dgpdelta)
        Y <- res.Y$dY
        glabel <- res.Y$glabel
        RRR[ttt,] <- estinumfactor(Y,glabel)$rvechat
        if(ttt %% 10 ==  0 ){cat("...")}
        if(ttt %% 100 == 0 ){cat(ttt,"\n") }
      }
      ppprrr <- colMeans(RRR == rvec)
      pr0 <- ppprrr[1]; pr1 <- ppprrr[2]; pr2 <- ppprrr[3]
      pr3 <- ppprrr[4]
      result1 <- cbind(dgpN,dgpPm,dgpdelta,pr0,pr1,pr2,pr3)
      result <- rbind(result,result1)
      cat('T =', dgpN, ', p =', dgpPm, 'delta =', dgpdelta, '\n', sep=' ')
      cat('=> Prob = (', ppprrr,')', '\n', sep=' ')
    }
  }
}
write.csv(result,"Numfactor.csv")



############################################################
# estimation for the factors and loadings
############################################################
dgpNs <- c(200,400) # the sample size
dgpM <- 3 # the number of groups
dgpR0 <- 2 # the number of factor in each group
dgpPs <- c(50,75,100) # the length of segments in each group
dgpdeltas <- c(0.3,0.7) # the strength of weak factor
tot.sim <- 1000
result <- NULL
for (nnn in 1:length(dgpNs)) {
  dgpN <- dgpNs[nnn]
  for (ppp in 1:length(dgpPs)) {
    dgpPm <- dgpPs[ppp]
    dgpP <- dgpPm*dgpM # the dimension of observation
    rvec <- rep(dgpR0,1+dgpM) # the sequence of the number of factor
    
    for (ddd in 1:length(dgpdeltas)) {
      dgpdelta <- dgpdeltas[ddd]
      PAerror1 <- PGerror1 <- PBerror1 <- PFerror1 <- NULL
      PAerror2 <- PBerror2 <- NULL
      PGFerror <- NULL
      
      for (ttt in 1:tot.sim) {
        res.Y <- dgpf(dgpN, dgpPm, dgpR0, dgpM, dgpdelta)
        Y <- res.Y$dY; glabel <- res.Y$glabel
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
        
        
        # estimation of ZhangPanYaoZhou
        result2 <- LFMest(Y,5,c(dgpR0,dgpR0*(dgpM+1)))
        Ahat2 <- result2$A; PAhat2 <- result2$PA
        Bhat2 <- result1$B; PBhat2 <- result2$PB
        
        PAerror2[ttt] <- compError(dPA,PAhat2,type = "F")
        PBerror2[ttt] <- compError(dPAB,PBhat2,type = "F")
        
        
        # estimation of BaiNg
        result3 <- StandardPCA(Y,sum(rvec))
        PGFhat <- result3$PGhat
        
        PGFerror[ttt] <- compError(dPGF,PGFhat,type = "F")
        
        
        if(ttt %% 10 ==  0 ){cat("...")}
        if(ttt %% 100 == 0 ){cat(ttt,"\n") }
      }
      oPAmean <- mean(PAerror1); oPAsd <- sd(PAerror1)
      oPGmean <- mean(PGerror1); oPGsd <- sd(PGerror1)
      oPBmean <- mean(PBerror1); oPBsd <- sd(PBerror1)
      oPFmean <- mean(PFerror1); oPFsd <- sd(PFerror1)
      
      zpyzPAzmean <- mean(PAerror2); zpyzPAzsd <- sd(PAerror2)
      zpyzPBzmean <- mean(PBerror2); zpyzPBzsd <- sd(PBerror2)
      
      bnPGFmean <- mean(PGFerror); bnPGFsd <- sd(PGFerror)
      
      result0 <- cbind(dgpN,dgpPm,dgpdelta,
                       oPAmean,oPAsd,oPGmean,oPGsd,
                       oPBmean,oPBsd,oPFmean,oPFsd,
                       zpyzPAzmean,zpyzPAzsd,zpyzPBzmean,zpyzPBzsd,
                       bnPGFmean,bnPGFsd)
      result <- rbind(result,result0)
      write.csv(result,"FactorEstimation.csv")
    }
  }
}

############################################################
# plot
############################################################
dgpN <- 100 # the sample size
dgpM <- 3 # the number of groups
dgpR0 <- 2 # the number of factor in each group
dgpPm <- 100 # the length of segments in each group
dgpdeltas <- seq(0.1,0.9,0.1) # the strength of weak factor
tot.sim <- 500
result <- NULL

dgpP <- dgpPm*dgpM # the dimension of observation
rvec <- rep(dgpR0,1+dgpM) # the sequence of the number of factor

for (ddd in 1:length(dgpdeltas)) {
  dgpdelta <- dgpdeltas[ddd]
  PAerror1 <- PGerror1 <- PBerror1 <- PFerror1 <- NULL
  PAerror2 <- PBerror2 <- NULL
  PGFerror <- NULL
  
  for (ttt in 1:tot.sim) {
    res.Y <- dgpf(dgpN, dgpPm, dgpR0, dgpM, dgpdelta)
    Y <- res.Y$dY; glabel <- res.Y$glabel
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
    
    # estimation of ZhangPanYaoZhou
    result2 <- LFMest(Y,5,c(dgpR0,dgpR0*(dgpM+1)))
    Ahat2 <- result2$A; PAhat2 <- result2$PA
    Bhat2 <- result1$B; PBhat2 <- result2$PB
    
    PAerror2[ttt] <- compError(dPA,PAhat2,type = "F")
    PBerror2[ttt] <- compError(dPAB,PBhat2,type = "F")

    # estimation of BaiNg
    result3 <- StandardPCA(Y,sum(rvec))
    PGFhat <- result3$PGhat
    
    PGFerror[ttt] <- compError(dPGF,PGFhat,type = "F")
    
    if(ttt %% 10 ==  0 ){cat("...")}
    if(ttt %% 100 == 0 ){cat(ttt,"\n")}
    }
  
  oPAmean <- mean(PAerror1); oPGmean <- mean(PGerror1)
  oPBmean <- mean(PBerror1); oPFmean <- mean(PFerror1)
  
  zpyzPAzmean <- mean(PAerror2); zpyzPBzmean <- mean(PBerror2)

  bnPGFmean <- mean(PGFerror)
  
  result0 <- cbind(dgpdelta,
                   oPAmean,oPGmean,
                   oPBmean,oPFmean,
                   zpyzPAzmean,zpyzPBzmean,
                   bnPGFmean)
  result <- rbind(result,result0)
}

library(ggplot2)
ggplot(result,aes(x=dgpdelta)) + 
  geom_line(aes(y=oPGmean, color = "R(G)"),size = 0.8) +
  geom_line(aes(y=oPFmean, color = "R(F)"),size = 0.8) +
  geom_line(aes(y=bnPGFmean, color = "R(GF)"),size = 0.8)+
  scale_color_manual(name = "",
    values = c("R(G)" = "#FF0000",
               "R(F)" = "#008000",
               "R(GF)" = "#0000FF"),
    labels = c(expression({R^2}(widehat(G))), 
               expression({R^2}(widehat(F))), 
               expression({R^2}(widehat((G*","*F)))))
  ) +
  labs(x = expression(delta), y = "Error") +
  scale_x_continuous(limits=c(0.1,0.9), breaks=seq(0.1,0.9,0.1)) +
  theme(legend.position=c(0.2,0.9),legend.background=element_blank())+
  theme(legend.key.size=unit(12,'pt'))


ggplot(result,aes(x=dgpdelta)) + 
  geom_line(aes(y=oPAmean, color = "R(A)"),size = 0.8) +
  geom_line(aes(y=oPBmean, color = "R(B)"),size = 0.8) +
  geom_line(aes(y=zpyzPAzmean, color = "R(zA)"),size = 0.8)+
  geom_line(aes(y=zpyzPBzmean, color = "R(zB)"),size = 0.8)+
  scale_color_manual(
    name = "",
    values = c("R(A)" = "#FF0000",
               "R(B)" = "#008000",
               "R(zA)" = "#0000FF",
               "R(zB)" = "#800080"),
    labels = c(expression({R^2}(widehat(A))),
               expression({R^2}(widehat(B))), 
               expression({zR^2}(widehat(A))),
               expression({zR^2}(widehat(B))))
  ) +
  labs(x = expression(delta), y = "Error") +
  scale_x_continuous(limits=c(0.1,0.9), breaks=seq(0.1,0.9,0.1)) +
  theme(legend.position=c(0.2,0.9),legend.background=element_blank())+
  theme(legend.key.size=unit(10,'pt'))


