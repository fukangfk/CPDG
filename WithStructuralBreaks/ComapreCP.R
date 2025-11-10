# rm(list = ls())
source('MainFunction.R')
source('DBH.R')
source('BKW.R')
library(ggplot2)

dgpN <- 200
dgpM <- 3
dgpK <- as.integer(c(0.3,0.7)*dgpN)
dgpPm <- 50
dgpP <- rep(dgpPm, dgpM)
glabel <- rep(1:dgpM, each = dgpPm)
dgpdelta <- 0.6
tot.sim <- 200
scenario <- 1
Khats1 <- Khats4 <- Khats5 <- matrix(rep(0,2*tot.sim),ncol = 2)
Khats2 <- Khats3 <- numeric(tot.sim)
for (ttt in 1:tot.sim) {
  if (scenario == 1){
    dgpR <- cbind(c(2,2,2,2),c(2,2,2,2))
    res.Y <- dgpf(dgpN, dgpP, dgpR, dgpM, dgpK, dgpdelta)
    rvec <- rowSums(dgpR)
  } else {
    if(scenario == 2){
      dgpR <- cbind(c(2,2,2,2),c(2,2,2,2))
      dgpC <- matrix(c(1,0,2,1),nrow = 2, byrow = TRUE)
      dgpD <- matrix(c(1,2,0,1),nrow = 2, byrow = TRUE)
      res.Y <- dgpf2(dgpN, dgpPm, dgpR, dgpM, dgpK, dgpdelta, dgpC, dgpD)
      rvec <- dgpR[,1]
    } else {
      dgpR <- cbind(c(4,4,4,4),c(4,4,4,4))
      dgpC <- matrix(c(1,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0),nrow = 4, byrow = TRUE)
      dgpD <- matrix(c(1,2,0,0,0,1,0,0,0,0,0,0,0,0,0,0),nrow = 4, byrow = TRUE)
      res.Y <- dgpf2(dgpN, dgpPm, dgpR, dgpM, dgpK, dgpdelta, dgpC, dgpD)
      rvec <- dgpR[,1]
    }
  } 

  Y <- res.Y$dY
  CP.result1 <- QMLE.CP(Y,glabel,rvec)
  CP.result2 <- DBHQML.CP(Y,1)
  # CP.result3 <- BKWLS.CP(Y,1)
  # CP.result4 <- BKWLS.CP(Y,2)
  CP.result5 <- DBHQML.CP(Y,2)
  
  Khats1[ttt,] <- CP.result1$Khat
  Khats2[ttt] <- CP.result2$Khat
  # Khats3[ttt] <- CP.result3$Khat
  # Khats4[ttt,] <- CP.result4$Khat
  Khats5[ttt,] <- CP.result5$Khat
  if(ttt %% 10 ==  0 ){cat("...")}
  if(ttt %% 100 == 0 ){cat(ttt,"\n") }
}
Khats1
Khats2
Khats3
Khats4
Khats5

p <- paste("64(III)")

Khats1 <- as.data.frame(Khats1)
colnames(Khats1) <- c('k1','k2')

o <- paste0("GCPD",p,".png",sep="")
png(filename = o,width = 500,height = 500,
    units = "px",bg = "white",res = 100)
ggplot(Khats1,aes(y=..density..)) +
  geom_histogram(aes(x = k1, fill = 'k1'), binwidth = 1) +
  geom_histogram(aes(x = k2, fill = 'k2'), binwidth = 1) +
  xlab('T') + ylab('Frequency') + 
  scale_fill_manual(name = '', values = c('k1' = 'skyblue2', 'k2' = 'orangered1'),
                    labels = c(expression(widehat(k)),expression(widehat(k)^g))) + 
  scale_x_continuous(limits=c(60,140), breaks=seq(60,140,10)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) + 
  theme(legend.position=c(0.9,0.9),legend.background=element_blank())+
  theme(legend.key.size=unit(10,'pt'))
dev.off()

############################################################
# k_0=k_0^g
o <- paste0("GCPDc",p,".png",sep="")
png(filename = o,width = 500,height = 500,
    units = "px",bg = "white",res = 100)
ggplot(Khats1,aes(y=..density..)) +
  geom_histogram(aes(x = k1), binwidth = 1, fill = 'skyblue2') +
  xlab('T') + ylab('Frequency') + 
  scale_x_continuous(limits=c(60,140), breaks=seq(60,140,10)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2))
dev.off()

o <- paste0("GCPDg",p,".png",sep="")
png(filename = o,width = 500,height = 500,
    units = "px",bg = "white",res = 100)
ggplot(Khats1,aes(y=..density..)) +
  geom_histogram(aes(x = k2), binwidth = 1, fill = 'orangered1') +
  xlab('T') + ylab('Frequency') + 
  scale_x_continuous(limits=c(60,140), breaks=seq(60,140,10)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2))
dev.off()
############################################################


Khats2 <- as.data.frame(Khats2)
colnames(Khats2) <- c('k1')

o <- paste0("DBH",p,".png",sep="")
png(filename = o,width = 500,height = 500,
    units = "px",bg = "white",res = 100)
ggplot(Khats2,aes(y=..density..)) +
  geom_histogram(aes(x = k1),fill="skyblue2",colour="skyblue2",binwidth = 1) + 
  xlab('T') + ylab('Frequency') + 
  scale_x_continuous(limits=c(60,140), breaks=seq(60,140,10)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2))
dev.off()

Khats3 <- as.data.frame(Khats3)
colnames(Khats3) <- c('k1')

o <- paste0("BKW1",p,".png",sep="")
png(filename = o,width = 500,height = 500,
    units = "px",bg = "white",res = 100)
ggplot(Khats3,aes(y=..density..)) +
  geom_histogram(aes(x = k1),fill="skyblue2",colour="skyblue2",binwidth = 1) + 
  xlab('T') + ylab('Frequency') + 
  scale_x_continuous(limits=c(60,140), breaks=seq(60,140,10)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2))
dev.off()

Khats4 <- as.data.frame(Khats4)
colnames(Khats4) <- c('k1','k2')

o <- paste0("BKW2",p,".png",sep="")
png(filename = o,width = 500,height = 500,
    units = "px",bg = "white",res = 100)
ggplot(Khats4,aes(y=..density..)) +
  geom_histogram(aes(x = k1, fill = 'k1'), binwidth = 1) +
  geom_histogram(aes(x = k2, fill = 'k2'), binwidth = 1) +
  xlab('T') + ylab('Frequency') + 
  scale_fill_manual(name = '', values = c('k1' = 'skyblue2', 'k2' = 'orangered1'),
                    labels = c(expression(widehat(k)[1]),expression(widehat(k)[2]))) + 
  scale_x_continuous(limits=c(60,140), breaks=seq(60,140,10)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) + 
  theme(legend.position=c(0.9,0.9),legend.background=element_blank())+
  theme(legend.key.size=unit(10,'pt'))
dev.off()


Khats5 <- as.data.frame(Khats5)
colnames(Khats5) <- c('k1','k2')

o <- paste0("DHB2",p,".png",sep="")
png(filename = o,width = 500,height = 500,
    units = "px",bg = "white",res = 100)
ggplot(Khats5,aes(y=..density..)) +
  geom_histogram(aes(x = k1, fill = 'k1'), binwidth = 1) +
  geom_histogram(aes(x = k2, fill = 'k2'), binwidth = 1) +
  xlab('T') + ylab('Frequency') + 
  scale_fill_manual(name = '', values = c('k1' = 'skyblue2', 'k2' = 'orangered1'),
                    labels = c(expression(widehat(k)[1]),expression(widehat(k)[2]))) + 
  scale_x_continuous(limits=c(60,140), breaks=seq(60,140,10)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) + 
  theme(legend.position=c(0.9,0.9),legend.background=element_blank())+
  theme(legend.key.size=unit(10,'pt'))
dev.off()