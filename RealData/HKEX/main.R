rm(list = ls())
library(ggplot2)
source("MainFunction.R")
source("BKW.R")
source("DBH.R")
load("HKstock.Rdata")
hsdm <- HKstock$hsdm
glabel <- HKstock$label[,2]
dates <- HKstock$dates

label <- dplyr::case_when(
  glabel %in% c(1, 3, 4, 9, 11) ~ 1,
  glabel %in% c(2, 5, 6, 8, 10) ~ 2,
  glabel %in% c(7, 12) ~ 3,
)
table(label)

date1 <- which(dates == "2020/1/2")
date2 <- which(dates == "2024/1/2")
hsdm0 <- hsdm[date1:date2,]

dN <- nrow(hsdm0); dp <- ncol(hsdm0)

res.R3 <- estinumfactor(hsdm0,label)
ppp <- 3
Rmat <- res.R3$R
Rratio <- data.frame(j=3:20,R=Rmat[Rmat[,1]==ppp,][3:20,2])
ggplot(data = Rratio, aes(x=j,y=R)) +
  geom_line() + 
  geom_point(shape=21,size=2,fill="white") + 
  xlab(expression(j)) + 
  #ylab(expression(hat(R)[j])) +
  ylab(expression(hat(R)[j]^{(3)})) +
  geom_text(aes(label=j),nudge_y=0.02,alpha=1)+
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12)
  )


# r0hat <- c(2)
# r1hat <- c(6,4,2)
# r2hat <- c(8,6,4,1)
# r3hat <- c(6,3,1)
# Rcomb <- as.matrix(expand.grid(r0hat,r1hat,r2hat,r3hat))
# which(rowSums(Rcomb) == 13)
# Rcomb[which(rowSums(Rcomb) == 13),]

rvec <- c(2,4,6,3) # dates[c(736,383)]
# rvec <- c(2,4,4,3) # dates[c(736,383)]
# rvec <- c(2,2,6,3) # dates[c(736,383)]
# rvec <- c(2,4,6,1) # dates[c(736,383)]
# rvec <- c(2,2,8,1) # dates[c(736,383)]
sum(rvec)

QMLE.CP(hsdm0,label,rvec) # dates[c(736,383)]
DBHQML.CP(hsdm0,1) # dates[c(736)]
BKWLS.CP(hsdm0,1) # dates[c(736)]


# split subsample

Y1 <- hsdm0[seq(1,383,1),]
Y2 <- hsdm0[seq(383,736,1),]
Y3 <- hsdm0[seq(737,dN,1),]


res.R33 <- estinumfactor(Y2,label)
ppp <- 0
Rmat <- res.R33$R
Rratio <- data.frame(j=3:20,R=Rmat[Rmat[,1]==ppp,][3:20,2])
ggplot(data = Rratio, aes(x=j,y=R)) +
  geom_line() + geom_point(shape=21,size=2,fill="white") + 
  xlab(expression(j)) + ylab(expression(hat(R)^{(1)})) +
  geom_text(aes(label=j),nudge_y=0.02,alpha=1)

r0hat <- c(1)
r1hat <- c(4,1)
r2hat <- c(4,1)
r3hat <- c(5,3)
Rcomb <- as.matrix(expand.grid(r0hat,r1hat,r2hat,r3hat))
which(rowSums(Rcomb) == 11)
Rcomb[which(rowSums(Rcomb) == 11),]

rvec1 <- c(2,6,4,1) 
sum(rvec1)
rvec2 <- c(2,4,1,5) 
sum(rvec2)
rvec3 <- c(1,4,1,5) 
sum(rvec3)

res1 <- gwest(Y1,label,rvec1)
res2 <- gwest(Y2,label,rvec2)
res3 <- gwest(Y3,label,rvec3)
A1hat <- res1$Ahat
A2hat <- res2$Ahat
A3hat <- res3$Ahat
B1hat <- res1$Bhat
B2hat <- res2$Bhat
B3hat <- res3$Bhat

dd1 <- A1hat%*%(solve(t(A1hat)%*%A1hat))%*%t(A1hat)
dd2 <- A2hat%*%(solve(t(A2hat)%*%A2hat))%*%t(A2hat)
dd3 <- A3hat%*%(solve(t(A3hat)%*%A3hat))%*%t(A3hat)

norm(dd1-dd2,'F')
norm(dd2-dd3,'F')


dd1 <- B1hat%*%(solve(t(B1hat)%*%B1hat))%*%t(B1hat)
dd2 <- B2hat%*%(solve(t(B2hat)%*%B2hat))%*%t(B2hat)
dd3 <- B3hat%*%(solve(t(B3hat)%*%B3hat))%*%t(B3hat)


norm(dd1-dd2,'F')
norm(dd2-dd3,'F')

#########################################################
source("Function/MainFunction.R")
source("Function/BKW.R")
source("Function/DBH.R")
source("Function/ZPYZ.R")
load("SP500stock.Rdata")

hsdm <- sp500stock$hsdm
glabel <- sp500stock$label[,1]
dates <- sp500stock$dates

res.R <- estifactor(hsdm,5)
Rratio <- data.frame(j=2:20,R=res.R$R[2:20])
ggplot(data = Rratio, aes(x=j,y=R)) +
  geom_line() + geom_point(shape=21,size=2,fill="white") + 
  xlab(expression(j)) + ylab(expression(hat(R)[j])) +
  scale_y_continuous(limits = c(1,2.2)) +
  geom_text(aes(label=j),nudge_y=0.05,alpha=1)
rv <- c(1,16)
dN <- nrow(hsdm); dp <- ncol(hsdm)
tau <- 1-1/log(dN)
res.AB <- LFMest(hsdm,5,rv)
Bhat <- res.AB$B
res.nogroup <- nocluster(Bhat,2)
id.group <- res.nogroup$ind.cluster1
id.nogroup <- res.nogroup$ind.cluster0
num.group <- res.nogroup$numcluster1
num.nogroup <- res.nogroup$numcluster0
MMM <- 3
B0hat <- Bhat[-id.nogroup,]

hsdm0 <- hsdm[,id.group] # 1258 x 466 data matrix
label <- clusterppp(B0hat,MMM) # group information
dN <- nrow(hsdm0); dp <- ncol(hsdm0)

res.R3 <- estinumfactor(hsdm0,label)
ppp <- 3
Rmat <- res.R3$R
Rratio <- data.frame(j=2:20,R=Rmat[Rmat[,1]==ppp,][2:20,2])
ggplot(data = Rratio, aes(x=j,y=R)) +
  geom_line() + geom_point(shape=21,size=2,fill="white") + 
  xlab(expression(j)) + ylab(expression(hat(R)^{(1)})) +
  geom_text(aes(label=j),nudge_y=0.02,alpha=1)

r0hat <- c(1)
r1hat <- c(1,3,6,9)
r2hat <- c(1,4,6,9,11,13)
r3hat <- c(3,5,10,12)
Rcomb <- as.matrix(expand.grid(r0hat,r1hat,r2hat,r3hat))
which(rowSums(Rcomb) == 16) # 11 14 17 28 31 37 53
Rcomb[which(rowSums(Rcomb) == 16),]

# rvec <- c(1,6,6,3) # dates[c(384,612)]
# rvec <- c(1,3,9,3) # dates[c(384,711)]
# rvec <- c(1,1,11,3) # dates[c(384,711)]
# rvec <- c(1,9,1,5) # dates[c(384,613)]
# rvec <- c(1,6,4,5) # dates[c(384,612)]
# rvec <- c(1,1,9,5) # dates[c(384,711)]
# rvec <- c(1,1,4,10) # dates[c(384,590)]
sum(rvec)

QMLE.CP(hsdm0,label,rvec) # dates[c(736,383)]
DBHQML.CP(hsdm0,1) # dates[c(384)]
DBHQML.CP(hsdm0,2) # dates[c(384,779)]
BKWLS.CP(hsdm0,1) # dates[c(384)]
BKWLS.CP(hsdm0,2) # dates[c(384,779)]

#########################################################
# group is unknown
source("Function/MainFunction.R")
source("Function/BKW.R")
source("Function/DBH.R")
source("Function/ZPYZ.R")
load("SP500stock.Rdata")

hsdm <- sp500stock$hsdm
glabel <- sp500stock$label[,1]
dates <- sp500stock$dates

res.R <- estifactor(hsdm,5)
Rratio <- data.frame(j=2:20,R=res.R$R[2:20])
ggplot(data = Rratio, aes(x=j,y=R)) +
  geom_line() + geom_point(shape=21,size=2,fill="white") + 
  xlab(expression(j)) + ylab(expression(hat(R)[j])) +
  scale_y_continuous(limits = c(1,2.2)) +
  geom_text(aes(label=j),nudge_y=0.05,alpha=1)
rv <- c(1,16)
dN <- nrow(hsdm); dp <- ncol(hsdm)

rvec <- c(1,15)
label <- rep(1,dp)

QMLE.CP(hsdm,label,rvec) # dates[c(384,475)]
DBHQML.CP(hsdm,1) # dates[c(384)]
DBHQML.CP(hsdm,2) # dates[c(384,779)]
BKWLS.CP(hsdm,1) # dates[c(384)]
BKWLS.CP(hsdm,2) # dates[c(384,779)]