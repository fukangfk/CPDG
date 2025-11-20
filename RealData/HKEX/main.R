rm(list = ls())
library(ggplot2)
source("MainFunction.R")
source("BKW.R")
source("DBH.R")
load("HKstock.Rdata")
hsdm0 <- HKstock$hsdm
glabel <- HKstock$label[,2]
dates <- HKstock$dates

label <- dplyr::case_when(
  glabel %in% c(1, 3, 4, 9, 11) ~ 1,
  glabel %in% c(2, 5, 6, 8, 10) ~ 2,
  glabel %in% c(7, 12) ~ 3,
)
table(label)

dN <- nrow(hsdm0); dN
dp <- ncol(hsdm0); dp

res.R3 <- estinumfactor(hsdm0,label)
ppp <- 1
Rmat <- res.R3$R
Rratio <- data.frame(j=1:20,R=Rmat[Rmat[,1]==ppp,][1:20,2])
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

QMLE.CP(hsdm0,label,rvec) # dates0[c(736,383)]
DBHQML.CP(hsdm0,1) # dates0[c(736)]
BKWLS.CP(hsdm0,2) # dates0[c(307,503)]



# split full sample
Y1 <- hsdm0[1:307,]
Y2 <- hsdm0[308:503,]
Y3 <- hsdm0[504:dN,]


res.R33 <- estinumfactor(Y3,label)
ppp <- 3
Rmat <- res.R33$R
Rratio <- data.frame(j=4:20,R=Rmat[Rmat[,1]==ppp,][4:20,2])
ggplot(data = Rratio, aes(x=j,y=R)) +
  geom_line() + geom_point(shape=21,size=2,fill="white") + 
  xlab(expression(j)) + ylab(expression(hat(R)^{(1)})) +
  geom_text(aes(label=j),nudge_y=0.02,alpha=1)


# r0hat <- c(2)
# r1hat <- c(2,7)
# r2hat <- c(3,6)
# r3hat <- c(3,7)
# Rcomb <- as.matrix(expand.grid(r0hat,r1hat,r2hat,r3hat))
# which(rowSums(Rcomb) == 14)
# Rcomb[which(rowSums(Rcomb) == 14),]

rvec1 <- c(2,6,4,3)
sum(rvec1)
rvec2 <- c(2,8,2,3)
sum(rvec2)
rvec3 <- c(2,3,3,7) 
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

