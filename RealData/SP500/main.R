library(ggplot2)
#########################################################
source("MainFunction.R")
load("SP500stock.Rdata")


hsdm <- sp500stock$hsdm
dates <- sp500stock$dates
hsdm0 <- hsdm
dates0 <- dates

# date1 <- which(dates == "2015-12-31")
# date2 <- which(dates == "2019-12-30")
# hsdm0 <- hsdm[date1:date2,]
# dates0 <- dates[date1:date2]

dN <- nrow(hsdm0); dN
dp <- ncol(hsdm0); dp

res.R <- estinumfactor(hsdm0)
Rratio <- data.frame(j=2:20,R=res.R$R[2:20,2])
ggplot(data = Rratio, aes(x=j,y=R)) +
  geom_line() + geom_point(shape=21,size=2,fill="white") + 
  xlab(expression(j)) + ylab(expression(hat(R)[j])) +
  scale_y_continuous(limits = c(1,1.75)) +
  geom_text(aes(label=j),nudge_y=0.05,alpha=1)

rvec <- c(1,12)
QMLE.CP(hsdm0,rvec) # dates0[c(384,483)]


# split full sample
Y1 <- hsdm0[1:384,]
Y2 <- hsdm0[385:483,]
Y3 <- hsdm0[485:dN,]


res.R33 <- estinumfactor(Y3)
Rmat <- res.R33$R
Rratio <- data.frame(j=2:20,R=Rmat[2:20,2])
ggplot(data = Rratio, aes(x=j,y=R)) +
  geom_line() + geom_point(shape=21,size=2,fill="white") + 
  xlab(expression(j)) + ylab(expression(hat(R))) +
  geom_text(aes(label=j),nudge_y=0.02,alpha=1)


rvec1 <- c(1,8)
rvec2 <- c(2,7)
rvec3 <- c(2,10)

res1 <- gwest(Y1,rvec1)
res2 <- gwest(Y2,rvec2)
res3 <- gwest(Y3,rvec3)
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



