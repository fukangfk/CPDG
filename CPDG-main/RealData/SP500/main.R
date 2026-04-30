library(ggplot2)
#########################################################
source("MainFunction.R")
load("SP500stock.Rdata")


hsdm <- sp500stock$hsdm
dates <- sp500stock$dates


dN <- nrow(hsdm); dN
dp <- ncol(hsdm); dp

res.R <- estinumfactor(hsdm)
Rratio <- data.frame(j=2:20,R=res.R$R[2:20,2])
ggplot(data = Rratio, aes(x=j,y=R)) +
  geom_line() + geom_point(shape=21,size=2,fill="white") + 
  xlab(expression(j)) + ylab(expression(hat(R)[j])) +
  scale_y_continuous(limits = c(1,1.75)) +
  geom_text(aes(label=j),nudge_y=0.05,alpha=1)

rvec <- c(1,12)
QMLE.CP(hsdm,rvec) # dates[c(384,483)]




