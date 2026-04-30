rm(list = ls())
library(ggplot2)
source("MainFunction.R")
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

QMLE.CP(hsdm0,label,rvec)



