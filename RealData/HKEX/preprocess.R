rm(list = ls())
hs <- read.csv("data_HK.csv")
sectorinfo <- read.csv("HK_sectors.csv",header = TRUE)


start.date <- which(hs[,1] == "2019/12/9")
end.date <- which(hs[,1] == "2023/12/8")

hs0 <- hs[c(start.date:end.date), -c(1)];
date0 <- hs[c(start.date:end.date),1]

dT0 <- nrow(hs0); dp0 <- ncol(hs0)
# calculate the daily log-returns of theses stocks
hs1 <- log(hs0[1:(dT0-1),]) - log(hs0[2:dT0,])


hsdm <- scale(hs1,center = TRUE, scale = FALSE)
label <- data.frame(index=sectorinfo$Stock,sector=sectorinfo$label)

HKstock <- list(hsdm = hsdm,
                label = label,
                dates = date0[-1])
save(HKstock,file = "HKstock.Rdata")




