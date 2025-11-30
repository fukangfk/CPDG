rm(list = ls())
hs <- read.csv("data_HK.csv")
sectorinfo <- read.csv("HK_sectors.csv",header = TRUE)



hs0 <- hs[, -c(1)];dates <- hs[,1]

dT0 <- nrow(hs0); dp0 <- ncol(hs0)
# calculate the daily log-returns of theses stocks
hs1 <- log(hs0[1:(dT0-1),]) - log(hs0[2:dT0,])

hsdm <- scale(hs1,center = TRUE, scale = FALSE)
label <- data.frame(index=sectorinfo$Stock,sector=sectorinfo$label)

start.date <- which(hs[,1] == "2019/12/9")
end.date <- which(hs[,1] == "2023/12/8")
HKstock <- list(hsdm = hsdm[c(start.date:end.date),],
                label = label,
                dates = dates[c((start.date+1):end.date)])
save(HKstock,file = "HKstock.Rdata")





