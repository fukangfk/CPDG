# rm(list = ls())
source('MainFunction.R')
############################################################
tot.sim <- 1000
# Scenario 1
dgpNs <- c(100,200)
dgpM <- 3
dgpR <- cbind(c(2,2,2,2),c(2,2,2,2))
dgpPs <- c(25, 50)
dgpdeltas <- c(0.6,0.8)
Mat <- NULL
for (nnn in 1:length(dgpNs)) {
  dgpN <- dgpNs[nnn]
  Ks <- matrix(c(0.4,0.5,0.6,
                 0.6,0.5,0.4),nrow = 3,ncol = 2)*dgpN
  for (ddd in 1:length(dgpdeltas)) {
    dgpdelta <- dgpdeltas[ddd]
    for (ppp in 1:length(dgpPs)) {
      dgpPm <- dgpPs[ppp]
      dgpPd <- rep(dgpPm, dgpM)
      glabel <- rep(1:dgpM, each = dgpPm)
      
      for (kkk in 1:nrow(Ks)) {
        dgpK <- Ks[kkk,]
        Khats1 <- matrix(rep(0,2*tot.sim),ncol = 2)
        for (ttt in 1:tot.sim) {
          res.Y <- dgpf(dgpN, dgpPd, dgpR, dgpM, dgpK, dgpdelta)
          Y <- res.Y$dY
          rvec <- rowSums(dgpR)
          CP.result <- QMLE.CP(Y,glabel,rvec)
          Khats1[ttt,] <- CP.result$Khat
        }
        
        # MAE
        MAE1 <- sum(abs(Khats1[,1] - dgpK[1]))/tot.sim
        MAE2 <- sum(abs(Khats1[,2] - dgpK[2]))/tot.sim
        # RMSE
        RMSE1 <- sqrt(sum((Khats1[,1] - dgpK[1])^2)/tot.sim)
        RMSE2 <- sqrt(sum((Khats1[,2] - dgpK[2])^2)/tot.sim)
        
        K1 <- dgpK[1]; K2 <- dgpK[2]
        Mat1 <- cbind(dgpN,dgpPm,dgpdelta,K1,K2,MAE1,RMSE1,MAE2,RMSE2)
        Mat <- rbind(Mat,Mat1)
        o <- paste0("result1.csv", sep="")
        write.csv(Mat, o)
      }
    }
  }
}

# Scenario 2
dgpNs <- c(100,200)
dgpM <- 3
dgpR <- cbind(c(2,2,2,2),c(2,2,2,2))
dgpPs <- c(25, 50)
dgpdeltas <- c(0.6,0.8)
dgpC <- matrix(c(1,0,2,1),nrow = 2, byrow = TRUE)
dgpD <- matrix(c(1,2,0,1),nrow = 2, byrow = TRUE)
Mat <- NULL
for (nnn in 1:length(dgpNs)) {
  dgpN <- dgpNs[nnn]
  Ks <- matrix(c(0.4,0.5,0.6,
                 0.6,0.5,0.4),nrow = 3,ncol = 2)*dgpN
  for (ddd in 1:length(dgpdeltas)) {
    dgpdelta <- dgpdeltas[ddd]
    for (ppp in 1:length(dgpPs)) {
      dgpPm <- dgpPs[ppp]
      dgpPd <- rep(dgpPm, dgpM)
      glabel <- rep(1:dgpM, each = dgpPm)
      
      for (kkk in 1:nrow(Ks)) {
        dgpK <- Ks[kkk,]
        Khats1 <- matrix(rep(0,2*tot.sim),ncol = 2)
        for (ttt in 1:tot.sim) {
          res.Y <- dgpf2(dgpN, dgpPd, dgpR, dgpM, dgpK, dgpdelta, dgpC, dgpD)
          Y <- res.Y$dY
          rvec <- dgpR[,1]
          CP.result <- QMLE.CP(Y,glabel,rvec)
          Khats1[ttt,] <- CP.result$Khat
        }
        
        # MAE
        MAE1 <- sum(abs(Khats1[,1] - dgpK[1]))/tot.sim
        MAE2 <- sum(abs(Khats1[,2] - dgpK[2]))/tot.sim
        # RMSE
        RMSE1 <- sqrt(sum((Khats1[,1] - dgpK[1])^2)/tot.sim)
        RMSE2 <- sqrt(sum((Khats1[,2] - dgpK[2])^2)/tot.sim)
        
        K1 <- dgpK[1]; K2 <- dgpK[2]
        Mat1 <- cbind(dgpN,dgpPm,dgpdelta,K1,K2,MAE1,RMSE1,MAE2,RMSE2)
        Mat <- rbind(Mat,Mat1)
        o <- paste0("result2.csv", sep="")
        write.csv(Mat, o)
      }
    }
  }
}

# Scenario 3
dgpNs <- c(100,200)
dgpM <- 3
dgpR <- cbind(c(4,4,4,4),c(4,4,4,4))
dgpPs <- c(25, 50)
dgpdeltas <- c(0.6,0.8)
dgpC <- matrix(c(1,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0),
               nrow = 4, byrow = TRUE)
dgpD <- matrix(c(1,2,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
               nrow = 4, byrow = TRUE)
Mat <- NULL
for (nnn in 1:length(dgpNs)) {
  dgpN <- dgpNs[nnn]
  Ks <- matrix(c(0.4,0.5,0.6,
                 0.6,0.5,0.4),nrow = 3,ncol = 2)*dgpN
  for (ddd in 1:length(dgpdeltas)) {
    dgpdelta <- dgpdeltas[ddd]
    for (ppp in 1:length(dgpPs)) {
      dgpPm <- dgpPs[ppp]
      dgpPd <- rep(dgpPm, dgpM)
      glabel <- rep(1:dgpM, each = dgpPm)
      
      for (kkk in 1:nrow(Ks)) {
        dgpK <- Ks[kkk,]
        Khats1 <- matrix(rep(0,2*tot.sim),ncol = 2)
        for (ttt in 1:tot.sim) {
          res.Y <- dgpf2(dgpN, dgpPd, dgpR, dgpM, dgpK, dgpdelta, dgpC, dgpD)
          Y <- res.Y$dY
          rvec <- dgpR[,1]
          CP.result <- QMLE.CP(Y,glabel,rvec)
          Khats1[ttt,] <- CP.result$Khat
        }
        
        # MAE
        MAE1 <- sum(abs(Khats1[,1] - dgpK[1]))/tot.sim
        MAE2 <- sum(abs(Khats1[,2] - dgpK[2]))/tot.sim
        # RMSE
        RMSE1 <- sqrt(sum((Khats1[,1] - dgpK[1])^2)/tot.sim)
        RMSE2 <- sqrt(sum((Khats1[,2] - dgpK[2])^2)/tot.sim)
        
        K1 <- dgpK[1]; K2 <- dgpK[2]
        Mat1 <- cbind(dgpN,dgpPm,dgpdelta,K1,K2,MAE1,RMSE1,MAE2,RMSE2)
        Mat <- rbind(Mat,Mat1)
        o <- paste0("result3.csv", sep="")
        write.csv(Mat, o)
      }
    }
  }
}

