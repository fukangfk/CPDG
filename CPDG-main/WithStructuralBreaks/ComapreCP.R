# rm(list = ls())
source('MainFunction.R')
source('DBH.R')
source('BKW.R')
library(ggplot2)
library(ggpattern)
library(dplyr)

dgpN <- 200
dgpM <- 3
dgpK <- as.integer(c(0.6,0.4)*dgpN)
dgpPm <- 50
dgpPd <- rep(dgpPm, dgpM)
glabel <- rep(1:dgpM, each = dgpPm)
dgpdelta <- 0.8
tot.sim <- 100
scenario <- 2
Khats1 <- Khats4 <- Khats5 <- matrix(rep(0,2*tot.sim),ncol = 2)
Khats2 <- Khats3 <- numeric(tot.sim)
for (ttt in 1:tot.sim) {
  if (scenario == 1){
    dgpR <- cbind(c(2,2,2,2),c(2,2,2,2))
    res.Y <- dgpf(dgpN, dgpPd, dgpR, dgpM, dgpK, dgpdelta)
    # rvec <- rowSums(dgpR)
    
    Y <- res.Y$dY
    rvec <- estinumfactor(Y,glabel)$rvechat
    while(sum(rvec>0)!=4){
      res.Y <- dgpf(dgpN, dgpPd, dgpR, dgpM, dgpK, dgpdelta)
      Y <- res.Y$dY
      rvec <- estinumfactor(Y,glabel)$rvechat
    }
  } else {
    if(scenario == 2){
      dgpR <- cbind(c(2,2,2,2),c(2,2,2,2))
      dgpC <- matrix(c(0.5,1,0,0.5),nrow = 2, byrow = TRUE)
      dgpD <- matrix(c(0.5,0,1,0.5),nrow = 2, byrow = TRUE)
      res.Y <- dgpf2(dgpN, dgpPd, dgpR, dgpM, dgpK, dgpdelta, dgpC, dgpD)
      # rvec <- dgpR[,1]
      
      Y <- res.Y$dY
      rvec <- estinumfactor(Y,glabel)$rvechat
      while(sum(rvec>0)!=4){
        res.Y <- dgpf2(dgpN, dgpPd, dgpR, dgpM, dgpK, dgpdelta, dgpC, dgpD)
        Y <- res.Y$dY
        rvec <- estinumfactor(Y,glabel)$rvechat
      }
    } else {
      dgpR <- cbind(c(4,4,4,4),c(4,4,4,4))
      dgpC <- matrix(c(0.5,1,0,0, 0,0.5,0,0, 0,0,0,0, 0,0,0,0),nrow = 4, byrow = TRUE)
      dgpD <- matrix(c(0.5,0,0,0, 1,0.5,0,0, 0,0,0,0, 0,0,0,0),nrow = 4, byrow = TRUE)
      res.Y <- dgpf2(dgpN, dgpPd, dgpR, dgpM, dgpK, dgpdelta, dgpC, dgpD)
      # rvec <- dgpR[,1]
      
      Y <- res.Y$dY
      rvec <- estinumfactor(Y,glabel)$rvechat
      while(sum(rvec>0)!=4 || rvec[1]<4){
        res.Y <- dgpf2(dgpN, dgpPd, dgpR, dgpM, dgpK, dgpdelta, dgpC, dgpD)
        Y <- res.Y$dY
        rvec <- estinumfactor(Y,glabel)$rvechat
      }
    }
  }
  
  
  CP.result1 <- QMLE.CP(Y,glabel)
  CP.result2 <- DBHQML.CP(Y,1)
  # CP.result3 <- BKWLS.CP(Y,1)
  CP.result4 <- BKWLS.CP(Y,2)
  # CP.result5 <- DBHQML.CP(Y,2)
  
  Khats1[ttt,] <- CP.result1$Khat
  Khats2[ttt] <- CP.result2$Khat
  # Khats3[ttt] <- CP.result3$Khat
  Khats4[ttt,] <- CP.result4$Khat
  # Khats5[ttt,] <- CP.result5$Khat
  if(ttt %% 10 ==  0 ){cat("...",ttt)}
  if(ttt %% 100 == 0 ){cat("\n") }
}


data <- data.frame(
  value = c(Khats1[,1], Khats1[,2],
            Khats2, 
            Khats4[,1], Khats4[,2]),
  group = factor(
    rep(c("global", "group", "DBH", "BKW1", "BKW2"), each = 100),
    levels = c("global", "group", "DBH", "BKW1", "BKW2"))
)
up_groups   <- c("global", "group")
down_groups <- c("DBH", "BKW1", "BKW2")

# ----------------- 手动计算直方图 -----------------
bins <- 15
hist_data <- data %>%
  group_by(group) %>%
  do({
    h <- hist(.$value, breaks = seq(50, 150, length.out = 16), plot = FALSE)
    data.frame(
      group = unique(.$group),
      x = h$mids,
      y = h$counts,
      y_dir = ifelse(unique(.$group) %in% up_groups, 1, -1)
    )
  }) %>%
  mutate(y = y * y_dir)

# ----------------- pattern/color/fill 设置 -----------------
pattern_values <- c("stripe", "circle", "none", "stripe", "circle")
line_types <- c("dashed", "solid", "dotdash", "dotted", "longdash")
fill_values <- c(
  "global" = "#1f77b4",  # 蓝色
  "group"  = "#aec7e8",  # 浅蓝色
  "DBH"    = "#ff7f0e",  # 橙色
  "BKW1"   = "#ffbb78",  # 浅橙色
  "BKW2"   = "#d62728"   # 红色
)
color_values <- c(
  "global" = "#1f55a0",  # 深蓝色
  "group"  = "#7ea6d6",  # 中蓝色
  "DBH"    = "#e67300",  # 深橙色
  "BKW1"   = "#e6994c",  # 中橙色
  "BKW2"   = "#b31a1a"   # 深红色
)

# ----------------- 图例标签 -----------------
legend_labels <- c(
  expression("FLZZG:" ~ widehat(k)),
  expression("FLZZG:" ~ widehat(k)^g),
  expression("DBH:" ~ widehat(k)),
  expression("BKW2:" ~ widehat(k)[1]),
  expression("BKW2:" ~ widehat(k)[2])
)


# ----------------- 绘图 -----------------
ggplot(hist_data, aes(x = x, y = y, pattern = group, 
                      fill = group, color = group, 
                      linetype = group)) +
  geom_col_pattern(
    width = (max(data$value)-min(data$value))/bins * 0.9,
    alpha = 0.3,
    pattern_density = 0.3,
    pattern_spacing = 0.04,
    pattern_fill = "white",
    pattern_colour = "black",
    size = 1.5,
    position = "identity"
  ) +
  scale_pattern_manual(values = pattern_values,
                       labels = legend_labels,
                       name = NULL) +
  scale_fill_manual(values = fill_values,
                    labels = legend_labels,
                    name = NULL) +
  scale_color_manual(values = color_values,
                     labels = legend_labels,
                     name = NULL) +
  scale_linetype_manual(values = line_types,
                        labels = legend_labels,
                        name = NULL) +
  coord_cartesian(xlim = c(50,150), ylim = c(-100,100)) +
  scale_x_continuous(breaks = seq(50,150,10)) +
  scale_y_continuous(breaks = seq(-100,100,20), labels = abs) +
  labs(x = "T", y = "Frequency", title = "") +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1), 
    plot.title = element_text(hjust = 0.5, size = 16)
  )


