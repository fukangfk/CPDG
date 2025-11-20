library(ggplot2)
library(ggpattern)
library(dplyr)


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
ggplot(hist_data, 
       aes(x = x, 
           y = y, 
           pattern = group,
           fill = group,
           color = group,
           linetype = group)
       ) +
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
  # ----------------- 坐标轴 -----------------
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

