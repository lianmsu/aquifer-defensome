# =============================================================================
# aquifer-defensome -- supplementary figures / Mobile resistance genes and the defensome
# =============================================================================
# Regression between MRG content and DS number, with the MazEF panel.
#
# Manuscript      : Fig. S13
# Main outputs    : 1_arg_dss_regresion.pdf, 1_MazEF_mrg3.pdf, 2_MazEF_mrg2.pdf
# Original script : R codes/Fig.S13.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA, the root of the data
# directory (its layout mirrors the one used for the study). The default is the
# current working directory; set the environment variable to point somewhere
# else instead of editing this file:
#
#     Sys.setenv(GTR_DATA = "/path/to/defensome")
#
# Scripts of the peer-review revision additionally read from GTR_REV (the
# revision working copy, default: the current directory).
# =============================================================================
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)

#对mags中的arg进行分析

rm(list=ls()) # Clear workspace
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)

#1. 读取数据
setwd(file.path(GTR_DATA, "groundwater/10mag_arg"))
magarg <- read.delim("mrg_gene_out1_type.csv")


#2. 处理数据
#mags中arg的数据
library(dplyr)
rownames(magarg) <- paste0(magarg[,1], ".fasta")
magarg <- magarg[,-1]
magarg_sum <- magarg %>% rowSums() %>% as.data.frame() %>% setNames("argsum")

#mags中dss的数据读取和处理
setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))

setwd(file.path(GTR_DATA, "groundwater/02去冗余信息"))#读取所有27578个mags的数据
allmags <- read.delim2('Wdb_ANI99.csv',sep = ',',
                       quote = "",header = T
) 

setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))#读取含有dss的mags的数据
dss <- read.delim2('ds_counts_matrix.csv',sep = ',',row.names = 1,
                   quote = "",header = T
) 

not_in_dss <- allmags$genome[!allmags$genome %in% rownames(dss)]# 找出 allmags 中在 dss 行名中不存在的 genome 值#将不含dss的mags也合并到一起

new_rows <- data.frame(
  genome = not_in_dss,  # 重复 genome 值
  stringsAsFactors = FALSE)# 使用不在 dss 行名中的 genome 值创建新的数据框

col_names <- setdiff(names(dss), "genome")  # 获取 dss 中除了 genome 之外的列名# 为 new_rows 添加与 dss 相同的列，并将这些列的值填充为 0
new_rows[col_names] <- lapply(col_names, function(col_name) rep(0, nrow(dss)))
rownames(new_rows) <- new_rows$genome
new_rows <- new_rows[,-1]

dss <- rbind(dss,new_rows)#合并new_rows与dss

magdss_sum <- dss %>% rowSums() %>% as.data.frame() %>% setNames("dsssum")#求和每个mags中的dss
magdss_sum <- magdss_sum %>%
  filter(rownames(.) %in% rownames(magarg_sum))

magdss_sum <- data.frame(genome = rownames(dss) , dsssum = dss$MazEF)#分别对五个主要的DS进行探究：
rownames(magdss_sum) <- magdss_sum$genome
magdss_sum <- select(magdss_sum, dsssum)

magdss_sum <- select(dss,-RM, -SoFIC, -CasFinder, -AbiE, -MazEF)
magdss_sum <- magdss_sum %>% rowSums() %>% as.data.frame() %>% setNames("dsssum")#求和每个mags中的dss
magdss_sum <- magdss_sum %>%
  filter(rownames(.) %in% rownames(magarg_sum))

magdss_sum <- select(dss,RM, SoFIC, CasFinder, AbiE, MazEF)
magdss_sum <- magdss_sum %>% rowSums() %>% as.data.frame() %>% setNames("dsssum")#求和每个mags中的dss
magdss_sum <- magdss_sum %>%
  filter(rownames(.) %in% rownames(magarg_sum))

#mags中dgs的数据
setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))

setwd(file.path(GTR_DATA, "groundwater/02去冗余信息"))#读取所有27578个mags的数据
allmags <- read.delim2('Wdb_ANI99.csv',sep = ',',
                       quote = "",header = T
) 

setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))#读取含有dss的mags的数据
dgs <- read.delim2('dg_counts_matrix.csv',sep = ',',row.names = 1,
                   quote = "",header = T
) 

not_in_dgs <- allmags$genome[!allmags$genome %in% rownames(dgs)]# 找出 allmags 中在 dss 行名中不存在的 genome 值#将不含dss的mags也合并到一起

new_rows <- data.frame(
  genome = not_in_dgs,  # 重复 genome 值
  stringsAsFactors = FALSE)# 使用不在 dss 行名中的 genome 值创建新的数据框

col_names <- setdiff(names(dgs), "genome")  # 获取 dss 中除了 genome 之外的列名# 为 new_rows 添加与 dss 相同的列，并将这些列的值填充为 0
new_rows[col_names] <- lapply(col_names, function(col_name) rep(0, nrow(dgs)))
rownames(new_rows) <- new_rows$genome
new_rows <- new_rows[,-1]

dgs <- rbind(dgs,new_rows)#合并new_rows与dss

magdgs_sum <- dgs %>% rowSums() %>% as.data.frame() %>% setNames("dgssum")#求和每个mags中的dss
magdgs_sum <- magdgs_sum %>%
  filter(rownames(.) %in% rownames(magarg_sum))



#3. 统计分析
#dss与arg的回归分析
# 按行名合并两个数据框
magdss_sum$genome <- rownames(magdss_sum)
magarg_sum$genome <- rownames(magarg_sum)
merge_dss_arg <- merge(magdss_sum,magarg_sum,by = 'genome')

# 计算相关系数
correlation <- cor(merge_dss_arg$dsssum, merge_dss_arg$argsum)
print(correlation)

# 进行线性回归分析
model <- lm(argsum ~ dsssum, data = merge_dss_arg)
# 查看回归模型的摘要
summary(model)


#dgs与arg的回归分析
magdgs_sum$genome <- rownames(magdgs_sum)
merge_dss_arg_dgs <- merge(merge_dss_arg,magdgs_sum,by = 'genome')




#4. 可视化
#dss与arg
library(ggplot2)
# 绘制散点图和回归线
p1 <- ggplot(merge_dss_arg, aes(x = dsssum, y = argsum)) +
  geom_point() +  # 添加散点图
  geom_smooth(method = "nls", se = TRUE) +  # 添加线性回归线，不显示置信区间
  labs(title = NULL, x = "DSs number (per MAG)", y = "ARGs number (per MAG)")+  # 添加标题和坐标轴标签
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.8,0.7),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p1
setwd(file.path(GTR_DATA, "groundwater/10mag_arg"))
ggsave("1_arg_dss_regresion.pdf", plot = p1, width = 5.5, height = 4.5, units = "in")

#参考华大海洋nature的图
# Read data
df.count <- merge_dss_arg %>%
  mutate(dsssum = factor(dsssum, levels = 0:83)) %>%
  mutate(group = case_when(dsssum == 1 ~ "group1", TRUE ~ "group2"))

df.msd <- merge_dss_arg %>%
  group_by(dsssum) %>%
  summarise_at(vars(argsum), list(max = max, mean = mean, sd = sd, median = median, skew = moments::skewness)) %>% 
  as.data.frame() %>%
  mutate(dsssum = as.numeric(dsssum))
df.msd

# Visualize

library(scales)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)

set.seed(12345)

p1 <- ggplot() +
  geom_violin(data = df.count, aes(x = dsssum, y = argsum, fill = dsssum),
              position = position_dodge(0.5), linewidth = 0.5,drop = F) +
  geom_jitter(data = df.count, aes(x = dsssum, y = argsum), size = 0.01, color = "grey") +
 # scale_x_discrete(breaks = 0:25, labels = c(0:25, ">=25")) +
  scale_y_continuous(trans = log2_trans(), breaks = c(1,2,4,8,16,32,64), labels = c(1,2,4,8,16,32,"≥64")) +
  scale_fill_viridis_d(option = 'D')+
# scale_fill_manual(values = c("1" = "#FFEFCC", "2" = "#FECCCC", "3" = "#99CCCC", "4" = "#CCCCFF", "5" = "#3399CC"))+
  labs(title = "", x = "MazEF system count", y = "MRG count", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
      #  axis.line.x = element_line(linewidth = 0.5),
      #  axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(angle = 90,size = 10,hjust = 0.5,vjust = 0.5),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p1
setwd(file.path(GTR_DATA, "groundwater/10mag_arg"))

cairo_pdf(width = 3,height = 3,filename = "2_MazEF_mrg.pdf")
p1
dev.off()

p2 <- ggscatter(df.msd, x = "dsssum", y = "max", add = "reg.line", add.params = list(color = "blue", fill = NA, linetype = "dashed")) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 2, label.y = 60) +
  stat_regline_equation(label.x = 2, label.y = 70) +
  labs(title = "", x = "MazEF system count", y = "Max MRG count", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
       # axis.line.x = element_line(linewidth = 0.5),
       # axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p2

setwd(file.path(GTR_DATA, "groundwater/10mag_arg"))
pdf("2_MazEF_mrg2.pdf", width = 3, height = 3, onefile = F)
p2
#ggdraw() +
#  draw_plot(ggplotify::as.ggplot(p1), x = 0, y = 0, width = 1, height = 1) +
#  draw_plot(ggplotify::as.ggplot(p2), x = 0.57, y = 0.605, width = 0.41, height = 0.38)
dev.off()

#width = 5.4, height = 3.8

##############################################################曲线拟合##########
# Fit the reversed logistic function using nls()
fit <- nls(max ~ a / (1 + exp(b * (dsssum - c))), data = df.msd, 
           start = list(a = 100, b = 0.5, c = 10))  # Provide initial values for a, b, and c
# Try new starting values based on the range of dsssum and max in your data
# fit <- nls(max ~ a / (1 + exp(b * (dsssum - c))), data = df.msd, 
#            start = list(a = max(df.msd$max), b = -1, c = median(df.msd$dsssum)))  # Refined start values

fit <- nls(max ~ a / (1 + exp(b * (dsssum - c))), data = df.msd, 
           start = list(a = 100, b = 2, c = 2))


# Check the summary of the fitted model
summary(fit)
a <- coef(fit)["a"]
b <- coef(fit)["b"]
c <- coef(fit)["c"]

# 假设你已经拟合了 nls 模型 'fit'
fitted_values <- predict(fit)  # 拟合值（预测值）
actual_values <- df.msd$max    # 实际值
residuals <- actual_values - fitted_values  # 残差
SST <- sum((actual_values - mean(actual_values))^2)  # 总平方和
SSE <- sum(residuals^2)  # 残差平方和
R_squared <- 1 - (SSE / SST)
print(paste("R-squared:", R_squared))

# Plotting using ggplot2
p3 <- ggplot(df.msd, aes(x = dsssum, y = max)) +
  geom_point() +  # Plot the points
  stat_function(fun = function(x) coef(fit)["a"] / (1 + exp(coef(fit)["b"] * (x - coef(fit)["c"]))), 
                color = "blue", linetype = "dashed", size =1) +  # Plot the fitted reversed logistic curve
  labs(title = "", x = "MazEF count", y = "Max MRG count", fill = "") +
  # 使用 annotate 手动标注公式
  annotate("text", x = Inf, y = Inf, label = sprintf("y = %.2f / (1 + exp(%.2f * (x - %.2f)))", a, b, c), 
           parse = FALSE, hjust = 1, vjust = 1.5)+
  # 标注 R²
  annotate("text", x = Inf, y = Inf, label = sprintf("R² = %.2f", R_squared), hjust = 1.1, vjust = 3.5)+
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        # axis.line.x = element_line(linewidth = 0.5),
        # axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p3
setwd(file.path(GTR_DATA, "groundwater/10mag_arg"))
pdf("1_MazEF_mrg3.pdf", width = 3, height = 3, onefile = F)
p3
#ggdraw() +
#  draw_plot(ggplotify::as.ggplot(p1), x = 0, y = 0, width = 1, height = 1) +
#  draw_plot(ggplotify::as.ggplot(p2), x = 0.57, y = 0.605, width = 0.41, height = 0.38)
dev.off()

########################批量绘制#################
library(ggplot2)
library(scales)
library(dplyr)
library(ggpubr)
library(cowplot)
plot_logistic_fit <- function(magdss_sum, magarg_sum) {
  # 确保数据框有行名
  magdss_sum$genome <- rownames(magdss_sum)
  magarg_sum$genome <- rownames(magarg_sum)
  
  # 合并两个数据框
  merge_dss_arg <- merge(magdss_sum, magarg_sum, by = 'genome')
  
  # 数据处理
  df.count <- merge_dss_arg %>%
    mutate(dsssum = factor(dsssum, levels = 0:83)) %>%
    mutate(group = case_when(dsssum == 1 ~ "group1", TRUE ~ "group2"))
  
  df.msd <- merge_dss_arg %>%
    group_by(dsssum) %>%
    summarise_at(vars(argsum), list(max = max, mean = mean, sd = sd, median = median, skew = moments::skewness)) %>%
    as.data.frame() %>%
    mutate(dsssum = as.numeric(dsssum))
  
  # 图1: Violin 图
  pall <- ggplot() +
    geom_violin(data = df.count, aes(x = dsssum, y = argsum, fill = dsssum),
                position = position_dodge(0.5), linewidth = 0.5, drop = F) +
    #geom_jitter(data = df.count, aes(x = dsssum, y = argsum), size = 0.01, color = "grey") +
    scale_y_continuous(trans = log2_trans(), breaks = c(1, 2, 4, 8, 16, 32, 64), labels = c(1, 2, 4, 8, 16, 32, "≥64")) +
    scale_fill_viridis_d(option = 'D') +
    labs(x = "MazEF system count", y = "MRG count", fill = "") +
    theme_bw() +
    theme(#plot.margin = unit(c(0.2, 0.2, 0.2, 0), "cm"),
          axis.title.x = element_text(size = 12, face = "plain"),
          axis.title.y = element_text(size = 12, face = "plain"),
          axis.text.x = element_text(angle = 90, size = 10, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 10),
          legend.position = "none")
  
  # Logistic 模型拟合
  fit <- nls(max ~ a / (1 + exp(b * (dsssum - c))), data = df.msd, 
             start = list(a = 100, b = 2, c = 2))
  
  # 提取拟合参数
  a <- coef(fit)["a"]
  b <- coef(fit)["b"]
  c <- coef(fit)["c"]
  
  # 计算 R²
  fitted_values <- predict(fit)
  actual_values <- df.msd$max
  residuals <- actual_values - fitted_values
  SST <- sum((actual_values - mean(actual_values))^2)
  SSE <- sum(residuals^2)
  R_squared <- 1 - (SSE / SST)
  
  # 图2: Logistic 拟合曲线
  pall2 <- ggplot(df.msd, aes(x = dsssum, y = max)) +
    geom_point() +  # 绘制散点
    stat_function(fun = function(x) a / (1 + exp(b * (x - c))), color = "blue", linetype = "dashed", size = 1) +  # 拟合曲线
    labs(x = "MazEF count", y = "Max MRG count") +
    annotate("text", x = Inf, y = Inf, label = sprintf("y = %.2f / (1 + exp(%.2f * (x - %.2f)))", a, b, c), 
             parse = FALSE, hjust = 1, vjust = 1.5) +
    annotate("text", x = Inf, y = Inf, label = sprintf("R² = %.2f", R_squared), hjust = 1.1, vjust = 3.5) +
    theme_bw() +
    theme(#plot.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.title.x = element_text(size = 12, face = "plain"),
          axis.title.y = element_text(size = 12, face = "plain"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.position = "none")
  
  # 返回两个图
  return(list(violin_plot = pall, logistic_plot = pall2))
}

magdss_sum <- dss %>% rowSums() %>% as.data.frame() %>% setNames("dsssum")#求和每个mags中的dss
magdss_sum <- magdss_sum %>%
  filter(rownames(.) %in% rownames(magarg_sum))
plots_all <- plot_logistic_fit(magdss_sum, magarg_sum)
plots_all$violin_plot
plots_all$logistic_plot

magdss_sum <- select(dss,-RM, -SoFIC, -CasFinder, -AbiE, -MazEF)
magdss_sum <- magdss_sum %>% rowSums() %>% as.data.frame() %>% setNames("dsssum")#求和每个mags中的dss
magdss_sum <- magdss_sum %>%
  filter(rownames(.) %in% rownames(magarg_sum))
plots_pan <- plot_logistic_fit(magdss_sum, magarg_sum)
plots_pan$violin_plot
plots_pan$logistic_plot

magdss_sum <- select(dss,RM, SoFIC, CasFinder, AbiE, MazEF)
magdss_sum <- magdss_sum %>% rowSums() %>% as.data.frame() %>% setNames("dsssum")#求和每个mags中的dss
magdss_sum <- magdss_sum %>%
  filter(rownames(.) %in% rownames(magarg_sum))
plots_corner <- plot_logistic_fit(magdss_sum, magarg_sum)
plots_corner$violin_plot
plots_corner$logistic_plot

magdss_sum <- data.frame(genome = rownames(dss) , dsssum = dss$MazEF)#分别对五个主要的DS进行探究：
rownames(magdss_sum) <- magdss_sum$genome
magdss_sum <- select(magdss_sum, dsssum)
plots_MazEF <- plot_logistic_fit(magdss_sum, magarg_sum)
plots_MazEF$violin_plot
plots_MazEF$logistic_plot

magdss_sum <- data.frame(genome = rownames(dss) , dsssum = dss$RM)#分别对五个主要的DS进行探究：
rownames(magdss_sum) <- magdss_sum$genome
magdss_sum <- select(magdss_sum, dsssum)
plots_RM <- plot_logistic_fit(magdss_sum, magarg_sum)
plots_RM$violin_plot
plots_RM$logistic_plot

magdss_sum <- data.frame(genome = rownames(dss) , dsssum = dss$SoFIC)#分别对五个主要的DS进行探究：
rownames(magdss_sum) <- magdss_sum$genome
magdss_sum <- select(magdss_sum, dsssum)
plots_SoFIC <- plot_logistic_fit(magdss_sum, magarg_sum)
plots_SoFIC$violin_plot
plots_SoFIC$logistic_plot

magdss_sum <- data.frame(genome = rownames(dss) , dsssum = dss$Cas)#分别对五个主要的DS进行探究：
rownames(magdss_sum) <- magdss_sum$genome
magdss_sum <- select(magdss_sum, dsssum)
plots_Cas <- plot_logistic_fit(magdss_sum, magarg_sum)
plots_Cas$violin_plot
plots_Cas$logistic_plot

magdss_sum <- data.frame(genome = rownames(dss) , dsssum = dss$AbiE)#分别对五个主要的DS进行探究：
rownames(magdss_sum) <- magdss_sum$genome
magdss_sum <- select(magdss_sum, dsssum)
plot_logistic_fit <- function(magdss_sum, magarg_sum) {
  # 确保数据框有行名
  magdss_sum$genome <- rownames(magdss_sum)
  magarg_sum$genome <- rownames(magarg_sum)
  
  # 合并两个数据框
  merge_dss_arg <- merge(magdss_sum, magarg_sum, by = 'genome')
  
  # 数据处理
  df.count <- merge_dss_arg %>%
    mutate(dsssum = factor(dsssum, levels = 0:83)) %>%
    mutate(group = case_when(dsssum == 1 ~ "group1", TRUE ~ "group2"))
  
  df.msd <- merge_dss_arg %>%
    group_by(dsssum) %>%
    summarise_at(vars(argsum), list(max = max, mean = mean, sd = sd, median = median, skew = moments::skewness)) %>%
    as.data.frame() %>%
    mutate(dsssum = as.numeric(dsssum))
  
  # 图1: Violin 图
  pall <- ggplot() +
    geom_violin(data = df.count, aes(x = dsssum, y = argsum, fill = dsssum),
                position = position_dodge(0.5), linewidth = 0.5, drop = F) +
    #geom_jitter(data = df.count, aes(x = dsssum, y = argsum), size = 0.01, color = "grey") +
    scale_y_continuous(trans = log2_trans(), breaks = c(1, 2, 4, 8, 16, 32, 64), labels = c(1, 2, 4, 8, 16, 32, "≥64")) +
    scale_fill_viridis_d(option = 'D') +
    labs(x = "MazEF system count", y = "MRG count", fill = "") +
    theme_bw() +
    theme(#plot.margin = unit(c(0.2, 0.2, 0.2, 0), "cm"),
      axis.title.x = element_text(size = 12, face = "plain"),
      axis.title.y = element_text(size = 12, face = "plain"),
      axis.text.x = element_text(angle = 90, size = 10, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(size = 10),
      legend.position = "none")
  
  # # Logistic 模型拟合
  # fit <- nls(max ~ a / (1 + exp(b * (dsssum - c))), data = df.msd, 
  #            start = list(a = 120, b = 0.9, c = 0.2),
  #            control = nls.control(maxiter = 10000, minFactor = 1e-10)
  #            )
  # 
  # # 提取拟合参数
  # a <- coef(fit)["a"]
  # b <- coef(fit)["b"]
  # c <- coef(fit)["c"]
  # 
  # # 计算 R²
  # fitted_values <- predict(fit)
  # actual_values <- df.msd$max
  # residuals <- actual_values - fitted_values
  # SST <- sum((actual_values - mean(actual_values))^2)
  # SSE <- sum(residuals^2)
  # R_squared <- 1 - (SSE / SST)
  
  # 图2: Logistic 拟合曲线
  pall2 <- ggscatter(df.msd, x = "dsssum", y = "max", add = "reg.line", add.params = list(color = "blue", fill = NA, linetype = "dashed")) +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 2, label.y = 60) +
    stat_regline_equation(label.x = 2, label.y = 70) +
    labs(title = "", x = "MazEF system count", y = "Max MRG count", fill = "") +
    labs(x = "MazEF count", y = "Max MRG count") +
    # annotate("text", x = Inf, y = Inf, label = sprintf("y = %.2f / (1 + exp(%.2f * (x - %.2f)))", a, b, c), 
    #          parse = FALSE, hjust = 1, vjust = 1.5) +
    # annotate("text", x = Inf, y = Inf, label = sprintf("R² = %.2f", R_squared), hjust = 1.1, vjust = 3.5) +
    theme_bw() +
    theme(#plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.title.x = element_text(size = 12, face = "plain"),
      axis.title.y = element_text(size = 12, face = "plain"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "none")
  
  
  
  # 返回两个图
  return(list(violin_plot = pall, logistic_plot = pall2))
}
plots_AbiE <- plot_logistic_fit(magdss_sum, magarg_sum)
plots_AbiE$violin_plot
plots_AbiE$logistic_plot

library(Cairo)
library(cowplot)

row1 <- plot_grid(plots_all$violin_plot, plots_all$logistic_plot,
               
               labels = LETTERS[1:2], nrow = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
row1
row2 <- plot_grid(
                  plots_corner$violin_plot, plots_corner$logistic_plot, plots_pan$violin_plot, plots_pan$logistic_plot,
                  labels = LETTERS[3:6], nrow = 1
                  #rel_heights = c(1, 1, 1.5),
                  # label_x = 0.09,            # 标签x坐标（1为最右侧）
                  #  label_y = 0.9             # 标签y坐标（1为最上方）
)
row2
row3 <- plot_grid(plots_RM$violin_plot,plots_SoFIC$violin_plot,plots_Cas$violin_plot,plots_AbiE$violin_plot,plots_MazEF$violin_plot,
                  plots_RM$logistic_plot,plots_SoFIC$logistic_plot,plots_Cas$logistic_plot,plots_AbiE$logistic_plot,plots_MazEF$logistic_plot,
                  ncol = 5)
row3

# 组合所有行
final_plot <- plot_grid(row1, row2, row3, nrow = 3, rel_heights = c(1, 1, 1.5))
final_plot

setwd(file.path(GTR_DATA, "groundwater/10mag_arg"))
cairo_pdf(width = 10,height = 8,filename = "2_ds_mrg2.pdf")
final_plot
dev.off()

##################################################################
#dgs与arg
ggplot(merge_dss_arg_dgs, aes(x = dgssum, y = argsum)) +
  geom_point() +  # 添加散点图
  geom_smooth(method = "lm", se = TRUE,formula = 'y ~ 1/x') +  # 添加线性回归线，不显示置信区间
  labs(title = "Linear Regression Model", x = "X", y = "Y")  # 添加标题和坐标轴标签
#没有明显的规律。得到的可能信息：DGs与ARGs的关系是互相排斥的。
#它们虽然会给细菌带来生存优势，但是当其中一个数量很多时，另一个数量就不会很多，以减小细菌的负荷。


#增加物种信息进一步探究
setwd(file.path(GTR_DATA, "groundwater/02去冗余信息"))
mags_tax <- read.delim2('gtdbtk.bac120.summary.tsv')
# 假设 df 是你的数据框，MAGs 是其中的一列
mags_tax$user_genome <- paste(mags_tax$user_genome, ".fasta", sep="")
mags_tax <- mags_tax[mags_tax$user_genome %in% rownames(dss),]

#处理物种信息
taxonoxy <- function(tax){
  tax <- tax[,1:2]
  max_cols <- max(sapply(strsplit(tax$classification, ";"), length))# 确定最大的拆分列数
  col_names <- paste0("col", 1:max_cols)# 生成列名
  library(tidyr)
  taxonomy <- tax %>%
    separate(classification, into = col_names, sep = ";", fill = "right")# 拆分字符串到列
}

mags_tax <- taxonoxy(mags_tax)

merge_dss_arg_dgs_tax <- merge(merge_dss_arg_dgs,mags_tax,by.x='genome',by.y = 'user_genome')

p2 <- ggplot(merge_dss_arg_dgs_tax, aes(x = dsssum, y = argsum)) +
  geom_point() +  # 添加散点图
  geom_smooth(method = "lm", se = TRUE) +  # 添加线性回归线，不显示置信区间
  labs(title = NULL, x = "DSs number (per MAG)", y = "ARGs number (per MAG)")+  # 添加标题和坐标轴标签
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.8,0.7),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p2
