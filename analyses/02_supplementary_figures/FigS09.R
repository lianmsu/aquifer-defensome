# =============================================================================
# aquifer-defensome -- supplementary figures / Plaque assay experiment
# =============================================================================
# Quantification and plotting of the three plaque-assay replicates (wet-lab
# validation).
#
# Manuscript      : Fig. S9
# Main outputs    : plaque_count_plot_2.pdf
# Original script : R codes/Fig.S9.R
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
GTR_REV <- normalizePath(Sys.getenv("GTR_REV", unset = "."), mustWork = FALSE)


library(dplyr)
library(reshape2)

setwd(file.path(GTR_REV, "01defensome/NC_revision/2experiments/plaque_results_20250508"))

plaque1 <- readxl::read_excel("plaque_1.xlsx") %>% as.data.frame()
# 2. 将第一列设为行名，并删除该列
rownames(plaque1) <- plaque1[[1]]  # 用第一列设置行名
colnames(plaque1)[1] <- "DS"
#plaque1 <- plaque1[-1]            # 移除数据框中的第一列
plaque1 <- melt(plaque1)

plaque2 <- readxl::read_excel("plaque_2.xlsx") %>% as.data.frame()
# 2. 将第一列设为行名，并删除该列
rownames(plaque2) <- plaque2[[1]]  # 用第一列设置行名
colnames(plaque2)[1] <- "DS"
#plaque1 <- plaque1[-1]            # 移除数据框中的第一列
plaque2 <- melt(plaque2)

plaque3 <- readxl::read_excel("plaque_3.xlsx") %>% as.data.frame()
# 2. 将第一列设为行名，并删除该列
rownames(plaque3) <- plaque3[[1]]  # 用第一列设置行名
colnames(plaque3)[1] <- "DS"
#plaque1 <- plaque1[-1]            # 移除数据框中的第一列
plaque3 <- melt(plaque3)

plaque_res <- rbind(plaque1,plaque2,plaque3)

plaque_res$PFU <- log10(plaque_res$value*200)

write.csv(plaque_res,"figs9b.csv")


library(ggpubr)
ggplot(plaque_res, aes(x = DS, y = PFU, fill = variable)) +
  geom_boxplot() +
  stat_compare_means(aes(group = variable), 
                     method = "t.test", 
                     label = "p.signif") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(ggplot2)
library(dplyr)

# 1. 准备数据
plaque_res <- plaque_res %>%
  mutate(
    # 设置DS的因子水平（指定顺序）
    DS = factor(DS, levels = c("pQE60", "RM_IIG_close", "RM_IIG_distant", 
                               "SoFIC_close", "SoFIC_distant", "pETDuet1", "CBASS_I")),
    # 创建分面变量
    facet_group = ifelse(DS %in% c("pETDuet1", "CBASS_I"), "Group2", "Group1")
  )


ggplot(plaque_res, aes(x = DS, y = PFU, fill = variable)) +
  geom_boxplot() +
  # 分面设置（Group1 和 Group2 分开）
  facet_grid(. ~ facet_group, scales = "free_x", space = "free_x") +
  # 美化图形
  labs(
    title = "Plaque Assay Results by DS and Variable",
    x = "DS",
    y = "PFU (log10)",
    fill = "Variable"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转 x 轴标签
    strip.background = element_blank(),  # 隐藏分面背景
    strip.text = element_blank()  # 隐藏分面标签（因为图形已清晰）
  ) +
  scale_fill_brewer(palette = "Set2")  # 使用更佳的颜色


library(dplyr)
library(ggplot2)

# 设定 DS 的顺序
ds_order <- c("pQE60", "RM_IIG_close", "RM_IIG_distant", 
              "SoFIC_close", "SoFIC_distant", "pETDuet1", "CBASS_I")

# 创建分面列：横向分面（Group1/Group2）和纵向分面（T1/T7）
plaque_res <- plaque_res %>%
  mutate(
    DS = factor(DS, levels = ds_order),  # 强制 DS 按指定顺序排列
    facet_horizontal = ifelse(DS %in% c("pETDuet1", "CBASS_I"), "Group2", "Group1"),  # 横向分面
    facet_vertical = variable  # 纵向分面（T1/T7）
  ) %>%
  filter(facet_vertical %in% c("T1", "T7"))  # 仅保留 T1 和 T7（按需调整）

plaque_count_plot <- 
ggplot(plaque_res, aes(x = DS, y = PFU, fill = variable)) +
  geom_boxplot(width = 0.7, alpha = 0.8) +  # 调整箱线宽度和透明度
  # 叠加数据点（抖动避免重叠）
  geom_jitter(
    width = 0.2,      # 水平抖动范围
    height = 0,       # 垂直不抖动（保持y值精确）
    aes(color = variable),  # 按variable着色
    size = 2,         # 点大小
    alpha = 0.8       # 透明度
  ) +
  # 双向分面：纵向 ~ 横向
  facet_grid(
    rows = vars(facet_vertical),  # 纵向：T1 在上，T7 在下
    cols = vars(facet_horizontal),  # 横向：Group1（前5个DS）左，Group2（后2个DS）右
    scales = "free_x",  # x 轴独立缩放
    space = "free_x"    # 避免空白占位
  ) +
  # 美化图形
  labs(
    title = "",
    x = "",
    y = "PFU/ml (log10)",
    fill = "Variable"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # 旋转 x 轴标签
    strip.background = element_rect(fill = "white"),  # 分面标签背景
    strip.text = element_text(face = "bold"),        # 分面标签加粗
    panel.spacing = unit(1, "lines"),                 # 增加分面间距
    legend.position = "na"
  ) +
  scale_fill_brewer(palette = "Set2")+  # 使用更佳的颜色
  scale_color_brewer(palette = "Set2") 
# scale_fill_manual(values = c("T1" = "#4E79A7", "T7" = "#F28E2B"))  # 自定义颜色
plaque_count_plot


setwd(file.path(GTR_REV, "01defensome/NC_revision/2experiments/plaque_results_20250508"))
ggsave("plaque_count_plot_2.pdf",plaque_count_plot ,width = 6,height = 5, units = "in")



wilcox_results  <- plaque_res %>% group_by(variable) %>%  # 按DS分组（可选：如果需要在不同DS下分别检验）
  wilcox_test(
  PFU ~ DS,        # 公式：PFU按variable分组比较
  paired = TRUE,         # 配对检验
 # id = batch,            # 配对依据（同一batch内的样本配对）
  p.adjust.method = "BH" # Benjamini-Hochberg校正
)


library(rstatix)
# Perform paired t-tests within each variable group
t_results <- plaque_res %>%
  group_by(variable) %>%
  t_test(
    PFU ~ DS,
    paired = TRUE,
    detailed = TRUE,  # Get more detailed output
    p.adjust.method = "BH"
  ) %>%
  add_significance()

# View results
print(t_results, n = Inf)



c1 <- c(30.5,26.5,25)
c2 <- c(11.5,10,9.5)
res <- t.test(c1, c2, paired = T)  # 要求c1和c2长度相同

pQE60 <- plaque_res[plaque_res$DS=="pQE60" & plaque_res$variable=="T1",]$PFU
SoFIC_close <- plaque_res[plaque_res$DS=="SoFIC_close" & plaque_res$variable=="T1",]$PFU
SoFIC_distant <- plaque_res[plaque_res$DS=="SoFIC_distant" & plaque_res$variable=="T1",]$PFU
RM_IIG_close <- plaque_res[plaque_res$DS=="RM_IIG_close" & plaque_res$variable=="T1",]$PFU
RM_IIG_distant <- plaque_res[plaque_res$DS=="RM_IIG_distant" & plaque_res$variable=="T1",]$PFU
pETDuet1 <- plaque_res[plaque_res$DS=="pETDuet1" & plaque_res$variable=="T1",]$PFU
CBASS_I <- plaque_res[plaque_res$DS=="CBASS_I" & plaque_res$variable=="T1",]$PFU


# t.test(pQE60, SoFIC_close, paired = T) #0.01344
# t.test(pQE60, SoFIC_distant, paired = T) #0.2568
# t.test(pQE60, RM_IIG_close, paired = T) #6.468e-06
# t.test(pQE60, RM_IIG_distant, paired = T) #0.1835
# t.test(pETDuet1, CBASS_I, paired = T) #0.08596

