# =============================================================================
# aquifer-defensome -- main figures / Defensome of CPR, Pseudomonadota and other phyla
# =============================================================================
# DS number across CPR / other phyla / Pseudomonadota, the phylum-level DS
# frequency heatmap and the per-phylum comparisons.
#
# Manuscript      : Fig. 2
# Main outputs    : 3phyla_dss.pdf
# Original script : R codes/Fig2.R
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


#Fig2A,B####################
p1 <- ggplot(ds_density_tax, aes(x = CPR, y = number)) +
  geom_violin(trim = FALSE, aes(fill = CPR),show.legend = F) + # 绘制小提琴图
  geom_boxplot(aes(fill = CPR),width = 0.2,outliers = FALSE, staplewidth = 0.5) +  # 隐藏箱线图中的异常值点
  #geom_jitter(width = 0.2, alpha = 0.7)+  # 显示异常值点
  # stat_compare_means(comparisons = list(c("1", "2"),c("SS", "WA"),c("SS", "WS"),
  #                                       c("SA", "WA"),c("SA", "WS"),
  #                                       c("WA", "WS")
  # ),
  # method = "wilcox.test", test = "p.adjust",label = "p.signif")+
  stat_pvalue_manual(md_p_val1,label = "p.adj.signif",hide.ns=T,
                     tip.length = 0,label.size = 5,color="black")+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  #stat_summary(fun=mean, geom="line", aes(group=CPR), col="black")+
  scale_fill_manual(values = c("CPR"="#FF7F0E","Other phyla" = "#ADD9EE","Pseudomonadota"="#7892C7"
  )) +
  #ylim(0,15)+
  #scale_y_continuous(limits = c(0,10))+
  coord_cartesian(ylim = c(0,15))+
  labs(x=NULL,y='DSs number (per MAG)')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p1

md_p_val2 <- ds_density_tax %>% 
  wilcox_test(density ~ CPR) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="CPR",dodge=0.5)
md_p_val2$p.adj

md_p_val2$y.position <- c(8,9,10)
p2 <- ggplot(ds_density_tax, aes(x = CPR, y = density)) +
  geom_violin(trim = FALSE, aes(fill = CPR),show.legend = F) + # 绘制小提琴图
  geom_boxplot(aes(fill = CPR),width = 0.2,outliers = FALSE, staplewidth = 0.5) +  # 隐藏箱线图中的异常值点
  #geom_jitter(width = 0.2, alpha = 0.7)+  # 显示异常值点
  # stat_compare_means(comparisons = list(c("1", "2"),c("SS", "WA"),c("SS", "WS"),
  #                                       c("SA", "WA"),c("SA", "WS"),
  #                                       c("WA", "WS")
  # ),
  # method = "wilcox.test", test = "p.adjust",label = "p.signif")+
  stat_pvalue_manual(md_p_val2,label = "p.adj.signif",hide.ns=T,
                     tip.length = 0,label.size = 5,color="black")+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  #stat_summary(fun=mean, geom="line", aes(group=CPR), col="black")+
  scale_fill_manual(values = c("CPR"="#FF7F0E","Other phyla" = "#ADD9EE","Pseudomonadota"="#7892C7"
  )) +
  #ylim(0,15)+
  #scale_y_continuous(limits = c(0,10))+
  coord_cartesian(ylim = c(0,10))+
  labs(x=NULL,y='DSs density (per MAG per Mb)')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p2

#Fig2C,D,E##############
dsr_phayla_density <- ggplot(dsr, aes(x = col2, y = density_mean)) +
  geom_bar(stat = "identity", fill = "#43A3EF") +
  geom_errorbar(
    aes(ymin = density_mean - 1.96*density_se, ymax = density_mean + 1.96*density_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = density_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 3  # 点大小
  )+
  geom_hline(yintercept = mean(dsr$density_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DSs richness (per MAG per Mb)") +
  scale_y_continuous(limits = c(0,1.5),expand = c(0,0))+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # 横轴标签旋转
    #axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dsr_phayla_density

dsr_phayla_richness <- ggplot(dsr, aes(x = col2, y = richness_mean)) +
  geom_bar(stat = "identity", fill = "#EF767B") +
  geom_errorbar(
    aes(ymin = richness_mean - 1.96*richness_se, ymax = richness_mean + 1.96*richness_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = richness_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 3  # 点大小
  )+
  geom_hline(yintercept = mean(dsr$richness_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DSs richness (per MAG)") +
  scale_y_continuous(limits = c(0,4.5),expand = c(0,0), label = function(x) sprintf("%.1f", x) )+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # 横轴标签旋转
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dsr_phayla_richness

dsr_phayla_length <- ggplot(dsr, aes(x = col2, y = length_Mb_mean)) +
  geom_bar(stat = "identity", fill = "#F9E7A7") +
  geom_errorbar(
    aes(ymin = length_Mb_mean - 1.96*length_Mb_se, ymax = length_Mb_mean + 1.96*length_Mb_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = length_Mb_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 3  # 点大小
  )+
  geom_hline(yintercept = mean(dsr$length_Mb_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "MAG size (Mb)") +
  scale_y_continuous(limits = c(0,6.5),expand = c(0,0), label = function(x) sprintf("%.1f", x) )+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # 横轴标签旋转
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dsr_phayla_length

library(cowplot)
g <- plot_grid(dsr_phayla_length, dsr_phayla_richness,dsr_phayla_density,
               labels = LETTERS[3:5], nrow = 3,rel_heights = c(1, 1, 1.5),
               label_x = 0.09,            # 标签x坐标（1为最右侧）
               label_y = 0.9             # 标签y坐标（1为最上方）
)
g

#Fig2F###################
#1.1读取MAGs上的DSs数据
setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
#1.读取数据表###
#读取所有27578个mags的数据
setwd(file.path(GTR_DATA, "groundwater/02去冗余信息"))
allmags <- read.delim2('Wdb_ANI99.csv',sep = ',',
                       quote = "",header = T
) 

#读取含有dss的mags的数据
setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
dss <- read.delim2('ds_counts_matrix.csv',sep = ',',row.names = 1,
                   quote = "",header = T,check.names = F
) 
#将不含dss的mags也合并到一起
# 找出 allmags 中在 dss 行名中不存在的 genome 值
not_in_dss <- allmags$genome[!allmags$genome %in% rownames(dss)]
# 使用不在 dss 行名中的 genome 值创建新的数据框
new_rows <- data.frame(
  genome = not_in_dss,  # 重复 genome 值
  stringsAsFactors = FALSE
)
# 为 new_rows 添加与 dss 相同的列，并将这些列的值填充为 0
col_names <- setdiff(names(dss), "genome")  # 获取 dss 中除了 genome 之外的列名
new_rows[col_names] <- lapply(col_names, function(col_name) rep(0, nrow(dss)))
rownames(new_rows) <- new_rows$genome
new_rows <- new_rows[,-1]
#合并new_rows与dss
dss <- rbind(dss,new_rows)


#1.2读取MAGs的物种注释信息
setwd(file.path(GTR_DATA, "groundwater/02去冗余信息"))

# mags_tax <- readxl::read_xlsx('bac.gtdb.r207.xlsx')
# # 假设 df 是你的数据框，MAGs 是其中的一列
# mags_tax$MAGs <- paste(mags_tax$MAGs, ".fasta", sep="")
# mags_tax <- mags_tax[mags_tax$MAGs %in% rownames(dss),]
# 
# #处理物种信息
# taxonoxy <- function(tax){
#   tax <- tax[,1:2]
#   max_cols <- max(sapply(strsplit(tax$Taxonomy, ";"), length))# 确定最大的拆分列数
#   col_names <- paste0("col", 1:max_cols)# 生成列名
#   library(tidyr)
#   taxonomy <- tax %>%
#     separate(Taxonomy, into = col_names, sep = ";", fill = "right")# 拆分字符串到列
# }
# 
# mags_tax <- taxonoxy(mags_tax)

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



#1.3合并物种信息与DSs信息
dss$genome <- rownames(dss)
tax_dss <- merge(dss,mags_tax,by.x = 'genome',by.y = 'user_genome')

#在19847个MAGs中找到

#2. 绘制各个门MAGs的DSs分布热图###

#删除未注释到门的行,剩下27564个MAGs
tax_dss <- tax_dss[tax_dss$col2 != 'p__',]

#2.1 数据处理
library(dplyr)
#对tax_dss的col2列进行处理，删除_后面的字母
# 假设 tax_dss 是你的数据框，col2 是你想要处理的列
library(stringr)
# 使用 sub 函数替换掉第三个 _ 及其后面的内容
tax_dss$col2 <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", tax_dss$col2)

#统计col2列每个元素的数量
element_counts <- tax_dss %>%
  count(col2)
#setwd(file.path(GTR_DATA, "groundwater/01样品信息"))
#write.csv(element_counts,"phyla_count.csv")

phyla100 <- element_counts[element_counts$n>=100,]

#分组统计
result <- function(df) {
  output <- df %>%
    group_by(col2) %>%
    summarise(across(where(is.numeric), ~mean(. != 0, na.rm = TRUE)))
  return(output)
}#这个函数和下面的函数的计算结果相同

result <- function(df) {
  output <- df %>%
    group_by(col2) %>%
    summarise(across(where(is.numeric), ~sum(. != 0, na.rm = TRUE) / sum(!is.na(.), na.rm = TRUE)))
  return(output)
}

result_data <- result(tax_dss) %>% as.data.frame()



#2.2 过滤数据
library(dplyr)
#筛选>100的mags
result_data <- result_data[result_data$col2 %in% phyla100$col2,]

rownames(result_data) <- result_data$col2
result_data <- result_data[,-1]

rownames(result_data) <- gsub("^p__", "", rownames(result_data))

# 假设 df 是你的数据框
# 计算每一行非零数值的比例
# df <- result_data[,-1] %>%
#   mutate(zero_proportion = rowSums(. != 0) / ncol(.)) %>%
#   # 过滤出非零数值比例大于1%的行
#   filter(rowSums(. != 0) / ncol(.) > 0.1)

#过滤列#
#删除全为0的列
result_data <- result_data[, !colSums(result_data) == 0]
#筛选主要的DSs
# 计算每一列非0数字所占的百分比并输出到新数据框
percent_nonzero <- colSums(dss != 0) / nrow(dss) * 100
new_df <- data.frame(Column = names(percent_nonzero), Percentage = percent_nonzero)
# 按照值的大小对数据集进行降序排列
ds_persentage <- new_df[order(-new_df$Percentage), ]
ds_persentage <- ds_persentage[ds_persentage$Percentage>=1,]
result_data <- result_data[,colnames(result_data) %in% ds_persentage$Column]
# # 计算每列非零值的比例
# non_zero_ratios <- sapply(result_data, function(x) sum(x != 0) / length(x))
# # 筛选出比例大于1%的列
# columns_to_keep <- names(non_zero_ratios)[non_zero_ratios > 0.1]
# # 只保留这些列
# df_filtered <- result_data[, columns_to_keep]

#2.3 绘制
#绘制热图
library(pheatmap)

# 设置颜色
#colors <- colorRampPalette(c('white', colorRampPalette('#F9B14F', '#8F2731')))(1000)
colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p <- pheatmap(t(result_data),
              #color = colorRampPalette(colors)(1000),
              #color = colorRampPalette(c('white','lightpink','red'))(1000),
              color = colors,
              cluster_cols = F,
              cluster_rows = F,
)
setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
write.csv(result_data,"Fig2F.csv",quote = F)
ggplot2::ggsave("3phyla_dss.pdf", plot = p, width = 4.5, height = 8, units = "in")
