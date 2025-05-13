#对

#1. 读取数据####
#1.1读取MAGs上的DSs数据
setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
#1.读取数据表###
#读取所有27578个mags的数据
setwd('D:/wenjian/research/defensome/groundwater/02去冗余信息')
allmags <- read.delim2('Wdb_ANI99.csv',sep = ',',
                       quote = "",header = T
) 

#读取含有dss的mags的数据
setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
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
setwd('D:/wenjian/research/defensome/groundwater/02去冗余信息')

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

#2. 绘制各个门MAGs的DSs分布热图####

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
#setwd('D:/wenjian/research/defensome/groundwater/01样品信息')
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
setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
write.csv(result_data,"Fig2F.csv",quote = F)
ggplot2::ggsave("3phyla_dss.pdf", plot = p, width = 4.5, height = 8, units = "in")


#3. 注释MAGs属级别的DSs分布热图####
#删除未注释到门的行,剩余21356个MAGs
tax_dss <- tax_dss[tax_dss$col6 != 'g__',]

#3.1 数据处理
library(dplyr)
#对tax_dss的col2列进行处理，删除_后面的字母
# 假设 tax_dss 是你的数据框，col2 是你想要处理的列
library(stringr)
# 使用 sub 函数替换掉第三个 _ 及其后面的内容
tax_dss$col6 <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", tax_dss$col6)


#分组统计
result <- function(df) {
  output <- df %>%
    group_by(col6) %>%
    summarise(across(where(is.numeric), ~mean(. != 0, na.rm = TRUE)))
  return(output)
}
result_data <- result(tax_dss) %>% as.data.frame()

#统计col2列每个元素的数量
genus_counts <- tax_dss %>%
  count(col6)
# setwd('D:/wenjian/research/defensome/groundwater/01样品信息')
# write.csv(genus_counts,"genera_count.csv")

genera50 <- genus_counts[genus_counts$n>=50,]
#3.2 过滤数据
#过滤行
#筛选>50的mags,有63个
result_data <- result_data[result_data$col6 %in% genera50$col6,]

rownames(result_data) <- result_data$col6
result_data <- result_data[,-1]

rownames(result_data) <- gsub("^g__", "", rownames(result_data))

#过滤列
#删除全为0的列
result_data <- result_data[, !colSums(result_data) == 0]
#提取主要的DSs
result_data <- result_data[,colnames(result_data) %in% ds_persentage$Column]

#3.3 绘制
#绘制热图
library(pheatmap)

# 创建颜色梯度，0为白色，其余为蓝色到红色的渐变
colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p2 <- pheatmap(t(result_data),
              #color = colorRampPalette(colors)(1000),
              #color = colorRampPalette(c('white','lightpink','red'))(1000),
              color = colors,
              cluster_cols = F,
              cluster_rows = F,
              )
p2
setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
ggplot2::ggsave("3genera_dss.pdf", plot = p2, width = 9, height = 5, units = "in")


#4. 主要细菌类别中DSs及其richness####
#21个主要门中DSs的多样性####
#基于代码第128行得到的result_data进行后续分析
row_sums_result <- apply(result_data, 1, function(row) sum(row > 0))
# Acidobacteriota    Actinomycetota         Bacillota      Bacteroidota  Bdellovibrionota 
# 80                90               103               104                77 
# Campylobacterota     Chloroflexota   Cyanobacteriota  Desulfobacterota   Elusimicrobiota 
# 84                64                62               108                44 
# Gemmatimonadota Methylomirabilota       Myxococcota      Nitrospinota      Nitrospirota 
# 34                54                68                46                80 
# Omnitrophota   Patescibacteria   Planctomycetota    Pseudomonadota          UBA10199 
# 43                87                83               124                31 
# Verrucomicrobiota 
# 82 
#绘制简单的柱状图
library(ggplot2)
row_sums_result_df <- as.data.frame(row_sums_result)
names(row_sums_result_df) <- c("Count")
dss_richness <- ggplot(row_sums_result_df, aes(x = rownames(row_sums_result_df), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DSs richness") +
  scale_y_continuous(limits = c(0,130),expand = c(0,0))+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.8,0.7),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
  
dss_richness
setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
ggsave("3phyla_dss_richness.pdf", plot = dss_richness, width = 5, height = 5, units = "in")

#21个门中DSs richness在每个MAGs中的平均值####
magdsr <- data.frame(genome=rownames(dss),dsr=rowSums(dss[,-140] > 0))

setwd('D:/wenjian/research/defensome/groundwater/03mag_info')
magsize <- read.csv('genome_stats.csv')
dsr_size <- merge(magdsr,magsize,by.x = 'genome',by.y = 'genome')

dsr_density <- data.frame(genome=dsr_size$genome,
                         density=dsr_size$dsr/dsr_size$size*1000000,
                         richness=dsr_size$dsr,
                         length_Mb=dsr_size$size/1000000
                         )

dsr_density <- merge(dsr_density,mags_tax,by.x = 'genome',by.y = 'user_genome')

# 使用 sub 函数替换掉第三个 _ 及其后面的内容
dsr_density$col2 <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", dsr_density$col2)

# 使用mutate和if_else增加新的CPR列
dsr_density <- dsr_density %>%
  mutate(CPR = case_when(col2 == "p__Patescibacteria" ~ "CPR",
                         col2 == "p__Pseudomonadota" ~ "Pseudomonadota",
                         TRUE ~ "Other phyla"
  ))
##################################
#计算显著性差异
library(rstatix)
library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means

md_p_val3 <- dsr_density %>% 
  wilcox_test(richness ~ CPR) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="CPR",dodge=0.5)
md_p_val3$p.adj

md_p_val3$y.position <- c(13,14,15)
p3 <- ggplot(dsr_density, aes(x = CPR, y = richness)) +
  geom_violin(trim = FALSE, aes(fill = CPR),show.legend = F) + # 绘制小提琴图
  geom_boxplot(aes(fill = CPR),width = 0.2,outliers = FALSE, staplewidth = 0.5) +  # 隐藏箱线图中的异常值点
  #geom_jitter(width = 0.2, alpha = 0.7)+  # 显示异常值点
  # stat_compare_means(comparisons = list(c("1", "2"),c("SS", "WA"),c("SS", "WS"),
  #                                       c("SA", "WA"),c("SA", "WS"),
  #                                       c("WA", "WS")
  # ),
  # method = "wilcox.test", test = "p.adjust",label = "p.signif")+
  stat_pvalue_manual(md_p_val3,label = "p.adj.signif",hide.ns=T,
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
p3
##########################

# result3 <- function(df) {
#   output <- df %>%
#     group_by(col2) %>%
#     summarise(across(where(is.numeric), ~mean(., na.rm = TRUE)))
#   return(output)
# }

result3 <- function(df) {
  output <- df %>%
    group_by(col2) %>%
    summarise(
      across(
        where(is.numeric),
        list(
          mean = ~mean(., na.rm = TRUE),
          sd = ~sd(., na.rm = TRUE),  # 添加标准差的计算
          se = ~sd(., na.rm = TRUE) / sqrt(sum(!is.na(.)))  # 计算标准误
        )
      )
    )
  return(output)
}
dsr <- result3(dsr_density) %>% as.data.frame()

dsr <- dsr[dsr$col2 %in% phyla100$col2,]

#绘制简单的柱状图
library(ggplot2)
dsr$col2 <- gsub('^p__','',dsr$col2)

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

# ggarrange(dsr_phayla_length, dsr_phayla_richness,dsr_phayla_density,
#           ncol=1, nrow=3, common.legend = TRUE, legend="bottom"
#           ) 

setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
write.csv(dsr,"Fig2CDE.csv",quote = F)
ggsave(plot = g,'3_dss_phyla_3.pdf',width = 10,height = 6.18,units = 'in')


#63个主要属中DSs的多样性####
#基于第202行得到的result_data进行下面分析
row_sums_result <- apply(result_data, 1, function(row) sum(row > 0))
# 2-02-FULL-66-22    2-12-FULL-35-15         Acidovorax      Acinetobacter       Algoriphagus 
# 21                 54                 52                 77                 39 
# Aliarcobacter       Alkanindiges      Aquabacterium       Arthrobacter           Azonexus 
# 59                 15                 65                 19                 68 
# Bacteriovorax      Brevundimonas          C7867-001        Cypionkella        Daejeonella 
# 47                 48                 20                 25                 29 
# Flavobacterium         GWA2-66-18        Gallionella        Gemmobacter       Giesbergeria 
# 74                 23                 79                 18                 43 
# Hydrogenophaga        Hylemonella           JACRCG01         Legionella      Limnohabitans 
# 59                 39                 54                 43                 36 
# Methylocystis   Methylomirabilis      Methylotenera      Mycobacterium       Nocardioides 
# 40                 34                 47                 22                 22 
# Noviherbaspirillum    Novosphingobium       Paludibacter   Pararheinheimera      Perlucidibaca 
# 32                 70                 53                 30                 35 
# Phenylobacterium       Planktophila        Polaromonas   Polynucleobacter Pseudobdellovibrio 
# 36                 12                 49                 16                 40 
# Pseudomonas        Ramlibacter         Rhodoferax          Rhodoluna       Rugosibacter 
# 62                 40                 72                  8                 38 
# Sediminibacterium      Sideroxyarcus        Sphingobium     Sphingorhabdus       Sulfuricella 
# 46                 45                 48                 40                 45 
# Sulfuricurvum       Sulfurimonas   Sulfurospirillum        Tabrizicola         Thiovibrio 
# 60                 49                 48                 33                 44 
# UBA2031            UBA2093            UBA2212            UBA3006            UBA6024 
# 29                  3                 35                 17                 45 
# UBA9959      Undibacterium            Zwartia 
# 19                 59                 19 
#绘制简单的柱状图
library(ggplot2)
row_sums_result_df <- as.data.frame(row_sums_result)
names(row_sums_result_df) <- c("Count")
dss_richness_genera <- ggplot(row_sums_result_df, aes(x = rownames(row_sums_result_df), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DSs richness") +
  scale_y_continuous(limits = c(0,80),expand = c(0,0))+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.8,0.7),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )

dss_richness_genera
setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
ggsave("3genera_dss_richness.pdf", plot = dss_richness_genera, width = 8, height = 5, units = "in")


library(cowplot)
g1 <- plot_grid(p1, p2,labels = LETTERS[1:2], nrow = 1)
g1
all <- plot_grid(g1, dss_richness,labels = c('','C'), nrow = 2)
all

ggsave(plot = all,'3CPR_pseu_men.pdf',width = 10,height = 6.18,units = 'in')
library(export)
graph2ppt(file="3.pptx", width=10, height=6.18)






#5. 查看CPR细菌门中DSs的特征####
#4.1 CPR的MAGs中DSs number density与其他门的对比
# 计算每个MAGs中的DSs密度
setwd('D:/wenjian/research/defensome/groundwater/03mag_info')
magsize <- read.csv('genome_stats.csv')
# 在R中，可以使用sub函数去掉数据框某一列字符串中的.fa后缀
# 假设数据框为df，列名为filename，可以使用以下代码：
magds <- data.frame(genome=rownames(dss),dss=rowSums(dss[,-140]))

ds_size <- merge(magds,magsize,by.x = 'genome',by.y = 'genome')

ds_density <- data.frame(genome=ds_size$genome,
                         density=ds_size$dss/ds_size$size*1000000,
                         number=ds_size$dss,
                         length_Mb=ds_size$size/1000000
                         )
# 在R中，可以使用floor函数将某一列的数字取小数点前的整数，并形成新的一列
# 假设数据框为df，列名为original_column，可以使用以下代码：
ds_density$magsize <- floor(ds_density$length_Mb)

#为ds_density 添加物种信息

ds_density_tax <- merge(ds_density,mags_tax,by.x = 'genome',by.y = 'user_genome')

#绘制箱线图比较不同门之间DSs number density的差异
# 假设你的数据框叫做ds_density_tax
# 首先加载dplyr包，它提供了mutate和if_else函数
library(dplyr)

# 使用mutate和if_else增加新的CPR列
ds_density_tax <- ds_density_tax %>%
  mutate(CPR = case_when(col2 == "p__Patescibacteria" ~ "CPR",
                         col2 == "p__Pseudomonadota" ~ "Pseudomonadota",
                         TRUE ~ "Other phyla"
                         ))

#write.csv(ds_density_tax,"ds_density_tax_r220.csv",row.names = F)


#计算显著性差异
library(rstatix)
library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means

md_p_val1 <- ds_density_tax %>% 
  wilcox_test(number ~ CPR) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="CPR",dodge=0.5)
md_p_val1$p.adj

md_p_val1$y.position <- c(13,14,15)
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


setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
write.csv(ds_density_tax,"Fig2AB.csv",quote = F)
ggsave("3CPR_DSs_number.pdf", plot = p2, width = 5.5, height = 3.8, units = "in")


#统计CPR与non-CPR的数量
sum(ds_density_tax$CPR=='CPR') #3004
sum(ds_density_tax$CPR!='CPR') #24574


#6. 统计CPR门里的各个纲的特点#######
setwd("D:/wenjian/research/defensome/groundwater/03mag_info")
ds_tax <- read.delim('ds_density_tax_r220.csv',sep = ',')
library(dplyr)
CPR <- ds_tax[ds_tax$CPR=='CPR',]
length(unique(CPR$col3))
#统计col3列每个元素的数量
element_counts <- CPR %>%
  count(col3)
#setwd('D:/wenjian/research/defensome/groundwater/01样品信息')
#write.csv(element_counts,"phyla_count.csv")

CPRc50 <- element_counts[element_counts$n>20,]
CPR50 <- CPR[CPR$col3 %in% CPRc50$col3,]

library(ggpubr)
library(rstatix)
md_p_val4 <- CPR50 %>% 
  wilcox_test(density ~ col3) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="col3",dodge=0.5)
md_p_val4$p.adj

md_p_val4$y.position <- c(13,14,15)
library(ggplot2)
p4 <- ggplot(CPR50, aes(x = col3, y = density)) +
  #geom_violin(trim = FALSE, aes(fill = CPR),show.legend = F) + # 绘制小提琴图
  geom_boxplot(aes(fill = CPR),width = 0.2,outliers = FALSE, staplewidth = 0.5) +  # 隐藏箱线图中的异常值点
  geom_jitter(width = 0.2, alpha = 0.7)+  # 显示异常值点
  # stat_compare_means(comparisons = list(c("1", "2"),c("SS", "WA"),c("SS", "WS"),
  #                                       c("SA", "WA"),c("SA", "WS"),
  #                                       c("WA", "WS")
  # ),
  # method = "wilcox.test", test = "p.adjust",label = "p.signif")+
  stat_pvalue_manual(md_p_val4,label = "p.adj.signif",hide.ns=T,
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
p4


result3 <- function(df) {
  output <- df %>%
    group_by(col3) %>%
    summarise(
      across(
        where(is.numeric),
        list(
          mean = ~mean(., na.rm = TRUE),
          sd = ~sd(., na.rm = TRUE),  # 添加标准差的计算
          se = ~sd(., na.rm = TRUE) / sqrt(sum(!is.na(.)))  # 计算标准误
        )
      )
    )
  return(output)
}
CPR50_res <- result3(CPR50)
# CPR50_res$col3 <- factor(CPR50_res$col3,levels = c(  "c__Paceibacteria" ,   "c__ABY1" , "c__JAEDAM01" ,       
#                                                      "c__Microgenomatia",   "c__Gracilibacteria",  "c__Saccharimonadia" ,
#                                                      "c__UBA1384",          "c__Andersenbacteria", "c__WWE3",            
#                                                      "c__Doudnabacteria" ))
CPR50_res$col3 <- factor(CPR50_res$col3,levels = c(  "c__Paceibacteria" , "c__Microgenomatia", "c__UBA1384", 
                                                      "c__Andersenbacteria","c__JAEDAM01" ,   "c__ABY1" ,  
                                                       "c__Doudnabacteria" , "c__Gracilibacteria", "c__WWE3",  
                                                         "c__Saccharimonadia"
                                                      ))

CPR50_density <- ggplot(CPR50_res, aes(x = col3, y = density_mean)) +
  geom_bar(stat = "identity", fill = "#43A3EF",width = 0.5) +
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
  geom_hline(yintercept = mean(CPR50_res$density_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DS density (per MAG per Mb)") +
  scale_y_continuous(limits = c(0,2.5))+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    
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
CPR50_density

setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
write.csv(CPR50_res,"Fig3E_3.csv",quote = F)
ggsave(plot = CPR50_density,'3_dss_CPR.pdf',width = 5,height = 4,units = 'in')


#CPR的各个属级别DSs分布特征热图####

CPR50$col3 <- gsub("^c__", "", CPR50$col3)
CPR20 <- tax_dss[tax_dss$col3 %in% CPRc50$col3,]
CPR20$col3 <- gsub("^c__", "", CPR20$col3)

liuxinglv <- function(df) {
  output <- df %>%
    group_by(col3) %>%
    summarise(across(where(is.numeric), ~sum(. != 0, na.rm = TRUE) / sum(!is.na(.), na.rm = TRUE)))
  return(output)
}

CPR20_liuxinglv <- liuxinglv(CPR20) %>% as.data.frame()
rownames(CPR20_liuxinglv) <- CPR20_liuxinglv$col3
CPR20_liuxinglv <- CPR20_liuxinglv[,-1]
#过滤数据
#筛选31个主要defense families
#筛选主要的DSs
# 计算每一列非0数字所占的百分比并输出到新数据框
percent_nonzero <- colSums(dss != 0) / nrow(dss) * 100
new_df <- data.frame(Column = names(percent_nonzero), Percentage = percent_nonzero)
# 按照值的大小对数据集进行降序排列
ds_persentage <- new_df[order(-new_df$Percentage), ]
ds_persentage <- ds_persentage[ds_persentage$Percentage>=1,]
result_data <- CPR20_liuxinglv[,colnames(CPR20_liuxinglv) %in% ds_persentage$Column]
 

#绘制热图
library(pheatmap)

# 设置颜色
#colors <- colorRampPalette(c('white', colorRampPalette('#F9B14F', '#8F2731')))(1000)
colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p <- pheatmap(t(result_data),
              #color = colorRampPalette(colors)(1000),
              #color = colorRampPalette(c('white','lightpink','red'))(1000),
              color = colors,
              cluster_cols = T,
              cluster_rows = F,
)
p

setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
write.csv(result_data,"Fig3D.csv",quote = F)
ggplot2::ggsave("3CPR_clade_dss.pdf", plot = p, width = 4, height = 6, units = "in")

#CPR的clade中DS的richness####
CPR_clade <- tax_dss[tax_dss$col2=="p__Patescibacteria",]

# 使用group_by和summarise函数进行分组求和
library(dplyr)
CPR_clade_res <- CPR_clade %>%
  group_by(col3) %>%
  summarise(across(.cols = where(is.numeric), sum, na.rm = TRUE))

CPR_clade_res$richness <- rowSums(CPR_clade_res[,-1]!=0)

CPR_clade_res_richness <- CPR_clade_res[,c(1,141)]

CPR_clade_num <- CPR_clade %>% count(col3)
CPR_clade_res_richness_20 <- CPR_clade_res_richness[CPR_clade_res_richness$col3 %in% CPR_clade_num[CPR_clade_num$n>=20,]$col3,]

#绘制CPR的各个属中DS的richness
library(ggplot2)
# 首先，我们需要确保col3是因子，并且按照richness的值排序

CPR_clade_res_richness_20$col3 <- factor(CPR_clade_res_richness_20$col3, 
    levels = unique(CPR_clade_res_richness_20$col3[order(CPR_clade_res_richness_20$richness,decreasing = T)]))
CPR_clade_res_richness_20$col3 <- factor(CPR_clade_res_richness_20$col3,levels = c(  "c__Paceibacteria" , "c__Microgenomatia", "c__UBA1384", 
                                                     "c__Andersenbacteria","c__JAEDAM01" ,   "c__ABY1" ,  
                                                     "c__Doudnabacteria" ,"c__Gracilibacteria",  "c__WWE3",  
                                                      "c__Saccharimonadia"
                                                     ))
CPR_clade_DSrichness_20 <- ggplot(CPR_clade_res_richness_20,aes(x=col3,y=richness))+
  geom_bar(stat = 'identity',fill = "#FEDDCE", width = 0.5)+
  labs(x=NULL,y='Defense family number')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
CPR_clade_DSrichness_20
setwd("D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense")
write.csv(CPR_clade_res_richness_20,"Fig3E_1.csv",quote = F)
ggsave("3_CPR_clade_DSrichness.pdf", plot = CPR_clade_DSrichness_20, width = 3.6, height = 3.8, units = "in")

#CPR各属DS的richness+ DS density + DS number + DS richness per MAG#########
CPR_clade$DSrichness <- rowSums(CPR_clade[,c(2:140)])

CPR20_clade_dsrichnessdensity <- CPR_clade[,c(1,143,148)]
CPR20_clade_dsrichnessdensity <- 
  CPR20_clade_dsrichnessdensity[CPR20_clade_dsrichnessdensity$col3 %in% CPR_clade_num[CPR_clade_num$n>=20,]$col3,]

result3 <- function(df) {
  output <- df %>%
    group_by(col3) %>%
    summarise(
      across(
        where(is.numeric),
        list(
          mean = ~mean(., na.rm = TRUE),
          sd = ~sd(., na.rm = TRUE),  # 添加标准差的计算
          se = ~sd(., na.rm = TRUE) / sqrt(sum(!is.na(.)))  # 计算标准误
        )
      )
    )
  return(output)
}
CPR20_res <- result3(CPR20_clade_dsrichnessdensity)

CPR20_res$col3 <- factor(CPR20_res$col3,levels = c(   "c__Paceibacteria" , "c__Microgenomatia", "c__UBA1384", 
                                                      "c__Andersenbacteria","c__JAEDAM01" ,   "c__ABY1" ,  
                                                      "c__Doudnabacteria" , "c__Gracilibacteria",  "c__WWE3",  
                                                      "c__Saccharimonadia"
                                                      ))


CPR20_density <- ggplot(CPR20_res, aes(x = col3, y = DSrichness_mean)) +
  geom_bar(stat = "identity", fill = "#d4e09b",width = 0.5) +
  geom_errorbar(
    aes(ymin = DSrichness_mean - 1.96*DSrichness_se, ymax = DSrichness_mean + 1.96*DSrichness_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = DSrichness_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 3  # 点大小
  )+
  geom_hline(yintercept = mean(CPR20_res$DSrichness_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DS richness (per MAG)") +
  scale_y_continuous(limits = c(0,2.8),labels = scales::number_format(accuracy = 0.1))+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
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
CPR20_density

setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
write.csv(CPR20_res,"Fig3E_2.csv",quote = F)
ggsave(plot = CPR50_density,'3_CPR_DSrichnessperMAG.pdf',width = 5,height = 4,units = 'in')


#合并上面几个图，共用横坐标
#CPR50_res
#CPR_clade_res_richness_20
#CPR20_res
CPR_ds_all <- merge(CPR50_res,CPR_clade_res_richness_20,by = 'col3')
CPR_ds_all <- merge(CPR_ds_all,CPR20_res)

library(Cairo)
library(cowplot)
cairo_pdf("3_CPR_clade_all2.pdf",width = 2.9,height = 6.78)
g <- plot_grid(CPR50_density,CPR20_density,CPR_clade_DSrichness_20,
               #labels = LETTERS[1:3],
               nrow = 3,
               rel_heights = c(1, 1, 1.2)
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g

setwd("D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense")
dev.off()




#7. DHR分布图#########
#每个门的DSs总数除以该门MAGs数量
setwd("D:/wenjian/research/defensome/groundwater/03mag_info")
ds_tax <- read.delim('ds_density_tax_r220.csv',sep = ',',check.names = FALSE)

DHR <- ds_tax %>% group_by(col2) %>% 
  summarise(mean_number = mean(number))
DHR$col2 <- gsub('^p__','',DHR$col2)

DHR$phy_type <- 'Others'
DHR[DHR$col2=='Patescibacteria',]$phy_type <- 'CPR'
DHR[DHR$col2=='Pseudomonadota',]$phy_type <- 'Pseudomonadota'
DHR <- DHR[DHR$col2!='' & DHR$mean_number>0,]

library(ggpubr)
dotp <- ggdotchart(DHR,x = 'col2', y = 'mean_number',
                 color = 'phy_type',                                # Color by groups
                  palette = c("#00AFBB","grey", "#E7B800", "#FC4E07"), # Custom color palette
                 sorting = "descending",                       # Sort value in descending order
                 add = "segments",                             # Add segments from y = 0 to dots
                 add.params = list(color = 'phy_type', size = 0.5), # Change segment color and size
                 #  group = "cyl",                                # Order by groups
                 dot.size = 1,                                 # Large dot size
                 #  label = round(dfm$mpg_z,1),                        # Add mpg values as dot labels
                 font.label = list(color = "white", size = 9,
                                   vjust = 0.5),               # Adjust label parameters
                 ggtheme = theme_pubr()                        # ggplot2 theme
) +
  labs(x=NULL,y='DSs per MAG')+
  theme(
    panel.background = element_blank(),  # 背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    
    # 设置边框粗细一致
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 边界线颜色和粗细
    panel.margin = unit(0.5, "lines"),  # 设置边框的粗细，与panel.border的size参数一致
    
    # 设置坐标轴标题和标签
    axis.title = element_text(face = "bold", size = 10, color = "grey20"),  # 坐标轴标题样式
    axis.text = element_text(size = 8, color = "grey20"),  # 坐标轴标签样式
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 坐标轴刻度线样式
    axis.line = element_blank()
  )+
  theme(
    legend.title = element_blank(),       # 删除图例标题
    legend.position = c(0.9, 0.9),  # 设置图例在图的右上角，使用1表示x和y坐标的最大值
    legend.justification = c(1, 1),  # 将图例对齐到右上角
    legend.direction = "horizontal",  # 设置图例水平排列
    #legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 图例标题样式
    legend.text = element_text(size = 10, color = "grey20"),  # 图例文本样式
    legend.background = element_blank(),  # 图例背景透明
    legend.box.background = element_blank(),
    #legend.box.background = element_rect(color = "grey80", fill = NA, size = 0.5),  # 图例边框样式，fill设置为NA以确保背景透明
    legend.key = element_blank()  # 去掉图例键的背景
  )
dotp
write.csv(DHR,"figs9a.csv")
#每个门的DSs数量除以MAGs大小之后再除以MAGs数量（density的平均值）
DHR2 <- ds_tax %>% group_by(col2) %>% 
  summarise(mean_density = mean(density))
DHR2$col2 <- gsub('^p__','',DHR2$col2)

DHR2$phy_type <- 'Others'
DHR2[DHR2$col2=='Patescibacteria',]$phy_type <- 'CPR'
DHR2[DHR2$col2=='Pseudomonadota',]$phy_type <- 'Pseudomonadota'
DHR2 <- DHR2[DHR2$col2!='' & DHR2$mean_density>0,]

library(ggpubr)
dotp2 <- ggdotchart(DHR2,x = 'col2', y = 'mean_density',
                   color = 'phy_type',                                # Color by groups
                   palette = c("#00AFBB","grey", "#E7B800", "#FC4E07"), # Custom color palette
                   sorting = "descending",                       # Sort value in descending order
                   add = "segments",                             # Add segments from y = 0 to dots
                   add.params = list(color = 'phy_type', size = 0.5), # Change segment color and size
                   #  group = "cyl",                                # Order by groups
                   dot.size = 1,                                 # Large dot size
                   #  label = round(dfm$mpg_z,1),                        # Add mpg values as dot labels
                   font.label = list(color = "white", size = 9,
                                     vjust = 0.5),               # Adjust label parameters
                   ggtheme = theme_pubr()                        # ggplot2 theme
) +
  labs(x=NULL,y='DSs density')+
  theme(
    panel.background = element_blank(),  # 背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    
    # 设置边框粗细一致
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 边界线颜色和粗细
    panel.margin = unit(0.5, "lines"),  # 设置边框的粗细，与panel.border的size参数一致
    
    # 设置坐标轴标题和标签
    axis.title = element_text(face = "bold", size = 10, color = "grey20"),  # 坐标轴标题样式
    axis.text = element_text(size = 8, color = "grey20"),  # 坐标轴标签样式
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 坐标轴刻度线样式
    axis.line = element_blank()
  )+
  theme(
    legend.title = element_blank(),       # 删除图例标题
    legend.position = c(0.9, 0.9),  # 设置图例在图的右上角，使用1表示x和y坐标的最大值
    legend.justification = c(1, 1),  # 将图例对齐到右上角
    legend.direction = "horizontal",  # 设置图例水平排列
    #legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 图例标题样式
    legend.text = element_text(size = 10, color = "grey20"),  # 图例文本样式
    legend.background = element_blank(),  # 图例背景透明
    legend.box.background = element_blank(),
    #legend.box.background = element_rect(color = "grey80", fill = NA, size = 0.5),  # 图例边框样式，fill设置为NA以确保背景透明
    legend.key = element_blank()  # 去掉图例键的背景
  )
dotp2
write.csv(DHR2,"figs9b.csv")
library(cowplot)
g <- plot_grid(dotp, dotp2,
               labels = LETTERS[1:2], nrow = 2
               #rel_heights = c(1, 1, 1.5),
              # label_x = 0.09,            # 标签x坐标（1为最右侧）
             #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g

setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
ggsave(plot = g,'3_DHR.pdf',width = 12,height = 6.18,units = 'in')


#8. 查看NCBI分类物种#####
setwd("D:/wenjian/research/defensome/groundwater/03mag_info")
ncbi_tax <- read.delim2("gtdbtk.bac120.summary_ncbi2.txt",)

#处理物种信息
ncbitaxonoxy <- function(tax){
  tax <- tax[,1:3]
  max_cols <- max(sapply(strsplit(tax$`Majority.vote.NCBI.classification`, ";"), length))# 确定最大的拆分列数
  col_names <- paste0("col", 1:max_cols)# 生成列名
  library(tidyr)
  taxonomy <- tax %>%
    separate(`Majority.vote.NCBI.classification`, into = col_names, sep = ";", fill = "right")# 拆分字符串到列
}

ncbi_tax <- ncbitaxonoxy(ncbi_tax)
bins <- rownames(dss)
bins <- sub("*.fasta","",bins)
ncbi_2 <- ncbi_tax[ncbi_tax$Genome.ID %in% bins,]


#9. 计算每个门的MAG size，DS number，DS density，defense family（diversity）。








