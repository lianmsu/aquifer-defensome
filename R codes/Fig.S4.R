
#1. 读取数据及其基本处理#############################################

#DSs读取
setwd("D:/wenjian/research/defensome/groundwater/archaea/02defense")
ar_defense <- read.csv("archaea_ds_counts_matrix.csv",check.names = F)
ar_defense$genome <- gsub('*.fna','',ar_defense$genome)

#50_10_MAG,checkm1
setwd("D:/wenjian/research/defensome/groundwater/archaea/01MAG_info")
MAG_5010 <- readxl::read_xlsx("2269_metainfo.xlsx")

#70_10_MAG,checkm1
setwd("D:/wenjian/research/defensome/groundwater/archaea/01MAG_info")
MAG_7010 <- readxl::read_xlsx("zhong_70_10.xlsx")

#MAG_size_info,checkm2
setwd("D:/wenjian/research/defensome/groundwater/archaea/01MAG_info")
MAG_checkm2 <- read.csv("quality_report_checkm2.tsv",sep = '\t',check.names = F)

#合并不含DS的MAG
not_in <- MAG_5010$binname[!MAG_5010$binname %in% ar_defense$genome]
# 使用不在 ar_defense 行名中的 genome 值创建新的数据框
new_rows <- data.frame(
  genome = not_in,  # 重复 genome 值
  stringsAsFactors = FALSE)
# 为 new_rows 添加与 ar_defense 相同的列，并将这些列的值填充为 0
col_names <- setdiff(names(ar_defense), "genome")  # 获取 dss 中除了 genome 之外的列名
new_rows[col_names] <- lapply(col_names, function(col_name) rep(0, nrow(ar_defense)))
rownames(new_rows) <- new_rows$genome
#合并new_rows与dss
ar_defense <- rbind(ar_defense,new_rows)

#提取70_10的MAGs
ar_defense_7010 <- ar_defense[ar_defense$genome %in% MAG_7010$binname,]
#setwd("D:/wenjian/research/defensome/groundwater/archaea/02defense")
#write.csv(file = 'ar_defense_7010.csv', ar_defense_7010,row.names = F)

#2. 古菌的DSs分布###############################################################

#2.1数据预处理
rownames(ar_defense_7010) <- ar_defense_7010$genome
ar_defense_7010 <- ar_defense_7010[,-1]

# 使用rowSums()对每行求和
row_sums <- rowSums(ar_defense_7010)
# 将求和结果转换为向量
sum_vector <- as.vector(row_sums)
row_names <- rownames(ar_defense_7010)
# 创建一个新的数据框，包含求和结果和行名
ar_7010_num_dss <- data.frame(Sum = sum_vector, row.names = row_names)

# 首先，确保"Sum"列的值是数值类型
ar_7010_num_dss$Sum <- as.numeric(ar_7010_num_dss$Sum)
# 使用table()函数统计Sum列每个数值的出现次数
sum_counts <- table(ar_7010_num_dss$Sum)
# 将统计结果转换为数据框
sum_counts_df <- as.data.frame(sum_counts)
# 重命名列
colnames(sum_counts_df) <- c("Value", "Frequency")
# 确保Value列是因子类型
sum_counts_df$Value <- as.factor(sum_counts_df$Value)

#2.2绘制柱状图
library(ggplot2)
# 使用geom_bar()绘制柱状图，指定stat为"identity"
p1 <- ggplot(sum_counts_df,aes(x=Value,y=Frequency))+
  geom_bar(stat = "identity",position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,fill="#E9AF89" # "#DFA4BA" "#E9AF89" "#81D0D6" "#B3DFC4"       
  )+
  geom_text(aes(label = Frequency), angle = 90, hjust = -0.1,
            position = position_dodge(width = 0.9), 
            # vjust = -0.50, # 调整文本位置，避免与柱子重叠
            color = "black", # 文本颜色，根据柱子颜色调整
            size = 3.5)+ # 文本大小
  scale_y_continuous(limits = c(0,600),expand = c(0,0))+
  labs( x="Number of DSs",y="Number of MAGs") +
  theme(axis.text.x=element_text(angle = 0,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(color="black"),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.line.x=element_line(color="black"),
        axis.line.y.left = element_line(color="grey30"),
        axis.line.y.right = element_line(color="grey30"),
        axis.line.x.bottom = element_line(color="grey30"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))
p1  
write.csv(sum_counts_df,"FigS4B.csv")
# 导出PDF
ggsave("1_p1_DSs_distribution.pdf", plot = p1, width = 5.5, height = 3.5, units = "in")

#3. Defense family组成##########################################################

#3.1数据预处理
# 计算每一列非0数字所占的百分比并输出到新数据框
percent_nonzero <- colSums(ar_defense_7010 != 0) / nrow(ar_defense_7010) * 100
new_df <- data.frame(Column = names(percent_nonzero), Percentage = percent_nonzero)
# 按照值的大小对数据集进行降序排列
ds_persentage <- new_df[order(-new_df$Percentage), ]

#3.2绘制柱状图
library(ggplot2)
p2 <- ggplot(new_df[new_df$Percentage>=1,]) +
  geom_bar(aes(x = reorder(Column,-Percentage), y = Percentage),
           stat = "identity", position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,
           fill="#E9AF89" # '#81D0D6'  "#E9AF89",sa"#9CB9D7",ss  "#DFA4BA",ws"#E9AF89",sa
  )+
  scale_y_continuous(limits = c(0,50),expand = c(0,0))+
  labs( x=NULL,y="Proportion present in 1453 MAGs (%)") +
  theme(axis.text.x=element_text(angle = 90,vjust=0.5,hjust=0.9,color="black"),
        axis.text.y=element_text(color="black"),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        text = element_text(family = 'sans'),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.line.x=element_line(color="black"),
        axis.line.y.left = element_line(color="grey30"),
        axis.line.y.right = element_line(color="grey30"),
        axis.line.x.bottom = element_line(color="grey30"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))
p2
write.csv(new_df,"FigS4A.csv")
ggsave("1_p2_DSs_percentage.pdf", plot = p2, width = 5, height = 3.5, units = "in")

#4. treemap展示DSs在MAGs中的平均数量############################################

#4.1数据预处理

#读取数据库中effector类型
setwd('D:/wenjian/research/defensome/groundwater/01defensefinder_info')
allds <- read.delim2('df_List_Systems_lpw.csv',sep = ',')
allds$Effector[allds$Effector == ""] <- "Unknown"
allds$title <- gsub("CasFinder", "Cas", allds$title)

#计算每个DSs在所有MAGs中的平均number
library(dplyr)
# 假设ar_defense_7010是你要处理的数据框,计算每列均值
means_df <- as.data.frame(colMeans(ar_defense_7010))
means_df$dsstype <- rownames(means_df)
#合并，注意：Cas与Casfinder，两个不对应
num_average <- merge(means_df,allds,by.x = 'dsstype',by.y = 'title')
colnames(num_average)[2] <- 'num_means'

#4.2 绘图treemap
library(treemap)
library(RColorBrewer)
library(treemapify)
# 添加多个分组信息
p3 <- ggplot(num_average, aes(area = num_means, fill = Effector,
                              label = dsstype,
                              subgroup = Effector)) +
  geom_treemap(colour='black') +
  geom_treemap_subgroup_border() +
  geom_treemap_subgroup_text(place = "centre", alpha = 0.5, 
                             colour = "white", fontface = "italic",family='Times') +
  geom_treemap_text(colour = "black", place = "topleft", reflow = T,family='Times') +
  scale_fill_manual(values = c("#FF6384", "#36A2EB", "#FFCE56", "#4BC0C0", "#C0C0C0"))+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_blank()
        )+
  guides(fill = guide_legend(nrow = 2))
p3
write.csv(num_average,"FigS4C.csv")
setwd("D:/wenjian/research/defensome/groundwater/archaea/02defense")
ggsave(plot = p3,"1_p3_mean_dsnum_treemap.pdf",width = 4, height = 4, units = "in")

#5. DPANN与其他门相比的特点(MAG size,DS richness,DS number,DS density)##########

#5.1 数据预处理

md_archaea <- merge(MAG_7010,MAG_checkm2,by.x = 'binname',by.y = 'Name')
ar_de_sum <- data.frame(ds_sum = rowSums(ar_defense_7010[,-1]),
                        genome = rownames(ar_defense_7010),
                        richness = apply(ar_defense_7010[,-1], 1, function(row) sum(row > 0))
                        )
md_archaea <- merge(md_archaea,ar_de_sum,by.x = 'binname',by.y = 'genome')
md_archaea$ds_density <- md_archaea$ds_sum/md_archaea$size*1000000

#5.2 统计比较
library(rstatix)

md_p_val1 <- md_archaea %>% 
  wilcox_test(ds_sum ~ superphylum) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="superphylum",dodge=0.5)
md_p_val1$p.adj

md_p_val2 <- md_archaea %>% 
  wilcox_test(ds_density ~ superphylum) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="superphylum",dodge=0.5)
md_p_val2$p.adj

md_p_val3 <- md_archaea %>% 
  wilcox_test(richness ~ superphylum) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="superphylum",dodge=0.5)
md_p_val3$p.adj

#5.3 绘制箱线+小提琴图
library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means

md_p_val1$y.position <- c(13,14,15)
p4 <- ggplot(md_archaea, aes(x = superphylum, y = ds_sum)) +
  geom_violin(trim = FALSE, aes(fill = superphylum),show.legend = F) + # 绘制小提琴图
  geom_boxplot(aes(fill = superphylum),width = 0.2,outliers = FALSE, staplewidth = 0.5) +  # 隐藏箱线图中的异常值点
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
  scale_fill_manual(values = c("DPANN"="#FF7F0E","Euryarchaeota" = "#ADD9EE","TACK"="#7892C7"
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
    axis.title = element_text( size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p4
write.csv(md_archaea,"FigS4EFG.csv")

p5 <- ggplot(md_archaea, aes(x = superphylum, y = ds_density)) +
  geom_violin(trim = FALSE, aes(fill = superphylum),show.legend = F) + # 绘制小提琴图
  geom_boxplot(aes(fill = superphylum),width = 0.2,outliers = FALSE, staplewidth = 0.5) +  # 隐藏箱线图中的异常值点
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
  scale_fill_manual(values = c("DPANN"="#FF7F0E","Euryarchaeota" = "#ADD9EE","TACK"="#7892C7"
  )) +
  #ylim(0,15)+
  #scale_y_continuous(limits = c(0,10))+
  coord_cartesian(ylim = c(0,15))+
  labs(x=NULL,y='DSs density (per MAG per Mb)')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p5

p6 <- ggplot(md_archaea, aes(x = superphylum, y = richness)) +
  geom_violin(trim = FALSE, aes(fill = superphylum),show.legend = F) + # 绘制小提琴图
  geom_boxplot(aes(fill = superphylum),width = 0.2,outliers = FALSE, staplewidth = 0.5) +  # 隐藏箱线图中的异常值点
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
  scale_fill_manual(values = c("DPANN"="#FF7F0E","Euryarchaeota" = "#ADD9EE","TACK"="#7892C7"
  )) +
  #ylim(0,15)+
  #scale_y_continuous(limits = c(0,10))+
  coord_cartesian(ylim = c(0,15))+
  labs(x=NULL,y='Defense families (per MAG)')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p6

library(cowplot)
g1 <- plot_grid(p4,p5,p6,
               labels = LETTERS[1:3], nrow = 1
               # rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               # label_y = 0.9             # 标签y坐标（1为最上方）
)
g1

setwd("D:/wenjian/research/defensome/groundwater/archaea/02defense")
ggsave("1_g1_superphylum_compare.pdf", plot = g1, width = 9, height = 3, units = "in")

#6. 去除不含DSs的MAG后，再看一下。###############################################

#6.1 统计检验
md_p_val4 <- md_archaea[md_archaea$ds_sum>0,] %>% 
  wilcox_test(ds_sum ~ superphylum) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="superphylum",dodge=0.5)
md_p_val4$p.adj

md_p_val5 <- md_archaea[md_archaea$ds_sum>0,] %>% 
  wilcox_test(ds_density ~ superphylum) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="superphylum",dodge=0.5)
md_p_val5$p.adj

md_p_val6 <- md_archaea[md_archaea$ds_sum>0,] %>% 
  wilcox_test(richness ~ superphylum) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="superphylum",dodge=0.5)
md_p_val6$p.adj

#6.2 绘制箱线+小提琴图
library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means

md_p_val4$y.position <- c(13,14,15)
p7 <- ggplot(md_archaea[md_archaea$ds_sum>0,], aes(x = superphylum, y = ds_sum)) +
  geom_violin(trim = FALSE, aes(fill = superphylum),show.legend = F) + # 绘制小提琴图
  geom_boxplot(aes(fill = superphylum),width = 0.2,outliers = FALSE, staplewidth = 0.5) +  # 隐藏箱线图中的异常值点
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
  scale_fill_manual(values = c("DPANN"="#FF7F0E","Euryarchaeota" = "#ADD9EE","TACK"="#7892C7"
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
    axis.title = element_text( size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p7

p8 <- ggplot(md_archaea[md_archaea$ds_sum>0,], aes(x = superphylum, y = ds_density)) +
  geom_violin(trim = FALSE, aes(fill = superphylum),show.legend = F) + # 绘制小提琴图
  geom_boxplot(aes(fill = superphylum),width = 0.2,outliers = FALSE, staplewidth = 0.5) +  # 隐藏箱线图中的异常值点
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
  scale_fill_manual(values = c("DPANN"="#FF7F0E","Euryarchaeota" = "#ADD9EE","TACK"="#7892C7"
  )) +
  #ylim(0,15)+
  #scale_y_continuous(limits = c(0,10))+
  coord_cartesian(ylim = c(0,15))+
  labs(x=NULL,y='DSs density (per MAG per Mb)')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p8

p9 <- ggplot(md_archaea[md_archaea$ds_sum>0,], aes(x = superphylum, y = richness)) +
  geom_violin(trim = FALSE, aes(fill = superphylum),show.legend = F) + # 绘制小提琴图
  geom_boxplot(aes(fill = superphylum),width = 0.2,outliers = FALSE, staplewidth = 0.5) +  # 隐藏箱线图中的异常值点
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
  scale_fill_manual(values = c("DPANN"="#FF7F0E","Euryarchaeota" = "#ADD9EE","TACK"="#7892C7"
  )) +
  #ylim(0,15)+
  #scale_y_continuous(limits = c(0,10))+
  coord_cartesian(ylim = c(0,15))+
  labs(x=NULL,y='Defense families (per MAG)')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p9

library(cowplot)
g1 <- plot_grid(p4,p5,p6,
                labels = LETTERS[1:3], nrow = 1
                # rel_heights = c(1, 1, 1.5),
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                # label_y = 0.9             # 标签y坐标（1为最上方）
)
g1

setwd("D:/wenjian/research/defensome/groundwater/archaea/02defense")
ggsave("1_g1_superphylum_compare.pdf", plot = g1, width = 9, height = 3, units = "in")

#7. 古菌中各个门中的DS流行率特点#################################################

#7.1 数据预处理
md_archaea <- merge(MAG_7010,MAG_checkm2,by.x = 'binname',by.y = 'Name')
ar_de_sum <- data.frame(ds_sum = rowSums(ar_defense_7010[,-1]),
                        genome = rownames(ar_defense_7010),
                        richness = apply(ar_defense_7010[,-1], 1, function(row) sum(row > 0))
)
md_archaea <- merge(md_archaea,ar_de_sum,by.x = 'binname',by.y = 'genome')
md_archaea$ds_density <- md_archaea$ds_sum/md_archaea$size*1000000

#处理物种信息
taxonoxy <- function(tax){
  tax <- tax[,c(1,4)]
  max_cols <- max(sapply(strsplit(tax$`GTDB classification v220`, ";"), length))# 确定最大的拆分列数
  col_names <- paste0("col", 1:max_cols)# 生成列名
  library(tidyr)
  taxonomy <- tax %>%
    separate(`GTDB classification v220`, into = col_names, sep = ";", fill = "right")# 拆分字符串到列
}
ar_tax <- taxonoxy(md_archaea)
md_archaea <- merge(md_archaea,ar_tax,by.x = 'binname',by.y = 'binname')

setwd("D:/wenjian/research/defensome/groundwater/archaea/01MAG_info")
write.csv( md_archaea,"md_archaea.csv", quote = F)

ar_defense_7010$genome <- rownames(ar_defense_7010)
tax_dss <- merge(ar_defense_7010,ar_tax,by.x = 'genome',by.y = 'binname')

#7.2 统计处理

#统计col2列每个元素的数量(统计每个门的MAG数量)
library(dplyr)
phylum_counts <- md_archaea %>%
  count(col2)

#分组统计
liuxinglv <- function(df) {
  output <- df %>%
    group_by(col2) %>%
    summarise(across(where(is.numeric), ~sum(. != 0, na.rm = TRUE) / sum(!is.na(.), na.rm = TRUE)))
  return(output)
}
phyla_liuxing <- liuxinglv(tax_dss) %>% as.data.frame()
rownames(phyla_liuxing) <- phyla_liuxing$col2 %>% gsub('p__*','',.)
phyla_liuxing <- phyla_liuxing[,-1]

#筛选主要defense families,在所有MAG中出现率大于0.5%
#筛选主要的DSs
# 计算每一列非0数字所占的百分比并输出到新数据框
percent_nonzero <- colSums(ar_defense_7010[,-1] != 0) / nrow(ar_defense_7010) * 100
new_df <- data.frame(Column = names(percent_nonzero), Percentage = percent_nonzero)
# 按照值的大小对数据集进行降序排列
ds_persentage <- new_df[order(-new_df$Percentage), ]
ds_persentage <- ds_persentage[ds_persentage$Percentage>=0.5,]

phyla_liuxing <- phyla_liuxing[,colnames(phyla_liuxing)%in%ds_persentage$Column]

#7.3 绘制热图
#绘制热图
library(pheatmap)

#设置门的超门
# super_phyla <- select(md_archaea,col2,superphylum) %>% unique() %>% as.data.frame()
# super_phyla$col2 <- gsub('p__*','', super_phyla$col2) 
# rownames(super_phyla) <- super_phyla$col2 
# super_phyla <- super_phyla[,-1]

# 设置颜色
#colors <- colorRampPalette(c('white', colorRampPalette('#F9B14F', '#8F2731')))(1000)
colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p10 <- pheatmap(t(phyla_liuxing),
              #color = colorRampPalette(colors)(1000),
              #color = colorRampPalette(c('white','lightpink','red'))(1000),
              color = colors,
              cluster_cols = T,
              cluster_rows = F,
              #annotation_col = super_phyla[,-1],
              fontsize_row = 10
              )
p10
setwd("D:/wenjian/research/defensome/groundwater/archaea/02defense")
write.csv(phyla_liuxing,"FigS4H.csv")
library(ggplot2)
ggsave("1_p10_DPANN_dss_liuxinglv.pdf", plot = p10, width = 4, height = 6, units = "in")

#只绘制DPANN的门中的DS流行率！

tax_dss_dpann <- tax_dss[tax_dss$genome %in% md_archaea[md_archaea$superphylum=='DPANN',]$binname,]
dpann_phyla <-tax_dss_dpann %>% count(col2)
phyla_liuxing2 <- liuxinglv(tax_dss_dpann) %>% as.data.frame()
rownames(phyla_liuxing2) <- phyla_liuxing2$col2 %>% gsub('p__*','',.)
phyla_liuxing2 <- phyla_liuxing2[,-1]

#筛选主要defense families,在所有MAG中出现率大于0.5%
#筛选主要的DSs
# 计算每一列非0数字所占的百分比并输出到新数据框
percent_nonzero <- colSums(ar_defense_7010[,-1] != 0) / nrow(ar_defense_7010) * 100
new_df <- data.frame(Column = names(percent_nonzero), Percentage = percent_nonzero)
# 按照值的大小对数据集进行降序排列
ds_persentage <- new_df[order(-new_df$Percentage), ]
ds_persentage <- ds_persentage[ds_persentage$Percentage>=0.5,]

phyla_liuxing2 <- phyla_liuxing2[,colnames(phyla_liuxing2)%in%ds_persentage$Column]

#7.3 绘制热图
#绘制热图
library(pheatmap)

#设置门的超门
# super_phyla <- select(md_archaea,col2,superphylum) %>% unique() %>% as.data.frame()
# super_phyla$col2 <- gsub('p__*','', super_phyla$col2) 
# rownames(super_phyla) <- super_phyla$col2 
# super_phyla <- super_phyla[,-1]

# 设置颜色
#colors <- colorRampPalette(c('white', colorRampPalette('#F9B14F', '#8F2731')))(1000)
colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p11 <- pheatmap(t(phyla_liuxing2),
                #color = colorRampPalette(colors)(1000),
                #color = colorRampPalette(c('white','lightpink','red'))(1000),
                color = colors,
                cluster_cols = T,
                cluster_rows = F,
                #annotation_col = super_phyla[,-1],
                fontsize_row = 10
)

setwd("D:/wenjian/research/defensome/groundwater/archaea/02defense")
write.csv(phyla_liuxing2,"Fig3F.csv",quote = F)
library(ggplot2)
ggsave("1_p11_DPANN_dss_liuxinglv.pdf", plot = p11, width = 4, height = 6, units = "in")

#7.4绘制DPANN中各个门的DS richness####
# 使用group_by和summarise函数进行分组求和
library(dplyr)
dpann_phyla_res <- tax_dss_dpann %>%
  group_by(col2) %>%
  summarise(across(.cols = where(is.numeric), sum, na.rm = TRUE))

dpann_phyla_res$richness <- rowSums(dpann_phyla_res[,-1]!=0)

dpann_phyla_res_richness <- dpann_phyla_res[,c(1,83)]

dpann_phyla_res_richness_10 <- dpann_phyla_res_richness[dpann_phyla_res_richness$col2 %in% dpann_phyla[dpann_phyla$n>=10,]$col2,]

#绘制DPANN的各个门中DS的richness
library(ggplot2)
# 首先，我们需要确保col3是因子，并且按照richness的值排序

dpann_phyla_res_richness_10$col2 <- factor(dpann_phyla_res_richness_10$col2, 
                                         levels = unique(dpann_phyla_res_richness_10$col2[order(dpann_phyla_res_richness_10$richness,decreasing = T)]))
dpann_phyla_res_richness_10$col2 <- factor(dpann_phyla_res_richness_10$col2,levels = c("p__Nanoarchaeota","p__Aenigmatarchaeota","p__B1Sed10-29",
                                                                                       "p__Micrarchaeota", "p__Iainarchaeota", "p__Undinarchaeota",   
                                                                                       "p__Altiarchaeota", "p__EX4484-52","p__SpSt-1190"
                                                                                       ))

dpann_phyla_DSrichness_10 <- ggplot(dpann_phyla_res_richness_10,aes(x=col2,y=richness))+
  geom_bar(stat = 'identity',fill = "#CDF0F2",width = 0.5)+
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
dpann_phyla_DSrichness_10
setwd("D:/wenjian/research/defensome/groundwater/archaea/02defense")
write.csv(dpann_phyla_res_richness_10,"Fig3G_1.csv",quote = F)
ggsave("1_dpann_phyla_DSrichness_10.pdf", plot = dpann_phyla_DSrichness_10, width = 3.5, height = 3.8, units = "in")

#7.5 绘制古菌DPANN内各个门中，DS richness per MAG#####
tax_dss_dpann <- tax_dss[tax_dss$genome %in% md_archaea[md_archaea$superphylum=='DPANN',]$binname,]
dpann_phyla <-tax_dss_dpann %>% count(col2)

tax_dss_dpann$DSrichness <- rowSums(tax_dss_dpann[,c(2:82)])
dpann_phyla_richness <- tax_dss_dpann[,c(1,84,90)]

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
dpann_res <- result3(dpann_phyla_richness)

dpann_res$col2 <- factor(dpann_res$col2, levels = c("p__Nanoarchaeota","p__Aenigmatarchaeota","p__B1Sed10-29",
                                                    "p__Micrarchaeota", "p__Iainarchaeota", "p__Undinarchaeota",   
                                               "p__Altiarchaeota", "p__EX4484-52","p__SpSt-1190"
                                               ))

dpann_density <- ggplot(dpann_res, aes(x = col2, y = DSrichness_mean)) +
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
  geom_hline(yintercept = mean(dpann_res$DSrichness_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DS richness(per MAG)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
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
dpann_density

setwd('D:/wenjian/research/defensome/groundwater/archaea/02defense')
write.csv(dpann_res,"Fig3G_2.csv",quote = F)
ggsave(plot = dpann_density,'1_dpann_ds_richnessperMAG_phyla_10.pdf',width = 4,height = 4,units = 'in')


#8. 古菌DPANN的各个门中的DSs density（per MAG per Mb）##########################

#8.1处理数据，来自7.1
dpann <- md_archaea[md_archaea$superphylum=="DPANN",]
#统计col2列每个元素的数量
element_counts <- dpann %>%
  count(col2.x)

dpann50 <- element_counts[element_counts$n>50,]
dpann50 <- dpann[dpann$col2 %in% dpann50$col2,]

dpann10 <- element_counts[element_counts$n>=10,]
dpann10 <- dpann[dpann$col2.x %in% dpann10$col2.x,]

#8.2统计分析
library(rstatix)
md_p_val4 <- dpann10 %>% 
  wilcox_test(ds_density ~ col2) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="col2",dodge=0.5)
md_p_val4$p.adj

#8.3绘制箱线图
library(ggpubr)
p11 <- ggplot(dpann10, aes(x = col2, y = ds_density)) +
  #geom_violin(trim = FALSE, aes(fill = CPR),show.legend = F) + # 绘制小提琴图
  geom_boxplot(aes(fill = col2),width = 0.2,outliers = FALSE, staplewidth = 0.5) +  # 隐藏箱线图中的异常值点
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
p11

#8.4柱状图，dpann各个门中的DSs density
result3 <- function(df) {
  output <- df %>%
    group_by(col2.x) %>%
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
dpann_density_res <- result3(dpann10)

dpann_density_res$col2 <- factor(dpann_density_res$col2.x, levels = c("p__Nanoarchaeota","p__Aenigmatarchaeota","p__B1Sed10-29",
                                                                    "p__Micrarchaeota", "p__Iainarchaeota", "p__Undinarchaeota",   
                                                                    "p__Altiarchaeota", "p__EX4484-52","p__SpSt-1190"
                                                                    ))
dpann_ds_density <- ggplot(dpann_density_res, aes(x = col2, y = ds_density_mean)) +
  geom_bar(stat = "identity", fill = "#43A3EF",width = 0.5) +
  geom_errorbar(
    aes(ymin = ds_density_mean - 1.96*ds_density_se, ymax = ds_density_mean + 1.96*ds_density_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = ds_density_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 3  # 点大小
  )+
  geom_hline(yintercept = mean(dpann_density_res$ds_density_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DS density (per MAG per Mb)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
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
dpann_ds_density

setwd('D:/wenjian/research/defensome/groundwater/archaea/02defense')
write.csv(dpann_density_res,"Fig3G_3.csv",quote = F)
ggsave(plot = dpann_density,'1_dpann_ds_density_phyla_10.pdf',width = 4,height = 4,units = 'in')


#########################
#合并绘制三个DPANN各门的DS

library(Cairo)
library(cowplot)
cairo_pdf("1_DPANN_phyla_all2.pdf",width = 2.9,height = 6.78)
g <- plot_grid(dpann_ds_density,dpann_density,dpann_phyla_DSrichness_10,
               #labels = LETTERS[1:3],
               nrow = 3,
               rel_heights = c(1, 1, 1.2)
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g

setwd("D:/wenjian/research/defensome/groundwater/archaea/02defense")
dev.off()


#9. 读取DGs数据####

#DGs读取
setwd("D:/wenjian/research/defensome/groundwater/archaea/02defense")
ar_dg <- read.csv("archaea_dg_counts_matrix.csv",check.names = F)
ar_dg$genome <- gsub("*.fna","",ar_dg$genome)

#50_10_MAG,checkm1
setwd("D:/wenjian/research/defensome/groundwater/archaea/01MAG_info")
MAG_5010 <- readxl::read_xlsx("2269_metainfo.xlsx")

#70_10_MAG,checkm1
setwd("D:/wenjian/research/defensome/groundwater/archaea/01MAG_info")
MAG_7010 <- readxl::read_xlsx("zhong_70_10.xlsx")
#write.table(MAG_7010,file="mag_7010_zhong.txt",quote = F,sep = '\t',row.names = F)
#MAG_size_info,checkm2
#setwd("D:/wenjian/research/defensome/groundwater/archaea/01MAG_info")
#MAG_checkm2 <- read.csv("quality_report_checkm2.tsv",sep = '\t',check.names = F)

#合并不含DS的MAG
not_in <- MAG_5010$binname[!MAG_5010$binname %in% ar_dg$genome]
# 使用不在 ar_defense 行名中的 genome 值创建新的数据框
new_rows <- data.frame(
  genome = not_in,  # 重复 genome 值
  stringsAsFactors = FALSE)
# 为 new_rows 添加与 ar_defense 相同的列，并将这些列的值填充为 0
col_names <- setdiff(names(ar_dg), "genome")  # 获取 dss 中除了 genome 之外的列名
new_rows[col_names] <- lapply(col_names, function(col_name) rep(0, nrow(ar_dg)))
rownames(new_rows) <- new_rows$genome
#合并new_rows与dss
ar_dg <- rbind(ar_dg,new_rows)

#提取70_10的MAGs
ar_dg_7010 <- ar_dg[ar_dg$genome %in% MAG_7010$binname,]
sum(rowSums(ar_dg_7010[,-1])!=0) #956个MAGs中找到353种DGs，共5427个。

#10. 古菌的DGs数量分布####

#10.1数据预处理
rownames(ar_dg_7010) <- ar_dg_7010$genome
ar_dg_7010 <- ar_dg_7010[,-1]

# 使用rowSums()对每行求和
row_sums <- rowSums(ar_dg_7010)
# 将求和结果转换为向量
sum_vector <- as.vector(row_sums)
row_names <- rownames(ar_dg_7010)
# 创建一个新的数据框，包含求和结果和行名
ar_7010_num_dgs <- data.frame(Sum = sum_vector, row.names = row_names)

# 首先，确保"Sum"列的值是数值类型
ar_7010_num_dgs$Sum <- as.numeric(ar_7010_num_dgs$Sum)
# 使用table()函数统计Sum列每个数值的出现次数
sum_counts <- table(ar_7010_num_dgs$Sum)
# 将统计结果转换为数据框
sum_counts_df <- as.data.frame(sum_counts)
# 重命名列
colnames(sum_counts_df) <- c("Value", "Frequency")
# 确保Value列是因子类型
sum_counts_df$Value <- as.factor(sum_counts_df$Value)

#"#DFA4BA" "#E9AF89"  "#B3DFC4"
#10.2绘制柱状图
library(ggplot2)
# 使用geom_bar()绘制柱状图，指定stat为"identity"
dg_distribution <- ggplot(sum_counts_df,aes(x=Value,y=Frequency))+
  geom_bar(stat = "identity",position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,fill="#E9AF89"#DFA4BA #E9AF89  #B3DFC4       
  )+
  geom_text(aes(label = Frequency), angle = 90, hjust = -0.1,
            position = position_dodge(width = 0.9), 
            # vjust = -0.50, # 调整文本位置，避免与柱子重叠
            color = "black", # 文本颜色，根据柱子颜色调整
            size = 3.5)+ # 文本大小
  scale_y_continuous(limits = c(0,600),expand = c(0,0))+
  labs( x="Number of DGs",y="Number of MAGs") +
  theme(axis.text.x=element_text(angle = 90,vjust=0.5,hjust=0.5,color="black"),
        axis.text.y=element_text(color="black"),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.line.x=element_line(color="black"),
        axis.line.y.left = element_line(color="grey30"),
        axis.line.y.right = element_line(color="grey30"),
        axis.line.x.bottom = element_line(color="grey30"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))
dg_distribution  
# 导出PDF
setwd("D:/wenjian/research/defensome/groundwater/archaea/02defense")
write.csv(sum_counts_df,"FigS4D1.csv")
ggsave("1_DGs_distribution.pdf", plot = dg_distribution, width = 5.5, height = 4.5, units = "in")

#11. 古菌中DGs的出现率####

#11.1数据预处理
# 计算每一列非0数字所占的百分比并输出到新数据框
percent_nonzero <- colSums(ar_dg_7010 != 0) / nrow(ar_dg_7010) * 100
new_df <- data.frame(Column = names(percent_nonzero), Percentage = percent_nonzero)
# 按照值的大小对数据集进行降序排列
ds_persentage <- new_df[order(-new_df$Percentage), ]

#11.2绘制柱状图
library(ggplot2)
dg_chuxian <- ggplot(new_df[new_df$Percentage>=1,]) +
  geom_bar(aes(x = reorder(Column,-Percentage), y = Percentage),
           stat = "identity", position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,
           fill='#E9AF89'#E9AF89,sa#9CB9D7,ss  #DFA4BA,ws#E9AF89,sa
  )+
  scale_y_continuous(limits = c(0,30),expand = c(0,0))+
  labs( x=NULL,y="Proportion present in 1453 MAGs (%)") +
  theme(axis.text.x=element_text(angle = 90,vjust=0.5,hjust=0.9,color="black"),
        axis.text.y=element_text(color="black"),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        text = element_text(family = 'sans'),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.line.x=element_line(color="black"),
        axis.line.y.left = element_line(color="grey30"),
        axis.line.y.right = element_line(color="grey30"),
        axis.line.x.bottom = element_line(color="grey30"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))
dg_chuxian
write.csv(new_df,"FigS4D2.csv")
ggsave("1_DG_chuxian.pdf", plot = dg_chuxian, width = 6, height = 5, units = "in")

#12. 古菌中基因组大小统计，DSs数量、密度、多样性与MAG大小的关系#################

#12.1数据预处理
md_archaea <- merge(MAG_7010,MAG_checkm2,by.x = 'binname',by.y = 'Name')
ar_de_sum <- data.frame(ds_sum = rowSums(ar_defense_7010[,-1]),
                        genome = ar_defense_7010$genome,
                        richness = apply(ar_defense_7010[,-1], 1, function(row) sum(row > 0))
)
md_archaea <- merge(md_archaea,ar_de_sum,by.x = 'binname',by.y = 'genome')
md_archaea$ds_density <- md_archaea$ds_sum/md_archaea$size*1000000

md_archaea$magsize <- floor(md_archaea$size/1000000)


#12.2统计
library(dplyr)
df_summary <- md_archaea %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(ds_density),
    sd_value = sd(ds_density)
  )
df_summary2 <- md_archaea %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(ds_sum),
    sd_value = sd(ds_sum)
  )
df_summary3 <- md_archaea %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(richness),
    sd_value = sd(ds_sum)
  )

#计算两列数字的相关性及显著性
cor.test(df_summary$magsize, df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(df_summary2$magsize, df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(df_summary3$magsize, df_summary3$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)

#几个变量对DGs的逐步回归分析



#12.3绘图
p12_dsdensity <- ggplot(df_summary, aes(x = magsize, y = mean_value)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.1) +
  labs(x = 'MAG size (Mb)', y = 'DGs density (per Mb)') +
  annotate("Text",label='ρ = 0.59, P = 0.038',x=3,y=6.5,family='Times')+
  #scale_color_manual(values = c("SA" = "#E9AF89", "SS" = "#9CB9D7", "WA" = "#B3DFC4", "WS" = "#DFA4BA")) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),  # 横轴标签旋转
    text = element_text(family = "Times New Roman"),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p12_dsdensity
ggsave("s1_p3_DGs_density_magsize.pdf", plot = p3, width = 5.4, height = 3.42, units = "in")

p12_dsnum <- ggplot(df_summary2, aes(x = magsize, y = mean_value)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.1) +
  labs(x = 'MAG size (Mb)', y = 'DGs density (per Mb)') +
  annotate("Text",label='ρ = 0.59, P = 0.038',x=3,y=6.5,family='Times')+
  #scale_color_manual(values = c("SA" = "#E9AF89", "SS" = "#9CB9D7", "WA" = "#B3DFC4", "WS" = "#DFA4BA")) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),  # 横轴标签旋转
    text = element_text(family = "Times New Roman"),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p12_dsnum

p12_dsrich <- ggplot(df_summary3, aes(x = magsize, y = mean_value)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.1) +
  labs(x = 'MAG size (Mb)', y = 'DGs density (per Mb)') +
  annotate("Text",label='ρ = 0.59, P = 0.038',x=3,y=6.5,family='Times')+
  #scale_color_manual(values = c("SA" = "#E9AF89", "SS" = "#9CB9D7", "WA" = "#B3DFC4", "WS" = "#DFA4BA")) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),  # 横轴标签旋转
    text = element_text(family = "Times New Roman"),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p12_dsrich


#13. DHR分布图#########
#每个门的DSs总数除以该门MAGs数量
setwd("D:/wenjian/research/defensome/groundwater/archaea/01mag_info")
md_archaea <- read.csv("md_archaea.csv")

DHR <- md_archaea %>% group_by(col2) %>% 
  summarise(mean_number = mean(ds_sum))
#DHR$col2 <- gsub('^p__','',DHR$col2)

# 使用match函数找到md_archaea中col2列在DHR中col2列的对应位置
matches <- match(DHR$col2, md_archaea$col2)
# 然后，使用这些匹配的位置来将md_archaea中的superphylum值赋给DHR的phy_type列
DHR$phy_type <- md_archaea$superphylum[matches]

DHR$phy_type <- ifelse(DHR$phy_type == 'DPANN', 'DPANN', 'Others')

DHR <- DHR[DHR$col2!='' & DHR$mean_number>0,]

library(ggpubr)
dotp <- ggdotchart(DHR,x = 'col2', y = 'mean_number',
                   color = 'phy_type',                                # Color by groups
                   palette = c("#FC4E07","grey", "#E7B800", "#FC4E07"), # Custom color palette
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
  labs(x=NULL,y='DS number (per MAG)')+
  theme(
    panel.background = element_blank(),  # 背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    
    # 设置边框粗细一致
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 边界线颜色和粗细
    panel.spacing = unit(0.5, "lines"),  # 设置边框的粗细，与panel.border的size参数一致
    
    # 设置坐标轴标题和标签
    axis.title = element_text(face = "bold", size = 10, color = "grey20"),  # 坐标轴标题样式
    axis.text = element_text(size = 8, color = "grey20"),  # 坐标轴标签样式
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 坐标轴刻度线样式
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
write.csv(DHR,"FigS9C.csv")

#每个门的DSs数量除以MAGs大小之后再除以MAGs数量（density的平均值）
DHR2 <- md_archaea %>% group_by(col2) %>% 
  summarise(mean_density = mean(ds_density))
#DHR2$col2 <- gsub('^p__','',DHR2$col2)

# 使用match函数找到md_archaea中col2列在DHR中col2列的对应位置
matches <- match(DHR2$col2, md_archaea$col2)
# 然后，使用这些匹配的位置来将md_archaea中的superphylum值赋给DHR的phy_type列
DHR2$phy_type <- md_archaea$superphylum[matches]

DHR2$phy_type <- ifelse(DHR2$phy_type == 'DPANN', 'DPANN', 'Others')

DHR2 <- DHR2[DHR2$col2!='' & DHR2$mean_density>0,]

library(ggpubr)
dotp2 <- ggdotchart(DHR2,x = 'col2', y = 'mean_density',
                    color = 'phy_type',                                # Color by groups
                    palette = c("#FC4E07","grey", "#E7B800", "#FC4E07"), # Custom color palette
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
  labs(x=NULL,y='DS density (per MAG per Mb)')+
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
write.csv(DHR2,"FigS9D.csv")

library(cowplot)
g <- plot_grid(dotp, dotp2,
               labels = LETTERS[3:4], nrow = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g

setwd('D:/wenjian/research/defensome/groundwater/archaea/02defense')
ggsave(plot = g,'1_p13_DHR.pdf',width = 12,height = 3.09,units = 'in')

