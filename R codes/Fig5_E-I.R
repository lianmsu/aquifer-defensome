#DIs
#DIs的判定条件，只满足DGs大于等于5以及DSs大于等于3不行。需要考虑基因之间的位置！


library(dplyr)
library(stringr)


#1. 读取数据####
#1.1 DIs数据
setwd("D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense")
dis <- read.delim('defense_islands_summary4.csv',header = T,check.names = F)
#提取MAGs名称
dis <- dis %>%
  mutate(genome = str_extract(Contig, "^(.*?fasta)"))

#提取样品名称
dis <- dis %>%
  mutate(sample = str_extract(Contig, "^[^_]*")) 

#1.2 每个MAGs对应的物种信息
setwd('D:/wenjian/research/defensome/groundwater/02去冗余信息')
mags_tax <- read.delim2('gtdbtk.bac120.summary.tsv')
# 假设 df 是你的数据框，MAGs 是其中的一列
mags_tax$user_genome <- paste(mags_tax$user_genome, ".fasta", sep="")
mags_tax <- mags_tax[mags_tax$user_genome %in% dis$genome,]

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

#1.3合并物种信息与DIs信息
dis <- merge(dis,mags_tax,by.x = 'genome',by.y = 'user_genome')

#1.4 读取每个样点的信息
setwd('D:/wenjian/research/defensome/groundwater')
env <- readxl::read_xlsx('GW_DATA_lpw.xlsx')

#1.5 将env信息合并到dis
dis <- merge(dis,env,by.x = "sample","sample")


#2. 将CPR的DIs 的基因数量与其他门比较####

# 使用mutate和if_else增加新的CPR列
# dis_tax <- dis %>%
#   mutate(CPR = if_else(col2 == "p__Patescibacteria", "CPR", "non-CPR"))

dis_tax <- dis %>%
  mutate(CPR = case_when(col2 == "p__Patescibacteria" ~ "CPR",
                         col2 == "p__Pseudomonadota" ~ "Pseudomonadota",
                         TRUE ~ "Other phyla"
  ))




#比较CPR与Pseudomonadota,non-CPR的差异####

sum(dis_tax$CPR=="CPR") #n=25
sum(dis_tax$CPR!="CPR") #n=637
sum(dis_tax$CPR=="Other phyla") #n=340
sum(dis_tax$CPR!="Pseudomonadota") #n=365


#计算显著性差异
library(rstatix)
md_p_val1 <- dis_tax %>% 
  wilcox_test(`Total Gene Count` ~ CPR) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="CPR",dodge=0.5)
md_p_val1$p.adj #6.090e-05 2.115e-06 3.030e-01

#2.1 绘图####
library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means

md_p_val1$y.position <- c(13,14,15)
p1 <- ggplot(dis_tax, aes(x = CPR, y = `Total Gene Count`)) +
  geom_violin(trim = FALSE, aes(fill = CPR),show.legend = F) + # 绘制小提琴图
  geom_boxplot(aes(fill = CPR),width = 0.2,outliers = FALSE, staplewidth = 0.5) +  # 隐藏箱线图中的异常值点
  #geom_jitter(width = 0.2, alpha = 0.7)+  # 显示异常值点
  # stat_compare_means(comparisons = list(c("1", "2"),c("SS", "WA"),c("SS", "WS"),
  #                                       c("SA", "WA"),c("SA", "WS"),
  #                                       c("WA", "WS")
  # ),
  # method = "wilcox.test", test = "p.adjust",label = "p.signif")+
  stat_pvalue_manual(md_p_val1,label = "p.adj.signif",hide.ns=F,
                     tip.length = 0.05,label.size = 5,color="black"
                     )+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  #stat_summary(fun=mean, geom="line", aes(group=CPR), col="black")+
  scale_fill_manual(values = c("CPR"="#FF7F0E","Other phyla" = "#ADD9EE","Pseudomonadota"="#7892C7"
  )) +
  #ylim(0,15)+
  #scale_y_continuous(limits = c(0,10))+
  coord_cartesian(ylim = c(0,15))+
  labs(x=NULL,y='DGs number in DIs')+
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

setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
write.csv(dis_tax,"fig5f")
ggsave("8CPR_DIs_genenumber.pdf", plot = p1, width = 5.5, height = 3.8, units = "in")

#比较基因岛总长度的差异
dis_tax$allgene <- dis_tax$`Last Record Num` - dis_tax$`First Record Num`+1
dis_tax$dgbili <- dis_tax$`Total Gene Count`/dis_tax$allgene

#
md_p_val2 <- dis_tax %>% 
  wilcox_test(allgene ~ CPR) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="CPR",dodge=0.5)
md_p_val2$p.adj #0.0007470 0.0002907 0.9600000
#
md_p_val3 <- dis_tax %>% 
  wilcox_test(dgbili ~ CPR) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="CPR",dodge=0.5)
md_p_val3$p.adj #0.213 0.192 1.000

md_p_val2$y.position <- c(26,28,30)
p2 <- ggplot(dis_tax, aes(x = CPR, y = allgene)) +
  geom_violin(trim = FALSE, aes(fill = CPR),show.legend = F) + # 绘制小提琴图
  geom_boxplot(aes(fill = CPR),width = 0.2,outliers = FALSE, staplewidth = 0.5) +  # 隐藏箱线图中的异常值点
  #geom_jitter(width = 0.2, alpha = 0.7)+  # 显示异常值点
  # stat_compare_means(comparisons = list(c("1", "2"),c("SS", "WA"),c("SS", "WS"),
  #                                       c("SA", "WA"),c("SA", "WS"),
  #                                       c("WA", "WS")
  # ),
  # method = "wilcox.test", test = "p.adjust",label = "p.signif")+
  stat_pvalue_manual(md_p_val2,label = "p.adj.signif",hide.ns=F,
                     tip.length = 0.05,label.size = 5,color="black"
  )+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  #stat_summary(fun=mean, geom="line", aes(group=CPR), col="black")+
  scale_fill_manual(values = c("CPR"="#FF7F0E","Other phyla" = "#ADD9EE","Pseudomonadota"="#7892C7"
  )) +
  #ylim(0,15)+
  #scale_y_continuous(limits = c(0,10))+
  coord_cartesian(ylim = c(0,30))+
  labs(x=NULL,y='DI length (genes)')+
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
write.csv(dis_tax,"fig5g.csv")
ggsave("8DI-length-(genes).pdf", plot = p2, width = 5.5, height = 3.8, units = "in")

md_p_val3$y.position <- c(1.6,1.8,2.0)
p3 <- ggplot(dis_tax, aes(x = CPR, y = dgbili)) +
  geom_violin(trim = FALSE, aes(fill = CPR),show.legend = F) + # 绘制小提琴图
  geom_boxplot(aes(fill = CPR),width = 0.2,outliers = FALSE, staplewidth = 0.5) +  # 隐藏箱线图中的异常值点
  #geom_jitter(width = 0.2, alpha = 0.7)+  # 显示异常值点
  # stat_compare_means(comparisons = list(c("1", "2"),c("SS", "WA"),c("SS", "WS"),
  #                                       c("SA", "WA"),c("SA", "WS"),
  #                                       c("WA", "WS")
  # ),
  # method = "wilcox.test", test = "p.adjust",label = "p.signif")+
  stat_pvalue_manual(md_p_val3,label = "p.adj.signif",hide.ns=F,
                     tip.length = 0.05,label.size = 5,color="black"
  )+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  #stat_summary(fun=mean, geom="line", aes(group=CPR), col="black")+
  scale_fill_manual(values = c("CPR"="#FF7F0E","Other phyla" = "#ADD9EE","Pseudomonadota"="#7892C7"
  )) +
  #ylim(0,15)+
  #scale_y_continuous(limits = c(0,10))+
  coord_cartesian(ylim = c(0,2))+
  labs(x=NULL,y='DG proportion in DIs')+
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
write.csv(dis_tax,"fig5h.csv")
ggsave("8DGs-proportion-in-DIs.pdf", plot = p3, width = 5.5, height = 3.8, units = "in")

library(cowplot)
g1 <- plot_grid(p1, p2,p3,labels = LETTERS[1:3], nrow = 1)
g1

ggsave(plot = g1,'8CPR_pseu_men_new.pdf',width = 7.9,height = 3.2,units = 'in')
library(export)
graph2ppt(file="8.pptx", width=10, height=6.18)

#3. 比较不同地质环境区的情况####

#3.1 基本统计

sum(dis_tax$Q=="1") #n=133
sum(dis_tax$Q=="2") #n=285
sum(dis_tax$Q=="3") #n=30
sum(dis_tax$Q=="4") #n=122
sum(dis_tax$Q=="5") #n=36
sum(dis_tax$Q=="6") #n=24
sum(dis_tax$Q=="7") #n=32

#计算显著性差异
library(rstatix)
md_p_val1 <- dis_tax %>% 
  wilcox_test(Total.Gene.Count ~ Q) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="Q",dodge=0.5)
md_p_val1$p.adj #ns

#转换为因子类型
dis_tax$Q <- as.character(dis_tax$Q)
dis_tax$Q <- ordered(dis_tax$Q,
                             levels = c("1", "2", "3","4","5","6","7")) # 替换为实际的水平名称
#3.2 绘图####
library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means


p2 <- ggplot(dis_tax, aes(x = Q, y = Total.Gene.Count)) +
  geom_boxplot(aes(fill = Q),outliers = F) +  # 隐藏箱线图中的异常值点
  #geom_jitter(width = 0.2, alpha = 0.7)+  # 显示异常值点
  # stat_compare_means(comparisons = list(c("1", "2"),c("SS", "WA"),c("SS", "WS"),
  #                                       c("SA", "WA"),c("SA", "WS"),
  #                                       c("WA", "WS")
  # ),
  # method = "wilcox.test", test = "p.adjust",label = "p.signif")+
  stat_pvalue_manual(md_p_val1,label = "p.adj.signif",hide.ns=T,
                     tip.length = 0.05,label.size = 5,color="black",
                     y.position = 15
  )+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  stat_summary(fun=mean, geom="line", aes(group=Q), col="black")+
  scale_fill_manual(values = c("1"="#FCB2AF","2" = "#9BDFDF",
                               "3" = "#FFE2CE", "4" = "#C4D8E9",
                               "5"="#BEBCDF","6"="#FB8C62","7"="#FFF3CA"
  )) +
  #ylim(0,15)+
  #scale_y_continuous(limits = c(0,10))+
  coord_cartesian(ylim = c(0,15))+
  labs(x=NULL,y='DI length (genes)')+
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

setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
ggsave("8dizhihuanjing_DIs_genenumber.pdf", plot = p2, width = 5.5, height = 3.8, units = "in")


#4. 比较DI中DGs占的比例####
#4.1 基本统计
#
dis$allgene <- dis$Last.Record.Num-dis$First.Record.Num+1

dis$dgbili <- dis$Total.Gene.Count/dis$allgene

mean(dis$dgbili) #0.6781

#计算所有DIs中，各类DSs的数量所占比例
alldss <- colSums(dis[, 9:106]) %>% as.data.frame() %>% setNames("dssnumber")
sum(dis$allgene) #7040
sum(alldss) #2091
# 添加新行
#alldss <- rbind(alldss, data.frame(dssgene = 2444, stringsAsFactors = FALSE))

# 给新行命名
#rownames(alldss)[nrow(alldss)] <- "NA"


# 筛选出dssnumber 小于50的行，并将它们合并为一行
alldss_filtered <- alldss %>%
  filter(dssnumber < 50) %>%
  summarise(combined_dssgene = sum(dssnumber)) %>% setNames("dssnumber")

# 将合并后的行添加到原始数据框的末尾
alldss <- rbind( alldss,alldss_filtered)
rownames(alldss)[nrow(alldss)] <- "Others"
alldss$ds <- rownames(alldss)

#绘制圈图
library(ggplot2)
library(ggpmisc)

# 绘制中空饼图
mycolor=c("#7294D4","#FB9A99","#FDBF6F","#A6CEE3","#1F78B4","#CAB2D6",
          "#FF7F00","#0B775E",
          "#33A02C","#6A3D9A","#B2DF8A")
pie <- ggplot(alldss[alldss$dssnumber>50,], aes(x = 2, y = dssnumber, fill = ds)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  # 使用 geom_bar 模拟饼图
  coord_polar(theta = "y") +  # 将坐标转换为极坐标
  xlim(0.5,2.5)+
  scale_fill_manual(values = mycolor) +
  geom_text( aes(label = dssnumber), 
             position = position_stack(vjust = 0.5),
             color = "white", # 文本颜色
             size = 5) + # 文本大小
  guides(fill = guide_legend(title = "DSs"))+
  theme_void() +  # 移除背景元素
  theme(legend.position = "none",  # 移除默认图例
        axis.text.x = element_blank(),
        axis.ticks = element_blank())  # 移除坐标轴标签和刻度
pie
# 创建自定义图例
legend_plot <- ggplot(alldss[alldss$dssnumber>50,], aes(x = 2, y = dssnumber, fill = ds)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mycolor) +
  theme_void() +
  guides(fill = guide_legend(title = ""))

# 提取图例
legend <- get_legend(legend_plot)

# 创建组合图
combined_plot <- ggdraw() +
  draw_plot(pie, 0, 0, 1, 1) +
  draw_grob(legend, 0.4, 0.4, 0.2, 0.25)  # 调整图例的位置和大小
combined_plot

setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
write.csv(alldss[alldss$dssnumber>50,],"fig5e.csv")
ggsave(plot = combined_plot,"8DIs_dss_pie.pdf",width = 7, height = 7, units = "in")

ggsave(plot = g1,'8CPR_pseu_men.pdf',width = 10,height = 6.18,units = 'in')

graph2ppt(file="8.pptx", width=10, height=6.18)

#5. 对DSs进行相关分析####
alldss <- dis[, 9:106]

# 计算相关系数矩阵，只计算数值型列
cor_matrix <- cor(alldss[sapply(alldss, is.numeric)])
library(corrplot)
library(ggplot2)
env.cor <- round(cor(alldss, method = "spearman"),3) # round(),对输出结果取小数点前三位,计算pearson相关性系数
env.cor
library(ggcorrplot)
env.p <-round(cor_pmat(alldss,method = "spearman"),3) # pearson系数
env.p

par(mfrow=c(1,1))
plot.new()
p1<- corrplot(corr =env.cor, p.mat = env.p,method = "square",
              type="upper",diag=FALSE,add=TRUE,
              tl.pos="n", tl.cex=0.5, tl.col="black",tl.srt = 0,tl.offset=0.5,
              insig="label_sig",sig.level = c(.001, .01, .05),
              pch.cex = 0.8,pch.col = "red") # tl.*参数依次设置变量标签的位置为左侧和右侧，字符大小为1，颜色为黑色，倾斜45度，距离热图边界0.5。pch*参数依次设置显著性标签大小为0.8，颜色为红色。
p2<-corrplot(corr = env.cor, method = "number",
             type="lower",
             tl.pos = "tp",tl.cex=0.5, tl.col="black",tl.srt = 90,tl.offset=0.5,
             cl.pos = "n",
             number.digits = 2,number.cex = 0.5,number.font = NULL
) # number*参数依次设置保留相关性系数3位小数点，字体大小为0.7，字体使用默认，只能设置par中的字体参数。后面会讲修改par字体参数。

corrplot(corr =env.cor, p.mat = env.p,method = "square",
         type="upper",diag=FALSE,#add=TRUE,
         tl.pos="d", tl.cex=0.5, tl.col="black",tl.srt = 0,tl.offset=0.5,
         insig="label_sig",sig.level = c(.001, .01, .05),
         pch.cex = 0.8,pch.col = "red")






#6. 构建多个DSs之间的网络图####
library(igraph)
library(Hmisc) #提供rcorr函数
library(ggplot2)
library(ggpubr)
library(ggsignif)
alldss <- as.matrix(alldss)

#计算相关性系数；
sp.cor<- rcorr(alldss,type="spearman") #这里计算的是所有物种之间的
#提取r、p值矩阵；
r.cor<-sp.cor$r
p.cor<-sp.cor$P
#使用Benjamini-Hochberg("FDR-BH")法进行多重检验校正；
p.adj <- p.adjust(p.cor, method="BH")
#确定物种间存在相互作用关系的阈值，将相关性R矩阵内不符合的数据转换为0
#r.cor[p.cor>=0.05|abs(r.cor)<0.4] = 0
r.cor[p.cor>=0.05] = 0
#对角线处的1不计
diag(r.cor) <- 0

#使用邻接矩阵（即相关系数矩阵）创建网络；
net_allshu<-graph.adjacency(r.cor,weight=T,mode="undirected")
#去掉冗余的边（multiple edges、loop edges）；
net_allshu<-simplify(net_allshu)
#提取权重
df_weight = E(net_allshu)$weight
# 设定边的宽度，这里我们将相关系数与边宽进行关联
E(net_allshu)$width = abs(df_weight)*5
#生成网络图的结点标签（OTU id）和degree属性；degree表示一个点的中心性
V(net_allshu)$label <- V(net_allshu)$name
V(net_allshu)$degree <- degree(net_allshu)

V(net_allshu)$type <- c(rep('aroma',times=119),rep('arg',times=67),
                        rep('mrg',times=47)
)

#查看网络图的对象结构;
print(net_allshu)

#将网络图导出为"graphml"、"gml"格式，方便导入Gephi中使用；
setwd("D:/wenjian/research/defensome/groundwater/04mag_DI_net")
write_graph(net_allshu, "DIsds_onlyp0.05.graphml", format="graphml")


#7.统计检验防御岛中不同defense families之间的共定位####
setwd("D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense")

#7.1 读取1. 中的数据
di_suzhu <- table(mags_tax$col2)

#7.2 卡方检验及绘图，这里可以比较防御岛与非防御岛上DSs的O/E
di_ds <- dis[,c(9:106)] %>% colSums()
di_ds_chi <-  chisq.test(di_ds)
di_ds_oe <- di_ds_chi$observed/di_ds_chi$expected %>% as.data.frame() %>%
  setNames("O/E")

setwd("D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense")
non_di_ds <- read.csv("ds_counts_matrix.csv",check.names = F)
di_ds <- di_ds %>% as.data.frame() %>% setNames('DIs')
non_di_ds <- non_di_ds[!(non_di_ds$genome %in% dis$genome), ]
no_di_ds <- non_di_ds[,-1] %>% colSums() %>% as.data.frame() %>% setNames('non DIs')

di_ds <- merge(di_ds,no_di_ds,by = 0)
rownames(di_ds) <- di_ds$Row.names
di_ds <- di_ds[,-1]
di_ds_chi <- chisq.test(di_ds)
di_ds_oe <- di_ds_chi$observed/di_ds_chi$expected %>% as.data.frame() %>%
  setNames(c("DIs","non DIs"))

#绘制热图
library(ComplexHeatmap)
library(circlize)  # 用于 colorRamp2 函数

# 使用subset函数选择满足条件的行
filtered_rows <- subset(di_ds_oe, `O/E` > 0.1)
# 定义颜色映射函数
cairo_pdf(width = 2.5,height = 8,filename = "8_di_ds_oe2.pdf")
col_fun <- colorRamp2(c(0,1, 2), c("white","#B6B8D8", "#FF7F0E"))  # 粉色到绿色的渐变
p_di_ds_oe <- Heatmap(di_ds_oe,
                       col = col_fun,
                       cluster_rows = F,
                       cluster_columns = F,
                       name = 'O/E',
                       row_title = NULL,
                      # row_split = colnames(di_ds_oe),
                       row_names_side = "left", # 行标题显示在左侧
                      row_names_gp = gpar(fontsize = 6),
                      # column_names_gp = gpar(fontsize = 6),
                       column_names_rot = 30,
                       rect_gp = gpar(col = "white", lwd = 0.5),    ##热图矩阵线设置
                       heatmap_legend_param = list(
                         title = "O/E", at = c(0,1, 2), 
                         labels = c("0","1", "≥2")),
)
p_di_ds_oe
setwd("D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense")
dev.off()



cairo_pdf(width = 12,height = 2.3,filename = "8_di_ds_oe3.pdf")
col_fun <- colorRamp2(c(0,1, 2), c("white","#B6B8D8", "#FF7F0E"))  # 粉色到绿色的渐变
p_di_ds_oe <- Heatmap(t(di_ds_oe),
                      col = col_fun,
                      cluster_rows = F,
                      cluster_columns = T,
                      name = 'O/E',
                      row_title = NULL,
                       row_split = colnames(di_ds_oe),
                      row_names_side = "left", # 行标题显示在左侧
                     # row_names_gp = gpar(fontsize = 10),
                       column_names_gp = gpar(fontsize = 8),
                      column_names_rot = 60,
                      rect_gp = gpar(col = "white", lwd = 0.5),    ##热图矩阵线设置
                      heatmap_legend_param = list(
                        title = "O/E", at = c(0,1, 2), 
                        labels = c("0","1", "≥2")),
)
p_di_ds_oe
setwd("D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense")
dev.off()

write.csv(di_ds_oe,"fig5i.csv")


# library(ggplot2)
# ggsave("8_di_ds_oe.pdf",p_di_ds_oe,width = 1.5,height = 7,units = 'in')
# library(ggpubr)
# ggexport(filename="8_di_ds_oe.pdf",p_di_ds_oe,width = 1.5,height = 7,units = 'in')

#7.3防御岛上不同defense family之间的共定位？
#根据fisher检验的原理，自己构造防御岛 非防御岛\n RM Cas，这种2x2表格 计算OR值？
#读取所有DSs信息
dss <- read.csv("ds_counts_matrix.csv",check.names = F)
dss_all_sum <- dss[,-1] %>% colSums() %>% as.data.frame() %>% setNames("DS")
dss_all_sum$DF <- rownames(dss_all_sum)
di_ds <- dis[,c(9:106)] %>% colSums() %>% as.data.frame() %>% setNames("DI_DS")
di_ds$DF <- rownames(di_ds)

allmag_di_ds <- merge(dss_all_sum,di_ds,by ="DF" ,all = T)
allmag_di_ds[is.na(allmag_di_ds)] <- 0

allmag_di_ds$`non DIs` <- allmag_di_ds$DS - allmag_di_ds$DI_DS
rownames(allmag_di_ds) <- allmag_di_ds$DF
f_mag_di_ds <- select(allmag_di_ds,DI_DS,`non DIs`)
f_mag_di_ds <- f_mag_di_ds[f_mag_di_ds$DI_DS>0,]

#一次性对所有数据进行fisher，超内存
library(rstatix)
df_fisher <- fisher.test(f_mag_di_ds,alternative = "two.sided",
            workspace = 1e+9) #超内存
df_fisher <- fisher.test(f_mag_di_ds[c(4,5),],alternative = "two.sided",
                         workspace = 1e+9)

calculate_fisher_tests <- function(f_mag_di_ds) {
  # 初始化一个空的数据框来存储结果
  results_df <- data.frame(
    Row1 = character(),
    Row2 = character(),
    Odds_Ratio = numeric(),
    P_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 获取f_mag_di_ds的行数
  n_rows <- nrow(f_mag_di_ds)
  
  # 对每两行进行循环
  for (i in 1:(n_rows)) {
    for (j in (i):n_rows) {
      # 提取当前行和下一行的行名
      row_names <- c(rownames(f_mag_di_ds)[i], rownames(f_mag_di_ds)[j])
      
      # 提取当前行和下一行的数据
      data_to_test <- f_mag_di_ds[c(i, j), ]
      
      # 确保data_to_test是一个2x2矩阵
      if (ncol(data_to_test) == 2) {
        # 执行Fisher精确检验
        fisher_result <- fisher.test(data_to_test, alternative = "two.sided", workspace = 1e+9)
        
        # 将行名、odds ratio和p.value添加到结果数据框中
        results_df <- rbind(results_df, data.frame(
          Row1 = row_names[1],
          Row2 = row_names[2],
          Odds_Ratio = fisher_result$estimate[1],
          P_Value = fisher_result$p.value,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(results_df)
# 示例调用
# 假设f_mag_di_ds是一个两列的数据框
# f_mag_di_ds <- data.frame(
#   A = c(1, 2, 3, 4),
#   B = c(5, 6, 7, 8),
#   row.names = c("Row1", "Row2", "Row3", "Row4")
# )
# results <- calculate_fisher_tests(f_mag_di_ds)
# print(results)
}

results_df <- calculate_fisher_tests(f_mag_di_ds)


library(tidyr)
tb_od = pivot_wider(results_df[,1:3], names_from = c('Row2'), values_from = 'Odds_Ratio') %>% as.data.frame()
#kable(tb_wide_new)
rownames(tb_od) <- tb_od$Row1
tb_od <- tb_od[,-1]
tb_od[is.na(tb_od)] <- 0

try <- tb_od+t(tb_od)
diag(try) <- 1

library(ComplexHeatmap)
library(circlize)  # 用于 colorRamp2 函数
col_fun <- colorRamp2(c(-2,2,6), c("green", "white", "red"))
DS_or <- Heatmap(try,
        col = col_fun,
        name = 'OR',
        cluster_rows = T,
        cluster_columns = T,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(
          title = "OR", at = c(-2,2,6), 
          labels = c("≤-2","2", "≥6")),
         )

library(ggplot2)
setwd("D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense")
cairo_pdf(width = 8,height = 8,filename = "8_DS_OR.PDF")
DS_or
dev.off()



######绘制p-value
library(tidyr)
tb_p = pivot_wider(results_df[,c(1,2,4)], names_from = c('Row2'), values_from = 'P_Value') %>% as.data.frame()
#kable(tb_wide_new)
rownames(tb_p) <- tb_p$Row1
tb_p <- tb_p[,-1]
tb_p[is.na(tb_p)] <- 0

p_value <- tb_p+t(tb_p)
diag(p_value) <- 1

library(ComplexHeatmap)
library(circlize)  # 用于 colorRamp2 函数
col_fun2 <- colorRamp2(c(0,0.05,0.1), c("blue", "white", "pink"))
DS_p_value <- Heatmap(p_value,
        col = col_fun2,
        name = 'P-value',
        cluster_rows = T,
        cluster_columns = T,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(
          title = "P-value", at = c(0,0.05,0.1), 
          labels = c("0","0.05", "≥0.1")),
        )

setwd("D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense")
cairo_pdf(width = 8,height = 8,filename = "8_DS_P-value.pdf")
DS_p_value
dev.off()




