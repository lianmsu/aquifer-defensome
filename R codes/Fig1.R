
#Fig1A########################

library(ggplot2)
p <- ggplot(new_df[new_df$Percentage>=1,]) +
  geom_bar(aes(x = reorder(Column,-Percentage), y = Percentage),
           stat = "identity", position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,
           fill='#81D0D6'#E9AF89,sa#9CB9D7,ss  #DFA4BA,ws#E9AF89,sa
  )+
  scale_y_continuous(limits = c(0,50),expand = c(0,0))+
  labs( x=NULL,y="Proportion present in 27578 MAGs (%)") +
  theme(axis.text.x=element_text(angle = 90,vjust=0.5,hjust=0.9,color="black"),
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
p
ggsave("DSs_percentage.pdf", plot = p, width = 10, height = 8, units = "in")


#Fig1B and Fig1C########################
library(reshape2)
library(dplyr)
library(ggplot2)

setwd("D:/wenjian/research/defensome/groundwater/88defense_finder_database")
DG_md <- read.table("uniq_DG.txt",header = T)


DG_md$novel <- (DG_md$novel_DG)/DG_md$DGs *100
DG_md$unnovel <- (DG_md$DGs - DG_md$novel_DG)/DG_md$DGs *100

df_plot <- select(DG_md,source,novel,unnovel)
df_plot <- melt(df_plot)

df_plot$variable <- factor(df_plot$variable,levels = c('unnovel','novel'))

colorsss <- c(
  "#EFEFEF",  # 青色
  "#66CCFF",  # 橙色
  # "#F5E320",  # 红色
  # "#CC281C",  # 绿色
  "#D1CB7E",  # 紫色
  "#5AAAB2",  # 棕色
  "#E77F5A",  # 粉色
  "#A0E6C5",  # 灰色
  "#70A47F",  # 橄榄色
  "#1f77b4",  # 蓝色
  "#D68112",  # 浅蓝色
  "#C9A79E",  # 浅橙色
  "#98df8a",  # 浅绿色
  #  "#ff9896",  # 浅红色
  "#6C6B6F",  # 浅紫色
  "#CED2DD",  # 浅棕色
  #  "#f7b6d2",  # 浅粉色
  "#4E616A",  # 浅橄榄色
  #  "#9edae5",  # 浅青色
  "#ED6764",  # 深蓝色
  "#74A0A4"   # 深绿色
)

dg_novel <-  
  ggplot(df_plot, aes(x = source, y = value, fill = variable)) +
  geom_bar(stat = "identity", width = 0.8) +  # 使用实际值绘制柱状图
  scale_y_continuous(limits = c(0, 100)) +  # 设置 y 轴的范围和刻度
  #scale_y_continuous(limits=c(0,100)) +  # 将y轴标签转换为百分比
  scale_fill_manual(values = colorsss)+
  # scale_x_discrete(limits=c('Cd11','Cr11','Cu11','phOH11','PNP11','PAP11','CK11',
  #                           'Cd21','Cr21','Cu21','phOH21','PNP21','PAP21','CK21',
  #                           'Seed'))+
  #coord_flip() +  # 翻转坐标轴
  #scale_y_continuous(limits=c(0,100),expand = c(0,0)) +  # 将y轴标签转换为百分比
  labs(x = "", y = "Relative proportion (%)", fill = "") +  # 设置轴标签和图例标题
  
  #  geom_vline(xintercept = c(7.5,14.5),linetype=2)+
  guides(fill=guide_legend(title = "",nrow =1,position = 'bottom'))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "black", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'bottom',
    #legend.position = c(0.2,0.2),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "black"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "black"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    #legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dg_novel

setwd("D:/wenjian/research/defensome/groundwater/88defense_finder_database")
ggsave("dg_novel.pdf",dg_novel,width = 2.5, height = 4.7, units = 'in')


gw_dg_phyla <- data.frame(type=c('groundwater','RefSeq'),
                          value=c(152,48)
)
colorsss <- c(
  "#FF8F8F",  # 青色
  "#EFEFEF"# 橙色
)
dg_phyla <-  
  ggplot(gw_dg_phyla, aes(x = type, y = value, fill = type)) +
  geom_bar(stat = "identity", width = 0.8) +  # 使用实际值绘制柱状图
  #scale_y_continuous(limits = c(0, 100)) +  # 设置 y 轴的范围和刻度
  #scale_y_continuous(limits=c(0,100)) +  # 将y轴标签转换为百分比
  scale_fill_manual(values = colorsss)+
  # scale_x_discrete(limits=c('Cd11','Cr11','Cu11','phOH11','PNP11','PAP11','CK11',
  #                           'Cd21','Cr21','Cu21','phOH21','PNP21','PAP21','CK21',
  #                           'Seed'))+
  #coord_flip() +  # 翻转坐标轴
  #scale_y_continuous(limits=c(0,100),expand = c(0,0)) +  # 将y轴标签转换为百分比
  labs(x = "", y = "Phylum number", fill = "") +  # 设置轴标签和图例标题
  
  #  geom_vline(xintercept = c(7.5,14.5),linetype=2)+
  guides(fill=guide_legend(title = "",nrow =1,position = 'bottom'))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "black", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'bottom',
    #legend.position = c(0.2,0.2),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "black"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "black"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    #legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dg_phyla
setwd("D:/wenjian/research/defensome/groundwater/88defense_finder_database")
ggsave("dg_phyla.pdf",dg_phyla,width = 2.5, height = 4.7, units = 'in')

#Fig1D####################
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

#2. 计算每个MAGs的DSs数量
# 使用rowSums()对每行求和
row_sums <- rowSums(dss)

# 将求和结果转换为向量
sum_vector <- as.vector(row_sums)
row_names <- rownames(dss)

# 创建一个新的数据框，包含求和结果和行名
num_dss <- data.frame(Sum = sum_vector, row.names = row_names)

# 假设new_df是你的数据框，其中包含一个名为"Sum"的列

# 首先，确保"Sum"列的值是数值类型
num_dss$Sum <- as.numeric(num_dss$Sum)
# 使用table()函数统计Sum列每个数值的出现次数
sum_counts <- table(num_dss$Sum)
# 将统计结果转换为数据框
sum_counts_df <- as.data.frame(sum_counts)
# 重命名列
colnames(sum_counts_df) <- c("Value", "Frequency")

# 确保Value列是因子类型
sum_counts_df$Value <- as.factor(sum_counts_df$Value)

#3.绘制柱状图
library(ggplot2)
# 使用geom_bar()绘制柱状图，指定stat为"identity"
p <- ggplot(sum_counts_df,aes(x=Value,y=Frequency))+
  geom_bar(stat = "identity",position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,fill="#81D0D6"#DFA4BA #E9AF89  #B3DFC4       
  )+
  geom_text(aes(label = Frequency), angle = 60, hjust = -0,vjust = -0.5,
            position = position_dodge(width = 0.9), 
            # vjust = -0.50, # 调整文本位置，避免与柱子重叠
            color = "black", # 文本颜色，根据柱子颜色调整
            size = 3.5)+ # 文本大小
  scale_y_continuous(limits = c(0,8500),expand = c(0,0))+
  labs( x="Number of DSs",y="Number of MAGs") +
  theme(axis.text.x=element_text(angle = 45,vjust=0.5,hjust=0.5,color="black"),
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
p  
# 导出PDF
ggsave("1_DSs_distribution.pdf", plot = p, width = 10, height = 8, units = "in")

#Fig1E################################
library(ggplot2)
#1. 读取
setwd("D:/wenjian/research/defensome/RefSeq")
refseq_ds <- read.csv("df-RefSeq.csv")


compare <- read.csv("DS_compare2.txt",sep = '\t') #地下水数据目前只放了细菌的结果
compare$DS_mean <- compare$DSs_number/compare$MAG_number

mytheme <-   theme(
  panel.background = element_blank(),  # 设置背景为透明
  panel.grid.major = element_blank(),  # 去掉主要网格线
  panel.grid.minor = element_blank(),  # 去掉次要网格线
  panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
  axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
  axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
  axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
  axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
  legend.position = 'none',
  #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
  legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
  legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
  legend.background = element_blank(),  # 设置图例背景为透明
  legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
  legend.key = element_blank()  # 去掉图例键的背景
)

#2. 绘制横着的柱状图
# 绘制横条形图
p_compare <- ggplot(compare, aes(x = reorder(env,DS_mean) ,y = DS_mean)) +
  geom_bar(stat = "identity", aes(fill = DS_mean),width = 0.5) +  # 使用stat = "identity"因为我们的值已经是计数或比例
  coord_flip() +  # 翻转坐标轴，使条形图水平
  scale_fill_gradient(low = "#FFF1EB",
                      high = "#ACE0F9")+
  #scale_fill_viridis_c(option = "inferno") + ggtitle("'inferno'")+
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "black")+
  geom_text(aes(label = round(DS_mean,2)),vjust = 0.5,hjust = 1.1,size = 5)+
  labs(x = "", y = "DSs number (per MAG)", title = "")+  # 设置标签和标题
  mytheme
p_compare
setwd("D:/wenjian/research/defensome/RefSeq")
ggsave("p_compare2.pdf",p_compare,width = 4.3,height = 3.7,units = 'in')

#Fig1F######################
#1.
setwd('D:/wenjian/research/defensome/groundwater/01defensefinder_info')
allds <- read.delim2('df_List_Systems_lpw.csv',sep = ',')

#2.
# 检查Effector列中的空字符串并替换为"Unknown"
allds$Effector[allds$Effector == ""] <- "Unknown"

#3.绘制饼图
# 绘制饼图
library(ggplot2)

allds$Effector <- as.factor(allds$Effector)

# 定义颜色
my_colors <- c("#FF6384", "#36A2EB", "#FFCE56", "#4BC0C0", "#C0C0C0")


p <- ggplot(allds, aes(x = "", fill = Effector)) + 
  geom_bar(width = 1, color = "white") + 
  coord_polar("y", start = 0) +
  scale_fill_manual(values = my_colors) +
  geom_text(stat = "count", aes(label = ..count..), 
            position = position_stack(vjust = 0.5),
            color = "white", # 文本颜色
            size = 5) + # 文本大小
  theme_void()  # 移除所有非数据元素
p
ggsave(plot = p,"DSseffectortype.pdf",width = 5, height = 3.5, units = "in")

#MAGs中的DSs数量统计结果与DSs effector类型结合，treemap##
#1. 读取DSs数据
#读取所有27578个mags的数据
setwd('D:/wenjian/research/defensome/groundwater/02去冗余信息')
allmags <- read.delim2('Wdb_ANI99.csv',sep = ',',
                       quote = "",header = T
) 

#读取含有dss的mags的数据
setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
dss <- read.delim2('ds_counts_matrix.csv',sep = ',',row.names = 1,
                   quote = "",header = T
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

#2. 计算每个DSs在所有MAGs中的平均number
library(dplyr)
# 假设dss是你要处理的数据框,计算每列均值
means_df <- as.data.frame(colMeans(dss))
means_df$dsstype <- rownames(means_df)

#将means_df dsstype列的.改为-
means_df$dsstype <- gsub("\\.", "-", means_df$dsstype)

#合并
num_average <- merge(means_df,allds,by.x = 'dsstype',by.y = 'title')
colnames(num_average)[2] <- 'num_means'

# 找出vector1中在vector2中不存在的字符串
#unique_to_vector1 <- means_df$dsstype[!means_df$dsstype %in% num_average$dsstype]

library(treemap)
library(RColorBrewer)
library(treemapify)


# 添加多个分组信息
p1 <- ggplot(num_average, aes(area = num_means, fill = Effector,
                              label = dsstype,
                              subgroup = Effector)) +
  geom_treemap(colour='black') +
  geom_treemap_subgroup_border() +
  geom_treemap_subgroup_text(place = "centre", alpha = 0.5, 
                             colour = "white", fontface = "italic",family='Times') +
  geom_treemap_text(colour = "black", place = "topleft", reflow = T,family='Times') +
  scale_fill_manual(values = c("#FF6384", "#36A2EB", "#FFCE56", "#4BC0C0", "#C0C0C0"))
p1

setwd('D:/wenjian/research/defensome/groundwater/01defensefinder_info')

write.csv(num_average,"Fig1F.csv",quote = F)
ggsave(plot = p1,"01mean_dsnum_treemap.pdf",width = 8, height = 5.5, units = "in")


#计算未知effector type DS的比例
sum(num_average[num_average$Effector == "Nucleic acid degrading",]$num_means)
sum(num_average[num_average$Effector == "Nucleic acid degrading",]$num_means)/sum(num_average$num_means)

sum(num_average[num_average$Effector == "Unknown",]$num_means) #1.2
sum(num_average[num_average$Effector == "Unknown",]$num_means)/sum(num_average$num_means) #36.9%

#Fig1G################################
setwd('D:/wenjian/research/defensome/groundwater/03mag_info')
magsize <- read.csv('genome_stats.csv')
# 在R中，可以使用sub函数去掉数据框某一列字符串中的.fa后缀
# 假设数据框为df，列名为filename，可以使用以下代码：
magds <- data.frame(genome=rownames(dss),dss=rowSums(dss))

ds_size <- merge(magds,magsize,by.x = 'genome',by.y = 'genome')

ds_density <- data.frame(genome=ds_size$genome,
                         density=ds_size$dss/ds_size$size*1000000,
                         number=ds_size$dss,
                         length_Mb=ds_size$size/1000000
)
# 在R中，可以使用floor函数将某一列的数字取小数点前的整数，并形成新的一列
# 假设数据框为df，列名为original_column，可以使用以下代码：
ds_density$magsize <- floor(ds_density$length_Mb)

# 在R中，可以使用dplyr库计算数据框中各个分组对应数据的平均值和标准偏差
# 假设有数据集df，包含group和value列，可以使用以下代码：
library(dplyr)
df_summary <- ds_density %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(density),
    sd_value = sd(density)
  )
df_summary2 <- ds_density %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(number),
    sd_value = sd(number)
  )

#计算两列数字的相关性及显著性
cor.test(df_summary$magsize, df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(df_summary2$magsize, df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)



library(ggplot2)
ggplot(df_summary2, aes(x = magsize, y = mean_value)) +
  geom_point(color='skyblue')+
  geom_line(color='skyblue')+
  geom_errorbar(aes(ymin = mean_value - sd_value,
                    ymax = mean_value + sd_value),
                width = 0.1,
                color='skyblue'
  )+
  labs(x='MAG size (Mb)',y='DSs density')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "right",  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )

p <- ggplot(df_summary, aes(x = magsize, y = mean_value)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.1) +
  labs(x = 'MAG size (Mb)', y = 'DSs density (per Mb)') +
  #scale_color_manual(values = c("SA" = "#E9AF89", "SS" = "#9CB9D7", "WA" = "#B3DFC4", "WS" = "#DFA4BA")) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p
write.csv(df_summary,"Fig1G_2.csv")
ggsave("DSs_density_magsize.pdf", plot = p, width = 5.4, height = 3.42, units = "in")

#Fig1G_3
#统计各个大小MAG的数量并绘制直方图
#
setwd('D:/wenjian/research/defensome/groundwater/03mag_info')
size_info <- read.delim2('magsize_fenlei.txt')


# 设定新的水平顺序
new_levels <- size_info$MAG_size.Mb.
size_info$MAG_size.Mb. <- factor(size_info$MAG_size.Mb. , levels = new_levels)

library(ggplot2)
p <- ggplot(size_info,aes(x=MAG_size.Mb.)) +
  geom_bar(aes(x = order(MAG_size.Mb.),y=Number), color = "white",
           stat = "identity", position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,
           fill='skyblue'
  )+
  scale_x_discrete(limits = size_info$MAG_size.Mb.)+
  scale_y_continuous(limits = c(0,10000),expand = c(0,0))+
  labs(x = 'MAG size (Mb)', y = 'Number') +
  #scale_color_manual(values = c("SA" = "#E9AF89", "SS" = "#9CB9D7", "WA" = "#B3DFC4", "WS" = "#DFA4BA")) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p
write.csv(size_info,"Fig1G_1.csv")
ggsave("DSs_magsize_number2.pdf", plot = p, width = 5.4, height = 3.42, units = "in")





