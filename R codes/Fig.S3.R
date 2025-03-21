
#1 MAG number vs DGs number###############
setwd('D:/wenjian/research/defensome/groundwater/s10mag_dg')
dgs <- read.delim('dg_counts_matrix.csv',sep = ',',check.names = F)
rownames(dgs) <- dgs$genome

#读取所有27578个mags的数据
setwd('D:/wenjian/research/defensome/groundwater/02去冗余信息')
allmags <- read.delim2('Wdb_ANI99.csv',sep = ',',
                       quote = "",header = T
) 

#将不含dgs的mags也合并到一起
# 找出 allmags 中在 dgs 行名中不存在的 genome 值
not_in_dgs <- allmags$genome[!allmags$genome %in% rownames(dgs)]
# 使用不在 dss 行名中的 genome 值创建新的数据框
new_rows <- data.frame(
  genome = not_in_dgs,  # 重复 genome 值
  stringsAsFactors = FALSE
)
# 为 new_rows 添加与 dss 相同的列，并将这些列的值填充为 0
col_names <- setdiff(names(dgs), "genome")  # 获取 dss 中除了 genome 之外的列名
new_rows[col_names] <- lapply(col_names, function(col_name) rep(0, nrow(dgs)))
rownames(new_rows) <- new_rows$genome

#合并new_rows与dss
dgs <- rbind(dgs,new_rows)

#2. 计算每个MAGs的DSs数量
# 使用rowSums()对每行求和
dgs <- dgs[,-1]
row_sums <- rowSums(dgs)

# 将求和结果转换为向量
sum_vector <- as.vector(row_sums)
row_names <- rownames(dgs)

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
p <- ggplot(sum_counts_df[sum_counts_df$Frequency>120,],aes(x=Value,y=Frequency))+
  geom_bar(stat = "identity",position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,fill="#81D0D6"#DFA4BA #E9AF89  #B3DFC4       
  )+
  geom_text(aes(label = Frequency), angle = 90, hjust = -0.1,
            position = position_dodge(width = 0.9), 
            # vjust = -0.50, # 调整文本位置，避免与柱子重叠
            color = "black", # 文本颜色，根据柱子颜色调整
            size = 3.5)+ # 文本大小
  scale_y_continuous(limits = c(0,10000),expand = c(0,0))+
  labs( x="Number of DGs",y="Number of MAGs") +
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
p  
write.csv(sum_counts_df,"FigS3B.csv",quote = F)
# 导出PDF
setwd("D:/wenjian/research/defensome/groundwater/s10mag_dg")
ggsave("s1_DGs_distribution.pdf", plot = p, width = 10, height = 8, units = "in")

#2 #######################
#2.统计每种DGs在所有MAGs中出现的百分比
# 计算每一列非0数字所占的百分比并输出到新数据框
percent_nonzero <- colSums(dgs != 0) / nrow(dgs) * 100
new_df <- data.frame(Column = names(percent_nonzero), Percentage = percent_nonzero)
# 按照值的大小对数据集进行降序排列
dg_persentage <- new_df[order(-new_df$Percentage), ]

#3.绘制柱状图
library(ggplot2)
p2 <- ggplot(new_df[new_df$Percentage>=1,]) +
  geom_bar(aes(x = reorder(Column,-Percentage), y = Percentage),
           stat = "identity", position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,
           fill='#81D0D6'#E9AF89,sa#9CB9D7,ss  #DFA4BA,ws#E9AF89,sa
  )+
  scale_y_continuous(limits = c(0,30),expand = c(0,0))+
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
p2
write.csv(new_df,"FigS3A.csv",quote = F)
ggsave("s1_DGs_percentage.pdf", plot = p2, width = 10, height = 8, units = "in")

#可选theme
# theme(
#   panel.background = element_blank(),  # 设置背景为透明
#   panel.grid.major = element_blank(),  # 去掉主要网格线
#   panel.grid.minor = element_blank(),  # 去掉次要网格线
#   #panel.border = element_blank(),  # 去掉边界线
#   panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
#   axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
#   axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
#   axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
#   axis.text.x = element_text(angle = 15,vjust = 0.85,hjust = 0.75), #横轴旋转
#   axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
#   legend.position = "right",  # 设置图例位置在右侧
#   legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
#   legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
#   legend.background = element_blank(),  # 设置图例背景为透明
#   legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
#   legend.key = element_blank()  # 去掉图例键的背景
# )

# #4.绘制百分比柱状圈图
# # 确保ggplot2包已安装并加载
# if (!require(ggplot2)) install.packages("ggplot2")
# library(ggplot2)
# 绘制环状柱状图
ggplot(new_df[new_df$Percentage>=1,]) +
  geom_bar(aes(x = Column, y = Percentage),
           stat = "identity", position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,
           fill='skyblue'
  ) +
  coord_polar()+
  scale_y_continuous(limits = c(-50,100),expand = c(0,0)) +
  guides(size="none")+
  labs( x=NULL,y=NULL) +
  theme( panel.background = element_blank(),
         plot.background = element_blank(),
         #  axis.text = element_blank(),
         #  legend.text = element_markdown(size = 8, face = "bold"),
         legend.title = element_text(vjust = 0.5,hjust=0.5,face="bold"),
         legend.position = c(0.5,0.5),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.ticks = element_blank(),
         plot.margin = margin(0,0,0,0))

#5.计算DSs的密度####
setwd('D:/wenjian/research/defensome/groundwater/03mag_info')
magsize <- read.csv('genome_stats.csv')
# 在R中，可以使用sub函数去掉数据框某一列字符串中的.fa后缀
# 假设数据框为df，列名为filename，可以使用以下代码：
magdg <- data.frame(genome=rownames(dgs),dgs=rowSums(dgs))

dg_size <- merge(magdg,magsize,by.x = 'genome',by.y = 'genome')

dg_density <- data.frame(genome=dg_size$genome,
                         density=dg_size$dgs/dg_size$size*1000000,
                         number=dg_size$dgs,
                         length_Mb=dg_size$size/1000000
)
# 在R中，可以使用floor函数将某一列的数字取小数点前的整数，并形成新的一列
# 假设数据框为df，列名为original_column，可以使用以下代码：
dg_density$magsize <- floor(dg_density$length_Mb)

# 在R中，可以使用dplyr库计算数据框中各个分组对应数据的平均值和标准偏差
# 假设有数据集df，包含group和value列，可以使用以下代码：
library(dplyr)
df_summary <- dg_density %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(density),
    sd_value = sd(density)
  )
df_summary2 <- dg_density %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(number),
    sd_value = sd(number)
  )

#计算两列数字的相关性及显著性
cor.test(df_summary$magsize, df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(df_summary2$magsize, df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)



library(ggplot2)

p3 <- ggplot(df_summary, aes(x = magsize, y = mean_value)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.1) +
  labs(x = 'MAG size (Mb)', y = 'DGs density (per Mb)') +
  annotate("Text",label='ρ = 0.59, P = 0.038',x=5,y=6.5,family='Times')+
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
p3
write.csv(df_summary,"FigS3C.csv")
ggsave("s1_p3_DGs_density_magsize.pdf", plot = p3, width = 5.4, height = 3.42, units = "in")

p4 <- ggplot(df_summary2, aes(x = magsize, y = mean_value)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.1) +
  labs(x = 'MAG size (Mb)', y = 'DGs number') +
  annotate("Text",label='ρ = 0.97, P < 2.2e-16',x=5,y=50,family='Times')+
  #scale_color_manual(values = c("SA" = "#E9AF89", "SS" = "#9CB9D7", "WA" = "#B3DFC4", "WS" = "#DFA4BA")) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),  # 横轴标签旋转
    text = element_text(family = "Times"),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p4
write.csv(df_summary2,"FigS3D.csv")
ggsave("s1_p4_DGs_number_magsize.pdf", plot = p4, width = 5.4, height = 3.42, units = "in")

library(cowplot)
g <- plot_grid(p3, p4,
               labels = LETTERS[2:3], nrow = 2
               #rel_heights = c(1, 1, 1.5),
              # label_x = 0.09,            # 标签x坐标（1为最右侧）
              # label_y = 0.9             # 标签y坐标（1为最上方）
)
g

setwd("D:/wenjian/research/defensome/groundwater/s1mag_dg")
ggsave("s1_DGnumber-density_magsize.pdf", plot = g, width = 5.4, height = 6.84, units = "in")



