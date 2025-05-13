#探究绘制多个特别关注的微生物门中DS的分布特征




#1. 读取通用数据

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

# 计算每个MAGs的DSs数量
# 使用rowSums()对每行求和
row_sums <- rowSums(dss)

# 将求和结果转换为向量
sum_vector <- as.vector(row_sums)
row_names <- rownames(dss)

# 创建一个新的数据框，包含求和结果和行名
num_dss <- data.frame(Sum = sum_vector, row.names = row_names)

#读取物种注释数据
setwd('D:/wenjian/research/defensome/groundwater/02去冗余信息')

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

#合并物种数据与DS number数据
num_dss$genome <- rownames(num_dss)

num_dss_tax <- merge(num_dss,mags_tax,by.x = "genome",by.y = "user_genome")

# 假设new_df是你的数据框，其中包含一个名为"Sum"的列

# 首先，确保"Sum"列的值是数值类型
num_dss_tax$Sum <- as.numeric(num_dss$Sum)

#2.  分门绘图#####
#p__Chlamydiota（衣原体）；p__Campylobacterota(厌氧自氧);p__Myxococcota(黏细菌); 
#p__Pseudomonadota(假单胞);p__Nitrospinota ; p__Nitrospirota ;p__Cyanobacteria

# 使用table()函数统计Sum列每个数值的出现次数
sum_counts <- table(num_dss_tax[num_dss_tax$col2=="p__Myxococcota",]$Sum)
# 将统计结果转换为数据框
sum_counts_df <- as.data.frame(sum_counts)
# 重命名列
colnames(sum_counts_df) <- c("Value", "Frequency")

# 确保Value列是因子类型
sum_counts_df$Value <- as.factor(sum_counts_df$Value)

sum_counts_df$Frequency2 <- sum_counts_df$Frequency/sum(sum_counts_df$Frequency)*100


#3.绘制柱状图
library(ggplot2)


# 计算参考线的位置
ref_line_value <- sum(num_dss_tax[num_dss_tax$col2=="p__Myxococcota",]$Sum)/sum(num_dss_tax$col2=="p__Myxococcota")

genome_num <- sum(num_dss_tax$col2=="p__Myxococcota")
genome_annotation <- paste("p__Myxococcota\n",genome_num,"genomes")

# 使用geom_bar()绘制柱状图，指定stat为"identity"
p_Myxococcota <- ggplot(sum_counts_df,aes(x=Value,y=Frequency2))+
  geom_bar(stat = "identity",position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,fill="#3d5a80"#DFA4BA #E9AF89  #B3DFC4       
  )+
  # geom_text(aes(label = Frequency), angle = 0, hjust = 0.5,vjust = -0.5,
  #           position = position_dodge(width = 0.9), 
  #           # vjust = -0.50, # 调整文本位置，避免与柱子重叠
  #           color = "black", # 文本颜色，根据柱子颜色调整
  #           size = 4)+ # 文本大小
  annotate("text", x = Inf, y = Inf, label = genome_annotation, hjust = 1.1,
           vjust = 1.1, size = 5, color = "black", parse = FALSE)+
  scale_y_continuous(expand = c(0,0))+
 # scale_y_continuous(limits = c(0, max(sum_counts_df$Frequency2) + 2))+ # 设置y轴的范围，留出顶部空间
  labs( x="Number of DSs",y="Percent genomes in phylum") +
  # geom_text(x = ref_line_value, y = max(sum_counts_df$Frequency), label = "Myxococcota",
  #           hjust = 1, vjust = 1, color = "black", size = 5) +  # 添加自定义文本
  geom_vline(xintercept = ref_line_value+1, linetype = "dashed", color = "red", linewidth = 0.5) +  # 添加参考线
  geom_vline(xintercept = 3.2+1, linetype = "dashed", color = "grey", size = 0.5) +
  theme(axis.text.x=element_text(angle = 90,vjust=0.5,hjust=0.5,color="black",size=10),
        axis.text.y=element_text(color="black",size = 10),
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
p_Myxococcota  
setwd('D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense')
write.csv(sum_counts_df,"figs11e.csv")
# 导出PDF
#setwd("D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense")
#ggsave("11_DSs_distribution_p__Myxococcota.pdf", plot = p, width = 3.86, height = 3.46, units = "in")

#################################p__Chlamydiota###################################################
# 使用table()函数统计Sum列每个数值的出现次数
sum_counts <- table(num_dss_tax[num_dss_tax$col2=="p__Chlamydiota",]$Sum)
# 将统计结果转换为数据框
sum_counts_df <- as.data.frame(sum_counts)
# 重命名列
colnames(sum_counts_df) <- c("Value", "Frequency")

# 确保Value列是因子类型
sum_counts_df$Value <- as.factor(sum_counts_df$Value)

sum_counts_df$Frequency2 <- sum_counts_df$Frequency/sum(sum_counts_df$Frequency)*100


#3.绘制柱状图
library(ggplot2)


# 计算参考线的位置
ref_line_value <- sum(num_dss_tax[num_dss_tax$col2=="p__Chlamydiota",]$Sum)/sum(num_dss_tax$col2=="p__Chlamydiota")

genome_num <- sum(num_dss_tax$col2=="p__Chlamydiota")
genome_annotation <- paste("p__Chlamydiota\n",genome_num,"genomes")

# 使用geom_bar()绘制柱状图，指定stat为"identity"
p_Chlamydiota <- ggplot(sum_counts_df,aes(x=Value,y=Frequency2))+
  geom_bar(stat = "identity",position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,fill="#98c1d9"#DFA4BA #E9AF89  #B3DFC4       
  )+
  # geom_text(aes(label = Frequency), angle = 0, hjust = 0.5,vjust = -0.5,
  #           position = position_dodge(width = 0.9), 
  #           # vjust = -0.50, # 调整文本位置，避免与柱子重叠
  #           color = "black", # 文本颜色，根据柱子颜色调整
  #           size = 4)+ # 文本大小
  annotate("text", x = Inf, y = Inf, label = genome_annotation, hjust = 1.1,
           vjust = 1.1, size = 5, color = "black", parse = FALSE)+
  scale_y_continuous(expand = c(0,0))+
  # scale_y_continuous(limits = c(0, max(sum_counts_df$Frequency2) + 2))+ # 设置y轴的范围，留出顶部空间
  labs( x="Number of DSs",y="Percent genomes in phylum") +
  # geom_text(x = ref_line_value, y = max(sum_counts_df$Frequency), label = "Myxococcota",
  #           hjust = 1, vjust = 1, color = "black", size = 5) +  # 添加自定义文本
  geom_vline(xintercept = ref_line_value+1, linetype = "dashed", color = "red", linewidth = 0.5) +  # 添加参考线
  geom_vline(xintercept = 3.2+1, linetype = "dashed", color = "grey", size = 0.5) +
  theme(axis.text.x=element_text(angle = 90,vjust=0.5,hjust=0.5,color="black",size=10),
        axis.text.y=element_text(color="black",size = 10),
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
p_Chlamydiota
write.csv(sum_counts_df,"figs11f.csv")
##############################p__Patescibacteria####################
# 使用table()函数统计Sum列每个数值的出现次数
sum_counts <- table(num_dss_tax[num_dss_tax$col2=="p__Patescibacteria",]$Sum)
# 将统计结果转换为数据框
sum_counts_df <- as.data.frame(sum_counts)
# 重命名列
colnames(sum_counts_df) <- c("Value", "Frequency")

# 确保Value列是因子类型
sum_counts_df$Value <- as.factor(sum_counts_df$Value)

sum_counts_df$Frequency2 <- sum_counts_df$Frequency/sum(sum_counts_df$Frequency)*100


#3.绘制柱状图
library(ggplot2)


# 计算参考线的位置
ref_line_value <- sum(num_dss_tax[num_dss_tax$col2=="p__Patescibacteria",]$Sum)/sum(num_dss_tax$col2=="p__Patescibacteria")

genome_num <- sum(num_dss_tax$col2=="p__Patescibacteria")
genome_annotation <- paste("p__Patescibacteria\n",genome_num,"genomes")

# 使用geom_bar()绘制柱状图，指定stat为"identity"
p_Patescibacteria <- ggplot(sum_counts_df,aes(x=Value,y=Frequency2))+
  geom_bar(stat = "identity",position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,fill="#e09f3e"#DFA4BA #E9AF89  #B3DFC4        
  )+
  # geom_text(aes(label = Frequency), angle = 0, hjust = 0.5,vjust = -0.5,
  #           position = position_dodge(width = 0.9), 
  #           # vjust = -0.50, # 调整文本位置，避免与柱子重叠
  #           color = "black", # 文本颜色，根据柱子颜色调整
  #           size = 4)+ # 文本大小
  annotate("text", x = Inf, y = Inf, label = genome_annotation, hjust = 1.1,
           vjust = 1.1, size = 5, color = "black", parse = FALSE)+
  scale_y_continuous(expand = c(0,0))+
  # scale_y_continuous(limits = c(0, max(sum_counts_df$Frequency2) + 2))+ # 设置y轴的范围，留出顶部空间
  labs( x="Number of DSs",y="Percent genomes in phylum") +
  # geom_text(x = ref_line_value, y = max(sum_counts_df$Frequency), label = "Myxococcota",
  #           hjust = 1, vjust = 1, color = "black", size = 5) +  # 添加自定义文本
  geom_vline(xintercept = ref_line_value+1, linetype = "dashed", color = "red", linewidth = 0.5) +  # 添加参考线
  geom_vline(xintercept = 3.2+1, linetype = "dashed", color = "grey", size = 0.5) +
  theme(axis.text.x=element_text(angle = 90,vjust=0.5,hjust=0.5,color="black",size=10),
        axis.text.y=element_text(color="black",size = 10),
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
p_Patescibacteria
write.csv(sum_counts_df,"figs11b.csv")
###
#setwd("D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense")
#ggsave("11_DSs_distribution_p__Patescibacteria.pdf", plot = p, width = 5.8, height = 3.9, units = "in")

##################################p__Campylobacterota###########################
# 使用table()函数统计Sum列每个数值的出现次数
sum_counts <- table(num_dss_tax[num_dss_tax$col2=="p__Campylobacterota",]$Sum)
# 将统计结果转换为数据框
sum_counts_df <- as.data.frame(sum_counts)
# 重命名列
colnames(sum_counts_df) <- c("Value", "Frequency")

# 确保Value列是因子类型
sum_counts_df$Value <- as.factor(sum_counts_df$Value)

sum_counts_df$Frequency2 <- sum_counts_df$Frequency/sum(sum_counts_df$Frequency)*100


#3.绘制柱状图
library(ggplot2)


# 计算参考线的位置
ref_line_value <- sum(num_dss_tax[num_dss_tax$col2=="p__Campylobacterota",]$Sum)/sum(num_dss_tax$col2=="p__Campylobacterota")

genome_num <- sum(num_dss_tax$col2=="p__Campylobacterota")
genome_annotation <- paste("p__Campylobacterota\n",genome_num,"genomes")

# 使用geom_bar()绘制柱状图，指定stat为"identity"
p_Campylobacterota <- ggplot(sum_counts_df,aes(x=Value,y=Frequency2))+
  geom_bar(stat = "identity",position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,fill="#74c69d"#DFA4BA #E9AF89  #B3DFC4          
  )+
  # geom_text(aes(label = Frequency), angle = 0, hjust = 0.5,vjust = -0.5,
  #           position = position_dodge(width = 0.9), 
  #           # vjust = -0.50, # 调整文本位置，避免与柱子重叠
  #           color = "black", # 文本颜色，根据柱子颜色调整
  #           size = 4)+ # 文本大小
  annotate("text", x = Inf, y = Inf, label = genome_annotation, hjust = 1.1,
           vjust = 1.1, size = 5, color = "black", parse = FALSE)+
  scale_y_continuous(expand = c(0,0))+
  # scale_y_continuous(limits = c(0, max(sum_counts_df$Frequency2) + 2))+ # 设置y轴的范围，留出顶部空间
  labs( x="Number of DSs",y="Percent genomes in phylum") +
  # geom_text(x = ref_line_value, y = max(sum_counts_df$Frequency), label = "Myxococcota",
  #           hjust = 1, vjust = 1, color = "black", size = 5) +  # 添加自定义文本
  geom_vline(xintercept = ref_line_value+1, linetype = "dashed", color = "red", linewidth = 0.5) +  # 添加参考线
  geom_vline(xintercept = 3.2+1, linetype = "dashed", color = "grey", size = 0.5) +
  theme(axis.text.x=element_text(angle = 90,vjust=0.5,hjust=0.5,color="black",size=10),
        axis.text.y=element_text(color="black",size = 10),
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
p_Campylobacterota
write.csv(sum_counts_df,"figs11a.csv")
# 导出PDF
#setwd("D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense")
#ggsave("11_DSs_distribution_p__Campylobacterota.pdf", plot = p, width = 3.86, height = 3.46, units = "in")

#########p__Nitrospirota############
# 使用table()函数统计Sum列每个数值的出现次数
sum_counts <- table(num_dss_tax[num_dss_tax$col2=="p__Nitrospirota",]$Sum)
# 将统计结果转换为数据框
sum_counts_df <- as.data.frame(sum_counts)
# 重命名列
colnames(sum_counts_df) <- c("Value", "Frequency")

# 确保Value列是因子类型
sum_counts_df$Value <- as.factor(sum_counts_df$Value)

sum_counts_df$Frequency2 <- sum_counts_df$Frequency/sum(sum_counts_df$Frequency)*100


#3.绘制柱状图
library(ggplot2)


# 计算参考线的位置
ref_line_value <- sum(num_dss_tax[num_dss_tax$col2=="p__Nitrospirota",]$Sum)/sum(num_dss_tax$col2=="p__Nitrospirota")

genome_num <- sum(num_dss_tax$col2=="p__Nitrospirota")
genome_num
genome_annotation <- paste("p__Nitrospirota\n",genome_num,"genomes")

# 使用geom_bar()绘制柱状图，指定stat为"identity"
p_Nitrospirota <- ggplot(sum_counts_df,aes(x=Value,y=Frequency2))+
  geom_bar(stat = "identity",position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,fill="#f3e9d2"#DFA4BA #E9AF89  #B3DFC4          
  )+
  # geom_text(aes(label = Frequency), angle = 0, hjust = 0.5,vjust = -0.5,
  #           position = position_dodge(width = 0.9), 
  #           # vjust = -0.50, # 调整文本位置，避免与柱子重叠
  #           color = "black", # 文本颜色，根据柱子颜色调整
  #           size = 4)+ # 文本大小
  annotate("text", x = Inf, y = Inf, label = genome_annotation, hjust = 1.1,
           vjust = 1.1, size = 5, color = "black", parse = FALSE)+
  scale_y_continuous(expand = c(0,0))+
  # scale_y_continuous(limits = c(0, max(sum_counts_df$Frequency2) + 2))+ # 设置y轴的范围，留出顶部空间
  labs( x="Number of DSs",y="Percent genomes in phylum") +
  # geom_text(x = ref_line_value, y = max(sum_counts_df$Frequency), label = "Myxococcota",
  #           hjust = 1, vjust = 1, color = "black", size = 5) +  # 添加自定义文本
  geom_vline(xintercept = ref_line_value+1, linetype = "dashed", color = "red", linewidth = 0.5) +  # 添加参考线
  geom_vline(xintercept = 3.2+1, linetype = "dashed", color = "grey", size = 0.5) +
  theme(axis.text.x=element_text(angle = 90,vjust=0.5,hjust=0.5,color="black",size=10),
        axis.text.y=element_text(color="black",size = 10),
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
p_Nitrospirota
write.csv(sum_counts_df,"figs11c.csv")
#######################p__Nitrospinota#########################
# 使用table()函数统计Sum列每个数值的出现次数
sum_counts <- table(num_dss_tax[num_dss_tax$col2=="p__Nitrospinota",]$Sum)
# 将统计结果转换为数据框
sum_counts_df <- as.data.frame(sum_counts)
# 重命名列
colnames(sum_counts_df) <- c("Value", "Frequency")

# 确保Value列是因子类型
sum_counts_df$Value <- as.factor(sum_counts_df$Value)

sum_counts_df$Frequency2 <- sum_counts_df$Frequency/sum(sum_counts_df$Frequency)*100


#3.绘制柱状图
library(ggplot2)


# 计算参考线的位置
ref_line_value <- sum(num_dss_tax[num_dss_tax$col2=="p__Nitrospinota",]$Sum)/sum(num_dss_tax$col2=="p__Nitrospinota")

genome_num <- sum(num_dss_tax$col2=="p__Nitrospinota")
genome_num
genome_annotation <- paste("p__Nitrospinota\n",genome_num,"genomes")

# 使用geom_bar()绘制柱状图，指定stat为"identity"
p_Nitrospinota <- ggplot(sum_counts_df,aes(x=Value,y=Frequency2))+
  geom_bar(stat = "identity",position = "dodge",
           show.legend = T,alpha = .9,linewidth = 0.5,fill="#b8b8ff"#DFA4BA #E9AF89  #B3DFC4          
  )+
  # geom_text(aes(label = Frequency), angle = 0, hjust = 0.5,vjust = -0.5,
  #           position = position_dodge(width = 0.9), 
  #           # vjust = -0.50, # 调整文本位置，避免与柱子重叠
  #           color = "black", # 文本颜色，根据柱子颜色调整
  #           size = 4)+ # 文本大小
  annotate("text", x = Inf, y = Inf, label = genome_annotation, hjust = 1.1,
           vjust = 1.1, size = 5, color = "black", parse = FALSE)+
  scale_y_continuous(expand = c(0,0))+
  # scale_y_continuous(limits = c(0, max(sum_counts_df$Frequency2) + 2))+ # 设置y轴的范围，留出顶部空间
  labs( x="Number of DSs",y="Percent genomes in phylum") +
  # geom_text(x = ref_line_value, y = max(sum_counts_df$Frequency), label = "Myxococcota",
  #           hjust = 1, vjust = 1, color = "black", size = 5) +  # 添加自定义文本
  geom_vline(xintercept = ref_line_value+1, linetype = "dashed", color = "red", linewidth = 0.5) +  # 添加参考线
  geom_vline(xintercept = 3.2+1, linetype = "dashed", color = "grey", size = 0.5) +
  theme(axis.text.x=element_text(angle = 90,vjust=0.5,hjust=0.5,color="black",size=10),
        axis.text.y=element_text(color="black",size = 10),
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
p_Nitrospinota
write.csv(sum_counts_df,"figs11d.csv")

library(Cairo)
library(cowplot)

g_6phyla <- plot_grid(p_Campylobacterota,p_Patescibacteria, p_Nitrospirota,
                  p_Nitrospinota,p_Myxococcota,p_Chlamydiota,
                  labels = LETTERS[1:6], nrow = 2
                  #rel_heights = c(1, 1, 1.5),
                  # label_x = 0.09,            # 标签x坐标（1为最右侧）
                  #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g_6phyla

setwd("D:/wenjian/research/defensome/groundwater/02ANI99MAGs_defense")
cairo_pdf(width = 10.3,height = 5.5,filename = "11_ds_6phyla.pdf")
g_6phyla
dev.off()



