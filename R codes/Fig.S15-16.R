# 安装并加载必要的包
# install.packages("ggplot2")
# install.packages("gggenes")
# install.packages("circlize")
library(ggplot2)
library(gggenes)
library(circlize)

#单个contig上####
#1. 读取数据
setwd("D:/wenjian/research/defensome/groundwater/31gene_circle")
virus_contig <- read.delim2("virus_3anti.gff",sep = '\t',header = F)
virus_contig <- virus_contig[-1,]
virus_contig$orientation <- ifelse(virus_contig$V7 == "-", 0, 1)

virus_contig$gene <- 'others'
virus_contig$gene[c(17,18,21,25)] <- c('Anti_CRISPR','Anti_CRISPR','Anti_CRISPR','Anti_Thoeris')

#2. 可视化(线性)
virus_3anti <- ggplot(virus_contig, aes(xmin = V4, xmax = V5, y = V3, fill = gene,
                         forward = orientation)) +
  scale_fill_manual(c(others="grey","#B6B8D8"))+
  geom_gene_arrow() +
  facet_wrap(~ V3, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")+
  theme_genes()
virus_3anti
setwd("D:/wenjian/research/defensome/groundwater/31gene_circle")
ggsave("virus_3anti.pdf",virus_3anti,width = 12,height = 1.2,units = 'in')

#prophage的contig与其对应的MAG的一起，进行比较####
#1. 读取数据
setwd("D:/wenjian/research/defensome/groundwater/31gene_circle")
virus_contig <- read.delim2("Ji15_bin.49_virus_anti.txt",sep = '\t',header = F)
virus_contig <- virus_contig[-1,]
virus_contig$orientation <- ifelse(virus_contig$V7 == "-", 0, 1)
virus_contig$gene <- 'others'
virus_contig$gene[c(21,23,73)] <- c('Anti_CRISPR','Anti_CRISPR','Anti_RM')
virus_contig$type <- 'prophage contig'

mag_contig <- read.delim2("Ji15_bin.49_defense.txt",sep = '\t',header = F)
mag_contig <- mag_contig[-1,]
mag_contig$orientation <- ifelse(mag_contig$V7 == "-", 0, 1)
mag_contig$gene <- 'others'
mag_contig$gene[c(2,4)] <- c('RM','RM')
mag_contig$type <- 'MAG contig'

all <- rbind(virus_contig,mag_contig)
#2. 可视化(线性)
p_mag_virus <- ggplot(all, aes(xmin = V4, xmax = V5, y = V3, fill = gene,
                                        forward = orientation)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_gene_arrow() +
  facet_wrap(~ type, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")+
  theme_genes()
p_mag_virus
setwd("D:/wenjian/research/defensome/groundwater/31gene_circle")
ggsave("p_mag_virus.pdf",p_mag_virus,width = 12,height = 1.2,units = 'in')


#十个prophage中含有大于等于3个anti-gene，对它们进行可视化#########
setwd("D:/wenjian/research/defensome/groundwater/31gene_circle")

prophage_anti_coloca <- read.delim2("coloca3_anti.gff",header = F)
prophage_anti_coloca$orientation <- ifelse(prophage_anti_coloca$V7 == "-", 0, 1)
library(dplyr)
library(stringr)
# 提取V9列中的ID部分
prophage_anti_coloca <- prophage_anti_coloca %>%
  mutate(ID_part = str_extract(V9, "ID=[^;]+")) %>%
  mutate(ID_suffix = str_extract(ID_part, "_\\d+")) %>%
  mutate(ORF = paste0(V1, ID_suffix)) %>%
  select(-ID_part, -ID_suffix)

gene_anno <- read.delim2("coloca3.txt",header = F)
colnames(gene_anno) <- c('sample','MAG','contig','ORF',"type",'subtype','gene')
prophage_anti_coloca <- merge(prophage_anti_coloca,gene_anno,by.x = 'ORF',by.y = 'ORF',all = T)

prophage_anti_coloca$type[is.na(prophage_anti_coloca$type)] <- 'Others'

#2. 可视化(线性)
library(ggplot2)
library(gggenes)
p_prophage_anti_colo <- ggplot(prophage_anti_coloca,
                               aes(xmin = V4, xmax = V5, 
                                   y = V1, fill = type,
                                 
                               forward = orientation)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_feature(
    data = na.omit(prophage_anti_coloca),
    aes(x = V4, y = V1, forward = T)
  ) +
  geom_feature_label(
    data = na.omit(prophage_anti_coloca),
    aes(x = V4, y = V1, label = gene, forward = T)
  ) +
  geom_gene_arrow() +
  facet_wrap(~ V1, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")+
  theme_genes()+
  theme(axis.title.y = element_blank(),# 隐藏y轴标题
        axis.text.y = element_blank()   # 隐藏y轴刻度标签
        )  

p_prophage_anti_colo
setwd("D:/wenjian/research/defensome/groundwater/31gene_circle")
ggsave("p_prophage_anti_colo2.pdf",p_prophage_anti_colo,width = 12,height = 12,units = 'in')

#挑选特定的contig进行可视化
prophage_anti_coloca_4 <- prophage_anti_coloca[prophage_anti_coloca$V1 %in% unique(prophage_anti_coloca$V1)[c(1,2,10)],]
prophage_anti_coloca_4 <- prophage_anti_coloca_4[prophage_anti_coloca_4$V4<20000,]

p_prophage_anti_colo4 <- ggplot(prophage_anti_coloca_4,
                               aes(xmin = V4, xmax = V5, 
                                   y = V1, fill = type,
                                   
                                   forward = orientation)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_feature(
    data = na.omit(prophage_anti_coloca_4),
    aes(x = V4, y = V1, forward = T)
  ) +
  geom_feature_label(
    data = na.omit(prophage_anti_coloca_4),
    aes(x = V4, y = V1, label = gene, forward = T)
  ) +
  geom_gene_arrow() +
  facet_wrap(~ V1, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")+  #Others的颜色是"#FDB462"
  # scale_fill_manual(values = c('grey','grey','grey','grey',
  #                              'grey','grey','grey'))+
  theme_genes()+
  theme(axis.title.y = element_blank(),# 隐藏y轴标题
        axis.text.y = element_blank()   # 隐藏y轴刻度标签
  )+
  labs(y = "",x=NULL)+
  guides(fill=guide_legend(title = "Genes",nrow =1,position = 'bottom'))

p_prophage_anti_colo4
setwd("D:/wenjian/research/defensome/groundwater/31gene_circle")
write.csv(prophage_anti_coloca_4,"Fig6A.csv",quote = F)
ggsave("p_prophage_anti_colo4.pdf",p_prophage_anti_colo4,width = 6.8,height = 4.3,units = 'in')


setwd("D:/wenjian/research/defensome/groundwater/31gene_circle/eggnog_10prophage_anticolocation")
eggnog <- read.delim("out.emapper.annotations",check.names = F)
eggnog <- select(eggnog,query,Description,Preferred_name,PFAMs)
prophage_anti_coloca_4 <- merge(prophage_anti_coloca_4,eggnog,
                                by.x = 'ORF',by.y = 'query',all = T)
prophage_anti_coloca_4 <- subset(prophage_anti_coloca_4, !is.na(type))


#绘制virus上defense与anti-defense共定位####
setwd("D:/wenjian/research/defensome/groundwater/31gene_circle/virus_ds_anti_both")

virus_65contig <- read.delim("65virus.gff",check.names = F,header = F)
virus_65contig$orientation <- ifelse(virus_65contig$V7 == "-", 0, 1)
# 提取V9列中的ID部分
virus_65contig <- virus_65contig %>%
  mutate(ID_part = str_extract(V9, "ID=[^;]+")) %>%
  mutate(ID_suffix = str_extract(ID_part, "_\\d+")) %>%
  mutate(ORF = paste0(V1, ID_suffix)) %>%
  select(-ID_part, -ID_suffix)

anti_gene <- read.delim("65virus_anti_annotation.txt",header = F,check.names = F)
colnames(anti_gene) <- c('sample','MAG','contig','ORF',"type",'subtype','gene')
defense_gene <- read.delim("65virus_DG_annotation.txt",header = F,check.names = F)
defense_gene <- defense_gene[,c(1,2,3)]
colnames(defense_gene) <- c('contig','ORF',"gene")
defense_gene$type <- 'Defense gene'

virus_65contig <- merge(virus_65contig,anti_gene,by = 'ORF',all = T)
virus_65contig <- merge(virus_65contig,defense_gene,by = 'ORF',all = T)



# 合并gene列
library(tidyr)
final_df <- unite(virus_65contig, "gene", gene.x, gene.y, sep = ";", na.rm = TRUE)
# 合并type列
final_df <- unite(final_df,"type",type.x,type.y,sep = ';',na.rm = TRUE)

final_df$type[final_df$type == ""] <- "Others"
final_df$gene[final_df$gene == ''] <- NA

virus_3contig <- final_df[final_df$V1 %in% unique(final_df$V1)[c(27,31,21)],]
virus_3contig <- virus_3contig[virus_3contig$V4<100000,]

p_virus_3contig <- ggplot(virus_3contig,
                                aes(xmin = V4, xmax = V5, 
                                    y = V1, fill = type,
                                    
                                    forward = orientation)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_feature(
    data = na.omit(virus_3contig[,c(-12:-14,-16,-18)]),
    aes(x = V4, y = V1, forward = T)
  ) +
  geom_feature_label(
    data = na.omit(virus_3contig[,c(-12:-14,-16,-18)]),
    aes(x = V4, y = V1, label = gene, forward = T)
  ) +
  geom_gene_arrow() +
  facet_wrap(~ V1, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")+  #Others的颜色是"#FDB462"
  # scale_fill_manual(values = c('grey','grey','grey','grey',
  #                              'grey','grey','grey'))+
  theme_genes()+
  theme(axis.title.y = element_blank(),# 隐藏y轴标题
        axis.text.y = element_blank()   # 隐藏y轴刻度标签
  )+
  labs(y = "",x=NULL)+
  guides(fill=guide_legend(title = "Genes",nrow =1,position = 'bottom'))
p_virus_3contig

setwd("D:/wenjian/research/defensome/groundwater/31gene_circle/virus_ds_anti_both")
ggsave("p_virus_3contig.pdf",p_virus_3contig,width = 6.8,height = 4.3,units = 'in')

#绘制acr与cas互作的两个contig######
setwd("D:/wenjian/research/defensome/groundwater/31gene_circle/acr_cas_two_contig")

acrcascontig <- read.delim("two_special_contig.gff",check.names = F,header = F)
acrcascontig$orientation <- ifelse(acrcascontig$V7 == "-", 0, 1)
# 提取V9列中的ID部分
acrcascontig <- acrcascontig %>%
  mutate(ID_part = str_extract(V9, "ID=[^;]+")) %>%
  mutate(ID_suffix = str_extract(ID_part, "_\\d+")) %>%
  mutate(ORF = paste0(V1, ID_suffix)) %>%
  select(-ID_part, -ID_suffix)

anti_gene <- read.delim("twocontig_anti_annotation.txt",header = F,check.names = F)
anti_gene <- anti_gene[,c(1,2,3)]
colnames(anti_gene) <- c('contig','ORF',"gene")
anti_gene$type <- 'Anti-defense gene'
#colnames(anti_gene) <- c('sample','MAG','contig','ORF',"type",'subtype','gene')
defense_gene <- read.delim("twocontig_DG_annotation.txt",header = F,check.names = F)
defense_gene <- defense_gene[,c(1,2,3)]
colnames(defense_gene) <- c('contig','ORF',"gene")
defense_gene$type <- 'Defense gene'

acrcascontig <- merge(acrcascontig,anti_gene,by = 'ORF',all = T)
acrcascontig <- merge(acrcascontig,defense_gene,by = 'ORF',all = T)



# 合并gene列
library(tidyr)
final_df <- unite(acrcascontig, "gene", gene.x, gene.y, sep = ";", na.rm = TRUE)
# 合并type列
final_df <- unite(final_df,"type",type.x,type.y,sep = ';',na.rm = TRUE)

final_df$type[final_df$type == ""] <- "Others"
final_df$gene[final_df$gene == ''] <- NA

# virus_3contig <- final_df[final_df$V1 %in% unique(final_df$V1)[c(27,31,21)],]
# virus_3contig <- virus_3contig[virus_3contig$V4<100000,]

#添加spacer和protospacer
protospacer <- final_df[50,]
protospacer$V4 <- 14571
protospacer$V5 <- 14600
protospacer$gene <- 'Protospacer'
protospacer$type <- 'Protospacer'

spacer <- final_df[1,]
spacer$V4 <- 126
spacer$V5 <- 689
spacer$gene <- 'CRISPR'
spacer$type <- 'CRISPR'

final_df <- rbind(final_df,spacer,protospacer)
# 假设数据框名为 final_df
final_df$type <- factor(final_df$type, levels = c("Anti-defense gene",
                                                  "CRISPR", "Defense gene",
                                                  "Protospacer","Others"
                                                  ))

p_acrcascontig <- ggplot(final_df,
                          aes(xmin = V4, xmax = V5, 
                              y = V1, fill = type,
                              
                              forward = orientation)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_feature(
    data = final_df[!is.na(final_df$gene),],
    aes(x = V4, y = V1, forward = T)
  ) +
  geom_feature_label(
    data = final_df,
    aes(x = V4, y = V1, label = gene, forward = T)
  ) +
  geom_gene_arrow() +  #data = 
  facet_wrap(~ V1, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")+  #Others的颜色是"#FDB462"
  # scale_fill_manual(values = c('grey','grey','grey','grey',
  #                              'grey','grey','grey'))+
  theme_genes()+
  theme(axis.title.y = element_blank(),# 隐藏y轴标题
        axis.text.y = element_blank()   # 隐藏y轴刻度标签
  )+
  labs(y = "",x=NULL)+
  guides(fill=guide_legend(title = "Genes",nrow =1,position = 'bottom'))
p_acrcascontig

setwd("D:/wenjian/research/defensome/groundwater/31gene_circle/acr_cas_two_contig")
write.csv(final_df,"Fig6C.csv",quote = F)
ggsave("p_acrcascontig.pdf",p_acrcascontig,width = 6.8,height = 3,units = 'in')

#绘制protospacer的结构,9repeats+8spacers

# 初始化参数
start_pos <- 126  # CRISPR的起始位置
repeat_length <- 36  # repeat的长度
spacer_length <- 30  # spacer的长度
num_repeats <- 9  # 重复的次数

# 创建一个空的数据框来存储结果
crispr_structure <- data.frame(
  type = character(),  # 类型：repeat 或 spacer
  start = integer(),   # 起始位置
  end = integer(),     # 终止位置
  stringsAsFactors = FALSE
)

# 生成repeat和spacer的起始和终止位置
for (i in 1:num_repeats) {
  # 添加repeat
  repeat_start <- start_pos
  repeat_end <- repeat_start + repeat_length - 1
  crispr_structure <- rbind(crispr_structure, data.frame(type = "repeat", start = repeat_start, end = repeat_end))
  
  # 更新start_pos以定位下一个spacer的起始位置
  start_pos <- repeat_end + 1
  
  # 如果不是最后一个repeat，添加spacer
  if (i < num_repeats) {
    spacer_start <- start_pos
    spacer_end <- spacer_start + spacer_length - 1
    crispr_structure <- rbind(crispr_structure, data.frame(type = "spacer", start = spacer_start, end = spacer_end))
    
    # 更新start_pos以定位下一个repeat的起始位置
    start_pos <- spacer_end + 1
  }
}

# 输出数据框
print(crispr_structure)
crispr_structure$type2 <- "CRISPR"
p_crispr <- ggplot(crispr_structure,
                         aes(xmin = start, xmax = end, y=type2,
                              fill = type,ymin = 0, ymax = 1
                             )) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  # geom_feature(
  #   data = final_df[!is.na(final_df$gene),],
  #   aes(x = V4, y = V1, forward = T)
  # ) +
  # geom_feature_label(
  #   data = final_df,
  #   aes(x = V4, y = V1, label = gene, forward = T)
  # ) +
  geom_rect()+
  #geom_gene_arrow() +  #data = 
  #facet_wrap(~ V1, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")+  #Others的颜色是"#FDB462"
  # scale_fill_manual(values = c('grey','grey','grey','grey',
  #                              'grey','grey','grey'))+
  theme_genes()+
  theme(axis.title.y = element_blank(),# 隐藏y轴标题
        axis.text.y = element_blank()   # 隐藏y轴刻度标签
  )+
  labs(y = "",x=NULL)+
  guides(fill=guide_legend(title = "",nrow =1,position = 'bottom'))
p_crispr

setwd("D:/wenjian/research/defensome/groundwater/31gene_circle/acr_cas_two_contig")
ggsave("p_crispr.pdf",p_crispr,width = 3.5,height = 1.5,units = 'in')


#绘制含有T4SS的防御岛######
setwd("D:/wenjian/research/defensome/groundwater/31gene_circle/DIs_with-T4SS")

acrcascontig <- read.delim("di119.gff",check.names = F,header = F)
acrcascontig$orientation <- ifelse(acrcascontig$V7 == "-", 0, 1)
# 提取V9列中的ID部分
acrcascontig <- acrcascontig %>%
  mutate(ID_part = str_extract(V9, "ID=[^;]+")) %>%
  mutate(ID_suffix = str_extract(ID_part, "_\\d+")) %>%
  mutate(ORF = paste0(V1, ID_suffix)) %>%
  select(-ID_part, -ID_suffix)

cdd_gene <- read.delim("DIs_nonDG_gene_CDD.txt",header = T,check.names = F)
# 通用模式提取 ">" 后面的字符串
cdd_gene$ORF <- sub(".*>\\s*([^ #]+).*", "\\1", cdd_gene$Query)

cdd_gene <- cdd_gene[,c("ORF","Short name")]
colnames(cdd_gene) <- c('ORF',"gene")
cdd_gene$type <- 'T4SS'

#colnames(anti_gene) <- c('sample','MAG','contig','ORF',"type",'subtype','gene')
defense_gene <- read.delim("di119_DG_annotation.txt",header = F,check.names = F)
defense_gene <- defense_gene[,c(1,2,3)]
colnames(defense_gene) <- c('contig','ORF',"gene")
defense_gene$type <- 'Defense gene'

acrcascontig <- merge(acrcascontig,cdd_gene,by = 'ORF',all = T)
acrcascontig <- merge(acrcascontig,defense_gene,by = 'ORF',all = T)



# 合并gene列
library(tidyr)
final_df <- unite(acrcascontig, "gene", gene.x, gene.y, sep = ";", na.rm = TRUE)
# 合并type列
final_df <- unite(final_df,"type",type.x,type.y,sep = ';',na.rm = TRUE)

final_df$type[final_df$type == ""] <- "Others"
final_df$gene[final_df$gene == ''] <- NA

final_df <- final_df[!is.na(final_df$V2), ]

# 假设数据框名为 final_df
final_df$type <- factor(final_df$type, levels = c("Defense gene",
                                                  "T4SS","Others"
))

custom_colors <- c("Defense gene" = "#7779AE", 
                   "T4SS" = "#E8C07C", 
                   "type3" = "green", 
                   "Others" = "#91C8ED")  # 包含 "others" 的自定义颜色
p_acrcascontig <- ggplot(final_df,
                         aes(xmin = V4, xmax = V5, 
                             y = V1, fill = type,
                             
                             forward = orientation)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_feature(
    data = final_df[!is.na(final_df$gene),],
    aes(x = V4, y = V1, forward = T)
  ) +
  geom_feature_label(
    data = final_df,
    aes(x = V4, y = V1, label = gene, forward = T)
  ) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  size = 0.1,
                  color = NA,
                  arrowhead_width = unit(1, "mm")) +  #data = 
  scale_fill_manual(values = custom_colors) +  # 使用自定义颜色
#  facet_wrap(~ V1, scales = "free", ncol = 1) +
 # scale_fill_brewer(palette = "Set3")+  #Others的颜色是"#FDB462"
  # scale_fill_manual(values = c('grey','grey','grey','grey',
  #                              'grey','grey','grey'))+
  theme_genes()+
  theme(axis.title.y = element_blank(),# 隐藏y轴标题
        axis.text.y = element_blank()   # 隐藏y轴刻度标签
  )+
  labs(y = "",x=NULL)+
  guides(fill=guide_legend(title = "Genes",nrow =1,position = 'bottom'))
p_acrcascontig

setwd("D:/wenjian/research/defensome/groundwater/31gene_circle/DIs_with-T4SS")
ggsave("p_DI-with-T4SS.pdf",p_acrcascontig,width = 6.8,height = 3,units = 'in')


#绘制含有TxSS的防御岛######
setwd("D:/wenjian/research/defensome/groundwater/31gene_circle/DIs_with-T4SS")

acrcascontig <- read.delim("allDIs_t.gff",check.names = F,header = F)
acrcascontig$orientation <- ifelse(acrcascontig$V7 == "-", 0, 1)
# 提取V9列中的ID部分
acrcascontig <- acrcascontig %>%
  mutate(ID_part = str_extract(V9, "ID=[^;]+")) %>%
  mutate(ID_suffix = str_extract(ID_part, "_\\d+")) %>%
  mutate(ORF = paste0(V1, ID_suffix)) %>%
  select(-ID_part, -ID_suffix)

cdd_gene <- read.delim("DIs_nonDG_gene_CDD_TxSS.txt",header = T,check.names = F)
# 通用模式提取 ">" 后面的字符串
cdd_gene$ORF <- sub(".*>\\s*([^ #]+).*", "\\1", cdd_gene$Query)

cdd_gene <- cdd_gene[,c("ORF","Short name","TxSS")]
colnames(cdd_gene) <- c('ORF',"gene","TxSS")
cdd_gene$type <- 'TxSS'

#colnames(anti_gene) <- c('sample','MAG','contig','ORF',"type",'subtype','gene')
defense_gene <- read.delim("allDIs_DG_annotation.txt",header = F,check.names = F)
defense_gene <- defense_gene[,c(1,2,3)]
colnames(defense_gene) <- c('contig','ORF',"gene")
defense_gene$type <- 'Defense gene'

acrcascontig <- merge(acrcascontig,cdd_gene,by = 'ORF',all = T)
acrcascontig <- merge(acrcascontig,defense_gene,by = 'ORF',all = T)



# 合并gene列
library(tidyr)
final_df <- unite(acrcascontig, "gene", gene.x, gene.y, sep = ";", na.rm = TRUE)
# 合并type列
final_df <- unite(final_df,"type",type.x,type.y,sep = ';',na.rm = TRUE)

final_df$type[final_df$type == ""] <- "Others"
final_df$gene[final_df$gene == ''] <- NA

#final_df <- final_df[!is.na(final_df$V2), ]

# 假设数据框名为 final_df
final_df$type <- factor(final_df$type, levels = c("Defense gene",
                                                  "TxSS","Others"
))

custom_colors <- c("Defense gene" = "#7779AE", 
                   "TxSS" = "#E8C07C", 
                   "type3" = "green", 
                   "Others" = "#91C8ED")  # 包含 "others" 的自定义颜色


# 提取最后一个下划线前的部分
cdd_contig <- sub("(_[^_]*)$", "", cdd_gene$ORF)

final_df <- final_df[final_df$V1 %in% cdd_contig,]

library(dplyr)
library(stringr)
# 提取 ORF 列末尾的数字并转换为数值
final_df2 <- final_df %>%
  mutate(orf_number = as.numeric(str_extract(ORF, "\\d+$"))) %>%
  arrange(V1, orf_number) %>%  # 按 V1 和末尾数字排序
  group_by(V1) %>%
  mutate(
    # 确定要保留的行（gene 不为 NA）
    keep = ifelse(!is.na(gene), TRUE, FALSE),
    # 找到要删除的两端行的索引
    to_remove = cumsum(keep) == 0 & !keep
  ) %>%
  filter(!(to_remove & !keep)) %>%  # 删除两端行
  ungroup() %>%
  select(-orf_number, -keep, -to_remove)  # 清理临时列


final_df2 <- final_df %>%
  mutate(orf_number = as.numeric(str_extract(ORF, "\\d+$"))) %>%
  arrange(V1, orf_number) %>%  # 按 V1 和末尾数字排序
  group_by(V1) %>%
  mutate(
    # 找到 gene 不为 NA 的行的索引
    is_gene_present = !is.na(gene),
    # 获取 gene 不为 NA 的行的索引
    gene_indices = which(is_gene_present),
    # 确定需要保留的行的范围
    to_keep = !(row_number() %in% c(gene_indices[1], gene_indices[length(gene_indices)]))
  ) %>%
  filter(to_keep | !is_gene_present) %>%  # 保留需要的行
  ungroup() %>%
  select(-orf_number, -is_gene_present, -gene_indices, -to_keep)  # 清理临时列


final_df2 <- final_df %>%
  mutate(orf_number = as.numeric(str_extract(ORF, "\\d+$"))) %>%
  arrange(V1, orf_number) %>%  # 按 V1 和末尾数字排序
  group_by(V1) %>%
  mutate(
    is_gene_present = !is.na(gene),
    gene_indices = list(which(is_gene_present)),  # 使用 list() 存储
    total_genes = length(gene_indices[[1]])  # 获取总基因数
  ) %>%
  ungroup() %>%
  mutate(
    to_keep = !((row_number() %in% c(gene_indices[[1]][1], gene_indices[[1]][total_genes])) | 
                  (row_number() %in% gene_indices[[1]][1:3]) | 
                  (row_number() %in% gene_indices[[1]][(total_genes - 2):total_genes]))
  ) %>%
  filter(to_keep | !is.na(gene)) %>%  # 保留需要的行
  select(-orf_number, -is_gene_present, -gene_indices, -total_genes, -to_keep)  # 清理临时列

final_df2 <- final_df %>%
  mutate(orf_number = as.numeric(str_extract(ORF, "\\d+$"))) %>%
  arrange(V1, orf_number) %>%  # 按 V1 和末尾数字排序
  group_by(V1) %>%
  mutate(
    # 找到 gene 不为 NA 的行的索引
    is_gene_present = !is.na(gene),
    # 获取 gene 不为 NA 的行的索引，并将其放入列表
    gene_indices = list(which(is_gene_present))
  ) %>%
  ungroup() %>%
  mutate(
    # 提取每组的基因索引
    first_index = sapply(gene_indices, function(x) x[1]) -2,
    last_index = sapply(gene_indices, function(x) x[length(x)] +2),
    # 确定需要保留的行的范围
    to_keep = !(row_number() %in% c(first_index, last_index))
  ) %>%
  filter(to_keep | !is.na(gene)) %>%  # 保留需要的行
  select(-orf_number, -is_gene_present, -gene_indices, -first_index, -last_index, -to_keep)  # 清理临时列




p_acrcascontig <- ggplot(final_df,
                         aes(xmin = V4, xmax = V5, 
                             y = V1, fill = type,
                             
                             forward = orientation)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_feature(
    data = final_df[!is.na(final_df$gene),],
    aes(x = V4, y = V1, forward = T,fill = type)
  ) +
  geom_feature_label(
    data = final_df[!is.na(final_df$gene),],
    aes(x = V4, y = V1, label = gene, forward = T,color = type)
  ) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  size = 0.1,
                  color = NA,  #控制基因箭头的边框
                  arrowhead_width = unit(1, "mm")) +  #data = 
  scale_fill_manual(values = custom_colors) +  # 使用自定义颜色
  scale_color_manual(values = custom_colors) +  # 添加文字颜色的手动映射
  facet_wrap(~ V1, scales = "free", ncol = 1) +
  # scale_fill_brewer(palette = "Set3")+  #Others的颜色是"#FDB462"
  # scale_fill_manual(values = c('grey','grey','grey','grey',
  #                              'grey','grey','grey'))+
  theme_genes()+
  theme(axis.title.y = element_blank(),# 隐藏y轴标题
        axis.text.y = element_blank()   # 隐藏y轴刻度标签
  )+
  labs(y = "",x=NULL)+
  guides(fill=guide_legend(title = "Genes",nrow =1,position = 'bottom'))
p_acrcascontig

setwd("D:/wenjian/research/defensome/groundwater/31gene_circle/DIs_with-T4SS")
ggsave("p_allDI-with-TxSS.pdf",p_acrcascontig,width = 20,height = 50,units = 'in',
       limitsize = FALSE)








final_df2 <- final_df %>%
  mutate(orf_number = as.numeric(str_extract(ORF, "\\d+$"))) %>%
  arrange(V1, orf_number) %>%  # 按 V1 和末尾数字排序
  group_by(V1) %>%
  mutate(
    is_gene_present = !is.na(gene),
    # 获取 gene 不为 NA 的行的索引
    gene_indices = list(which(is_gene_present))
  ) %>%
  ungroup() %>%
  mutate(
    # 提取每组的第一个和最后一个基因索引
    first_index = sapply(gene_indices, function(x) x[1]),
    last_index = sapply(gene_indices, function(x) x[length(x)]),
    # 确定需要保留的行
    to_keep = row_number() %in% c(first_index ：last_index) | !is.na(gene)
  ) %>%
  filter(to_keep) %>%  # 仅保留需要的行
  select(-orf_number, -is_gene_present, -gene_indices, -first_index, -last_index, -to_keep)  # 清理临时列

final_df2 <- final_df %>%
  mutate(orf_number = as.numeric(str_extract(ORF, "\\d+$"))) %>%
  arrange(V1, orf_number) %>%  # 按 V1 和末尾数字排序
  group_by(V1) %>%
  mutate(
    is_gene_present = !is.na(gene),
    # 获取 gene 不为 NA 的行的索引
    gene_indices = list(which(is_gene_present))
  ) %>%
  mutate(
    # 提取每组的第一个和最后一个基因索引
    first_index = sapply(gene_indices, function(x) x[1]),
    last_index = sapply(gene_indices, function(x) x[length(x)]),
    # 标记要保留的行
    to_keep = (row_number() >= first_index & row_number() <= last_index) | !is.na(gene)
  ) %>%
  ungroup() %>%
  filter(to_keep) %>%  # 仅保留需要的行
  select(-orf_number, -is_gene_present, -gene_indices, -first_index, -last_index, -to_keep)  # 清理临时列


p_acrcascontig2 <- ggplot(final_df2,
                         aes(xmin = V4, xmax = V5, 
                             y = V1, fill = type,
                             
                             forward = orientation)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_feature(
    data = final_df2[!is.na(final_df2$gene),],
    aes(x = V4, y = V1, forward = T,fill = type)
  ) +
  geom_feature_label(
    data = final_df2[!is.na(final_df2$gene),],
    aes(x = V4, y = V1, label = gene, forward = T,color = type)
  ) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  size = 0.1,
                  color = NA,  #控制基因箭头的边框
                  arrowhead_width = unit(1, "mm")) +  #data = 
  scale_fill_manual(values = custom_colors) +  # 使用自定义颜色
  scale_color_manual(values = custom_colors) +  # 添加文字颜色的手动映射
  facet_wrap(~ V1, scales = "free", ncol = 1) +
  # scale_fill_brewer(palette = "Set3")+  #Others的颜色是"#FDB462"
  # scale_fill_manual(values = c('grey','grey','grey','grey',
  #                              'grey','grey','grey'))+
  theme_genes()+
  theme(axis.title.y = element_blank(),# 隐藏y轴标题
        axis.text.y = element_blank()   # 隐藏y轴刻度标签
  )+
  labs(y = "",x=NULL)+
  guides(fill=guide_legend(title = "Genes",nrow =1,position = 'bottom'))
p_acrcascontig2

setwd("D:/wenjian/research/defensome/groundwater/31gene_circle/DIs_with-T4SS")
ggsave("p_allDI-with-TxSS_2.pdf",p_acrcascontig2,width = 6,height = 8,units = 'in',
       limitsize = FALSE)








