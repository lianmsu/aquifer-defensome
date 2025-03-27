#不含MGE的MAG，8763个 与  不含DSs的7731个MAG之间的重叠？没有一个重叠MAG

#加载库
library(dplyr)
library(reshape2)
#1. 读取prophage, plasmid，基本的统计分析####
#1.1 prophage和plasmid的读取
setwd("D:/pku/01groundwater/01defensome/NC_revision/1additional_analyses/ICEs")

#prophage <- read.delim2('virus_taxonomy_count.tsv')
#rownames(prophage) <- prophage[,1]
#sum(prophage[,-1]) #40302个
#prophage_sum <- rowSums(prophage[,-1]) %>% as.data.frame() %>% setNames('prophagesum')

all_phage <- read.delim("all_prophage_contig.tsv")
all_phage$genome <- gsub("\\.fasta.*","",all_phage$seq_name)
all_phage$prophage <- ifelse(all_phage$topology == "Provirus","Prophage","Phage")
library(dplyr)
library(tidyr)
# 根据 genome 列分组，并统计每组中 prophage 列的数量
virus_sum <- all_phage %>%
  group_by(genome,prophage) %>%
  summarize(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = prophage, values_from = count, values_fill = list(count = 0))

plasmid <- read.delim2('plasmid_taxonomy_count.tsv')
rownames(plasmid) <- plasmid[,1]
sum(plasmid[,-1]) #163897
plasmid_sum <- rowSums(plasmid[,-1]) %>% as.data.frame() %>% setNames("plasmidsum")

#intersection <- intersect(prophage$Prefix, plasmid$Prefix)# 取交集

#integrons
all_integrons <- read.delim("complete_integron_contig.id",header = F)
colnames(all_integrons) <- "contig_id"
all_integrons$genome <- gsub("\\.fasta.*","",all_integrons$contig_id)
library(stringr)
all_integrons$length <- str_extract(all_integrons$contig_id,"(?<=length_)\\d+")
all_integrons$length <- as.numeric(all_integrons$length)
integron_sum <- all_integrons %>%
  group_by(genome) %>%
  summarize(integron_num = n())

#ICE/IME
all_ICE <- read.delim("all_ICE_contigid.txt",header = F)
colnames(all_ICE) <- "contig_id"
all_ICE$genome <- gsub("\\_contig.*","",all_ICE$contig_id)
all_ICE$genome <- gsub(">","",all_ICE$genome)
# 使用正则表达式提取 start 和 end 位置
all_ICE$start <- as.numeric(gsub(".*_(\\d+)-(\\d+)$", "\\1", all_ICE$contig_id))
all_ICE$end <- as.numeric(gsub(".*_(\\d+)-(\\d+)$", "\\2", all_ICE$contig_id))
# 计算长度并存储到新列 length 中
all_ICE$length <- all_ICE$end - all_ICE$start + 1
ICE_sum <- all_ICE %>%
  group_by(genome) %>%
  summarize(ICE_num = n())


#1.2 MAGs物种信息的读取
setwd('D:/wenjian/research/defensome/groundwater/02去冗余信息')
mags_tax <- read.delim2('gtdbtk.bac120.summary.tsv')
# 假设 df 是你的数据框，MAGs 是其中的一列
mags_tax$user_genome <- paste(mags_tax$user_genome, ".fasta", sep="")
#mags_tax <- mags_tax[mags_tax$user_genome %in% rownames(dss),]

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






# #统计既不含质粒，又不含prophage的MAG，为0
# allmags$genome <- sub('.fasta' , '',allmags$genome)
# no_plasmid <- allmags$genome[!allmags$genome %in% plasmid$Prefix]
# no_plasmid_prophage <- no_plasmid[!no_plasmid %in% prophage$Prefix]
# 
# no_plasmid_prophage_dss <- no_plasmid_prophage[no_plasmid_prophage %in% not_in_dss]

#2. 统计质粒和噬菌体在不同门中的数量分布####
#将物种信息合并到质粒与噬菌体的数量表中
plasmid_sum$genome <- rownames(plasmid_sum)

mags_tax$user_genome <- sub('.fasta','',mags_tax$user_genome)
plasmid_sum <- merge(plasmid_sum,mags_tax,by.x = 'genome',by.y = 'user_genome')

virus_sum <- merge(virus_sum,mags_tax,by.x = 'genome',by.y = 'user_genome')

integron_sum <- merge(integron_sum,mags_tax,by.x = 'genome',by.y = 'user_genome')
ICE_sum <- merge(ICE_sum,mags_tax,by.x = 'genome',by.y = 'user_genome')

#对门级别的物种注释进行修正
plasmid_sum$col2 <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", plasmid_sum$col2)
plasmid_sum$col2 <- gsub("^p__", "", plasmid_sum$col2)

virus_sum$col2 <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", virus_sum$col2)
virus_sum$col2 <- gsub("^p__", "", virus_sum$col2)

integron_sum$col2 <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", integron_sum$col2)
integron_sum$col2 <- gsub("^p__", "", integron_sum$col2)

ICE_sum$col2 <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", ICE_sum$col2)
ICE_sum$col2 <- gsub("^p__", "", ICE_sum$col2)

#绘制质粒、噬菌体来自门的分布，桑基图; 
plasmid_sum$plasmid <- 'Plasmid'
#prophage_sum$prophage <- 'Prophage'

plasmid_md <- select(plasmid_sum,plasmidsum,col2,plasmid)
# 使用dplyr的group_by和summarise函数进行操作
plasmid_md <- plasmid_md %>%
  group_by(col2) %>%
  summarise(plasmidsum = sum(plasmidsum, na.rm = TRUE))
plasmid_md$plasmid <- 'Plasmid'

phage_md <- select(virus_sum,Phage,col2)
prophage_md <- select(virus_sum,Prophage,col2)
# 使用dplyr的group_by和summarise函数进行操作
phage_md <- phage_md %>%
  group_by(col2) %>%
  summarise(phagesum = sum(Phage, na.rm = TRUE))
phage_md$phage <- 'Phage'

prophage_md <- prophage_md %>%
  group_by(col2) %>%
  summarise(prophagesum = sum(Prophage, na.rm = TRUE))
prophage_md$prophage <- 'Prophage'

integron_md <- select(integron_sum,integron_num,col2)
integron_md <- integron_md %>%
  group_by(col2) %>%
  summarize(integronsum = sum(integron_num,na.rm = T))
integron_md$integron <- "Integron"

ICE_md <- select(ICE_sum,ICE_num,col2)
ICE_md <- ICE_md %>%
  group_by(col2) %>%
  summarize(ICEsum = sum(ICE_num,na.rm = T))
ICE_md$ICE <- "ICEs/IMEs"

#####桑基图##############
# 筛选出数量小于500的行，并将它们合并为一行
plasmid_filtered <- plasmid_md %>%
  filter(plasmidsum < 500) %>%
  summarise(combined_plasum = sum(plasmidsum)) %>% setNames("plasmidsum")
plasmid_filtered$col2 <- 'Other phyla'
plasmid_filtered$plasmid <- 'Plasmid'
plasmid_filtered <- select(plasmid_filtered,col2,plasmidsum,plasmid)

plasmid_md <- plasmid_md[plasmid_md$plasmidsum>500,]
plasmid_md <- rbind(plasmid_md,plasmid_filtered)
plasmid_md <- select(plasmid_md,plasmid,col2,plasmidsum)
colnames(plasmid_md) <- c('source','target','value')
colnames(plasmid_md) <- c("subclass_label","predict","count")

# 筛选出数量小于500的行，并将它们合并为一行
phage_filtered <- phage_md %>%
  filter(phagesum < 500) %>%
  summarise(combined_prosum = sum(phagesum)) %>% setNames("phagesum")
phage_filtered$col2 <- 'Other phyla'
phage_filtered$phage <- 'Phage'
phage_filtered <- select(phage_filtered,col2,phagesum,phage)

phage_md <- phage_md[phage_md$phagesum>500,]
phage_md <- rbind(phage_md,phage_filtered)
phage_md <- select(phage_md,phage,col2,phagesum)
colnames(phage_md) <- c('source','target','value')
colnames(phage_md) <- c("subclass_label","predict","count")

# 筛选出数量小于45的行，并将它们合并为一行
prophage_filtered <- prophage_md %>%
  filter(prophagesum < 45) %>%
  summarise(combined_prosum = sum(prophagesum)) %>% setNames("prophagesum")
prophage_filtered$col2 <- 'Other phyla'
prophage_filtered$prophage <- 'Prophage'
prophage_filtered <- select(prophage_filtered,col2,prophagesum,prophage)

prophage_md <- prophage_md[prophage_md$prophagesum>=45,]
prophage_md <- rbind(prophage_md,prophage_filtered)
prophage_md <- select(prophage_md,prophage,col2,prophagesum)
colnames(prophage_md) <- c('source','target','value')
colnames(prophage_md) <- c("subclass_label","predict","count")

# Integron
integron_filtered <- integron_md %>%
  filter(integronsum < 90) %>%
  summarise(combined_prosum = sum(integronsum)) %>% setNames("integronsum")
integron_filtered$col2 <- 'Other phyla'
integron_filtered$integron <- 'Integron'
integron_filtered <- select(integron_filtered,col2,integronsum,integron)

integron_md <- integron_md[integron_md$integronsum>=90,]
integron_md <- rbind(integron_md,integron_filtered)
integron_md <- select(integron_md,integron,col2,integronsum)
colnames(integron_md) <- c('source','target','value')
colnames(integron_md) <- c("subclass_label","predict","count")

#ICE
ICE_filtered <- ICE_md %>%
  filter(ICEsum < 40) %>%
  summarise(combined_prosum = sum(ICEsum)) %>% setNames("ICEsum")
ICE_filtered$col2 <- 'Other phyla'
ICE_filtered$ICE <- 'ICEs/IMEs'
ICE_filtered <- select(ICE_filtered,col2,ICEsum,ICE)

ICE_md <- ICE_md[ICE_md$ICEsum>=40,]
ICE_md <- rbind(ICE_md,ICE_filtered)
ICE_md <- select(ICE_md,ICE,col2,ICEsum)
colnames(ICE_md) <- c('source','target','value')
colnames(ICE_md) <- c("subclass_label","predict","count")

library(tidyverse)
library(ggsankeyfier)
df_plasmid <- plasmid_md %>% 
  pivot_stages_longer(.,
                      stages_from = c("subclass_label","predict"),
                      values_from = "count")

df_sankey <- rbind(plasmid_md,phage_md,prophage_md,integron_md,ICE_md)
df_sankey <-  df_sankey %>% 
  pivot_stages_longer(.,
                      stages_from = c("subclass_label","predict"),
                      values_from = "count")

mge_phyla_sankey <- 
  ggplot(data= df_sankey,aes(x = stage, y = count, group = node,
                             edge_id = edge_id,connector = connector)) +
  geom_sankeyedge(aes(fill = node),
                  position = position_sankey(order = "ascending",
                                             v_space ="auto",
                                             width = 0.05),
                  show.legend = F)+
  geom_sankeynode(aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",
                                             v_space ="auto",
                                             width = 0.05),
                  show.legend = F) +
  scale_fill_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                               "Integron"="#e8a598", "ICEs/IMEs"="#9f7e69",
                               "Pseudomonadota" = "#FCDABA",
                               "Patescibacteria" = "#00BA38",
                               "Phage" = "#ED6764")) +  # 手动设置每个节点的颜色
  scale_color_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                                "Integron"="#e8a598", "ICEs/IMEs"="#9f7e69",
                                "Pseudomonadota" = "#FCDABA", 
                                "Patescibacteria" = "#00BA38",
                                "Phage" = "#ED6764")) + # 手动设置每个节点的颜色
  coord_cartesian(clip = "off")+
  geom_text(data=df_sankey %>% filter(connector=="from"),
            aes(label = node),
            stat = "sankeynode",
            position = position_sankey(width = 0.05,
                                       v_space ="auto",
                                       order="ascending",
                                       nudge_x=-0.05),
            hjust=1,cex=2.8,color="black")+
  geom_text(data=df_sankey %>% filter(connector=="to"),aes(label = node),
            stat = "sankeynode",
            position = position_sankey(width = 0.05,
                                       v_space ="auto",
                                       order="ascending",
                                       nudge_x=0.05),
            hjust=0,cex=2.8,color="black")+
  geom_text(data=df_sankey %>% filter(connector=="to"),aes(label = count),
            stat = "sankeynode",
            position = position_sankey(width = 0.05,
                                       v_space ="auto",
                                       order="ascending",
                                       nudge_x=0.5),
            hjust=0,cex=2.8,color="black")+
  theme(panel.background = element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0,3,0,0,unit = "cm"))
mge_phyla_sankey
setwd("D:/wenjian/research/defensome/groundwater/05MGE")
setwd("D:/pku/01groundwater/01defensome/NC_revision/1additional_analyses/ICEs")
write.csv(df_sankey,"fig5a.csv")
#ggsave("1_mge_phyla_sankey_new.pdf",mge_phyla_sankey,width = 5.56,height = 4.39,units = "in")


############################不合并低丰度MGE的门,计算每个门中plasmid-prophage比例####
#重新读取2.部分数据(150行之前)

plasmid_md$col2[1] <-"Unassigned"
prophage_md$col2[1] <- "Unassigned"
phage_md$col2[1] <- "Unassigned"
integron_md$col2[1] <- "Unassigned"

plasmid_md <- select(plasmid_md,plasmid,col2,plasmidsum)
prophage_md <- select(prophage_md,prophage,col2,prophagesum)
phage_md <- select(phage_md,phage,col2,phagesum)
integron_md <- select(integron_md,integron,col2,integronsum)
ICE_md <- select(ICE_md,ICE,col2,ICEsum)

colnames(plasmid_md) <- c('mge_plasmid','phyla','plasmidsum')
colnames(prophage_md) <- c('mge_prophage','phyla','prophagesum')
colnames(phage_md) <- c('mge_phage','phyla','phagesum')
colnames(integron_md) <- c('mge_integron','phyla','integronsum')
colnames(ICE_md) <- c('mge_ICE','phyla','ICEsum')

#plasmid_prophage <- merge(plasmid_md,prophage_md,by = "phyla", all = T)
all_phyla_mge <- plasmid_md %>%
  full_join(prophage_md, by = "phyla") %>%
  full_join(phage_md, by = "phyla") %>%
  full_join(integron_md, by = "phyla") %>%
  full_join(ICE_md, by = "phyla")

all_phyla_mge <- all_phyla_mge %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)))

all_phyla_mge$mge_plasmid <- 'Plasmid'
all_phyla_mge$mge_prophage <- 'Prophage'
all_phyla_mge$mge_phage <- 'Phage'
all_phyla_mge$mge_integron <- 'Integron'
all_phyla_mge$mge_ICE <- 'ICE/IME'

all_phyla_mge$mgesum <- all_phyla_mge$plasmidsum + all_phyla_mge$prophagesum+
  all_phyla_mge$phagesum+all_phyla_mge$integronsum+all_phyla_mge$ICEsum
  
#all_phyla_mge$plasmidbiprophage <- all_phyla_mge$plasmidsum / all_phyla_mge$prophagesum

all_phyla_mge <- all_phyla_mge %>% filter(., phagesum>500 | plasmidsum>500)
all_phyla_mge$plasmidbili <- all_phyla_mge$plasmidsum/all_phyla_mge$mgesum *100
all_phyla_mge$prophagebili <- all_phyla_mge$prophagesum/all_phyla_mge$mgesum *100
all_phyla_mge$phagebili <- all_phyla_mge$phagesum/all_phyla_mge$mgesum *100
all_phyla_mge$integronbili <- all_phyla_mge$integronsum/all_phyla_mge$mgesum *100
all_phyla_mge$ICEbili <- all_phyla_mge$ICEsum/all_phyla_mge$mgesum *100


df_plot <- select(all_phyla_mge,phyla,
                  mge_plasmid,mge_prophage,mge_phage,mge_integron,mge_ICE,
                  plasmidbili,prophagebili,phagebili,integronbili,ICEbili)
df_plot <- melt(df_plot)

# 假设 df_plot 是你的数据框
# 筛选 variable == "prophagebili" 的行
phagebili_rows <- df_plot[df_plot$variable == "phagebili", ]

# 按照 value 列对 phyla 列排序
ordered_phyla <- phagebili_rows$phyla[order(phagebili_rows$value)]

# 重新定义 phyla 列为因子，并按照 ordered_phyla 的顺序排序
df_plot$phyla <- factor(df_plot$phyla, levels = ordered_phyla)

# 打印结果查看
print(df_plot)


library(ggplot2)
colorsss <- c(
  "#CED2DD",  # 青色
  "#ED6764",  # 橙色
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

# 绘制百分比堆叠柱状图 
mge_plasmid_prophage_bili<-  
  ggplot(df_plot, aes(x = phyla, y = value, fill = variable)) +
  geom_bar(stat = "identity", width = 0.8) +  # 使用实际值绘制柱状图
  #scale_y_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0, 100, by = 20)) +  # 设置 y 轴的范围和刻度
  #scale_y_continuous(limits=c(0,100)) +  # 将y轴标签转换为百分比
  #scale_fill_manual(values = colorsss)+
  scale_fill_manual(values = c("plasmidbili" = "#619CFF", "prophagebili"="#fa9500",
                               "integronbili"="#e8a598", "ICEbili"="#9f7e69",
                               "phagebili" = "#ED6764")) +  # 手动设置每个节点的颜色
  # scale_x_discrete(limits=c('Cd11','Cr11','Cu11','phOH11','PNP11','PAP11','CK11',
  #                           'Cd21','Cr21','Cu21','phOH21','PNP21','PAP21','CK21',
  #                           'Seed'))+
  #coord_flip() +  # 翻转坐标轴
  #scale_y_continuous(limits=c(0,100),expand = c(0,0)) +  # 将y轴标签转换为百分比
  labs(x = "", y = "Relative proportion (%)", fill = "") +  # 设置轴标签和图例标题
  theme_minimal()+  # 使用简洁主题
  #  geom_vline(xintercept = c(7.5,14.5),linetype=2)+
  guides(fill=guide_legend(title = "",nrow =1))+
  theme(axis.text.x=element_text(angle = 90,vjust=0.5,hjust=0.5,color="black",size =10),
        axis.text.y=element_text(color="black",size = 10),
        plot.background = element_rect(fill="white"), 
        panel.background = element_rect(fill="white"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0,"cm"),
        legend.position = 'top',
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"),
        axis.line.x=element_line(color="black"),
        axis.line.y.left = element_line(color="grey30"),
        axis.line.y.right = element_line(color="grey30"),
        axis.line.x.bottom = element_line(color="grey30"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank())+
  guides(y = guide_axis(minor.ticks = TRUE))
mge_plasmid_prophage_bili

setwd("D:/pku/01groundwater/01defensome/NC_revision/1additional_analyses/ICEs")
write.csv(df_plot,"fig5b.csv")
#ggsave("1mge_plasmid_prophage_bili_new.pdf", mge_plasmid_prophage_bili, width =5 ,height =4 ,units = 'in' )
#############################不合并门


#3. 读取所有MAGs中的dss信息####
#读取数据表###
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


#4. MAGs中DSs的密度####
setwd('D:/wenjian/research/defensome/groundwater/03mag_info')
magsize <- read.csv('genome_stats.csv')
magsize$genome <- sub('.fasta','',magsize$genome)
# 在R中，可以使用sub函数去掉数据框某一列字符串中的.fa后缀
# 假设数据框为df，列名为filename，可以使用以下代码：
magds <- data.frame(genome=rownames(dss),dss=rowSums(dss))

ds_size <- merge(magds,magsize,by.x = 'genome',by.y = 'genome')
ds_size$genome <- sub('.fasta', '',ds_size$genome)
ds_density <- data.frame(genome=ds_size$genome,
                         density=ds_size$dss/ds_size$size*1000000,
                         number=ds_size$dss,
                         length_Mb=ds_size$size/1000000
)


# #5. prophage、plasmid数量与MAGs大小的关系####
# prophage_sum$genome <- rownames(prophage_sum)
# magsize$genome <- sub('.fasta', '',magsize$genome)
# prophage_sum_size <- merge(prophage_sum,magsize,
#                            by = 'genome')
# 
# #检验prophage数量与MAG size的关系
# rand<-rnorm(12427,sd=1e-10)
# prophage_sum_size$prophagesum <- prophage_sum_size$prophagesum+rand
# prophage_sum_size$size <- prophage_sum_size$size+rand
# cor.test(prophage_sum_size$prophagesum,
#          prophage_sum_size$size,
#          alternative = "two.sided",
#          method = "spearman",
#          conf.level = 0.95
# )#ρ0.157，p<2.2e-16
# 
# 
# plasmid_sum$genome <- rownames(plasmid_sum)
# plasmid_sum_size <- merge(plasmid_sum,magsize,
#                           by = 'genome'
# )
# 
# #检验plasmid数量与MAG size的关系
# cor.test(plasmid_sum_size$plasmidsum,
#          plasmid_sum_size$size,
#          alternative = "two.sided",
#          method = "spearman",conf.level = 0.95)#ρ0.268 p<2.2e-16
# plot(plasmid_sum_size$plasmidsum,
#      plasmid_sum_size$size)

#6. prophage、plasmid、染色体上DSs的密度####
#可能不是目前这种算法。应该是DGs/DSs数量除以（所在的prophage的长度，所在plasmid的长度，所在染色体contig的长度）

# #prophage上的DSs密度
# rownames(prophage_dss) <- prophage_dss$genome
# prophage_dss <- prophage_dss[,-1]
# prophage_dss_sum <- rowSums(prophage_dss) %>% as.data.frame() %>% setNames('prophage_DSs')
# prophage_dss_sum$genome <- rownames(prophage_dss_sum)
# prophage_dss_sum$genome <- sub('.fasta', '',prophage_dss_sum$genome)
# 
# pv_dss_density <- merge(prophage_dss_sum,magsize,by = 'genome')
# pv_dss_density$pvdensity <- pv_dss_density$prophage_DSs/pv_dss_density$size*1000000
# pv_dss_density$type <- 'prophage'  
# mean(pv_dss_density$pvdensity)#0.5346812
# 
# #plasmid上的DSs密度
# rownames(plasmid_dss) <- plasmid_dss$genome
# plasmid_dss <- plasmid_dss[,-1]
# plasmid_dss_sum <- rowSums(plasmid_dss) %>% as.data.frame() %>% setNames('plasmid_DSs')
# plasmid_dss_sum$genome <- rownames(plasmid_dss_sum)
# plasmid_dss_sum$genome <- sub('.fasta', '',plasmid_dss_sum$genome)
# 
# pd_dss_density <- merge(plasmid_dss_sum,magsize,by = 'genome')
# pd_dss_density$pddensity <- pd_dss_density$plasmid_DSs/pd_dss_density$size*1000000
# pd_dss_density$type <- 'plasmid'  
# mean(pd_dss_density$pddensity)#0.5601555
# 
# #染色体上的DSs密度，分别尝试不计入/计入那7731个不含DSs的MAGs
# # 使用dplyr的left_join进行匹配
# ds_size_updated <- ds_size %>%
#   left_join(pd_dss_density[,c(1,2)], by = "genome") %>%
#   mutate(dss = ifelse(!is.na(plasmid_DSs), dss - plasmid_DSs, dss))
# 
# ds_size_updated <- ds_size_updated %>%
#   left_join(pv_dss_density[,c(1,2)], by = "genome") %>%
#   mutate(dss = ifelse(!is.na(prophage_DSs), dss - prophage_DSs, dss))
# 
# ds_size_updated <- ds_size_updated[,c(1,2,3,4,5)]
# ds_size_updated$dsdensity <- ds_size_updated$dss/ds_size_updated$size*1000000
# ds_size_updated$type <- 'Chromosome'
# 
# #合并
# colnames(pd_dss_density) <-
#   colnames(pv_dss_density) <-
#   colnames(ds_size_updated) <- 
#   c('genome','dsnum','magsize','N50','N90','dsdensity','type')
# 
# all <- rbind(ds_size_updated,pv_dss_density,pd_dss_density)
# mean(all[all$type=='Chromosome',]$dsdensity)
# 
# #统计分析比较
# library(rstatix)
# md_p_val <- all %>% 
#   wilcox_test(dsdensity ~ type) %>%
#   adjust_pvalue(p.col="p",method="bonferroni") %>%
#   add_significance(p.col="p.adj") %>% 
#   add_xy_position(x="region",dodge=0.5)
# 
# library(ggplot2)
# 
# ggplot(all[all$dsnum>0,],aes(x=type,y=dsdensity))+
#   geom_boxplot()


# #7. 探究prophage plasmid数量与DSs数量关系
# setwd('D:/wenjian/research/defensome/groundwater/03mag_info')
# tax <- read.csv('ds_density_tax_r220.csv')
# 
# prophage_sum_size$genome <- paste(prophage_sum_size$genome,'.fasta',sep = "")
# prophage_ds_mag <- merge(prophage_sum_size,tax,by = 'genome')
# 
# cor.test(prophage_ds_mag$prophagesum,
#          prophage_ds_mag$number,
#          alternative = "two.sided",
#          method = "spearman",
#          conf.level = 0.95
# )#ρ0.1，p<2.2e-16
# plot(prophage_ds_mag$size,
#      prophage_ds_mag$number)
# 
# plasmid_sum_size$genome <- paste(plasmid_sum_size$genome,'.fasta',sep = '')
# plasmid_ds_mag <- merge(plasmid_sum_size,tax,by = 'genome')
# 
# cor.test(plasmid_ds_mag$plasmidsum,
#          plasmid_ds_mag$number,
#          alternative = "two.sided",
#          method = "spearman",
#          conf.level = 0.95
# )#ρ=-0.00308，p=0.7071。用这个来说明MGE的影响小?
# plot(plasmid_ds_mag$plasmidsum,
#      plasmid_ds_mag$number)

#7. 计算染色体、质粒和病毒上的DSs density（DSs数量/对应的总长度）###############

#7.1读取数据，预处理

#Integron的contig总长度信息
setwd("D:/pku/01groundwater/01defensome/NC_revision/1additional_analyses/ICEs")
all_integrons <- read.delim("complete_integron_contig.id",header = F)
colnames(all_integrons) <- "contig_id"
all_integrons$genome <- gsub("\\.fasta.*","",all_integrons$contig_id)
library(stringr)
all_integrons$length <- str_extract(all_integrons$contig_id,"(?<=length_)\\d+")
all_integrons$length <- as.numeric(all_integrons$length)
integron_contig_size <- all_integrons %>% group_by(genome) %>% 
  summarise(integron_size = sum(length))

#integron的DS信息
integron_ds <- read.delim("all_integron_dss2.txt")
# 提取从开头到 cov_ 及其后面的数字部分
integron_ds$contig_id <- str_extract(integron_ds$sys_id, "^.*cov_\\d+(\\.\\d+)?")
integron_ds$genome <- gsub(".fasta","",integron_ds$genome)
integron_ds_sum <- integron_ds %>% 
  group_by(genome) %>% 
  summarize(Integron_DS_num = sum(n()))

#ICE的contig总长度信息
all_ICE <- read.delim("all_ICE_contigid.txt",header = F)
colnames(all_ICE) <- "contig_id"
all_ICE$genome <- gsub("\\_contig.*","",all_ICE$contig_id)
all_ICE$genome <- gsub(">","",all_ICE$genome)
# 使用正则表达式提取 start 和 end 位置
all_ICE$start <- as.numeric(gsub(".*_(\\d+)-(\\d+)$", "\\1", all_ICE$contig_id))
all_ICE$end <- as.numeric(gsub(".*_(\\d+)-(\\d+)$", "\\2", all_ICE$contig_id))
# 计算长度并存储到新列 length 中
all_ICE$length <- all_ICE$end - all_ICE$start + 1
ICE_contig_size <- all_ICE %>%
  group_by(genome) %>%
  summarize(ICE_size = sum(length))

#ICE的DS信息
ICE_ds <- read.delim("all_ICE_dss2.txt")
# 提取从开头到 cov_ 及其后面的数字部分
#ICE_ds$contig_id <- str_extract(ICE_ds$sys_id, "^.*cov_\\d+(\\.\\d+)?")
ICE_ds$genome <- gsub("_contig","",ICE_ds$genome)
ICE_ds_sum <- ICE_ds %>% 
  group_by(genome) %>% 
  summarize(ICE_DS_num = sum(n()))

#prophage总长度信息
all_prophage <- read.delim("all_prophage.tsv",header = T)
all_prophage$genome <- gsub("\\.fasta.*","",all_prophage$source_seq)
prophage_contig_size <- all_prophage %>%
  group_by(genome) %>%
  summarize(prophage_size = sum(length))


#prophage的DS信息
prophage_ds <- read.delim("all_prophage_dss2.txt")
# 提取从开头到 cov_ 及其后面的数字部分
#ICE_ds$contig_id <- str_extract(ICE_ds$sys_id, "^.*cov_\\d+(\\.\\d+)?")
prophage_ds$genome <- gsub(".fasta","",prophage_ds$genome)
prophage_ds_sum <- prophage_ds %>% 
  group_by(genome) %>% 
  summarize(prophage_DS_num = sum(n()))



#质粒、病毒的contig总长度信息
setwd("D:/pku/01groundwater/01defensome/NC_revision/1additional_analyses/ICEs")
virus_contig <- read.delim2("all_prophage_contig.tsv")
plasmid_contig <- read.delim2("all_plasmid_contig.tsv")

virus_contig$genome <- sub(".fasta.*$", "", virus_contig$seq_name)
plasmid_contig$genome <- sub(".fasta.*$", "", plasmid_contig$seq_name)

library(dplyr)
virus_contig_size <- virus_contig %>% group_by(genome) %>% 
  summarise(virus_size = sum(length))
plasmid_contig_size <- plasmid_contig %>% group_by(genome) %>% 
  summarise(plasmid_size = sum(length))




#MAG的总长度信息
setwd("D:/wenjian/research/defensome/groundwater/03mag_info")
mag_ds <- read.csv("ds_density_tax_r220.csv")
mag_ds$genome <- sub(".fasta.*$", "", mag_ds$genome)
mag_ds$mag_size <- mag_ds$length_Mb*1000000

#合并
md <- merge(mag_ds,virus_contig_size,by = "genome",all = T)
md <- mag_ds %>%
  full_join(virus_contig_size, by = "genome") %>%
  full_join(plasmid_contig_size, by = "genome") %>%
  full_join(integron_contig_size, by = "genome") %>%
  full_join(ICE_contig_size, by = "genome") %>%
  full_join(prophage_contig_size, by = "genome")
md$phage_size <- md$virus_size - md$prophage_size
md$phage_size[md$phage_size < 0] <- 0
  
md$virus_size[is.na(md$virus_size)] <- 0
md$plasmid_size[is.na(md$plasmid_size)] <- 0
md$integron_size[is.na(md$integron_size)] <- 0
md$ICE_size[is.na(md$ICE_size)] <- 0
md$prophage_size[is.na(md$prophage_size)] <- 0

#计算染色体的size
md$chromosome <- md$mag_size-(md$virus_size+md$plasmid_size+md$integron_size+md$ICE_size)

#读取质粒、病毒上的DSs数量
setwd("D:/pku/01groundwater/01defensome/NC_revision/1additional_analyses/ICEs")

virus_dss <- read.delim2('all_virus_dss_counts_matrix.csv',sep = ',',check.names = F) 
virus_dss$genome <- gsub("*.fasta","",virus_dss$genome)
sum(virus_dss[,-1]) #1120个MAGs中，1774个DSs，106个DS families

plasmid_dss <- read.delim2('all_plasmid_dss_counts_matrix.csv',sep = ',',check.names = F)
plasmid_dss$genome <- gsub("*.fasta","",plasmid_dss$genome)
#sum(plasmid_dss[,-1]) #3278个MAGs中，6070个DSs

#sum(prophage_dss$genome %in% plasmid_dss$genome)#280个MAG中同时含有prophage DS和plasmid DS
#sum(prophage$Prefix %in% plasmid$Prefix)#8501个MAG中同时含有prophage和plasmid

virus_dss_sum <- data.frame(genome = virus_dss$genome,
                               virus_dss_number = rowSums(virus_dss[,-1])
)
plasmid_dss_sum <- data.frame(genome = plasmid_dss$genome,
                              Plasmid_dss_number = rowSums(plasmid_dss[,-1])
)

md <- md %>%
  full_join(virus_dss_sum,by = "genome") %>%
  full_join(plasmid_dss_sum,by = "genome") %>%
  full_join(integron_ds_sum,by = "genome") %>%
  full_join(ICE_ds_sum,by = "genome") %>%
  full_join(prophage_ds_sum,by = "genome")

md$virus_dss_number[is.na(md$virus_dss_number)] <- 0
md$Plasmid_dss_number[is.na(md$Plasmid_dss_number)] <- 0
md$Integron_DS_num[is.na(md$Integron_DS_num)] <- 0
md$ICE_DS_num[is.na(md$ICE_DS_num)] <- 0
md$prophage_DS_num[is.na(md$prophage_DS_num)] <- 0
md$phage_DS_num <- md$virus_dss_number - md$prophage_DS_num
md$phage_DS_num[md$phage_DS_num < 0] <- 0
md$chromosome_dss_number <- md$number -
  md$virus_dss_number -
  md$Plasmid_dss_number -
  md$Integron_DS_num -
  md$ICE_DS_num
md$chromosome_dss_number[md$chromosome_dss_number < 0] <- 0

#7.2计算各种MGE上的DSs密度

#prophage
prophage_ds_density <- select(md,prophage_size,prophage_DS_num)
prophage_ds_density <- prophage_ds_density[prophage_ds_density$prophage_size>0,]
prophage_ds_density$ds_density <- prophage_ds_density$prophage_DS_num/prophage_ds_density$prophage_size*1000
prophage_ds_density$location <- "Prophage"
colnames(prophage_ds_density)[1:2] <- c("size","dsnumber")

#phage上
phage_ds_density <- select(md,phage_size,phage_DS_num)
phage_ds_density <- phage_ds_density[phage_ds_density$phage_size>0,]
phage_ds_density$ds_density <- phage_ds_density$phage_DS_num/phage_ds_density$phage_size*1000
phage_ds_density$location <- "Phage"
colnames(phage_ds_density)[1:2] <- c("size","dsnumber")

#质粒上
plasmid_ds_density <- select(md,plasmid_size,Plasmid_dss_number)
plasmid_ds_density <- plasmid_ds_density[plasmid_ds_density$plasmid_size>0,]
plasmid_ds_density$ds_density <- plasmid_ds_density$Plasmid_dss_number/plasmid_ds_density$plasmid_size*1000
plasmid_ds_density$location <- "Plasmid"
colnames(plasmid_ds_density)[1:2] <- c("size","dsnumber")

#integron
integron_ds_density <- select(md,integron_size,Integron_DS_num)
integron_ds_density <- integron_ds_density[integron_ds_density$integron_size>0,]
integron_ds_density$ds_density <- integron_ds_density$Integron_DS_num/integron_ds_density$integron_size*1000
integron_ds_density$location <- "Integron"
colnames(integron_ds_density)[1:2] <- c("size","dsnumber")

#ICE
ICE_ds_density <- select(md,ICE_size,ICE_DS_num)
ICE_ds_density <- ICE_ds_density[ICE_ds_density$ICE_size>0,]
ICE_ds_density$ds_density <- ICE_ds_density$ICE_DS_num/ICE_ds_density$ICE_size*1000
ICE_ds_density$location <- "ICEs/IMEs"
colnames(ICE_ds_density)[1:2] <- c("size","dsnumber")

#染色体上
chromosome_ds_density <- select(md,chromosome,chromosome_dss_number)
chromosome_ds_density$ds_density <- chromosome_ds_density$chromosome_dss_number/chromosome_ds_density$chromosome*1000
chromosome_ds_density$location <- "Chromosome"
colnames(chromosome_ds_density)[1:2] <- c("size","dsnumber")

mge_ds_density <- rbind(prophage_ds_density,phage_ds_density,
                        plasmid_ds_density,integron_ds_density,ICE_ds_density,
                        chromosome_ds_density)

#7.3 统计分析
library(rstatix)
md_p_val <- mge_ds_density[mge_ds_density$dsnumber>0,] %>% 
  wilcox_test(ds_density ~ location) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="location",dodge=0.5)

library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means

#md_p_val$y.position <- c(0.4,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.8,0.84,0.88,0.92,0.96)
md_p_val$y.position <- seq(from = 0.3, to = 0.6, length.out = 15)
p_mge_ds_density <- 
  ggplot(mge_ds_density[mge_ds_density$dsnumber>0,],aes(x=location,y=ds_density))+
  geom_jitter(aes(group=1,color = location),width = 0.2, alpha = 0.3)+  # 显示异常值点
  geom_boxplot(aes(fill = location),width = 0.2,outlier.shape = NA, staplewidth = 0.5,linewidth=0.2) +  # 隐藏箱线图中的异常值点
  geom_violin(trim = T, aes(fill = location),show.legend = F,linewidth=0.2) + # 绘制小提琴图
  # stat_compare_means(comparisons = list(c("1", "2"),c("SS", "WA"),c("SS", "WS"),
  #                                       c("SA", "WA"),c("SA", "WS"),
  #                                       c("WA", "WS")
  # ),
  # method = "wilcox.test", test = "p.adjust",label = "p.signif")+
  stat_pvalue_manual(md_p_val,label = "p.adj.signif",hide.ns=T,
                     tip.length = 0,label.size = 5,color="black")+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,size=1.2,show.legend = F)+
  #stat_summary(fun=mean, geom="line", aes(group=CPR), col="black")+
  scale_fill_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                               "Integron"="#e8a598", "ICEs/IMEs"="#00BA38",
                               "Chromosome" = "#FCDABA",
                               "Patescibacteria" = "#9f7e69",
                               "Phage" = "#ED6764")) +
  scale_color_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                               "Integron"="#e8a598", "ICEs/IMEs"="#00BA38",
                               "Chromosome" = "#FCDABA",
                               "Phage" = "#ED6764")) +
  #ylim(0,15)+
  #scale_y_continuous(limits = c(0,10))+
  coord_cartesian(ylim = c(0,0.6))+
  labs(x=NULL,y='DS density (per kb)')+
  mytheme
p_mge_ds_density
setwd("D:/pku/01groundwater/01defensome/NC_revision/1additional_analyses/ICEs")
write.csv(mge_ds_density[mge_ds_density$dsnumber>0,], "fig5c.csv")

ggsave("1_mge_ds_density_new.pdf",p_mge_ds_density,width = 4.2,height = 3.5,units = "in")



mytheme <-  theme(
  panel.background = element_blank(),  # 设置背景为透明
  panel.grid.major = element_blank(),  # 去掉主要网格线
  panel.grid.minor = element_blank(),  # 去掉次要网格线
  panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
  axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
  axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
  axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
  axis.ticks = element_line(color = "black", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
  legend.position = 'none',
  #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
  legend.title = element_text( size = 12, color = "black"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
  legend.text = element_text(size = 10, color = "black"),  # 设置图例文本字体大小为10，颜色为深灰色
  legend.background = element_blank(),  # 设置图例背景为透明
  legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
  legend.key = element_blank()  # 去掉图例键的背景
)

#8.O/E值来自卡方检验；OR值来自fisher检验。二者类似。###################
library(dplyr)
##计算卡方值以及期望和观测值
#x <- chisq.test(matrix(c(80,300,200,220),ncol = 2))
#Roe <- x$observed / x$expected
#Roe
###T细胞在Cancer组OR值
#fisher.test(matrix(c(80,300,200,220),ncol = 2))

#8.1数据准备，每个MGE的各个denfese family的数量
#读取质粒、病毒上的DSs数量
setwd("D:/pku/01groundwater/01defensome/NC_revision/1additional_analyses/ICEs")

virus_dss <- read.delim2('all_virus_dss_counts_matrix.csv',sep = ',',check.names = F) 
virus_dss$genome <- gsub("*.fasta","",virus_dss$genome)
sum(virus_dss[,-1]) #1120个MAGs中，1774个DSs，106个DS families

plasmid_dss <- read.delim2('all_plasmid_dss_counts_matrix.csv',sep = ',',check.names = F)
plasmid_dss$genome <- gsub("*.fasta","",plasmid_dss$genome)

prophage_dss <- read.delim2("all_prophage_dss_counts_matrix.csv",sep = ',',check.names = F)
prophage_dss$genome <- gsub("*.fasta","",prophage_dss$genome)
prophage_dss <- prophage_dss %>%
  mutate(across(-1, as.numeric))

integron_dss <- read.delim2("all_integron_dss_counts_matrix.csv",sep = ',',check.names = F)
integron_dss$genome <- gsub("*.fasta","",integron_dss$genome)
integron_dss <- integron_dss %>%
  mutate(across(-1, as.numeric))

ICE_dss <- read.delim2("all_ICE_dss_counts_matrix.csv",sep = ',',check.names = F)
ICE_dss$genome <- gsub("*_contig","",ICE_dss$genome)
ICE_dss <- ICE_dss %>%
  mutate(across(-1, as.numeric))

#MAGs上所有dss
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
new_rows[col_names] <- lapply(col_names, function(col_name) rep(0, nrow(new_rows)))
rownames(new_rows) <- new_rows$genome
new_rows <- new_rows[,-1]
#合并new_rows与dss
dss <- rbind(dss,new_rows)

#把几个mge的ds family合并
ds_family_plasmid <- colSums(plasmid_dss[,-1]) %>% as.data.frame() %>% setNames("Plasmid")
ds_family_plasmid$dsfamily <- rownames(ds_family_plasmid)

ds_family_virus <- colSums(virus_dss[,-1]) %>% as.data.frame() %>% setNames("Virus")
ds_family_virus$dsfamily <- rownames(ds_family_virus)

ds_family_prophage <- colSums(prophage_dss[,-1]) %>% as.data.frame() %>% setNames("Prophage")
ds_family_prophage$dsfamily <- rownames(ds_family_prophage)
ds_family_prophage$dsfamily[ds_family_prophage$dsfamily=="Cas"] <- "CasFinder"

ds_family_integron <- colSums(integron_dss[,-1]) %>% as.data.frame() %>% setNames("Integron")
ds_family_integron$dsfamily <- rownames(ds_family_integron)

ds_family_ICE <- colSums(ICE_dss[,-1]) %>% as.data.frame() %>% setNames("ICEs/IMEs")
ds_family_ICE$dsfamily <- rownames(ds_family_ICE)

ds_family_all <- colSums(dss) %>% as.data.frame() %>% setNames("all")
ds_family_all$dsfamily <- rownames(ds_family_all)

ds_family_bind <- ds_family_all %>%
  full_join(ds_family_virus,by = 'dsfamily') %>%
  full_join(ds_family_plasmid,by = 'dsfamily') %>%
  full_join(ds_family_prophage,by = 'dsfamily') %>%
  full_join(ds_family_integron,by = 'dsfamily') %>%
  full_join(ds_family_ICE,by = 'dsfamily')

ds_family_bind[is.na(ds_family_bind)] <- 0

ds_family_bind <- ds_family_bind %>%
  rowwise() %>%
  mutate(Virus = ifelse(Virus == 0, Prophage, Virus))
# 将 all 列中为 0 的部分替换为其他列（除了第二列）的和
ds_family_bind <- ds_family_bind %>%
  rowwise() %>%
  mutate(all = ifelse(all == 0, sum(c_across(c(-2,-5))), all))



ds_family_bind$Chromosome <- ds_family_bind$all - 
  ds_family_bind$Plasmid - 
  ds_family_bind$Virus -
  ds_family_bind$Integron-
  ds_family_bind$`ICEs/IMEs`

ds_family_bind$Phage <- ds_family_bind$Virus - ds_family_bind$Prophage

ds_family_final <- select(ds_family_bind,dsfamily,Plasmid,Phage,
                          Prophage,Chromosome,Integron,`ICEs/IMEs`) %>% as.data.frame()

rownames(ds_family_final) <- ds_family_final$dsfamily
ds_family_final <- ds_family_final[,-1]
ds_family_final <- ds_family_final[ds_family_final$Chromosome>=0,]

#8.2 计算O/E，即卡方检验

oe <- chisq.test(ds_family_final,simulate.p.value = T)

oe_bi <- oe$observed/oe$expected %>% as.data.frame()
oe_bi <- na.omit(oe_bi)
oe_bi <- oe_bi[oe_bi$Chromosome!=0,]
oe_bi <- oe_bi[oe_bi$Plasmid!=0 & oe_bi$Phage!=0,]
#oe_bi <- oe_bi[oe_bi$Plasmid!=0 & oe_bi$Prophage!=0 & oe_bi$Integron!=0 & oe_bi$`ICEs/IMEs`!=0,]

#or <- fisher.test(ds_family_final[,-3],workspace = 1e+010)


#绘制热图
library(ComplexHeatmap)
library(circlize)  # 用于 colorRamp2 函数

# 定义颜色映射函数
col_fun <- colorRamp2(c(0,1, 2), c("white","#B6B8D8", "#FF7F0E"))  # 粉色到绿色的渐变
setwd("D:/pku/01groundwater/01defensome/NC_revision/1additional_analyses/ICEs")
cairo_pdf(width = 12,height = 2.8,filename = "1_mag_ds_oe2_new2.pdf")
p_mge_ds_oe <- Heatmap(t(oe_bi),
                       col = col_fun,
                       cluster_rows = F,
                       cluster_columns = T,
                       name = 'O/E',
                       row_title = NULL,
                       row_split = colnames(oe_bi),
                       row_names_side = "left", # 行标题显示在左侧
                       column_names_gp = gpar(fontsize = 8),
                       column_names_rot = 60,
                       rect_gp = gpar(col = "white", lwd = 0.5),    ##热图矩阵线设置
                       heatmap_legend_param = list(
                         title = "O/E", at = c(0,1, 2), 
                         labels = c("0","1", "≥2")),
)
p_mge_ds_oe
dev.off()
setwd("D:/pku/01groundwater/01defensome/NC_revision/1additional_analyses/ICEs")
write.csv(oe_bi,"fig5d.csv")

draw(p_mge_ds_oe)
library(ggplot2)
ggsave("1_mag_ds_oe.pdf",p_mge_ds_oe,width = 13,height = 2,units = 'in')
library(ggpubr)
ggexport("1_mag_ds_oe.pdf",p_mge_ds_oe,width = 13,height = 2,units = 'in')
