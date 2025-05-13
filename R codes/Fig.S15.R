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

plasmid_md$col2[1] <-"Unassigned"
prophage_md$col2[1] <- "Unassigned"
phage_md$col2[1] <- "Unassigned"
integron_md$col2[1] <- "Unassigned"

#####桑基图##############
library(tidyverse)
library(ggsankeyfier) 
library(MetBrewer)
df_ICE <- ICE_md %>% 
  pivot_stages_longer(.,stages_from = c("ICE", "col2"),
                      values_from = "ICEsum")

ICE_sankey <- 
  ggplot(data=df_ICE,aes(x = stage,y =ICEsum,group = node,
                      edge_id = edge_id,connector = connector))+
  # 绘制第 1，2 层级
  geom_sankeyedge(aes(fill = node),
                  position = position_sankey(order ="ascending",v_space="auto",
                                             width = 0.01))+
  geom_sankeynode(aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df_ICE %>% filter(connector=="from"),
            aes(label = node),stat = "sankeynode",
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
            hjust=0,size=4,fontface="plain",color="grey20")+
  geom_text(data=df_ICE %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3,vjust=0.5,color="grey20",fontface="plain")+
  geom_text(data=df_ICE %>% filter(connector=="to"),aes(label = ICEsum),
            stat = "sankeynode",
            position = position_sankey(width = 0.05,
                                       v_space ="auto",
                                       order="ascending",
                                       nudge_x=0.5),
            hjust=0,cex=2.8,color="black")+
  coord_cartesian(clip="off")+
  scale_x_discrete(position = "top")+
  scale_fill_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                               "Integron"="#e8a598", "ICEs/IMEs"="#9f7e69",
                               "Pseudomonadota" = "#FCDABA",
                               "Patescibacteria" = "#00BA38",
                               "Phage" = "#ED6764")) + 
  scale_color_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                                "Integron"="#e8a598", "ICEs/IMEs"="#9f7e69",
                                "Pseudomonadota" = "#FCDABA", 
                                "Patescibacteria" = "#00BA38",
                                "Phage" = "#ED6764")) +
  # scale_fill_manual(values = met.brewer("Nizami")) +
  # scale_color_manual(values = met.brewer("Nizami")) +
  theme_void()+
  theme(plot.margin = margin(0,0,0,0,unit = "cm"),
        axis.text.x=element_text(color="black",face="plain",size=12,
                                 margin = margin(b=0.05,unit = "cm")),
        legend.position="none")
ICE_sankey



# 筛选出数量小于10的行，并将它们合并为一行
plasmid_filtered <- plasmid_md %>%
  filter(plasmidsum < 700) %>%
  summarise(combined_plasum = sum(plasmidsum)) %>% setNames("plasmidsum")
plasmid_filtered$col2 <- 'Other phyla'
plasmid_filtered$plasmid <- 'Plasmid'
plasmid_filtered <- select(plasmid_filtered,col2,plasmidsum,plasmid)

plasmid_md <- plasmid_md[plasmid_md$plasmidsum>=700,]
plasmid_md <- rbind(plasmid_md,plasmid_filtered)
plasmid_md <- select(plasmid_md,plasmid,col2,plasmidsum)
df_plasmid <-  plasmid_md %>% 
  pivot_stages_longer(.,stages_from = c("plasmid", "col2"),
                      values_from = "plasmidsum")
Plasmid_sankey <- 
  ggplot(data=df_plasmid,aes(x = stage,y = plasmidsum,group = node,
                         edge_id = edge_id,connector = connector))+
  # 绘制第 1，2 层级
  geom_sankeyedge(aes(fill = node),
                  position = position_sankey(order ="ascending",v_space="auto",
                                             width = 0.01))+
  geom_sankeynode(aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df_plasmid %>% filter(connector=="from"),
            aes(label = node),stat = "sankeynode",
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
            hjust=0,size=4,fontface="plain",color="grey20")+
  geom_text(data=df_plasmid %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3,vjust=0.5,color="grey20",fontface="plain")+
  geom_text(data=df_plasmid %>% filter(connector=="to"),aes(label = plasmidsum),
            stat = "sankeynode",
            position = position_sankey(width = 0.05,
                                       v_space ="auto",
                                       order="ascending",
                                       nudge_x=0.5),
            hjust=0,cex=2.8,color="black")+
  coord_cartesian(clip="off")+
  scale_x_discrete(position = "top")+
  scale_fill_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                               "Integron"="#e8a598", "ICEs/IMEs"="#9f7e69",
                               "Pseudomonadota" = "#FCDABA",
                               "Patescibacteria" = "#00BA38",
                               "Phage" = "#ED6764")) + 
  scale_color_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                                "Integron"="#e8a598", "ICEs/IMEs"="#9f7e69",
                                "Pseudomonadota" = "#FCDABA", 
                                "Patescibacteria" = "#00BA38",
                                "Phage" = "#ED6764")) +
  # scale_fill_manual(values = met.brewer("Nizami")) +
  # scale_color_manual(values = met.brewer("Nizami")) +
  theme_void()+
  theme(plot.margin = margin(0,0,0,0,unit = "cm"),
        axis.text.x=element_text(color="black",face="plain",size=12,
                                 margin = margin(b=0.05,unit = "cm")),
        legend.position="none")
Plasmid_sankey

####Phage#
# 筛选出数量小于500的行，并将它们合并为一行
phage_filtered <- phage_md %>%
  filter(phagesum < 520) %>%
  summarise(combined_prosum = sum(phagesum)) %>% setNames("phagesum")
phage_filtered$col2 <- 'Other phyla'
phage_filtered$phage <- 'Phage'
phage_filtered <- select(phage_filtered,col2,phagesum,phage)

phage_md <- phage_md[phage_md$phagesum>=520,]
phage_md <- rbind(phage_md,phage_filtered)
phage_md <- select(phage_md,phage,col2,phagesum)

df_phage <-  phage_md %>% 
  pivot_stages_longer(.,stages_from = c("phage", "col2"),
                      values_from = "phagesum")
phage_sankey <- 
  ggplot(data=df_phage,aes(x = stage,y = phagesum,group = node,
                             edge_id = edge_id,connector = connector))+
  # 绘制第 1，2 层级
  geom_sankeyedge(aes(fill = node),
                  position = position_sankey(order ="ascending",v_space="auto",
                                             width = 0.01))+
  geom_sankeynode(aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df_phage %>% filter(connector=="from"),
            aes(label = node),stat = "sankeynode",
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
            hjust=0,size=4,fontface="plain",color="grey20")+
  geom_text(data=df_phage %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3,vjust=0.5,color="grey20",fontface="plain")+
  geom_text(data=df_phage %>% filter(connector=="to"),aes(label = phagesum),
            stat = "sankeynode",
            position = position_sankey(width = 0.05,
                                       v_space ="auto",
                                       order="ascending",
                                       nudge_x=0.5),
            hjust=0,cex=2.8,color="black")+
  coord_cartesian(clip="off")+
  scale_x_discrete(position = "top")+
  scale_fill_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                               "Integron"="#e8a598", "ICEs/IMEs"="#9f7e69",
                               "Pseudomonadota" = "#FCDABA",
                               "Patescibacteria" = "#00BA38",
                               "Phage" = "#ED6764")) + 
  scale_color_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                                "Integron"="#e8a598", "ICEs/IMEs"="#9f7e69",
                                "Pseudomonadota" = "#FCDABA", 
                                "Patescibacteria" = "#00BA38",
                                "Phage" = "#ED6764")) +
  # scale_fill_manual(values = met.brewer("Nizami")) +
  # scale_color_manual(values = met.brewer("Nizami")) +
  theme_void()+
  theme(plot.margin = margin(0,0,0,0,unit = "cm"),
        axis.text.x=element_text(color="black",face="plain",size=12,
                                 margin = margin(b=0.05,unit = "cm")),
        legend.position="none")
phage_sankey

#Prophage
# 筛选出数量10的行，并将它们合并为一行
prophage_filtered <- prophage_md %>%
  filter(prophagesum < 50) %>%
  summarise(combined_prosum = sum(prophagesum)) %>% setNames("prophagesum")
prophage_filtered$col2 <- 'Other phyla'
prophage_filtered$prophage <- 'Prophage'
prophage_filtered <- select(prophage_filtered,col2,prophagesum,prophage)

prophage_md <- prophage_md[prophage_md$prophagesum>=50,]
prophage_md <- rbind(prophage_md,prophage_filtered)
prophage_md <- select(prophage_md,prophage,col2,prophagesum)

df_prophage <-  prophage_md %>% 
  pivot_stages_longer(.,stages_from = c("prophage", "col2"),
                      values_from = "prophagesum")
prophage_sankey <- 
  ggplot(data=df_prophage,aes(x = stage,y = prophagesum,group = node,
                           edge_id = edge_id,connector = connector))+
  # 绘制第 1，2 层级
  geom_sankeyedge(aes(fill = node),
                  position = position_sankey(order ="ascending",v_space="auto",
                                             width = 0.01))+
  geom_sankeynode(aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df_prophage %>% filter(connector=="from"),
            aes(label = node),stat = "sankeynode",
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
            hjust=0,size=4,fontface="plain",color="grey20")+
  geom_text(data=df_prophage %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3,vjust=0.5,color="grey20",fontface="plain")+
  geom_text(data=df_prophage %>% filter(connector=="to"),aes(label = prophagesum),
            stat = "sankeynode",
            position = position_sankey(width = 0.05,
                                       v_space ="auto",
                                       order="ascending",
                                       nudge_x=0.5),
            hjust=0,cex=2.8,color="black")+
  coord_cartesian(clip="off")+
  scale_x_discrete(position = "top")+
  scale_fill_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                               "Integron"="#e8a598", "ICEs/IMEs"="#9f7e69",
                               "Pseudomonadota" = "#FCDABA",
                               "Patescibacteria" = "#00BA38",
                               "Phage" = "#ED6764")) + 
  scale_color_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                                "Integron"="#e8a598", "ICEs/IMEs"="#9f7e69",
                                "Pseudomonadota" = "#FCDABA", 
                                "Patescibacteria" = "#00BA38",
                                "Phage" = "#ED6764")) +
  # scale_fill_manual(values = met.brewer("Nizami")) +
  # scale_color_manual(values = met.brewer("Nizami")) +
  theme_void()+
  theme(plot.margin = margin(0,0,0,0,unit = "cm"),
        axis.text.x=element_text(color="black",face="plain",size=12,
                                 margin = margin(b=0.05,unit = "cm")),
        legend.position="none")
prophage_sankey


# Integron
integron_filtered <- integron_md %>%
  filter(integronsum < 45) %>%
  summarise(combined_prosum = sum(integronsum)) %>% setNames("integronsum")
integron_filtered$col2 <- 'Other phyla'
integron_filtered$integron <- 'Integron'
integron_filtered <- select(integron_filtered,col2,integronsum,integron)

integron_md <- integron_md[integron_md$integronsum>=45,]
integron_md <- rbind(integron_md,integron_filtered)
integron_md <- select(integron_md,integron,col2,integronsum)

df_integron <-  integron_md %>% 
  pivot_stages_longer(.,stages_from = c("integron", "col2"),
                      values_from = "integronsum")
integron_sankey <- 
  ggplot(data=df_integron,aes(x = stage,y = integronsum,group = node,
                              edge_id = edge_id,connector = connector))+
  # 绘制第 1，2 层级
  geom_sankeyedge(aes(fill = node),
                  position = position_sankey(order ="ascending",v_space="auto",
                                             width = 0.01))+
  geom_sankeynode(aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df_integron %>% filter(connector=="from"),
            aes(label = node),stat = "sankeynode",
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
            hjust=0,size=4,fontface="plain",color="grey20")+
  geom_text(data=df_integron %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3,vjust=0.5,color="grey20",fontface="plain")+
  geom_text(data=df_integron %>% filter(connector=="to"),aes(label = integronsum),
            stat = "sankeynode",
            position = position_sankey(width = 0.05,
                                       v_space ="auto",
                                       order="ascending",
                                       nudge_x=0.5),
            hjust=0,cex=2.8,color="black")+
  coord_cartesian(clip="off")+
  scale_x_discrete(position = "top")+
#  scale_fill_manual(values = met.brewer("Nizami")) +
  scale_fill_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                               "Integron"="#e8a598", "ICEs/IMEs"="#9f7e69",
                               "Pseudomonadota" = "#FCDABA",
                               "Patescibacteria" = "#00BA38",
                               "Phage" = "#ED6764")) + 
  scale_color_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                                "Integron"="#e8a598", "ICEs/IMEs"="#9f7e69",
                                "Pseudomonadota" = "#FCDABA", 
                                "Patescibacteria" = "#00BA38",
                                "Phage" = "#ED6764")) +
  # scale_color_manual(values = met.brewer("Nizami")) +
  theme_void()+
  theme(plot.margin = margin(0,0,0,0,unit = "cm"),
        axis.text.x=element_text(color="black",face="plain",size=12,
                                 margin = margin(b=0.05,unit = "cm")),
        legend.position="none")
integron_sankey


#ggsave("1_mge_phyla_sankey_new.pdf",mge_phyla_sankey,width = 5.56,height = 4.39,units = "in")

setwd("D:/pku/01groundwater/01defensome/NC_revision/1additional_analyses/ICEs")

write.csv(df_plasmid,"figs15a.csv")
write.csv(df_phage,"figs15b.csv")
write.csv(df_prophage,"figs15c.csv")
write.csv(df_integron,"figs15d.csv")
write.csv(df_ICE,"figs15e.csv")

library(Cairo)
library(cowplot)
#cairo_pdf(filename = "MGE_hosts.pdf" , width = 9,height = 6.5)
g1 <- plot_grid(Plasmid_sankey,phage_sankey,prophage_sankey,integron_sankey,
                ICE_sankey,
               labels = LETTERS[1:5], nrow = 2,ncol = 3
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
               )
g1



#setwd("D:/pku/01groundwater/03eu_pre/11eukdetect")
dev.off()



