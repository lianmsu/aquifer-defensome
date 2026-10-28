# =============================================================================
# aquifer-defensome -- supplementary figures / MAG abundance by niche breadth
# =============================================================================
# MAG abundance grouped by niche breadth and the stacked abundance
# composition per zone.
#
# Manuscript      : Fig. S2b
# Main outputs    : abs_stack.pdf, abs_stack_zone2.pdf
# Original script : R codes/Fig.S2_B.R
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

#对MAG按照生态位宽度分组？

#
setwd(file.path(GTR_DATA, "groundwater/s12mag_abundance"))

#mag丰度数据
mag_abun <- read.csv('bac.rpkm.tsv',sep = '\t')

#物种组成数据
tax <- read.delim('ds_density_tax_r220.csv',sep = ',')

#样品分区数据
setwd(file.path(GTR_DATA, "groundwater"))
gwdata <- readxl::read_xlsx('GW_DATA_lpw.xlsx')
colnames(gwdata) <- c('zhandian','dizhi','meta','q','jing_type',"chengya","zone")

tax$genome <- sub('.fasta','',tax$genome)
all <- merge(tax,mag_abun,by.x='genome',by.y='Genome')

spe <- all[-c(1:13)]
tax <- all[c(6:12)]
colnames(tax) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
#tax$Phylum <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", tax$Phylum)

# 2.1.1 物种组成数据按照门进行汇总
## spe和tax数据表中物种排序一致
library(tidyverse)
phy2 <- spe %>%
  group_by(tax$Phylum) %>% # 使用tax中的门水平进行分类
  summarise_all(sum) %>%
  rename(Phylum = `tax$Phylum`) %>%
  gather(key="samples",value = "abun",-Phylum) %>% # 数据形式转换：“宽”转“长”
  left_join(gwdata,by=c("samples"="meta")) %>%
  select(samples,zone,Phylum,abun) %>%
  group_by(samples,zone,Phylum) %>% # 求均值
  summarise_all(mean)

dim(phy)
head(phy)

# 2.1.2 颜色
library(ggsci)
col=pal_d3("category20")(20)
col2 = pal_d3("category20",alpha = 0.5)(20)
mypal=c(col,col2[-8])

# 2.1.3  物种组成堆叠柱形图-绝对丰度
library(ggplot2)
setwd(file.path(GTR_DATA, "groundwater/s12mag_abundance"))
pdf("abs_stack.pdf",width = 12,height = 10,family="Times")
ggplot()+
  geom_bar(data=phy2[phy2$abun>10,],
           aes(x=samples,
               weight=abun,
               fill=reorder(Phylum,-abun)),
           position = "stack",# 百分比堆叠图position = "fill"
           width=1)+
  facet_grid(.~zone)+
  scale_fill_manual(values = mypal[-18])+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  guides(fill=guide_legend(title = "Phylum",ncol = 2))+
  labs(y = "Absolute abundance",x=NULL)+
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", 
                                  size = 12, color = "black"),
        axis.text.x = element_blank(),  # 隐藏横轴标签
        axis.text.y = element_text(face = "bold", 
                                   size = 10, color = "black"),
        axis.ticks.x = element_blank(),  # 隐藏横轴刻度线
        strip.text.x = element_text(face = "bold", 
                                    size = 12, color = "black"),
        panel.grid = element_blank(),
        legend.title = element_text(face = "bold", 
                                    size = 12, color = "black"),
        legend.text = element_text(face = "bold", 
                                   size = 10, color = "black"))
dev.off()

# 对 spe 数据框的每一列进行归一化，使得每列求和为100
spe_normalized <- as.data.frame(lapply(spe, function(col) {
  col / sum(col, na.rm = TRUE) * 100
}))


phy <- spe_normalized %>%
  group_by(tax$Phylum) %>% # 使用tax中的门水平进行分类
  summarise_all(sum) %>%
  rename(Phylum = `tax$Phylum`) %>%
  gather(key="samples",value = "abun",-Phylum) %>% # 数据形式转换：“宽”转“长”
  left_join(gwdata,by=c("samples"="meta")) %>%
  select(zone,Phylum,abun) %>%
  group_by(zone,Phylum) %>% # 求均值
  summarise_all(mean)

# 将 abun 列值小于 1 的行的 Phylum 列改为 "Others"
phy$Phylum[phy$abun < 1] <- "Others"

library(dplyr)

# 合并 Phylum 列为 "Others" 且 zone 列相同的行
phy_merged <- phy %>%
  group_by(zone, Phylum) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>% 
  ungroup()

# 保留非 Others 的行，并合并 Others 行
phy_final <- phy %>%
  filter(Phylum != "Others") %>%
  bind_rows(phy_merged %>% filter(Phylum == "Others"))

dim(phy_final)
head(phy_final)

# 2.1.2 颜色
library(ggsci)
col=pal_d3("category20")(20)
col2 = pal_d3("category20",alpha = 0.5)(20)
mypal=c(col,col2[-8])

# 2.1.3  物种组成堆叠柱形图-绝对丰度
library(ggplot2)
setwd(file.path(GTR_DATA, "groundwater/s12mag_abundance"))
pdf("abs_stack_zone2.pdf",width = 8.5,height = 7)

phy_final$zone_roman <- as.roman(phy_final$zone) %>% as.character()

ggplot()+
  geom_bar(data=phy_final,
           aes(x=zone_roman,
               weight=abun,
               fill=reorder(Phylum,-abun)),
           position = "stack",# 百分比堆叠图position = "fill"
           width=0.5)+
  #facet_grid(.~zone)+
  scale_fill_manual(values = mypal[-18])+
  #scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  guides(fill=guide_legend(title = "Phylum",nrow =4))+
  labs(y = "Relative abundance (%)",x=NULL)+
  theme(legend.position = "bottom",legend.direction = "horizontal",legend.box = "bottom",
        axis.title = element_text(face = "bold", 
                                  size = 12, colour = "black"),
        #axis.text.x = element_blank(),  # 隐藏横轴标签
        axis.text.y = element_text(
                                   size = 10, color = "black"),
        strip.text.x = element_text(face = "bold", 
                                    size = 12, color = "black"),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),  # 设置图例项的大小
        #legend.text = element_text(size = 8),  # 设置图例文本的大小
        legend.text = element_text(face = "bold", 
                                   size = 8, color = "black"))

dev.off()

