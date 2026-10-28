# =============================================================================
# aquifer-defensome -- supplementary figures / VHR and genome size
# =============================================================================
# Prophage number per MAG vs MAG size, VHR vs genome size, and the features
# of 'supercarrier' MAGs.
#
# Manuscript      : Fig. S7
# Main outputs    : 1_prophage_MAGsize.pdf
# Original script : R codes/Fig.S7.R
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
GTR_REV <- normalizePath(Sys.getenv("GTR_REV", unset = "."), mustWork = FALSE)


#1每个MAG size与其中找到的phage的数量的关系。
library(dplyr)
#1. 读取数据###########
setwd(file.path(GTR_DATA, "groundwater/4VHR_DHR"))
VHR <- readxl::read_xlsx("VHR.xlsx")

#1.2 每个MAG中的质粒、病毒、DSs、DGs、及物种信息
setwd(file.path(GTR_DATA, "groundwater/03mag_info"))
mag_info <- read.delim2("ds_density_tax_r220.csv",sep = ',')
mag_info$genome <- gsub("*.fasta","",mag_info$genome)

setwd(file.path(GTR_DATA, "groundwater/05MGE"))
prophage <- read.delim2('virus_taxonomy_count.tsv')
rownames(prophage) <- prophage[,1]
sum(prophage[,-1]) #40302个
prophage_sum <- rowSums(prophage[,-1]) %>% as.data.frame() %>% setNames('prophagesum')
prophage_sum$genome <- rownames(prophage_sum)

plasmid <- read.delim2('plasmid_taxonomy_count.tsv')
rownames(plasmid) <- plasmid[,1]
sum(plasmid[,-1]) #163897
plasmid_sum <- rowSums(plasmid[,-1]) %>% as.data.frame() %>% setNames("plasmidsum")
plasmid_sum$genome <- rownames(plasmid_sum)

mag_ds_mge <- merge(mag_info,prophage_sum,by = 'genome', all = T)
mag_ds_mge <- merge(mag_ds_mge,plasmid_sum,by = 'genome', all = T)
mag_ds_mge[is.na(mag_ds_mge)] <- 0

DHR <- mag_ds_mge %>% group_by(col2) %>% 
  summarise(across(where(is.numeric),mean))

mytheme <-  theme(
  panel.background = element_blank(),  # 设置背景为透明
  panel.grid.major = element_blank(),  # 去掉主要网格线
  panel.grid.minor = element_blank(),  # 去掉次要网格线
  panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
  axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
  axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
  axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
  axis.ticks = element_line(color = "black", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
  legend.position = 'none',
  #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
  legend.title = element_text( size = 12, color = "black"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
  legend.text = element_text(size = 10, color = "black"),  # 设置图例文本字体大小为10，颜色为深灰色
  legend.background = element_blank(),  # 设置图例背景为透明
  legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
  legend.key = element_blank()  # 去掉图例键的背景
)

DHR <- DHR[DHR$prophagesum < 7.5,]
cor.test(DHR$prophagesum,DHR$magsize,method = "spearman")
phyla_gomomesize_phage <-
  ggplot(DHR,aes(x = magsize, y = prophagesum, color = "cyan3"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 7, y = 5, 
           label = "Spearman: ρ = 0.18 P = 0.04", vjust = 1.5, hjust = 1)+
  labs(x='MAG size (Mb)',y='Prophage number (per MAG)')+
  scale_color_manual(values = c( "#e6beae",
                                 "indianred1"
  )) +  # 手动设置每个category的颜色
  mytheme
phyla_gomomesize_phage

write.csv(DHR,"figs7b.csv")

#genomesize与phage/plasmid/all MGE的关系;对genomesize进行分组。
setwd(file.path(GTR_REV, "01defensome/NC_revision/1additional_analyses/genomesize_phage"))

#mag_ds_mge
#cor.test(mag_ds_mge$length_Mb,
#         mag_ds_mge$plasmidsum,method = "spearman")

#plot(mag_ds_mge$length_Mb,
#     mag_ds_mge$plasmidsum,)

df_summary2 <- mag_ds_mge %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(prophagesum),
    sd_value = sd(prophagesum)
  )
df_summary2$sd_value <- df_summary2$sd_value/10
#df_summary2$sd_value[7] <- 0
cor.test(df_summary2$magsize, df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)

cor.test(df_summary2[df_summary2$magsize<11,]$magsize, 
         df_summary2[df_summary2$magsize<11,]$mean_value,
         alternative = "two.sided",method = "spearman",conf.level = 0.95)

library(ggplot2)
size_phage <- ggplot(df_summary2, aes(x = magsize, y = mean_value)) +
  geom_point() +
  geom_line(linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.1,linewidth =0.3) +
  labs(x = 'MAG size (Mb)', y = 'Prophage number (per MAG)') +
  annotate("text", x = 10, y = 6, 
           label = "Spearman: ρ = 0.88 P = 0.0007", vjust = 1.5, hjust = 1)+
  #scale_color_manual(values = c("SA" = "#E9AF89", "SS" = "#9CB9D7", "WA" = "#B3DFC4", "WS" = "#DFA4BA")) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
size_phage
write.csv(df_summary2,"figs7a.csv")

#2VHR与MAG size关系。##########
setwd(file.path(GTR_REV, "01defensome/NC_revision/1additional_analyses/genomesize_phage"))
VHR <- read.csv("VHR-DHR.csv",header = T)

cor.test(VHR$logVHR, VHR$magsize, method = "spearman")
cor.test(VHR$VHR, VHR$magsize, method = "spearman")


mytheme <-  theme(
  panel.background = element_blank(),  # 设置背景为透明
  panel.grid.major = element_blank(),  # 去掉主要网格线
  panel.grid.minor = element_blank(),  # 去掉次要网格线
  panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
  axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
  axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
  axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
  axis.ticks = element_line(color = "black", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
  legend.position = 'none',
  #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
  legend.title = element_text( size = 12, color = "black"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
  legend.text = element_text(size = 10, color = "black"),  # 设置图例文本字体大小为10，颜色为深灰色
  legend.background = element_blank(),  # 设置图例背景为透明
  legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
  legend.key = element_blank()  # 去掉图例键的背景
)

logVHR_magsize <-
  ggplot(VHR,aes(x = logVHR, y = magsize, color = "cyan3"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 3.5, y = 5, 
           label = "Spearman: ρ = 0.17 P = 0.13", vjust = 1.5, hjust = 1)+
  labs(x='MAG size (Mb)',y='logVHR')+
  scale_color_manual(values = c( "#e6beae",
                                 "indianred1"
  )) +  # 手动设置每个category的颜色
  mytheme
logVHR_magsize
write.csv(VHR,"figs7c.csv")

DHR <- mag_ds_mge %>% group_by(col2) %>% 
  summarise(across(where(is.numeric),mean))
DHR <- DHR[DHR$prophagesum < 7.5,]
cor.test(DHR$prophagesum, DHR$number, method = "spearman")

DHR$density <- DHR$number/DHR$magsize
cor.test(DHR$prophagesum, DHR$density, method = "spearman")

p1 <- ggplot(DHR,aes(x = prophagesum, y = number, color = "cyan3"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 5, y = 20, 
           label = "Spearman: ρ = 0.32 P = 0.0003", vjust = 1.5, hjust = 1)+
  labs(x='Prophage number (per MAG)',y='DS number (per MAG)')+
  scale_color_manual(values = c( "cyan3",
                                 "indianred1"
  )) +  # 手动设置每个category的颜色
  mytheme
p1
write.csv(DHR,"figs7def.csv")

p2 <- ggplot(DHR,aes(x = prophagesum, y = density, color = "indianred1"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 5, y = 5, 
           label = "Spearman: ρ = 0.26 P = 0.003", vjust = 1.5, hjust = 1)+
  labs(x='Prophage number (per MAG)',y='DS density (per MAG)')+
  scale_color_manual(values = c( 
    "indianred1","cyan3"
  )) +  # 手动设置每个category的颜色
  mytheme
p2

#读取所有27578个mags的数据
setwd(file.path(GTR_DATA, "groundwater/02去冗余信息"))
allmags <- read.delim2('Wdb_ANI99.csv',sep = ',',
                       quote = "",header = T
) 

#读取含有dss的mags的数据
setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
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

# 使用apply函数计算每行大于0的数字的数量
ds_richness <- apply(dss, 1, function(x) sum(x > 0)) %>% as.data.frame() %>% setNames("ds_richness")
ds_richness$genome <- rownames(ds_richness)
ds_richness$genome <- gsub(".fasta","",ds_richness$genome)

mag_ds_mge <- merge(mag_ds_mge,ds_richness,by.x = "genome",by.y = "genome")

#统计
mag_ds_mge$density <- mag_ds_mge$density %>% as.numeric()

DHR <- mag_ds_mge %>% group_by(col2) %>% 
  summarise(across(where(is.numeric),mean))

DHR <- DHR[DHR$prophagesum < 7.5,]

cor.test(DHR$prophagesum, DHR$ds_richness, method = "spearman")
plot(DHR$plasmidsum, DHR$number)
p5 <- ggplot(DHR,aes(x = prophagesum, y = ds_richness, color = "#9c89b8"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 5, y = 8, 
           label = "Spearman: ρ = 0.27 P = 0.002", vjust = 1.5, hjust = 1)+
  scale_y_continuous(breaks = seq(0, max(DHR[-c(40,131),]$ds_richness), by = 2))+
  labs(x='Prophage number (per MAG)',y='DS richness (per MAG)')+
  scale_color_manual(values = c( "#9c89b8",
                                 "indianred1"
  )) +  # 手动设置每个category的颜色
  mytheme
p5

library(Cairo)
library(cowplot)
cairo_pdf("1_prophage_MAGsize.pdf",width = 8,height = 5)
g <- plot_grid(size_phage,phyla_gomomesize_phage,logVHR_magsize,p1,p2,p5,
               labels = LETTERS[1:6], nrow = 2
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_REV, "01defensome/NC_revision/1additional_analyses/genomesize_phage"))
dev.off()



#3supercarrier是否有更多的phage，OGT等#####################






