# =============================================================================
# aquifer-defensome -- main figures / VHR/DHR and the defensome
# =============================================================================
# Relationship between virus-to-host ratio (VHR), prophage/plasmid load and
# DS number, plus DPANN and CPR clade characteristics.
#
# Manuscript      : Fig. 3
# Main outputs    : VHR_DHR.pdf, 3CPR_clade_dss.pdf, 3_CPR_clade_DSrichness.pdf,
#                   1_dpann_phyla_DSrichness_10.pdf, 1_p11_DPANN_dss_liuxinglv.pdf
# Original script : R codes/Fig3.R
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


##Fig3ABC###############################
#VHR与DHR的关系;每个MAG中prophage number-DS number；plasmid number；
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

#2. 统计分析MAG上mge数量与DS数量之间关系，MAG上没有规律#########
mag_ds_mge_no0 <- mag_ds_mge[mag_ds_mge$number>0,]

# 进行Spearman相关性的假设检验
spearman_test <- cor.test(mag_ds_mge$number, mag_ds_mge$prophagesum, method = "spearman")
print(spearman_test)
plot(mag_ds_mge$number, mag_ds_mge$prophagesum)

spearman_test <- cor.test(mag_ds_mge$number, mag_ds_mge$plasmidsum, method = "spearman")
print(spearman_test)
plot(mag_ds_mge$number, mag_ds_mge$plasmidsum)

# 进行Spearman相关性的假设检验
spearman_test <- cor.test(mag_ds_mge_no0$number, mag_ds_mge_no0$prophagesum, method = "spearman")
print(spearman_test)
plot(mag_ds_mge_no0$number, mag_ds_mge_no0$prophagesum)

#3. 在门水平分析mge与DS数量关系,对各个门中的MAG计算DS数量平均，DS密度平均########

#3.1 统计
mag_ds_mge$density <- mag_ds_mge$density %>% as.numeric()
DHR <- mag_ds_mge %>% group_by(col2) %>% 
  summarise(across(where(is.numeric),mean))


# 对本身进行Spearman相关性的假设检验
cor.test(DHR$prophagesum, DHR$number, method = "spearman")
plot(DHR$prophagesum, DHR$number)

cor.test(DHR$prophagesum, DHR$density, method = "spearman")
plot(DHR$prophagesum, DHR$density)

cor.test(DHR$plasmidsum, DHR$number, method = "spearman")
plot(DHR$plasmidsum, DHR$number)

cor.test(DHR$plasmidsum, DHR$density, method = "spearman")
plot(DHR$plasmidsum, DHR$density)


#绘图
library(ggplot2)
p1 <- ggplot(DHR[-c(40,131),],aes(x = prophagesum, y = number, color = "cyan3"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 10, y = 20, 
           label = "Spearman: ρ = 0.29 P = 0.0008", vjust = 1.5, hjust = 1)+
  labs(x='Prophage number (per MAG)',y='DS number (per MAG)')+
  scale_color_manual(values = c( "cyan3",
                                 "indianred1"
  )) +  # 手动设置每个category的颜色
  mytheme
p1
p2 <- ggplot(DHR[-c(40,131),],aes(x = prophagesum, y = density, color = "indianred1"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 10, y = 4, 
           label = "Spearman: ρ = 0.28 P = 0.001", vjust = 1.5, hjust = 1)+
  labs(x='Prophage number (per MAG)',y='DS density (per MAG)')+
  scale_color_manual(values = c( 
    "indianred1","cyan3"
  )) +  # 手动设置每个category的颜色
  mytheme
p2
p3 <- ggplot(DHR[-c(40,131),],aes(x = plasmidsum, y = number, color = "cyan3"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 15, y = 20, 
           label = "Spearman: ρ = 0.38 P = 5.38e-06", vjust = 1.5, hjust = 1)+
  labs(x='Plasmid number (per MAG)',y='DS number (per MAG)')+
  scale_color_manual(values = c( "cyan3",
                                 "indianred1"
  )) +  # 手动设置每个category的颜色
  mytheme
p3
p4 <- ggplot(DHR[-c(40,131),],aes(x = plasmidsum, y = density, color = "indianred1"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 15, y = 4, 
           label = "Spearman: ρ = 0.22 P = 0.01", vjust = 1.5, hjust = 1)+
  labs(x='Plasmid number (per MAG)',y='DS density (per MAG)')+
  scale_color_manual(values = c( 
    "indianred1","cyan3"
  )) +  # 手动设置每个category的颜色
  mytheme
p4
p1
p2
p3
p4



#对VHR进行Spearman相关性的假设检验
HR <- merge(VHR,DHR,by.x = "phylum",by.y = "col2")
HR$logVHR <- log10(HR$VHR)
HR$lognumber <- log10(HR$number)

#setwd(file.path(GTR_DATA, "groundwater/4VHR_DHR"))
#write.csv(HR,"VHR-DHR.csv",quote = F)

cor.test(HR$VHR,HR$plasmidsum,method = "spearman")
plot(HR$logVHR,HR$density)

cor.test(HR$logVHR,HR$number,method = "spearman")
plot(HR$logVHR,HR$number)

cor.test(HR$logVHR,HR$density,method = "spearman")
plot(HR$logVHR,HR$density)

#3.2 绘图
library(ggplot2)

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

logVHR_DSnumber <-
  ggplot(HR,aes(x = logVHR, y = number, color = "cyan3"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 3.5, y = 14, 
           label = "Spearman: ρ = 0.3 P = 0.006", vjust = 1.5, hjust = 1)+
  labs(x='log10VHR',y='DS number (per MAG)')+
  scale_color_manual(values = c( "cyan3",
                                 "indianred1"
  )) +  # 手动设置每个category的颜色
  mytheme
logVHR_DSnumber

logVHR_DSdensity <-
  ggplot(HR,aes(x = logVHR, y = density, color = "indianred1"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 3.5, y = 4, 
           label = "Spearman: ρ = 0.22 P = 0.05", vjust = 1.5, hjust = 1)+
  labs(x='log10VHR',y='DS density (per MAG)')+
  scale_color_manual(values = c( 
    "indianred1","cyan3"
  )) +  # 手动设置每个category的颜色
  mytheme
logVHR_DSdensity


library(Cairo)
library(cowplot)
cairo_pdf(width = 6,height = 8)
g <- plot_grid(logVHR_DSnumber,logVHR_DSdensity,p1, p2,p3,p4,
               labels = LETTERS[1:6], nrow = 3
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "groundwater/4VHR_DHR"))
dev.off()

ggsave("VHR_DHR.pdf",g,width = 6,height = 8,units = 'in')

library(ggpubr)
ggexport(filename = "VHR_DHR.pdf",g,width = 6,height = 8,units = 'in')


#计算每个门中，DS richness的总数#############

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

cor.test(DHR$prophagesum, DHR$ds_richness, method = "spearman")
plot(DHR$plasmidsum, DHR$number)
p5 <- ggplot(DHR[-c(40,131),],aes(x = prophagesum, y = ds_richness, color = "#9c89b8"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 10, y = 8, 
           label = "Spearman: ρ = 0.26 P = 0.002", vjust = 1.5, hjust = 1)+
  scale_y_continuous(breaks = seq(0, max(DHR[-c(40,131),]$ds_richness), by = 2))+
  labs(x='Prophage number (per MAG)',y='DS richness (per MAG)')+
  scale_color_manual(values = c( "#9c89b8",
                                 "indianred1"
  )) +  # 手动设置每个category的颜色
  mytheme
p5

cor.test(DHR$plasmidsum, DHR$ds_richness, method = "spearman")
plot(DHR$plasmidsum, DHR$number)
p6 <- ggplot(DHR[-c(40,131),],aes(x = plasmidsum, y = ds_richness, color = "#9c89b8"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 15, y = 8.5, 
           label = "Spearman: ρ = 0.36 P = 1.778e-05", vjust = 1.5, hjust = 1)+
  scale_y_continuous(breaks = seq(0, max(DHR[-c(40,131),]$ds_richness), by = 2))+  # 设置纵轴刻度线为一位数字
  labs(x='Plasmid number (per MAG)',y='DS richness (per MAG)')+
  scale_color_manual(values = c( 
    "#9c89b8","cyan3"
  )) +  # 手动设置每个category的颜色
  mytheme
p6

HR <- merge(VHR,DHR,by.x = "phylum",by.y = "col2")
HR$logVHR <- log10(HR$VHR)
HR$lognumber <- log10(HR$number)


cor.test(HR$logVHR,HR$ds_richness,method = "spearman")
plot(HR$logVHR,HR$ds_richness)
logVHR_DSrichness <-
  ggplot(HR,aes(x = logVHR, y = ds_richness, color = "#9c89b8"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 3.5, y = 7, 
           label = "Spearman: ρ = 0.32 P = 0.004", vjust = 1.5, hjust = 1)+
  labs(x='log10VHR',y='DSs richness (per MAG)')+
  scale_color_manual(values = c( 
    "#9c89b8","cyan3"
  )) +  # 手动设置每个category的颜色
  mytheme
logVHR_DSrichness

dss$genome <- rownames(dss)
dss$genome <- gsub(".fasta","",dss$genome)
dss_tax <- merge(dss,mag_info,by.x = "genome",by.y = "genome")

dss_tax_phy <- dss_tax %>% group_by(col2) %>% 
  summarise(across(where(is.numeric),sum)) %>% as.data.frame()
rownames(dss_tax_phy) <- dss_tax_phy$col2

dss_tax_phy_dsrich <- apply(dss_tax_phy[,c(2:140)], 1, function(x) sum(x > 0)) %>% as.data.frame() %>% setNames("ds_richnesssum")
dss_tax_phy_dsrich$phyla <- rownames(dss_tax_phy_dsrich)

HR <- merge(HR,dss_tax_phy_dsrich,by.x = "phylum",by.y = "phyla")

cor.test(HR$logVHR,HR$ds_richnesssum,method = "spearman")
plot(HR$logVHR,HR$ds_richness)
logVHR_DSrichnesssum <-
  ggplot(HR,aes(x = logVHR, y = ds_richnesssum, color = "#9c89b8"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 3.5, y = 4, 
           label = "Spearman: ρ = 0.32 P = 0.004", vjust = 1.5, hjust = 1)+
  labs(x='log10VHR',y='DS richness (per MAG)')+
  scale_color_manual(values = c( 
    "#9c89b8","cyan3"
  )) +  # 手动设置每个category的颜色
  mytheme
logVHR_DSrichnesssum


library(Cairo)
library(cowplot)
cairo_pdf(width = 8,height = 8)
g <- plot_grid(logVHR_DSnumber,logVHR_DSdensity,logVHR_DSrichness,p1, p2,p5,p3,p4,p6,
               labels = LETTERS[1:9], nrow = 3
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "groundwater/4VHR_DHR"))
dev.off()

write.csv(HR,"Fig3ABC.csv",quote = F)

##Fig3D###################
#绘制热图
library(pheatmap)

# 设置颜色
#colors <- colorRampPalette(c('white', colorRampPalette('#F9B14F', '#8F2731')))(1000)
colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p <- pheatmap(t(result_data),
              #color = colorRampPalette(colors)(1000),
              #color = colorRampPalette(c('white','lightpink','red'))(1000),
              color = colors,
              cluster_cols = T,
              cluster_rows = F,
)
p

setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
write.csv(result_data,"Fig3D.csv",quote = F)
ggplot2::ggsave("3CPR_clade_dss.pdf", plot = p, width = 4, height = 6, units = "in")


##Fig3E###################################3


CPR_clade_DSrichness_20 <- ggplot(CPR_clade_res_richness_20,aes(x=col3,y=richness))+
  geom_bar(stat = 'identity',fill = "#FEDDCE", width = 0.5)+
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
CPR_clade_DSrichness_20
setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
write.csv(CPR_clade_res_richness_20,"Fig3E_1.csv",quote = F)
ggsave("3_CPR_clade_DSrichness.pdf", plot = CPR_clade_DSrichness_20, width = 3.6, height = 3.8, units = "in")

CPR20_density <- ggplot(CPR20_res, aes(x = col3, y = DSrichness_mean)) +
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
  geom_hline(yintercept = mean(CPR20_res$DSrichness_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DS richness (per MAG)") +
  scale_y_continuous(limits = c(0,2.8),labels = scales::number_format(accuracy = 0.1))+
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
CPR20_density

setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
write.csv(CPR20_res,"Fig3E_2.csv",quote = F)
ggsave(plot = CPR50_density,'3_CPR_DSrichnessperMAG.pdf',width = 5,height = 4,units = 'in')


CPR50_density <- ggplot(CPR50_res, aes(x = col3, y = density_mean)) +
  geom_bar(stat = "identity", fill = "#43A3EF",width = 0.5) +
  geom_errorbar(
    aes(ymin = density_mean - 1.96*density_se, ymax = density_mean + 1.96*density_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = density_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 3  # 点大小
  )+
  geom_hline(yintercept = mean(CPR50_res$density_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DS density (per MAG per Mb)") +
  scale_y_continuous(limits = c(0,2.5))+
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
CPR50_density

setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
ggsave(plot = CPR50_density,'3_dss_CPR.pdf',width = 5,height = 4,units = 'in')


##Fig3F###############################
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

setwd(file.path(GTR_DATA, "groundwater/archaea/02defense"))
write.csv(phyla_liuxing2,"Fig3F.csv",quote = F)
library(ggplot2)
ggsave("1_p11_DPANN_dss_liuxinglv.pdf", plot = p11, width = 4, height = 6, units = "in")

##Fig3G#############################
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
setwd(file.path(GTR_DATA, "groundwater/archaea/02defense"))
write.csv(dpann_phyla_res_richness_10,"Fig3G_1.csv",quote = F)
ggsave("1_dpann_phyla_DSrichness_10.pdf", plot = dpann_phyla_DSrichness_10, width = 3.5, height = 3.8, units = "in")


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

setwd(file.path(GTR_DATA, "groundwater/archaea/02defense"))
write.csv(dpann_res,"Fig3G_2.csv",quote = F)
ggsave(plot = dpann_density,'1_dpann_ds_richnessperMAG_phyla_10.pdf',width = 4,height = 4,units = 'in')


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

setwd(file.path(GTR_DATA, "groundwater/archaea/02defense"))
write.csv(dpann_density_res,"Fig3G_3.csv",quote = F)
ggsave(plot = dpann_density,'1_dpann_ds_density_phyla_10.pdf',width = 4,height = 4,units = 'in')




