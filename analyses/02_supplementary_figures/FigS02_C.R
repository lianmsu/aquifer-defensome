# =============================================================================
# aquifer-defensome -- supplementary figures / DS density across geological settings
# =============================================================================
# DS density by geological environment zone, well type, confined/unconfined
# aquifer and well depth; upset plot of defense families. Superseded by
# FigS08.R (which adds the source-data export).
#
# Manuscript      : Fig. S2c
# Main outputs    : 4_geographic_analyses.pdf, 4_log10DSs_density_jing-type2.pdf,
#                   4_p_DSs_log10density_dizhihuanjing2.pdf, DSs_density_chengya2.pdf
# Original script : R codes/Fig.S2_C.R
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

#前面这些是DSs density####

#！！！在绘图时，scale_y_continuous(limits = c(0,10))以及ylim(0,5)这种设置y轴，
#会改变绘制出来的图形结果。需要使用coord_cartesian(ylim = c(0,5))来放大局部而不改变数据使用。

#数据读取####
#1. 读取env数据
setwd(file.path(GTR_DATA, "groundwater/02去冗余信息"))
env <- readxl::read_xlsx('env.xlsx',sheet = 3)

setwd(file.path(GTR_DATA, "groundwater"))
gwdata <- readxl::read_xlsx('GW_DATA_lpw.xlsx')
colnames(gwdata) <- c('zhandian','dizhi','meta','q','jing_type',"chengya","zone")
gwdata$zone2 <- as.roman(gwdata$zone)


#2. 读取MAGs中的DSs数据
setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
#1.读取数据表###
#读取所有27578个mags的数据
setwd(file.path(GTR_DATA, "groundwater/02去冗余信息"))
allmags <- read.delim2('Wdb_ANI99.csv',sep = ',',
                       quote = "",header = T
) 

#读取含有dss的mags的数据
setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
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

#比较不同地质环境分区的MAGs中的DSs density####

#3. 计算每个MAGs中的DSs密度
setwd(file.path(GTR_DATA, "groundwater/03mag_info"))
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
# library(dplyr)
# df_summary <- ds_density %>%
#   group_by(magsize) %>%
#   summarise(
#     mean_value = mean(density),
#     sd_value = sd(density)
#   )
# df_summary2 <- ds_density %>%
#   group_by(magsize) %>%
#   summarise(
#     mean_value = mean(number),
#     sd_value = sd(number)
#   )


#4. 为ds_density添加地址环境分区列
#4.1 
library(stringr)
library(dplyr)
ds_density <- ds_density %>%
  mutate(ID = str_extract(genome, "^[^_]+"))

#4.2
ds_density <- ds_density %>%
  mutate(region = if_else(!is.na(match(ID, gwdata$meta)), 
                          gwdata$zone[match(ID, gwdata$meta)], NA))
#转换为因子类型
ds_density$region <- as.character(ds_density$region)
ds_density$region <- ordered(ds_density$region,
                             levels = c("1", "2", "3","4","5","6","7")) # 替换为实际的水平名称
ds_density <- na.omit(ds_density)

sum(ds_density$region=='1') #4809
sum(ds_density$region=='2') #12233
sum(ds_density$region=='3') #2966
sum(ds_density$region=='4') #3769
sum(ds_density$region=='5') #1659
sum(ds_density$region=='6') #1027
sum(ds_density$region=='7') #1115

#5. 绘制箱线+小提琴+散点复合图
library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means

#计算显著性差异
library(rstatix)
md_p_val1 <- ds_density %>% 
  wilcox_test(density ~ region) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="region",dodge=0.5)

# ds_density$log10density <- log10(ds_density$density+1)
# md_p_val2 <- ds_density %>% 
#   wilcox_test(log10density ~ region) %>%
#   adjust_pvalue(p.col="p",method="bonferroni") %>%
#   add_significance(p.col="p.adj") %>% 
#   add_xy_position(x="region",dodge=0.5)

md_p_val1$y.position <- seq(from = 0.6, to = 0.6 + (0.03 * (21 - 1)), by = 0.03)

p <- ggplot(ds_density, aes(x = region, y = log10(density+1))) +
  geom_boxplot(aes(fill = region),outliers = F) +  # 隐藏箱线图中的异常值点
  #geom_jitter(width = 0.2, alpha = 0.7)+  # 显示异常值点
  # stat_compare_means(comparisons = list(c("1", "2"),c("SS", "WA"),c("SS", "WS"),
  #                                       c("SA", "WA"),c("SA", "WS"),
  #                                       c("WA", "WS")
  # ),
  # method = "wilcox.test", test = "p.adjust",label = "p.signif")+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  stat_pvalue_manual(md_p_val1,label = "p.adj.signif",hide.ns=T,
                     tip.length = 0,label.size = 5,color="black")+
  scale_fill_manual(values = c("1"="#FCB2AF","2" = "#9BDFDF",
                               "3" = "#FFE2CE", "4" = "#C4D8E9",
                               "5"="#BEBCDF","6"="#FB8C62","7"="#FFF3CA"
                               )) +
  #ylim(0,15)+
  #scale_y_continuous(limits = c(0,10))+
  #coord_cartesian(ylim = c(0,5))+
  labs(x=NULL,y='log10 (DS density + 1)')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "black", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "plain", size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "black", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.8,0.7),  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "black"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "black"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p
setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
ggsave("4_p_DSs_log10density_dizhihuanjing2.pdf", plot = p, width = 6, height = 5, units = "in")

#按照新建改建井分类####
#读取分类数据
setwd(file.path(GTR_DATA, "groundwater/01样品信息"))
gwdata <- readxl::read_xlsx('GW_DATA.xlsx')
colnames(gwdata) <- c('zhandian','dizhi','16s','meta','zone','jing_type','chengya')

#
# ds_density <- ds_density %>%
#   mutate(region = if_else(!is.na(match(ID, gwdata$meta)), 
#                           gwdata$zone[match(ID, gwdata$meta)], NA))
ds_density <- ds_density %>%
  mutate(jingtype = if_else(!is.na(match(ID, gwdata$meta)), 
                          gwdata$jing_type[match(ID, gwdata$meta)], NA))
# ds_density <- ds_density %>%
#   mutate(chengya = if_else(!is.na(match(ID, gwdata$meta)), 
#                           gwdata$chengya[match(ID, gwdata$meta)], NA))

#删除jingtype为NA和0的行
ds_density <- na.omit(ds_density)
ds_density <- ds_density[ds_density$jingtype != 0,]

md_p_val2 <- ds_density %>% 
  wilcox_test(density ~ jingtype) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="jingtype",dodge=0.5)

md_p_val2$y.position <- c(0.8,0.9,1.0)

p2 <- ggplot(ds_density, aes(x = jingtype, y = log10(density+1))) +
  geom_boxplot(aes(fill = jingtype),outliers = F) +  # 隐藏箱线图中的异常值点
  #geom_jitter(width = 0.2, alpha = 0.7)+  # 显示异常值点
  # stat_compare_means(comparisons = list(c("1", "2"),c("SS", "WA"),c("SS", "WS"),
  #                                       c("SA", "WA"),c("SA", "WS"),
  #                                       c("WA", "WS")
  # ),
  # method = "wilcox.test", test = "p.adjust",label = "p.signif")+
  stat_pvalue_manual(md_p_val2,label = "p.adj.signif",hide.ns=T,
                     tip.length = 0,label.size = 5,color="black")+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
 # stat_summary(fun=mean, geom="line", aes(group=jingtype), col="black")+
  scale_fill_manual(values = c("R"="#FDBE6F","NC" = "#1F78B4","NP" = "#A5CEE2"   #"R"="#CFD99D","NC" = "#8F85C2","NP" = "#00A4B1"
  )) +
  #ylim(0,5)+
 # coord_cartesian(ylim = c(0,5))+
  labs(x=NULL,y='log10 (DS density + 1)')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.8,0.7),  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p2
setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
ggsave("4_log10DSs_density_jing-type2.pdf", plot = p2, width = 4, height = 5, units = "in")

#按照承压/潜水来比较MAGs中DSs density####
#读取分类数据
setwd(file.path(GTR_DATA, "groundwater/01样品信息"))
gwdata <- readxl::read_xlsx('GW_DATA.xlsx')
colnames(gwdata) <- c('zhandian','dizhi','16s','meta','zone','jing_type','chengya')

ds_density <- ds_density %>%
  mutate(chengya = if_else(!is.na(match(ID, gwdata$meta)),
                          gwdata$chengya[match(ID, gwdata$meta)], NA))

#只提取chengya为10和20的样本
ds_density <- na.omit(ds_density)
ds_density_2 <- ds_density[ds_density$chengya == 10 | ds_density$chengya == 20,]
#将chengya列的20替换为confined groundwater，10替换为phreatic water
ds_density_2 <- ds_density_2 %>%
  mutate(chengya = if_else(chengya == 20, "confined groundwater",
                           if_else(chengya == 10, "phreatic water", as.character(chengya))))

#计算分组统计量
md_p_val3 <- ds_density_2 %>% 
  wilcox_test(density ~ chengya) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="chengya",dodge=0.5)

#绘图

p3 <- ggplot(ds_density_2, aes(x = chengya, y = log10(density+1))) +
  geom_boxplot(aes(fill = chengya),outliers = F) +  # 隐藏箱线图中的异常值点
  #geom_jitter(width = 0.2, alpha = 0.7)+  # 显示异常值点
  # stat_compare_means(comparisons = list(c("1", "2"),c("SS", "WA"),c("SS", "WS"),
  #                                       c("SA", "WA"),c("SA", "WS"),
  #                                       c("WA", "WS")
  # ),
  # method = "wilcox.test", test = "p.adjust",label = "p.signif")+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,show.legend = F)+
  stat_pvalue_manual(md_p_val3,label = "p.adj.signif",hide.ns=F,
                     tip.length = 0,label.size = 5,color="black")+
  scale_fill_manual(values = c("confined groundwater"="#7892C7","phreatic water" = "#ADD9EE"
  )) +
  #ylim(0,5)+
  coord_cartesian(ylim = c(0,1.2))+
  labs(x=NULL,y='log10 (DS density + 1)')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.8,0.7),  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p3

setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
ggsave("DSs_density_chengya2.pdf", plot = p3, width = 5.5, height = 3.8, units = "in")


#分析井深的影响#######################
setwd(file.path(GTR_DATA, "groundwater/01样品信息"))
gwdata <- readxl::read_xlsx('GW_DATA.xlsx')
colnames(gwdata) <- c('zhandian','dizhi','16s','meta','zone','jing_type','chengya')

ds_density <- data.frame(genome=ds_size$genome,
                         density=ds_size$dss/ds_size$size*1000000,
                         number=ds_size$dss,
                         length_Mb=ds_size$size/1000000
)

ds_density <- ds_density %>%
  mutate(ID = str_extract(genome, "^[^_]+"))


depth <- env[,c(1,10)]
colnames(depth) <- c('ID',"Depth")
depth$Depth[depth$ID == "YU28"] <- 600
ds_density_depth <- merge(ds_density,depth,by = "ID")

cor.test(ds_density_depth$Depth,
         ds_density_depth$number,method = "spearman")
plot(ds_density_depth$Depth,
     ds_density_depth$number)

cor.test(ds_density_depth$Depth,
         ds_density_depth$density,method = "spearman")
plot(ds_density_depth$Depth,
     ds_density_depth$density)



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

depth_DSnumber <-
  ggplot(ds_density_depth,aes(x = Depth, y = number, color = "#809BCE"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 500, y = 60, 
           label = "Spearman: ρ = 0.02 P = 0.02", vjust = 1.5, hjust = 1)+
  labs(x='Well depth',y='DS number (per MAG)')+
  scale_color_manual(values = c( "#809BCE",
                                 "indianred1"
  )) +  # 手动设置每个category的颜色
  mytheme
depth_DSnumber

depth_DSdensity <-
  ggplot(ds_density_depth,aes(x = Depth, y = density, color = "#B8DFD2"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 500, y = 20, 
           label = "Spearman: ρ = 0.00 P = 0.58", vjust = 1.5, hjust = 1)+
  labs(x='Well depth',y='DS density (per MAG per Mb)')+
  scale_color_manual(values = c( "#B8DFD2",
                                 "indianred1"
  )) +  # 手动设置每个category的颜色
  mytheme
depth_DSdensity

setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
library(Cairo)
library(cowplot)
cairo_pdf("4_geographic_analyses.pdf",width = 6,height = 8)

g1 <- plot_grid(p2,p3,
               labels = LETTERS[2:3], nrow = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1
g2 <- plot_grid(depth_DSnumber,depth_DSdensity,
               labels = LETTERS[4:5], nrow = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g2


g <- plot_grid(p, g1,g2,
               labels = LETTERS[1:6], nrow = 3
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g

setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
dev.off()




#由于每个MAGs的DSs number与DSs density规律可能不一样。下面绘制DSs number####
#直接对上面代码进行修改，把density改为number

#对七大分区的defense families 绘制upset图######
dss$genome <- rownames(dss)
dss$site <- sapply(strsplit(dss$genome,"_"),`[`,1)
dss_zone <- merge(dss,gwdata,by.x = 'site',by.y = 'meta')

library(dplyr)
summarized_df <- dss_zone %>% group_by(zone2) %>%
  summarise(across(where(is.numeric),sum,na.rm=TRUE))%>% as.data.frame()

summarized_df <- select(summarized_df,-q,-chengya,-zone)
rownames(summarized_df) <- summarized_df$zone2
summarized_df <- summarized_df[,-1]

#计算七大分区DSs number总数
total_ds_zone <- data.frame(zone=rownames(summarized_df),
                            total_ds=rowSums(summarized_df)
                            )
library(ggplot2)
zone_ds <- ggplot(total_ds_zone,aes(x=reorder(zone,total_ds), y=total_ds,fill=zone)) +  #reorder调整排列顺序
  geom_bar(stat = "identity",  color="black", width = 0.5, position = position_dodge(0.5))+#调整条形宽度以及条形距离 )+ 
  coord_flip()+
  scale_fill_manual(values = c("I"="#FCB2AF","II"="#9BDFDF",
                               "III"="#FFE2CE","IV"= "#C4D8E9",
                               "V"= "#BEBCDF","VI"="#FB8C62","VII"="#FFF3CA"
  )) +
  geom_text(aes(label=total_ds), hjust=1, color="grey20")+
  xlab("Zone") + ylab("Number of DSs")+#可以通过color、size等自行调整标签属性
  theme(
    text = element_text(family = 'sans'),
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text( size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.8,0.7),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
zone_ds
ggexport(filename = "4_bar_zone_ds_sum.pdf", plot = zone_ds, width = 5, height = 3, units = "in")


t_ds_zone <- t(summarized_df) %>% as.data.frame()
t_ds_zone <- t_ds_zone %>% mutate(across(everything(),~if_else(. > 0,1,0)))

library(UpSetR)
plotzone_family <-  upset(t_ds_zone, sets = colnames(t_ds_zone), 
        sets.bar.color = c("II"="#9BDFDF","I"="#FCB2AF","IV"= "#C4D8E9",
                           "III"="#FFE2CE",
                           "V"= "#BEBCDF","VI"="#FB8C62","VII"="#FFF3CA"
        ),
        main.bar.color = "lightblue",
        order.by = "freq",
        point.size = 1.5,
        line.size = 0.5,
        text.scale = c(1.5,1.5,1,1,1.5,1.5),
        sets.x.label = 'Defense families',
        empty.intersections = "on")
plotzone_family

setwd(file.path(GTR_DATA, "groundwater/02ANI99MAGs_defense"))
#ggsave("4_upset_zone_ds_families.pdf", plot = plotzone_family, width = 6.5, height = 6, units = "in")
ggexport(filename = "4_upset_zone_ds_families.pdf", plot = plotzone_family, width = 8, height = 6, units = "in")

write.csv(total_ds_zone,"FigS2C_1.csv",quote = F)
write.csv(t_ds_zone,"FigS2C_2.csv",quote = F)
