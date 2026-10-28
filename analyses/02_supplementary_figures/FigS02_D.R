# =============================================================================
# aquifer-defensome -- supplementary figures / Defense families across zones
# =============================================================================
# Upset plot of the defense families found in the seven aquifer zones
# (archaeal MAGs).
#
# Manuscript      : Fig. S2d
# Main outputs    : 4_upset_zone_ds_families.pdf
# Original script : R codes/Fig.S2_D.R
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



#1. 读取数据####
#1. 读取env数据
setwd(file.path(GTR_DATA, "groundwater/02去冗余信息"))
env <- readxl::read_xlsx('env.xlsx',sheet = 3)

setwd(file.path(GTR_DATA, "groundwater"))
gwdata <- readxl::read_xlsx('GW_DATA_lpw.xlsx')
colnames(gwdata) <- c('zhandian','dizhi','meta','q','jing_type',"chengya","zone")
gwdata$zone2 <- as.roman(gwdata$zone)

#2. 读取MAGs中的DSs数据
setwd(file.path(GTR_DATA, "groundwater/archaea/02defense"))
ar_defense_7010 <- read.csv(file = 'ar_defense_7010.csv')



#对七大分区的defense families 绘制upset图######
ar_defense_7010$site <- sapply(strsplit(ar_defense_7010$genome,"_"),`[`,1)
dss_zone <- merge(ar_defense_7010,gwdata,by.x = 'site',by.y = 'meta')

summarized_df <- dss_zone %>% group_by(zone2) %>%
  summarise(across(where(is.numeric),sum,na.rm=TRUE))%>% as.data.frame()

summarized_df <- select(summarized_df,-q,-chengya,-zone)
rownames(summarized_df) <- summarized_df$zone2
summarized_df <- summarized_df[,-1]

#计算七大分区DSs number总数
total_ds_zone <- data.frame(zone=rownames(summarized_df),
                            total_ds=rowSums(summarized_df)
)

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
setwd(file.path(GTR_DATA, "groundwater/archaea/02defense"))
library(ggpubr)
ggexport(filename = "4_ar_bar_zone_ds_sum.pdf", plot = zone_ds, width = 5, height = 3, units = "in")


t_ds_zone <- t(summarized_df) %>% as.data.frame()
t_ds_zone <- t_ds_zone %>% mutate(across(everything(),~if_else(. > 0,1,0)))

library(UpSetR)
plotzone_family <-  upset(t_ds_zone, sets = colnames(t_ds_zone), 
                          sets.bar.color = c("II"="#9BDFDF","I"="#FCB2AF","IV"= "#C4D8E9",
                                             "III"="#FFE2CE",
                                             "V"= "#BEBCDF","VI"="#FB8C62","VII"="#FFF3CA"
                          ),
                          main.bar.color = "#EBB694",
                          order.by = "freq",
                          point.size = 1.5,
                          line.size = 0.5,
                          text.scale = c(1.5,1.5,1,1,1.5,1.5),
                          sets.x.label = 'Defense families',
                          empty.intersections = "on")
plotzone_family

setwd(file.path(GTR_DATA, "groundwater/archaea/02defense"))
#ggsave("4_upset_zone_ds_families.pdf", plot = plotzone_family, width = 6.5, height = 6, units = "in")
ggexport(filename = "4_ar_upset_zone_ds_families.pdf", plot = plotzone_family, width = 8, height = 6, units = "in")

write.csv(total_ds_zone,"FigS2D_1.csv",quote = F)
write.csv(t_ds_zone,"FigS2D_2.csv",quote = F)
