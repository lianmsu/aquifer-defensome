# =============================================================================
# aquifer-defensome -- supplementary figures / Function of non-defense genes in defense islands
# =============================================================================
# COG composition of the eggNOG annotation of the non-defense genes found
# inside defense islands.
#
# Manuscript      : Fig. S16
# Main outputs    : COG composition figure
# Original script : R codes/Fig.S16.R
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

#探究DIs上非DGs基因的功能组成
library(dplyr)
library(ggplot2)


#1. 读取数据
setwd(file.path(GTR_DATA, "groundwater/04mag_DI_gene"))
eggnog <- readxl::read_excel('out.emapper.annotations.xlsx')
eggnog <- eggnog[1:1867,]

#2. 探究COG组成####
# 使用 dplyr 包的 count 函数来统计 COG_category 列的各个字符串的数量
category_counts <- eggnog %>%
  count(COG_category)

# # 首先，过滤出 n 列小于20的行
# filtered_counts2 <- category_counts %>%
#   filter(n < 10)
# other <- data.frame(COG_category='Others',
#                     n=sum(filtered_counts2$n)
#                     )
# # 首先，过滤出 n 列大于20的行
# filtered_counts1 <- category_counts %>%
#   filter(n > 10)
# 
# category_counts <- rbind(filtered_counts1,other)
# category_counts$COG_category[1] <- 'Unknown'

library(ggpubr)
p1 <- ggdotchart(category_counts,x = 'COG_category', y = 'n',
           color = "#00AFBB",                                # Color by groups
          # palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
         #  group = "cyl",                                # Order by groups
           dot.size = 5,                                 # Large dot size
         #  label = round(dfm$mpg_z,1),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 9,
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           ) +
      labs(x=NULL,y='Gene number')+
  annotate("text", x = Inf, y = Inf, label = "COG annotation", hjust = 1.1,
           vjust = 2, size = 5, color = "black", parse = FALSE)+
      theme(
        panel.background = element_blank(),  # 设置背景为透明
        panel.grid.major = element_blank(),  # 去掉主要网格线
        panel.grid.minor = element_blank(),  # 去掉次要网格线
        panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
        axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
        axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
        axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
        axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
        legend.position = 'none',
        #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
        legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
        legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
        legend.background = element_blank(),  # 设置图例背景为透明
        legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
        legend.key = element_blank()  # 去掉图例键的背景
      )
p1      
setwd(file.path(GTR_DATA, "groundwater/04mag_DI_gene"))
ggsave(plot=p1,'1_DI_nonDGgene_COG.pdf',width = 10,height = 3,units = 'in')

#3. KO, PFAM等的组成
ko_counts <- eggnog %>%
  count(KEGG_ko)

pfam <-  eggnog %>%
  count(PFAMs)

p2 <- ggdotchart(ko_counts[ko_counts$n>=10 & ko_counts$n<1000,],x = 'KEGG_ko', y = 'n',
                 color = "#FC4E07",                                # Color by groups
                 # palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
                 sorting = "descending",                       # Sort value in descending order
                 add = "segments",                             # Add segments from y = 0 to dots
                 add.params = list(color = "lightgray", size = 2), # Change segment color and size
                 #  group = "cyl",                                # Order by groups
                 dot.size = 5,                                 # Large dot size
                 #  label = round(dfm$mpg_z,1),                        # Add mpg values as dot labels
                 font.label = list(color = "white", size = 9,
                                   vjust = 0.5)               # Adjust label parameters
                
) +
  annotate("text", x = Inf, y = Inf, label = "KEGG annotation", hjust = 1.1,
           vjust = 2, size = 5, color = "black", parse = FALSE)+
  labs(x=NULL,y='Gene number')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p2

p3 <- ggdotchart(pfam[pfam$n>=10 & pfam$n<400,],x = 'PFAMs', y = 'n',
                 color = "#E7B800",                                # Color by groups
                 # palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
                 sorting = "descending",                       # Sort value in descending order
                 add = "segments",                             # Add segments from y = 0 to dots
                 add.params = list(color = "lightgray", size = 2), # Change segment color and size
                 #  group = "cyl",                                # Order by groups
                 dot.size = 5,                                 # Large dot size
                 #  label = round(dfm$mpg_z,1),                        # Add mpg values as dot labels
                 font.label = list(color = "white", size = 9,
                                   vjust = 0.5),               # Adjust label parameters
                 ggtheme = theme_pubr()                        # ggplot2 theme
) +
  labs(x=NULL,y='Gene number')+
  annotate("text", x = Inf, y = Inf, label = "PFAM annotation", hjust = 1.1,
           vjust = 2, size = 5, color = "black", parse = FALSE)+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
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
p3

library(cowplot)
g1 <- plot_grid(p2,p3,labels = LETTERS[2:3], nrow = 1,rel_heights = c(1, 2))
g1

g2 <- plot_grid(p1,g1,labels = LETTERS[1], nrow = 2,rel_heights = c(1, 2))
g2

setwd(file.path(GTR_DATA, "groundwater/04mag_DI_gene"))
ggsave(plot = g2,'1_DI_nonDGgene_all.pdf',width = 6.7,height = 6.7,units = 'in')





