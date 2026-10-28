# =============================================================================
# aquifer-defensome -- main figures / Arms-race figure assembly
# =============================================================================
# Assembly fragment of Fig. 6: draws the prophage/anti-defense gene maps
# from the objects built by FigS17-18.R.
#
# Manuscript      : Fig. 6
# Main outputs    : p_prophage_anti_colo4.pdf, p_acrcascontig.pdf (objects built in
#                   FigS17-18.R)
# Original script : R codes/Fig6.R
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


##Fig6A###########################
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
setwd(file.path(GTR_DATA, "groundwater/31gene_circle"))
write.csv(prophage_anti_coloca_4,"Fig6A.csv",quote = F)
ggsave("p_prophage_anti_colo4.pdf",p_prophage_anti_colo4,width = 6.8,height = 4.3,units = 'in')

##Fig6C################################
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

setwd(file.path(GTR_DATA, "groundwater/31gene_circle/acr_cas_two_contig"))
write.csv(final_df,"Fig6C.csv",quote = F)
ggsave("p_acrcascontig.pdf",p_acrcascontig,width = 6.8,height = 3,units = 'in')

