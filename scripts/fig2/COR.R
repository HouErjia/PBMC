devtools::install_github("Github-Yilei/ggcor")
library(vegan)
library(dplyr)
library(ggcor)
library(ggplot2)
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2_add_Hou"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
BA=read.csv("./BA.csv",row.names = 1,check.names = F)
ENV=read.csv("./ENV.csv",row.names = 1,check.names = F)

spec=list()
for(i in 1:length(BA)){
  spec[[names(BA)[i]]]=i
}

df_mantel <- mantel_test(BA, ENV,
                         mantel.fun = 'mantel',
                         spec.dist.method = 'bray',
                         env.dist.method = 'euclidean',
                         spec.select = spec
                         # spec.select = list(`Primary BAs` = c(1,2,4),
                         #                    `Secondary BAs` = c(3,5:8)
                                            # )
                         )#将群落数据按组进行分开
df_mantel <- df_mantel %>%
  mutate(df_r = cut(r, breaks = c(-Inf, 0.25,0.5, Inf),
                    labels = c("< 0.25", "0.25 - 0.5",  ">= 0.5")),#定义Mantel的R值范围标签
         df_p = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("<=0.01", "0.01 - 0.05", ">= 0.05")))#定义Mantel的P值范围标签
quickcor(ENV,method = "spearman", type = "upper", cor.test = T, cluster.type = "all") +#环境因子之间的相关性热图
  geom_square() +#相关性显示形式
  geom_mark(r = NA,sig.thres = 0.05, size = 3.5, colour = "black")+#显著性标签
  scale_fill_gradient2( high = '#4CB847', mid = 'white',low = '#204E90') + #颜色设置
  anno_link(df_mantel, aes(color = df_p,
                           size = df_r))+
  scale_size_manual(values = c(0.5, 1, 1.5))+#连线粗细设置
  scale_color_manual(values = c("#ff5454","#ffa459","#59b2ff"))+#线条颜色设置
  guides(fill = guide_colorbar(title = "correlation", order = 1),#图例相关设置
         size = guide_legend(title = "Mantel's r",order = 2),
         color = guide_legend(title = "Mantel's p", order = 3),
         linetype = "none")+
  theme(legend.margin = margin(0,0,0,0))
ggsave(file.path(fig_out_dir, "cor.pdf"),units = "cm",width = 20,height = 18)

