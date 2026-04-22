# devtools::install_github("Github-Yilei/ggcor")
library(vegan)
library(dplyr)
library(ggcor)
library(ggplot2)
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2_add_Hou"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
# 1. 读取数据
ENV_all <- read_excel("fig4F_ASV_相对丰度.xlsx", sheet = "相对丰度_转置")
# 选取第1列（Sample）和第2~11列、第32~41列
# 注意：R 中列索引从1开始，第2~11列对应索引 2:11，第32~41列对应索引 32:41
ENV_all <- ENV_all  %>%select(1, 2:11, 32:41)
ENV_all <- as.data.frame(ENV_all)
rownames(ENV_all) <- ENV_all$Sample
ENV_all <- ENV_all[, -1]

BA <- read.csv("fig4F_selected_genes_TE_transposed.csv", row.names = 1, check.names = FALSE)
BA[is.na(BA)] <- 0

# 2. 筛选共同样本
common_samples <- intersect(rownames(BA), rownames(ENV_all))
cat("共同样本:", common_samples, "\n")

# 3. 在共同样本中过滤全0列
ENV_common <- ENV_all[common_samples, ]
# 只保留在共同样本中非全0的列
non_zero_cols_common <- which(colSums(ENV_common) > 0)
ENV <- ENV_common[, non_zero_cols_common]

BA <- BA[common_samples, ]

cat("过滤后ENV数据维度:", dim(ENV), "\n")
cat("过滤后BA数据维度:", dim(BA), "\n")

# 4. 检查数据质量
cat("ENV中零值比例:", sum(ENV == 0) / (nrow(ENV) * ncol(ENV)), "\n")
cat("BA中零值比例:", sum(BA == 0) / (nrow(BA) * ncol(BA)), "\n")

spec=list()
for(i in 1:length(BA)){
  spec[[names(BA)[i]]]=i
}
library(vegan)
ENV_hell <- decostand(ENV, method = "hellinger")

df_mantel <- mantel_test(BA, ENV_hell,
                         mantel.fun = 'mantel',
                         spec.dist.method = 'euclidean',
                         env.dist.method = 'euclidean',
                         spec.select = spec)
# df_mantel <- mantel_test(BA, ENV,
#                          mantel.fun = 'mantel',
#                          spec.dist.method = 'bray',
#                          env.dist.method = 'euclidean',
#                          spec.select = spec
#                          # spec.select = list(`Primary BAs` = c(1,2,4),
                         #                    `Secondary BAs` = c(3,5:8)
                                            # )
                         # )#将群落数据按组进行分开
df_mantel <- df_mantel %>%
  mutate(df_r = cut(r, breaks = c(-Inf, 0.25,0.5, Inf),
                    labels = c("< 0.25", "0.25 - 0.5",  ">= 0.5")),#定义Mantel的R值范围标签
         df_p = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("<=0.01", "0.01 - 0.05", ">= 0.05")))#定义Mantel的P值范围标签


quickcor(ENV,method = "spearman", type = "upper", cor.test = T, cluster.type = "all") +#环境因子之间的相关性热图
  geom_tile() +#相关性显示形式
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
  theme(legend.margin = margin(0,0,0,0),
        axis.text.x.top = element_blank())
ggsave(file.path(fig_out_dir, "cor.pdf"),units = "cm",width = 20,height = 20)

