rm(list=ls())
library(tidyverse)
library(pheatmap) 
library(RColorBrewer)
mytheme_noborder <- theme_bw() + 
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = 9, colour = 'black'),
    axis.text = element_text(size = 9, colour = 'black'),
    axis.ticks = element_line(size = 0.4, colour = 'black'),
    axis.ticks.length = unit(0.1, 'cm'),
    axis.line.x = element_line(size = 0.4, colour = 'black'),
    axis.line.y = element_line(size = 0.4, colour = 'black')
  )

mytheme <- theme_bw() + 
  theme(
    panel.border = element_rect(size = 0.4, colour = 'black'),
    panel.grid = element_blank(),
    axis.title = element_text(size = 10, colour = 'black'),
    axis.text = element_text(size = 10, colour = 'black'),
    axis.ticks = element_line(size = 0.4, colour = 'black'),
    axis.ticks.length = unit(0.1, 'cm')
  )
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig5")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig5"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
input_file <- "fig5c_average_heatmap.csv"  
df <- read_csv(input_file)  

df_sorted <- df %>% 
  mutate(
    third_base = substr(codon, 3, 3),  # 提取第3位碱基
    category = ifelse(third_base %in% c("A", "T"), "AU3", "GC3"),  # 分类
    sort_order = ifelse(category == "AU3", 1, 2)  # AU3 排第1组，GC3 排第2组
  ) %>% 
  arrange(sort_order, codon)  # 先按类别排序，再按密码子排序

# 整理成热图矩阵（行=密码子，列=分组）
mat <- df_sorted %>% 
  column_to_rownames(var = "codon") %>%  # 密码子作为行名
  select(-c(third_base, category, sort_order)) %>%  # 移除辅助列
  as.matrix()  # 转为矩阵

# ----------------------
# 3. 行注释：标记 AU3/GC3
# ----------------------
anno_row <- df_sorted %>% 
  distinct(codon, category) %>%  # 去重，保留密码子-类别对应关系
  column_to_rownames(var = "codon")  # 行名=密码子

# ----------------------
# 4. 绘制热图（按 AU3/GC3 分组展示）
# ----------------------
# 自定义颜色（红蓝渐变）
color_palette <- colorRampPalette(brewer.pal(11, "RdBu"))(100)

# 绘制热图：行按 AU3→GC3 排序，列保持分组顺序
pdf(file.path(fig_out_dir, "CE_average_heatmap.pdf"), width =5.5/2.54, height = 19/2.54)
pheatmap(
  mat,
  color = color_palette,
  annotation_row = anno_row,  # 行注释（AU3/GC3）
  annotation_colors = list(
    category = c(AU3 ="#3182BD" , GC3 = "#E64B35")  # 类别颜色
  ),
  show_rownames = TRUE,  # 显示密码子行名
  show_colnames = TRUE,  # 显示分组列名
  cluster_rows = FALSE,  # 关闭行聚类（保持手动排序）
  cluster_cols = FALSE,  # 关闭列聚类（保持原顺序）
  gaps_row = sum(df_sorted$category == "AU3"),  # 在 AU3/GC3 之间加分隔线
  main = "Z-score Heatmap (AU3/GC3 Grouped)")
dev.off()


