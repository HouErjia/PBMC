# ============================================================
# Translation Efficiency Heatmap of Immune-related Genes
# 按行 scale：每个基因内部比较 H 和 W 的 TE 高低
# 基因名逆时针旋转 90°
# ============================================================

rm(list = ls())


setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2_add_Hou"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
library(pheatmap)
library(gplots)
library(ggplot2)

# =========================
# 1. 绘图主题
# =========================
mytheme <- theme_bw() + 
  theme(
    panel.border = element_rect(linewidth = 0.4, colour = "black"),
    panel.grid = element_blank(),
    axis.title = element_text(size = 10, colour = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.length = unit(0.1, "cm")
  )

# =========================
# 2. 读取数据
# =========================
data <- read.csv(
  "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou/fig2f_selected_genes_data_byH.csv",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

cat("原始数据：\n")
print(head(data))

cat("\n数据结构：\n")
print(str(data))

# =========================
# 3. 只保留 TE 数值列
# =========================
# 第一列 X.1 是 Ensembl ID，不参与绘图
data_numeric <- data[, -1]

# 确保 H 和 W 是数值型
data_numeric$H <- as.numeric(data_numeric$H)
data_numeric$W <- as.numeric(data_numeric$W)

cat("\n用于绘图的原始 TE 矩阵：\n")
print(data_numeric)

# =========================
# 4. 按行 scale / 按基因 Z-score
# =========================
# 原始矩阵结构：
# 行 = Gene
# 列 = H / W
#
# scale(data_numeric) 是按列 scale：
#   比较 H 组内部不同基因之间的高低，
#   比较 W 组内部不同基因之间的高低。
#
# t(scale(t(data_numeric))) 是按行 scale：
#   比较同一个基因内部 H 和 W 的高低。
#
# 这一步适合展示：
#   同一个基因中，W 组 TE 是否高于 H 组。

data_row_zscore <- t(scale(t(data_numeric)))

# 转换为 data.frame，保留基因名和分组名
data_row_zscore <- as.data.frame(data_row_zscore)

cat("\n按行 scale 后的 TE Z-score 矩阵：\n")
print(data_row_zscore)

# =========================
# 5. 转置用于画图
# =========================
# pheatmap 中：
# 行 = H / W
# 列 = Gene
data_t <- t(data_row_zscore)

cat("\n最终用于热图的按行 scale 矩阵：\n")
print(data_t)

# =========================
# 6. 设置颜色
# =========================
# 蓝色 = 同一基因内部 TE 较低
# 白色 = 中间值
# 红色 = 同一基因内部 TE 较高
heatmap_colors <- colorRampPalette(
  c(
    "#002256",
    "#175391",
    "#96C7DF",
    "#D1E5F0",
    "white",
    "#F4A582",
    "#D45C4A",
    "#B2182B",
    "#67001F"
  )
)(100)

# =========================
# 7. 输出 PDF 热图
# =========================
pdf(
  "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou/fig2h_corrected_byH_rowZscore.pdf",
  width = 15 / 2.54,
  height = 10 / 2.54
)

pheatmap(
  data_t,
  color = heatmap_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  cellwidth = 15,
  cellheight = 15,
  angle_col = 90,
  main = "Translation Efficiency Heatmap of Immune-related Genes",
  border_color = "white",
  gaps_row = NULL,
  gaps_col = NULL
)

dev.off()

cat("\n✅ 按行 scale 热图已输出：D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou/fig2h_corrected_byH_rowZscore.pdf\n")


