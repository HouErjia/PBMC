library(scatterplot3d)
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
data <- read.csv("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2\fig2h_CCR5.csv", header = TRUE)

#修改
sample_colors <- ifelse(grepl("^W_", names(data)[-1]), "#CA3F4C", "#485998")
colors <- rep(sample_colors, each = nrow(data))

# 构建坐标轴数据
x <- rep(data$position, 4)  # X轴：基因组位置
y <- rep(1:4, each = nrow(data))  # Y轴：样本编号
z <- unlist(data[, -1])      # Z轴：reads数量（排除position列）

# 绘制3D柱状图
library(scatterplot3d)
pdf(file.path(fig_out_dir, "CCR5.pdf"), width = 12/2.54, height = 10/2.54, pointsize = 8)
s3d <- scatterplot3d(
  x = x, y = y, z = z,color = colors,
  type = "h", lwd = 3,pch = NA,
  xlab = "Genomic Position",ylab = "",zlab = "Reads Count",
  main = "",grid = FALSE,box = FALSE,
  y.ticklabs = c("W_RNAseq", "H_RNAseq", "H_Riboseq", "W_Riboseq"),  # 修正标签
  ylim = c(0.5, 4.5),
  axis = TRUE)
dev.off()

# 添加图例（按颜色分组）
legend("topright", 
       legend = c("W Samples (RNAseq & Riboseq)", "H Samples (RNAseq & Riboseq)"),
       fill = c("#CA3F4C", "#485998"), 
       title = "Sample Groups")

