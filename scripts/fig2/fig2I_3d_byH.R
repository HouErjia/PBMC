library(scatterplot3d)
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2_add_Hou"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
# 读取数据，确保正确读取
data <- read.csv("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou/CD93/CDS_ENSSSCG00000013593_coverage_avg.csv", header = TRUE)

# 检查数据结构和列名
cat("数据维度:", dim(data), "\n")
cat("列名:", names(data), "\n")
cat("前几行数据:\n")
print(head(data))

# 动态计算样本数量
num_samples <- ncol(data) - 1  # 减去Position列
cat("样本数量:", num_samples, "\n")

# 设置颜色 - 根据样本名称
sample_names <- names(data)[-1]
sample_colors <- ifelse(grepl("^W_", sample_names), "#CA3F4C", "#485998")
cat("样本名称:", sample_names, "\n")
cat("样本颜色:", sample_colors, "\n")

# 构建数据向量 - 确保使用正确的列名
x <- rep(data$Position, num_samples)  # 注意：Position首字母大写
y <- rep(1:num_samples, each = nrow(data))
z <- as.numeric(unlist(data[, -1]))  # 确保转换为数值型
colors <- rep(sample_colors, each = nrow(data))

# 验证数据长度
cat("数据长度验证:\n")
cat("x (位置) 长度:", length(x), "\n")
cat("y (样本) 长度:", length(y), "\n")
cat("z (读数) 长度:", length(z), "\n")
cat("colors 长度:", length(colors), "\n")

# 检查是否有NA值
cat("NA值检查:\n")
cat("x中的NA:", sum(is.na(x)), "\n")
cat("y中的NA:", sum(is.na(y)), "\n")
cat("z中的NA:", sum(is.na(z)), "\n")

# 如果长度不一致，停止执行
if(length(x) != length(y) | length(y) != length(z) | length(z) != length(colors)) {
  stop("错误：数据向量长度不一致！")
}

# 绘制3D柱状图
pdf(file.path(fig_out_dir, "CD93_corrected.pdf"), width = 12/2.54, height = 10/2.54, pointsize = 8)

s3d <- scatterplot3d(
  x = x, y = y, z = z, color = colors,
  type = "h", lwd = 2, pch = NA,
  xlab = "Genomic Position", ylab = "", zlab = "Reads Count",
  main = "CD93 Coverage", grid = FALSE, box = FALSE,
  y.ticklabs = sample_names,
  ylim = c(0.5, num_samples + 0.5),
  axis = TRUE)

# 添加图例
legend("topright", 
       legend = c("W Samples", "H Samples"),
       fill = c("#CA3F4C", "#485998"), 
       title = "Sample Groups",
       bty = "n")

dev.off()

cat("3D图形已成功保存为 CD93_3d.pdf\n")


