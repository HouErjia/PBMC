
# 读取两个表
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2_add_Hou"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
df1 <- read.csv("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou/CD93/CSV_H6_W6/CDS_ENSSSCG00000007116_coverage.csv", header = TRUE)
df2 <- read.csv("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou/CD93/CSV_H7_W10/CDS_ENSSSCG00000007116_coverage.csv", header = TRUE)

# 保留位置列
pos <- df1[, 1]

# 找到共有的数值列名（去掉位置列）
common_cols <- intersect(colnames(df1)[-1], colnames(df2)[-1])

# 取出数值部分
val1 <- df1[, common_cols]
val2 <- df2[, common_cols]

# 对应列取平均
avg_vals <- (val1 + val2) / 2

# 合并位置列与平均结果
merged_avg <- cbind(Position = pos, avg_vals)

# 保存为新 CSV 文件
write.csv(merged_avg, "CDS_ENSSSCG00000013593_coverage_avg.csv", row.names = FALSE)

# 查看前几行结果
head(merged_avg)



