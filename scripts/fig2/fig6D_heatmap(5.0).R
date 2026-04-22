# =========================
# 手动添加未映射菌名的完整代码
# =========================

# ---- 环境清理与加载包 ----
rm(list = ls())

setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2_add_Hou"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
library(tidyverse)
library(readxl)
library(pheatmap)
library(RColorBrewer)

# ---- 文件读取 ----
data <- read_excel("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou/fig6D_Lys(new).xlsx", 
                   sheet = "Sheet1", col_names = FALSE)

# ---- 提取分组信息（第1行）----
full_group_row <- as.character(data[1, ])
group_info <- c()
current_group <- ""
for (i in 2:length(full_group_row)) {
  if (!is.na(full_group_row[i]) && full_group_row[i] != "") {
    current_group <- full_group_row[i]
  }
  group_info <- c(group_info, current_group)
}

# ---- 提取样本名称（第2行）----
sample_names <- as.character(data[2, 2:ncol(data)])
sample_names <- sample_names[!is.na(sample_names) & sample_names != ""]
group_info <- group_info[1:length(sample_names)]

# ---- 提取表达矩阵 ----
expression_data <- data[3:nrow(data), 2:ncol(data)]
gene_names <- data[3:nrow(data), 1] %>% pull()
expression_df <- as.data.frame(expression_data)
rownames(expression_df) <- gene_names
colnames(expression_df) <- sample_names
expression_matrix <- as.matrix(expression_df)
mode(expression_matrix) <- "numeric"

# ---- 按分组顺序排列列 ----
group_order <- c("H PBMC", "W PBMC", "Healthy", "Weaning")
ordered_samples <- sample_names[order(factor(group_info, levels = group_order))]
expression_matrix <- expression_matrix[, ordered_samples]
group_info_ordered <- group_info[order(factor(group_info, levels = group_order))]

# ---- 创建完整的菌名映射（属名/种名斜体显示）----
complete_bacteria_mapping <- c(
  # 使用 sp. (单数)
  "ASV8" = "ASV8 (<i>Coprococcus</i> sp.)",
  "ASV1555" = "ASV1555 (<i>Enterococcus</i> sp.)",
  "ASV54" = "ASV54 (<i>Clostridium</i> sp.)",
  "ASV30" = "ASV30 (<i>Clostridium</i> sp.)",
  "ASV73" = "ASV73 (<i>Turicibacter</i> sp.)",
  "ASV9" = "ASV9 (<i>Romboutsia</i> sp.)",
  "ASV17" = "ASV17 (<i>Sharpea</i> sp.)",
  "ASV1834" = "ASV1834 (<i>Lachnospiraceae</i> sp.)",
  "ASV22" = "ASV22 (<i>Lactobacillus</i> sp.)",
  "ASV3402" = "ASV3402 (<i>Megasphaera</i> sp.)",
  "ASV1540" = "ASV1540 (<i>Megasphaera</i> sp.)",
  "ASV77" = "ASV77 (<i>Lactobacillus</i> sp.)",
  "ASV1853" = "ASV1853 (<i>Lactobacillus</i> sp.)",
  "ASV38" = "ASV38 (<i>Blautia</i> sp.)",
  "ASV48" = "ASV48 (<i>Blautia</i> sp.)",
  "ASV1515" = "ASV1515 (<i>Streptococcus</i> sp.)",
  "ASV2996" = "ASV2996 (<i>Fusobacterium</i> sp.)",
  "ASV2874" = "ASV2874 (<i>Fusobacterium</i> sp.)",
  "ASV167" = "ASV167 (<i>Collinsella</i> sp.)",
  "ASV2" = "ASV2 (<i>Lactobacillus</i> sp.)",
  "ASV318" = "ASV318 (<i>Peptostreptococcus</i> sp.)",
  "ASV2586" = "ASV2586 (<i>Campylobacter</i> sp.)",
  "ASV2860" = "ASV2860 (<i>Phascolarctobacterium</i> sp.)",
  "ASV172" = "ASV172 (<i>Mogibacterium</i> sp.)",
  "ASV3058" = "ASV3058 (<i>Allisonella</i> sp.)",
  "ASV1833" = "ASV1833 (<i>Enterococcus</i> sp.)",
  "ASV7" = "ASV7 (<i>Lactobacillus</i> sp.)",
  "ASV7036" = "ASV7036 (<i>Bacteroides</i> sp.)",
  "ASV5" = "ASV5 (<i>UCG-002</i> sp.)",
  
  # 使用 spp. (复数)
  "ASV45" = "ASV45 (<i>Escherichia–Shigella</i> spp.)",
  "ASV1517" = "ASV1517 (<i>Escherichia–Shigella</i> spp.)",
  
  # 具体种名
  "ASV1831" = "ASV1831 (<i>Eubacterium coprostanoligenes</i>)",
  "ASV21" = "ASV21 (<i>Clostridium scindens</i>)",
  "ASV1089" = "ASV1089 (<i>Fusobacterium gastrosuis</i>)",
  "ASV2839" = "ASV2839 (<i>Pasteurella aerogenes</i>)",
  "ASV3729" = "ASV3729 (<i>Cloacibacillus porcorum</i>)",
  
  # 特殊分类
  "ASV31" = "ASV31 (<i>Lactobacillus vaginalis</i> Bacillus)",
  "ASV1527" = "ASV1527 (<i>Lactobacillus vaginalis</i> Bacillus)",
  "ASV33" = "ASV33 (gut metagenome)",
  "ASV4" = "ASV4 (<i>Christensenellaceae</i> sp.)",
  "ASV28" = "ASV28 (<i>Christensenellaceae</i> sp.)",
  "ASV3787" = "ASV3787 (<i>UCG.002</i>)",
  "ASV55" = "ASV55 (<i>Anaerotignum lactatifermentans</i>)",
  
  # PBMC样本
  "H6_PBMC" = "H6_PBMC",
  "H7_PBMC" = "H7_PBMC",
  "W6_PBMC" = "W6_PBMC",
  "W10_PBMC" = "W10_PBMC"
)


# 创建菌名映射向量
bacteria_names <- complete_bacteria_mapping

# 检查是否还有缺失的样本
missing_samples <- setdiff(ordered_samples, names(bacteria_names))
if(length(missing_samples) > 0) {
  cat("警告: 以下样本在菌名映射中找不到，将使用原始名称:\n")
  cat(missing_samples, "\n")
  # 为缺失的样本添加默认映射
  for(sample in missing_samples) {
    bacteria_names[sample] <- sample
  }
} else {
  cat("所有样本都已成功映射!\n")
}

# 创建显示用的菌名标签
display_names <- bacteria_names[ordered_samples]

# ---- 创建行注释（样本分组）----
annotation_row <- data.frame(
  Group = factor(group_info_ordered, levels = group_order),
  row.names = ordered_samples
)

# ---- 分组颜色 ----
ann_colors <- list(
  Group = c(
    "H PBMC" = "#E64B35",   # 红色
    "W PBMC" = "#3182BD",   # 蓝色
    "Healthy" = "#4DAF4A",  # 绿色
    "Weaning" = "#FF7F00"   # 橙色
  )
)

# ---- 使用RdBu配色方案 ----
heatmap_colors <- colorRampPalette(brewer.pal(11, "RdBu"))(100)

# ---- 绘制临时热图以获取聚类信息 ----
temp_heatmap <- pheatmap(
  t(expression_matrix),
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  color = heatmap_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  fontsize = 7,
  fontsize_row = 7,
  fontsize_col = 7,
  treeheight_col = 15,
  main = "Glu",
  cellwidth = 10,
  cellheight = 8,
  silent = TRUE
)

# ---- 提取聚类结果并添加cluster注释 ----
col_clusters <- cutree(temp_heatmap$tree_col, k = 5)
cluster_annotation <- data.frame(
  Cluster = factor(col_clusters),
  row.names = names(col_clusters)
)

# 更新注释颜色，添加cluster颜色
ann_colors$Cluster <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")[1:length(unique(col_clusters))]
names(ann_colors$Cluster) <- as.character(1:length(unique(col_clusters)))

# ---- 绘制最终热图并输出PDF ----
pdf("fig6D_Val_transposed_heatmap_with_complete_bacteria_names.pdf", 
    width = 15/2.54, height = 19/2.54)  # 增加宽度以容纳更长的菌名

pheatmap(
  t(expression_matrix),
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  annotation_row = annotation_row,   # 样本分组注释
  annotation_col = cluster_annotation, # 基因聚类注释
  annotation_colors = ann_colors,
  color = heatmap_colors,
  breaks = seq(-2, 2, length.out = 101),
  labels_row = display_names,  # 使用完整菌名映射作为行标签
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  fontsize = 7,
  fontsize_row = 7,  # 可能需要调整字体大小
  fontsize_col = 7,
  treeheight_col = 15,
  main = "Lys",
  cellwidth = 10,
  cellheight = 8,
  margin = c(3, 5)  # 调整边距，确保菌名显示完整
)

dev.off()

# ---- 保存更新后的完整菌名映射 ----
updated_mapping <- data.frame(
  sample_id = names(bacteria_names),
  display_name = bacteria_names,
  stringsAsFactors = FALSE
)

write.csv(updated_mapping, "fig6D_updated_fig6D_fig6D_updated_complete_bacteria_mapping.csv", row.names = FALSE)

# ---- 输出统计信息 ----
cat("=== 数据统计 ===\n")
cat("密码子分类数量:", nrow(expression_matrix), "\n")
cat("样本数量:", ncol(expression_matrix), "\n")
cat("\n分组统计:\n")
print(table(group_info_ordered))
cat("\n聚类统计:\n")
print(table(col_clusters))
cat("\n分组顺序:", paste(group_order, collapse = " → "), "\n")
cat("\n热图已保存为: fig6D_Glu(3.0)_transposed_heatmap_with_complete_bacteria_names.pdf\n")
cat("更新后的菌名映射已保存为: updated_fig6D_fig6D_updated_complete_bacteria_mapping.csv\n")

# ---- 输出菌名映射信息 ----
cat("\n=== 菌名映射 ===\n")
cat("使用的菌名映射: 完整手动映射\n")
cat("总映射数量:", length(bacteria_names), "\n")
cat("当前样本使用的映射:\n")
for(i in 1:length(display_names)) {
  cat(ordered_samples[i], " → ", display_names[i], "\n")
}

# ---- 检查数据一致性 ----
cat("\n=== 数据一致性检查 ===\n")
cat("表达矩阵样本数:", ncol(expression_matrix), "\n")
cat("菌名映射样本数:", length(display_names), "\n")
cat("样本是否匹配:", identical(colnames(expression_matrix), names(display_names)), "\n")


