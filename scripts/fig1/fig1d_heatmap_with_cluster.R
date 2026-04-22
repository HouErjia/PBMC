rm(list = ls())
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig1")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig1"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
library(pheatmap)
library(DESeq2)

results <- read.csv("./fig1c_PBMC_rpkm_diff.csv", header = TRUE, row.names = 1)
sig_genes <- subset(results, results$adj.P.Val < 0.05 & abs(logFC) != 0)
expr_matrix <- read.csv("./fig1c_rpkm_values.csv", row.names = 1)
expr_matrix <- expr_matrix[rownames(sig_genes), ]
group_info <- data.frame(group = factor(rep(c("H", "W"), each = 5)), row.names = colnames(expr_matrix))

p1 <- pheatmap(as.matrix(expr_matrix), scale = "row", treeheight_row = 100,
              show_rownames = FALSE, show_colnames = TRUE, 
              annotation_col = group_info, color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
              cluster_rows = TRUE, cluster_cols = TRUE, 
              clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete",
              draw=FALSE) 

clusters <- cutree(as.hclust(p1$tree_row), k=5) 
sig_genes$cluster <- clusters
cluster_annotation <- data.frame(cluster = factor(sig_genes$cluster), row.names = rownames(sig_genes))

p2 <-pheatmap(as.matrix(expr_matrix), scale = "row", treeheight_row = 100,
         show_rownames = FALSE, show_colnames = TRUE, 
         annotation_col = group_info, annotation_row = cluster_annotation, # 新增：添加cluster标识作为行注释
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_rows = TRUE, cluster_cols = TRUE, 
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete")



