rm(list=ls())
library("RColorBrewer")
library("xlsx")
library(pheatmap)
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig3")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig3"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
data<-read.xlsx("./fig3_correlation_PBMC_diff.xlsx",1,header=T,row.names=1)
data=as.data.frame(data)

group <- data.frame(g = c(rep("down_H", 2), rep("down_W", 2), rep("up_H", 2), rep("up_W", 2)))
rownames(group) <- colnames(data)  # 确保行名与 data 的列名一致

annotation_colors <- list(g = c( down_H = "#4B5DA0", down_W = "#B53F48",up_H = "#4B5DA0", up_W = "#B53F48" ))

pdf(file.path(fig_out_dir, "PBMC_diff_column.pdf"), width = 24/2.54, height = 45/2.54)
pheatmap(
  data,
  scale = "column",
  annotation_col = group,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  show_colnames = TRUE,cluster_cols = TRUE,
  fontsize = 10,
  angle_row = 180,     
  cellwidth = 10,  
  cellheight = 10)
dev.off()

cor_matrix<-read.xlsx("./fig3_correlation_PBMC_diff.xlsx",1,header=T,row.names=1)
selected_codons <- c("CTT","CTC","CTA","CTG","TTA","TTG")
codon_selected_cols <- grep(paste(selected_codons, collapse = "|"),
                            colnames(cor_matrix),
                            ignore.case = TRUE)
selected_cor_matrix <- cor_matrix[, codon_selected_cols]
# 绘制特定密码子热图
library(pheatmap)
p1 <- pheatmap(
  selected_cor_matrix,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  show_rownames = TRUE,show_colnames = TRUE, cellwidth = 20, cellheight = 20 )


