rm(list=ls())
library(pheatmap)
library(gplots)
library(ggplot2)
library(xlsx)
mytheme<-theme_bw()+theme(panel.border=element_rect(size=0.4,colour='black'),
                          panel.grid=element_blank(),axis.title=element_text(size=10,colour='black'),
                          axis.text=element_text(size=10,colour='black'),
                          axis.ticks=element_line(size=0.4,colour='black'),axis.ticks.length=unit(0.1,'cm'))
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
data<-read.xlsx("./fig2f.xlsx",4,header=T,row.names=1)
data_t <- t(data)
data_log <- log2(data_t + 1)
pdf(file.path(fig_out_dir, "ribo_immue_related_genes_heatmap.pdf"),width=10/2.54,height = 5/2.54)
pheatmap(data_log,
  color = colorRampPalette(c( "#002256","#175391","#96C7DF", "#D1E5F0","white","#F4A582","#D45C4A", "#B2182B","#67001F"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE,  show_rownames = TRUE,
  show_colnames = TRUE, fontsize_row = 8, 
  fontsize_col = 8, cellwidth = 15, cellheight = 15, # 单元格宽度和高度
  main = "Ribo-seq Heatmap of Immune-related Genes ",
  border_color = "white", cell_spacing = c(0, 0), gaps_row = NULL,  gaps_col = NULL )
dev.off()


