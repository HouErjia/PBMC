rm(list=ls())
library(ggplot2)

#PBMC柱状图
rm(list=ls())
library(ggplot2)
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig1")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig1"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
file <- "./fig1d_GO_clusters.xlsx"
library(xlsx)
cluster1=read.xlsx(file,1)
cluster2=read.xlsx(file,2)
cluster3=read.xlsx(file,3)
cluster4=read.xlsx(file,4)
cluster5=read.xlsx(file,5)
cluster6=read.xlsx(file,6)
#以padj为横坐标
p1 <- ggplot(cluster1, aes(x = -log2(padj), y = reorder(Term, -log2(padj)))) +
  geom_bar(stat = "identity", fill = "gray", color = "black", linewidth = 0.8, width = 0.9) + 
  labs(title = "GO Pathways by p.adj",
       x = "p.adj (-log2 scale)",
       y = "GO Pathway") +
  theme_minimal() +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0)), limits = c(0, max(-log2(cluster1$padj)) + 0.1)) +  
  scale_y_discrete(expand = expansion(add = c(0.1, 0.1))) +  
  theme(axis.text.y = element_text(size = 15, hjust = 0.9),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.line.x = element_line(color = "black", linewidth = 1),
        axis.ticks.x = element_line(color = "black", linewidth = 1),
        axis.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())
#以gene ratio为横坐标
p1 <- ggplot(cluster4, aes(x = gene_ratio, y = reorder(Term, gene_ratio))) + 
  geom_bar(stat = "identity", fill = "gray", color = "black", linewidth = 0.8, width = 0.9) + 
  labs(title = "GO Pathways by Gene Ratio",
       x = "Gene Ratio", 
       y = "GO Pathway") +
  theme_minimal() +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0)), limits = c(0, max(cluster4$gene_ratio) + 0.01)) + 
  scale_y_discrete(expand = expansion(add = c(0.1, 0.1))) +  
  theme(axis.text.y = element_text(size = 15, hjust = 0.9),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.line.x = element_line(color = "black", linewidth = 1),
        axis.ticks.x = element_line(color = "black", linewidth = 1),
        axis.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())
p1
ggsave(file.path(fig_out_dir, "cluster4.pdf"),p1,units = "in",width=20/2.54,height = 8/2.54)


