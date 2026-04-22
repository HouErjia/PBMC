#### 去除批次效应
# 安装 install.packages("devtools")
rm(list=ls())
# devtools::install_github("zhangyuqing/sva-devel")
library(sva)
# 导入counts文件
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig1/RNAseq数据前处理")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig1/RNAseq数据前处理"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
count_data <- read.csv("./correlation_readcounts_PBMC.csv",header=T,row.names=1)
count_data <- as.matrix(count_data)
# 导入样本信息
metadata <- read.csv("./PBMC_group.csv")

## 去除批次效应
count_data_ComBat <- ComBat_seq(count_data, batch=metadata$group1, group = metadata$group)
write.csv(count_data_ComBat, file = "./combat_PBMC_readcounts.csv", row.names = TRUE)

library(edgeR)
group<-metadata[,2]
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(count_data_ComBat)
contrast.matrix<-makeContrasts("HvsW",levels=design)

fit <- lmFit(count_data_ComBat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, coef=1, n=Inf)
write.csv(results, file = "./differential_expression_results_HvsW2.csv")

group=data.frame(c(rep('H', 5),rep('W', 5)))
rownames(group) <- colnames(count_data_ComBat)
library(pheatmap)
significant_genes <- rownames(results)[results$adj.P.Val < 0.05 & abs(results$logFC) > 2 & results$B > 0]
cat("Number of significant genes:", length(significant_genes), "\n")
pdf(file.path(fig_out_dir, "PBMC_combat.pdf"),width=25/2.54,height = 28/2.54)
pheatmap_result <- pheatmap(edata$E[significant_genes, ], annotation_col = group, show_rownames = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_cols = FALSE)
dev.off()

row_clusters <- cutree(pheatmap_result$tree_row, k=3)
output <- data.frame(gene = significant_genes, row_clusters)
write.table(output, file="aovRNA_cluster.txt", quote = FALSE, sep="\t", row.names = FALSE)
cluster_files <- c("aovRNA_cluster1.txt", "aovRNA_cluster2.txt", "aovRNA_cluster3.txt", "aovRNA_cluster4.txt")
for (i in 1:3) {
  cluster_data <- subset(output, row_clusters == i)
  write.table(cluster_data, file=cluster_files[i], quote = FALSE, sep="\t", row.names = FALSE)
}
#绘制火山图
library(grid)
library(ggplot2)
mytheme<-theme_bw()+theme(panel.border=element_rect(size=0.4,colour='black'),
                          panel.grid=element_blank(),axis.title=element_text(size=9,colour='black'),
                          axis.text=element_text(size=9,colour='black'),
                          axis.ticks=element_line(size=0.4,colour='black'),
                          axis.ticks.length=unit(0.1,'cm'))

output_data <- read.csv("./differential_expression_results_HvsW.csv", row.names = 1)
output_data$flag <- ifelse(output_data$adj.P.Val < 0.05,
                           ifelse(output_data$logFC > 0, "Up", "Down"),
                           "NS")

volcano_plot <- ggplot(output_data, aes(x = logFC, y = -log10(adj.P.Val),, col = flag)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c('Up' = 'red', 'Down' = 'lightblue', 'NS' = 'grey')) +
  mytheme +
  theme(legend.title = element_blank(), strip.background = element_blank()) +
       xlab("Log2 Fold Change (logFC)") +
       ylab ("-log10(Adjusted P-value)")

print(volcano_plot)


