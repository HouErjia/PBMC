rm(list = ls())
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig1")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig1"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
library(limma)
library(edgeR)
rpkm_data<- read.csv("./fig1c_rpkm_values.csv", header = TRUE, row.names = 1)
log_rpkm <- log2(rpkm_data + 1)
group <- c(rep("H", 5), rep("W", 5))
design <- model.matrix(~0+factor(group))
colnames(design) <- levels(factor(group))

# 使用voom函数将RPKM转换为权重
v <- voom(log_rpkm, design, plot=TRUE)

# 拟合线性模型
fit <- lmFit(v, design)
contrast_matrix <- makeContrasts(H_vs_W = W - H, levels=design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, number=Inf, adjust="fdr")
head(results)
write.table(results,"./fig1c_PBMC_rpkm_diff.csv", sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE)

results<- read.csv("./fig1c_PBMC_rpkm_diff.csv", header = TRUE, row.names = 1)
library(ggplot2)
mytheme<-theme_bw()+theme(panel.border=element_rect(size=0.4,colour='black'),
                          panel.grid=element_blank(),axis.title=element_text(size=9,colour='black'),
                          axis.text=element_text(size=9,colour='black'),
                          axis.ticks=element_line(size=0.4,colour='black'),axis.ticks.length=unit(0.1,'cm'))
results$regulation <- ifelse(results$logFC > 0 & results$adj.P.Val < 0.05, "up",
                          ifelse(results$logFC < 0 & results$adj.P.Val < 0.05, "down", "NS"))
results$regulation <- factor(results$regulation, levels = c("NS", "down", "up"))# 将 regulation 转换为因子类型以用于颜色映射
table(results$regulation)
p1 <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val + 1e-299), colour = regulation)) +  # 添加一个小常数避免负无穷大
  geom_point(size = 0.4, alpha = 0.6,shape=16) +  # 设置点大小和透明度
  scale_colour_manual(values = c("NS" = "grey", "down" = "#4B5DA0", "up" = "#BF3F4A")) +
  mytheme +  
  theme(legend.title = element_blank(), strip.background = element_blank()) +
  xlab("log2 Fold change") +
  ylab("-log10 FDR")+
  xlim(-4, 4)
print(p1)
ggsave(file.path(fig_out_dir, "Volcano_PBMC_rpkm_diff.pdf"),units = "in",width=7/2.54,height=4/2.54)


