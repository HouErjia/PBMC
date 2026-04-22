rm(list=ls())
mytheme_noborder<-theme_bw()+theme(panel.border=element_blank(),
                                   panel.grid=element_blank(),
                                   axis.title=element_text(size=9,colour='black'),
                                   axis.text=element_text(size=9,colour='black'),
                                   axis.ticks=element_line(size=0.4,colour='black'),
                                   axis.ticks.length=unit(0.1,'cm'),
                                   axis.line.x=element_line(size=0.4,colour='black'),
                                   axis.line.y=element_line(size=0.4,colour='black'))
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig5")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig5"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
library(xlsx)
enc_data<-read.xlsx("./fig5f_results_right.xlsx",1,header=T,row.names=1)
Groups=data.frame(group=c(rep("H",20),rep("W",22)))
enc_data$Group <- Groups$group
enc_data$Species <- rownames(enc_data)
# 计算理论曲线
enc_data$Theoretical_ENC <- 2 + enc_data$GC3 + 29/(enc_data$GC3^2 + (1-enc_data$GC3)^2)

# 绘制ENC图
library(ggplot2)
library(ggrepel)

# x_min <- 0
# x_max <- max(enc_data$GC3, na.rm = TRUE)
# y_min <- 20
# y_max <- max(c(enc_data$ENC, enc_data$Theoretical_ENC), na.rm = TRUE)
pdf(file.path(fig_out_dir, "ENC_asv.pdf"), width = 8/2.54, height = 8/2.54)
p <- ggplot(enc_data, aes(x = GC3, y = ENC)) +
  geom_point(aes(color = Group), size = 4, alpha = 0.8, shape = 16, stroke = 0) + 
  geom_line(aes(y = Theoretical_ENC), linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = c("#4B5EA0", "#CA3F4C"), name = "Group") +
  geom_text_repel(aes(label = Species), size = 3, max.overlaps = 20) + 
  labs(
    subtitle = "42 Bacterial Species Analysis (20 Healthy, 22 Stress)",
    x = "GC3 Content",
    y = "ENC Value",
    caption = "Dashed line: Theoretical curve under neutral evolution"
  ) +
  mytheme_noborder +  # 替换为自定义主题
  theme(legend.position = "right") 
print(p)
dev.off()


library(FactoMineR)
library(ggplot2)

rscu_data <- readxl::read_excel("./RSCU.xlsx")

rownames(rscu_data) <- rscu_data$Codon
rscu_data$Codon <- NULL

# 转置数据（菌种为行，密码子为列）
species_rscu <- t(rscu_data)

# 执行对应分析
coa_result <- CA(species_rscu, graph = FALSE)

# 提取坐标
coa_coords <- as.data.frame(coa_result$row$coord[, 1:2])
colnames(coa_coords) <- c("Dim1", "Dim2")
coa_coords$Species <- rownames(coa_coords)
Groups <- data.frame( Species = colnames(rscu_data),   group = c(rep("H", 20), rep("W", 22)))
coa_coords <- merge(coa_coords, Groups, by = "Species")

# 计算方差解释比例
variance <- coa_result$eig[, 2]
dim1_var <- round(variance[1], 1)
dim2_var <- round(variance[2], 1)

# 绘制COA图
pdf(file.path(fig_out_dir, "ENC_COA.pdf"), width = 8/2.54, height = 8/2.54)
p2 <-ggplot(coa_coords, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = group), size = 4, alpha = 0.8, shape = 16, stroke = 0) + 
  stat_ellipse(aes(group = group, color = group), level = 0.95) +  
  scale_color_manual(values = c("H" = "#4B5EA0", "W" = "#CA3F4C")) + 
  labs(
    x = paste0("Dimension 1 (", dim1_var, "%)"),
    y = paste0("Dimension 2 (", dim2_var, "%)"), ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed")
print(p2)
dev.off()
library(corrplot)
library(gplots)

# 计算菌种间的RSCU相关性
library(FactoMineR)
library(corrplot)

rscu_data <- readxl::read_excel("./RSCU.xlsx")
rownames(rscu_data) <- rscu_data$Codon  # 设置行名为密码子
rscu_data$Codon <- NULL
species_rscu <- t(rscu_data)  # 转置后行=物种，列=密码子

# 检查转置后列名是否丢失（重要！）
if (is.null(colnames(species_rscu))) {
  colnames(species_rscu) <- rownames(rscu_data)
}

# 执行对应分析
coa_result <- CA(species_rscu, graph = FALSE)
coa_coords <- as.data.frame(coa_result$row$coord[, 1:2])
colnames(coa_coords) <- c("Dim1", "Dim2")
coa_coords$Species <- rownames(coa_coords)

# 添加分组信息
Groups <- data.frame(
  Species = rownames(species_rscu),  # 物种名在行
  group = c(rep("H", 20), rep("W", 22))
)
coa_coords <- merge(coa_coords, Groups, by = "Species")

# 创建相关系数结果表（关键修复！）
corr_results <- data.frame(
  Codon = colnames(species_rscu)  # 或 rownames(rscu_data)
)
stopifnot(nrow(corr_results) > 0) # 验证非空

# 计算与Dim1和Dim2的相关系数
corr_results$Corr_Dim1 <- apply(
  species_rscu, 2, 
  function(x) cor(x, coa_coords$Dim1, method = "pearson")
)
corr_results$Corr_Dim2 <- apply(
  species_rscu, 2, 
  function(x) cor(x, coa_coords$Dim2, method = "pearson"))

# 计算p值
corr_results$Pvalue_Dim1 <- apply(
  species_rscu, 2, 
  function(x) cor.test(x, coa_coords$Dim1, method = "pearson")$p.value
)

corr_results$Pvalue_Dim2 <- apply(
  species_rscu, 2, 
  function(x) cor.test(x, coa_coords$Dim2, method = "pearson")$p.value
)

# 添加显著性标记
corr_results$Sig_Dim1 <- ifelse(corr_results$Pvalue_Dim1 < 0.001, "***",
                                ifelse(corr_results$Pvalue_Dim1 < 0.01, "**",
                                       ifelse(corr_results$Pvalue_Dim1 < 0.05, "*", "")))
corr_results$Sig_Dim2 <- ifelse(corr_results$Pvalue_Dim2 < 0.001, "***",
                                ifelse(corr_results$Pvalue_Dim2 < 0.01, "**",
                                       ifelse(corr_results$Pvalue_Dim2 < 0.05, "*", "")))

# 结果排序（按与Dim1的相关性）
corr_results <- corr_results[order(-abs(corr_results$Corr_Dim1)), ]

# 输出前20个最相关的密码子
head(corr_results, 20)

# 可视化相关系数矩阵
corr_matrix <- cor(species_rscu, coa_coords[, c("Dim1", "Dim2")], method = "pearson")
colnames(corr_matrix) <- c("Dim1", "Dim2")

# 创建相关系数热图
pdf(file.path(fig_out_dir, "ENC_pearson.pdf"), width = 4/2.54, height = 10/2.54)
p3 <-corrplot(
  corr_matrix, 
  method = "color", 
  tl.col = "black", 
  tl.cex = 0.7,
  cl.ratio = 0.2,
  cl.align.text = "l",
  mar = c(0, 0, 4, 0))
mtext(paste0("Dim1 (", dim1_var, "% variance)"), side = 3, line = 0.5, cex = 0.8, font = 2)
mtext(paste0("Dim2 (", dim2_var, "% variance)"), side = 3, line = 1.5, cex = 0.8, font = 2)
print(p3)
dev.off()
# 保存结果到CSV文件
write.csv(corr_results, "RSCU_COA_Correlations.csv", row.names = FALSE)

