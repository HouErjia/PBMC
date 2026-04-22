rm(list=ls())
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig5")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig5"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
library(vegan)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(readxl)
library(compositions)
mytheme_noborder<-theme_bw()+theme(panel.border=element_blank(),panel.grid=element_blank(),axis.title=element_text(size=9,colour='black'),axis.text=element_text(size=9,colour='black'),axis.ticks=element_line(size=0.4,colour='black'),axis.ticks.length=unit(0.1,'cm'),axis.line.x=element_line(size=0.4,colour='black'),axis.line.y=element_line(size=0.4,colour='black'))
mytheme<-theme_bw()+theme(panel.border=element_rect(size=0.4,colour='black'),panel.grid=element_blank(),axis.title=element_text(size=9,colour='black'),axis.text=element_text(size=9,colour='black'),axis.ticks=element_line(size=0.4,colour='black'),axis.ticks.length=unit(0.1,'cm'))

jum <- read_excel("./asv_correlation_jum_codon.xlsx", sheet = 3, col_names = TRUE)
codon <- read_excel("./asv_correlation_jum_codon.xlsx", sheet = 5, col_names = TRUE)
jum=as.data.frame(jum)
codon=as.data.frame(codon)
#菌群丰度标准化：中心对数比变换 (CLR)
jum_clr <- clr(jum[, -1] + 1)   # 加1避免零值
rownames(jum_clr) <- jum$GCA
colnames(jum_clr) 

#密码子指标标准化：z-score标准化
codon_z <- scale(codon[, -1])
rownames(codon_z) <- codon$jum


# **1. 菌群丰度与样品的聚类分析**
# 计算样品间距离（基于菌群丰度）
distance_jum <- dist(t(jum_clr))
hc_jum <- hclust(distance_jum)

# 绘制聚类图并标记组别（前10腹泻，后10健康）
plot(hc_jum, main = "菌群丰度样品聚类（腹泻组 vs 健康组）", 
     sub = "红：腹泻 | 蓝：健康", cex = 0.8)
groups <- c(rep("#B53F48", 10), rep("#4B5DA0", 10))
order <- hc_jum$order
rect.hclust(hc_jum, k = 2, border = "black")
axis(1, at = 1:20, 
     labels = ifelse(1:20 %in% 1:10, "腹泻", "健康"), 
     col = groups[order], las = 2)


# **2. 密码子使用与菌种的聚类分析**
# 计算菌种间距离（基于密码子指标）
distance_codon <- dist(codon_z)
hc_codon <- hclust(distance_codon)

# 绘制菌种聚类图
plot(hc_codon, main = "菌种密码子使用模式聚类", 
     sub = "菌种间相似性", cex = 0.8)

# **3. 密码子指标特征聚类（密码子间相似性）**
distance_codon_features <- dist(t(codon_z))
hc_codon_features <- hclust(distance_codon_features)

plot(hc_codon_features, main = "密码子特征聚类（密码子间相似性）", 
     sub = "密码子指标的聚类", cex = 0.8)

# **4. 联合分析：菌群与密码子数据的PCA**
# 合并数据（菌种为行，菌群丰度+密码子指标为列）
merged_data <- cbind(jum_clr, codon_z)

# PCA分析
pca <- prcomp(merged_data, scale. = TRUE)
summary(pca)

# 绘制PCA图
plot(pca$x, 
     main = "PCA of 菌群丰度与密码子指标",
     pch = 19,
     col = c(rep("#B53F48", 10), rep("#4B5DA0", 10)),
     xlab = paste0("PC1 (", round(pca$sdev[1]^2/sum(pca$sdev^2)*100, 1), "%)"),
     ylab = paste0("PC2 (", round(pca$sdev[2]^2/sum(pca$sdev^2)*100, 1), "%)"))
legend("topright", 
       legend = c("腹泻组", "健康组"),
       col = c("#B53F48", "#4B5DA0"),
       pch = 19)


# **5. 热图可视化**
# 菌群丰度热图（样品×菌种）

jum_melt <- melt(jum_clr)
colnames(jum_melt)
jum_clr_matrix <- as.matrix(jum_clr)
rownames(jum_clr_matrix) <- rownames(jum_clr)
colnames(jum_clr_matrix) <- colnames(jum_clr)
custom_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

p <-pheatmap( jum_clr_matrix,color = custom_colors,scale = "row",
  cluster_rows = TRUE, cluster_cols = TRUE, main = "菌群丰度热图（样品×菌种）",
  show_rownames = TRUE, show_colnames = TRUE,fontsize = 8,legend = TRUE,legend_main = "标准化丰度")
ggsave(p,filename = file.path(fig_out_dir, "correlation_jum_fengdu.pdf"),units = "in",width=15/2.54,height=18/2.54)

# 密码子指标热图（菌种×密码子）
codon_melt <- melt(codon_z)
library(pheatmap)
p1 <-pheatmap(codon_z,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "密码子使用热图（聚类）",
         show_rownames = TRUE,
         show_colnames = TRUE)
ggsave(p1, filename = file.path(fig_out_dir, "correlation_asv_codon.pdf"), units = "in", width = 22/2.54, height = 18/2.54)

# 6. 关联性分析：菌群与密码子指标的相关性矩阵**
common_samples <- intersect(rownames(jum_clr), rownames(codon_z))
jum_clr_common <- jum_clr[common_samples, ]
codon_z_common <- codon_z[common_samples, ]
jum_mat <- as.matrix(jum_clr_common)
codon_mat <- as.matrix(codon_z_common)

block_size <- 50  # 减小分块大小以避免内存溢出
cor_matrix <- matrix(NA, 
                     nrow = ncol(jum_mat), 
                     ncol = ncol(codon_mat),
                     dimnames = list(colnames(jum_mat), colnames(codon_mat)))

for (i in seq(1, ncol(jum_mat), by = block_size)) {
  end_i <- min(i + block_size - 1, ncol(jum_mat))
  for (j in seq(1, ncol(codon_mat), by = block_size)) {
    end_j <- min(j + block_size - 1, ncol(codon_mat))
    sub_cor <- cor(
      jum_mat[, i:end_i], 
      codon_mat[, j:end_j], 
      method = "pearson",
      use = "pairwise.complete.obs"
    )
    cor_matrix[i:end_i, j:end_j] <- sub_cor
    gc() } }
# 绘制相关性热图
library("xlsx")
library(pheatmap)
p3 <- pheatmap(cor_matrix,
               color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
               main = "菌群丰度与密码子使用的相关性热图",
               show_rownames = TRUE,
               show_colnames = TRUE)
ggsave(p3, filename = file.path(fig_out_dir, "PBMC_correlation_jum_codon.pdf"), units = "in", width =25/2.54, height = 12/2.54)

saveRDS(cor_matrix, "./PBMC_correlation_jum_codon_matrix.rds")
write.csv(cor_matrix, "./PBMC_correlation_jum_codon_matrix.csv", row.names = TRUE)

# 绘制全密码子热图（直接保存为PDF）
#手动改表
cor_matrix <- read.csv("./PBMC_correlation_jum_codon_matrix.csv", row.names = 1)
pheatmap(
  cor_matrix_modified,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "菌群丰度与密码子使用的相关性热图",
  show_rownames = TRUE,
  show_colnames = TRUE,
  filename = file.path(fig_out_dir, "PBMC_correlation_jum_codon_full.pdf"),  # 直接保存为PDF
  width = 25/2.54, 
  height = 12/2.54
)

# 1. 挑选特定密码子
#手动改表
library("xlsx")
cor_matrix<-read.xlsx("./correlation_jum_codon_matrix.xls",1,header=T,row.names=1)
selected_codons <- c("CTT","CTC","CTA","CTG","TTA","TTG")
codon_selected_cols <- grep(paste(selected_codons, collapse = "|"), 
                            colnames(cor_matrix),
                            ignore.case = TRUE)
selected_cor_matrix <- cor_matrix[, codon_selected_cols]

# 绘制特定密码子热图
library(pheatmap)
p4 <- pheatmap(
  selected_cor_matrix,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  main = "Heat map of the correlation between Phe and flora abundance",
  show_rownames = TRUE,
  show_colnames = TRUE)
ggsave(p4, filename = file.path(fig_out_dir, "correlation_Leu_fengdu.pdf"), units = "in", width = 7/2.54, height = 11/2.54)

# 2. 密码子分类及计算平均值
codon_third_base <- substr(colnames(codon_mat), 3, 3)
group_AT <- codon_third_base %in% c("A", "T")
group_GC <- codon_third_base %in% c("G", "C")

cols_AT <- which(group_AT)
cols_GC <- which(group_GC)

mean_AT <- rowMeans(cor_matrix[, cols_AT, drop = FALSE], na.rm = TRUE)
mean_GC <- rowMeans(cor_matrix[, cols_GC, drop = FALSE], na.rm = TRUE)

codon_group_mean <- data.frame(
  Microbe = rownames(cor_matrix),
  `A/T Mean` = mean_AT,
  `G/C Mean` = mean_GC)

# 转换为长格式
codon_group_melt <- melt(codon_group_mean, id.vars = "Microbe", variable.name = "Group",value.name = "Correlation")

# 绘制箱线图
p5 <- ggplot(codon_group_melt, aes(x = Group, y = Correlation, fill = Group)) +
  geom_boxplot() +
  labs(title = "菌群丰度与两类密码子的平均相关性",
       x = "密码子第三位碱基类型",
       y = "平均相关系数") +
  mytheme_noborder
print(p5)
ggsave(p5, filename = file.path(fig_out_dir, "correlation_GC3mean_fengdu.pdf"), units = "in", width = 7.5/2.54, height = 6/2.54)
#密码子第三位碱基类型与平均相关系数之间相关性较弱，但G.C.Mean组可能略强于A.T.Mean组，需进一步统计验证

# 绘制分组热图
p6 <- pheatmap(
  cbind(mean_AT, mean_GC),
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_cols = FALSE)
ggsave(p6, filename = file.path(fig_out_dir, "correlation_GC3AT3_fengdu.pdf"), units = "in", width = 7/2.54, height = 12/2.54)

