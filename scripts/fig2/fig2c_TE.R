rm(list=ls())
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
rpkm_rna <- read.csv(file = "./fig2c_RNA_PBMC_rpkm.csv", row.names = 1, header = TRUE)
rpkm_ribo <- read.csv(file = "./fig2c_ribo_PBMC_rpkm.csv", row.names = 1, header = TRUE)
#先算平均
mean_rpkm_rna_H <- rowMeans(rpkm_rna[, 1:5])
mean_rpkm_rna_W <- rowMeans(rpkm_rna[, 6:10])
mean_rpkm_ribo_H <- rowMeans(rpkm_ribo[, 1:2])
mean_rpkm_ribo_W <- rowMeans(rpkm_ribo[, 3:4])
#再取交集，之后过滤
common_H <- intersect(names(mean_rpkm_rna_H), names(mean_rpkm_ribo_H))
common_W <- intersect(names(mean_rpkm_rna_W), names(mean_rpkm_ribo_W))
ficommon_H <- common_H[mean_rpkm_rna_H[common_H] >= 1 & mean_rpkm_ribo_H[common_H] >= 1]
ficommon_W <- common_W[mean_rpkm_rna_W[common_W] >= 1 & mean_rpkm_ribo_W[common_W] >= 1]

# 对过滤后的交集进行排序
sort_ficommon_H <- sort(mean_rpkm_rna_H[ficommon_H], decreasing = TRUE)
sort_ficommon_W <- sort(mean_rpkm_rna_W[ficommon_W], decreasing = TRUE)
fi_rna_H <- mean_rpkm_rna_H[ficommon_H]
fi_rna_W <- mean_rpkm_rna_W[ficommon_W]
fi_ribo_H <- mean_rpkm_ribo_H[ficommon_H]
fi_ribo_W <- mean_rpkm_ribo_W[ficommon_W]

TE_H <- fi_ribo_H / fi_rna_H 
TE_W <- fi_ribo_W / fi_rna_W 
TE_df_H <- data.frame(Transcript = names(sort_ficommon_H), TE_H = TE_H[names(sort_ficommon_H)])
TE_df_W <- data.frame(Transcript = names(sort_ficommon_W), TE_W = TE_W[names(sort_ficommon_W)])

# 合并两个数据框，找到共同的转录本
combined_TE_df <- merge(TE_df_H, TE_df_W, by = "Transcript", all = FALSE)

library(dplyr)
# 计算每个转录本在两组之间的log2 fold change和p值
test_results <- combined_TE_df %>%
  mutate(
    log2_fold_change = log2(TE_W / TE_H),
    regulation = case_when(
      log2_fold_change > 2 ~ "up",
      log2_fold_change < -1/2 ~ "down",
      TRUE ~ "unchanged"
    )
)

significant_transcripts <- test_results %>%
  filter(regulation != "unchanged")
# 输出结果到CSV文件
write.table(significant_transcripts, file = "./translation_efficiency_changes.csv", sep = ",", row.names = FALSE, quote = FALSE)


