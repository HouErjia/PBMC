setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2_add_Hou"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
# 读取刚才生成的宽格式数据
te_data <- read_csv("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou/translation_efficiency_auto_wide.csv")

# 定义要提取的基因列表
target_genes <- c("CD82", "CD93", "TNFSF13", "RELT", "CCL25", "CD24")

# 提取这些基因的数据
selected_genes <- te_data %>% 
  filter(external_gene_name %in% target_genes)

# 查看提取结果
print(selected_genes)

# 保存到新文件
write_csv(selected_genes, "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou/fig4F_selected_genes_TE.csv")

message("特定基因数据已保存到: D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou/fig4F_selected_genes_TE.csv")
#---- 读取数据-----
library(tibble)
data <- read.csv("fig4F_selected_genes_TE.csv")

# 转置数据
# 首先提取基因名作为列名，TE样本作为行名
transposed_data <- data %>%
  select(-Transcript) %>%  # 移除转录本ID列
  column_to_rownames("external_gene_name") %>%  # 将基因名设为行名
  t() %>%  # 转置矩阵
  as.data.frame()  # 转换回数据框

# 查看转置后的数据
print(transposed_data)

# 如果需要保存转置后的数据
write.csv(transposed_data, "fig4F_selected_genes_TE_transposed.csv", row.names = TRUE)

