setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2_add_Hou"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
# 加载必要的包
library(readxl)
library(dplyr)
library(tidyr)

# 读取Excel文件
df <- read_excel("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou/all菌（new）.xlsx", sheet = "alljun")

# 提取指定的ASV列，并按照指定顺序排列
selected_asvs <- c("ASV8", "ASV1831", "ASV54", "ASV30", "ASV73", "ASV9", "ASV17", 
                   "ASV4", "ASV28", "ASV3787", "ASV45", "ASV1834", "ASV1089", 
                   "ASV1527", "ASV38", "ASV48", "ASV2860", "ASV22", "ASV1515", 
                   "ASV2874", "ASV2", "ASV2839", "ASV3402", "ASV1540", "ASV167")

# 找到这些ASV在数据框中的列位置
asv_columns <- which(colnames(df) %in% selected_asvs)

# 提取allAA_AT3和allAA_GC3行
at3_row <- which(df$asv == "allAA_AT3")
gc3_row <- which(df$asv == "allAA_GC3")

# 创建转置的新表
new_table <- data.frame(
  row.names = c("allAA_AT3", "allAA_GC3")
)

# 按照指定顺序添加每个ASV列的数据
for(asv_name in selected_asvs) {
  if(asv_name %in% colnames(df)) {
    # 确保数据是数值类型
    at3_value <- as.numeric(df[at3_row, asv_name])
    gc3_value <- as.numeric(df[gc3_row, asv_name])
    
    new_table[[asv_name]] <- c(at3_value, gc3_value)
  }
}

# 查看结果
print("最终结果:")
print(new_table)

# 保存到新文件
write.csv(new_table, "Val_ASV_data.csv")


