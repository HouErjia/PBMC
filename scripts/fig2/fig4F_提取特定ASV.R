setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2_add_Hou"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
#挑选需要列
# 加载必要的包
library(readxl)
library(dplyr)
if (!require(openxlsx)) {
  install.packages("openxlsx")
  library(openxlsx)
}

# 用户提供的ASV列表（按顺序）
asv_list <- c("ASV9", "ASV13", "ASV4", "ASV44", "ASV3", "ASV40", "ASV30", "ASV24", 
              "ASV1555", "ASV51", "ASV5", "ASV1902", "ASV119", "ASV8", "ASV73", 
              "ASV1844", "ASV10377", "ASV54", "ASV3787", "ASV67", "ASV28", "ASV15", 
              "ASV1861", "ASV1834", "ASV354", "ASV1832", "ASV575", "ASV180", "ASV34", 
              "ASV1876", "ASV1527", "ASV2", "ASV2874", "ASV77", "ASV2996", "ASV1515", 
              "ASV36", "ASV1997", "ASV2894", "ASV3058", "ASV19", "ASV3064", "ASV2909", 
              "ASV1853", "ASV3027", "ASV2839", "ASV2860", "ASV1089", "ASV318", 
              "ASV1514", "ASV2586", "ASV422", "ASV3094", "ASV211", "ASV172", "ASV3073", 
              "ASV2699", "ASV3063", "ASV65", "ASV3059")

# 读取Excel文件数据
data <- read_excel("副本各肠段ASV-all.xlsx", sheet = "Sheet1")

# 提取指定ASV并按顺序排列
result <- data %>% 
  filter(asv %in% asv_list) %>% 
  mutate(asv = factor(asv, levels = asv_list)) %>% 
  arrange(asv) %>% 
  mutate(asv = as.character(asv))  # 将因子转回字符（可选）

# 查看结果
print(result)
 
#保存结果
write.csv(result, "fig4F_ASV.csv", row.names = FALSE)


# 读取数据
data <- read.csv("fig4F_ASV.csv")

# 直接筛选并重新排列列
new_data <- data %>%
  select(
    asv, 
    # H开头的CaeC列
    matches("^H.*CaeC"),
    # W开头的CaeC列  
    matches("^W.*CaeC")
  )

# 查看结果
head(new_data)

# 保存结果
write.csv(new_data, "CaeC筛选数据.csv", row.names = FALSE)


# 更详细的版本（包含百分比格式）
library(dplyr)
library(openxlsx)


#整理顺序及求相对丰度
# 读取数据
data <- read.csv("CaeC筛选数据.csv")

# 提取和排序列名
H_columns <- grep("^H_", names(data), value = TRUE)
W_columns <- grep("^W_", names(data), value = TRUE)

sort_columns_by_number <- function(col_names) {
  numbers <- as.numeric(gsub("^[HW]_(\\d+).*", "\\1", col_names))
  col_names[order(numbers)]
}

H_columns_sorted <- sort_columns_by_number(H_columns)
W_columns_sorted <- sort_columns_by_number(W_columns)

# 重新排列数据
data_sorted <- data %>%
  select(asv, all_of(H_columns_sorted), all_of(W_columns_sorted))

# 计算相对丰度（百分比格式）
numeric_columns <- setdiff(names(data_sorted), "asv")
row_sums <- rowSums(data_sorted[, numeric_columns], na.rm = TRUE)

data_relative <- data_sorted
data_relative[, numeric_columns] <- sapply(data_sorted[, numeric_columns], function(x) {
  round(ifelse(row_sums > 0, (x / row_sums) * 100, 0), 4)  # 转换为百分比，保留4位小数
})

# 创建Excel文件
wb <- createWorkbook()

# 添加三个sheet
addWorksheet(wb, "排序后原始数据")
addWorksheet(wb, "相对丰度_百分比")
addWorksheet(wb, "汇总信息")

writeData(wb, "排序后原始数据", data_sorted)
writeData(wb, "相对丰度_百分比", data_relative)

# 添加汇总信息
summary_info <- data.frame(
  项目 = c("总ASV数量", "H组样本数", "W组样本数", "处理时间"),
  数值 = c(nrow(data_sorted), length(H_columns_sorted), length(W_columns_sorted), format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeData(wb, "汇总信息", summary_info)

# 保存文件
saveWorkbook(wb, "CaeC数据_完整分析.xlsx", overwrite = TRUE)

cat("分析完成！已生成包含3个sheet的Excel文件。\n")


#——————数据转置——————
# 安装必要的包（如果尚未安装）
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

# 加载包
library(readxl)
library(openxlsx)  # 用于写入Excel文件

# 读取原始Excel文件中的两个sheet
# 假设文件路径为：
file_path <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou/fig4F_CaeC数据_完整分析.xlsx"

# 读取绝对丰度sheet（根据您的数据，可能是"排序后原始数据"）
abs_abundance <- read_excel(file_path, sheet = "排序后原始数据")
# 读取相对丰度sheet（根据您的数据，是"相对丰度_百分比"）
rel_abundance <- read_excel(file_path, sheet = "相对丰度_百分比")

# 转置函数
transpose_df <- function(df) {
  # 第一列作为行名
  rownames <- df[[1]]
  # 转置数据部分
  t_df <- as.data.frame(t(df[-1]))
  # 设置列名
  colnames(t_df) <- rownames
  # 添加样本名列
  t_df <- cbind(Sample = rownames(t_df), t_df)
  rownames(t_df) <- NULL
  return(t_df)
}

# 转置两个表
abs_abundance_t <- transpose_df(abs_abundance)
rel_abundance_t <- transpose_df(rel_abundance)

# 查看转置后的数据结构
head(abs_abundance_t[, 1:5])  # 显示前5列
head(rel_abundance_t[, 1:5])

# 可选：将转置后的数据保存到新的Excel文件
output_file <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou/fig4F_CaeC_transposed.xlsx"
wb <- createWorkbook()

addWorksheet(wb, "绝对丰度_转置")
writeData(wb, "绝对丰度_转置", abs_abundance_t)

addWorksheet(wb, "相对丰度_转置")
writeData(wb, "相对丰度_转置", rel_abundance_t)

saveWorkbook(wb, output_file, overwrite = TRUE)

