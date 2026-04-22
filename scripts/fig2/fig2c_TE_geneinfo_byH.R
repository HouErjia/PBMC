setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
# 转录本ID批量转换为基因信息的R脚本
# 需要先安装必要的包

# 安装必要的包（如果尚未安装）
if (!require("biomaRt")) {
  install.packages("biomaRt")
}
if (!require("dplyr")) {
  install.packages("dplyr")
}
if (!require("readr")) {
  install.packages("readr")
}

# 加载包
library(biomaRt)
library(dplyr)
library(readr)

# 设置Ensembl连接
setup_ensembl_connection <- function() {
  tryCatch({
    # 尝试连接猪的Ensembl数据库
    ensembl <- useEnsembl(biomart = "genes", dataset = "sscrofa_gene_ensembl")
    return(ensembl)
  }, error = function(e) {
    # 如果连接失败，尝试其他镜像
    message("主镜像连接失败，尝试使用备用镜像...")
    tryCatch({
      ensembl <- useEnsembl(biomart = "genes", 
                           dataset = "sscrofa_gene_ensembl",
                           mirror = "useast")
      return(ensembl)
    }, error = function(e2) {
      message("所有镜像连接失败，请检查网络连接")
      return(NULL)
    })
  })
}

# 批量转换转录本ID到基因信息
transcript_to_gene_batch <- function(transcript_ids, mart) {
  if (is.null(mart)) {
    stop("Ensembl连接失败，无法进行转换")
  }
  
  # 分批处理，避免请求过大
  batch_size <- 500
  n_batches <- ceiling(length(transcript_ids) / batch_size)
  all_results <- data.frame()
  
  for (i in 1:n_batches) {
    start_idx <- (i-1) * batch_size + 1
    end_idx <- min(i * batch_size, length(transcript_ids))
    batch_ids <- transcript_ids[start_idx:end_idx]
    
    message(sprintf("正在处理第 %d/%d 批 (%d 个转录本)...", i, n_batches, length(batch_ids)))
    
    tryCatch({
      # 获取基因信息
      batch_results <- getBM(
        attributes = c(
          "ensembl_transcript_id", 
          "ensembl_gene_id", 
          "external_gene_name",
          "gene_biotype",
          "description"
        ),
        filters = "ensembl_transcript_id",
        values = batch_ids,
        mart = mart
      )
      
      all_results <- bind_rows(all_results, batch_results)
      
      # 添加延迟，避免请求过于频繁
      if (i < n_batches) {
        Sys.sleep(1)
      }
      
    }, error = function(e) {
      message(sprintf("第 %d 批处理失败: %s", i, e$message))
    })
  }
  
  return(all_results)
}

# 主函数
main <- function() {
  # 读取CSV文件
  message("正在读取CSV文件...")
  data <- read_csv("./fig2_translation_efficiency_changes.csv")
  
  # 提取唯一的转录本ID
  transcript_ids <- unique(data$Transcript)
  message(sprintf("找到 %d 个唯一的转录本ID", length(transcript_ids)))
  
  # 建立Ensembl连接
  message("正在连接Ensembl数据库...")
  ensembl_mart <- setup_ensembl_connection()
  
  if (is.null(ensembl_mart)) {
    stop("无法建立Ensembl连接，请检查网络设置")
  }
  
  # 批量转换
  message("开始批量转换转录本ID到基因信息...")
  gene_info <- transcript_to_gene_batch(transcript_ids, ensembl_mart)
  
  # 检查转换结果
  if (nrow(gene_info) == 0) {
    stop("没有获取到任何基因信息")
  }
  
  message(sprintf("成功转换了 %d/%d 个转录本", 
                 nrow(gene_info), length(transcript_ids)))
  
  # 重命名列以便合并
  colnames(gene_info)[1] <- "Transcript"
  
  # 合并回原数据
  message("正在合并结果到原数据...")
  final_data <- data %>%
    left_join(gene_info, by = "Transcript")
  
  # 保存结果
  output_file <- "./fig2c_translation_efficiency_changes_with_gene_info.csv"
  write_csv(final_data, output_file)
  
  # 生成转换统计报告
  converted_count <- sum(!is.na(final_data$ensembl_gene_id))
  conversion_rate <- converted_count / nrow(final_data) * 100
  
  message("\n=== 转换完成 ===")
  message(sprintf("输入转录本数量: %d", length(transcript_ids)))
  message(sprintf("成功转换数量: %d", nrow(gene_info)))
  message(sprintf("最终数据转换率: %.2f%%", conversion_rate))
  message(sprintf("结果已保存到: %s", output_file))
  
  # 显示前几行结果预览
  message("\n结果预览:")
  print(head(final_data[, c("Transcript", "ensembl_gene_id", "external_gene_name", "regulation")]))
  
  return(final_data)
}

# 执行主函数
result <- main()

