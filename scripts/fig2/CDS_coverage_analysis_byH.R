# ============================================
# 📊 CDS Coverage Analysis Pipeline
# Author: ChatGPT
# Version: 1.0
# ============================================

# ---- 1. 自动安装并加载依赖 ----
required_packages <- c("GenomicAlignments", "GenomicRanges", "dplyr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("GenomicAlignments", "GenomicRanges")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# ---- 2. 参数设置 ----
gene_id <- "ENSSSCG00000007116"   # 🧬 目标基因 Ensembl ID

# 🧪 BAM 文件列表（可替换）
#bam_files <- list(
 # H_RNAseq   = "D:/DATA/PBMC/fig2_RNAseq_bam/sorted_bam/RNAseq_H7.sorted.bam",
  #W_RNAseq   = "D:/DATA/PBMC/fig2_RNAseq_bam/sorted_bam/RNAseq_w10.sorted.bam",
  #H_Riboseq  = "D:/DATA/PBMC/new_fig2_ribo_bam/original_bam/sorted_bam/ribo_H7_PBMC.sorted.bam",
  #W_Riboseq  = "D:/DATA/PBMC/new_fig2_ribo_bam/original_bam/sorted_bam/ribo_w10_PBMC.sorted.bam"
#)
bam_files <- list(
H_RNAseq   = "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_raw/fig2_rna_bam/RNAseq_H6.sorted.bam",
W_RNAseq   = "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_raw/fig2_rna_bam/RNAseq_W6.sorted.bam",
H_Riboseq  = "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_raw/fig2_ribo_bam_new/ribo_H6_PBMC.sorted.bam",
W_Riboseq  = "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_raw/fig2_ribo_bam_new/ribo_W6_PBMC.sorted.bam"
)
# ---- 3. 输出路径设置（固定到指定目录下） ----
base_output_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou"
output_dir <- file.path(
  base_output_dir,
  paste0("CDS_Coverage_", gene_id, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2_add_Hou"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
cat("📁 输出目录:", normalizePath(getwd(), winslash = "/"), "\n")
# ---- 4. 主分析函数 ----
analyze_cds_coverage <- function(bam_files, gene_id) {
  results <- list()
  for (sample_name in names(bam_files)) {
    bam_file <- bam_files[[sample_name]]
    cat("\n🔍 正在处理样本:", sample_name, "\n")
    
    # 获取 BAM 中的序列信息
    bam_seqs <- scanBamHeader(bam_file)[[1]]$targets
    target_seq <- grep(gene_id, names(bam_seqs), value = TRUE)
    
    if (length(target_seq) == 0) {
      warning("❌ 样本 ", sample_name, " 中未找到基因 ", gene_id)
      next
    }
    
    seq_name <- target_seq[1]
    seq_length <- bam_seqs[seq_name]
    cat("✅ 找到目标序列:", seq_name, " (长度:", seq_length, "bp)\n")
    
    # 提取该 CDS 的比对并计算覆盖度
    param <- ScanBamParam(
      what = c("pos", "qwidth", "strand"),
      which = GRanges(seq_name, IRanges(1, seq_length))
    )
    aln <- readGAlignments(bam_file, param = param)
    cov <- coverage(aln)[[seq_name]]
    
    results[[sample_name]] <- list(
      sequence_name = seq_name,
      length = seq_length,
      coverage = as.numeric(cov),
      total_reads = length(aln)
    )
  }
  return(results)
}

# ---- 5. 运行分析 ----
cds_results <- analyze_cds_coverage(bam_files, gene_id)

# ---- 6. 可视化覆盖度（折线图）----
plot_cds_coverage <- function(results) {
  par(mfrow = c(2, 2))
  for (sample in names(results)) {
    if (!is.null(results[[sample]])) {
      plot(
        results[[sample]]$coverage,
        type = "l",
        main = paste0(sample, " - ", results[[sample]]$sequence_name),
        xlab = "CDS位置 (bp)",
        ylab = "覆盖度",
        col = "blue",
        lwd = 1.5
      )
      grid()
    }
  }
}
cat("\n📈 绘制覆盖度曲线...\n")
plot_cds_coverage(cds_results)

# ---- 7. 创建覆盖度矩阵 CSV ----
create_coverage_dataframe <- function(results) {
  max_length <- max(sapply(results, function(x) length(x$coverage)))
  df <- data.frame(position = 1:max_length)
  for (sample_name in names(results)) {
    cov_vec <- results[[sample_name]]$coverage
    if (length(cov_vec) < max_length) {
      cov_vec <- c(cov_vec, rep(NA, max_length - length(cov_vec)))
    }
    df[[sample_name]] <- cov_vec
  }
  return(df)
}

coverage_df <- create_coverage_dataframe(cds_results)
write.csv(coverage_df, paste0("CDS_", gene_id, "_coverage.csv"), row.names = FALSE)

# ---- 8. 统计信息 ----
create_stats_dataframe <- function(results) {
  data.frame(
    Sample = names(results),
    Sequence_Name = sapply(results, `[[`, "sequence_name"),
    CDS_Length = sapply(results, `[[`, "length"),
    Total_Reads = sapply(results, `[[`, "total_reads"),
    Mean_Coverage = sapply(results, function(x) mean(x$coverage, na.rm = TRUE)),
    Median_Coverage = sapply(results, function(x) median(x$coverage, na.rm = TRUE)),
    Max_Coverage = sapply(results, function(x) max(x$coverage, na.rm = TRUE)),
    Coverage_Percentage_Above_1 = sapply(results, function(x) mean(x$coverage >= 1, na.rm = TRUE) * 100)
  )
}

stats_df <- create_stats_dataframe(cds_results)
write.csv(stats_df, paste0("CDS_", gene_id, "_stats.csv"), row.names = FALSE)

# ---- 9. 归一化覆盖度 ----
create_normalized_coverage <- function(results) {
  cov_df <- create_coverage_dataframe(results)
  norm_df <- data.frame(position = cov_df$position)
  for (sample_name in names(results)) {
    max_cov <- max(cov_df[[sample_name]], na.rm = TRUE)
    norm_df[[sample_name]] <- ifelse(max_cov > 0, cov_df[[sample_name]] / max_cov, cov_df[[sample_name]])
  }
  return(norm_df)
}

normalized_df <- create_normalized_coverage(cds_results)
write.csv(normalized_df, paste0("CDS_", gene_id, "_normalized_coverage.csv"), row.names = FALSE)

# ---- 10. 摘要报告 ----
summary_txt <- paste(
  "📊 CDS覆盖度分析摘要",
  paste("目标基因:", gene_id),
  paste("样本数:", length(cds_results)),
  paste("生成时间:", Sys.time()),
  "\n📈 统计信息表：",
  capture.output(print(stats_df)),
  sep = "\n"
)
writeLines(summary_txt, paste0("CDS_", gene_id, "_summary.txt"))

cat("\n✅ 全部分析完成！\n")
cat("📁 结果保存于:", normalizePath(output_dir), "\n")



