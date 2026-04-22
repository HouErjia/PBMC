setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2_add_Hou"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
# ============================================================
# 翻译效率 (Translation Efficiency, TE) 自动计算脚本
# 自动检测样本列 (H*, W*)，计算 TE = Ribo_RPKM / RNA_RPKM
# ============================================================

# ---- 安装与加载包 ----
packages <- c("biomaRt", "dplyr", "readr", "tidyr", "stringr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
# ---- 确保加载正确的 dplyr ----
# 指定使用 dplyr::select 防止函数冲突
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
# ---- 建立 Ensembl 连接 ----
setup_ensembl_connection <- function() {
  tryCatch({
    useEnsembl(biomart = "genes", dataset = "sscrofa_gene_ensembl")
  }, error = function(e) {
    message("主镜像连接失败，尝试备用镜像...")
    tryCatch({
      useEnsembl(biomart = "genes", dataset = "sscrofa_gene_ensembl", mirror = "useast")
    }, error = function(e2) {
      stop("所有镜像连接失败，请检查网络连接。")
    })
  })
}

# ---- 批量转换转录本ID → 基因信息 ----
transcript_to_gene_batch <- function(transcript_ids, mart) {
  batch_size <- 500
  n_batches <- ceiling(length(transcript_ids) / batch_size)
  all_results <- data.frame()
  
  for (i in seq_len(n_batches)) {
    start <- (i - 1) * batch_size + 1
    end <- min(i * batch_size, length(transcript_ids))
    ids <- transcript_ids[start:end]
    message(sprintf("正在处理第 %d/%d 批 (%d 个转录本)...", i, n_batches, length(ids)))
    
    tryCatch({
      res <- getBM(
        attributes = c(
          "ensembl_transcript_id",
          "ensembl_gene_id",
          "external_gene_name",
          "gene_biotype",
          "description"
        ),
        filters = "ensembl_transcript_id",
        values = ids,
        mart = mart
      )
      all_results <- bind_rows(all_results, res)
      if (i < n_batches) Sys.sleep(1)
    }, error = function(e) {
      message(sprintf("第 %d 批失败: %s", i, e$message))
    })
  }
  return(all_results)
}

# ---- 主函数 ----
calculate_TE_auto <- function(rna_file = "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2/fig2c_RNA_PBMC_rpkm.csv",
                              ribo_file = "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2/fig2c_ribo_PBMC_rpkm.csv",
                              output_prefix = "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou/translation_efficiency_auto") {
  
  message("=== 读取数据 ===")
  rna <- read_csv(rna_file)
  ribo <- read_csv(ribo_file)
  
  # 确保第一列是 Transcript ID
  colnames(rna)[1] <- "Transcript"
  colnames(ribo)[1] <- "Transcript"
  
  # 自动检测样本列（以 H 或 W 开头）
  rna_samples <- colnames(rna)[str_detect(colnames(rna), "^(H|W)")]
  ribo_samples <- colnames(ribo)[str_detect(colnames(ribo), "^(H|W)")]
  
  shared_samples <- intersect(rna_samples, ribo_samples)
  message(sprintf("检测到共有样本列: %s", paste(shared_samples, collapse = ", ")))
  
  # 提取样本列 + 转录本列
  rna_sel <- rna %>% select(any_of(c("Transcript", shared_samples)))
  ribo_sel <- ribo %>% select(any_of(c("Transcript", shared_samples)))
  
  # 取共有转录本
  shared_ids <- intersect(rna_sel$Transcript, ribo_sel$Transcript)
  message(sprintf("共有转录本数量: %d", length(shared_ids)))
  
  rna_sel <- rna_sel %>% filter(Transcript %in% shared_ids)
  ribo_sel <- ribo_sel %>% filter(Transcript %in% shared_ids)
  
  # 转为长格式
  rna_long <- rna_sel %>%
    pivot_longer(-Transcript, names_to = "Sample", values_to = "RNA_RPKM")
  ribo_long <- ribo_sel %>%
    pivot_longer(-Transcript, names_to = "Sample", values_to = "Ribo_RPKM")
  
  # 合并并计算 TE
  combined <- inner_join(rna_long, ribo_long, by = c("Transcript", "Sample"))
  combined <- combined %>%
    mutate(
      TE = ifelse(RNA_RPKM > 0, Ribo_RPKM / RNA_RPKM, NA),
      Log2_TE = ifelse(!is.na(TE) & TE > 0, log2(TE), NA)
    )
  
  # ---- 转换转录本ID到基因 ----
  message("=== Ensembl ID 转换 ===")
  mart <- setup_ensembl_connection()
  gene_info <- transcript_to_gene_batch(unique(combined$Transcript), mart)
  
  combined <- left_join(combined,
                        gene_info,
                        by = c("Transcript" = "ensembl_transcript_id"))
  
  # ---- 宽格式输出 ----
  te_wide <- combined %>%
    select(Transcript, external_gene_name, Sample, TE) %>%
    pivot_wider(names_from = Sample, values_from = TE, names_prefix = "TE_")
  
  # ---- 保存结果 ----
  long_file <- paste0(output_prefix, "_long.csv")
  wide_file <- paste0(output_prefix, "_wide.csv")
  write_csv(combined, long_file)
  write_csv(te_wide, wide_file)
  
  message("\n=== 翻译效率计算完成 ===")
  message(sprintf("保存文件:\n  %s\n  %s", long_file, wide_file))
  
  print(head(te_wide))
  return(list(long = combined, wide = te_wide))
}

# ---- 执行 ----
result <- calculate_TE_auto()

