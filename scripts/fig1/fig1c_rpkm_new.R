rm(list=ls())
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig1")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig1"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
library(limma)
library(edgeR)
library(rtracklayer)
library(parallel)

# 导入GTF文件（请将注释文件放入 data_raw/reference/）
gtf_path <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_raw/reference/Sus_scrofa.Sscrofa11.1.113.gtf.gz"
gtf <-as.data.frame(x = rtracklayer::import(con = gtf_path))
table(gtf$type)

# 提取所有外显子并按基因名称分组
exon <- gtf[gtf$type == "exon", c("start", "end", "gene_name","transcript_id")]
exon_by_transcript_id <- split(exon, exon$transcript_id)

# 创建集群进行并行处理
num_cores <- floor(0.5 * detectCores())
cl <- makeCluster(num_cores)

# 并行计算每个基因的总长度
gene_length <- parLapply(cl = cl, X = exon_by_transcript_id, fun = function(x) {
  sum(x$end - x$start + 1) # 计算每个外显子长度并求和
})

stopCluster(cl) # 停止集群

gene_length <- data.frame(transcript_id = names(gene_length), length = as.numeric(gene_length))
write.table(x = gene_length, file = "./Sus_scrofa.Sscrofa11.1.11_length.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# 读取基因长度和读取计数
gene_length <- read.table(file = "./Sus_scrofa.Sscrofa11.1.11_length.tsv", header = TRUE, row.names = 1)
counts <- read.csv(file = "./newdata_transcripts.csv", row.names = 1, header = TRUE, check.names = FALSE)

# 筛选掉所有零表达的基因
counts <- counts[rowSums(counts) > 0, ]

# 取交集特征
intersect_features <- intersect(rownames(counts), rownames(gene_length))

# 筛选匹配的基因
counts <- counts[intersect_features, ]
gene_length <- gene_length[intersect_features, , drop = FALSE]

# 计算RPKM
total_counts <- colSums(counts)
gene_length_kb <- gene_length[, "length"] / 1000
rpkm <- t(t(counts) / gene_length_kb) / (total_counts / 1e6)

# 输出结果
write.table(rpkm, file = "./rpkm_values_new.csv", sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE)



