# PBMC Paper Package

该目录按论文复现优先的方式整理，当前结构基于 D:\PBMC交接\PBMC_sqluo 主目录，并合并了 D:\PBMC交接 下的补充资料。

## 目录角色
- data_raw/: 原始输入与参考文件（BAM/BAM 索引、原始统计输出、原始校验文件）
- scripts/: 分析脚本（.R、.pl、.ipynb 等）
- data_processed/: 中间结果与衍生表格
- metadata/: 样本与分组元数据
- env/: 环境信息（sessionInfo.txt）

## 外部补充来源
- D:\PBMC交接\fig2_ribo_bam_new -> data_raw/fig2_ribo_bam_new/
- D:\PBMC交接\fig2_add_Hou -> scripts/fig2_add_Hou/、data_processed/fig2_add_Hou/、results/figures/fig2_add_Hou/

## 文件统计（MANIFEST.tsv，更新于 2026-04-22）
- data_processed: 168
- data_raw: 67
- scripts: 52
- metadata: 2
- env: 1

## 追溯说明
- MANIFEST.tsv 记录逐文件的 category/source/destination/size_bytes。
- data_raw/checksums/md5_raw.txt 保存原始数据校验信息。
- README.md 与 MANIFEST.tsv 本身不纳入 MANIFEST.tsv 统计。
