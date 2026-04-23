# PBMC Paper Package

该目录按论文复现优先的方式整理，本版结构统计基于 D:\PBMC交接\PBMC_sqluo\PBMC_paper_package\PBMC 当前文件树自动生成。

## 目录角色
- data_raw/: 原始输入与参考文件（BAM/BAM 索引、原始统计输出、原始校验文件）
- scripts/: 分析脚本（.R、.pl、.ipynb 等）
- data_processed/: 中间结果与衍生表格
- metadata/: 样本与分组元数据
- env/: 环境信息（sessionInfo.txt）

## 文件统计（MANIFEST.tsv，更新于 2026-04-23）
- tracked files total: 139
- data_processed: 97
- data_raw: 5
- scripts: 34
- metadata: 2
- env: 1

## 追溯说明
- MANIFEST.tsv 记录逐文件的 category/source/destination/size_bytes。
- 本版 MANIFEST 基于 PBMC 当前文件结构生成，source 与 destination 均为仓库内实际绝对路径。
- README.md 与 MANIFEST.tsv 本身不纳入 MANIFEST.tsv 统计。
