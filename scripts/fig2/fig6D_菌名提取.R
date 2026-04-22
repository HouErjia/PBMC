setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2_add_Hou")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2_add_Hou"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
# =========================
# 修正菌名映射提取问题
# =========================

# ---- 重新提取并处理数据 ----
# 使用你已有的bacteria_mapping数据

# 首先移除重复项并确保所有42个ASV都被包含
unique_mapping <- bacteria_mapping %>%
  distinct(asv_id, .keep_all = TRUE)

# 检查是否所有42个ASV都在
cat("原始ASV数量:", nrow(bacteria_mapping), "\n")
cat("去重后ASV数量:", nrow(unique_mapping), "\n")
cat("缺失的ASV:", setdiff(bacteria_mapping$asv_id, unique_mapping$asv_id), "\n")

# ---- 改进的菌名解析逻辑 ----
formatted_mapping <- unique_mapping %>%
  mutate(
    # 提取属名
    genus = str_match(full_name, "g__([^:_]+)")[,2],
    
    # 提取种名
    species = str_match(full_name, "s__([^:_]+)")[,2],
    
    # 提取科名（对于ASV1834这种情况）
    family = str_match(full_name, "f__([^:_]+)")[,2],
    
    # 改进的显示名称生成逻辑
    display_name = case_when(
      # 情况1: 有具体种名（如Pasteurella_aerogenes）
      !is.na(species) & !str_detect(species, "unclassified|uncultured") ~ 
        paste0(asv_id, " (", str_replace(species, "_", " "), ")"),
      
      # 情况2: 只有属名，使用spp.表示复数
      !is.na(genus) & is.na(species) ~ 
        paste0(asv_id, " (", genus, " spp.)"),
      
      # 情况3: 未分类的属，使用sp.表示单数
      !is.na(genus) & str_detect(species, "unclassified") ~ 
        paste0(asv_id, " (", genus, " sp.)"),
      
      # 情况4: 未培养的属，使用sp.表示单数
      !is.na(genus) & str_detect(species, "uncultured") ~ 
        paste0(asv_id, " (", genus, " sp.)"),
      
      # 情况5: 只有科名，没有属名（如ASV1834）
      !is.na(family) & is.na(genus) ~ 
        paste0(asv_id, " (", family, " sp.)"),
      
      # 默认情况：显示ASV编号
      TRUE ~ asv_id
    )
  )

# ---- 特殊处理一些已知的菌名 ----
# 手动修正一些特殊的菌名格式
special_cases <- c(
  "ASV55" = "ASV55 (Anaerotignum lactatifermentans)",
  "ASV31" = "ASV31 (Lactobacillus vaginalis Bacillus)", 
  "ASV3787" = "ASV3787 (UCG.002)",
  "ASV45" = "ASV45 (Escherichia Shigella)",
  "ASV1834" = "ASV1834 (Lachnospiraceae sp.)"  # 添加ASV1834的特殊处理
)

# 应用特殊处理
for(asv_id in names(special_cases)) {
  if(asv_id %in% formatted_mapping$asv_id) {
    formatted_mapping$display_name[formatted_mapping$asv_id == asv_id] <- special_cases[asv_id]
  }
}

# ---- 创建最终的菌名映射向量 ----
final_bacteria_names <- setNames(formatted_mapping$display_name, formatted_mapping$asv_id)

# 添加PBMC样本
final_bacteria_names <- c(
  final_bacteria_names,
  "W6_PBMC" = "W6_PBMC",
  "W10_PBMC" = "W10_PBMC", 
  "H7_PBMC" = "H7_PBMC",
  "H6_PBMC" = "H6_PBMC"
)

# ---- 输出结果 ----
cat("=== 修正后的菌名映射 ===\n")
cat("总样本数量:", length(final_bacteria_names), "\n\n")

for(i in 1:length(final_bacteria_names)) {
  cat(names(final_bacteria_names)[i], " → ", final_bacteria_names[i], "\n")
}

# ---- 检查是否所有原始ASV都在最终结果中 ----
missing_asvs <- setdiff(bacteria_mapping$asv_id, names(final_bacteria_names))
if(length(missing_asvs) > 0) {
  cat("\n警告: 以下ASV未包含在最终结果中:\n")
  cat(missing_asvs, "\n")
}

#--- 对于菌名查漏补缺---
# =========================
# 找出图2中有而图1没有的菌名并创建完整表格
# =========================

# ---- 图1：已生成的菌名映射 ----
figure1_names <- c(
  "ASV8" = "ASV8 (Coprococcus sp.)",
  "ASV1555" = "ASV1555 (Enterococcus sp.)",
  "ASV1831" = "ASV1831 (Eubacterium)",
  "ASV54" = "ASV54 (Clostridium sp.)",
  "ASV21" = "ASV21 (Clostridium)",
  "ASV30" = "ASV30 (Clostridium sp.)",
  "ASV73" = "ASV73 (Turicibacter sp.)",
  "ASV9" = "ASV9 (Romboutsia sp.)",
  "ASV17" = "ASV17 (Sharpea sp.)",
  "ASV5" = "ASV5 (UCG-002 sp.)",
  "ASV4" = "ASV4 (Christensenellaceae sp.)",
  "ASV28" = "ASV28 (Christensenellaceae sp.)",
  "ASV3787" = "ASV3787 (UCG.002)",
  "ASV1834" = "ASV1834 (Lachnospiraceae sp.)",
  "ASV22" = "ASV22 (Lactobacillus sp.)",
  "ASV3402" = "ASV3402 (Megasphaera sp.)",
  "ASV1540" = "ASV1540 (Megasphaera sp.)",
  "ASV77" = "ASV77 (Lactobacillus sp.)",
  "ASV1853" = "ASV1853 (Lactobacillus sp.)",
  "ASV1089" = "ASV1089 (Fusobacterium)",
  "ASV38" = "ASV38 (Blautia sp.)",
  "ASV48" = "ASV48 (Blautia sp.)",
  "ASV1515" = "ASV1515 (Streptococcus sp.)",
  "ASV2996" = "ASV2996 (Fusobacterium sp.)",
  "ASV2874" = "ASV2874 (Fusobacterium sp.)",
  "ASV167" = "ASV167 (Collinsella sp.)",
  "ASV2" = "ASV2 (Lactobacillus sp.)",
  "ASV318" = "ASV318 (Peptostreptococcus sp.)",
  "ASV2586" = "ASV2586 (Campylobacter sp.)",
  "ASV2860" = "ASV2860 (Phascolarctobacterium sp.)",
  "ASV2839" = "ASV2839 (Pasteurella)",
  "ASV1527" = "ASV1527 (Lactobacillus)",
  "ASV172" = "ASV172 (Mogibacterium sp.)",
  "ASV3058" = "ASV3058 (Allisonella sp.)",
  "W6_PBMC" = "W6_PBMC",
  "W10_PBMC" = "W10_PBMC", 
  "H7_PBMC" = "H7_PBMC",
  "H6_PBMC" = "H6_PBMC"
)

# ---- 图2：原始菌名列表 ----
figure2_raw <- c(
  "ASV8:s__unclassified_g__Coprococcus",
  "ASV1555:s__unclassified_g__Enterococcus",
  "ASV1831:s__Eubacterium_coprostanoligenes",
  "ASV54:s__unclassified_g__Clostridium_sensu_stricto_1",
  "ASV21:s__Clostridium_scindens",
  "ASV30:s__unclassified_g__Clostridium_sensu_stricto_1",
  "ASV73:s__unclassified_g__Turicibacter",
  "ASV9:s__unclassified_g__Romboutsia",
  "ASV17:s__uncultured_bacterium_g__Sharpea",
  "ASV5:s__unclassified_g__UCG-002",
  "ASV4:s__uncultured_spirochete_g__Christensenellaceae_R-7_group",
  "ASV28:s__uncultured_prokaryote_g__Christensenellaceae_R-7_group",
  "ASV3787:s__uncultured_organism_g__UCG-002",
  "ASV1834:s__unclassified_f__Lachnospiraceae",
  "ASV22:s__unclassified_g__Lactobacillus",
  "ASV3402:s__unclassified_g__Megasphaera",
  "ASV1540:s__unclassified_g__Megasphaera",
  "ASV77:s__unclassified_g__Lactobacillus",
  "ASV1853:s__unclassified_g__Lactobacillus",
  "ASV1089:s__Fusobacterium_gastrosuis",
  "ASV38:s__unclassified_g__Blautia",
  "ASV48:s__unclassified_g__Blautia",
  "ASV1515:s__unclassified_g__Streptococcus",
  "ASV2996:s__unclassified_g__Fusobacterium",
  "ASV2874:s__unclassified_g__Fusobacterium",
  "ASV167:s__unclassified_g__Collinsella",
  "ASV2:s__unclassified_g__Lactobacillus",
  "ASV3402:s__unclassified_g__Megasphaera",
  "ASV2874:s__unclassified_g__Fusobacterium",
  "ASV2996:s__unclassified_g__Fusobacterium",
  "ASV318:s__unclassified_g__Peptostreptococcus",
  "ASV2586:s__unclassified_g__Campylobacter",
  "ASV1515:s__unclassified_g__Streptococcus",
  "ASV5:s__unclassified_g__UCG-002",
  "ASV2860:s__unclassified_g__Phaseolaretobacterium",
  "ASV48:s__unclassified_g__Blautia",
  "ASV38:s__unclassified_g__Blautia",
  "ASV2839:s__Pasteurella_aerogenes",
  "ASV1527:s__Lactobacillus_vaginalis_g__Bacillus",
  "ASV172:s__uncultured_bacterium_g__Mogibacterium",
  "ASV1089:s__Fusobacterium_gastrosuis",
  "ASV3058:s__uncultured_bacterium_g__Allisonella"
)

# ---- 处理图2数据，提取ASV编号 ----
figure2_asvs <- unique(str_extract(figure2_raw, "^[^:]+"))

# ---- 找出图2中有而图1没有的ASV ----
missing_in_figure1 <- setdiff(figure2_asvs, names(figure1_names))

cat("图2中有而图1没有的ASV编号:\n")
print(missing_in_figure1)

# ---- 为缺失的ASV创建格式化名称 ----
# 首先创建图2的完整映射
figure2_mapping <- data.frame(
  asv_id = str_extract(figure2_raw, "^[^:]+"),
  full_name = figure2_raw,
  stringsAsFactors = FALSE
) %>%
  distinct(asv_id, .keep_all = TRUE)  # 去重，每个ASV只保留一次

# 为缺失的ASV生成格式化名称
missing_formatted <- figure2_mapping %>%
  filter(asv_id %in% missing_in_figure1) %>%
  mutate(
    # 提取属名
    genus = str_match(full_name, "g__([^:_]+)")[,2],
    
    # 提取种名
    species = str_match(full_name, "s__([^:_]+)")[,2],
    
    # 提取科名
    family = str_match(full_name, "f__([^:_]+)")[,2],
    
    # 生成显示名称
    display_name = case_when(
      # 情况1: 有具体种名
      !is.na(species) & !str_detect(species, "unclassified|uncultured") ~ 
        paste0(asv_id, " (", str_replace(species, "_", " "), ")"),
      
      # 情况2: 只有属名，使用spp.表示复数
      !is.na(genus) & is.na(species) ~ 
        paste0(asv_id, " (", genus, " spp.)"),
      
      # 情况3: 未分类的属，使用sp.表示单数
      !is.na(genus) & str_detect(species, "unclassified") ~ 
        paste0(asv_id, " (", genus, " sp.)"),
      
      # 情况4: 未培养的属，使用sp.表示单数
      !is.na(genus) & str_detect(species, "uncultured") ~ 
        paste0(asv_id, " (", genus, " sp.)"),
      
      # 情况5: 只有科名，没有属名
      !is.na(family) & is.na(genus) ~ 
        paste0(asv_id, " (", family, " sp.)"),
      
      # 默认情况：显示ASV编号
      TRUE ~ asv_id
    )
  )

# ---- 创建完整的新表格 ----
# 合并图1和缺失的ASV
complete_bacteria_names <- c(
  figure1_names,
  setNames(missing_formatted$display_name, missing_formatted$asv_id)
)

# ---- 输出结果 ----
cat("\n=== 完整的菌名映射表 ===\n")
cat("总样本数量:", length(complete_bacteria_names), "\n\n")

# 按ASV编号排序输出
sorted_names <- complete_bacteria_names[order(names(complete_bacteria_names))]
for(i in 1:length(sorted_names)) {
  cat(names(sorted_names)[i], " → ", sorted_names[i], "\n")
}

# ---- 保存完整表格 ----
complete_table <- data.frame(
  sample_id = names(complete_bacteria_names),
  display_name = complete_bacteria_names,
  stringsAsFactors = FALSE
) %>%
  arrange(sample_id)  # 按样本ID排序

write.csv(complete_table, "fig6D_fig6D_updated_complete_bacteria_mapping.csv", row.names = FALSE)

cat("\n完整菌名映射表已保存为: fig6D_fig6D_updated_complete_bacteria_mapping.csv\n")


