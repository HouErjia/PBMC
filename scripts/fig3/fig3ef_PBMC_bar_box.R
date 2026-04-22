rm(list=ls())
library(dplyr)
library(ggplot2)
library(readxl)
mytheme_noborder <- theme_bw() + 
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = 9, colour = 'black'),
    axis.text = element_text(size = 9, colour = 'black'),
    axis.ticks = element_line(size = 0.4, colour = 'black'),
    axis.ticks.length = unit(0.1, 'cm'),
    axis.line.x = element_line(size = 0.4, colour = 'black'),
    axis.line.y = element_line(size = 0.4, colour = 'black')
  )

mytheme <- theme_bw() + 
  theme(
    panel.border = element_rect(size = 0.4, colour = 'black'),
    panel.grid = element_blank(),
    axis.title = element_text(size = 10, colour = 'black'),
    axis.text = element_text(size = 10, colour = 'black'),
    axis.ticks = element_line(size = 0.4, colour = 'black'),
    axis.ticks.length = unit(0.1, 'cm')
  )
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig3")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig3"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
excel_file <- "PBMC_barplot.xlsx" 

df_healthy <- read_excel(excel_file, sheet = "frequency_H") %>%
  rename(codon = name, frequency = H) %>%  # 重命名列以匹配原代码
  mutate(group = "frequency_H")

df_diarrhea <- read_excel(excel_file, sheet = "frequency_W") %>%
  rename(codon = name, frequency = W) %>%  # 重命名列以匹配原代码
  mutate(group = "frequency_W")

#单张图片
df_combined <- full_join(df_healthy, df_diarrhea, by = "codon", suffix = c("_healthy", "_diarrhea")) %>%
  mutate( diff = frequency_diarrhea - frequency_healthy,
    color = ifelse(diff > 0, "Up", "Down"),
    third_base = substr(codon, 3, 3),
    category = ifelse(third_base %in% c("A", "T"), "AU3", "GC3")  )
# 绘制带正负值的柱状图
ggplot(df_combined, aes(x = codon, y = diff, fill = color)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  # 添加零线
  scale_fill_manual(values = c("Up" = "#E64B35", "Down" = "#3182BD")) +  # 设置正负色
  facet_wrap(~category, scales = "free_x", nrow = 1) +  # 按AU3/GC3分面
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    x = "Codon",
    y = "Difference (Diarrhea - Healthy)",
    title = "Codon Usage Difference Between Healthy and Diarrhea Groups"
  )

#组图
#健康
df_healthy <- read_excel(excel_file, sheet = "RSCU_H") %>% 
  rename(codon = name, frequency = H) %>%  
  mutate(
    third_base = substr(codon, 3, 3),  # 提取第3位碱基
    category = ifelse(third_base %in% c("A", "T"), "AU3", "GC3"),
    diff = frequency - 0.95  # 计算差异值
  ) %>%
  arrange(diff) %>%  # 按差异值从小到大排序
  mutate(codon = factor(codon, levels = codon)) 

pdf(file.path(fig_out_dir, "PBMC_barplot_healthy.pdf"), width = 18/2.54, height = 7/2.54)
p_healthy_bar <- ggplot(df_healthy, aes(x = codon, y = diff, fill = category)) +
  geom_bar(stat = "identity", width = 0.6) + 
  #geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.4) +  # 零线
  scale_fill_manual(values = c("AU3" = "#E64B35", "GC3" = "#3182BD")) +  # 原图配色
  mytheme +
  labs(x = "Codon",y = "Difference (Healthy - Baseline)",title = "Healthy Group - Codon Usage (Sorted by Diff)" ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "top")
print(p_healthy_bar)
dev.off()

# 绘制健康组 - 箱线图（保持原逻辑不变）
df_healthy_box <- df_healthy %>% 
  group_by(category) %>% 
  summarise(
    value = list(diff)  # 用差异值做箱线图
  ) %>% 
  unnest(value)
library(ggpubr)
df_healthy_box <- df_healthy %>% 
  group_by(category) %>% 
  summarise(value = list(diff)   ) %>% 
  unnest(value)

pdf(file.path(fig_out_dir, "PBMC_boxplot_healthy.pdf"), width =3.5/2.54, height = 4.5/2.54)
p_healthy_box <- ggplot(df_healthy_box, aes(x = category, y = value, fill = category)) +
  geom_boxplot(width = 0.6, outlier.shape = 21,outlier.size = 1, 
  outlier.fill = NA,   color = "black",  size = 0.3 ) +
  scale_fill_manual(values = c("AU3" = "#E64B35", "GC3" = "#3182BD")) + 
  stat_compare_means(
    method = "wilcox.test",  
    label = "p.format",    
    label.x = 1.5,  
    label.y = max(df_healthy_box$value) * 1.1  ) +
  mytheme +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 0.5),
    legend.position = "none")
print(p_healthy_box)
dev.off()

