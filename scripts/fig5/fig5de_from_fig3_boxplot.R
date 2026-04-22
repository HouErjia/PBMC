rm(list = ls())


setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig5")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig5"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
library(dplyr)
library(ggplot2)
library(readxl)
library(ggpubr)

mytheme_noborder <- theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = 9, colour = "black"),
    axis.text = element_text(size = 9, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    axis.line.x = element_line(linewidth = 0.4, colour = "black"),
    axis.line.y = element_line(linewidth = 0.4, colour = "black")
  )

mytheme <- theme_bw() +
  theme(
    panel.border = element_rect(linewidth = 0.4, colour = "black"),
    panel.grid = element_blank(),
    axis.title = element_text(size = 10, colour = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    plot.title = element_text(hjust = 0.5)
  )

input_file <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig5/fig5de.xlsx"
output_dir <- dirname(input_file)
plot_width_cm <- 5.6
plot_height_cm <- 7.2

resolve_sheet_name <- function(candidates, workbook_sheets) {
  matched <- candidates[candidates %in% workbook_sheets]
  if (length(matched) == 0) {
    stop(sprintf("None of the sheets were found: %s", paste(candidates, collapse = ", ")))
  }
  matched[1]
}

read_metric_data <- function(sheet_name) {
  metric_df <- as.data.frame(read_excel(input_file, sheet = sheet_name))
  if (!all(c("H", "W") %in% names(metric_df))) {
    stop(sprintf("Sheet %s must contain columns named H and W.", sheet_name))
  }

  data.frame(
    group = factor(
      c(rep("H", sum(!is.na(metric_df$H))), rep("W", sum(!is.na(metric_df$W)))),
      levels = c("H", "W")
    ),
    value = c(metric_df$H[!is.na(metric_df$H)], metric_df$W[!is.na(metric_df$W)])
  )
}

build_boxplot <- function(plot_df, fill_colors, y_label, title_text) {
  value_range <- diff(range(plot_df$value, na.rm = TRUE))
  signif_y <- max(plot_df$value, na.rm = TRUE) + value_range * 0.02

  ggplot(plot_df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = 21,
      outlier.size = 1,
      outlier.fill = NA,
      color = "black",
      linewidth = 0.3
    ) +
    scale_fill_manual(values = fill_colors) +
    stat_compare_means(
      comparisons = list(c("H", "W")),
      aes(label = after_stat(p.signif)),
      method = "wilcox.test",
      label.y = signif_y,
      tip.length = 0.02
    ) +
    mytheme +
    labs(x = NULL, y = y_label, title = title_text) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none"
    )
}

sheet_names <- excel_sheets(input_file)
gc3_sheet <- resolve_sheet_name(c("CG3", "GC3"), sheet_names)
cai_sheet <- resolve_sheet_name(c("CAI"), sheet_names)

gc3_df <- read_metric_data(gc3_sheet)
cai_df <- read_metric_data(cai_sheet)

group_colors <- c("H" = "#4B5EA0", "W" = "#CA3F4C")

gc3_plot <- build_boxplot(
  plot_df = gc3_df,
  fill_colors = group_colors,
  y_label = "GC3 Comparison",
  title_text = "GC3"
)

cai_plot <- build_boxplot(
  plot_df = cai_df,
  fill_colors = group_colors,
  y_label = "CAI Distribution",
  title_text = "CAI"
)

pdf(file.path(output_dir, "fig5de_GC3_from_fig3_boxplot.pdf"), width = plot_width_cm / 2.54, height = plot_height_cm / 2.54)
print(gc3_plot)
dev.off()

pdf(file.path(output_dir, "fig5de_CAI_from_fig3_boxplot.pdf"), width = plot_width_cm / 2.54, height = plot_height_cm / 2.54)
print(cai_plot)
dev.off()

stats_df <- data.frame(
  Metric = c("GC3", "CAI"),
  Sheet = c(gc3_sheet, cai_sheet),
  H_n = c(sum(gc3_df$group == "H"), sum(cai_df$group == "H")),
  W_n = c(sum(gc3_df$group == "W"), sum(cai_df$group == "W")),
  Wilcoxon_P = c(
    wilcox.test(value ~ group, data = gc3_df, exact = FALSE)$p.value,
    wilcox.test(value ~ group, data = cai_df, exact = FALSE)$p.value
  )
)

write.csv(
  stats_df,
  file = file.path(output_dir, "fig5de_from_fig3_boxplot_stats.csv"),
  row.names = FALSE
)

print(stats_df)

