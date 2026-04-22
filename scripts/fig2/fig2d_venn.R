rm(list = ls())
library(ggvenn)
library(tidyverse)
library(ggtext)
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig2")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig2"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
#导入基因名文件，A,B分别只含一列基因名
data_diff <- "./fig2d_venn.xlsx"
library(xlsx)
Ribo <- read.xlsx(data_diff, 1,head=TRUE)
RNAseq<- read.xlsx(data_diff, 2,head=TRUE)

#将文件合并为list
 list=list("Riboseq"=Ribo$transcripts, "RNAseq"=RNAseq$transcripts)
#出图
p1=ggvenn(list,show_percentage = F,show_elements = F,label_sep = ",",
            digits = 1,stroke_color = NA,
            fill_color = c("#dd7208","#fbca50","lightblue"),
            set_name_size = 5,
            set_name_color = c("#dd7208","#fbca50","lightblue"),
            text_size = 5)
p1
ggsave(file.path(fig_out_dir, "venn_input.pdf"),units = "in",width=10/2.54,height=10/2.54)
#导出交集
intersection=intersect(Ribo$transcripts,RNAseq$transcripts)
result <- data.frame(Intersection = intersection)
write.xlsx(result, file = "./venn_intersection.xlsx")
#导出去交集
difference <- setdiff(Ribo$transcripts,RNAseq$transcripts)
result <- data.frame(Difference = difference)
write.xlsx(result, file = "./venn_de_in.xlsx")


