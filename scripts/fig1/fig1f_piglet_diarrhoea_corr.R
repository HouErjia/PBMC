rm(list = ls())
library(psych)
library(tidyr)
library(dplyr)
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig1")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig1"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
data=read.csv("./piglet_diarrhoea.csv",row.names = 1)
df=data.frame(t(data[,-11]))
res=corr.test(df,df,use="pairwise",
              method = "spearman",adjust = "holm",alpha = 0.05)

cor_df <- as.data.frame(res$r)
cor_df[lower.tri(cor_df, diag = TRUE)] <- NA
cor_df$cluster=data$cluster[rownames(cor_df)%in%rownames(data)]

pvalue_df <- as.data.frame(res$p)
pvalue_df[lower.tri(pvalue_df, diag = TRUE)] <- NA
pvalue_df$cluster=data$cluster[rownames(pvalue_df)%in%rownames(data)]


# 将相关系数和p值转换为长格式
long_cor=cor_df %>%
  as_tibble(rownames = "factor_1")%>%
  pivot_longer(cols = -c(factor_1,cluster), names_to = "factor_2", values_to = "correlation",values_drop_na = T)

long_pvalue <- pvalue_df %>%
  as_tibble(rownames = "factor_1") %>%
  pivot_longer(cols = -c(factor_1,cluster), names_to = "factor_2", values_to = "p_value",values_drop_na = T)






# 合并到一个数据框
df <- bind_cols(long_cor, p_value=long_pvalue$p_value)
df$padj=p.adjust(df$p_value,method = "fdr")
df=df[df$padj<0.05 & abs(df$correlation)>0.4,]
# write.table(df,"./corr_piglet_diarrhoea.txt",row.names = F)
write.csv(df,"./corr_piglet_diarrhoea.csv",row.names = F)


