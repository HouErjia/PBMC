rm(list=ls())

library(ropls) 
library(ggforce)
library(ggplot2) 
library(ggprism) 
library(xlsx)

options(java.parameters = "-Xmx4g")
setwd("D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/data_processed/fig1")
fig_out_dir <- "D:/PBMC交接/PBMC_sqluo/PBMC_paper_package/results/figures/fig1"
dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)
#meta<-read.xlsx("./correlation_RNAseq_DC.xlsx",1,header=T,row.names=1)
meta<-read.csv("./fig1c_rpkm_values.csv",header=T,row.names=1)
meta=t(meta)
meta=as.data.frame(meta)
meta=meta[c(1:10),]
meta=meta[,which(colSums(meta)>0)]
meta=meta[which(rowSums(meta)>0),]

group=data.frame(group=c(rep("H",5),rep("W",5)))
row.names(group) <- row.names(meta)

df1_oplsda <- opls(meta, predI = 2, orthoI = 0,crossvalI = 10)

#提取坐标值
data <- as.data.frame(df1_oplsda@scoreMN)
o1 <- df1_oplsda@orthoScoreMN[,1]
data$o1 <- o1
data$group = group$group
data$samples = rownames(data)
#提取解释度
x_lab <- df1_oplsda@modelDF[1, "R2X"] * 100
y_lab <- df1_oplsda@modelDF[2, "R2X"] * 100
#绘图
col=c("H"="#4e62ab","W"="#d6404e")

p1 <- ggplot(data,aes(x=p1,y=p2,color=group,))+
  theme_bw()+#主题设置
  geom_point(size=2)+#绘制点图并设定大小
  theme(panel.grid = element_blank(),
  plot.title = element_text(size = 10, hjust = 0.5)) +
  geom_vline(xintercept = 0,lty="dashed",color="grey")+
  geom_hline(yintercept = 0,lty="dashed",color="grey")+#图中虚线
  labs(x=paste0("P1 (",x_lab,"%)"),
       y=paste0("P2(",y_lab,"%)"))+#将x、y轴标题改为贡献度
  stat_ellipse(data=data,
               geom = "polygon",level = 0.95,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(values = col) +#点的颜色设置
  scale_fill_manual(values = col)+
  theme(axis.title.x=element_text(size=10),#修改X轴标题文本
        axis.title.y=element_text(size=10,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=10),#修改x轴刻度标签文本
        axis.text.x=element_text(size=10),#修改y轴刻度标签文本
        panel.grid=element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5))+#隐藏网格线
  ggtitle("PCA")

p1
library(ggplot2)

ggsave(file.path(fig_out_dir, "RNAseq_PCA.pdf"),units = "cm",width = 10.5,height = 8)


