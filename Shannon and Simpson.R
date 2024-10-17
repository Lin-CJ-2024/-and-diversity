
#设置文件路径
setwd("") 
rm(list = ls())

#导入需要的包
library(reshape2)
library(dplyr)
library(stringr)
library(vegan)
library(ggplot2)
library(ggpubr)

#读取物种数据和分组数据
genus <- read.delim('otu.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
group <- read.delim("group.txt",header = T)

###计算多样性指数
#Shannon指数
Shannon <- diversity(genus, index = "shannon", MARGIN = 2, base = exp(1))   #MARGIN：1用列数据算，2用行数据算
#simpson指数
Simpson <- diversity(genus, index = "simpson", MARGIN = 2, base = exp(1))
#物种丰富度
Richness <- specnumber(genus, MARGIN = 2) 

#将以上三个指数先统计成表
index <- as.data.frame(cbind(Shannon, Simpson, Richness))

#计算obs，chao，ace指数
tgenus<-ceiling(as.data.frame(t(genus)))  #转置物种数据之后取整数（向上取整）
obs_chao_ace <- t(estimateR(tgenus))  #estimateR获取obs，chao，ace指数

#将obs，chao，ace指数与前面指数计算结果合并
index$Chao <- obs_chao_ace[,2]
index$Ace <- obs_chao_ace[,4]
index$obs <- obs_chao_ace[,1]

#计算Pielou指数
index$Pielou <- Shannon / log(Richness)  #这里使用自然对数ln，部分文献会使用log10和log2

#加一列sample以便后续匹配
index$sample<- c(rownames(index))

#合并分组信息与多样性指数
data <- merge(index,group,by = 'sample')

#导出alpha多样性指数
write.csv(data,"5-alpha data.csv")

###以Shannon指数为例画箱线图（即上文图一），需要其他指数替换Shannon即可
shapiro.test(Simpson)
#my_comparisons <- list(c("Da", "Sha"), c("Da", "Su"), c("Da", "Xi"), c("Sha", "Xi"), c("Sha", "Su"), c("Su", "Xi"))
Shannon<- ggplot(data,aes(x=group,y=Shannon))+  
  stat_boxplot(geom = "errorbar", linewidth=0.2)+  #添加误差线
  geom_boxplot(aes(fill=group),linewidth=0.2)+  #箱线图
  scale_fill_manual(values =c("#34bf49","#ff4c4c","#0099e5","#FFB90F","#9370DB"))+
  #stat_compare_means(method = "kruskal.test", label.y = 2,label.x = 0.7,label = "p.signif") +
  stat_compare_means(comparisons=list(c("M1", "M2"), c("M1", "M3"), c("M1", "M4"), c("M1", "M5"), c("M2", "M3"), c("M2", "M4"),
                                      c("M2","M5"),c("M3","M5"),c("M4","M5")),
                     method="wilcox.test",label = "p.signif",size=5,vjust = 0.5,hjust=1.3)+theme_bw()+
  theme(panel.grid = element_blank(),
    axis.text.x = element_text(size = 16,face = "plain",hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16,face = "plain",hjust = 0.5),
        axis.title.y = element_text(vjust=0.2,size=16,face="bold"),
        strip.text = element_text(size = 16),
        legend.position ="none")+ylab(NULL)+
  facet_grid(.~"Shannon")
Shannon
# Note:When drawing Simpson, replace "Shannon" with Simpson.

