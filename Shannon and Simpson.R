

setwd("") 
rm(list = ls())


library(reshape2)
library(dplyr)
library(stringr)
library(vegan)
library(ggplot2)
library(ggpubr)


genus <- read.delim('otu.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
group <- read.delim("group.txt",header = T)



Shannon <- diversity(genus, index = "shannon", MARGIN = 2, base = exp(1))   #MARGIN：1用列数据算，2用行数据算

Simpson <- diversity(genus, index = "simpson", MARGIN = 2, base = exp(1))

Richness <- specnumber(genus, MARGIN = 2) 


index <- as.data.frame(cbind(Shannon, Simpson, Richness))


tgenus<-ceiling(as.data.frame(t(genus))) 
obs_chao_ace <- t(estimateR(tgenus)) 


index$Chao <- obs_chao_ace[,2]
index$Ace <- obs_chao_ace[,4]
index$obs <- obs_chao_ace[,1]


index$Pielou <- Shannon / log(Richness) 


index$sample<- c(rownames(index))


data <- merge(index,group,by = 'sample')


write.csv(data,"5-alpha data.csv")


shapiro.test(Shannon)

Shannon<- ggplot(data,aes(x=group,y=Shannon))+  
  stat_boxplot(geom = "errorbar", linewidth=0.2)+ 
  geom_boxplot(aes(fill=group),linewidth=0.2)+ 
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

