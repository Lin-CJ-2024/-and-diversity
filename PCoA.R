rm(list=ls())#clear Global Environment
setwd('')


# install.packages("vegan")
# install.packages("ggplot2")

library(vegan)
library(ggplot2)


otu_raw <- read.table(file="data/otu.txt",sep="\t",header=T,check.names=FALSE ,row.names=1)

otu <- t(otu_raw)


otu.distance <- vegdist(otu)



otu.distance1 <- as.matrix(round(otu.distance,digits = 3))
write.table(otu.distance1,'output_data/matrix_bray_curtis.txt', row.names = T, sep = '\t', quote = F)

pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)



pc12 <- as.data.frame(pc12)

pc12$samples <- row.names(pc12)
head(pc12)


p <- ggplot(pc12,aes(x=V1, y=V2))+
  geom_point(size=3)+
  theme_bw()
p


group <- read.table("data/group.txt", sep='\t', header=T)

colnames(group) <- c("samples","group")

df <- merge(pc12,group,by="samples")
head(df)

color=c("#34bf49","#ff4c4c","#0099e5","#FFB90F","#9370DB")
p1 <-ggplot(data=df,aes(x=V1,y=V2,
                       color=group,shape=group))+
  theme_bw()+
  geom_point(size=1.8)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  #geom_text(aes(label=samples, y=V2+0.03,x=V1+0.03,  vjust=0),size=3.5)+
 # guides(color=guide_legend(title=NULL))+
  labs(x=paste0("PC1 ",pc[1],"%"),
       y=paste0("PC2 ",pc[2],"%"))+
  stat_ellipse(data=df,
               geom = "polygon",level=0.9,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(values = color) +
  scale_fill_manual(values = c("#34bf49","#ff4c4c","#0099e5","#FFB90F","#9370DB"))+
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank())
p1
