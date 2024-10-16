rm(list=ls())#clear Global Environment
setwd('')#设置工作路径

#安装所需R包
# install.packages("vegan")
# install.packages("ggplot2")
#加载包
library(vegan)#计算距离时需要的包
library(ggplot2)#绘图包

#读取数据，一般所需是数据行名为样本名、列名为OTUxxx的数据表
otu_raw <- read.table(file="data/otu.txt",sep="\t",header=T,check.names=FALSE ,row.names=1)
#由于排序分析函数所需数据格式原因，需要对数据进行转置
otu <- t(otu_raw)

#计算bray_curtis距离
otu.distance <- vegdist(otu)


#导出表格
otu.distance1 <- as.matrix(round(otu.distance,digits = 3))
write.table(otu.distance1,'output_data/matrix_bray_curtis.txt', row.names = T, sep = '\t', quote = F)
#pcoa分析
pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)#解释度

###绘图###
#pcl2原来是matrix,转化为data.frame
pc12 <- as.data.frame(pc12)
#给pc12添加samp1es变量
pc12$samples <- row.names(pc12)
head(pc12)

#绘图
p <- ggplot(pc12,aes(x=V1, y=V2))+#指定数据、X轴、Y轴
  geom_point(size=3)+#绘制点图并设定大小
  theme_bw()#主题
p

#读入分组文件
group <- read.table("data/group.txt", sep='\t', header=T)
#修改列名
colnames(group) <- c("samples","group")
#将绘图数据和分组合并
df <- merge(pc12,group,by="samples")
head(df)
#绘图
color=c("#34bf49","#ff4c4c","#0099e5","#FFB90F","#9370DB")#颜色变量
p1 <-ggplot(data=df,aes(x=V1,y=V2,
                       color=group,shape=group))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(size=1.8)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+#图中虚线
  #geom_text(aes(label=samples, y=V2+0.03,x=V1+0.03,  vjust=0),size=3.5)+#添加数据点的标签
 # guides(color=guide_legend(title=NULL))+#去除图例标题
  labs(x=paste0("PC1 ",pc[1],"%"),
       y=paste0("PC2 ",pc[2],"%"))+#将x、y轴标题改为贡献度
  stat_ellipse(data=df,
               geom = "polygon",level=0.9,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(values = color) +#点的颜色设置
  scale_fill_manual(values = c("#34bf49","#ff4c4c","#0099e5","#FFB90F","#9370DB"))+
  theme(axis.title.x=element_text(size=12),#修改X轴标题文本
        axis.title.y=element_text(size=12,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=10),#修改x轴刻度标签文本
        axis.text.x=element_text(size=10),#修改y轴刻度标签文本
        panel.grid=element_blank())#隐藏网格线
p1
