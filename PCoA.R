rm(list=ls())#clear Global Environment
setwd('C:/Users/user/Desktop/α-diversity/马铃薯 BOLD/扩增子测序数据下游分析合集')#设置工作路径

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

##########PCoA图添加箱线图######
#在PCA图的x和y轴添加箱线图，可以实现进一步展示组间差异
colnames(df) <- c("samples", "PC1", "PC2", "group")#修改数据列名
#加载包，对组间进行统计检验以及组合图的拼接
library(ggpubr)
library(ggsignif)
# 绘制y轴为PC2值的分组箱线图
p2 <- ggplot(df,aes(x=group,y=PC2))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,linewidth=0.5)+#添加误差线,注意位置，放到最后则这条线不会被箱体覆盖
  geom_boxplot(aes(fill=group), #绘制箱线图函数
               outlier.colour="white",size=0.5)+#异常点去除
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(color = "white"),#坐标轴的线设为显示
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),#关闭刻度
        legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values=c("#34bf49","#ff4c4c","#0099e5","#FFB90F","#9370DB"))+#指定颜色
  geom_signif(comparisons = list(c("M1", "M2"), c("M1", "M3"), c("M1", "M4"), c("M1", "M5"), c("M2", "M3"), c("M2", "M4"),
                                 c("M2","M5"),c("M3","M5"),c("M4","M5")),# 设置需要比较的组
              map_signif_level = T, #是否使用星号显示
              test = t.test, ##计算方法
              y_position = c(1.1,1.2,1.3,1.5,1.6,1.7,1.8,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4),#图中横线位置 设置
              tip_length = c(c(0,0),
                             c(0,0),
                             c(0,0),
                             c(0,0),
                             c(0,0),
                             c(0,0)),#横线下方的竖线设置
              size=0.8,color="black")
p2
# 绘制y轴为PC1值的分组箱线图
p3 <- ggplot(df,aes(x=group,y=PC1))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,linewidth=0.5)+#添加误差线,注意位置，放到最后则这条线不会被箱体覆盖
  coord_flip()+
  geom_boxplot(aes(fill=group), #绘制箱线图函数
               outlier.colour="white",size=0.5)+#异常点去除
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(color = "white"),#坐标轴的线设为显示
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),#关闭刻度
        legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values=c("#34bf49","#ff4c4c","#0099e5","#FFB90F"))+#指定颜色
  geom_signif(comparisons = list(c("A", "B"), c("A", "C"), c("A", "D"), c("B", "C"), c("B", "D"), c("C", "D")),# 设置需要比较的组
              map_signif_level = T, #是否使用星号显示
              test = t.test, ##计算方法
              y_position = c(0.48,0.55,0.62,0.73,0.80,0.48),#图中横线位置 设置
              tip_length = c(c(0,0),
                             c(0,0),
                             c(0,0),
                             c(0,0),
                             c(0,0),
                             c(0,0)),#横线下方的竖线设置
              size=0.8,color="black")
p3
# ggpubr::ggarrange()函数对图进行拼接
ggarrange(p3, NULL, p1, p2, widths = c(5,2), heights = c(2,5), align = "hv")
#保存图片
ggsave('fig/PCoA.pdf',width=12,height = 10)
