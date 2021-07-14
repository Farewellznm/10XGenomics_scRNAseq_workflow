##########################################
# based on 4_DimensionaalityReductionClustering 
# using ggplot2 for tsne clustering
library(ggplot2)
phe<-data.frame(rownames(sce@meta.data),cluster = sce@meta.data$seurat_clusters)
tsne_pos = Embeddings(sce,"tsne")
dat=cbind(tsne_pos,phe)
head(dat)

p=ggplot(dat,aes(x=tSNE_1,y=tSNE_2,color=cluster))+geom_point(size= 1)
p=p+stat_ellipse(data=dat,aes(x=tSNE_1,y=tSNE_2,fill=cluster,color=cluster),
                 geom = "polygon",alpha=0.2,level=0.9,type="t",linetype =2,show.legend = F)+coord_fixed()
print(p)

theme= theme(panel.grid =element_blank()) +   ## 删去网格
  theme(panel.border = element_blank(),panel.background = element_blank()) +   ## 删去外层边框
  theme(axis.line = element_line(size=1, colour = "black"))
p+theme+guides(colour = guide_legend(override.aes = list(size=5)))
