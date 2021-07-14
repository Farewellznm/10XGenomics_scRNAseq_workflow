##########################################
#3_CellQc 细胞质控
#接下来的分析都是基于Seurat对象
rm(list = ls())
load(file = "result/sce.big.merge.ls_12.Rdata")
raw_sce<-sce.big
#索引线粒体基因
index_mt<-grepl("^mt",rownames(raw_sce),ignore.case = T)
mt_gene<-rownames(raw_sce)[index_mt]
#索引核糖体基因
index_RP<-grepl("^RP[SL]",rownames(raw_sce),ignore.case = T)
RP_gene<-rownames(raw_sce)[index_RP]

#质控前绘图

#线粒体基因百分比
raw_sce[["percent_mt"]] <- PercentageFeatureSet(raw_sce, pattern = "^MT-")
#核糖体基因百分比,这里展示计算的原理
C<-GetAssayData(object = raw_sce, slot = "counts")
percent_RP <- Matrix::colSums(C[RP_gene,])/Matrix::colSums(C)*100
raw_sce <- AddMetaData(raw_sce, percent_RP, col.name = "percent_RP")

plot1 <- FeatureScatter(raw_sce, feature1 = "nCount_RNA", feature2 = "percent_mt")
plot2 <- FeatureScatter(raw_sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
VlnPlot(raw_sce, features = c("percent_RP", "percent_mt"), ncol = 2)
VlnPlot(raw_sce, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#细胞质控
raw_sce1 <- subset(raw_sce, subset = nFeature_RNA > 200  & nCount_RNA > 1000 & percent_mt < 20)
#质控后绘图
VlnPlot(raw_sce1,features = c("nFeature_RNA","nCount_RNA"),ncol =2)
save(raw_sce1,file = "result/sce.Rdata")

