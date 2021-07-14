##########################################
# 4.降维，聚类，可视化
# 降维用PCA，排除共线性，获得贡献率排名靠前的特征
# 聚类用Louvain algorithm 社群发现算法，Louvain 算 法又叫 graph-based 聚类算法，是目前单细胞测序数据分析中细胞分类算法中常用的算
#法，该算法在效率和效果上都表现比较好
# 可视化使用UMAP和t-SNE
rm(list = ls())
load("result/sce.Rdata")
sce<-raw_sce1
#标准化
sce<-NormalizeData(sce,normalization.method = "LogNormalize", scale.factor = 10000)
GetAssay(sce,assay = "RNA")
#前2000个高变feature RNA,vst法，特征方差越大越是高变基因
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
#中心化
sce <- ScaleData(sce)
#主成分分析PCA
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce))
pdf("result/top2000VariableFeaturePlot.pdf")
VariableFeaturePlot(sce)
dev.off()
#top18 PCA
pdf("result/DimHeatmap.pdf")
DimHeatmap(sce, dims = 1:18, cells = 200, balanced = TRUE)
dev.off()
pdf("result/PCA_Elbow.pdf")
ElbowPlot(sce)
dev.off()
# 这里判断为top15
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.6)
table(sce@meta.data$RNA_snn_res.0.6)
#聚类成了8个cluster
#使用t—SNE进行可视化
set.seed(2021)
sce <- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE)
pdf("result/tsne.pdf")
DimPlot(sce,reduction = "tsne",label=T)
dev.off()
#不同样本的细胞聚类情况
table(sce@meta.data$seurat_clusters,sce@meta.data$orig.ident)
pdf("result/tsne_by_orig.ident.pdf")
DimPlot(sce,reduction = "tsne",label = T,split.by = "orig.ident")
dev.off()
#优化的可视化可以见4-1
##########################################
#寻找每个cluster的高变代表基因，并选取前5个，进行可视化
# Marker 基因分析
table(sce@meta.data$seurat_clusters)
#寻找高变基因
p <- list()
for (i in unique(sce@meta.data$seurat_clusters)){
  markers_df <- FindMarkers(object = sce, ident.1 = i, min.pct = 0.25)
  print(x = head(markers_df))
  markers_genes =  rownames(head(x = markers_df, n = 4))
  p1 <- VlnPlot(object = sce, features =markers_genes,log =T,ncol = 2)
  p[[i]][[1]] <- p1
  p2 <- FeaturePlot(object = sce, features=markers_genes,ncol = 2)
  p[[i]][[2]] <- p2
  
}
p[[1]][[1]]
p[[2]][[1]]
p[[3]][[2]]

sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25,
                              thresh.use = 0.25)
DT::datatable(sce.markers)
#热图可视化每个cluster的marker基因表达差异
#热图可视化每个cluster的marker基因表达差异
# dplyr包提供管道函数 %>%
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10,avg_log2FC)
DoHeatmap(sce,top10$gene,size=0.5)+NoLegend()
save(sce,sce.markers,file = 'result/sce_output_merge.Rdata')
