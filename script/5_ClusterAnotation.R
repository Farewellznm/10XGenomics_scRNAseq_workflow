##########################################
# 5.ClusterAnotation 细胞类型注释
# 有很多方法，可以用 singleR 自动化注释
load("result/sce_output_merge.Rdata")
library(SingleR)
library(Seurat)
ref.se <- celldex::HumanPrimaryCellAtlasData()
sce_for_SingleR <- GetAssayData(sce, slot="data")
clusters <- sce@meta.data$seurat_clusters
pred.hesc <- SingleR(test = sce_for_SingleR, ref = ref.se,
                     labels = ref.se$label.fine, clusters = clusters,
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
table(pred.hesc$labels)
plotScoreHeatmap(pred.hesc)
#免疫基因注释
ref.se_1<-celldex::ImmGenData()
pred.hesc_1<- SingleR(test = sce_for_SingleR, ref = ref.se_1,
                     labels = ref.se_1$label.fine, clusters = clusters,
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
table(pred.hesc_1$labels)
plotScoreHeatmap(pred.hesc_1)

#  29 fine cell types
ref.se_2<-celldex::MonacoImmuneData()
pred.hesc_2<- SingleR(test = sce_for_SingleR, ref = ref.se_2,
                      labels = ref.se_2$label.fine, clusters = clusters,
                      assay.type.test = "logcounts", assay.type.ref = "logcounts")
table(pred.hesc_2$labels)
plotScoreHeatmap(pred.hesc_2)

# tsne可视化

celltype<-data.frame(ClusterID =rownames(pred.hesc_2),celltype = pred.hesc_2$labels,stringsAsFactors = F)
#如下为sce对象注释细胞cluster鉴定结果。
sce@meta.data$celltype = "NA"
#先新增列celltype，值均为NA，然后利用下一行代码循环填充
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
DimPlot(sce, group.by="celltype", label=T, label.size= 4, reduction='tsne')
# 可以用ggplot2 来绘制注释好的细胞分群图，详见5-1