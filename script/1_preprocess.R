#载入的包
library(tidyverse)
library(Seurat)
#10X Genomics 单细胞测序标准下游分析流程
#########################################
# 1_PREPROCESS,进行文件的分组
fs=list.files('./data/GSE135927_RAW/','^GSM')
fs
# Split up a string into pieces;simplify = T 返回矩阵
samples=str_split(fs,'_',simplify = T)[,1]
lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  #为每个样本创建子文件夹
  folder=paste0("./data/GSE135927_RAW/", str_split(y[1],'_',simplify = T)[,1])
  dir.create(folder,recursive = T)
  #重命名文件，并移动到相应的子文件夹里
  # file.rename(from,to),将A文件移动到另外文件夹中并且重新命名为B
  file.rename(paste0("./data/GSE135927_RAW/",y[1]),file.path(folder,"barcodes.tsv.gz"))
  file.rename(paste0("./data/GSE135927_RAW/",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("./data/GSE135927_RAW/",y[3]),file.path(folder,"matrix.mtx.gz"))
})
samples<-list.files("data/GSE135927_RAW/")
samples
