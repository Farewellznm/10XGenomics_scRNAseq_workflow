#########################################
#2_CreatSeuratObj 创建Seurat对象
scelist<-lapply(samples,function(pro){
  folder<-file.path("data/GSE135927_RAW/",pro)
  CreateSeuratObject(counts = Read10X(folder),project = pro)
})
#合并两组Seurat对象
sce.big <- merge(scelist[[1]],scelist[[2]],project = "ls_12",add.cell.ids =samples)
sce.big
table(sce.big$orig.ident)
#保存Seurat对象
save(sce.big,file = "result/sce.big.merge.ls_12.Rdata")
