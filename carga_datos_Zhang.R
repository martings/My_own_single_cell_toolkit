library(Seurat)
setwd("/media/maxpower/DatosII/Fededa_Lab/")

pbmc<-ReadMtx(mtx = "GSE169246_TNBC_RNA.counts.mtx.gz", features = "GSE169246_TNBC_RNA.feature.tsv.gz",
                cells = "GSE169246_TNBC_RNA.barcode.tsv.gz", feature.column = 1)

pbmc <- CreateSeuratObject(counts = pbmc,project ="Zhang_21")

pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(object = pbmc)

pbmc <- RunPCA(object = pbmc, verbose = T)
pbmc <- RunUMAP(object = pbmc, dims = 1:30)


plot <- DimPlot(pbmc, reduction = "umap", group.by = "orig.ident" )
plot

pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.1)

plot <- DimPlot(pbmc, reduction = "umap" )
plot

saveRDS(pbmc,"Zhang_21_whole_cells_RNA.rds")

