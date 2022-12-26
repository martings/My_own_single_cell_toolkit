setwd("/media/maxpower/DatosII/Fededa_Lab/")

library(Seurat)
library(tidyverse)
pbmc<-readRDS("Zhang_21_whole_cells_RNA.rds")


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5_x_clus<-pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write_csv(x = top5_x_clus, file = "top5_clusters_por_integracion.csv")
write_csv(x = pbmc.markers, file = "AllMarkers_por_clusters.csv")

setwd("Resultados/")
#Este paso se puede obviar porque ya viene normalizado por cellranger count
#Features plot umap
for ( i in levels(top5_x_clus$cluster)) {
  nombre<-paste0("Cluster_",i)
  genes<-subset(top5_x_clus, subset = cluster%in%i)[,7]
  plot1<-FeaturePlot(pbmc, reduction = "umap",features = genes$gene,  min.cutoff = "q10", max.cutoff = "q90" )
  plot1 + labs(title = nombre)
  plot1
  ggsave(paste0(nombre,"_umap_hd_presentacion.png"), plot = plot1,height = 8.27 ,width = 8.27)
}

dot<-DotPlot(pbmc, features = top5_x_clus$gene)+ RotatedAxis()
ggsave("Dotplot_top5_Total.pdf", plot = dot,height = 8.27 ,width = 18.27)
#Al pedo, no se ve
#VlnPlot(pbmc.integrated, features = top5_x_clus$gene)+ RotatedAxis()

for ( i in levels(top5_x_clus$cluster)) {
  nombre<-paste0("Cluster_",i)
  genes<-subset(top5_x_clus, subset = cluster%in%i)[,7]
  plot1<-VlnPlot(pbmc, features = genes$gene)
  plot1 + labs(title = nombre)
  plot1
  ggsave(paste0("VlnPlot_",nombre,"_cluster_top5_umap.png"), plot = plot1,height = 8.27 ,width =  11.69)
}


plot1<-DimPlot(pbmc, reduction = "umap", label = T, label.size = 6)
ggsave("Zhang_21_whole_cells_RNA_UMAP.png" ,plot = plot1 ,height = 8.27, width =11.69)

#Marcadores B
B_Naive<-c("CD19", "MS4A1", "FCER2", "TCL1A")
B_Foliculares<-c("CD19", "MS4A1", "NEIL1", "RGS13", "MEF2B", "BCL6")
B_Mem<-c("CD19", "MS4A1","CD27", "TNFRSF13B")
B_Plas<-  c("CD79A","MZB1","IGHG1","XBP1")
  


plot1<-FeaturePlot(pbmc, features = B_Naive, min.cutoff = "p10", raster = F)
ggsave("B_Naive_Zhang_21_whole_cells_RNA_UMAP.png" ,plot = plot1 ,height = 8.27, width =11.69)


plot1<-FeaturePlot(pbmc, features = B_Foliculares, min.cutoff = "p10", raster = F)
ggsave("B_Foliculares_Zhang_21_whole_cells_RNA_UMAP.png" ,plot = plot1 ,height = 8.27, width =11.69)

plot1<-FeaturePlot(pbmc, features = B_Mem, min.cutoff = "p10", raster = F)
ggsave("B_Mem_Zhang_21_whole_cells_RNA_UMAP.png" ,plot = plot1 ,height = 8.27, width =11.69)

plot1<-FeaturePlot(pbmc, features = B_Plas, min.cutoff = "p10", raster = F)
ggsave("B_Plas_Zhang_21_whole_cells_RNA_UMAP.png" ,plot = plot1 ,height = 8.27, width =11.69)

#Subset dataset por cluster para quedarme solo con las celulas B
#Creo por ahora que sertian los cluster 7, 3, 6


pbmc_sub<-subset(pbmc, idents = c("3", "6", "7"))

pbmc_sub<-SCTransform(pbmc_sub)
pbmc_sub<-RunPCA(pbmc_sub, dims= 1:20)
pbmc_sub<-RunUMAP(pbmc_sub, dims= 1:20)
pbmc_sub<-FindNeighbors(pbmc_sub, dims= 1:20)
pbmc_sub<-FindClusters(pbmc_sub, resolution = 0.1)
DimPlot(pbmc_sub, reduction = "umap")
DimPlot(pbmc_sub, reduction = "umap")


#Marcadores B
B_Naive<-c("CD19", "MS4A1", "FCER2", "TCL1A")
B_Foliculares<-c("CD19", "MS4A1", "NEIL1", "RGS13", "MEF2B", "BCL6")
B_Mem<-c("CD19", "MS4A1","CD27", "TNFRSF13B")
B_Plas<-  c("CD79A","MZB1","IGHG1","XBP1")



plot1<-FeaturePlot(pbmc_sub, features = B_Naive, min.cutoff = "p10", raster = F)
ggsave("Subset_2_B_Naive_Zhang_21_whole_cells_RNA_UMAP.png" ,plot = plot1 ,height = 8.27, width =11.69)


plot1<-FeaturePlot(pbmc_sub, features = B_Foliculares, min.cutoff = "p10", raster = F)
ggsave("Subset_2_B_Foliculares_Zhang_21_whole_cells_RNA_UMAP.png" ,plot = plot1 ,height = 8.27, width =11.69)

plot1<-FeaturePlot(pbmc_sub, features = B_Mem, min.cutoff = "p10", raster = F)
ggsave("Subset_2_B_Mem_Zhang_21_whole_cells_RNA_UMAP.png" ,plot = plot1 ,height = 8.27, width =11.69)

plot1<-FeaturePlot(pbmc_sub, features = B_Plas, min.cutoff = "p10", raster = F)
ggsave("Subset_2_B_Plas_Zhang_21_whole_cells_RNA_UMAP.png" ,plot = plot1 ,height = 8.27, width =11.69)

saveRDS(pbmc_sub, "Subcluster_3_6_7_B_Cells_y_otras.rds")



pbmc_sub<-subset(pbmc, idents = c("3","7"))

pbmc_sub<-SCTransform(pbmc_sub)
pbmc_sub<-RunPCA(pbmc_sub, dims= 1:20)
pbmc_sub<-RunUMAP(pbmc_sub, dims= 1:20)
pbmc_sub<-FindNeighbors(pbmc_sub, dims= 1:20)
pbmc_sub<-FindClusters(pbmc_sub, resolution = 0.1)
DimPlot(pbmc_sub, reduction = "umap")
DimPlot(pbmc_sub, reduction = "umap")


#Marcadores B
B_Naive<-c("CD19", "MS4A1", "FCER2", "TCL1A")
B_Foliculares<-c("CD19", "MS4A1", "NEIL1", "RGS13", "MEF2B", "BCL6")
B_Mem<-c("CD19", "MS4A1","CD27", "TNFRSF13B")
B_Plas<-  c("CD79A","MZB1","IGHG1","XBP1")



plot1<-FeaturePlot(pbmc_sub, features = B_Naive, min.cutoff = "p10", raster = F)
ggsave("Subset_3_B_Naive_Zhang_21_whole_cells_RNA_UMAP.png" ,plot = plot1 ,height = 8.27, width =11.69)


plot1<-FeaturePlot(pbmc_sub, features = B_Foliculares, min.cutoff = "p10", raster = F)
ggsave("Subset_3_B_Foliculares_Zhang_21_whole_cells_RNA_UMAP.png" ,plot = plot1 ,height = 8.27, width =11.69)

plot1<-FeaturePlot(pbmc_sub, features = B_Mem, min.cutoff = "p10", raster = F)
ggsave("Subset_3_B_Mem_Zhang_21_whole_cells_RNA_UMAP.png" ,plot = plot1 ,height = 8.27, width =11.69)

plot1<-FeaturePlot(pbmc_sub, features = B_Plas, min.cutoff = "p10", raster = F)
ggsave("Subset_3_B_Plas_Zhang_21_whole_cells_RNA_UMAP.png" ,plot = plot1 ,height = 8.27, width =11.69)

saveRDS(pbmc_sub, "Subcluster_3_y_7_B_Cells.rds")

