setwd("/media/maxpower/DatosII/Fededa_Lab/")
library(Seurat)
library(tidyverse)

###Todas las celulas
pbmc<-readRDS("Zhang_21_whole_cells_RNA.rds")

#Subset dataset por cluster para quedarme solo con las celulas B

pbmc_sub<-subset(pbmc, idents = c("3", "6", "7"))
rm(pbmc)
pbmc_sub<-SCTransform(pbmc_sub)
pbmc_sub<-RunPCA(pbmc_sub, dims= 1:20)
pbmc_sub<-RunUMAP(pbmc_sub, dims= 1:20)
pbmc_sub<-FindNeighbors(pbmc_sub, dims= 1:20)
pbmc_sub<-FindClusters(pbmc_sub, resolution = 0.1)
plot1<-DimPlot(pbmc_sub, reduction = "umap")
ggsave("UMAP_Subcluster_3_6_7_B_Cells_20_Dim.png" ,plot = plot1 ,height = 8.27, width =11.69)

saveRDS(pbmc_sub, "Subcluster_3_6_7_B_Cells_20_Dim.rds")


###Todas las celulas
pbmc<-readRDS("Zhang_21_whole_cells_RNA.rds")
#Subset dataset por cluster para quedarme solo con las celulas T

pbmc_sub<-subset(pbmc, idents = c("0", "1", "2"))
rm(pbmc)
saveRDS(pbmc_sub, file = "subset_Cluster_0_1_2.rds")
gc()
pbmc_sub<-SCTransform(pbmc_sub)
gc()
pbmc_sub<-RunPCA(pbmc_sub, dims= 1:20)
pbmc_sub<-RunUMAP(pbmc_sub, dims= 1:20)
pbmc_sub<-FindNeighbors(pbmc_sub, dims= 1:20)
pbmc_sub<-FindClusters(pbmc_sub, resolution = 0.1)
DimPlot(pbmc_sub, reduction = "umap")
plot1<-DimPlot(pbmc_sub, reduction = "umap")
ggsave("UMAP_Subcluster_0_1_2_T_Cells_20_Dim.png" ,plot = plot1 ,height = 8.27, width =11.69)
saveRDS(pbmc_sub, "Subcluster_0_1_2_T_Cells_20_Dim.rds")




###Todas las celulas
pbmc<-readRDS("Zhang_21_whole_cells_RNA.rds")

#Subset dataset por cluster para quedarme solo con los macrofagos
pbmc_sub<-subset(pbmc, idents = c("4", "5"))
rm(pbmc)
pbmc_sub<-SCTransform(pbmc_sub)
pbmc_sub<-RunPCA(pbmc_sub, dims= 1:20)
pbmc_sub<-RunUMAP(pbmc_sub, dims= 1:20)
pbmc_sub<-FindNeighbors(pbmc_sub, dims= 1:20)
pbmc_sub<-FindClusters(pbmc_sub, resolution = 0.1)
plot1<-DimPlot(pbmc_sub, reduction = "umap")
ggsave("UMAP_Subcluster_4_5_Macrofagos_20_Dim.png" ,plot = plot1 ,height = 8.27, width =11.69)
saveRDS(pbmc_sub, "Subcluster_4_5_Macrofagos_20_Dim.rds")


