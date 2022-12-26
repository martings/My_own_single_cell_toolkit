#########Instalacion requerida de instalar monocle3##########
#sudo apt-get install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")
#https://rpubs.com/mahima_bose/Seurat_and_Monocle3_p


BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                       'HDF5Array', 'terra', 'ggrastr'), force = T)

install.packages("devtools")
library(devtools)

devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("satijalab/seurat-wrappers")
######### Cargo paquetes y datos pre analizados de Seurat
library(monocle3)
library(Seurat)
#library(SeuratDisk)
#library(SeuratWrappers)

setwd("/media/megas/0123-4567/Fededa_Lab/")

pbmc<-readRDS("RDS/Subcluster_4_5_Macrofagos_20_Dim.rds")

#########Pre subset si quiero hacer un zoom in en Seurat previo#####
######En # estan los codigos necesarios para hacer un subset, ya sea por Cell_Type_pre_establecido
#levels(pbmc$Cell_Type)
#pbmc_sub<-subset(pbmc, ident = c("Basal_Basal", "Basal_Myoepithelial", "Stem", "Inmuno", "Luminal_HormSensProg", "Luminal_HormSensInt", "Luminal_HormSensDif"))
#rm(pbmc)
#pbmc_sub[["RNA"]] <- pbmc_sub[["integrated"]]

##Normalizacion y dimensiones reduccionales del subset hecho
#pbmc_sub <- SCTransform(pbmc_sub)
#pbmc_sub <- RunPCA(pbmc_sub)
#pbmc_sub <- RunUMAP(pbmc_sub, dims = 1:20)
#pbmc_sub <- FindNeighbors(pbmc_sub, dims = 1:20)
#pbmc_sub <- FindClusters(pbmc_sub, resolution = 0.1 )
#DimPlot(pbmc_sub, reduction = "umap", group.by = c("Newcluster", "Cell_Type"))
#DimPlot(pbmc_sub, reduction = "umap")

#DefaultAssay(pbmc_sub) <- "RNA"
#rm(pbmc)

#########Transformacion de los objetos de seurat a monocle3#####
try<-pbmc
rm(pbmc)
gc()
gene_annotation <- as.data.frame(rownames(try@reductions[["pca"]]@feature.loadings),
                                 row.names = rownames(try@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"


cell_metadata <- as.data.frame(try@assays[["RNA"]]@counts@Dimnames[[2]],
                               row.names = try@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"


New_matrix <- try@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(try@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

list_cluster <- try@active.ident
names(list_cluster) <- try@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-try@reductions[["umap"]]@cell.embeddings

#El PCA
gene_loadings<-try@reductions[["pca"]]@feature.loadings
cds_from_seurat@reduce_dim_aux@listData[["PCA"]]<-gene_loadings
rm(try)
gc()

##########Ya transformado a monocle3 hago las trayectorias#####
cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = F,verbose = T)

plot_cells(cds_from_seurat, 
           color_cells_by = 'cluster',
           label_groups_by_cluster=TRUE,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=2)

cds_from_seurat <-order_cells(cds_from_seurat)

genes_a_graph<- c("IFI30", "KRT17", "KRT9")
my_genes <- row.names(subset(fData(cds_from_seurat), gene_short_name %in% genes_a_graph)) 
cds_subset <- cds_from_seurat_one_branch[my_genes,]
plot_genes_in_pseudotime(cds_subset)


plot_cells(cds_from_seurat, 
           color_cells_by = 'pseudotime',
           label_groups_by_cluster=TRUE,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=2)


##########Analisis de branching si se quisiera######
cds_from_seurat_one_branch<-choose_graph_segments(cds_from_seurat)
cds_from_seurat_one_branch<-preprocess_cds(cds_from_seurat_one_branch)
cds_from_seurat_one_branch<-reduce_dimension(cds_from_seurat_one_branch, reduction_method = "UMAP")
cds_from_seurat_one_branch<-cluster_cells(cds_from_seurat_one_branch, reduction_method = "UMAP")
cds_from_seurat_one_branch <- learn_graph(cds_from_seurat_one_branch, use_partition = T,verbose = T)

cds_from_seurat_one_branch <-order_cells(cds_from_seurat_one_branch)

plot_cells(cds_from_seurat_one_branch, 
           color_cells_by = 'pseudotime',
           label_groups_by_cluster=TRUE,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=2)



saveRDS(pbmc_sub, "Subset_.rds")
saveRDS(cds_from_seurat, "RDS/CDS_Monocle_Macrofagos.rds")
